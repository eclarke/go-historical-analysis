# Copyright 2012 by Erik Clarke. All rights reserved.

"""Classes to represent records for NCBI's Gene Expression Omnibus (GEO).

http://www.ncbi.nlm.nih.gov/geo/
"""


import re
from collections import defaultdict
from copy import deepcopy
import numpy
from numpy import array
import scipy.stats as stats
from statsmodels.stats import multitest


def _truncate(string, trunc=20):
    length = len(string)
    if length <= trunc:
        return string
    if trunc < 8:
        return '[...]'
    half = (trunc // 2) - 3

    p1 = ''.join([x for x in string[:half]])
    p2 = ''.join([x for x in string[length - half:]])
    return '[...]'.join((p1, p2))


def _pad(string, pad_to=20, justify='left'):
    length = len(string)
    if length >= pad_to:
        return string
    spacing = pad_to - length
    if justify == 'right':
        spaces = ''.join([' ' for x in xrange(spacing)])
        return spaces + string
    if justify == 'center':
        spaces_half = ''.join([' ' for x in xrange(spacing / 2)])
        if spacing % 2 == 0:
            return string.join((spaces_half, spaces_half))
        else:
            return string.join((spaces_half, spaces_half + ' '))
    if justify == 'left':
        spaces = ''.join([' ' for x in xrange(spacing)])
        return string + spaces


class SOFTRecord(object):
    """Base class for SOFT-format file representations (GDS, GSE, GSM, and GPL).

    SOFT file format structure:
    - Line starts indicate line types:
        '^': entity indicator line (PLATFORM, SAMPLE, SERIES, DATASET)
        '!': entity attribute line
        '#': data table header description line
        n/a: data table row

    - Line format:
        [line-type char] [label] = [value]

    Attributes:
        id:      entity id (GEO accession)
        type:    entity type
        meta:    metadata information for the record {key:value}
        source:  location of actual SOFT file (optional)
    """
    def __init__(self, type, id):
        self.id = id
        self.type = type
        self.source = None
        self.meta = {}

    def __repr__(self):
        return self.id

    def full(self):
        """Print all information associated with this SOFTRecord."""
        self.print_metadata()

    def print_metadata(self):
        """Prints the metadata associated with this data."""
        print("*** %s: %s ***" % (self.type, self.id))
        for entry in self.meta:
            print("%s:\t%s" % (entry, self.meta[entry]))


class Record(SOFTRecord):
    """Represents two basic GEO record types: samples (GSM), and platforms (GPL).

    GEO Datasets and Series are special records that have their own classes.

    Attributes:
        id:      entity id (GEO accession)
        type:    entity type
        meta:    metadata information for the record {key:value}
        table:   list of rows, split by tabs. First row is the header.
                 [header, row1, row2, ...]
        columns: column names and descriptions.
                 {colname: description}
    """

    def __init__(self, type, id):
        super(Record, self).__init__(type, id)
        self.table = []
        self.columns = {}

    def print_table(self, length=20, max_width=80, col_width=20):
        """Prints a nicely-formatted overview of the data table.

        Arguments: (set any to < 0 to disable)
        length      number of rows to print (omits middle)
        max_width:  maximum width of the table (omits ending columns if needed)
        col_width:  width of each column (trims middle if needed)
        """
        num_rows = len(self.table)
        printed_ellipses = False
        for c, row in enumerate(self.table):
            # if we're within the first or last part of the table
            # (or displaying all rows)
            if (length < 0) or (c < length // 2) or (c > num_rows - length // 2):
                line = ''
                i = 0

                def fill(x):
                    _pad(_truncate(x, trunc=col_width), pad_to=col_width)

                # set the column width if width is not -1
                entry = fill(row[i]) if col_width > 0 else row[i]
                # build the row, minding the max_width (if not -1)
                while ((i < len(row)) and ((max_width < 0)
                       or (i < len(row) and len(line) + len(entry) < max_width))):
                    entry = fill(row[i]) if col_width > 0 else row[i]
                    line = "%s\t%s" % (line, entry)
                    i += 1
                if i < len(row) - 1:
                    line += "\t [...]"
                print line
            elif not printed_ellipses:
                print("\t[...]")
                printed_ellipses = True

    def print_columns(self):
        """Prints the column information associated with the data table."""
        for i,col in enumerate(self.table[0]):
            print("[%d] %s:\t%s" % (i, col, self.columns[col]['description']))

    def full(self):
        self.print_metadata()
        self.print_table()
        self.print_columns()


class Dataset(Record):
    """Represents a GEO Dataset (GDS) record.

    Attributes:
        id:      entity id (GEO accession)
        type:    entity type
        meta:    metadata information for the record {key:value}
        table:   list of rows, split by tabs. First row is the header.
                 [header, row1, row2, ...]
        columns: column names, descriptions, and factor subsets.
                 {colname: {description:'', factor1:'', factor2:'', ...}}
        factors: factors and their subsets
                 {factor1: [subset1, subset2,...], factor2:[...], ...}
        _matrix: data table as a numpy matrix for numerical analysis
    """
    def __init__(self, id):
        super(Dataset, self).__init__("DATASET", id)
        self.factors = defaultdict(dict)
        self._platform = None
        self._matrix = None

    def print_columns(self, factor=None, subset=None):
        """Readable display of data columns, with associated description,
        factor, and subset.

        Arguments:
            factor:     restrict to only display this factor's samples
            subset:    restrict to only display this subset's samples
        """
        if subset and not factor:
            factors = [factor for factor in self.factors
                        if (subset in self.factors[factor])]
        else:
            factors = sorted(self.factors.keys()) if not factor else [factor]
        print '\t'.join(['NAME', 'DESCRIPTION'] +
                        [factor.upper() for factor in factors])
        # we reference the header, as it's in order and dicts are not
        for col_name in self.table[0]:
            column = self.columns[col_name]
            if col_name.startswith('GSM'):
                print '\t'.join([col_name, column['description']] +
                                [column[factor] for factor in factors])
            else:
                print '\t'.join((col_name, column['description']))

    def matrix(self, refresh=False, omitNulls=True, nullVal='null'):
        """Returns a numpy matrix of values built from the dataset's table.

        This takes some time to return the first time it is run, but results are
        cached for later use.

        Arguments:
            refresh:    rebuild the matrix (if data values changed)
                        [default = False]
        """
        # memoize the matrix (as it takes some time to generate)
        if self._matrix != None and not refresh:
            return self._matrix

        def tofloat(x):
            try:
                return float(x)
            except ValueError:
                return float('nan')
        if omitNulls:
            print "converting to numpy array and omitting probes containing " + nullVal
        header = self.table[0]
        samples = [x for x in header if re.match(r'GSM\d+$', x)]
        end_samples = header.index(samples[len(samples) - 1])
        matrix = []
        for row in self.table[1:]:
            values = row[2:end_samples + 1]
            if omitNulls and nullVal in values:
                continue
            matrix.append([tofloat(x) for x in values])
        self._matrix = array(matrix)
        return self._matrix

    def to_numeric(self):
        return NumericDataset(self)


class Series(SOFTRecord):
    """Represents a GEO Series (GSE) record.

    Attributes:
        id:      entity id (GEO accession)
        type:    entity type
        meta:    metadata information for the record {key:value}
        platforms:  list of platform records associated with series
                    [GPL0, GPL1, ...]
        samples: list of sample records associated with the series
                 [GSM0, GSM1, ...]
    """

    def __init__(self, id):
        super(Series, self).__init__("SERIES", id)
        self.platforms = []
        self.samples = []


class NumericDataset(SOFTRecord):
    """Represents a special form of a GEO Dataset that can be more easily
    be used in statistical and numerical analysis.

    Attributes:
        meta:       the original metadata information from the dataset
        matrix:     a numpy matrix of the probe values, omitting non-sample cols
        probes:     a numpy array of the probes, corresponding to rows in the
                    matrix
        header:     a numpy array of the header row from original table
        factors:    the original factor information from the dataset
    """

    def __init__(self, dataset=None):
        super(NumericDataset, self).__init__("NUMERIC DATASET", dataset.id)
        self.matrix = deepcopy(dataset.matrix())
        self.probes = array([x[:2] for x in dataset.table[1:]])
        self.header = array(dataset.table[0])
        self.factors = deepcopy(dataset.factors) if dataset else None
        self.meta = deepcopy(dataset.meta) if dataset else None
        self._log2xformed = False
        self._filtered = False
        # heuristic for determining if already log2 transformed:
        print "checking if matrix has been normalized..."
        if dataset.meta['value_type'] == 'count':
            print "matrix not normalized, taking the log2 of each value..."
            self.log2xform()
        else:
            print "matrix data has been modified, doing no further normalization..." 
            

        """ Old algorithm for checking if normalized
        sm = numpy.sort(self.matrix, axis=None)
        l = len(sm) / 100
        a = sm[::-1][l]
        if a > 100:  # reverse and then choose 1st percentile value (top hundredth)
            print "value of 1st percentile matrix value above 100 (%f); log2 transforming" % a
            self.log2xform()
        """

    def log2xform(self):
        """Returns this dataset with the binary log applied to each value in
        the data matrix.

        This operation is performed in-place for performance reasons.
        """
        if not self.log2xformed():
            self.matrix = numpy.log2(self.matrix)
            self._log2xformed = True
        return self

    def log2xformed(self):
        """Returns True if the matrix has been log2 transformed since being
        imported as a NumericDataset."""
        return self._log2xformed

    def filter(self, fn=numpy.median):
        """Filters probes from the data matrix that have maximum values below
        the value returned from `fn(self.matrix)`. This fn is by default the
        matrix median (numpy.median).

        Returns the same dataset with the rows and probes below this value
        removed, so as to allow chaining of methods (such as
        dataset.filter().log2xform()).

        This operation is performed in-place for performance reasons.

        Arguments:
            fn:     a function that returns a single value for a matrix input.
        """
        if not self.filtered():
            level = fn(self.matrix)
            above = array([max(x) > level for x in self.matrix])
            print "Filter: removing %d/%d probes." % (len([x for x in above if not x]), len(self.matrix))
            self.matrix = self.matrix[above, :]
            self.probes = self.probes[above, :]
            self._filtered = True
        return self

    def filtered(self):
        return self._filtered

    def diffexpressed(self, _subset, _factor, fdr_limit, verbose=True):
        """Returns an array of probes that are differentially expressed according
        to the following method:

        1)  Perform an independent t-test on the probe values for the specified
            subset against the probe values for the non-subset samples.
        2)  To correct for multiple testing errors, calculate the
            Benjamini-Hochberg FDR q-value for each p-value.
        3)  Filter probes where the q-value is above the cutoff.

        Arguments:
            _subset:    the subset to test for expressed genes
            _factor:    the factor the subset belongs to
            fdr_limit:  the FDR q-value representing the upper limit for results
        """
        if not self.filtered():
            print("Warning: Finding differentially expressed genes on an unfiltered matrix may fail. Run dataset.filter().")

        matrix = self.matrix
        probes = self.probes
        samples = self.factors[_factor][_subset]

        inA = array([x in samples for x in self.header[2:]])
        A = numpy.transpose(matrix[:, inA])
        B = numpy.transpose(matrix[:, numpy.invert(inA)])

        t, pvals = stats.ttest_ind(A, B)
        rejected, qvals = multitest.fdrcorrection(pvals, alpha=qval_limit)

        # probe values are [probe_name, entrez_id] form (hence x[0])
        diffexp = [x[0] for i, x in enumerate(probes) if qvals[i] < qval_limit]
        if verbose:
            print("%d samples, %d differentially expressed genes in %s: %s" % (len([x for x in inA if x]), len(diffexp), _factor, _subset))
        return diffexp


    def diffexpressed_alt(self, _subset, _factor, pval_cutoff, d_avg_cutoff, verbose=True, more=False):
        """Returns an array of probes that are differentially expressed according
        to the following method:

        1)  Perform an independent t-test on the probe values for the specified
            subset against the probe values for the non-subset samples.
        2)  Measure the magnitude of difference between the mean value for
            each probe across the subset samples and the mean value for the
            non-subset samples.
        3)  Filter results from 1) and 2) based on specified cutoffs
        4)  Return the probes that were significant in both measures.

        Consider the effects of multiple comparisons when selecting the p-value
        cutoff.

        Arguments:
            _subset:    the subset to test for expressed genes
            _factor:    the factor the subset belongs to
            pval_cutoff:    the maximum p-value to consider significant
            d_avg_cutoff:   the minimum difference in magnitude
        """

        matrix = self.matrix
        probes = self.probes
        samples = self.factors[_factor][_subset]
        nprobes = range(len(matrix))

        inA = array([x in samples for x in self.header[2:]])
        A = numpy.transpose(matrix[:, inA])
        B = numpy.transpose(matrix[:, numpy.invert(inA)])
        mA = numpy.mean(A, axis=0)
        mB = numpy.mean(B, axis=0)

        t, pvals = stats.ttest_ind(A, B)
        # boolean arrays (T if significant, F otherwise) for the cutoffs
        sig_pvals = [x < pval_cutoff for x in pvals]
        sig_diffs = [abs(mA[i] - mB[i]) > d_avg_cutoff for i in nprobes]
        sig_union = array([sig_pvals[i] and sig_diffs[i] for i in nprobes])
        diffexp = probes[sig_union, :]
        if verbose:
            print("%s samples, %s differentially expressed genes in %s: %s" % (len([x for x in inA if x]), len(diffexp), _factor, _subset))
        if more:
            return diffexp, pvals, (A, B, mA, mB), (sig_pvals, sig_diffs, sig_union)
        else:
            return [x[0] for x in diffexp]

        
