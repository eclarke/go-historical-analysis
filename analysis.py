import MySQLdb as my
import csv
import math
from collections import defaultdict
from contextlib import closing
from optparse import OptionParser
from ConfigParser import ConfigParser
import itertools

TOP_SQL = """
SELECT year, {table}.dataset, {table}.subset, {table}.term, {table}.pval, qval
        FROM {table}, 
            (SELECT dataset, subset, term, pval 
                FROM {table} WHERE year = %s AND dataset like %s AND subset like %s AND ontology like %s
                ORDER BY pval LIMIT 10) as top
        WHERE {table}.dataset like top.dataset 
            AND {table}.subset like top.subset
            AND {table}.term like top.term
        ORDER BY {table}.dataset, {table}.subset, {table}.term, year;
"""

RAND_SQL = """
SELECT year, {table}.dataset, {table}.subset, {table}.term, {table}.pval, qval
        FROM {table}, 
            (SELECT dataset, subset, term, pval 
                FROM {table} WHERE year = %s AND dataset like %s AND subset like %s AND ontology like %s AND term not 'n/a'
                ORDER BY rand() LIMIT 10) as top
        WHERE {table}.dataset like top.dataset 
            AND {table}.subset like top.subset
            AND {table}.term like top.term
        ORDER BY {table}.dataset, {table}.subset, {table}.term, year;
"""

SIGSET_SQL = """
SELECT dataset, subset 
    FROM {table} WHERE pval < 0.05
    GROUP BY subset;
"""


IN2012_SQL = """
CALL in2012(%s,%s,%s,%s);
"""


RANKTERMS_SQL = """
SELECT goid, term FROM {table} 
    WHERE dataset like %s and year like %s and subset like %s and ontology like %s
    ORDER BY pval
"""


def tofloat(x):
    """Safely converts input to a Float or None if impossible."""
    try:
        return float(x)
    except ValueError:
        return None


def mean(seq):
    """Returns arithmetic mean."""
    if seq:
        return sum(seq) / len(seq)
    else:
        return 0


def convergence(resultset):
    """Calculates the "convergence", or average slope, of a set of score vectors.

    A large negative number represents overall decrease in values, while a small number
    (close to zero) represents little to no change."""

    rs = []
    for row in resultset:
        rs.append([tofloat(x) for x in row if tofloat(x) and tofloat(x) < 1])
    t = [x[-1] - x[0] for x in rs if x]
    print t
    return mean(t)


def setDelta(resultset):
    deltas = delta(resultset)
    d = math.log10(mean(deltas)) if deltas else None
    return d
    

def delta(resultset):
    deltas = []
    for term, group in itertools.groupby(resultset, lambda x: x[3]):
        _pvals = []
        for year in group:
            _pvals.append(year[4])
        d = _pvals[0] / _pvals[-1]
        deltas.append(d)
    if len(deltas) == 10:
        return deltas
    else:
        return None


def getIn2012All(ds, ss, onto, conn, table='results'):
    percentages = []
    for year in range(2004,2013):
        with closing(conn.cursor()) as c:
            c.execute(IN2012_SQL, (ds,ss,onto,year))
            percentages.append(c.fetchone())
    return percentages


def getTopOverTime(ds, subset, ontology, conn, table='results', year=2012):
    """Returns the top terms for a given year across all years in which they
    were recorded. The returned result set is in the form (year, dataset,
    subset, term, pval, qval (if available))"""
    with closing(conn.cursor()) as c:
        c.execute(TOP_SQL.format(table=table), (year, ds, subset, ontology))
        return c.fetchall()


def getRandOverTime(ds, subset, ontology, conn, table='results', year=2012):
    """Returns the top terms for a given year across all years in which they
    were recorded. The returned result set is in the form (year, dataset,
    subset, term, pval, qval (if available))"""
    with closing(conn.cursor()) as c:
        c.execute(RAND_SQL.format(table=table), (year, ds, subset, ontology))
        return c.fetchall()
   

def getRankedTerms(ds, ss, onto, year, conn, table='results'):
    with closing(conn.cursor()) as c:
        c.execute(RANKTERMS_SQL.format(table=table), (ds, year, ss, onto))
        return c.fetchall()


def findRankOverTime(ds, ss, onto, root_year, conn, table='results', top=10):
    root_terms = getRankedTerms(ds, ss, onto, root_year, conn, table)
    size = len(root_terms)
    top_terms = root_terms[:top]
    terms = {}
    pctiles = defaultdict(list)
    for year in range(2004,2013):
        print "Getting all terms for year %d" % year
        terms[year] = root_terms if year == root_year else getRankedTerms(ds, ss, onto, year, conn, table)
        year_size = len(terms[year])
        for term,goid in top_terms:
            rank = terms[year].index((term,goid)) if (term,goid) in terms[year] else None
            pctile = rank/year_size if rank else 0
            pctiles[year].append((term,pctile,rank,year_size))
    return pctiles
 

def organize(resultset, header=True):
    """Organizes a resultset from getTopOverTime into a table where the headers
    are years (2004-2012) and the rows are the scores for each term over the 
    years."""
    res = defaultdict(list)
    for row in resultset:
        res[row[3]].append((row[0], row[-2:]))

    rows = [['',]+range(2004,2013)] if header else []
    
    for key,val in res.iteritems():
        tmp = [key,]
        for i in xrange(2004, 2013):
            _score = [x[1][0] for x in val if x[0] == i]
            score = _score[0] if _score else ''
            tmp.append(score)
        rows.append(tmp)
    return rows


def getSigDatasets(conn, table='results'):
    """Returns datasets+subsets that have terms with pvals less than 0.05"""
    with closing(conn.cursor()) as c:
        print("Fetching significant subsets...")
        c.execute(SIGSET_SQL.format(table=table))
        res = c.fetchall()
        print("Found %d significant subsets..." % len(res))
        return res


def calculateAllDeltas(year, ontology, conn, table='results', rand=False, expanded=False):
    sigds = getSigDatasets(conn)
    resultsets = []
    for ds, ss in sigds:
        print "Fetching top terms from %s for %s::%s... " % (year, ds, ss),
        if rand:
            rs = getRandOverTime(ds, ss, ontology, conn, table, year)
        else:
            rs = getTopOverTime(ds, ss, ontology, conn, table, year)
        d = delta(rs) if expanded else setDelta(rs)
        print d
        if d:
            resultsets.append([ds, ss, d, year])
    return resultsets


def main():
    parser = OptionParser(usage="%prog [options] year outfile")
    parser.add_option("--expanded", action='store_true', dest='expanded', default=False)
    parser.add_option("--rand", action='store_true', dest='rand', default=False)
    parser.add_option("--config", action='store', dest='config', default='configs/settings.cfg')
    parser.add_option("-o", action='store', dest='ontology', type='choice', choices=['MF','CC','BP'], default='BP')
    parser.add_option("--table", action='store', dest='table')
    parser.add_option("--db", action='store', dest='db')
    opts, args = parser.parse_args()
    
    config = ConfigParser()
    config.read(opts.config)
    
    if not args:
        parser.print_help()
        parser.error('Missing required arguments')

    year, outf = args

    ontology = opts.ontology
    table = opts.table if opts.table else config.get('MySQL', 'table')
    db = opts.db if opts.db else config.get('MySQL','db')
    user = config.get('MySQL', 'user')
    passwd = config.get('MySQL','pass')
    host = config.get('MySQL', 'host')
    rand = opts.rand

    conn = my.connect(host, user, passwd, db)
    results = calculateAllDeltas(year, ontology, conn, table, rand, opts.expanded)
    print "Average delta: ", mean([x[2] for x in results])
    with open(outf, 'wt') as out:
        writer = csv.writer(out)

        if opts.expanded:
            writer.writerow(['DATASET','SUBSET']+range(1,11)+['ROOT_YEAR',])
            for row in results:
                row = row[0:2]+list(row[2])+[row[3],]
                writer.writerow(row)
        else:
            [writer.writerow(row) for row in results]


if __name__ == '__main__':
    main()
