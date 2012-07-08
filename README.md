go-historical-analysis
======================

Pipeline for the mass analysis of GEO datasets using historical versions of the Gene Ontology annotations.

Dependencies
------------
0. Python >= 2.6
1. numpy 
2. scipy (be sure to run scipy.test('full') and ensure it passes!)
3. statsmodels (for the multiple testing FDR)
4. MySQLdb (in pip as mysql-python) [optional- can use SQLite but do not run concurrent processes!]


Usage
------

1.  Creation of the flattened ontology files:

    We used the CVS access point provided by the GO Consortium to find historical versions of the
    GO structure (.owl file) and UniProtKB human gene annotations. These loosely follow releases
    in mid-May of each year going back until 2004. 

    The .owl file was then read in using the accompanying JAR file (go-flattener.jar) outputted
    in a tab-delimited format as follows:
    > [GO ID] [PARENT TERM 1] [PARENT TERM 2] ...
    where each parent term is defined by having an "is a" or "part of" relationship with the
    child. Consequently, terms with no parents are either obsolete or top-level terms and exist
    only on their line. The accompanying Annotations.py file contains functions to read in these
    flat files and return them as a Python dictionary object.

    In conjunction with the annotation file corresponding to that date, we create an annotation
    object that is represented as follows:
    > {'meta':
    >	'year':20XX,
    >  'anno':
    >	'GO:0001234':
    >	    {'name':'name of GO term',
    >	     'genes': [list of UniProt ids associated with this term or this term's children]}
    > 	'GO:0005678':
    >	    ...}
    > }
    This is stored as a JSON object, which can be found in data/goa-[year].json
	
2.  Selection of the datasets:
    
    We selected over 200 datasets from the NCBI GDS database using the following query:
    > ((("count"[Sample Value Type] OR "transformed count"[Sample Value Type])
    >   AND "disease state"[Subset Variable Type] AND "homo sapiens"[Organism])
    >  NOT "time"[Subset Variable Type]) 

3.  Pre-filter dataset probes:
    
    We 
    
3.  Detecting differentially-expressed genes:

    After preliminary filtering 

   
   




