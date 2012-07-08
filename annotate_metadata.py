import MySQLdb as my
import urllib
import urllib2

import __init__ as geo

from contextlib import closing

from ConfigParser import ConfigParser
from Bio import Entrez

config = ConfigParser()
config.read('settings.cfg')

# MySQL settings
MYUSER = config.get('MySQL', 'user')
MYHOST = config.get('MySQL', 'host')
MYPASS = config.get('MySQL', 'pass')
MYDB = config.get('MySQL', 'db')

# Entrez settings
Entrez.email = config.get('Entrez', 'email')

# Annotator settings
API_KEY = config.get('Annotator', 'apikey')
ANNO_EMAIL = config.get('Annotator', 'email')
ANNO_URL = 'http://rest.bioontology.org/obs/annotator'

DB = None

def dbconnect():
    return my.connect(MYHOST, MYUSER, MYPASS, MYDB)


def get_datasets(pval=0.01):
    db = dbconnect()
    with closing(db.cursor()) as c:
        c.execute('SELECT distinct dataset FROM results WHERE pval < %s', (pval,))
        datasets = [x[0] for x in c.fetchall()]
    db.close()
    return datasets


def insert_metadata(datasets):
    fetch = lambda x: geo.fetch(x, destdir='data/')
    with closing(dbconnect()) as db:
        with closing(db.cursor()) as c:
            for ds in datasets:
                dataset = fetch(ds)
                did = dataset.id
                desc = dataset.meta['description']
                pmid = dataset.meta['pubmed_id'] if 'pubmed_id' in dataset.meta else None
                print "fetching abstract..."
                abstract = Entrez.efetch(id=pmid, db='pubmed', rettype='abstract', retmode='text').read() if pmid else None
                print "inserting metadata..."
                c.execute("REPLACE INTO metadata VALUES(%s, %s, %s, %s);", (did, desc, pmid, abstract))
                print "fetching annotations..."
                desc_annos = annotate(desc)
                
                abst_annos = annotate(abstract.replace('\n',' ')) if abstract else None
                print "inserting annotations...",
                desc_sql = "REPLACE INTO annotations(dataset, goid, term, source) values('%s', %%s, %%s, '%s')" % (did, 'description')
                c.executemany(desc_sql, desc_annos)
                if abstract:
                    abst_sql = "REPLACE INTO annotations(dataset, goid, term, source) values('%s', %%s, %%s, '%s')" % (did, 'abstract')
                    c.executemany(abst_sql, abst_annos)
                db.commit()
                print "done.\n"

            

def annotate(text):
    GO_local_id = 46440
    GO_virtual_id = 1070

    params = {
        'longestOnly':'false',
        'wholeWordOnly':'true',
        'withContext':'true',
        'filterNumber':'true', 
        'stopWords':'',
        'withDefaultStopWords':'true', 
        'isStopWordsCaseSenstive':'false', 
        'minTermSize':'1', 
        'scored':'true',  
        'withSynonyms':'true', 
        'ontologiesToExpand':'',   
        'ontologiesToKeepInResult': GO_virtual_id,   
        'isVirtualOntologyId':'true', 
        'semanticTypes':'',
        'levelMax':'1',
        'mappingTypes':'null', 
        'textToAnnotate': text, 
        'format':'tabDelimited',
        'apikey':API_KEY,
        }

    submitUrl = '%s/submit/%s' % (ANNO_URL, ANNO_EMAIL)
    pd = urllib.urlencode(params)
    fh = urllib2.urlopen(submitUrl, pd)
    results = [x.strip('\n').split('\t') for x in fh.readlines()]
    terms = list(set([(row[1].split('/')[1], row[2]) for row in results]))
    return terms
    


