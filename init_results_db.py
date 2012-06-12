#!/usr/bin/env python

""" Sets up results.db SQLite database for enrichment analysis. """

import sqlite3
import os

schema = """
drop table if exists results;
create table results (
    _id     text    primary key,
    goid    text,
    term    text,
    pval    number,
    dataset text,
    factor  text,
    subset  text,
    year    number
);

drop table if exists datasets;
create table datasets (
    published   number,
    id          text    primary key
);
"""

db_filename = 'results.db'
db_exists = os.path.exists(db_filename)

def initdb(force_reinit=False):
    if db_exists and not force_reinit:
        print("Database exists. To reinitialize db, add 'force' to the command line args.")
        return

    with sqlite3.connect(db_filename) as conn:
        print("Initializing database with schema:")
        print schema
        conn.executescript(schema)

if __name__ == '__main__':
    import sys
    force_reinit = 'force' in sys.argv
    initdb(force_reinit)

