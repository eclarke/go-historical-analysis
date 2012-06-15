#!/usr/bin/env python

""" Sets up results.db SQLite database for enrichment analysis. """

import sqlite3
import os

with open('results_db_schema.sql') as schema_file:
    schema = schema_file.read()

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

