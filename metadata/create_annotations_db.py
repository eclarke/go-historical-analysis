import sqlite3

annocommands = ["DROP TABLE if exists {table_name};",
"""CREATE TABLE {table_name} (
       db, 
       id, 
       sym, 
       qualifier, 
       goid, 
       dbref, 
       evcode, 
       withorfrom, 
       aspect, 
       name, 
       synonyms, 
       type, 
       taxon, 
       date, 
       assignedby, 
       direct,
       unique(id,goid)
);""",
"CREATE INDEX {table_name}_goid_idx ON {table_name} (goid);"]

ontocommands = ["DROP TABLE if exists {table_name};",
"""CREATE TABLE {table_name} (term, ancestor, unique(term, ancestor));""",
"CREATE INDEX {table_name}_term_idx ON {table_name} (term);"]


def main(db, _range, opt):
    assert opt in ['anno','go']
    c = sqlite3.connect(db)
    for year in _range:
        tbl_prefix = opt
        name = tbl_prefix + str(year)
        commands = ontocommands if opt == 'go' else annocommands
        for sql in commands:
            c.execute(sql.format(table_name=name))
            c.commit()
        print "Created table for year ", year
    c.close()


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 4:
        print "Usage: python create_annotations_db.py <start year> <end year> <db name> <anno OR go>"
        print "       anno: create annotation tables (drops any existing tables)"
        print "       go: create gene ontology tables (term:ancestor) (drops any existing tables)"
    else:
        _range = range(int(sys.argv[1]), int(sys.argv[2])+1)
        main(sys.argv[3], _range, sys.argv[4])
