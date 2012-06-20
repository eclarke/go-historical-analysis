
drop table if exists results;
create table results (
       _id   	     char(32),
       _subid	     char(32),
       ontology	     char(2),
       goid	     char(10),
       term	     text,
       pval	     double,
       qval	     double,
       dataset	     char(7),
       factor	     text,
       subset	     text,
       year	     year(4),
       num_annos     int,
       num_genes     int,
       anno_max	     int,
       anno_min	     int,
       primary key (_subid),
       index (_id),
       index (dataset)
);

-- after everything has finished processing, 
-- recommend creating indexes around pval, qval, etc