
-- drop table if exists results;
create table results (
       _id   	     varchar(32),
       ontology	     varchar(2),
       goid	     varchar(10),
       term	     text,
       pval	     double,
       qval	     double,
       dataset	     varchar(7),
       factor	     text,
       subset	     text,
       year	     year(4),
       num_annos     int,
       num_genes     int,
       anno_max	     int,
       anno_min	     int,
       primary key (_id),
       index (dataset)
);

-- after everything has finished processing, 
-- recommend creating indexes around pval, qval, etc