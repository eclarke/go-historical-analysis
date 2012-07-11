
drop table if exists results;
create table results (
       ontology	     char(2),
       goid	     char(10),
       term	     text,
       pval	     double,
       qval	     double,
       dataset	     char(7),
       factor	     text,
       subset	     text(50),
       year	     year(4),
       num_annos     int,
       num_genes     int,
       anno_max	     int,
       anno_min	     int,
       min_depth     int,
       max_depth     int,
       filter_similar boolean,
       filter_size    boolean,
       filter_depth   boolean,
       shuffled	     float,
       primary key (dataset, subset(50), year, ontology, goid, filter_similar, filter_size, filter_depth)
);

-- after everything has finished processing, 
-- recommend creating indexes around pval, qval, etc