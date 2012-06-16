drop table if exists results;
create table results (
       _id   	     text    primary key,
       ontology	     text,
       goid	     text,
       term	     text,
       pval	     number,
       qval	     number,
       dataset	     text,
       factor	     text,
       subset	     text,
       year	     number,
       num_annos     number,
       num_genes     number,
       anno_max	     number,
       anno_min	     number
);