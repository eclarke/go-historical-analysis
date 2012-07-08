
drop table if exists metadata;
create table metadata (
       dataset	      char(7),
       description    text,
       pmid	      text,
       abstract	      text,
       primary key (dataset),
       index (dataset)
);

drop table if exists annotations;
create table annotations (
       dataset		 char(7),
       goid		 char(10),
       term		 text,
       source		 text,
       constraint unique (dataset,goid),
       index (dataset)
);