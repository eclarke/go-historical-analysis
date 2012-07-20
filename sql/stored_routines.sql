DELIMITER //

DROP PROCEDURE IF EXISTS SelectTopForYear //
CREATE PROCEDURE SelectTopForYear(IN _ds CHAR(10), IN _sub TEXT, IN _onto CHAR(2), IN _year YEAR(4), IN _max INT)
BEGIN
    SELECT * FROM results WHERE dataset like _ds AND subset like _sub AND ontology like _onto AND year like _year
    ORDER BY pval limit _max;
END//
    

DROP PROCEDURE IF EXISTS SelectTopShuffled //
CREATE PROCEDURE SelectTopShuffled (IN inpct float)
BEGIN
    SELECT *
        FROM shuffled WHERE pct = inpct 
        ORDER BY dataset, subset, term, pval
        LIMIT 10;
END //

DROP PROCEDURE IF EXISTS SelectTopAcrossAll //
CREATE PROCEDURE SelectTopAcrossAll (IN indset CHAR(10), IN insub TEXT, IN _onto CHAR(2))
BEGIN
    SELECT year, results.dataset, results.subset, results.term, results.pval, qval
        FROM results, 
            (SELECT dataset, subset, term, pval 
                FROM results WHERE year = 2012 AND dataset like indset AND subset like insub AND ontology like _onto
                ORDER BY pval LIMIT 10) as top
        WHERE results.dataset like top.dataset 
            AND results.subset like top.subset
            AND results.term like top.term
        ORDER BY results.dataset, results.subset, results.term, year;
END //

DROP PROCEDURE IF EXISTS SelectTopShuffledAcrossAll //
CREATE PROCEDURE SelectTopShuffledAcrossAll (IN indset CHAR(10))
BEGIN
    SELECT year, results.dataset, results.subset, results.term, results.pval, qval
        FROM results, 
            (SELECT dataset, subset, term, pval 
                FROM shuffled WHERE year = 2012 
		     AND dataset like indset 
		     AND pct = inpct
                ORDER BY pval LIMIT 10) as top
        WHERE results.dataset like top.dataset 
            AND results.subset like top.subset
            AND results.term like top.term
        ORDER BY results.dataset, results.subset, results.term, year;
END //

DELIMITER ;
