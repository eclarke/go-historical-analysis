DELIMITER //
CREATE PROCEDURE SelectTopShuffled (IN inpct float)
BEGIN
    SELECT dataset, subset, term, pval, qval
        FROM shuffled WHERE pct == inpct 
        ORDER BY dataset, subset, term, pval
        LIMIT 10;
END //

CREATE PROCEDURE SelectTopAcrossAll ()
BEGIN
    SELECT year, dataset, subset, term, pval, qval
        FROM results, 
            (SELECT dataset, subset, term, pval 
                FROM results WHERE year == 2012 
                ORDER BY pval LIMIT 10) as top
        WHERE dataset like top.dataset 
            AND subset like top.subset
            AND term like top.term
        ORDER BY dataset, subset, term, year
        LIMIT 10
END //
DELIMITER ;
