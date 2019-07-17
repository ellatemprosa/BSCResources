*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter3/ncon1pr.sas                  *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: Sample size N to provide a desired precision for a        *;
*          confidence interval for a single probability.             *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
*
*  P IS THE PROBABILITY;
*
*  A = ERROR CONFIDENCE LEVEL;
*
*  E THE DESIRED LEVEL OF PRECISION AS A FRACTION OF P, I.E.
     THE RESULTING N WILL PROVIDE
     1 - A CONFIDENCE THAT THE TRUE PROBABILITY IS +/- E*P;
 
title1 'n for confidence interval for a single probability';
DATA ONE;
INPUT P A E;
CARDS;
.5   .10 0.2
.45  .10 0.2
.55  .10 0.2
;
DATA TWO; SET ONE;
ERROR=E*P;
ZA = PROBIT(A/2);
N = (ZA/ERROR)**2;
N = N*(P*(1-P));
PROC PRINT; VAR P E ERROR  A ZA  N ;
RUN;
