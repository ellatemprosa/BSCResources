*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter3/nconprob.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: N for a confidence interval for the difference between    *;
*          probabilities for 2 groups.                               *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
*
*  P1 AND P2 ARE THE PROBABILITIES;
*         P2 = D*P1;
*  Q1 AND Q2 ARE THE SAMPLE FRACTIONS, = 0.5 FOR EQUAL SAMPLES;
*  A = ERROR CONFIDENCE LEVEL;
*  E THE DESIRED LEVEL OF PRECISION AS A FRACTION OF P1, I.E. IF
     P1 = .4 AND P2 = .3 AND E = .2, THEN THE RESULTING N WILL PROVIDE
     1 - A CONFIDENCE THAT THE TRUE DIFFERENCE IS +/- E*P1 = .08;
 
title 'N for confidence interval for the difference between probabilities for 2 groups';
 
DATA ONE;
INPUT P1 D  Q1 Q2 A E;
CARDS;
.5    .8   .5 .5 .10 0.2
.5   1.0   .5 .5 .10 0.2
.45   .8   .5 .5 .10 0.2
.45  1.0   .5 .5 .10 0.2
.55   .8   .5 .5 .10 0.2
.55  1.0   .5 .5 .10 0.2
;
DATA TWO; SET ONE;
P2 = P1*D;
ERROR=E*P1;
ZA = PROBIT(A/2);
N = (ZA/ERROR)**2;
N = N*((P1*(1-P1)/Q1) + (P2*(1-P2)/Q2));
N1 = N*Q1; N2 = N*Q2;
PROC PRINT; VAR P1 D P2 E ERROR  A ZA Q1 Q2 N N1 N2;
RUN;
