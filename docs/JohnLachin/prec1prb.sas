*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter3/prec1prb.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: The precision for a confidence interval for a single      *;
*          probability evaluated over the range 0.5-1.0 with plots   *;
*          for a given sample size.                                  *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
*
*  P IS THE PROBABILITY, evaluated over the range 0.5 - 1.0;
*
*  A = ERROR CONFIDENCE LEVEL (ALPHA);
*
*  N IS THE SAMPLE SIZE
*
*  THE PROGRAM SOLVES FOR THE RESULTING PRECISION EXPRESSED AS
*
*  ERROR = THE WIDTH OF THE CONFIDENCE INTERVAL
*
*  E = THE LEVEL OF PRECISION AS A FRACTION OF P, I.E.
       THE WIDTH IS  +/- E*P;
*
*  Error and E are plotted as a function of P for the specified
   values of A and N;
 
title 'Precision for confidence interval for a single probability';
 
DATA ONE;
INPUT A N;
CARDS;
  .10  60
  .10  180
;
DATA TWO; SET ONE;
P=0.475;
 
DO I = 1 TO 20;
  P = P + 0.025;
  S = (P*(1-P))/N;
  ZA = - PROBIT(A/2);
  ERROR=ZA * SQRT(S);
  E = ERROR/P;
  OUTPUT;
END;
 
PROC SORT; BY N;
PROC PRINT; VAR P A ZA  N ERROR E;
PROC PLOT; PLOT ERROR*P = '*'
                E*P = 'e'
     / VAXIS = 0 TO .24 BY 0.02
       HAXIS = 0.5 TO 1.0 BY 0.05
       OVERLAY;
     BY N;
PROC PLOT; PLOT ERROR*P = '*'
     / VAXIS = 0 TO .24 BY 0.02
       HAXIS = 0.5 TO 1.0 BY 0.05
       OVERLAY;
     BY N;
PROC PLOT; PLOT E*P = 'e'
     / VAXIS = 0 TO .24 BY 0.02
       HAXIS = 0.5 TO 1.0 BY 0.05
       OVERLAY;
     BY N;
RUN;
