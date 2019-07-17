*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter3/npwrprob.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: Sample size based on the power function for the test of   *;
*          2 proportions.                                            *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
*
*  P1 AND P2 ARE THE PROBABILITIES;
*         P2 = D*P1;
*  Q1 AND Q2 ARE THE SAMPLE FRACTIONS, = 0.5 FOR EQUAL SAMPLES;
*  A = TYPE I ERROR LEVEL ONE OR TWO-SIDED, E.G. INPUT 0.025 FOR A
     TWO-SIDED 0.05 TEST
*  POWER = POWER;
;
 
title 'Sample Size Based On Power For Test Of 2 Proportions';
DATA ONE;
INPUT P1 D  Q1 Q2 A POWER;
CARDS;
0.40 0.70 0.5 0.5 0.025 0.90
;
 
** Problem 3.3 (incomplete);
0.2 0.05  .5  .5  0.05 0.80
0.2 0.10  .5  .5  0.05 0.80
0.2 0.05  .5  .5  0.05 0.90
0.2 0.10  .5  .5  0.05 0.90
;
 
DATA TWO; SET ONE;
P2 = P1*D;
DELTA=P1-P2;
ZA = PROBIT(1-A);
ZB = PROBIT(POWER);
PBAR = Q1*P1 + Q2*P2;
SIG0 = PBAR*(1-PBAR)*( (1/Q1) + (1/Q2) );
SIG0 = SQRT(SIG0);
SIG1 = ( P1*(1-P1)/Q1) + ( P2*(1-P2)/Q2);
SIG1 = SQRT(SIG1);
N = ( (ZA*SIG0) + (ZB*SIG1) ) / DELTA;
N = N**2;
N1 = Q1 * N;
N2 = Q2 * N;
PROC PRINT; VAR P1 D P2 DELTA PBAR A ZA ZB POWER N Q1 N1 Q2 N2;
RUN;
