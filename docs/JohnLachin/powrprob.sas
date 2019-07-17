*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter3/powrprob.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: Power for the test of 2 proportions over a specified      *;
*          range of control group probabilities and relative risks.  *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
*
*  P1L AND P1U ARE THE RANGES OF CONTROL PROBABILITIES;
*  DL AND DU ARE THE RANGES OF RELATIVE RISKS WITH TREATMENT SUCH THAT
*         P2 = D*P1;
*  PINC IS THE SIZE OF THE INCREMENTS OVER WHICH P1 IS EVALUATED
*  DINC IS THE SIZE OF THE INCREMENTS OVER WHICH D  IS EVALUATED
*  Q1 AND Q2 ARE THE SAMPLE FRACTIONS, = 0.5 FOR EQUAL SAMPLES;
*  A = TYPE I ERROR LEVEL ONE OR TWO-SIDED, E.G. INPUT 0.025 FOR A
     TWO-SIDED 0.05 TEST
*  N = TOTAL SAMPLE SIZE
;
 
DATA ONE;
INPUT P1L P1U PINC DL DU DINC Q1 N A;
CARDS;
 
0.02 0.03 0.0025 1.5 3.5 0.25 0.5 609 0.025
;
DATA TWO; SET ONE;
  ZA = PROBIT(1-A);
  Q2 = 1 - Q1;
  N1 = Q1 * N;
  N2 = Q2 * N;
 
DO P1 = P1L TO P1U BY PINC;
 
DO D = DL TO DU BY DINC;
  P2 = D*P1;
  DELTA=ABS(P1-P2);
  PBAR = Q1*P1 + Q2*P2;
  SIG0 = PBAR*(1-PBAR)*( (1/Q1) + (1/Q2) );
  SIG0 = SQRT(SIG0);
  SIG1 = ( P1*(1-P1)/Q1) + ( P2*(1-P2)/Q2);
  SIG1 = SQRT(SIG1);
  ZB = ( (DELTA*SQRT(N)) - (ZA*SIG0) ) / SIG1;
  POWER = PROBNORM(ZB);
OUTPUT;
END;
END;
 
PROC PRINT; VAR P1 D P2 DELTA A N Q1 ZB POWER;
RUN;
