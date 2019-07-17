*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter3/ssnprob.sas                  *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: Sample size based on power for the test of 2 proportions  *;
*          over a specified range of control group probabilities     *;
*          and relative risks.                                       *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
*
*  P1L AND P1U ARE THE RANGES OF CONTROL PROBABILITIES;
*  RL AND RU ARE THE RANGES OF relative risks R=p2/p1 such that
*         P2 = R*P1;
*  PINC IS THE SIZE OF THE INCREMENTS OVER WHICH P1 IS EVALUATED
*  RINC IS THE SIZE OF THE INCREMENTS OVER WHICH R  IS EVALUATED
*  Q1 AND Q2 ARE THE SAMPLE FRACTIONS, = 0.5 FOR EQUAL SAMPLES;
*  A = TYPE I ERROR LEVEL ONE OR TWO-SIDED, E.G. INPUT 0.025 FOR A
     TWO-SIDED 0.05 TEST
*  B = type 2 error
*  N = TOTAL SAMPLE SIZE
;
 
DATA ONE;
INPUT P1L P1U PINC rL rU rINC Q1 A B; *** note P1u > p1l and ru > rl;
CARDS;
0.22 0.34 0.04 0.60 0.67 0.07 0.5 0.0125 0.2
0.22 0.34 0.04 0.60 0.67 0.07 0.5 0.0125 0.1
;
 
DATA TWO; SET ONE;
  ZA = PROBIT(1-A);
  Zb = PROBIT(1-b);
  Q2 = 1 - Q1;
 
DO P1 = P1L TO P1U BY PINC;
 
DO r = rL TO rU BY rINC;
  P2 = r*P1;
  DELTA=ABS(P1-P2);
  PBAR = Q1*P1 + Q2*P2;
  SIG0 = PBAR*(1-PBAR)*( (1/Q1) + (1/Q2) );
  SIG0 = SQRT(SIG0);
  SIG1 = ( P1*(1-P1)/Q1) + ( P2*(1-P2)/Q2);
  SIG1 = SQRT(SIG1);
  n = ( (ZA*SIG0) + (ZB*SIG1)) / (DELTA);
  n = n**2;
  N1 = Q1 * N;
  N2 = Q2 * N;
OUTPUT;
END;
END;
 
PROC PRINT; VAR P1 r P2 DELTA A B Q1 n;
RUN;
