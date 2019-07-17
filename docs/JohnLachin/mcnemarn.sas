*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter5/mcnemarn.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: Sample size calculations for McNemar's test for paired    *;
*          or matched 2x2 tables using the unconditional power       *;
*          function, as illustrated in Example 5.9.                  *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
*  OR = ODDS RATIO = P12/P21
*  P12 AND P21 ARE THE DISCORDANT PROBABILITIES SATISFYING OR
        P12 = P21 * OR
*  Alpha = TYPE I ERROR LEVEL ONE MOR TWO-SIDED, E.G. INPUT 0.025 FOR A
     TWO-SIDED 0.05 TEST
*  Beta = TYPE II ERROR LEVEL, 1 - B = POWER
*
*  THE PROGRAM COMPUTES
*
*  N = TOTAL SAMPLE SIZE, USING the multinomial EXPRESSION
;
 
%macro doit;
DATA TWO; SET ONE;
P12 = P21 * OR;
D=P12 - P21;
Pd = P12 + P21;
ZA = PROBIT(1-Alpha);
ZB = PROBIT(1-Beta);
PBAR = Pd/2;
 
***  NULL VARIANCE;
SIG0U = SQRT(PD);
NUM0 = ZA*SIG0U;
 
***  MULTINOMIAL UNCONDITIONAL VARIANCE AND COMPUTATION;
SIG1U = PD - (d**2);
SIG1U = SQRT(SIG1U);
NUM1 = ZB*SIG1U;
NU =  (NUM0 + NUM1) / D;
NU = NU**2;
proc print;
run;
%mend;
 
data one; input OR p21 alpha beta;
cards;
2.0 0.125 0.025 0.10
1.8 0.25 0.025 0.10
;
 
%doit;
 
run;
