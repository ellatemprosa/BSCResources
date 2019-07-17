*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter5/mcnemupw.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: Power for McNemar's test for paired or matched 2x2        *;
*          tables using the unconditional power function. This       *;
*          calculation is not illustrated in the book.               *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
*
* UNCONDITIONAL CALCULATIONS FOR Power
* FOR TEST OF 2 MATCHED MOR PAIRED PROPORTIONS;
 
*  OR = ODDS RATIO = P12/P21
*  P12 AND P21 ARE THE DISCORDANT PROBABILITIES SATISFYING OR
        P12 = P21 * OR
*  Alpha = TYPE I ERROR LEVEL ONE MOR TWO-SIDED, E.G. INPUT 0.025 FOR A
     TWO-SIDED 0.05 TEST
*  N = TOTAL SAMPLE SIZE
*
*  THE PROGRAM COMPUTES the associated level of power
*
*  Beta = TYPE II ERROR LEVEL, 1 - B = POWER, USING the multinomial EXPRESSION
;
 
 
%macro doit;
DATA TWO; SET ONE;
P12 = P21 * OR;
D=P12 - P21;
Pd = P12 + P21;
ZA = PROBIT(1-Alpha);
PBAR = Pd/2;
 
***  NULL VARIANCE;
SIG0U = SQRT(PD);
NUM0 = ZA*SIG0U;
 
***  MULTINOMIAL UNCONDITIONAL VARIANCE AND COMPUTATION;
SIG1U = PD - (d**2);
SIG1U = SQRT(SIG1U);
zb = (sqrt(N)*abs(d) - num0)/sig1u;
power = probnorm(zb);
proc print;
run;
%mend;
 
data one; input OR p21 alpha N;
cards;
1.8 0.25 0.025 180
1.8 0.25 0.025 150
1.8 0.25 0.025 125
;
 
%doit;
 
run;
