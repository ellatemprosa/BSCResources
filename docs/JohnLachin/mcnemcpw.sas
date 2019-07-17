*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter5/mcnemcpw.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: Power for McNemar's test for paired or matched 2x2        *;
*          tables using the conditional power function, as           *;
*          illustrated in Example 5.10.                              *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
*
*  M = # discordant pairs
*  OR = DISEASE ODDS RATIO = P12/P21
*  Alpha = TYPE I ERROR LEVEL ONE MOR TWO-SIDED, E.G. INPUT 0.025 FOR A
     TWO-SIDED 0.05 TEST
*
*  THE PROGRAM COMPUTES
*     Power USING the conditional power EXPRESSION
;
 
%macro doit;
DATA TWO; SET ONE;
ZA = PROBIT(1-Alpha);
 
SIG0C = OR+1;
SIG1C = 2*SQRT(OR);
zb = ( (sqrt(M)*abs(or-1)) - ZA*sig0c ) / sig1c;
power = PROBnorm(zb);
proc print;
run;
%mend;
 
data one; input M OR alpha;
cards;
80 2.0 0.025
60 2.0 0.025
;
 
%doit;
 
run;
