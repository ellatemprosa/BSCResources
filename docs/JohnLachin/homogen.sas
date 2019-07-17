*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter4/homogen.sas                  *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: A PROC IML macro that computes the MVLE and a Wald test   *;
*          based on the MVLE and its estimated variance for a        *;
*          parameter estimate within independent strata, and also    *;
*          computes the contrast Wald test of the hypothesis of      *;
*          homogeneity. Here the parameter estimates could be any    *;
*          quantity. The estimate and is variance within each        *;
*          stratum are input using "cards".                          *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
%macro assoc;
data thetas; set tests (keep = Theta1-Theta&Levels);
data SEs; set tests (keep = SE1- SE&Levels);
 
proc iml;
reset print;  *** for test only;
use thetas;
read all into thetas;
use SEs;
read all into SEs;
var = diag( SEs##2);
j = j(1, &Levels, 1);
j1 = j(&Levels-1, 1, -1);
j2 = i(&Levels-1);
c =  j1 ¦¦ j2 ;
ivar = inv(var);
 
*** association test;
WT = j*ivar / (j*ivar*j`);
Theta = WT*thetas`;
SETheta = Sqrt(WT*var*WT`);
Z = Theta / SETheta;
P = 2 * (1 - probnorm(abs(Z)));
 
*** homogeneity test;
homog = (thetas*c`)* inv(c*var*c`) * (c*thetas`);
ph = 1 - Probchi(homog, &Levels-1);
 
out = Theta ¦¦ SETheta ¦¦ Z ¦¦ P ¦¦ ph;
create dataout from out;
Append from out;
 
quit;
 
data tests; if _N_ = 1 then set dataout; set tests;
rename col1 = Theta
       col2 = SETheta
       col3 = Z
       col4 = P
       col5 = ph;
 
proc print;
%mend assoc;
 
data one; input i thetak SEk;
cards;
1   0.0650  0.1741
2  -0.2134  0.2829
3  -0.3499  0.1742
4  -0.1695  0.2471
;
 
 
%let Levels = 4;
 
data two; set one end=eof;
drop i thetak SEk;
array Thetas (&Levels) Theta1-Theta&Levels;
array SEs    (&Levels) SE1- SE&Levels;
retain Theta1-Theta&Levels SE1- SE&Levels;
Thetas(i) = thetak;
SEs(i) = SEk;
if eof then output;
 
data tests; set two;
 
 %assoc;
 
run;
