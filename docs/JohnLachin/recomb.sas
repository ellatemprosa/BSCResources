*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter6/recomb.sas                   *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: calls the macro %Newton to obtain an iterative estimate   *;
*          of the recombination fraction for Fisher's maize data     *;
*          presented in Example 6.6. Newton-Raphson and Fisher       *;
*          scoring are illustrated, each starting from the moment    *;
*          estimate starting value and also from the null value.     *;
*          This generates the computations described in Tables 6.1   *;
*          and 6.2.                                                  *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
*** The macro can be accessed either by "submitting" a job to sas before submitting this
    job, or by statements such as the following (remove the * first);
 
* filename newton 'c:\sasjobs\bioweb\chapter6\newton.sas';
* %include newton;
 
*** Define the homogeneous equation to be solved, the solution is at yofx=0;
 
%macro fofx(theta);
u=(&x1/(2+&theta)) - ((&x2+&x3)/(1-&theta)) + (&x4/&theta);
logl=(&x1*log((2+&theta)/4)) + ((&x2+&x3)*log((1-&theta)/4)) + (&x4*log(&theta/4));
yofx=u;
%mend fofx;
 
*** Newton-Raphson using hessian to determine the step size, dyofx;
%macro dfofx(theta);
du=(&x1/(2+&theta)**2) + ((&x2+&x3)/(1-&theta)**2) + (&x4/&theta**2);
dyofx=-du;
%mend dfofx;
 
data one;
x0=0.05704611;
%let x1=1997;
%let x2= 906;
%let x3= 904;
%let x4= 32;
%let error=0.000000001;
%let maxitn=15;
 
title1 'Fisher recombination fraction data';
title2 'Newton-Raphson iteration starting from moment estimator';
%newton;
 
proc print; var xlast ylast dylast xroot;
format xlast ylast dylast xroot 15.8;
run;
 
data one;
x0=0.25;
%let x1=1997;
%let x2= 906;
%let x3= 904;
%let x4= 32;
%let error=0.000000001;
%let maxitn=15;
 
title1 'Fisher recombination fraction data';
title2 'Newton-Raphson iteration starting from the null value';
%newton;
 
proc print; var xlast ylast dylast xroot;
format xlast ylast dylast xroot 15.8;
run;
 
*** Fisher scoring using negative information to determine the step size, dyofx;
 
%macro dfofx(theta);
N=&x1+&x2+&x3+&x4;
Inf=(N/4)*((1/(2+&theta)) + (2/(1-&theta)) + (1/&theta));
dyofx=-Inf;
%mend dfofx;
 
data one;
x0=0.05704611;
%let x1=1997;
%let x2= 906;
%let x3= 904;
%let x4= 32;
%let error=0.000000001;
%let maxitn=15;
 
title2 'Fisher scoring starting from moment estimator';
%newton;
 
proc print; var xlast ylast dylast xroot logl;
format xlast ylast dylast xroot logl 15.8;
run;
 
data one;
x0=0.25;
%let x1=1997;
%let x2= 906;
%let x3= 904;
%let x4= 32;
%let error=0.000000001;
%let maxitn=15;
 
title2 'Fisher scoring starting from the null value';
%newton;
 
proc print; var xlast ylast dylast xroot logl;
format xlast ylast dylast xroot logl 15.8;
run;
