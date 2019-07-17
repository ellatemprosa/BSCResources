*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter7/ulcrlogi.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: fits the interaction model to the ulcer clinical trial    *;
*          data presented in Example 7.11.                           *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
options nodate formdlim='-'; *** ls=65;
 
data one;
input k n1 d1 n2 d2;
cards;
 1  42  16  47  20
 2  12   9   9   4
 3  46  28  44  16
*****;
Title1 'Example 4.1: Ulcer Clinical Trial';
 
data two; set one;
keep i j k f z2 z3;
z2=0; if k=2 then z2=1;
z3=0; if k=3 then z3=1;
d3=n1-d1;
d4=n2-d2;
i = 1; j = 1; f =d1; output;
i = 0; j = 1; f =d2; output;
i = 1; j = 0; f =d3; output;
i = 0; j = 0; f =d4; output;
 
data three; set two;
z2=0; if k=2 then z2=1;
z3=0; if k=3 then z3=1;
Iz2=I*z2;
iz3=I*z3;
n=1;
 
proc logistic descending;
  model j = i z2 z3 / rl;
weight f;
title2 'main effects model fit in Proc Logistic';
 
proc logistic descending;
  model j = i z2 z3 iz2 iz3 / rl;
weight f;
test iz2=iz3=0;
OUTPUT out=next pred=pred xbeta=xbeta;
title2 'interaction model fit in Proc Logistic';
 
proc sort; by i z2 z3;
data next; set next; by i z2 z3;
if first.z3;
odds = exp(xbeta);
 
proc print;
 
proc genmod data=three;
class i k;
model j/n = i k i*k
 / dist=binomial link=logit type3;
weight f;
title2 'interaction model fit in Proc Genmod';
 
run;
