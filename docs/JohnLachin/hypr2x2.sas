*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter2/hypr2x2.sas                  *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: computes hypergeometric probabilites for all possible     *;
*          2x2 tables with fixed margins, as a function of a         *;
*          specified odds ratio.                                     *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
%macro tablep;
  cp = probhypr(n, n1, m1, a, or);
  icp= 1 - cp;
  if a=al then p = cp;
  if a>al then p = cp - probhypr(n, n1, m1, a-1, or);
*** if p=. then p=cp; *** the first possible table;
%mend;
 
data one;
*** all possible 2 x 2 tables with fixed margins;
 
*** specify the margins of the 2x2 table;
 
n1=19; n2=10; m1=15; m2=14;
al=max(0, m1-n2); au=min(m1,n1);
n = n1 + n2;
 
do i = al to au by 1;
 b = m1-i;
 c = n1-i;
 d =m2 -b;
 a=i;
 
*** specify the odds ratio;
 or = 1;  **** null hypothesis;
** or = 0.0127; **** lower 95% ci for this example;
** or = 1.0952; **** upper 95% ci for this example;
%tablep;
 output;
end;
 
proc print; var N n1 m1 a b p cp icp;
title 'all possible tables with probabilities (p), cum. p (cp), and inverse cum. p (icp)';
 
run;
