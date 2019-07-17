*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter4/cmhtstsn.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: computes the sample size for the Cochran-Mantel-Haenszel  *;
*          test for stratified 2x2 tables and was used for the       *;
*          computations in Example 4.25 and Table 4.11. The program  *;
*          also computes the sample size from the unconditional      *;
*          marginal parameters, but without an adjustment for bias.  *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
options formdlim='-';
 
%macro witwal;
data two; set one end=eof;
retain v1 v0 e gv ap1 ap2 ar1 ar2 0;
p2 = p1 * or / (1 - p1 + (p1*or));
p  = (r * p1) + ((1 - r)* p2);
d = p2 - p1;
e = e + (l * r * (1-r) * d);
v1 = v1 + (l * r * (1-r) * ( ((1-r)*p1*(1-p1)) + (r*p2*(1-p2)) ) );
v0 = v0 + (l * r * (1-r) * p * (1-p));
**** averages for pooled analysis;
  ap1 = ap1 + (l * r * p1);
  ap2 = ap2 + (l * (1-r) * p2);
  ar1  = ar1  + (l * r);
  ar2  = ar2  + (l * (1-r));
if eof then do;
   za = probit(0.975);
   zb = probit(0.90);
   s0 = sqrt(v0);
   s1 = sqrt(v1);
   n = ( ( (za*s0) + (zb*s1) ) / e ) **2;
end;
title2 'stratified sample size computation';
proc print;
 
proc sort; by descending i;
 
data; set two;
if _n_ = 1;
p1 = ap1/ar1;
p2 = ap2/ar2;
q1 = ar1; q2 = ar2;
delta=p2-p1;
pbar = q1*p1 + q2*p2;
sig0 = pbar*(1-pbar)*( (1/q1) + (1/q2) );
sig0 = sqrt(sig0);
sig1 = ( p1*(1-p1)/q1) + ( p2*(1-p2)/q2);
sig1 = sqrt(sig1);
n = ( (za*sig0) + (zb*sig1) ) / delta;
n = n**2;
n1 = q1 * n;
n2 = q2 * n;
rr = p2/p1;
or = p2*(1-p1)/(p1*(1-p2));
proc print; var p1 p2 delta rr or pbar za zb n q1 n1 q2 n2;
title2 'sample size for pooled (unstratified) analysis';
 
%mend;
 
****   i = stratum
       l = lambda(i) = n(i) / n
       r = rho(i) = n(i,1) / n(i)
                    n(i,2) = (1 - r)*n(i)
       p1 = tau(i,1) = probability of event in group 1
       or = odds ratio in the i-th stratum = p2(1-p1) / p1(1-p2)
            p2 = p1*(or) / Ý1 - p1 + (p1*or)¨
               = p1 / Ý(1/or)(1-p1) + p1¨
****;
 
data one; input i l r p1 or;
cards;
1 0.15 0.3  0.75    4.0
2 0.50 0.6  0.70    3.2
3 0.35 0.45 0.65    1.8
;
title1 'Example 4.25, Table 4.11: heterogeneous and Unbalanced';
 
%witwal;
 
data one; input i l r p1 or;
cards;
1 0.15 0.5 0.75    4.0
2 0.50 0.5 0.70    3.2
3 0.35 0.5 0.65    1.8
;
title1 'Example 4.25, Table 4.11: heterogeneous and balanced';
 
%witwal;
 
data one; input i l r p1 or;
cards;
1 0.15 0.3  0.75    2.5224
2 0.50 0.6  0.70    2.5224
3 0.35 0.45 0.65    2.5224
;
title1 'Example 4.25, Table 4.11: homogeneous and Unbalanced';
 
%witwal;
 
data one; input i l r p1 or;
cards;
1 0.15 0.5 0.75    2.5224
2 0.50 0.5 0.70    2.5224
3 0.35 0.5 0.65    2.5224
;
title1 'Example 4.25, Table 4.11: homogeneous and balanced';
 
%witwal;
 
run;
 
run;
