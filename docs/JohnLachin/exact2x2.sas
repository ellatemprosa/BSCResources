*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter2/exact2x2.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: This program takes the exact confidence limits on the     *;
*          odds ratio from StatXact (or SAS) and then computes the   *;
*          corresponding limits on the probabilities, and the exact  *;
*          corresponding limits on the risk difference and the       *;
*          relative risk. For illustration, the exact computations   *;
*          in SAS PROC FREQ are also presented for the odds ratio,   *;
*          as well as the asymptotic limits for the other scales.    *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
%macro convert;
or=&orin;
term1=or*(n1+m1) + (n2-m1);
term2=sqrt((term1**2) - 4*or*m1*(n1*(or-1)));
term3=2*(n1*(or-1));
t=(term1-term2)/term3;
%mend;
 
data one; input a b c d orl oru;
n1=a+c;
n2=b+d;
m1=a+b;
m2=c+d;
N=n1+n2;
 
cards;
7 8 12 2 0.0127 1.0952
;
 
data two; set one;
%let orin=orl; %convert; t1l=t;
  t2l=(m1 - n1*t)/n2;
%let orin=oru; %convert; t1u=t;
  t2u=(m1 - n1*t)/n2;
rdl=t1l-t2l;
rdu=t1u-t2u;
rrl=t1l/t2l;
rru=t1u/t2u;
proc print;
 
run;
 
data two; set one;
i = 0; j = 0; f = d; output;
i = 0; j = 1; f = c; output;
i = 1; j = 0; f = b; output;
i = 1; j = 1; f = a; output;
 
proc freq; tables i*j / cmh riskdiff;
exact or;
weight f;
 
run;
