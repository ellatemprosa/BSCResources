*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter7/logibeta.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: For a given design matrix and sampling fractions, this    *;
*          routine computes the expected logistic model parameters   *;
*          (betas) using two approaches. The first is to use the     *;
*          iterative maximum likelihood solution using a large set   *;
*          of hypothetical observations. The other is the wls or     *;
*          gsk method described by Rochon (1989) based on a          *;
*          one-step wls computation. This is used for Example 7.9    *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
options formdlim='-';
 
%macro getbetas;
 
data two; set one;
p2 = (or*p1) / ( (1-p1) + (or*p1));
 
data three; set two (keep=i i2 i3 li ri p1 p2);
retain tn;
if _n_ = 1 then tn = 10000;
j = 0; k = 0; f = int(tn * li * ri     * (1-p1)); output;
j = 0; k = 1; f = int(tn * li * ri     *    p1 ); output;
j = 1; k = 0; f = int(tn * li * (1-ri) * (1-p2)); output;
j = 1; k = 1; f = int(tn * li * (1-ri) *    p2 ); output;
 
proc logistic descending; model k = i2 i3 j ;
     weight f;
title2 'Betas from a logistic model for marginal treatment with stratum effects';
 
*** This section computes the wls or gsk logit paramter estimates
    as in Rochon (1989);
 
%let ncell=6;  ************ specify # cells;
data ps; set two (keep = p1 p2 li ri) end=eof;
title2 'Betas from a wls solution as in Rochon (1989)';
 
keep ps1-ps&ncell lp1-lp&ncell qs1-qs&ncell;
retain qs1-qs&ncell 0;         **** array of sample fractions;
retain ps1-ps&ncell m 0;         **** array of ps;
retain lp1-lp&ncell 0;             **** array of logits;
array qs (&ncell) qs1-qs&ncell;
array ps (&ncell) ps1-ps&ncell;
array lps (&ncell) lp1-lp&ncell;
m=m+1;
qs(m)=li*ri;
ps(m)=p1;
lps(m)=log(p1)-log(1-p1);
m=m+1;
qs(m)=li*(1-ri);
ps(m)=p2;
lps(m)=log(p2)-log(1-p2);
 
if eof then output;
 
data q; set ps (keep=qs1-qs&ncell);
data p; set ps (keep=ps1-ps&ncell);
data lp; set ps (keep=lp1-lp&ncell);
 
proc iml;
reset print;  *** for test only;
use q;
read all into q;
use p;
read all into p;
use lp;
read all into lp;
p=p`;
lp = lp`;
   x ={ 1 0 0 0,
        1 0 0 1,
        1 1 0 0,
        1 1 0 1,
        1 0 1 0,
        1 0 1 1};
 
*** cell sample fractions;
   q = diag(q);
   sinv = diag(p#({1}-p));
   sinv = q*sinv;
   vbeta=inv(x`*sinv*x);
   beta= vbeta * (x`*sinv*lp);
   obeta=0;
do while (abs(obeta-beta) > {0.00000001});
 obeta=beta;
 nu=-x*beta;                          *** reweighted parameters;
   p=1/({1}+exp(nu));
   sinv = diag(p#({1}-p));
   sinv = q*sinv;
   vbeta=inv(x`*sinv*x);
   beta= vbeta * (x`*sinv*lp);
end;
run;
 
quit;
 
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
 
data one; input i li ri p1 or;
i2=0; if i=2 then i2=1;
i3=0; if i=3 then i3=1;
 
cards;
1 0.15 0.5 0.75    4.0
2 0.50 0.5 0.70    3.2
3 0.35 0.5 0.65    1.8
;
title1 'Example 7.9';
 
%getbetas;
 
run;
