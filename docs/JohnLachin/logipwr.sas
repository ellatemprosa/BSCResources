*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter7/logipwr.sas                  *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: computes the power of a Wald test in logistic regression  *;
*          for a given sample size as described in Example 7.9.      *;
*          note this version requires specification of desired       *;
*          non-centrality  parameter values for specified df and     *
*          power                                                     *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
proc iml;
 
reset print ;
 *** power of wald test for logistic regression;
 *** see sas SUGI supplemental users guide, chapter 29, p. 392;
 
print '*** Example 7.9';
 
   x ={ 1 0 0 0,  1 0 0 1,  1 1 0 0,  1 1 0 1,  1 0 1 0,  1 0 1 1};
   n ={200};
*** tx sample fractions;
    rho = {0.5 , 0.5};
*** strata sample fractions;
    qs = {0.15 , 0.5 , 0.35 };
    q = qs@rho;
    q = diag(q);
   b ={ 1.2412, -0.3158, -0.7680 , 0.9390 };
   logit =-x*b;
   p ={ 1}/ ({1} + exp(logit));
   sinv = diag(p#({1}-p));
   sinv = q*sinv;
   l ={ 0 0 0 1 };
   crit = cinv({0.95},{ 1});
   delta=(l*b)` * inv(l* inv(x`*sinv*x) *l`) * (l*b);
   nc = n * delta;
   power= 1 - probchi(crit,1,nc);
run;
