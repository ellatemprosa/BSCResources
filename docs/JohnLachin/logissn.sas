*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter7/logissn.sas                  *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: computes sample size for logistic regression based on     *;
*          the power of a Wald test for a given design matrix,       *;
*          sample fractions and set of coefficients. This is used    *;
*          for the computations described in Example 7.9. In         *;
*          addition to the 1 df test described in the Example, the   *;
*          determination of sample size for the model Wald test is   *;
*          also illustrated.                                         *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
*** uses coefficients for main-effects models in LogistBetas.sas;
 
proc iml;
reset print ;
 
print 'Sample Size Computations in Example 7.9';
 
*** -- This example uses the same design matrix, sample fractions and odds ratios
       described in my article on sample size presented in the Encyclopedia of Biostatistics.
       Note that the coefficients used here differ from those in the biostat encyclopedia.
       Those are wrong because i took the main effects from a saturated model with
       interactions, ignoring the interaction term.;
 
   x ={ 1 0 0 0,
        1 0 0 1,
        1 1 0 0,
        1 1 0 1,
        1 0 1 0,
        1 0 1 1};
 
*** treatment or exposure group sample fractions within each stratum;
    rho = {0.5 , 0.5};
*** strata sample fractions;
    qs = {0.15 , 0.5 , 0.35 };
    q = qs@rho;
    q = diag(q);
***b ={ 1.09861, -0.25131, -0.47957 , 0.99524 };  *** old wrong coefficients;
   b ={ 1.2412, -0.3158, -0.7680 , 0.9390 };
   logit =-x*b;
   p ={ 1}/ ({1} + exp(logit));
   sinv = diag(p#({1}-p));
   sinv = q*sinv;
 
print '3 df model test';
print 'Using exemplary data coefficients';
 
   l ={ 0 1 0 0,
        0 0 1 0,
        0 0 0 1};
   crit = cinv({0.95},{ 3});        ***** specify 1-alpha and df;
   nc= cnonct(crit,{ 3},{ 0.1});    ***** specify df and beta;
   vbeta=inv(x`*sinv*x);
   delta=(l*b)` * inv(l* vbeta *l`) * (l*b);
   n = nc/ delta;
 
print 'same using wls coefficients, one step';
   b ={ -1.091516,  -0.176635,  -0.52179,  0.9214117};
   logit =-x*b;
   p ={ 1}/ ({1} + exp(logit));
   sinv = diag(p#({1}-p));
   sinv = q*sinv;
   delta=(l*b)` * inv(l* vbeta *l`) * (l*b);
   crit = cinv({0.95},{ 1});        ***** specify 1-alpha and df;
   nc= cnonct(crit,{ 1},{ 0.1});    ***** specify df and beta;
   n = nc/ delta;
 
 
print '1 df test of treatment effect';
print 'Using exemplary data coefficients';
 
   b ={ 1.2412, -0.3158, -0.7680 , 0.9390 };
   logit =-x*b;
   p ={ 1}/ ({1} + exp(logit));
   sinv = diag(p#({1}-p));
   sinv = q*sinv;
   l ={ 0 0 0 1 };
   crit = cinv({0.95},{ 1});        ***** specify 1-alpha and df;
   nc= cnonct(crit,{ 1},{ 0.1});    ***** specify df and beta;
   vbeta=inv(x`*sinv*x);
   delta=(l*b)` * inv(l* vbeta *l`) * (l*b);
   n = nc/ delta;
 
print 'same using wls coefficients, one step';
   b ={1.2249787, -0.299133,  -0.754154, 0.9262417 };
   logit =-x*b;
   p ={ 1}/ ({1} + exp(logit));
   sinv = diag(p#({1}-p));
   sinv = q*sinv;
   delta=(l*b)` * inv(l* vbeta *l`) * (l*b);
   n = nc/ delta;
 
run;
quit;
 
run;
