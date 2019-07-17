*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter9/expwlogl.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: This program determines the level of power provided by a  *;
*          specified sample size for a test of the equality of the   *;
*          survival distributions for two groups under an exponential*;
*          model assuming uniform entry over the period (0,R) and    *;
*          continued follow-up up to T>R years in all subjects. This *;
*          program uses the log hazard ratio as the basis for the    *;
*          test, rather than difference in hazards. The power        *;
*          function is presented in equation (2.4) of Lachin and     *;
*          Foulkes (1986). The program also provides for losses to   *;
*          follow-up that are also exponentially distributed, using  *;
*          (4.1) of Lachin and Foulkes.                              *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
options ls=132;
**;
**    to run this program, the steps are:
**       1.  specify the values of the macro variables;
**       2.  create a data set called indata;
**       3.  invoke the macro %power;
**;
**    macro variables:;
**       title = a text description of the design
**       alpha = two-sided type i error probability;
**       qe = sample fractions for the experimental group,;
**            assumed the same for all strata;
**       qc = sample fractions for the control group,;
**            assumed the same for all strata, q**= 1.0 - &qe;
**       n = total n;
**;
**    data set variables, one record per study or scenario;
**       r = period of recruitment;
**       t = total period of study;
**       lc = hazard rate (lc) for the control group;
**       rr = relative hazard in treaed group, le = rr*lc;
**       ldc = hazard rate for losses to followup in the control group;
**       lde = hazard rate for losses to followup in the experimental group;
**;
**;
%macro  qx(l,rr,tt);
   qx = exp(-(&l*(&tt-&rr)));
%mend;
%macro  gm(l,z,tt);
   %qx(&l,&z,&tt);
   qxz = qx;
   %qx(&l,0.0,&tt);
   qx0 = qx;
   gm = (qxz - qx0) / &l;
%mend;
%macro  sgd(hh,l,z,tt);
   %gm(&hh,&z,&tt);
   sgd = 1 / ( ( (&z*&l) - (&l* gm) ) / (&hh*&z) );
%mend;
%macro  num(hh,l,rr,tt,q);
   %sgd(&hh,&l,&rr,&tt);
   num = sgd / &q;
%mend;
%macro  ev(z,tt,l,hh);
   %gm(&hh,&z,&tt);
   ev = (&l/(&z*&hh)) * (&z - gm);
%mend;
 
%macro power;
 
title 'Power for the ratio of exponential hazards with exponential losses';
 
data one; set indata;
      le = rr*lc;
      lh = le*&qe + lc*&qc;
      the = lh + lde;
      thc = lh + ldc;
      te = lde + le;
      tc = ldc + lc;
      %ev(r,t,le,te);
      pde = ev;
      %ev(r,t,lc,tc);
      pdc = ev;
      %ev(r,t,lde,te);
      ple = ev;
      %ev(r,t,ldc,tc);
      plc = ev;
 
data two; set one;
za=probit(1 - (&alpha/2));
      %num(the,lh,r,t,&qe);
      tnuma = num;
      %num(thc,lh,r,t,&qc);
      tnumb = num;
      snull = tnuma + tnumb;
      %num(te,le,r,t,&qe);
      tnuma = num;
      %num(tc,lc,r,t,&qc);
      tnumb = num;
      salt = tnuma + tnumb;
      nterm= sqrt(n)*abs(log(rr));
      zb =  ( nterm - (za * sqrt(snull)) ) / sqrt(salt) ;
      power=probnorm(zb);
      dc = pdc * n * &qc;
      de = pde * n * &qe;
      nlc = plc * n * &qc;
      nle = ple * n * &qe;
 
data; set two;
if _n_ = 1 then do;
alph = &alpha;
qfe = &qe; qfc = &qc;
 
file print;
put "&title";
put @5 'type 1 error (2-sided):'       @38 alph   10.5;
put @5 '     z-value (2-sided):'       @38 za     10.5;
put @5 'sample fractions in e and c'   @38 (qfe qfc) (10.5 10.5);
put /;
put
'   periods                   relative'
        '   --------------- e v e n t s --------------- '
        '  --------------- l o s s e s ---------------';
put
'recruit follow    n    power  hazard'
       '     hazards      probabilities    exp. number '
       '     hazards      probabilities    exp. number';
put  @37
       ' control treated control treated control treated'
       ' control treated control treated control treated';
end;
file print;
put (r t) (7.1 7.1) n 6.0 power 8.3 rr 8.4
    (lc le pdc pde dc de ldc lde plc ple nlc nle)
    (8.4 8.4 8.4 8.4 8.2 8.2 8.4 8.4 8.4 8.4 8.2 8.2);
%mend;
 
 run;
 
%let title = LNCS power in Example 9.10;
%let alpha = 0.10;  *** 2-sided, 0.05 one-sided;
%let qc=0.5; %let qe=0.5;
 
data indata; input n r t lc rr ldc lde;
cards;
125 4  6  0.30 0.5  0.0  0.0
125 4  6  0.30 0.6  0.0  0.0
125 4  6  0.30 0.666667  0.0  0.0
125 4  6  0.30 0.5  0.05  0.05
125 4  6  0.30 0.6  0.05  0.05
125 4  6  0.30 0.666667  0.05  0.05
;
 
%power;
 
run;
