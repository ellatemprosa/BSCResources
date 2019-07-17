*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter8/rate0.sas                    *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: The macro %rate0 that computes the variance of the mean   *;
*          rate within each group under the null hypothesis with an  *;
*          adjustment for over-dispersion within each group, and the *;
*          z test of the differences between rates.                  *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
options ls=80 formdlim='-';
 
***%rate0(allevts,events,grp_e,outrr);
 
***********************************************************;
 
*==============   initialze macro values  ================*;
   *-----------   variable title  -------------------*;
   %let tlt='dcct rates of severe hypoglycemia';
   *--------------   group          values  ------------*;
   %let grptxt1='standard';  **** coded as 0 *******;
   %let grptxt2='experimental';** coded as 1 *******;
   %let rrtxt='(ex/std)';    **** grp1/grp0  *******;
   *--------------   stratification values  ------------*;
   %let cattxt1='adolescents';**** coded as 0 *******;
   %let cattxt2='adults';     **** coded as 1 *******;
 
*******************************************************************
*
*  Specify :      indata= (data to be used);
*                 evt= (outcome event);
*                     must have variables in the data set:
                         i&evt = indicator, 1=yes
                         n&evt = number events
                         futime= years of exposure
*                 grp= (group);
*                 out= (output data set);
******************************************************************;
 
 
 
%macro rate0(indata,evt,grp,out);********************;
  title2 'rate of event per patient-time unit of follow-up';
 
*-------- get estimates ------------------*;
data rr; set &indata;
  proc sort; by &grp ;
 
   proc means n mean sum data=rr noprint;
      by &grp ;
      var n&evt;
      output out=bygrp1 sum=sevt mean=mevt
           stderr=vevt;
 
   proc means n mean sum data=rr noprint;
      by &grp ;
      var futime;
      output out=bygrp2 sum=sftime mean=mftime;
 
data bygrp_s;merge bygrp1 bygrp2 ; by &grp ;
 
   proc means n mean sum data=rr noprint;
      var n&evt;
      output out=bygrp01 sum=sevt mean=mevt
           stderr=vevt;
 
   proc means n mean sum data=rr noprint;
      var futime;
      output out=bygrp02 sum=sftime mean=mftime;
 
data bygrp_0; merge bygrp01 bygrp02 ;
 
*---------------------------------------------------;
 
data mhat; set bygrp_s; by &grp ;
  mhat=(sevt/sftime);
  keep &grp mhat;
 
data mhat0; set bygrp_0;
  mhat0=(sevt/sftime);
  keep mhat0;
 
data rr; if _n_=1 then set mhat0;
 set rr; by &grp;
  retain mhat0;
  x=n&evt;
  t=futime;
  tsq=t**2;
  mt=mhat0*t;
 
  a=((x-mt)**2/t-mhat0)/t;  *** Problem 8.3;
 
*----- Moment equation --------------*;
 
  proc means n mean sum data=rr noprint;
    by &grp;
    weight t;
    var a;
    output out=out1 mean=var_i;
 
*--------------------------------------------------*;
 
  proc means n mean sum data=rr noprint;
    by &grp ;
    var t tsq;
    output out=out2 sum=sumt sumtsq;
 
data out; if _n_=1 then set mhat0;
 merge out1 out2 mhat; by &grp ;
 
 retain mhat0;
    var_m=(var_i*sumtsq+mhat0*sumt)/(sumt**2);
    logm=log(mhat);
    v_logm = var_m/(mhat0**2);
    se_logm=sqrt(v_logm);
    logm_l=logm-1.96*se_logm;
    logm_u=logm+1.96*se_logm;
    m_l=exp(logm_l);
    m_u=exp(logm_u);
* proc print;
data out; set out;
  proc print;
  var &grp  var_i var_m
      logm v_logm se_logm logm_l logm_u
      mhat m_l m_u mhat0
      _freq_ ;
title3 'variance components and variance of mu-hat under the null, and 95% C.I.';
 
 proc transpose data=out out=mhat prefix=mhat;
    var mhat;
 
  proc transpose data=out out=var_m prefix=var_m;
    var var_m;
 
data mhat; set mhat;
  mdif=mhat2-mhat1;
data var_m; set var_m;
  vv=(var_m2+var_m1);
  sedif=sqrt(vv);
 
data usedata;
  merge mhat var_m;
 
  z = mdif/sedif;
  prob=2*(1-probnorm(abs(z)));
proc print; var mdif sedif z prob;
title3 'difference in rates, SE under H0 with over-dispersion, and the Z-test';
 
run;
%mend;******************************************************;
 
run;
