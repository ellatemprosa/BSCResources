*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter9/swl-gof.sas                  *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: assesses the PH model assumptions using tests of          *;
*          interaction with time, plots of the log(-log(survival))   *;
*          function within strata, and the Lin (1991) test of the    *;
*          PH assumption using the SAS macro %gofcox. This performs  *;
*          the computations described in Example 9.8 and generates   *;
*          the functions plotted in Figure 9.3.                      *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
filename gofcox '/jml/biostatmethods/chapter9/gofcox.sas';
 
%include gofcox;
 
*** The data file is containted in the datasets directory;
 
libname bioweb '/jml/biostatmethods/datasets';
 
Title1 'Lagakos Squamous Cell Carcinoma Data';
 
*** Define covariates in the model. Additional covariates are defined in the procs below;
 
%let ivars= age perfstat group;
 
data zero; input &ivars;
cards;
 0     0     0
;
*** data set with zero values for all covariates;
 
data one; set bioweb.lagakos;
fail = cause; if cause=2 then fail=1;
delta=fail;
 
group=1-treatmnt;   *** 1=A, 0=B;
 
proc phreg data=one outest=oneest;
  model time*delta(0) = &ivars / corrb covb rl;
  output out=two survival=s xbeta=xb stdxbeta=sexb;
  baseline out=three covariates=zero
           survival=s logsurv=lns loglogs=lnlns xbeta=xb stdxbeta=sexb;
title2 'treatment adjusted for age and performance status';
 
data two; set two;
   s0 = s**(exp(-xb));
title3 'risk scores and survival estimates';
proc sort; by xb;
***proc print data=two;
 
data threea threeb; set three;
if xb=0 then output threea;
else output threeb;
 
data threeb; set threeb;
proc print; var time &ivars xb sexb s lns lnlns;
proc plot;
  plot s*time;
  plot lns*time;
  plot lnlns*time;
title4 'plots of background functions for x = means';
 
data threea; set threea;
proc print; var time &ivars xb sexb s lns lnlns;
proc plot;
  plot s*time;
  plot lns*time;
  plot lnlns*time;
title4 'plots of background functions for x = zero';
 
proc phreg data=one;
  model time*delta(0) = &ivars aget;
  aget=age*log(time);                    *** define interaction covariate;
title2 'model with Cox test of non-proportionality for age';
 
proc phreg data=one;
  model time*delta(0) = &ivars perft;
  perft=perfstat*log(time);
title2 'model with Cox test of non-proportionality for perfstat';
 
proc phreg data=one;
  model time*delta(0) = &ivars groupt;
  groupt=group*log(time);
title2 'model with Cox test of non-proportionality for group';
 
proc rank data=one groups=2 out=onegrp; var &ivars;
     ranks agegrp perfgrp groupgrp;
*** divide distribution of age into two groups above and below median;
 
proc phreg data=onegrp;
  model time*delta(0) = perfstat group;
  strata agegrp;
  output out=two survival=s xbeta=xb stdxbeta=sexb;
title2 'PH model stratified by age group';
 
data two; set two;
   lnt=log(time);
   s0 = s**(exp(-xb));
   if agegrp=0 then s0l = s0;
   if agegrp=1 then s0h = s0;
   lnlns0l = log(-log(s0l));
   lnlns0h = log(-log(s0h));
 
proc print; var time agegrp perfstat group xb s0 s0l s0h lnlns0l lnlns0h;
title3 'risk scores and survival estimates';
proc plot;
  plot lnlns0l*lnt='l'
       lnlns0h*lnt='h'
       / overlay;
title3 'plots of log(-log(s0)) within age strata';
 
run;
 
%gofcox (data=one, time=time, EvtFlag=delta, vars=&ivars);
 
run;
