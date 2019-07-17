*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter8/rate.sas                     *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Authors: Shuping Lan                                               *;
*          John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: Two macros for the computation of rates and relative      *;
*          risks with an adjustment for over-dispersion. The macro   *;
*          %rate computes event rates and relative risks, with an    *;
*          adjustment for over-dispersion using the moment           *;
*          estimator for the dispersion variance. The macro          *;
*          %adjrate conducts a stratified analysis of relative       *;
*          risks over strata, computing the MVLE with the            *;
*          adjustment for overdispersion. The latter is not          *;
*          illustrated in the text. This program must be submitted   *;
*          before calling the macros. Instructions for use of the    *;
*          macros are provided in the this program.                  *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
options ls=80 nodate formdlim='-';
 
 
***********************************************************;
 
******* use of these macros as follows;
/*;
*****  copy this into your job and make the necessary changes;
 
*==============   initialze macro values  ================*;
 
* For macros %rate or %AdjRate;
   *-----------   variable title  -------------------*;
   %let tlt='dcct rates of severe hypoglycemia';
   *--------------   group          values  ------------*;
   %let grptxt1='standard';  **** coded as 0 *******;
   %let grptxt2='experimental';** coded as 1 *******;
   %let rrtxt='(ex/std)';    **** grp1/grp0  *******;
 
* For macro %AdjRate;
   *--------------   stratification values  ------------*;
   %let cattxt1='adolescents';**** coded as 0 *******;
   %let cattxt2='adults';     **** coded as 1 *******;
*
*  Specify the following in the macro call
*
*      indata= data set to be used;
*      evt= outcome event variable name fragment;
*           must have variables in the data set:
*            i&evt = indicator, 1=yes
*            n&evt = number events
*      futime= time of exposure in whatever units desired;
*      grp= group variable name, coded as 0,1;
*      out= output data set name;
*  For example, if the event name fragment is "event" then the data set
*      must contain the variables "ievent" and "nevent";
*
*  For an unstratified analysis, specify:
*
*    %rate(data,evt,grp,out);
*
*  For a stratified analysis, also specify
*      strata = the stratification variable, coded as 0,1;
*  Then specify;
*
*    %adjrate(data,evt,grp,strata,out);
*
*  This also includes the overall unstratified analysis;
*
* ****** end usage instructions;
*/;
******************************************************************;
 
 
%macro rate(indata,evt,grp,out);
  title2 'rate of event per patient-time unit of follow-up and relative risk (rr)';
 
data rr; set &indata;
rate=n&evt/futime;
 
PROC SORT; BY &grp; RUN;
proc univariate freq; var n&evt;
title3 'Distribution of number of events, overall and by group';
proc univariate freq; var n&evt; by &grp;
 
PROC MEANS N SUM MEAN uss css var; VAR FUtime;
title3 'Distribution of exposure time, overall and by group';
PROC MEANS N SUM MEAN uss css var; VAR FUtime; BY &grp;
 
PROC MEANS N SUM MEAN uss css var; VAR Rate; WEIGHT FUtime;
title3 'Mean rate, overall and by group';
PROC MEANS N SUM MEAN uss css var; VAR Rate; WEIGHT FUtime; BY &grp;
 
*-------- get estimates ------------------*;
data rr; set rr;
  proc sort; by &grp ;
 
   proc means n mean sum data=rr noprint;
      by &grp ;
      var n&evt;
      output out=bygrp1 sum=sevt mean=mevt
           stderr=vevt;
 
   proc means n mean sum data=rr noprint;
      by &grp ;
      var futime;
      output out=bygrp2 sum=sfutime mean=mfutime;
 
data bygrp; merge bygrp1 bygrp2 ; by &grp ;
 
*---------------------------------------------------;
data mhat; set bygrp; by &grp ;
  mhat=(sevt/sfutime);
  keep &grp sevt sfutime mhat;
 
proc print; var &grp sevt sfutime mhat;
title3 'number of events, total follow-up time and mean rate within each group';
 
data rr; merge rr mhat; by &grp ;
  x=n&evt;
  t=futime;
  tsq=t**2;
  Ex=mhat*t;
  a=((x-Ex)**2/t-mhat)/t;  *** Terms for dispersion variance as in Problem 8.3;
  sres = ((x - Ex)**2)/Ex;  *** for variance test;
 
  ressq=(rate-mhat)**2;
 
  proc means n mean sum data=rr noprint;
    by  &grp ;
    var sres;
    output out=outres sum=homogx2 n= nh;
 
data; set outres;
p=1 - probchi(homogx2, nh-1);
proc print; var &grp homogx2 nh p;
title3 'Cochran variance test of homogeneity within each group';
 
*----- use moment estimator of over-dispersion variance --------------*;
 
  proc means n mean sum data=rr noprint;
    by  &grp ;
    weight t;
    var a ressq;
    output out=out1 mean=var_i mr sum=sa wssres;
*--------------------------------------------------*;
 
  proc means n mean sum data=rr noprint;
    by &grp ;
    var t tsq;
    output out=out2 sum=sumt sumtsq;
 
data out; merge out1 out2 mhat; by &grp ;
 
  if var_i < 0 then var_i=0;  /* set to 0 */;
 
    var_m=(var_i*sumtsq+mhat*sumt)/(sumt**2);
    logm=log(mhat);
    v_logm = var_m/(mhat**2);
    se_logm=sqrt(v_logm);
    logm_l=logm-1.96*se_logm;
    logm_u=logm+1.96*se_logm;
    m_l=exp(logm_l);
    m_u=exp(logm_u);
 
data out; set out;
  proc print;
  var &grp  wssres var_i var_m
      logm se_logm logm_l logm_u
      mhat m_l m_u
      _freq_ ;
title3 'variance components, overdispersed variance and 95% confidence limits';
 
 proc transpose data=out out=logm prefix=logm;
    var logm;
 
  proc transpose data=out out=vlogm prefix=vlogm;
    var v_logm;
data logm; set logm;
  log_rr=logm2-logm1;
data vlogm; set vlogm;
  vv=(vlogm2+vlogm1);
  se_ln_rr=sqrt(vlogm2+vlogm1);
 
data usedata;
  merge logm vlogm;
 
  rr=(exp(log_rr));           **** relative risk (e/s)*****;
  rr_l=(exp(log_rr-1.96*se_ln_rr));
  rr_u=(exp(log_rr+1.96*se_ln_rr));
 
  z_rr = log_rr/se_ln_rr;
  prob=2*(1-probnorm(abs(z_rr)));
 
proc print; var log_rr rr se_ln_rr rr_l rr_u z_rr prob;
title3 'relative risks, unadjusted';
run;
%mend;******************************************************;
 
**+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++;
 
*******************************************************************
*      statistics for event rate and relative risk               *
* macro's  for rr (group 1/group 2 )  ************************
*                                                                 *
*      programmer:  shuping lan                                   *
*      date: 11/93                                                *
*            11/96 modified by jl to print var components
*                  and to use futime rather than fuday
*
*  this version gives rr by 2 groups and by 2 strata              *
*  specify :      indata= (data to be used);
*                 evt= (outcome event);
*                 grp= (group);
*                 cat= (strata);
*                 out= (output data set);
******************************************************************;
 
 
 
%macro adjrate(indata,evt,grp,cat,out);********************;
  title2 'rate of event and relative risk (rr) stratified by &cat';
 
 
*-------- get stratifies estimates ------------------*;
data rr; set &indata;
  if &cat=. then delete;
  rate=n&evt/futime;
  proc sort; by &grp &cat;
 
   proc summary data=rr;
      by &grp &cat;
      var i&evt n&evt;
      output out=bycat1 sum=s1-s2 mean=m1-m2
           stderr=v1-v2;
 
   proc summary data=rr vardef=weight;
      by &grp &cat;
      weight futime;
      var rate;
      output out=bycat2 sum=s3 mean=m3
           stderr=v3;
   proc summary data=rr;
      by &grp &cat;
      var futime;
      output out=bycat3 sum=s4 mean=m4;
 
data bycat_s;merge bycat1 bycat2 bycat3 ; by &grp &cat;
 
*---------------------------------------------------;
data temp; set bycat_s; by &grp &cat;
  sumx=s2;
  sumy=s4;
  mhat=(sumx/sumy);
  keep &grp &cat mhat;
 
data rr; merge rr temp; by &grp &cat;
  x=n&evt;
  t=futime;
  tsq=t**2;
  Ex=mhat*t;
 
  a=((x-Ex)**2/t-mhat)/t;
 
*----- use moment estimator --------------*;
 
  proc summary data=rr;
    by  &grp &cat;
    weight t;
    var a;
    output out=out1 mean=var_i;
 
*--------------------------------------------------*;
 
  proc summary data=rr;
    by &grp &cat;
    var t tsq mhat;
    output out=out2 sum=sumt sumtsq summhat
                    mean=m_t m_tsq mhat;   *** mean mhat is mhat **;
 
data out; merge out1 out2; by &grp &cat;
 
  if var_i < 0 then var_i=0;  /* set to 0 */;
 
    var_m=(var_i*sumtsq+mhat*sumt)/(sumt**2);
    logm=log(mhat);       /* logm is theda in jl's notes */;
    v_logm = var_m/(mhat**2);
    se_logm=sqrt(v_logm);
    logm_l=logm-1.96*se_logm;
    logm_u=logm+1.96*se_logm;
    m_l=exp(logm_l);
    m_u=exp(logm_u);
* proc print;
data out; set out;
 keep &grp &cat var_i var_m mhat m_l m_u logm v_logm _freq_
      se_logm logm_l logm_u;
  proc print;
  var &grp &cat
      var_i var_m
      logm se_logm logm_l logm_u
      mhat m_l m_u
      _freq_ ;
title3 'variance components, variance and 95% C.I. within group and strata';
 
 proc transpose data=out out=logm prefix=logm;
    var logm;
 
  proc transpose data=out out=vlogm prefix=vlogm;
    var v_logm;
data logm; set logm;
  dd1=logm3-logm1;
  dd2=logm4-logm2;
data vlogm; set vlogm;
  vv1=(vlogm3+vlogm1);
  vv2=(vlogm4+vlogm2);
  se_dd1=sqrt(vlogm3+vlogm1);
  se_dd2=sqrt(vlogm4+vlogm2);
 
data usedata;
  merge logm vlogm out;
  n=_freq_;
 
  *---- get adjusted rates weight by inverse variance ---------*;
  ivv1=1/vv1;
  ivv2=1/vv2;
  ww1=(1/vv1)/(1/vv1+1/vv2);
  ww2=(1/vv2)/(1/vv1+1/vv2);
  dd_adj=ww1*dd1+ww2*dd2;
  vv_adj=(ww1**2)*vv1+(ww2**2)*vv2;
  se_adj=sqrt(vv_adj);
 
  rr1=(exp(dd1));           **** relative risk *****;
  rr2=(exp(dd2));
  rr1_l=(exp(dd1-1.96*se_dd1));
  rr1_u=(exp(dd1+1.96*se_dd1));
  rr2_l=(exp(dd2-1.96*se_dd2));
  rr2_u=(exp(dd2+1.96*se_dd2));
 
  rr_adj=(exp(dd_adj));
  rra_l=(exp(dd_adj-1.96*se_adj));
  rra_u=(exp(dd_adj+1.96*se_adj));
 
  prob1=2*(1-probnorm(abs(dd1/se_dd1)));
  prob2=2*(1-probnorm(abs(dd2/se_dd2)));
  pp_adj=2*(1-probnorm(abs(dd_adj/se_adj)));
 
  keep &grp &cat mhat m_l m_u rr1 rr1_l rr1_u rr2 rr2_l rr2_u n
       prob1 prob2 dd1 dd2 vv1 vv2
       rr_adj rra_u rra_l pp_adj;
*-------- test for interaction ------------------------*;
 
  proc iml;
 
  start testi;
  use usedata;
  read point 1 var{dd1 dd2} into dd;
  read point 1 var{vv1 vv2} into va;
  vv=diag(va);
  cc={ -1 1 };
  k=ncol(cc);
  df=k-1;
  chisq=dd*(cc`)*inv(cc*vv*(cc`))*cc*(dd`);
  finish testi;
  run testi;
 
  create tests;
  append;
  data tests; set tests;
  pp=1-probchi(chisq,df);
 
*-------- add test results -------------------------*;
data strdata; merge usedata tests;
  proc sort; by &grp &cat;
* proc print;
 
******* done with stratified estimates ***************;
 
*----------- within tx group comparisons  -------------*;
*----------- and combined difference (stratified) ---------*;
data within; merge logm vlogm;
*proc print;
data within; set within;
  ff1=logm3-logm4;
  ff2=logm1-logm2;
  vvf1=(vlogm3+vlogm4);
  vvf2=(vlogm1+vlogm2);
  se_ff1=sqrt(vlogm3+vlogm4);
  se_ff2=sqrt(vlogm1+vlogm2);
  probf1=2*(1-probnorm(abs(ff1/se_ff1)));
  probf2=2*(1-probnorm(abs(ff2/se_ff2)));
  keep probf1 probf2 ;
*proc print;
************** done with within tx group comparison *******;
 
 
*-------- get pooled estimates ----------------------*;
data rr; set &indata;
  if &cat=. then delete;
  rate=n&evt/futime;
  proc sort; by &grp;
 
   proc summary data=rr;
      by &grp;
      var i&evt n&evt;
      output out=pool1 sum=s1-s2 mean=m1-m2
           stderr=v1-v2;
 
   proc summary data=rr vardef=weight;
      by &grp;
      weight futime;
      var rate;
      output out=pool2 sum=s3 mean=m3
           stderr=v3;
   proc summary data=rr;
      by &grp;
      var futime;
      output out=pool3 sum=s4 mean=m4;
data bycat_p; merge pool1 pool2 pool3 ; by &grp;
 
data temp; set bycat_p;by &grp;
  n=_freq_;
  sumx=s2;
  sumy=s4;
  mhat=(sumx/sumy);
  keep &grp mhat;
data rr; merge rr temp; by &grp;
  x=n&evt;
  t=futime;
  tsq=t**2;
  Ex=mhat*t;
 
  a=((x-Ex)**2/t-mhat)/t;
 
 
*----- use moment estmator --------------*;
 
  proc summary data=rr;
    by  &grp;
    weight t;
    var a;
    output out=out1 mean=var_i;
 
*--------------------------------------------------*;
 
  proc summary data=rr;
    by &grp ;
    var t tsq mhat;
    output out=out2 sum=sumt sumtsq summhat
                    mean=m_t m_tsq mhat;   *** mean mhat is mhat **;
 
data out; merge out1 out2; by &grp;
 
  if var_i < 0 then var_i=0;  /* set to 0 */;
 
    var_m=(var_i*sumtsq+mhat*sumt)/(sumt**2);
    logm=log(mhat);       /* logm is theda in jl's notes */;
    v_logm = var_m/(mhat**2);
    se_logm=sqrt(v_logm);
    logm_l=logm-1.96*se_logm;
    logm_u=logm+1.96*se_logm;
    m_l=exp(logm_l);
    m_u=exp(logm_u);
* proc print;
data out; set out;
 keep &grp var_i var_m mhat m_l m_u logm v_logm _freq_
      se_logm logm_l logm_u;
  proc print;
  var &grp
      var_i var_m
      logm se_logm logm_l logm_u
      mhat m_l m_u
      _freq_ ;
title3 'variance components, variance and 95% C.I. within group';
 
 proc transpose data=out out=logm prefix=logm;
    var logm;
 
  proc transpose data=out out=vlogm prefix=vlogm;
    var v_logm;
data logm; set logm;
  dd1=logm2-logm1;
data vlogm; set vlogm;
  vv1=(vlogm2+vlogm1);
  se_dd1=sqrt(vlogm2+vlogm1);
 
data usedata;
  merge logm vlogm out;
 
  n=_freq_;
 
  rr1=(exp(dd1));         ****** relative risk (e/s) *****;
  rr1_l=(exp(dd1-1.96*se_dd1));
  rr1_u=(exp(dd1+1.96*se_dd1));
 
  prob1=2*(1-probnorm(abs(dd1/se_dd1)));
 
  keep &grp mhat m_l m_u rr1 rr1_l rr1_u n
       prob1 dd1 vv1;
 
data poodata; set usedata; by &grp;
* proc print;
************ done with pooled data ******************;
 
 
*-------- get the original numbers------------------*;
data outs; merge strdata
                 bycat_s(keep=&grp &cat m1 s2 s4);
           by &grp &cat;
data outp; merge poodata
                 bycat_p(keep=&grp m1 s2 s4 );
           by &grp;
data outw; set within;
 
*-------- set them together ------------------------*;
data &out; set outs outp outw;
 
  cat="&cat";
  evt="&evt";
  pevt=m1; nevt=s2; fu=s4;
  array mmm   mhat m_l m_u;
    do over mmm;
    mmm=mmm * 100;   ******** rate per 100 pt yrs *******;
    end;
 
  keep cat evt &cat
       &grp n pevt nevt fu
       mhat m_l m_u
       rr1 rr1_l rr1_u
       rr2 rr2_l rr2_u prob1 prob2
       rr_adj rra_u rra_l pp_adj
       chisq df pp
       probf1 probf2;
 
data use;
   length cat $6 evt $4;
   set &out;
 
 
data rruse;
   length cat $6 evt $4;
   set &out;
   keep cat evt rr1 rr2 rr1_u rr1_l rr2_u rr2_l probf1 probf2
            rr_adj rra_l rra_u prob1 prob2 pp_adj chisq  pp;
proc sort;
   by evt;
data rruse2;
   set rruse (keep=evt rr2 rr2_l rr2_u prob2 pp rr_adj rra_l rra_u
                   pp_adj chisq);
   by evt;
   if rr2=. then delete;
data rruse3;
   set rruse (keep=evt probf1 probf2);
   by evt;
   if probf1=. then delete;
proc sql;
   create table work.join1 as
   select rruse.evt, rr1, rr1_u, rr1_l, prob1,
          rruse2.pp, rruse2.rr2, rruse2.rr2_l, rruse2.rr2_u,
          rruse2.prob2, rruse2.rr_adj, rruse2.rra_l, rruse2.rra_u,
          rruse2.pp_adj, rruse2.chisq
   from work.rruse, work.rruse2
   where rruse.evt=rruse2.evt;
proc sql;
   create table work.rrboth as
   select join1.evt, rr1, rr1_u, rr1_l, prob1,
          pp, rr2, rr2_l, rr2_u,
          prob2, rr_adj, rra_l, rra_u,
          pp_adj, chisq, probf1, probf2
   from work.join1, work.rruse3
   where join1.evt=rruse3.evt;
data &grp;
   set use (keep=evt &cat n pevt nevt mhat &grp);
   if &grp=1  ;
   rename n=n11
          pevt = p11
          nevt = e11
          mhat = r11;
data std;
   set use (keep=evt &cat n  pevt nevt mhat &grp);
   if &grp=0  ;
   rename n=n21
          pevt = p21
          nevt = e21
          mhat = r21;
data cols;
   merge &grp std;
   prtvar=_n_;
proc sort;
   by evt;
 
data prttab;
   merge cols rrboth;
   by evt;
   if mod(prtvar,3) = 0 and rr1=. then delete;
proc sort;
   by prtvar;
 
*****************************************************************;
data _null_;
   set prttab;
   p11=p11*100;
   p21=p21*100;
   file print header=h notitles;
 if prtvar=1 then
   put @2 &tlt/
       @5 &cattxt1   @18 n11 4. @23 p11 5.2 @31 e11 4. @45 r11 5.3
       @5            @61 n21 4. @66 p21 5.2 @74 e21 4. @88 r21 5.3
       @98 rr1 5.3 @108 '(' rr1_l 5.3 ', ' rr1_u 5.3 ')' @127 prob1 6.4 ;
   if prtvar=2 then put
       @5 &cattxt2  @18 n11 4. @23 p11 5.2 @31 e11 4. @45 r11 5.3
       @5            @61 n21 4. @66 p21 5.2 @74 e21 4. @88 r21 5.3
       @98 rr2 5.3 @108 '(' rr2_l 5.3 ', ' rr2_u 5.3 ')' @127 prob2 6.4/
       @85 'difference in rel. risk, p=' pp 6.4/;
   if prtvar=3 then put
       @5 'overall'  @18 n11 4. @23 p11 5.2 @31 e11 4. @45 r11 5.3
       @5            @61 n21 4. @66 p21 5.2 @74 e21 4. @88 r21 5.3
       @98 rr1 5.3 @108 '(' rr1_l 5.3 ', ' rr1_u 5.3 ')' @127 prob1 6.4/
       @5 'stratified-adjusted'
       @98 rr_adj 5.3 @108 '(' rra_l 5.3 ', ' rra_u 5.3 ')' @127 pp_adj
           6.4/
       @5 'within tx group' @40 'p = 'probf1 6.4 @83 'p = 'probf2 6.4//;
 
return;
h:
put //@30 &grptxt2       @73 &grptxt1   //
     @23 '% with' @31 'number of' @45 'rate/'
     @66 '% with' @74 'number of' @88 'rate/' @98 'relative risk'/
     @18 'n' @23 'event'  @32 'events'    @45 '100 py'
     @61 'n' @66 'event'  @75 'events'    @88 '100 py'
     @98 &rrtxt     @111 '95% c.i.' @127 'pvalue'/;
return;
run;
%mend;******************************************************;
 
**+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++;
 
 
run;
