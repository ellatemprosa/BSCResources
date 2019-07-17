*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter9/survanal.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: is a set of survival analysis macros that provides the    *;
*          basic computations of survival probabilities and hazard   *;
*          rates in each of two groups, either for data in           *;
*          continuous time or using the actuarial method, with       *;
*          generalized Mantel-Haenszel or G-rho tests for            *;
*          differences between groups. Instructions for use are      *;
*          given in the program.                                     *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
Macro TABLES(data=data set name, default is last):
 
     This macro reads an input data set constructed by the user which
     contains one observation per subject as follows:
         TIME: the time of the event or of censoring for that subject.
              These values must be less than 99999.
         DELTA: The indicator variable with values:
              1 if TIME is an event time for this subject
              0 if TIME is a censoring time for this subject
         GROUP: the treatement group number with values 1 or 2 for
              each subject.
 
     This macro then constructs a data set TABLE which contains counts
     of the numbers at risk, the numbers of events, and the numbers of
     censored observations in each group and all total for each distinct
     value of TIME.
 
Macro ACTUAR(width= interval width):
     This macro groups the time values into intervals of the desired
     width and then applies the actuarial adjustment to obtain the
     adjusted numbers at risk in each interval.
 
Macro TESTS: Arguments
 
     %let type='plimit' or  %let type='actuar'  *** CASE SENSITIVE ***
     plimit is the default:
 
     This macro computes and plots the survival curves and conducts the
     tests of significance.  The default option for the plot routine
     'plimit' calls for the product limit estimators to be plotted
     (a step function).  The 'actuar' option calls for only the interval
     associated values to be plotted.
 
     Confidence limits for the survival function are computed using the log, logit
     and log(-log) transformations.
 
     The hazard function is also computed.
 
     In addition to the family of weighted Mantel-Haenszel tests, pointwise tests
     of the difference between the survival functions is computed using the null
     variance of the group difference.
 
     The estimated difference between the survival probabilities, its
     the SE and confidence limits are computed.
 
     The estimated ratio of the cumulative hazards, or difference in
     log survival probabilities, its SE and confidence limits are computed.
 
     The estimated log survival odds ratio, its SE and confidence limits are computed.
 
     If the data is available as a table giving the numbers at risk,
     the numbers censored and the numbers of events, this macro can be
     called directly without having to first call the above macros. The input data
     set should contain the variables
         time = event time or end of interval for grouped time.
         d, d1, d2 = numbers of events total, in group1 and group 2 at that time.
         w, w1, w2 = numbers lost to follow-up in this interval (for actuarial
                     computations, not needed for Kaplan-Meier).
         n, n1, n2 = the numbers at risk for a Kaplan-Meier calculation, or the numbers
                     entering the interval for an actuarial computation.
**==================================================================**;
 
%macro tables(data=_last_);
** set up table input for calculation of product-limit life tables;
proc sort data=&data; by descending time;
data table; set &data; by descending time;
retain w w1 w2 d d1 d2 n n1 n2;
if _n_=1 then do;
   n = 0; n1 = 0; n2 = 0;
   end;
n = n+1;
if group=1 then n1=n1+1;
if group=2 then n2=n2+1;
if first.time then do;
   d = 0.0; d1 = 0.0; d2 = 0.0;
   w = 0.0; w1 = 0.0; w2 = 0.0;
   end;
if delta = 1 then do;
   d = d+1;
   if group=1 then d1=d1+1;
   if group=2 then d2=d2+1;
   end;
if delta = 0 then do;
   w = w+1;
   if group=1 then w1=w1+1;
   if group=2 then w2=w2+1;
   end;
if last.time then do;
   if time ne . then output;
   end;
proc sort; by time;
 
proc print; var time n d w n1 d1 w1 n2 d2 w2;
 title5
 'Table of basic data for each time value';
%mend;
 
 
**==================================================================**;
 
 
%macro actuar(width=);
** modify table input for calculation of actuarial life tables;
** over intervals of specified width;
data tempa; set table;
if time=0 then time=0.001;
interval = ceil(time/&width);
proc sort; by interval;
 
data tempb; set tempa; by interval;
keep time interval aw aw1 aw2 ad ad1 ad2 an an1 an2;
retain aw aw1 aw2 ad ad1 ad2 an an1 an2;
if first.interval then do;
   an = n; an1 = n1; an2 = n2;
   aw = 0; aw1 = 0; aw2 = 0;
   ad = 0; ad1 = 0; ad2 = 0;
   end;
aw = aw + w; aw1 = aw1 + w1; aw2 = aw2 + w2;
ad = ad + d; ad1 = ad1 + d1; ad2 = ad2 + d2;
if last.interval then do;
   time=interval*&width;
   output;
   end;
***proc print;
 
data table; set tempb;
keep interval time w w1 w2 d d1 d2 n n1 n2;
w = aw; d = ad; n = an;
w1 = aw1; d1 = ad1; n1 = an1;
w2 = aw2; d2 = ad2; n2 = an2;
n = n - (w/2);
n1 = n1 - (w1/2);
n2 = n2 - (w2/2);
proc print; var interval time w w1 w2 d d1 d2 n n1 n2;
 title5
'Table of basic data for each time value with actuarial adjustment';
%let type='actuar';
%mend;
 
 
**==================================================================**;
 
 
%macro tests;
** calculate survivor functions and rank tests;
 
data ztable ltimes; set table end=eof;
retain lastime1 lastime2 lastime;
if _n_ = 1 then do;
   lastime=0; lastime1=0; lastime2=0;
end;
elapsed=.; elapsed1=.; elapsed2=.;
if d>0 then do;
   elapsed = time - lastime;
   lastime=time;
end;
if d1>0 then do;
   elapsed1 = time - lastime1;
   lastime1=time;
end;
if d2>0 then do;
   elapsed2 = time - lastime2;
   lastime2=time;
end;
output ztable;
if eof then output ltimes;
 
data ztwo; if _n_=1 then set ltimes (keep = lastime1 lastime2 lastime);
   set ztable (keep = time n n1 n2 d d1 d2 elapsed elapsed1 elapsed2);
 
data two twol; set ztwo end=eof;
retain sp sp1 sp2 vlnsp vlnsp1 vlnsp2 vlnsp10 vlnsp20 z95 0;
retain td tpy td1 tpy1 td2 tpy2 0;
 
if _n_=1 then do;
   z95=probit(0.975);
   sp = 1.0; sp1 = 1.0; sp2 = 1.0;
   vlnsp = .; vlnsp1 = .; vlnsp2 = .;
   end;
if n > 0 then do;
   p = (n-d)/n;       q = d/n;
   if elapsed>0 then do;
      If "%lowcase(&type)"="plimit" then do;
        py = n*elapsed;
        td=max(td+d, td);
        tpy=max(tpy+py, tpy);
        if q > 0 then do;
          h=q/elapsed;               *** Nelson/Aalen hazard estimate;
          vh=(q*p/n)/(elapsed**2);
          seh=sqrt(vh);
        end;
      end;
      If "%lowcase(&type)"="actuar" then do;
        py=n*((1+p)*elapsed)/2;
        td=max(td+d, td);
        tpy=max(tpy+py, tpy);
        if q > 0 then do;
          h=2*q/((1+p)*elapsed);      *** actuarial hazard estimate;
          vh=((h**2)/d)*(1 - ((h*elapsed/2)**2));
          seh=sqrt(vh);
        end;
      end;
   end;
  if time le lastime then do;
   sp = sp*p;
  if 0<sp<1 then do;
   lnsp=log(sp);                       *** log (S) ;
   if vlnsp=. then vlnsp=0;
   vlnsp=vlnsp + (q/(n*p));              * variance of log survival;
   selnsp=sqrt(vlnsp);
   lnspu=lnsp+z95*selnsp;                * 95% confidence interval (CI);
   lnspl=lnsp-z95*selnsp;
   alnspu=exp(lnspu);                    * asymmetric CI limits on survival probability;
   alnspl=exp(lnspl);
   vsp = (sp**2) * vlnsp;              *** Greenwood variance estimate;
   sesp = sqrt(vsp);                     * and confidence interval;
   spu=sp+z95*sesp;
   spl=sp-z95*sesp;
   lnlnsp=log(-log(sp));               *** log(-log(S));
   vlnlnsp=(lnsp**(-2))*vlnsp;           * variance of log(-log(S));
   selnlnsp=sqrt(vlnlnsp);
   lnlnspu=lnlnsp+z95*selnlnsp;          * CI on log(-log(S);
   lnlnspl=lnlnsp-z95*selnlnsp;
   alnlnspu=exp(-exp(lnlnspl));          * asymmetric CI on S;
   alnlnspl=exp(-exp(lnlnspu));
   lno=log(sp/(1-sp));                 *** logit of S;
   vlno=((1-sp)**-2)*vlnsp;              * Variance of logit of S;
   if lno=. then vlno=.;
   selno=sqrt(vlno);
   lnou=lno+z95*selno;                   * CI on logit of S;
   lnol=lno-z95*selno;
   alnou=1/(1+exp(-lnou));               * asymmetric CI on S;
   alnol=1/(1+exp(-lnol));
   end;
  end;
  end;
  if time gt lastime then sp=.;
if n1 > 0 then do;                       *** Repeat above in Group 1;
   p1 = (n1-d1)/n1;   q1 = d1/n1;
   if elapsed1>0 then do;
      if "%lowcase(&type)" = 'plimit' then do;
        py1 = n1*elapsed1;
        td1=max(td1+d1, td1);
        tpy1=max(tpy1+py1, tpy1);
        if q1 > 0 then do;
          h1=q1/elapsed1;
          vh1=(q1*p1/n1)/(elapsed1**2);
          seh1=sqrt(vh1);
        end;
      end;
      if "%lowcase(&type)" = 'actuar' then do;
        py1=n1*((1+p1)*elapsed1)/2;
        td1=max(td1+d1, td1);
        tpy1=max(tpy1+py1, tpy1);
        if q1 > 0 then do;
          h1=2*q1/((1+p1)*elapsed1);
          vh1=((h1**2)/d1)*(1 - ((h1*elapsed1/2)**2));
          seh1=sqrt(vh1);
          end;
        end;
   end;
  if time le lastime1 then do;
   sp1 = sp1*p1;
  if 0<sp1 < 1 then do;
   lnsp1=log(sp1);                     *** log (S);
   if vlnsp1=. then vlnsp1=0;
   vlnsp1 = vlnsp1 + (q1/(n1*p1));       * variance;
   vlnsp10 = vlnsp10 + (q/(n1*p));       * variance estimated under H0;
   selnsp1=sqrt(vlnsp1);
   lnspu1=lnsp1+z95*selnsp1;
   lnspl1=lnsp1-z95*selnsp1;
   alnspu1=exp(lnspu1);
   alnspl1=exp(lnspl1);
   vsp1 = (sp1**2) * vlnsp1;           *** Greenwood estimate;
   sesp1 = sqrt(vsp1);
   vsp10 = (sp**2) * vlnsp10;            * estimate under H0;
   sesp10 = sqrt(vsp10);
   spu1=sp1+z95*sesp1;
   spl1=sp1-z95*sesp1;
   lnlnsp1=log(-log(sp1));             *** log(-log(S));
   vlnlnsp1=(lnsp1**(-2))*vlnsp1;
   selnlns1=sqrt(vlnlnsp1);
   lnlnspu1=lnlnsp1+z95*selnlns1;
   lnlnspl1=lnlnsp1-z95*selnlns1;
   allspu1=exp(-exp(lnlnspl1));
   allspl1=exp(-exp(lnlnspu1));
   lno1=log(sp1/(1-sp1));              *** logit of S;
   vlno1=((1-sp1)**-2)*vlnsp1;
   if lno1=. then vlno1=.;
   selno1=sqrt(vlno1);
   lnou1=lno1+z95*selno1;
   lnol1=lno1-z95*selno1;
   alnou1=1/(1+exp(-lnou1));
   alnol1=1/(1+exp(-lnol1));
   end;
  end;
  end;
  if time gt lastime1 then sp1=.;
if n2 > 0 then do;                    *** Repeat above in Group 2;
   p2 = (n2-d2)/n2;   q2 = d2/n2;
   if elapsed2>0 then do;
      if "%lowcase(&type)" = 'plimit' then do;
        py2 = n2*elapsed2;
        td2=max(td2+d2, td2);
        tpy2=max(tpy2+py2, tpy2);
        if q2 > 0 then do;
          h2=q2/elapsed2;
          vh2=(q2*p2/n2)/(elapsed2**2);
          seh2=sqrt(vh2);
        end;
      end;
      if "%lowcase(&type)" = 'actuar' then do;
        py2=n2*((1+p2)*elapsed2)/2;
        td2=max(td2+d2, td2);
        tpy2=max(tpy2+py2, tpy2);
        if q2 > 0 then do;
          h2=2*q2/((1+p2)*elapsed2);
          vh2=((h2**2)/d2)*(1 - ((h2*elapsed2/2)**2));
          seh2=sqrt(vh2);
        end;
      end;
   end;
  if time le lastime2 then do;
   sp2 = sp2*p2;
  if 0<sp2 < 1 then do;
   lnsp2=log(sp2);
   if vlnsp2=. then vlnsp2=0;
   vlnsp2 = vlnsp2 + (q2/(n2*p2));
   vlnsp20 = vlnsp20 + (q/(n2*p));
   selnsp2=sqrt(vlnsp2);
   lnspu2=lnsp2+z95*selnsp2;
   lnspl2=lnsp2-z95*selnsp2;
   alnspu2=exp(lnspu2);
   alnspl2=exp(lnspl2);
   vsp2 = (sp2**2) * vlnsp2;
   sesp2 = sqrt(vsp2);
   vsp20 = (sp**2) * vlnsp20;
   sesp20 = sqrt(vsp20);
   spu2=sp2+z95*sesp2;
   spl2=sp2-z95*sesp2;
   lnlnsp2=log(-log(sp2));
   vlnlnsp2=(lnsp2**(-2))*vlnsp2;
   selnlns2=sqrt(vlnlnsp2);
   lnlnspu2=lnlnsp2+z95*selnlns2;
   lnlnspl2=lnlnsp2-z95*selnlns2;
   allspu2=exp(-exp(lnlnspl2));
   allspl2=exp(-exp(lnlnspu2));
   lno2=log(sp2/(1-sp2));
   vlno2=((1-sp2)**-2)*vlnsp2;
   if lno2=. then vlno2=.;
   selno2=sqrt(vlno2);
   lnou2=lno2+z95*selno2;
   lnol2=lno2-z95*selno2;
   alnou2=1/(1+exp(-lnou2));
   alnol2=1/(1+exp(-lnol2));
   end;
  end;
  end;
  if time gt lastime2 then sp2=.;
output two;
 
if eof then output twol;
****
  data set twol contains the times of the last event in each group,
  the lastime values;
****;
 
data twol; set twol;
   th = td / tpy;
   th1 = td1 / tpy1;
   th2 = td2 / tpy2;
 
 title4 "Survival function estimates using the &type option";
 title5
'Proportion surviving each interval (p), ';
 title6
 'Lifetable estimate of survivors (sp), and s.e. (sesp), ';
 title7
 'Linearized hazard over interval (h) and s.e. (seh)';
proc print data=two;
 var time elapsed n d p sp sesp h py seh;
title8 'Life-table for combined groups';
proc print data=two;
 var time spl spu lnsp selnsp lnspl lnspu alnspl alnspu;
title9 'Pointwise 95% confidence limits using s and log(s)';
proc print data=two;
 var time lnlnsp selnlnsp lnlnspl lnlnspu alnlnspl alnlnspu;
title9 'Pointwise 95% confidence limits using log(log(s))';
proc print data=two;
 var time lno selno lnol lnou alnol alnou;
title9 'Pointwise 95% confidence limits using logit(s)';
 
proc print data=two;
 var time elapsed1 n1 d1 p1 sp1 sesp1 h1 py1 seh1;
title8 'Life-table for group 1';
proc print data=two;
 var time spu1 spl1 lnsp1 selnsp1 lnspl1 lnspu1 alnspl1 alnspu1;
title9 'Pointwise 95% confidence limits using s and log(s)';
proc print data=two;
 var time lnlnsp1 selnlns1 lnlnspl1 lnlnspu1 allspl1 allspu1;
title9 'Pointwise 95% confidence limits using log(log(s))';
proc print data=two;
 var time lno1 selno1 lnou1 lnou1 alnol1 alnou1;
title9 'Pointwise 95% confidence limits using logit(s)';
 
proc print data=two;
 var time elapsed2 n2 d2 p2 sp2 sesp2 h2 py2 seh2;
title8 'Life-table for group 2';
proc print data=two;
 var time spu2 spl2 lnsp2 selnsp2 lnspl2 lnspu2 alnspl2 alnspu2;
title9 'Pointwise 95% confidence limits using s and log(s)';
proc print data=two;
 var time lnlnsp2 selnlns2 lnlnspl2 lnlnspu2 allspl2 allspu2;
title9 'Pointwise 95% confidence limits using log(log(s))';
proc print data=two;
 var time lno2 selno2 lnou2 lnou2 alnol2 alnou2;
title9 'Pointwise 95% confidence limits using logit(s)';
 
proc print data=twol;
  var td tpy th td1 tpy1 th1 td2 tpy2 th2;
title8 'Linearized hazards over interval (0, t¨, combined and within each group';
 
data twox; set two;        *** Compute difference in Survival at each time;
if 0<sp1<1 and 0<sp2<1 then do;
  diff=sp1-sp2;                   *** difference in survival probabilities;
   sediff=sqrt(vsp1+vsp2);
   diffu=diff+z95*sediff;
   diffl=diff-z95*sediff;
   kmchi=(sp1-sp2)**2/(vsp10+vsp20);   *** TEST OF SIGNIFICANCE OF DIFFERENCE;
   kmp=1 - probchi(kmchi,1);
  lndiff=lnsp1-lnsp2;               *** log ratio of survival;
  sratio=exp(lndiff);               *** ratio sp1/sp2;
  selndiff=sqrt(vlnsp1+vlnsp2);
   lndiffu=lndiff+z95*selndiff;
   lndiffl=lndiff-z95*selndiff;
   sratiou=exp(lndiffu);
   sratiol=exp(lndiffl);
  lnlndiff=lnlnsp1-lnlnsp2;         *** difference in log(-log(S));
  chratio=exp(lnlndiff);            *** ratio of cumulative hazards;
  selnlndf=sqrt(vlnlnsp1+vlnlnsp2);
   lnlndfu=lnlndiff+z95*selnlndf;
   lnlndfl=lnlndiff-z95*selnlndf;
   chratiou=exp(lnlndfu);
   chratiol=exp(lnlndfl);
  lnor=lno1-lno2;                   *** difference in logit(S);
  survOR = exp(lnor);               *** survival odds ratio;
  selnor=sqrt(vlno1+vlno2);
   lnoru=lnor+z95*selnor;
   lnorl=lnor-z95*selnor;
   survORu=exp(lnoru);
   survORl=exp(lnorl);
end;
 
proc print; var time sp1 vsp1 sp2 vsp2 diff sediff diffl diffu;
title8 'Difference in survival probabilities, se and 95% confidence limits';
proc print; var time sp1 sesp10 n2 sp2 sesp20 kmchi kmp;
title8 'Pointwise tests of difference of Kaplan-Meier estimates using null variance';
proc print; var time lnsp1 vlnsp1 lnsp2 vlnsp2 lndiff selndiff lndiffl lndiffu
                sratio sratiol sratiou;
title8 'Difference in log survival probabilities, se, CI, and survival ratio and CI';
proc print; var time lnlnsp1 vlnlnsp1 lnlnsp2 vlnlnsp2 lnlndiff selnlndf lnlndfl lnlndfu
                chratio chratiol chratiou;
title8 'Difference in log(-log) survival probabilities, se, CI and cum. hazard ratio and CI';
proc print; var time lno1 vlno1 lno2 vlno2 lnor selnor lnorl lnoru survOR survORl survORu;
title8 'Difference in log survival odds, se, CI and survival odds ratio and CI';
 
********  Generate Plots of survival and hazard functions;
 
data twoa; set twol (keep=lastime lastime1 lastime2);
  time = -1; sp=1.0; sp1=1.0; sp2=1.0;
 
data twob; set twoa two (keep = time sp alnlnspu alnlnspl sp1 sp2);
retain lt lt1 lt2 lagtime osl osl1 osl2;
if _n_ = 1 then do;
   lagtime=-1;
   lt = lastime;
   lt1 = lastime1;
   lt2 = lastime2;
   osl=1.0;
   osll=1.0;
   oslu=1.0;
   osl1=1.0;
   osl2=1.0;
   end;
if _n_ > 1 then do;
  t=time;
  sl = sp;
  slu = alnlnspu;
  sll = alnlnspl;
  sl1 = sp1;
  sl2 = sp2;
  output;
  if "%lowcase(&type)" = 'plimit' or "%lowcase(&type)" ¬= 'actuar' then do;
    ta = lagtime + 1;
    if time > ta then do;
      sl = osl;
      slu = oslu;
      sll = osll;
      sl1 = osl1;
      sl2 = osl2;
      tb = time -1;
      do t = ta to tb;
         output;
      end;
    end;
    lagtime=time;
    osl = sp;
    oslu = slu;
    osll = sll;
    osl1 = sp1;
    osl2 = sp2;
  end;
end;
proc sort; by t;
 
data twob; set twob;
     s=sl;
     su = slu;
     sl = sll;
     s1=sl1;
     s2=sl2;
if t gt lt then delete;
if t gt lt1 then s1=.;
if t gt lt2 then s2=.;
 
***proc print; ***var t s s1 s2;
 
 proc plot nolegend;
      plot s*t = '*'
           su*t = '-'
           sl*t = '-'
    / overlay
      vaxis = 0 to 1.0 by .10
      vzero hzero
   ;
 title5
'The survival function in combined group with log(-log) limits';
 
 proc plot nolegend;
      plot s1*t = 'a'
           s2*t = 'b'
    / overlay
      vaxis = 0 to 1.0 by .10
      vzero hzero
   ;
 title5
'Plots of the survival functions in each group';
 
 proc plot nolegend data=two;
      plot h1*time = 'a'
           h2*time = 'b'
    / overlay
      vzero hzero
   ;
 title5
'Plots of the hazard functions in each group';
 
data one; set two
   (keep = n n1 n2 d d1 d2 sp);
 
title5 ' ';
 
****************** tests within time intervals or strata;
 
a=d1; b=d2; c=n1-d1; d=n2-d2;
n=n1+n2; dt=d1+d2;  p=dt/n;
o_e=d1-dt*n1/n;
v=n1*n2*dt*(n-dt)/(n*n*(n-1));   *** Mantel-Haenszel variance;
ls=o_e; lv = v;                  *** logrank test components;
gws = n*ls;  gwv=(n**2)*v;       *** Gehan Wilcoxon;
pws = sp*ls; pwv=(sp**2)*v;      *** Peto-Peto-Prentice Wilcoxon;
clw = sp*v;                *** numerator for correlation of logrank and PPP-W tests;
 
****proc print;
 
*** aggregated over time intervals or strata;
 
proc means sum noprint;
   var ls lv gws gwv pws pwv clw;
   output out=three
   sum= tls tlv tgws tgwv tpws tpwv tclw;
data; set three
  (keep=tls tlv tgws tgwv tpws tpwv tclw);
 
*** tests aggregated over time intervals or strata;
 
   logrank = (tls**2)/tlv;
   plr = 1 - probchi(logrank,1);
   gehan = (tgws**2)/tgwv;
   pgw = 1 - probchi(gehan,1);
   ppp_w = (tpws**2)/tpwv;
   ppw = 1 - probchi(ppp_w,1);
   clw = tclw/sqrt(tlv*tpwv);
   mert = ((sqrt(logrank) + sqrt(ppp_w))**2) / (2*(1+clw));
   pmert = 1 - probchi(mert,1);
 
*** Peto one-step estimates;
 
lnrr=(tls/tlv);
relrisk=exp(lnrr);
vlnrr=1/tlv;
lnrru = lnrr + 1.96*sqrt(vlnrr);
lnrrl = lnrr - 1.96*sqrt(vlnrr);
relrisku = exp(lnrru);
relriskl = exp(lnrrl);
 
lneor=(tpws/tpwv);
eventOR=exp(lneor);
vlneor=1/tpwv;
lneoru = lneor + 1.96*sqrt(vlneor);
lneorl = lneor - 1.96*sqrt(vlneor);
eventORu = exp(lneoru);
eventORl = exp(lneorl);
 
file print n=ps;
put @20 'Test Statistics';
put ' ';
put 'All tests use the Mantel-Haenszel variance';
put ' ';
put 'Logrank test:';
put ' ';
put @10 'Test statistic'
    @50  tls   12.5;
put @10 'Variance      '
    @50  tlv   12.5;
put ' ';
put @10 'Test value (chi-square)'
    @50  Logrank  12.5;
put @10 'p-value'
    @50  plr      12.5;
 
put ' ';
put 'Gehan Wilcoxon test';
put ' ';
put @10 'Test statistic'
    @50  tgws  12.5;
put @10 'Variance      '
    @50  tgwv  12.5;
put ' ';
put @10 'Test value (chi-square)'
    @50  gehan    12.5;
put @10 'p-value'
    @50  pgw      12.5;
 
put ' ';
put 'Peto-Peto-Prentice Wilcoxon test';
put ' ';
put @10 'Test statistic'
    @50  tpws  12.5;
put @10 'Variance      '
    @50  tpwv  12.5;
put ' ';
put @10 'Test value (chi-square)'
    @50  ppp_w    12.5;
put @10 'p-value'
    @50  ppw      12.5;
 
put ' ';
put 'Gastwirth maxi-min efficient robust test (MERT)';
put ' ';
put @10 'Covariance between the logrank and the';
put @10 '     Peto-Peto-Prentice Wilcoxon tests'
    @50  tclw  12.5;
put @10 'Correlation between tests'
    @50  clw   12.5;
put ' ';
put @10 'MERT test value (chi-square)'
    @50  mert     12.5;
put @10 'p-value'
    @50  pmert    12.5;
 
proc print; var lnrr vlnrr lnrrl lnrru relrisk relriskl relrisku;
 title5 'Peto relative risk estimate and confidence limits (for Logrank test)';
 
proc print; var lneor vlneor lneorl lneoru eventOR eventORl eventORu;
 title5 'Peto cumulative event odds ratio and confidence limits (for PPP-Wilcoxon test)';
 
%mend;
 
run;
