*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter9/KernelSm.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: This macro computes the kernel smoothed estimate of a     *;
*          hazard function or intensity function for a counting      *;
*          process based on possibly recurrent event time data.      *;
*          The macro uses the Epanechnikov kernel.  The estimate,    *;
*          its variance and the kernel are presented in equations    *;
*          (9.137)-(9.139).  The program allows for computations for *;
*          two groups of subjects.  See Usage below.                 *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*---------- ----------------------------------------------------------*;
*
*  Usage:
*
   %smooth(data= data set name or leave out to use the last as default,
      intime= either day, week, month or year,
      width = the band width in same units as intime,
      ratetime=day, month or year,
      rateK= multiplier, eg 100 for per 100 patient years,
      plottime=time asis for plots (day, month or year),
      width = the band width,
      axis = axis specifications for the plots such as
             haxis = 0 to 50 by 10 vaxis= l to 10 by 1
      );
*
*  The input data set must consist of a
*  single observation with 5 arrays of variables with a variable maxj
*  containing the number of elements (times)
*    t1-tmaxj =  Distinct times at which events are observed, t1=first
*       event time, t2=second, etc.
*    xe1-xemaxj = number of events observed at each event time in the
*       experimental group (1).
*    ye1-yemaxj = number of subjects at risk in the experimental group (1)
*       at each event time.
*    xc1-xcmaxj = numbers of events in the control group (2).
*    yc1-ycmaxj = numbers at risk in the control group (2).
*
*  See the program hypokrnl.sas for an example.
*
*  The input data set may be constructed using the program timessg.sas.
*;
 
%macro smooth(data=_last_,intime=,width=,ratetime=,rateK=,
      plottime=,axis =);
 
data inone;
 set &data;
 
call symput('ntimes',compress(put(maxj,8.0)));
 
data inone; set inone;
keep time fte vte ftc vtc;
 
  array t(&ntimes) t1-t&ntimes;
  array xe(&ntimes) xe1-xe&ntimes;
  array xc(&ntimes) xc1-xc&ntimes;
  array ye(&ntimes) ye1-ye&ntimes;
  array yc(&ntimes) yc1-yc&ntimes;
 
tmax=t(maxj);
 
band = &width;
 
qj=1;
 
do time = 0 to tmax;
  fte=0; ftc=0;
  vte=0; vtc=0;
  q=0; aq=0; bq=0;
  nqj=qj;
 
*** kernel smoothing;
 
 do j = qj to maxj;
   arg = (time - t(j))/band;
    if arg >= 1 then nqj=j;
    knl=0; aq=0; bq=0; gq=0;
    if -1  < arg < 1 then do;
      knl =  0.75 * (1 - arg**2);
      if time < band then do;       **** tail correction, abgk, 1993, p. 251;
        q=time/band;
        gq=0.75*( ((2/15 + (q**3)/3 - (q**5)/5)*(2/3 + q - (q**3)/3))
           - (((1-q**2)**4)/16) );
        aq=(2/15 + (q**3)/3 - (q**5)/5)/gq;
        bq=((1-q**2)**2)/(4*gq);
        knl=knl*(aq + bq*arg);
      end;
      if time > (tmax-band) then do;
        q=(tmax-time)/band;
        gq=0.75*( ((2/15 + (q**3)/3 - (q**5)/5)*(2/3 + q - (q**3)/3))
          - (((1-q**2)**4)/16) );
        aq=(2/15 + (q**3)/3 - (q**5)/5)/gq;
        bq=((1-q**2)**2)/(4*gq);
        knl=knl*(aq - bq*arg);
      end;
        if ye(j) > 0 then do;
          fte = fte + (knl * xe(j) / ye(j));
          vte = vte + ( knl*knl * xe(j)*(ye(j)-xe(j)) / (ye(j)**3) );
        end;
        if yc(j) > 0 then do;
          ftc = ftc + (knl * xc(j) / yc(j));
          vtc = vtc + ( knl*knl * xc(j)*(yc(j)-xc(j)) / (yc(j)**3) );
        end;
     end;
     if arg <= -1 then go to out1;
  end;
 
out1:;
qj = nqj;
 
fte = fte/band;
vte = vte/(band**2);
ftc = ftc/band;
vtc = vtc/(band**2);
 
output inone;
end;
 
data intwo; set inone;
 
seue = fte + (1.96*sqrt(vte));
sele = fte - (1.96*sqrt(vte));
seuc = ftc + (1.96*sqrt(vtc));
selc = ftc - (1.96*sqrt(vtc));
 
****  convert time and rates to years if time not in years;
 
If "%lowcase(&intime)"="day" then XT=365.25;
If "%lowcase(&intime)"="week" then XT=52;
If "%lowcase(&intime)"="month" then XT=12;
If "%lowcase(&intime)"="year" then XT=1;
  year = time / xt;
  fte = fte * xt;
  sele = sele * xt;
  seue = seue * xt;
  ftc = ftc * xt;
  selc = selc * xt;
  seuc = seuc * xt;
****  convert estimates to rates in desired time units;
If "%lowcase(&ratetime)"="day" then XT=365.25;
If "%lowcase(&ratetime)"="week" then XT=52;
If "%lowcase(&ratetime)"="month" then XT=12;
If "%lowcase(&ratetime)"="year" then XT=1;
  fte = fte * &rateK / xt;
  sele = sele * &rateK / xt;
  seue = seue * &rateK /xt;
  ftc = ftc * &rateK /xt;
  selc = selc * &rateK /xt;
  seuc = seuc * &rateK /xt;
****;
 
if sele < 0 then sele = 0;
if selc < 0 then selc = 0;
 
**** convert time to desired time units for plots;
If "%lowcase(&plottime)"="day" then XT=365.25;
If "%lowcase(&plottime)"="week" then XT=52;
If "%lowcase(&plottime)"="month" then XT=12;
If "%lowcase(&plottime)"="year" then XT=1;
time = year*xt;
****;
 
 proc plot;
  plot fte*time='e'
       ftc*time='c'
  /overlay &axis;
title3 "smoothed estimate of the intensity function -- all patients";
title4 "for group 1 (e) and group 2 (c)";
title5 "rate per &ratek &ratetime(s) versus &plottime";
 
 proc plot;
  plot fte*time='e' seue*time='-' sele*time='-'
       ftc*time='s' seuc*time='-' selc*time='-'
  /overlay &axis;
title6 'with 95% confidence band presented';
 
run;
%mend;
 
run;
