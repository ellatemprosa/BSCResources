*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter9/AG_tests.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: This macro computes the estimated cumulative intensity    *;
*          for possibly recurrent event times in two groups as       *;
*          presented in (9.132). These are plotted.
*          The program then computes the Aalen-Gill test statistics  *;
*          for possibly recurrent event times with an allowance for  *;
*          ties, as in (9.149)-(9.150). The logrank and Gehan-       *;
*          Wilcoxon tests are presented.                             *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
*
*  Usage:
*
*  %AGtests(data=input data set name, last is the default);
*
*
*  The input data set must consist of a single observation with 5 arrays
*  of variables and a variable maxj that specifies the number of elements
*  in these arrays:
*    t1-tmaxj =  Distinct times at which events are observed, t1=first
*       event time, t2=second, etc.
*    xe1-xemaxj = number of events observed at each event time in the
*       experimental group (1).
*    ye1-yemaxj = number of subjects at risk in the experimental group (1)
*       at each event time.
*    xc1-xcmaxj = numbers of events in the control group (2).
*    yc1-ycmaxj = numbers at risk in the control group (2).
*
*  The unit of time is arbitrary.
*
*  See the program hypotest.sas for an example.
*
*  The input data set may be constructed using the program timessg.sas.
*;
 
 
%macro AGtests(data=_last_);
 
*--------------------------------------------------------------------------------;
*  compute the estimated intensity and integrated intensit at each event time;
*--------------------------------------------------------------------------------;
 
data intimes;
 set &data;
 
call symput('ntime',compress(put(maxj,8.0)));
 
data tplot; set intimes;
 
  retain presre presrc presr epretime cpretime 0;
  keep re rc cre crc sre src time srey srcy sr sry;
 
  array xe(&ntime) xe1-xe&ntime;
  array xc(&ntime) xc1-xc&ntime;
  array ye(&ntime) ye1-ye&ntime;
  array yc(&ntime) yc1-yc&ntime;
  array t(&ntime) t1-t&ntime;
  do j=1 to maxj;
    time=t(j);
    re=.;
    if xe(j) > 0 then do;
      ltime = time - epretime;
      re= (xe(j)/ye(j)) / ltime; *** linearized intensities, not used;
      epretime=time;
      end;
    rc=.;
    if xc(j) > 0 then do;
      ltime = time - cpretime;
      rc= (xc(j)/yc(j)) / ltime;
      cpretime=time;
      end;
    sre=presre+(xe(j)/ye(j));
    src=presrc+(xc(j)/yc(j));
    cre=sre;
    crc=src;
    if xe(j)=0 then cre=.;
    if xc(j)=0 then crc=.;
    presre=sre;
    presrc=src;
    sr=presr+ ( (xc(j)+xe(j)) / (yc(j)+ye(j)) );
    presr=sr;
    srey=ye(j);
    srcy=yc(j);
    sry = ye(j) + yc(j);
    output;
    end;
/*;
proc print data=tplot;
  var re rc sre src sr time srey srcy sry;
 
** these tend to be poor plots and are not usually helpful;
proc plot data=tplot;
  plot  re*time='e'  rc*time='s'/overlay;
title2 'linearized intensity between successive event times';
*/;
 
proc plot data=tplot;
  plot cre*time='e' crc*time='c'/
       overlay;
title2 'cumulative intensity rates for group 2 (c) and group 1 (e)';
 
*--------------------------------------------------------------------------------;
* compute Aalen-Gill test statistics for possibly recurrent event times;
*--------------------------------------------------------------------------------;
 
data test;
 set intimes;
 
  keep maxj tlsum tlvar twsum twvar x2s x2w x2sp x2wp;
  array xe(&ntime) xe1-xe&ntime;
  array xc(&ntime) xc1-xc&ntime;
  array ye(&ntime) ye1-ye&ntime;
  array yc(&ntime) yc1-yc&ntime;
  array y(&ntime) y1-y&ntime;
 
  tlsum=0;
  tlvar=0;
  twsum=0;
  twvar=0;
 
  do j=1 to maxj;
    xt=xe(j)+xc(j);
    yt=ye(j)+yc(j);
    if ye(j)*yc(j) > 0 then do;
      tl= xe(j) - (xt*ye(j)/yt);
      tlsum = tlsum + tl;
      tw= yt*tl;
      twsum = twsum + tw;
      tlv = ye(j)*yc(j)*xt*(yt-xt);
      tlv = tlv/((yt**2)*(yt-1));
      tlvar = tlvar + tlv;
      twvar = twvar + (yt**2)*tlv;
    end;
  end;
 
  x2s=(tlsum**2)/tlvar;
  x2w=(twsum**2)/twvar;
  x2sp=1-probchi(x2s,1);
  x2wp=1-probchi(x2w,1);
 
proc print;
  var maxj tlsum tlvar x2s x2sp;
title2 'logrank test numerator (sum), variance, chi-square and p-value';
proc print;
  var maxj twsum twvar x2w x2wp;
title2 'Gehan Wilcoxon test numerator (sum), variance, chi-square and p-value';run;
 
%mend;
 
run;
