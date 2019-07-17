*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter9/hypotime.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: This program generate the data set hytimes that contains a*;
*          single observation with the numbers at risk and numbers of*;
*          events at each distinct event time. This data set is then *;
*          used with the macro smooth in the program kernelsm.sas to *;
*          generate kernel smoothed estimates of the intensity       *;
*          function over time, or the macro AGtests in AGtests.sas to*;
*          compute the estimates of the cumulative intensity         *;
*          functions and the Aalen-Gill tests for recurrent event    *;
*          processes. This program uses the input data set hyevents  *;
*          that has one observation for each event for each subject, *;
*          or a single observation if no events.                     *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
**********************************************************************;
 
* The output data set will contain a single observation with the following
* arrays of variables:
*
*  t1-t&ntime   =  the set of distinct ordered event times;
*  xe1-xe&ntime =  the number of events at each time in group 1;
*  xc1-xc&ntime =  the number at risk at each time in group 1;
*  ye1-ye&ntime =  the number of events at each time in group 2;
*  yc1-yc&ntime =  the number at risk at each time in group 2;
*  y1-y&ntime   =  the total number at risk at each time;
*  maxj = number of elements in the arrays ;
 
title "dcct -hypoglycemia analyses for book";
 
libname hypobook '/jml/biostatmethods/datasets/hypoglycemia';
 
data events; set hypobook.hyevents;
%let ntime=1565;
%let group=intgroup;
%let id = patient;
 
**********************************************************************;
 
proc sort data=events; by etime;
 
data rank; set events;
  retain rank pretime 0;
  keep retime;
  if etime eq . then do;
    retime=etime;
    end;
  if etime gt . then do;
    retime=rank;
    if etime gt pretime then do;
      rank=rank+1;
      retime=rank;
      pretime=etime;
      end;
    end;
 
data ranktime;
  merge events rank;
*------add the event time array to each record------------------------;
proc sort data=ranktime;
  by &id etime;
 
*------ dimension arrays at the number of distinct event times -------;
 
data one;
  set ranktime end=last;
  retain t1-t&ntime maxj 0;
  keep t1-t&ntime maxj;
  array t(&ntime) t1-t&ntime;
  if retime ne . then do;
    j=retime;
    maxj=max(maxj,j);
    t(j)=etime;
    end;
  if last then output;
 
data two; if _n_=1 then set one; set ranktime;
 
*------determine the x arrays and the y arrays------------------------;
 
data hypobook.hytimes;
 
  set two end=last;
  by &id etime;
  retain xe1-xe&ntime xc1-xc&ntime
         ye1-ye&ntime yc1-yc&ntime
         y1-y&ntime 0;
  keep xe1-xe&ntime xc1-xc&ntime
       ye1-ye&ntime yc1-yc&ntime
       y1-y&ntime t1-t&ntime maxj;
  array t(&ntime) t1-t&ntime;
  array xe(&ntime) xe1-xe&ntime;
  array xc(&ntime) xc1-xc&ntime;
  array ye(&ntime) ye1-ye&ntime;
  array yc(&ntime) yc1-yc&ntime;
  array y(&ntime) y1-y&ntime;
  if &group eq 1 then do;
    if retime gt . then do;
      j=retime;
      xe(j)=xe(j)+1;
      end;
    if last.&id then do;
      do j=1 to maxj;
        if ftime ge t(j) then do;
          ye(j)=ye(j)+1;
          y(j)=y(j)+1;
          end;
        end;
      end;
    end;
  if &group eq 2 then do;
    if retime gt . then do;
      j=retime;
      xc(j)=xc(j)+1;
      end;
    if last.&id then do;
      do j=1 to maxj;
        if ftime ge t(j) then do;
          yc(j)=yc(j)+1;
          y(j)=y(j)+1;
          end;
        end;
      end;
    end;
  if last then output;
 
proc print;
  var  t1-t&ntime;
proc print;
  var  xe1-xe&ntime;
proc print;
  var  xc1-xc&ntime;
proc print;
  var  ye1-ye&ntime;
proc print;
  var  yc1-yc&ntime;
proc print;
  var  y1-y&ntime;
 
run;
