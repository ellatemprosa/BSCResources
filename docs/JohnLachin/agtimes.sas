*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter9/agtimes.sas                  *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: is a more general macro that can be used to generate the  *;
*          appropriate data set with the arrays for event times,     *;
*          numbers of events and numbers at risk from any data set   *;
*          with recurrent events. In <I>FG-CGDtm.sas</I>, this       *;
*          program is applied to the Fleming and Harrington (1991)   *;
*          CGD data described in Problem 9.18. The macro requires    *;
*          that the data set contain a patient id variable, a group  *;
*          variable (1=experimental, 2=control), both of which are   *;
*          specified as macro variables, and additional variables    *;
*          to represent an event time (<I>etime</I>), the maximum    *;
*          follow-up time (<I>ftime</I>), and an indicator           *;
*          <I>delta</I> to indicate whether an event was observed    *;
*          at the current time (1) or the observation is right       *;
*          censored (0) at that time.                                *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
**********************************************************************;
 
* The input data set must contain the following three variables with at least
* one observation per subject, and one observation per event per subject.
*   an id number specified in the macro variable &id;;
*   group =1 if experimental or (+), 2 if control or (-);
*   time = time of an event or right censoring;
*   delta = 1 if an event at this time or zero if censored;
* If a subject had no events then the record would specify the censoring time.
* If a subject had k events and was later censored then the first k records
*   would specify the event times and the k+1th record the censoring time;
* If the subject was no longer at risk after the kth event time, then there
*   would only be k records as above;
*
*Usage:
*
* The macro variables below specify the input and output data set names, the
* variable representing treatment group coded as (1,2), the id variable
* designating each distinct subject.;
*
* Filename or libname statements are also required as appropriate;
*
* The required macro variables must be defined such as:
*
*    %let indata= cgddata;   **** the input data set;
*    %let outdata= biosbook.cgdtimes;  **** the output data set;
*    %let group=z1; *** group variable coded 1 or 2;
*    %let id=id;    *** subject ID variable;
*
* The output data set will contain a single observation with the following
* arrays of variables:
*
*  t1-t&ntimes   =  the set of distinct ordered event times;
*  xe1-xe&ntimes =  the number of events at each time in group 1;
*  xc1-xc&ntimes =  the number at risk at each time in group 1;
*  ye1-ye&ntimes =  the number of events at each time in group 2;
*  yc1-yc&ntimes =  the number at risk at each time in group 2;
*  y1-y&ntimes   =  the total number at risk at each time;
*  maxj = length of the arrays;
* ;
 
%macro timeset;
 
proc sort data=&indata; by &id time;
 
data events; set &indata;
event=delta;
etime=.; if event=1 then etime=time;
ftime=time;
 
keep &id &group event etime ftime;
 
proc means noprint; by &id;
 var event; output out=nevents sum=nevents;
proc means data=nevents;
 var nevents;
proc freq data=nevents; tables nevents;
title2 'distribution of number of events per subject';
 
data rfe; set events;
 
run;
 
*---------add total number of events per &id to each record-------;
proc sort data=rfe;
  by &id etime;
 
data events;
  merge nevents (keep=&id nevents) rfe;
  by &id;
  if first.&id then do;
    evnum=0;
    end;
  if event eq 1 then evnum+1;
run;
 
*--------summarize randomization and followup with events data--------;
proc sort data=events;
  by &id etime;
 
data kevents; set events;
proc sort; by etime;
 
data count; set kevents;
keep etime ftime;
 
data count; set count end=eof; by etime;
  retain count 0;
  keep count;
  if etime gt . then do;
  if last.etime then do;
      count=count+1;
    end; end;
  if eof then do;
     output;
     call symput('ntimes',compress(put(count,8.0)));
  end;
proc print;
title2 'number of distinct event times';
 
run;
 
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
  retain t1-t&ntimes maxj 0;
  keep t1-t&ntimes maxj;
  array t(&ntimes) t1-t&ntimes;
  if retime ne . then do;
    j=retime;
    maxj=max(maxj,j);
    t(j)=etime;
    end;
  if last then output;
 
data two; if _n_=1 then set one; set ranktime;
 
*------determine the x arrays and the y arrays------------------------;
 
data &outdata;
 
  set two end=last;
  by &id etime;
  retain xe1-xe&ntimes xc1-xc&ntimes
         ye1-ye&ntimes yc1-yc&ntimes
         y1-y&ntimes 0;
  keep xe1-xe&ntimes xc1-xc&ntimes
       ye1-ye&ntimes yc1-yc&ntimes
       y1-y&ntimes t1-t&ntimes maxj;
  array t(&ntimes) t1-t&ntimes;
  array xe(&ntimes) xe1-xe&ntimes;
  array xc(&ntimes) xc1-xc&ntimes;
  array ye(&ntimes) ye1-ye&ntimes;
  array yc(&ntimes) yc1-yc&ntimes;
  array y(&ntimes) y1-y&ntimes;
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
  var  t1-t&ntimes;
proc print;
  var  xe1-xe&ntimes;
proc print;
  var  xc1-xc&ntimes;
proc print;
  var  ye1-ye&ntimes;
proc print;
  var  yc1-yc&ntimes;
proc print;
  var  y1-y&ntimes;
 
 
%mend;
 
run;
