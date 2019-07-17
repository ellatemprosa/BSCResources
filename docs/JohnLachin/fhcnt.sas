*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter8/fhcnt.sas                    *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: reads the data file containing the detailed information   *;
*          on the recurrent infections from Fleming and Harrington   *;
*          (1991) and computes the numbers of events and exposure    *;
*          time for each individual. Additional computations can     *;
*          then be performed.                                        *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
title1 'CGD data from Fleming and Harrington (1991)';
 
filename fhcgd '/jml/biostatmethods/datasets/fh-cgd.dat';
 
data one;
infile fhcgd;
input id rdt idt z1 z2 z3 z4 z5 z6 z7 z8 z9 t1 t2 d s;
trt=0; if z1=1 then trt=1;
delta=0; if d=1 then delta=1;
time=t1;
inherit=z2-1;
 
proc sort; by id time;
 
data counts; set one; by id time;
if last.id;
if s=1 then nevents=delta;        *** one record for this subject, either last time an
                                      event or censoring;
if s > 1 then do;
   if delta=0 then nevents=s-1;   *** prior events but last time a censoring time;
   if delta=1 then nevents=s;     *** prior events, last time an event time;
end;
 
ievents=0; if nevents>1 then ievents=1;
futime=time;
 
proc print; var id z1 z2 z3 z4 z5 z6 z7 z8 z9 nevents futime;
title2 'data listing used for book data set counts';
 
proc freq; tables trt*nevents / nocol nopercent;
proc sort; by trt;
proc means; var futime; by trt;
 
run;
