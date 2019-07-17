*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter9/fh-cgdtm.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: uses the macro <I>%timeset</I> in <I>AGtimes.sas</I> to   *;
*          generate the data set <I>CGDtimes</I> in a format         *;
*          required to compute the kernel smoothed intensity         *;
*          estimates and the Aalen-Gill tests. This data set may     *;
*          then be used with the other macros and programs to        *;
*          perform additional analyses as specified in Problem       *;
*          9.18.                                                     *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
title1 'CGD data from Fleming and Harrington (1991)';
 
 
 Filename agtimes '/jml/biostatmethods/chapter9/AGtimes.sas';
 %include agtimes;
 
 Filename cgddata '/jml/biostatmethods/datasets/FH-CGD.dat';
 
 libname biosbook '/jml/biostatmethods/datasets';
 
*** required macro variables;
 
    %let indata= cgddata;   **** the input data set;
    %let outdata= biosbook.cgdtimes;  **** the output data set;
    %let group=z1;
    %let id=id;
 
data cgddata; infile cgddata;
input id rdt idt z1 z2 z3 z4 z5 z6 z7 z8 z9 t1 t2 d s;
delta=0; if d=1 then delta=1;
time=t1;
 
run;
 
%timeset;
 
run;
