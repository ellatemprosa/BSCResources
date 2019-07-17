*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter9/dcctneph.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: . This program uses the data set nephdata.dat with the    *;
*          macros provided in survanal.sas to perform the modified   *;
*          Kaplan-Meier lifetable analyses presented tables 9.5-5    *;
*          of Example 9.2 and to compute the tests presented in      *;
*          Example 9.5. The dataset contains additional variables    *;
*          not mentioned in the book and could be used for           *;
*          supplemental exercises.                                   *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
  filename survanal '/jml/biostatmethods/chapter9/survanal.sas';
 
  %include survanal;
 
*** The data file is containted in the datasets directory;
 
filename nephdata 'c:\sasjobs\bioweb\datasets\nephdata.dat';
 
data one; infile nephdata;
input patient int primary
etdpatb neur aer0 neph2flg neph2vis duration female age adult bcval5 hbael bmi
mhba1-mhba9;
 
if primary=0;  *** select secondary cohort only;
if aer0<40;  *** no microalbuminuria on entry;
 
data intimes; set one;
keep delta time group;
 
delta=neph2flg; *** indicator variable for microalbuminuria or worse at this visit;
time=neph2vis;  *** quarterly visit number (evaluations only done annually);
group=int+1; *** 1=conventional, 2=intensive;
if time > 36 then time = 40;
time=time/4;  *** convert quarterly visit number to year of followup;
 
title2 'Microalbuminuria in the DCCT Secondary Cohort';
 
proc sort; by descending time;
 
%TABLES;
%LET type=PLIMIT; %TESTS;
 
run;
