*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter9/swl-anal.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: uses the macro %survanal to perform survival analyses of  *;
*          the Lagakos squamous cell data in the complete cohort,    *;
*          and within subgroups defined by performance status. The   *;
*          analysis within the non-ambulatory subgroup is presented  *;
*          in Example 9.1 (Tables 9.1-3 and Figure 9.1) and Example  *;
*          9.4.                                                      *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
  filename survanal '/jml/biostatmethods/chapter9/survanal.sas';
 
  %include survanal;
 
*** The data file is containted in the datasets directory;
 
  libname bioweb '/jml/biostatmethods/datasets';
 
Title1 'Lagakos Squamous Cell Carcinoma Data';
 
** 194 patients with squamous cell carcinoma:
  perfstat:0=ambulatory, 1=non-ambulatory performance status
  treatmnt: 0=A, 1=B
  age in years
  time in weeks
  cause of failure: 0=censored, 1=local spread of disease, 2=metastatic spread,
  fail (defined below): spread of disease, either local or metastatic
;
 
data intwo; set bioweb.lagakos;
fail = cause; if cause=2 then fail=1;
delta=fail;
group=treatmnt+1;
 
data intimes; set intwo;
keep delta time group;
 
%tables;
%let type=plimit; %tests;
 
run;
 
data intimes; set intwo;
keep delta time group;
 
if perfstat=0;
 
title2 'subgroup with ambulatory status';
%tables;
%let type=plimit; %tests;
 
run;
 
data intimes; set intwo;
keep delta time group;
 
if perfstat=1;
 
title2 'subgroup with non-ambulatory status';
%tables;
%let type=plimit; %tests;
 
run;
