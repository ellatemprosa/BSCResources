*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter9/swl-ph.sas                   *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: conducts the various proportional hazards regression      *;
*          analyses presented in Example 9.6, including those with   *;
*          covariate interactions, nested effects, and stratified    *;
*          analyses.                                                 *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
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
 
data inone; set bioweb.lagakos;
fail = cause; if cause=2 then fail=1;
delta=fail;
 
group=1-treatmnt;   *** 1=A, 0=B;
 
ageperf=age*perfstat;    *** interaction variables;
agegroup=age*group;
perfgrp=perfstat*group;
 
data intimes; set inone;
 
proc phreg;
model time*delta(0) = group / risklimits;
title2 'groupment un-adjusted for covariates';
 
proc phreg;
model time*delta(0) = age perfstat / risklimits;
title2 'overall covariate effects';
 
proc phreg;
model time*delta(0) = age perfstat group / risklimits;
title2 'group adjusted for age and performance status';
 
proc phreg;
model time*delta(0) = age perfstat group ageperf agegroup perfgrp / risklimits;
title2 'coefficient interaction model';
 
proc phreg;
model time*delta(0) = age perfstat group0 group1 / risklimits;
group0= (1-perfstat)*group;
group1= perfstat*group;
title2 'group coefficient nested within perfstat model';
 
proc phreg;
model time*delta(0) = age group / risklimits;
strata perfstat;
title2 'stratified by performance status with homogeneous covariate effects';
 
proc phreg;
model time*delta(0) = age0 age1 group0 group1 / risklimits;
strata perfstat;
age0= (1-perfstat)*age;
age1= perfstat*age;
group0= (1-perfstat)*group;
group1= perfstat*group;
title2 'stratified by performance status with heterogeneous covariate effects';
 
proc phreg;
model time*delta(0) = age group0 group1 / risklimits;
strata perfstat;
group0= (1-perfstat)*group;
group1= perfstat*group;
title2 'stratified by performance status with heterogeneous group effect';
 
run;
