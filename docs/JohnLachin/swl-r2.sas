*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter9/swl-r2.sas                   *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: computes the Schemper and Kent-O'Quigley R-square         *;
*          measures of explained variation in the proportional       *;
*          hazards regression analysis of the Lagakos data           *;
*          presented in Example 9.6.                                 *;
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
 
data intimes; set inone;
 
%let nvar=3;
%let nvar2=9;
%let ivars= age perfstat group;
 
proc lifetest noprint outs=kmest;
time time*delta(0);
title2 'background survival function estimate';
 
proc phreg data=intimes OUTEST=modfit;
model time*delta(0) = &ivars / risklimits;
output out=PHsurv xbeta=xb survival=sPH;
title2 'PH model with treatment adjusted for age and performance status';
 
**************** Schemper measure v2 **************************;
 
***proc print data=modfit;
***proc print data=phsurv;
 
data kmest; set kmest;
if _censor_=0;
if survival=1 then delete;
if survival=. then delete;
 
data PHss; set Phsurv;
if delta=1;
 
proc sort data=kmest; by time;
proc sort data=PHss; by time;
 
%let ntimes=44;   *** specify number of distinct event times;
 
data kms; set kmest end=eof; by time;
keep s1-s&ntimes t1-t&ntimes ns;                *** ns = # event times;
retain s1-s&ntimes t1-t&ntimes ns 0;
array kms (&ntimes) s1-s&ntimes;
array ts (&ntimes) t1-t&ntimes;
if survival=. then survival=0;
if first.time then do;
  ns=ns+1;
  ts(ns)=time;
  kms(ns) = survival;
end;
if eof then output;
 
*** proc print;
 
data PHs; set PHss end=eof; by time;
keep PHs1-PHs&ntimes PHt1-PHt&ntimes nphs;
retain PHs1-PHs&ntimes PHt1-PHt&ntimes nphs 0;
array PHs (&ntimes) PHs1-PHs&ntimes;
array PHts (&ntimes) PHt1-PHt&ntimes;
if first.time then do;
  nphs=nphs+1;
  phts(nphs)=time;
  s0=sph**(exp(-xb));
  PHs(nphs) = s0;
end;
if eof then output;
 
*** proc print;
 
data ss; merge phs kms;
 
data survs; if _n_=1 then set ss;
 set PHsurv;
array PHso (&ntimes) PHs1-PHs&ntimes;
array kms (&ntimes) s1-s&ntimes;
array ts (&ntimes) t1-t&ntimes;
retain d0 d1 0;
 
************ Schemper's v2 as described by Schemper and Stare (1996), not Schemper (1990, 92);
 
if delta=1 then do;
  k=ns;
  di0=0; di1=0;
  do j = 1 to ns;
     sij = 1; if time <= ts(j) then sij=0;
     phs = phso(j)**(exp(xb));
     di1 = di1 + (sij - phs)**2;
     di0 = di0 + (sij - kms(j))**2;
  end;
end;
if delta=0 then do;
  k=0;
  di0=0; di1=0;
  do j = 1 to ns;
    sij = 1;
    if time >= ts(j) then do;
     k = k + 1;
     phs = phso(j)**(exp(xb));
     di1 = di1 + (sij - phs)**2;
     di0 = di0 + (sij - kms(j))**2;
    end;
  end;
end;
d1 = d1 + (di1/k);
d0 = d0 + (di0/k);
 
data v2; set survs end=eof;
if eof;
v2= 1 - d1/d0;
 
proc print; var d1 d0 v2;
title3 'Schemper''s measure v2 as described by Schemper and Stare (1996)';
 
 
******************  Kent-O'Quigley approximate measures *****************************;
 
proc corr data=intimes cov outp=cov noprint; var age perfstat group;
 
run;
 
DATA BETAS; set modfit;
ARRAY BETA {&NVAR} BETA1-BETA&NVAR;
ARRAY X {&NVAR} &IVARS;
if _TYPE_ = 'PARMS';
DO J = 1 TO &NVAR;
  BETA(J) = X(J);
END;
***PROC PRINT;
 
DATA COVMAT; set cov;
if _TYPE_ = 'COV';
 
DATA COVVEC; SET COVMAT END=EOF;
KEEP COVS1-COVS&NVAR2;
RETAIN COVS1-COVS&NVAR2;
ARRAY COVS {&NVAR, &NVAR};
ARRAY X {&NVAR} &IVARS;
I = _N_;
DO J = 1 TO &NVAR;
  COVS(I,J) = X(J);
END;
IF EOF THEN OUTPUT COVVEC;
***PROC PRINT;
 
data covx; set covvec (keep = covs1-covs&nvar2);
data betas; set betas (keep = beta1-beta&nvar);
 
title3 'Kent-O''Quigley approximate measures';
title4 'R2WA is the model R2, R2WJ is the diagonal matrix of covariate partial R2 values';
 
proc iml;
*** reset print;  *** for test only;
use betas;
read all into beta;
use covx;
read all into cov;
beta=beta`;
var = shape(cov,&nvar,&nvar);
A=beta`*VAR*beta;
R2WA = A/(A+1);   *** Model R2;
vari=inv(var);
ivj=diag(vari);
vj=inv(ivj);
betaj=diag(beta);
aj=betaj*vj*betaj;
aj1=diag(aj+1);
r2wj=aj*inv(aj1);  *** Partial R2 explained by each covariate;
 
print R2WA R2WJ;
 
quit;
 
RUN;
 
RUN;
