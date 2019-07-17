*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter9/swl-is.sas                   *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: computes the robust information sandwich estimate of the  *;
*          covariance matrix of the coefficient estimates for the    *;
*          proportional hazards regression analysis of the Lagakos   *;
*          data presented in Example 9.6.                            *;
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
 
TStart=0; TStop=time;
 
data intimes; set inone;
 
proc sort; by time delta;
 
%let nvar=3;     *** set dimensions for beta vector and covariance matrix;
%let nvar2=9;
%let ivars= age perfstat group;   *** list covariate names;
 
proc phreg covout outest=modfit;
model (TStart, TStop)*delta(0) = &ivars / risklimits COVB;
      output out=Out1 dfbeta=df1-df&nvar xbeta=xbeta/ order=data;
title2 'treatment group adjusted for age and performance status';
title3 'Information Sandwich Robust Estimate of Coefficient Covariance Matrix';
 
   run;
 
   proc iml;
      use Out1;
      read all var{df1 df2 df3} into x;
      v=x` * x;
      reset noname;
      vname={"X1","X2","X3"};
      print,"Estimated Covariance Matrix",, vÝcolname=vname
        rowname=vname format=15.10¨;
      create RCov from vÝcolname=vname rowname=vname¨;
      append from vÝrowname=vname¨;
 
      quit;
 
   run;
 
DATA BETAS; set modfit;
ARRAY BETA {&NVAR} BETA1-BETA&NVAR;
ARRAY X {&NVAR} &IVARS;
if _TYPE_ = 'PARMS';
DO J = 1 TO &NVAR;
  BETA(J) = X(J);
END;
**** PROC PRINT;
 
DATA COVMAT; set modfit;
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
**** PROC PRINT;
 
data covb; set covvec (keep = covs1-covs&nvar2);
data betas; set betas (keep = beta1-beta&nvar);
data dfbetas; set out1 (keep = df1-df&nvar);
 
proc iml;
*** reset print;  *** for test only;
use dfbetas;
read all into x;
covb=x` * x;                          *** robust information sandwich;
use betas;
read all into beta;
use covb;
read all into cov;
beta=beta`;
modcov = shape(cov,&nvar,&nvar);      *** model covariance, inverse information;
inf=inv(modcov);                      *** estimated information;
varb=diag(covb);                      *** robust variances and 95% C.I.;
seb=sqrt(varb);
j=shape(1,&nvar,1);
za=probit(0.975);
cl=beta` - za*j`*seb;
cu=beta` + za*j`*seb;
clrr=exp(cl);
curr=exp(cu);
Waldj=beta`*inv(varb)*(diag(beta));   *** robust Wald covariate tests;
Waldjp=1-probchi(Waldj,1);
betaÝ1,1¨=0;
ic=I(&nvar-1);
zv=shape(0,(&nvar-1),1);
c=(zv¦¦ic);
c=c`;
Waldmod=(c`*beta)`*inv(c`*covb*c)*(c`*beta);   *** Robust Wald model test;
Waldmodp=1-probchi(Waldmod,(&nvar-1));
 
print modcov inf;
title3 'model-based estimated information (Inf) and covariance matrix of estimates (modcov)';
 
print covb;
title3 'Wei-Lin robust information sandwich covariance matrix';
 
print beta seb cl cu clrr curr;
title3 'betas, robust SEs, robust 95% C.I. on betas, 95% C.I. on relative hazards (RR)';
 
print waldj waldjp waldmod waldmodp;
title3 'Wald robust test of covariate effects and p-value, and model test and p-value';
 
quit;
 
run;
 
**************************************************************************************
  Centered score statistic (wo) and robust covariance under null
**************************************************************************************;
 
data events; set intimes;
id=_n_;
if delta=1;
proc sort; by time id;
 
data events; set events; by time id;
if last.time then do;
  output events;
end;
 
proc transpose data=events prefix=t out=times;;
var time;
proc transpose data=events prefix=id out=idtimes;;
var id;
 
data times; merge times idtimes;
 
data inall; if _n_=1 then set times; set intimes;
id=_n_;
 
data terms; set inall end=eof;
%let ntimes=44;                 *** specify number of distinct event times;
%let nxb=132;                   *** p x ntimes;
%let ndata=194;                 *** number of observations;
 
keep xb1-xb&nxb nt1-nt&ntimes;
 
array x{&nvar} &ivars;
array times{&ntimes} t1-t&ntimes;
array xbar{&nvar,&ntimes} xb1-xb&nxb;
array nt{&ntimes} nt1-nt&ntimes;
retain xb1-xb&nxb 0;
retain nt1-nt&ntimes 0;
do l=1 to &ntimes;
   if time >= times(l) then do;
      nt(l)=nt(l)+1;
      do j = 1 to &nvar;
        xbar(j,l) = xbar(j,l) + x(j);
      end;
    end;
end;
 
if eof then do;
  do l=1 to &ntimes;
      do j = 1 to &nvar;
        xbar(j,l) = xbar(j,l)/nt(l);
      end;
  end;
  output;
end;
 
data scores;
 if _N_=1 then set terms (keep=xb1-xb&nxb);
 set inall end=eof;
keep nd uo1-uo&nvar;
array x{&nvar} &ivars;
array xbar{&nvar,&ntimes} xb1-xb&nxb;
array uo{&nvar} uo1-uo&nvar;
array times{&ntimes} t1-t&ntimes;
 
do j = 1 to &nvar;
 uo(j) = 0;
end;
 
nd=id;
 
if delta=1 then do;
  jtime=.;
  do l=1 to &ntimes;
    if time=times(l) then jtime=l;
  end;
  do j = 1 to &nvar;
   uo(j) = x(j) - xbar(j,jtime);
  end;
end;
 
proc sort; by nd;
 
data datan; set inall;
keep nd deltan timen jtimen xn1-xn&nvar;
array x{&nvar} &ivars;
array xn{&nvar} xn1-xn&nvar;
array times{&ntimes} t1-t&ntimes;
 
do j = 1 to &nvar;
   xn(j)=x(j);
end;
 
nd=id;
 
jtime=.;
if delta=1 then do l=1 to &ntimes;
  if time=times(l) then jtime=l;
end;
 
jtimen=jtime;
deltan=delta;
timen=time;
proc sort; by nd;
 
data doubln; set inall;
keep id nd time jtime delta &ivars t1-t&ntimes;
array x{&nvar} &ivars;
array times{&ntimes} t1-t&ntimes;
array idtime{&ntimes} id1-id&ntimes;
 
jtime=.;
if delta=1 then do l=1 to &ntimes;
  if time=times(l) then jtime=l;
end;
 
do nd =1 to &ndata;
  output;
end;
proc sort; by nd id;
 
data alln; merge doubln datan; by nd;
 
data eus;
if _N_=1 then set terms (keep=xb1-xb&nxb nt1-nt&ntimes);
 set alln end=eof; by nd;
keep nd eu1-eu&nvar;
retain eu1-eu&nvar;
array times{&ntimes} t1-t&ntimes;
array ntimes{&ntimes} nt1-nt&ntimes;
array xbar{&nvar,&ntimes} xb1-xb&nxb;
array x{&nvar} &ivars;
array xn{&nvar} xn1-xn&nvar;
array eu{&nvar} eu1-eu&nvar;
 
if first. nd then do j = 1 to &nvar;
 eu(j) = 0;
end;
 
if delta=1 and time <= timen then do;
  do j = 1 to &nvar;
    eu (j) = eu(j) + ( ( xn(j) - xbar(j,jtime) ) / ntimes(jtime) );
  end;
end;
if last.nd then output;
 
data ws; merge scores eus; by nd;
array uo{&nvar} uo1-uo&nvar;
array eu{&nvar} eu1-eu&nvar;
array w{&nvar} w1-w&nvar;
 
do j = 1 to &nvar;
  w(j) = uo(j) - eu(j);
end;
 
DATA THREE; set ws end=eof;
ARRAY w {&NVAR} w1-w&NVAR;
ARRAY wo {&NVAR} wo1-wo&NVAR;
retain wo1-wo&NVAR 0;
ARRAY wprodo {&NVAR, &NVAR} wprodo1-wprodo&NVAR2;
RETAIN wprodo1-wprodo&NVAR2 0;                        *** outer null score products;
 
do j = 1 to &nvar;
  wo(j) = wo(j) + w(j);
do k = 1 to &nvar;
  wprodo(j,k) = wprodo(j,k) + w(j)*w(k);
end;
end;
 
if eof then output;
 
data ws; set three (keep = wo1-wo&nvar);
data wprodos; set three (keep = wprodo1-wprodo&nvar2);
 
 
proc iml;
***reset print;  *** for test only;
use ws;
read all into wo;
wo = wo`;
use wprodos;
read all into Wprodo;
 
 wpo = shape(wprodo,&nvar,&nvar);
 score = wo`*inv(wpo)*wo;                  *** robust score test;
 scorep=1-probchi(score,(&nvar-1));
 
print 'Centered score vector (wo) and robust covariance (wpo) under the model null hypothesis'
  wo wpo;
print 'Robust model score test and p-value'
  score scorep;
 
quit;
 
RUN;
 
 
RUN;
