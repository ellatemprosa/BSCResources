*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter7/logirob.sas                  *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: computes the robust information sandwich estimate of the  *;
*          covariance matrix of the coefficient estimates and        *;
*          robust confidence intervals and Wald tests for the renal  *;
*          data as shown in Example 7.8. This also computes the      *;
*          score vector and the information sandwich estimate under  *;
*          the null hypothesis, and the robust efficient score       *;
*          test. First run Renal.sas to set up the renal data set.   *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
 
data one; set renal;
y= micro24;
 
data inall; set one;
if y ne .;
n=1;
subject = _n_;
 
%let nvar=6;    *** including the intercept;
%let nvar2=36;
%let ivars= int hbael yearsdm sbp female;
 
PROC Genmod data=inall;
   class subject;
   MODEL Y/n = &ivars / COVB link=logit dist=bin;
   repeated subject=subject / type=unstr covb;
TITLE2 'LOGISTIC REGRESSION MODEL through genmod with Robust Analysis';
 
PROC LOGISTIC data=inall DESCENDING COVOUT OUTEST=modfit;
   MODEL Y = &ivars / COVB;
   OUTPUT OUT=TWO XBETA=XBETA PRED=PRED;
TITLE2 'LOGISTIC REGRESSION MODEL and robust computations using IML';
 
DATA BETAS; set modfit;
ARRAY BETA {&NVAR} BETA1-BETA&NVAR;
ARRAY X {&NVAR} INTERCEP &IVARS;
if _TYPE_ = 'PARMS';
DO J = 1 TO &NVAR;
  BETA(J) = X(J);
END;
**PROC PRINT;
 
DATA COVMAT; set modfit;
if _TYPE_ = 'COV';
 
DATA COVVEC; SET COVMAT END=EOF;
KEEP COVS1-COVS&NVAR2;
RETAIN COVS1-COVS&NVAR2;
ARRAY COVS {&NVAR, &NVAR};
ARRAY X {&NVAR} INTERCEP &IVARS;
I = _N_;
DO J = 1 TO &NVAR;
  COVS(I,J) = X(J);
END;
IF EOF THEN OUTPUT COVVEC;
**PROC PRINT;
 
proc means data=inall n sum noprint; var y;
output out=four n=yn sum=ysum;
 
DATA THREE;
 if _n_=1 then set four (keep = yn ysum);
 SET TWO (KEEP = Y &ivars xbeta PRED) end=eof;
array x (&nvar) intercep &ivars;
retain predo;
RETAIN uprod1-uprod&NVAR2 0;                *** outer score products;
RETAIN uo1-uo&NVAR 0;                       *** scores under null hypothesis;
RETAIN uprdo1-uprdo&NVAR2 0;              *** outer null score products;
ARRAY uprod {&NVAR, &NVAR} uprod1-uprod&NVAR2;
ARRAY uo {&NVAR} uo1-uo&NVAR;
ARRAY uprdo {&NVAR, &NVAR} uprdo1-uprdo&NVAR2;
intercep=1;
if xbeta ne .;
if _n_=1 then do;
 predo=ysum/yn;
 do i = 1 to &nvar;
 do j = 1 to &nvar;
  uprod(i,j)=0;
 end;
 end;
end;
resid = (y - pred);
resido = (y - predo);
do i = 1 to &nvar;
  ui=x(i)*resid;
  uio=x(i)*resido;
  uo(i) = uo(i) + uio;
do j = 1 to &nvar;
  uj=x(j)*resid;
  uprod(i,j)=uprod(i,j) + ui*uj;
  ujo=x(j)*resido;
  uprdo(i,j)=uprdo(i,j) + uio*ujo;
end;
end;
 
if eof then output;
***proc print;
***title2
'Elements of the estimated information sandwich and the cov betas';
 
data uprods; set three (keep = uprod1-uprod&nvar2);
data uos; set three (keep = uo1-uo&nvar);
data uprdos; set three (keep = uprdo1-uprdo&nvar2);
data covb; set covvec (keep = covs1-covs&nvar2);
data betas; set betas (keep = beta1-beta&nvar);
 
proc iml;
reset print;
use uprods;
read all into uprod;  **** elements of J matrix;
use uos;
read all into uo;
uo = uo`;
use uprdos;
read all into uprdo;  **** elements of J matrix under H0;
use betas;
read all into beta;
use covb;
read all into cov;
beta=beta`;
var = shape(cov,&nvar,&nvar);
up = shape(uprod,&nvar,&nvar);
covb=var*up*var;                         *** information sandwich;
varb=diag(covb);
seb=sqrt(varb);
j=shape(1,&nvar,1);
za=probit(0.975);
cl=beta` - za*j`*seb;
cu=beta` + za*j`*seb;
clor=exp(cl);
cuor=exp(cu);
Waldj=beta`*inv(varb)*(diag(beta));      *** Wald tests for coefficients;
Waldjp=1-probchi(Waldj,1);
beta›1,1®=0;
ic=I(&nvar-1);
zv=shape(0,(&nvar-1),1);
c=(zv¶¶ic);
c=c`;
Waldmod=(c`*beta)`*inv(c`*covb*c)*(c`*beta);  *** Model robust Wald test;
Waldmodp=1-probchi(Waldmod,(&nvar-1));
 
upo = shape(uprdo,&nvar,&nvar);              *** Robust model score test;
uo›1,1®=0;
score = uo`*inv(upo)*uo;
scorep=1-probchi(score,(&nvar-1));
 
quit;
 
RUN;
 
RUN;
