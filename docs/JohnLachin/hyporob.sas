*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter8/hyporob.sas                  *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: fits the over-dispersed quasi-likelihood robust Poisson   *;
*          regression models described in Example 8.6, and the       *;
*          robust information sandwich analysis using PROC GENMOD    *;
*          that is described in Example 8.7. Additional computatons  *;
*          using a PROC IML program, as in LogiRob.sas, provide the  *;
*          robust sandwich estimate plus robust inferences.          *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
options ls=80 formdlim='-' ;
 
filename dccthypo '/jml/biostatmethods/datasets/hypoglycemia/dccthypo.dat';
 
****  specify the parameters for this model;
%let nvar=9;   *** # parameters, including intercept;
%let nvar2=81; *** same squared;
%let ivars= grp insulin duration female adult bcval5 hbael hxcoma;  ****model covariate list;
 
data one;
infile dccthypo;
input grp nevents fuday iu duration
 female adult bcval5 hbael hxcoma obweight;
 
  ievents = 0; if nevents > 0 then ievents=1;
  fuyears = fuday/365.25;
  lnyears = log(fuyears);
  insulin = iu/obweight;
  if grp=1 then group='Exp';  *** Exp = Intensive;
  if grp=0 then group='Std';  *** Std = Conventional;
 
d=nevents;  *** these are used later in the additional computations;
t=fuyears;
 
patient = _n_;
 
  Label grp='tx group (1=Int, 0=Conv)'
        group='treatment group: Int or Conv'
        nevents='# severe hypoglycemia episodes'
        ievents='any hypoglycemia episodes (1=Y, 0=N)'
        fuday='total days of follow-up in the study'
        fuyears='total years of follow-up in the study'
        lnyears='log years of follow-up'
        insulin='insulin units per kg weight'
        duration='months diabetes duration'
        female='female (1) or male (0)'
        adult='adult (1) or adolescent (0)'
        bcval5='C-peptide in pmol/mL'
        hbael='level of HbA1c at initial screening'
        hxcoma='Prior history of coma/seizure'
        obweight='body weight in kg.';
 
RUN;
 
PROC GENMOD DATA=one;
  class group;
  model d = group insulin duration female adult bcval5 hbael hxcoma
  / dist = poisson link = log offset = lnyears
    scale=2.606 pscale;
title2 'covariate adjusted treatment group effect';
title3 'Over-dispersed quasi-likelihood using scaled Pearson Chi-square';
RUN;
 
PROC GENMOD DATA=one;
  class group patient;
  model d = group insulin duration female adult bcval5 hbael hxcoma
  / dist = poisson
    link = log
    offset = lnyears;
    repeated subject=patient / type=unstr covb;
title3 'Robust Sandwich Compuations using PROC GENMOD';
RUN;
 
**** SAS code that will save estimates of the mean computed using   ;
**** the parameter estimates and values of the covariate vector     ;
**** in the input data set, together with 95% (by default) CI for   ;
**** the mean. The mean estimates and the CI's are not printed      ;
**** in the default listing file, but are saved in a SAS data set   ;
 
%global _disk_ ; %let _disk_=on   ;
%global _print_; %let _print_=off ;
 
title2 'covariate full model with no over-dispersion adjustments';
title3 'and additional computations of robust inferences';
PROC GENMOD DATA=one;
  MAKE 'OBSTATS' OUT=PRED  ;  ** Data set containing the mean estimates;
  MAKE 'PARMEST' OUT=PARAMS;  ** Data set containing the parameters estimates;
  MAKE 'COV'     OUT=COVB  ;  ** Data set containing the cov of param. est.;
  model d = grp insulin duration female adult bcval5 hbael hxcoma
  / dist = poisson
    link = log
    offset = lnyears
    covb OBSTATS
  ;
PROC CONTENTS DATA=PRED;
PROC print DATA=PARAMS;
PROC print DATA=COVB;
RUN;
 
data two; merge one pred;
resid = (d - pred);
proc univariate; var resid;
RUN;
 
 
DATA BETAS; set params end=eof;
keep BETA1-BETA&NVAR;
ARRAY BETA {&NVAR} BETA1-BETA&NVAR;
retain BETA1-BETA&NVAR;
if _n_ <= &nvar then beta(_n_) = estimate;
if eof then output;
PROC PRINT;
 
DATA COVVEC; SET COVB END=EOF;
KEEP COVS1-COVS&NVAR2;
RETAIN COVS1-COVS&NVAR2;
ARRAY COVS {&NVAR, &NVAR};
ARRAY X {&NVAR} prm1-prm&nvar;
I = _N_;
DO J = 1 TO &NVAR;
  COVS(I,J) = X(J);
END;
IF EOF THEN OUTPUT COVVEC;
PROC PRINT;
 
proc means data=one sum; var d t; output out=four sum=dsum tsum;
 
DATA THREE;
 if _n_=1 then set four (keep = dsum tsum);
 SET TWO (KEEP = d &ivars xbeta PRED) end=eof;
array x (&nvar) intercep &ivars;
retain predo;
RETAIN uprod1-uprod&NVAR2 0;                *** outer score products;
RETAIN uo1-uo&NVAR 0;                       *** scores under null hypothesis;
RETAIN uprodo1-uprodo&NVAR2 0;              *** outer null score products;
ARRAY uprod {&NVAR, &NVAR} uprod1-uprod&NVAR2;
ARRAY uo {&NVAR} uo1-uo&NVAR;
ARRAY uprodo {&NVAR, &NVAR} uprodo1-uprodo&NVAR2;
intercep=1;
if _n_=1 then do;
 predo=dsum/tsum;
 do i = 1 to &nvar;
 do j = 1 to &nvar;
  uprod(i,j)=0;
 end;
 end;
end;
resid = (d - pred);
resido = (d - predo);
do i = 1 to &nvar;
  ui=x(i)*resid;
  uio=x(i)*resido;
  uo(i) = uo(i) + uio;
do j = 1 to &nvar;
  uj=x(j)*resid;
  uprod(i,j)=uprod(i,j) + ui*uj;
  ujo=x(j)*resido;
  uprodo(i,j)=uprodo(i,j) + uio*ujo;
end;
end;
 
if eof then output;
proc print;
title2
'Elements of the estimated information sandwich and the cov betas';
 
data uprods; set three (keep = uprod1-uprod&nvar2);
data uos; set three (keep = uo1-uo&nvar);
data uprodos; set three (keep = uprodo1-uprodo&nvar2);
data covb; set covvec (keep = covs1-covs&nvar2);
data betas; set betas (keep = beta1-beta&nvar);
 
proc iml;
reset print;  *** for test only;
use uprods;
read all into uprod;  **** elements of J matrix;
use uos;
read all into uo;
uo = uo`;
use uprodos;
read all into uprodo;  **** elements of J matrix under H0;
use betas;
read all into beta;
use covb;
read all into cov;
beta=beta`;
var = shape(cov,&nvar,&nvar);
up = shape(uprod,&nvar,&nvar);
covb=var*up*var;                        *** information sandwich;
varb=diag(covb);
seb=sqrt(varb);
j=shape(1,&nvar,1);
za=probit(0.975);
cl=beta` - za*j`*seb;
cu=beta` + za*j`*seb;
clor=exp(cl);
cuor=exp(cu);
Waldj=beta`*inv(varb)*(diag(beta));     *** Wald tests for coefficients;
Waldjp=1-probchi(Waldj,1);
betaÝ1,1¨=0;
ic=I(&nvar-1);
zv=shape(0,(&nvar-1),1);
c=(zv¦¦ic);
c=c`;
Waldmod=(c`*beta)`*inv(c`*covb*c)*(c`*beta);  *** Model robust Wald test;
Waldmodp=1-probchi(Waldmod,(&nvar-1));
 
upo = shape(uprodo,&nvar,&nvar);              *** Robust model score test;
uoÝ1,1¨=0;
score = uo`*inv(upo)*uo;
scorep=1-probchi(score,(&nvar-1));
 
quit;
 
RUN;
 
RUN;
