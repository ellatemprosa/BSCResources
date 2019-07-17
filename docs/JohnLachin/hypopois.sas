*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter8/hypopois.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: fits the Poisson regression models using the SAS program  *;
*          shown in Table 8.2 that generates the output shown in     *;
*          tables 8.3, 8.4 and 8.5. The output generated differs     *;
*          slightly from that shown in the tables. This program      *;
*          uses the class variable for treatment group defined       *;
*          using the categories "exp" (experimental) versus "std"    *;
*          (standard). Since "std" has the higher alphabetical order,*;
*          this group is used as the reference category for the      *;
*          model effects shown in Tables 8.4 and 8.5. These group    *;
*          labels were later changed to Intensive and conventional,  *;
*          respectively. Thus, I edited the program output and       *;
*          manually changed "Exp" to "Int" and "Std" to "Conv".      *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
*** Note that the data set dccthypo is used not the data as shown in Table 8.1;
 
options ls=80 nodate;
 
filename dccthypo '/jml/biostatmethods/datasets/hypoglycemia/dccthypo.dat';
 
data allevts;
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
 
  Label grp='tx group (1=Int, 0=Conv)'
        group='treatment group: Int or Conv'
        nevents='# severe hypoglycemia episodes'
        ievents='any hypoglycemia episodes (1=Y, 0=N)'
        fuday='total days of follow-up in the study'
        fuyears='total years of follow-up in the study'
        lnyears='log years of follow-up'
        rate='rate of episodes per year of follow-up'
        insulin='insulin units per kg weight'
        duration='months diabetes duration'
        female='female (1) or male (0)'
        adult='adult (1) or adolescent (0)'
        bcval5='C-peptide in pmol/mL'
        hbael='level of HbA1c at initial screening'
        hxcoma='Prior history of coma/seizure'
        obweight='body weight in kg.';
 
data one; set allevts;
 
proc genmod;
  model nevents =
  / dist = poisson
    link = log
    offset = lnyears
    ;
TITLE1 'Poisson regression models of risk of hypoglycemia';
title2 'null model';
 
proc genmod;
  class group;
  model nevents = group
  / dist = poisson
    link = log
    offset = lnyears
    ;
title2 'unadjusted treatment group effect';
 
proc genmod; class group;
  model nevents = group
      insulin duration female
      adult bcval5 hbael hxcoma
  / dist = poisson
    link = log
    offset = lnyears
    covb
    ;
title2 'covariate adjusted treatment group effect';
 
RUN;
 
**** SAS code that will save estimates of the mean computed using   ;
**** the parameter estimates and values of the covariate vector     ;
**** in the input data set, together with 95% (by default) CI for   ;
**** the mean. The mean estimates and the CI's are not printed      ;
**** in the default listing file, but are saved in a SAS data set   ;
 
%global _disk_ ; %let _disk_=on   ;
%global _print_; %let _print_=off ;
 
PROC GENMOD DATA=one;
  MAKE 'OBSTATS' OUT=PRED  ;  ** Data set containing the mean estimates;
  class group;
  model nevents = group insulin duration female
      adult bcval5 hbael hxcoma
  / dist = poisson link = log offset = lnyears covb
        OBSTATS
  ;
RUN;
 
**proc print data=pred;
 
data two; merge one pred;
if _n_ < 50;
proc print;
 
RUN;
PROC CONTENTS DATA=PRED;
RUN;
 
proc means data=one sum noprint;
 var d t; output out=four sum=dsum tsum;
 
data two;
 if _n_=1 then set four (keep = dsum tsum);
 merge one pred;
resid = (d - pred);
pred0=(dsum/tsum)*t;
resid0 = (d - pred0);
regresid=(pred-pred0);
 
**proc univariate;
** var resid resido regresid;
 
data three; set two end=eof;
retain ssr ssr0 ssregr 0;
retain lnl0 lnl0p lnl lnlp lnlf lnlfp 0;
ssr=ssr + resid**2;
ssr0=ssr0 + resid0**2;
ssregr=ssregr + regresid**2;
lnl = lnl - pred + (d*log(pred));                    *** likelihood without constant;
lnlp = lnlp -pred + (d*log(pred)) - lgamma(d+1);     *** complete likelihood with constant;
lnl0 = lnl0 -pred0 + (d*log(pred0));                 *** null model likelihood w/o constant;
lnl0p = lnl0p -pred0 + (d*log(pred0)) - lgamma(d+1); *** null model L with constant;
dld=0; if d > 0 then dld = d*log(d);
lnlf = lnlf - d + dld;                               *** saturated model L without constant;
lnlfp = lnlfp -d + dld - lgamma(d+1);                *** saturated model L with constant;
Nobs=_n_;
if eof then output;
 
data three; set three;
 
R2resid = ssregr/ssr0;
**devmp = (-2*lnlp) - (-2*lnlfp); *** same result as next line;
devm = (-2*lnl) - (-2*lnlf);
**dev0p = (-2*lnl0p) - (-2*lnlfp);
dev0 = (-2*lnl0) - (-2*lnlf);
 
X2LR = dev0 - devm;
pvalue = 1 - probchi(x2lr, 8);   ************ specify the df;
R2RL = 1 - exp(-x2lr/nobs);
 
R2logL = (lnl0p - lnlp)/lnl0p;
R2Dev = (dev0 - devm)/dev0;
 
proc print;
 
proc print; var ssr ssr0 ssregr R2resid;
title3 'Sum of squares of model residuals, null model rediduals, regression and R2 residual';
 
proc print; var lnl0p lnlp R2logl;
title3 'Complete log likelihoods (including constant) and R2 log likelihood';
 
proc print; var dev0 devm r2dev;
title3 'null and full model deviances, R2 deviance explained';
 
proc print; var x2lr pvalue r2rl;
title3 'model likelihood test and p-value, Maddalas R2 LR';
 
RUN;
