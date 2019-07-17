*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter7/macklogi.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: uses the Mack et al. data described in Breslow and Day    *;
*          (1980) from a matched case control study. The data set    *;
*          is provided by Norman Breslow on his web site at the      *;
*          University of Washington Department of Biostatistics.     *;
*          The data set downloaded from his site in the fall of      *;
*          1999 is contained in the file                             *;
*          /jml/biostatmethods/datasets/bresmack.dat and the         *;
*          documentation provided by Dr. Breslow is provided in      *;
*          /jml/biostatmethods/datasets/bresmack.text                *;
*          This program computes additional variables used in        *;
*          Problem 7.18.                                             *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
*** this uses the data set from Norm Breslow at U. Washington;
 
filename mackdata '/jml/biostatmethods/datasets/bresmack.dat';
 
***
These are data from a study of the risk of endometrial cancer among
women who resided in a retirement community in Los Angeles.  Each of
the 63 cases was matched to 4 controls within one year of age of the
case, who had the same marital status, and who entered the community
at about the same time.
  The variables in the data set are:
  case: 1 if case, 0 if control, in sets of 5 (1+4)
  age:  in years (a matching variable)
  gbdx: history of gallbladder disease (1=yes, 0=no)
  hypert: history of hypertension, "
  obese: (1=y, 0=n, .=UNKNOWN)
  estrogen: history of estrogen use, "
  conjest: history of conjugated estrogen use, "
  othrest: hostory of non-conjugated estrogen use, "
        NOTE: these above 2 variables are derived from the other
              variables in the next data step below.
  DOSE: Dose of CONJUGATED estrogen use
      0 = 0
      1 = 0.3
      2 = 0.301-0.624
      3 = 0.625
      4 = 0.626-1.249
      5 = 1.25
      6 = 1.26-2.50
      9 = Unknown
  dosecode: dose group of CONJUGATED estrogen use
     0=none (history=no above)
     1=low  (0.1-0.299)
     2=average (0.3-0.625)
     3=high dose (0.626+)
     .=unknown
     Note that these categories are further groupings of Breslows Dose categories.
  est_dur: duration of conjugated estrogen use in months.
      Values > 96 are truncated at 96
  nonest: history of use of non-estrogen drug use (1=y,0=n)
;
 
DATA ONE;
infile mackdata;
INPUT caseset CASE age GBDX HYPERT OBESE ESTROGEN DOSE
 est_dur NONEST;
 
DATA two; SET one;
RETAIN CONTROLN 0;
 
IF CASE=1 THEN DO;
   CONTROLN=0;
END;
 
IF CASE=0 THEN CONTROLN=CONTROLN+1;
 
time = 2 - case;
 
conjest=dose;
 if dose > 1 then conjest=1;
 if dose = 9 then dose=.;
 
dosecode=dose;
  if dose > 3 then dosecode = 3;
 
if obese = . then obese = 0;  **** combine missing with not obese;
 
othrest = estrogen - conjest;
 
if dosecode ge 0 then do;
   dos0=0; dos1=0; dos2=0; dos3=0;
 end;
if dosecode = 0 then dos0=1;
if dosecode = 1 then dos1=1;
if dosecode = 2 then dos2=1;
if dosecode = 3 then dos3=1;
 
title1 'Analysis of Breslow-Day endometrial cancer matched-pairs data';
 
proc univariate freq; var  estrogen conjest dosecode est_dur
     dos0 dos1 dos2 dos3 gbdx hypert obese nonest;
 
RUN;
 
proc phreg data=two;
   model time*case(0) = gbdx hypert obese nonest / ties = discrete;
   strata caseset;
Title2 'Conditional logistic model fit through PHREG';
title3 'model with covariates only';
 
*** insert additional analyses here;
 
run;
 
data four; set two;
if case=1 or controln=1;
***proc print;
title2 'subset of data with 1 case and 1 control';
 
run;
 
data five; set;
retain xc1-xc12;
array xcase (12) xc1-xc12;
array xcont (12) xo1-xo12;
array xin   (12) estrogen conjest dosecode est_dur
     dos0 dos1 dos2 dos3 gbdx hypert obese nonest;
array xout  (12) xest xconjest xdose xest_dur
     xdos0 xdos1 xdos2 xdos3 xgbdx xhypert xobese xnonest;
if case=1 then do;
   do i = 1 to 12;
      xcase(i) = xin(i);
   end;
end;
 
if case=0 then do;
   do i = 1 to 12;
      xcont(i) = xin(i);
      xout(i)  = xcase(i) - xcont(i);
   end;
   y=1;
   output;
end;
 
data six; set five;
if _N_ < 15;
 
proc print; var caseset xest xconjest xdose xest_dur
     xdos0 xdos1 xdos2 xdos3 xgbdx xhypert xobese xnonest;
 
run;
 
proc logistic data=five;
   model y = xgbdx xhypert xobese xnonest / noint;
title2 'conditional logistic regression for 1:1 matching through LOGISTIC';
title3 'model with covariates only';
 
*** insert additional analyses here;
 
run;
 
run;
