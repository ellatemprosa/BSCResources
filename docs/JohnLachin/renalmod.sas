*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter7/renalmod.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: fits the logistic regression models with interactions     *;
*          for the renal data as shown in Example 7.12. The job      *;
*          renal.sas must be run to set-up the renal data before     *;
*          running this job.                                         *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
options formdlim='-';
 
data one; set renal;
i1=int*yearsdm;
i2=hbael*sbp;
n=1;
 
proc logistic descending;
model micro24 = int hbael yearsdm sbp female/ RL;
 
proc genmod data =one;
class int female;
model micro24/n = int hbael yearsdm sbp female
  / dist=binomial link=logit type3;
 
proc logistic descending;
model micro24 = int hbael yearsdm sbp female i1 i2 /covb RL;
 
proc genmod data =one;
class int female;
model micro24/n = int hbael yearsdm sbp female int*yearsdm hbael*sbp
  / dist=binomial link=logit type3;
 
run;
