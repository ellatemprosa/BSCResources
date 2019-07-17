*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter8/randrate.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Authors: Yvonne Sparling                                           *;
*          John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: is a supplemental program that uses a fixed point         *;
*          algorithm to compute convergent estimates of the          *;
*          over-dispersion parameters, the mean rates and their      *;
*          variances, both under the alternative hypothesis and      *;
*          also under the null hypothesis. The latter could be used  *;
*          as the basis for a Z-test of the difference between two   *;
*          groups. These computations are not described in the       *;
*          text, but follow from those described in Section 4.10.2.  *;
*          The program is applied to the DCCT data used in Example   *;
*          8.3.                                                      *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
*** Author: Yvonne Sparling;
 
*** Date: July 26, 1999;
 
 ******----- * ----- * ----- * ----- * ----- * ----- * ----- * ----- * ----- * ----- *********
   ITERATIVE PROCEDURE TO ESTIMATE LAMBDA & ITS VARIANCE FOR A POISSON MODEL
   USING THE MOMENT GENERATING ESTIMATOR FROM EXPRESSION 8.37 OF BIOSTATISTICAL METHODS BOOK.
   JULY 16, 1999
 ******----- * ----- * ----- * ----- * ----- * ----- * ----- * ----- * ----- * ----- *********;
 
options ps=60 ls=80 nodate nonumber;
 
filename dccthypo '/jml/biostatmethods/DataSets/hypoglycemia/dccthypo.dat';
 
data allevts;
infile dccthypo;
input patient grp_e nevents fuday iu duration
 female adult bcval5 hbael hxcoma obweight;
 
  ievents = 0; if nevents > 0 then ievents=1;
  fuyears = fuday/365.25;
  rate_1cy =  100*nevents/fuyears;
  revents = nevents/fuyears;
  Label hxcoma='Prior history of coma/seizure'
        fuyears='Total Years of Follow-up';
  lnyears = log(fuyears);
  insulin = iu/obweight;
  if grp_e=1 then group='experimental';
  if grp_e=0 then group='standard';
 
run;
 
        ******************************************
        * ESTIMATES UNDER ALTERNATIVE HYPOTHESIS *
        ******************************************;
 
 
data allevts (keep=d t r grp_e);
  set allevts (rename=(nevents=d fuyears=t revents=r));
run;
 
 
proc sort data=allevts; by grp_e;
proc summary data=allevts;
  by grp_e;
  var d t;
  output out=out0 sum=d0 t0;
run;
 
proc sort data=out0; by grp_e;
proc sort data=allevts; by grp_e; run;
data temp0;
  merge allevts out0; by grp_e;
        r0  = D0/T0;                    * Lambda estimate ;
        num = T*((r-r0)**2) - r0;       * Numerator for moment estimate for dispersion parameter;
run;
 
proc sort data=temp0; by grp_e;
proc summary data=temp0;
  by grp_e;
  var d t num ;
  output out=out0 sum=d0 t0 num0;
run;
 
proc sort data=out0; by grp_e; run;
data temp0;
  merge allevts out0; by grp_e;
        r0 = D0/T0;
       num = T*((r-r0)**2) - r0;
        s0 = num0 / t0;                 * Moment Estimate for dispersion parameter;
      num2 = (t**2)*s0 + r0*t;          * Numerator for estimate of variance of lambda estimate ;
run;
 
proc sort data=temp0; by grp_e;
proc summary data=temp0;
  by grp_e;
  var d t num num2;
  output out=out0 sum=d0 t0 num0 num20 ;
run;
 
data est0 (keep= r0 s0 varr0 varlnr0 j grp_e) ;
  set out0;
      r0 = d0 / t0;
      s0 = num0 / t0;
   varr0 = num20 / (t0**2);             * Estimated variance of estimated lambda ;
 varlnr0 = varr0 / (r0**2);             * Estimated variance for log scale ;
       j = 0;
  label r0='Lambda Est' s0='Dispersion Parameter Est'
        varr0='Est Var of Lambda Est' varlnr0='Est Var of Ln(Lambda)'
        grp_e='Group';
run;
 
title 'ESTIMATES FOR THE Jth STEP OF ITERATIVE METHOD';
title2 'UNDER THE ALTERNATIVE HYPOTHESIS';
proc print data=est0 noobs label;
var j grp_e r0 s0 varr0 varlnr0;
run;
 
 
%macro pois;
 
%LET CURR=0;
%LET PREV=-1;
%LET CONVERG=No;
 
%DO %UNTIL("&CONVERG"="Yes");
 
%LET CURR=%EVAL(&CURR+1);
%LET PREV=%EVAL(&PREV+1);
 
proc sort data=est&prev; by grp_e;
proc sort data=allevts; by grp_e; run;
 
data temp&curr;
  merge allevts  (keep = t d r grp_e)
        est&prev (keep = r&prev s&prev varr&prev varlnr&prev grp_e);
        by grp_e;
                vri = (r&prev/t) + s&prev;    * Var(ri);
                dw  = r/vri;                  * Weighted no. of events;
                tw  = 1/vri;                  * Weight;
run;
 
proc sort data=temp&curr; by grp_e;
proc summary data=temp&curr;
  by grp_e;
  var dw tw;
  output out=out1 sum=dw0 tw0;
run;
 
proc sort data=out1; by grp_e; run;
data temp&curr;
  merge temp&curr out1; by grp_e;
       r&curr = Dw0/Tw0;                    * Lambda estimate;
       num    = T*((r-r&curr)**2) - r&curr; * Numerator for moment estimate of dispersion param;
run;
 
proc sort data=temp&curr; by grp_e;
proc summary data=temp&curr;
  by grp_e;
  var dw tw num t;
  output out=out1 sum=dw0 tw0 num0 t0;
run;
 
proc sort data=out1; by grp_e; run;
data temp&curr;
  merge temp&curr out1; by grp_e;
      r&curr = Dw0 / Tw0;
      num    = T*((r-r&curr)**2) - r&curr;
      s&curr = num0 / t0;                 * Moment estimate for dispersion param ;
      num2   = (t**2)*s&curr + r&curr*t;  * Numerator for estimate of variance of lambda estimate;
run;
 
proc sort data=temp&curr; by grp_e;
proc summary data=temp&curr;
  by grp_e;
  var dw tw num num2 t;
  output out=out1 sum=dw0 tw0 num0 num20 t0;
run;
 
data est&curr (keep= r&curr s&curr varr&curr varlnr&curr j grp_e) ;
  set out1;
      r&curr = dw0 / tw0;
      s&curr = num0 / t0;
   varr&curr = num20 / (t0**2);             * Estimated variance of estimated lambda ;
 varlnr&curr = varr&curr / (r&curr**2);     * Estimated variance for log scale ;
           j = &curr;
  label r&curr='Lambda Est' s&curr='Dispersion Param Est'
        varr&curr='Est Var of Lambda Est' varlnr&curr='Est Var of Ln(Lambda)'
        grp_e='Group';
 
run;
 
proc sort data=est&curr; by grp_e;
proc sort data=est&prev; by grp_e; run;
 
%put &converg=;
 
data converg;
  merge est&curr est&prev (drop=j);  by grp_e;
        rdiff = abs(r&curr-r&prev);
        sdiff = abs(s&curr-s&prev);
        if (rdiff le 0.00001) & (sdiff le 0.00001) then call symput('converg','Yes');
        else call symput('converg','No');
run;
 
%IF "&CONVERG"="Yes" %THEN %DO;
proc print data=converg noobs label;
var j grp_e r&curr s&curr varr&curr varlnr&curr;
run;
%END;
 
%END;
 
%mend pois;
 
%pois;
 
 
 
 
 
        ***********************************
        * ESTIMATES UNDER NULL HYPOTHESIS *
        ***********************************;
 
data allevts;
  set allevts;
  sorter=1;
run;
 
proc summary data=allevts;
  var d t;
  output out=out0 sum=d0 t0;
run;
 
data out0; set out0; sorter=1; run;
 
proc sort data=out0; by sorter;
proc sort data=allevts; by sorter; run;
data temp0;
  merge allevts out0; by sorter;
        r0  = D0/T0;                    * Lambda estimate ;
        num = T*((r-r0)**2) - r0;       * Numerator for moment estimate for dispersion parameter;
run;
 
proc sort data=temp0; by grp_e;
proc summary data=temp0;
  by grp_e;
  var num t;
  output out=out0 sum=num0 den0;
run;
 
proc sort data=out0; by grp_e; run;
data temp0;
  merge out0 (keep=grp_e num0 den0)
        temp0; by grp_e;
        s0 = num0 / den0;                * Moment Estimate for dispersion parameter;
        num2 = (t**2)*s0 + r0*t;         * Numerator for estimate of variance of lambda estimate ;
run;
 
proc sort data=temp0; by grp_e;
proc summary data=temp0;
  by grp_e;
  var r0 s0 num2 t;
  output out=out0 sum=junk1 junk2 num20 t0 mean=r0 s0;
run;
 
data est0 (keep= r0 s0 varr0 varlnr0 j grp_e) ;
  set out0 (drop=junk1 junk2);
        varr0 = num20 / (t0**2);         * Estimated variance of estimated lambda ;
        varlnr0 = varr0 / (r0**2);       * Estimated variance for log scale ;
        j = 0;
        label r0='Lambda Est' s0='Dispersion Parameter Est'
              varr0='Est Var of Lambda Est' varlnr0='Est Var of Ln(Lambda)'
              grp_e='Group';
run;
 
title 'ESTIMATES FOR THE Jth STEP OF ITERATIVE METHOD';
title2 'UNDER THE NULL HYPOTHESIS';
proc print data=est0 noobs label;
var j grp_e r0 s0 varr0 varlnr0;
run;
 
 
%macro poisnull(curr,prev);
 
%LET CURR=0;
%LET PREV=-1;
%LET CONVERG=No;
 
%DO %UNTIL("&CONVERG"="Yes");
 
%LET CURR=%EVAL(&CURR+1);
%LET PREV=%EVAL(&PREV+1);
 
proc sort data=est&prev; by grp_e;
proc sort data=allevts; by grp_e; run;
 
data temp&curr;
  merge allevts  (keep = t d r grp_e sorter)
        est&prev (keep = r&prev s&prev varr&prev varlnr&prev grp_e);
        by grp_e;
                vri = (r&prev/t) + s&prev;    * Var(ri);
                dw  = r/vri;                  * Weighted no. of events;
                tw  = 1/vri;                  * Weight;
run;
 
proc summary data=temp&curr;
  var dw tw;
  output out=out1 sum=dw0 tw0;
run;
data out1; set out1; sorter=1; run;
 
proc sort data=temp&curr; by sorter;
proc sort data=out1; by sorter; run;
data temp&curr;
  merge temp&curr out1; by sorter;
       r&curr = Dw0/Tw0;                    * Lambda estimate;
       num    = T*((r-r&curr)**2) - r&curr; * Numerator for moment estimate of dispersion param;
run;
 
proc sort data=temp&curr; by grp_e; run;
proc summary data=temp&curr;
  by grp_e;
  var num t;
  output out=out1 sum=num0 t0;
run;
 
proc sort data=out1; by grp_e; run;
data temp&curr;
  merge temp&curr out1; by grp_e;
      s&curr = num0 / t0;                 * Moment estimate for dispersion param ;
      num2   = (t**2)*s&curr + r&curr*t;  * Numerator for estimate of variance of lambda estimate;
run;
 
proc sort data=temp&curr; by grp_e;
proc summary data=temp&curr;
  by grp_e;
  var s&curr r&curr num2 t;
  output out=out1 sum=junk1 junk2 num20 t0 mean=s&curr r&curr;
run;
 
data est&curr (keep= r&curr s&curr varr&curr varlnr&curr j grp_e) ;
  set out1;
   varr&curr = num20 / (t0**2);             * Estimated variance of estimated lambda ;
 varlnr&curr = varr&curr / (r&curr**2);     * Estimated variance for log scale ;
           j = &curr;
  label r&curr='Lambda Est' s&curr='Dispersion Param Est'
        varr&curr='Est Var of Lambda Est' varlnr&curr='Est Var of Ln(Lambda)'
        grp_e='Group';
run;
 
proc sort data=est&curr; by grp_e;
proc sort data=est&prev; by grp_e; run;
 
data converg;
  merge est&curr est&prev (drop=j);  by grp_e;
        sdiff = abs(s&curr-s&prev);
        rdiff = abs(r&curr-r&prev);
        if (rdiff le 0.00001) & (sdiff le 0.00001) then call symput('converg','Yes');
        else call symput('converg','No');
run;
 
%IF "&CONVERG"="Yes" %THEN %DO;
proc print data=converg noobs label;
var j grp_e r&curr s&curr varr&curr varlnr&curr;
run;
%END;
 
%END;
 
%mend;
 
%poisnull;
