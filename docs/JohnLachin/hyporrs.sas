*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter8/hyporrs.sas                  *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: The analysis of the DCCT rates of hypoglycemia as in      *;
*          Examples 8.1, 8.2 and 8.3 using the macros %rate and      *;
*          %rate0. The programs for these macros must be submitted   *;
*          before submitting this program, or %include statements    *;
*          must be added. This program also conducts a               *;
*          stratified-adjusted analysis, stratifying by adult        *;
*          versus adolescent, using the macro %adjrate. These        *;
*          results are not shown in text. The analysis uses the      *;
*          data set dccthypo.dat with the path                       *;
*          /jml/biostatmethods/datasets/hypoglycemia/dccthypo.dat.   *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
filename dccthypo '/jml/biostatmethods/datasets/hypoglycemia/dccthypo.dat';
 
***********************************************************;
 
*==============   Initialize Macro Values  ================*;
   *-----------   Variable title  -------------------*;
   %LET TLT='DCCT rates of severe hypoglycemia';
   *--------------   GROUP          VALUES  ------------*;
   %LET GRPTXT1='Standard';  **** CODED AS 0 *******;
   %LET GRPTXT2='Experimental';** CODED AS 1 *******;
   %LET RRTXT='(Ex/Std)';    **** GRP1/GRP0  *******;
   *--------------   Stratification values  ------------*;
   %LET CATTXT1='Adolescents';**** CODED AS 0 *******;
   %LET CATTXT2='Adults';     **** CODED AS 1 *******;
 
 
data allevts;
infile dccthypo;
input grp_e nevents fuday iu duration
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
 
  futime=fuyears;                       *** time in years;
  revents = nevents/fuyears;
 
data one; set;
if _n_ < 50;
proc print; var group nevents fuyears revents insulin duration female
      adult bcval5 hbael hxcoma;
Title1 'DCCT rates of severe hypoglycemia';
Title2 'Subset of observations';
 
%RATE(allevts,events,grp_e,outrr);
%RATE0(allevts,events,grp_e,outrr);
 
%AdjRate(allevts,events,grp_e,adult,outrr);
 
 
RUN;
