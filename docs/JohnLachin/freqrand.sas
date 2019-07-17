*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter4/freqrand.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: This is the principal macro %freqrand that computes the   *;
*          random effects model one-step adjusted estimates for      *;
*          multiple 2x2 tables. The macro starts with the            *;
*          computation of the variance component estimate for each   *;
*          scale, followed by the random effects model               *;
*          computations. The macro %freqanal must be invoked before  *;
*          this macro. See Ktab2x2.sas for an example.               *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
***** MACRO freqrand;  **** MACRO %freqanal must be called before this is called;
 
*** This macro computes the random effects model one-step adjusted estimates for
multiple 2x2 tables. The macro starts with the computation of the Cochran test
of homogeneity on each scale, with the variance component estimates, followed by
the random effects model computations.
**;
 
%MACRO freqrand;
 
TITLE2 'Random Effects ANALYSIS OF MULTIPLE 2 X 2 TABLES';
 
data seven; set seven;
  vblor = tmor*(chlor - dfh) / ((tmor**2) - tmorsq);
  vblrr = tmrr*(chlrr - dfh) / ((tmrr**2) - tmrrsq);
  vbrd = tmd*(chrd - dfh) / ((tmd**2) - tmdsq);
  vblor = max(0, vblor);
  vblrr = max(0, vblrr);
  vbrd = max(0, vbrd);
proc print; var vblor vblrr vbrd;
TITLE3 'Between stratum variance components';
 
data eight; if _n_=1 then set seven (keep = vblor vblrr vbrd);
 SET TWO;
TITLE3
 'random effect MVLE ESTIMATORS';
 
if stratum NE 'marginal';
tvrd = vrd + vbrd;
tvlrr = vlrr + vblrr;
tvlor = vlor + vblor;
rmD = 1/TVRD;
rmRR = 1/TVLRR;
rmOR = 1/TVLOR;     ** TErmS FOR rmVLE ADJUSTED ESTImATORS;
rAD = rmD*RD;
rALRR = rmRR*LRR;
rALOR = rmOR*LOR;
 
proc print;
 var stratum tvrd tvlrr tvlor rmd rmrr rmor rad ralrr ralor;
 
PROC MEANS SUM noprint;
   VAR rmD rAD rmRR rALRR
       rmOR rALOR;
   OUTPUT OUT=nine
   SUM= TrmD TrAD TrmRR TrALRR
       TrmOR TrALOR;
 
DATA nine; SET nine
  (keep =TrALOR TrmOR TrALRR TrmRR TrAD TrmD);
 
*** ESTIMATORS AGGREGATED OVER STRATA;
 
   rALOR = TrALOR/TrmOR;  raor = exp(ralor);
   TVALOR = 1/TrmOR;
   rALRR = TrALRR/TrmRR;  rarr = exp(ralrr);
   TVALRR=1/TrmRR;
   rAD = TrAD/TrmD;
   TVAD= 1/TrmD;
 
   U95rAD = rAD + 1.96*SQRT(tVAD);
   L95rAD = rAD - 1.96*SQRT(tVAD);
   U95rALRR = rALRR + 1.96*SQRT(tVALRR);
   L95rALRR = rALRR - 1.96*SQRT(tVALRR);
   U95rARR = EXP(U95rALRR);
   L95rARR = EXP(L95rALRR);
   U95rALOR = rALOR + 1.96*SQRT(tVALOR);
   L95rALOR = rALOR - 1.96*SQRT(tVALOR);
   U95rAOR = EXP(U95rALOR);
   L95rAOR = EXP(L95rALOR);
 
title3 'Random effect aggregate estimators';
PROC PRINT; VAR rAD tVAD L95rAD U95rAD;
 TITLE4 'RISK DIFFERENCE: ESTIMATE AND CONFIDENCE LIMITS';
PROC PRINT; VAR rARR rALRR tVALRR
 L95rALRR U95rALRR L95rARR U95rARR;
 TITLE4 'RELALTIVE RISK: ESTIMATE AND CONFIDENCE LIMITS';
PROC PRINT; VAR rAOR rALOR tVALOR
 L95rALOR U95rALOR L95rAOR U95rAOR;
 TITLE4 'ODDS RATIO: ESTIMATE AND CONFIDENCE LIMITS';
 
RUN;
 
%MEND;
 
RUN;
 
RUN;
