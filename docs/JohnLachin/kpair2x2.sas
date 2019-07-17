*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter5/kpair2x2.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: The macro %Kpaired for the stratified analysis of K       *;
*          paired or matched 2x2 tables. This is used in Examples    *;
*          5.12 and 5.13.                                            *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
%MACRO KPAIRED;
 
TITLE2 'Stratified-Adjusted Analysis of K matched or paired 2x2 tables';
 
DATA TWO; SET ONE;
 
*** TESTS WITHIN STRATA;
 
CHISQ=(f-g)**2/(f+g);  P = 1 - PROBCHI(CHISQ,1);
 
*** ESTIMATORS WITHIN STRATA;
 
N=e+f+g+h;
P11=e/N; P12=f/N; P21=g/N; P22=h/N;
nr1=e+f; nc1=e+g;
P1=nr1/N;  P2=nc1/N;
M=f+g;
za=probit(0.975);
 
**    LOG ODDS RATIO SCALE;
 
or=f/g; lor=log(or); vlor=m/(f*g);
U95Lor = Lor + 1.96*SQRT(VLor);
L95Lor = Lor - 1.96*SQRT(VLor);
U95or = EXP(U95Lor);
L95or = EXP(L95Lor);
      MOR = 1/Vlor;        ** WEIGHTS FOR AGGREGATE ESTIMATOR;
      ALOR = MOR*Lor;      ** CONTRIBUTION TO ESTIMATOR;
        MORsq = MOR**2;     ** terms for random effects model;
 
**    LOG RELATIVE RISK SCALE;
 
rr=P1/P2; Lrr=LOG(rr); vlrr=M/(nr1*nc1);
U95lrr = lrr + 1.96*SQRT(Vlrr);
L95lrr = lrr - 1.96*SQRT(Vlrr);
U95rr = EXP(U95lrr);
L95rr = EXP(L95lrr);
      MRR = 1/VlRR;        ** WEIGHTS FOR AGGREGATE ESTIMATOR;
      ALRR = MRR*LRR;      ** CONTRIBUTION TO ESTIMATOR;
        MRRsq = MRR**2;
 
PROC PRINT; VAR K e f g h Nr1 Nc1 p12 p21 P1 P2 chisq p;
 TITLE4
 'FREQUENCIES AND PROPorTIONS WITHIN the table and McNemars test';
PROC PRINT; VAR or Lor VLor
 L95Lor U95Lor L95or U95or;
 TITLE4 'conditional ODDS RATIOS, LOG ODDS AND CONFIDENCE LIMITS';
PROC PRINT; VAR rr Lrr VLrr
 L95Lrr U95Lrr L95rr U95rr;
 TITLE4 'Population averaged RR, log rr AND CONFIDENCE LIMITS';
 
PROC PRINT; var k lor vlor mor alor lrr vlrr mrr alrr;
  title4 'contributions to the MVLE';
 
*** AGGREGATED OVER STRATA;
 
PROC MEANS SUM noprint;
   VAR MRR ALRR
       MOR ALOR
       morsq mrrsq;
   OUTPUT OUT=three
   SUM=TMRR TALRR
       TMOR TALOR
       tmorsq tmrrsq;
 
DATA four; SET THREE;
 
*** ESTIMATORS AGGREGATED OVER STRATA;
 
   ALOR = TALOR/TMOR;
   VALOR = 1/TMOR;
aor=exp(alor);
U95Lor = ALOR + 1.96*SQRT(VALOR);
L95Lor = ALOR - 1.96*SQRT(VALOR);
U95or = EXP(U95Lor);
L95or = EXP(L95Lor);
 
   ALRR = TALRR/TMRR;
   VALRR=1/TMRR;
arr=exp(alrr);
U95lrr = Alrr + 1.96*SQRT(Valrr);
L95lrr = Alrr - 1.96*SQRT(Valrr);
U95rr = EXP(U95lrr);
L95rr = EXP(L95lrr);
 
title3 'MVLE analysis over strata';
PROC PRINT; VAR Aor ALOR VALOR
 L95Lor U95Lor L95or U95or;
 TITLE4 'conditional ODDS RATIOS, LOG ODDS AND CONFIDENCE LIMITS';
PROC PRINT; VAR arr aLrr VaLrr
 L95Lrr U95Lrr L95rr U95rr;
 TITLE4 'Population averaged RR, log rr AND CONFIDENCE LIMITS';
 
*** ESTIMATORS AGGREGATED OVER STRATA;
 
data six; if _n_=1 then set four (keep = alor alrr tmor tmrr tmorsq tmrrsq);
 set two (keep = lor lrr mor mrr);
data six; set six end=eof;
retain chlor chlrr dfh;
 
**** TESTS OF HOMOGENEITY (INTERACTION) USING MVLE;
 
if _n_ = 1 then do;
 chlor=0; chlrr=0; dfh=0;
end;
chlor = chlor + (mor*((lor-alor)**2));
chlrr = chlrr + (mrr*((lrr-alrr)**2));
dfh = dfh+1;
if eof then do;
  dfh = dfh -1;
  pchlor = 1 - probchi(chlor, dfh);
  pchlrr = 1 - probchi(chlrr, dfh);
  vblor = tmor*(chlor - dfh) / ((tmor**2) - tmorsq);
  vblrr = tmrr*(chlrr - dfh) / ((tmrr**2) - tmrrsq);
  vblor = max(0, vblor);
  vblrr = max(0, vblrr);
  output;
end;
proc print; var chlor pchlor chlrr pchlrr dfh
  vblor vblrr;
TITLE3
 'TESTS OF HOMOGENEITY USING MVLE, AND THE BETWEEN VARIANCE COMPONENTS';
 
data two; set one;
keep i k w;
i=12; w=f; output;
i=21; w=g; output;
proc freq; table i*k / all;
weight w;
title3 'contingency test of independence from 2*K table';
 
RUN;
%MEND;
 
 
DATA ONE;
INPUT K e f g h;
CARDS;
 1 10  6  5 11
 2 10 10  5 11
 3 10  6  1 11
 4 10  9  3 11
*****;
TITLE2 'DCCT pregnancy *** e and h are dummy numbers, ignore RR';
%KPAIRED;
 
DATA ONE;
INPUT K e f g h;
CARDS;
 1 10  1 3 11
 2 10  4 3 11
 3 10  4 5 11
;
 
title1 'H-L App 5 with 1 control, stratified by chk for case/control';
TITLE2 ' *** e and h are dummy numbers, ignore RR';
%KPAIRED;
 
 
DATA ONE;
INPUT K e f g h;
CARDS;
 1 10  20 21 11
 2 10  2   3 11
 3 10  11 14 11
 4 10  2   3 11
;
 
title1 'Categorical Analysis p 271, flu by vaccination stratified by lung disease';
TITLE2 '*** e and h are dummy numbers, ignore RR';
 
%KPAIRED;
 
run;
