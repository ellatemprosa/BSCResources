*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter4/randiter.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Authors: John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Yvonne Sparling                                           *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: The macro %random that performs a fixed-point iterative   *;
*          random effects model analysis of multiple 2x2 tables.     *;
*          The analysis is performed simultaneously for each of the  *;
*          three scales until convergence on all scales is reached.  *;
*          This is a stand-alone macro and does not require that     *;
*          freqanal or freqrand be run beforehand.                   *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
**---**---**---**---**
 
   MACRO: random
   date:   7/29/99
   author: Yvonne Sparling
 
   THE PROGRAM CALCULATES RANDOM EFFECTS MVLE ESTIMATES & THEIR VARIANCES USING AN
   ITERATIVE PROCEDURE.
 
   IF THE OVERDISPERSION PARAMETER (VBRD, VBLOR, VBLRR) IS 0, THEN THE FINAL ESTIMATES ARE THE
   ESTIMATES FROM THE PREVIOUS STEP. IN THIS CASE, THE LAST AND SECOND TO LAST ESTIMATES
   USUALLY WILL NOT BE WITHIN 0.00001 AS THE CONVERGENCE CRITERIA SPECIFIES.
 
**---**---**---**---****---**---**---**---****---**---**---**---****---**---**---**---****---**;
 
%MACRO random;
 
TITLE1 'ITERATIVE RANDOM EFFECTS ANALYSIS OF MULTIPLE 2 X 2 TABLES';
 
DATA ONE; SET ONE;
STRATUM = '        ';STRATUM = K;
 
DATA TWO; SET ONE;
 
TITLE3 'ESTIMATORS AND TESTS'
 ' WITHIN STRATA';
 
A=D1; B=D2; C=N1-D1; D=N2-D2;
N=N1+N2; DT=D1+D2;  P=DT/N;
M=1/(1/N1 + 1/N2);
P1=D1/N1;P2=D2/N2;
RD=P1-P2;
RR=P1/P2; LRR=LOG(RR);
OR=A*D/(B*C); LOR=LOG(OR);
 
VRD = (P1*(1-P1)/N1) + (P2*(1-P2)/N2);
VLRR = ((1-P1)/(N1*P1)) + ((1-P2)/(N2*P2));
VLOR = (1/(N1*P1*(1-P1))) + (1/(N2*P2*(1-P2)));
 
U95RD = RD + 1.96*SQRT(VRD);
L95RD = RD - 1.96*SQRT(VRD);
U95LRR = LRR + 1.96*SQRT(VLRR);
L95LRR = LRR - 1.96*SQRT(VLRR);
U95RR = EXP(U95LRR);
L95RR = EXP(L95LRR);
U95LOR = LOR + 1.96*SQRT(VLOR);
L95LOR = LOR - 1.96*SQRT(VLOR);
U95OR = EXP(U95LOR);
L95OR = EXP(L95LOR);
 
PROC PRINT; VAR STRATUM A B C D N1 N2 P1 P2 P;
 TITLE4
 'FREQUENCIES AND PROPORTIONS WITHIN STRATA';
PROC PRINT; VAR STRATUM
 RD VRD L95RD U95RD;
 TITLE4 'RISK DIFFERENCES AND CONFIDENCE LIMITS';
PROC PRINT; VAR STRATUM
 RR LRR VLRR
 L95LRR U95LRR L95RR U95RR;
 TITLE4 'RELATIVE RISK, LOG RR AND CONFIDENCE LIMITS';
PROC PRINT; VAR STRATUM
 OR LOR VLOR
 L95LOR U95LOR L95OR U95OR;
 TITLE4 'ODDS RATIOS, LOG ODDS AND CONFIDENCE LIMITS';
 
DATA THREE; SET TWO;
TITLE3
 'STRATIFIED-ADJUSTED MVLE ESTIMATORS ';
 
MD = 1/VRD;
MRR = 1/VLRR;
MOR = 1/VLOR;     ** TERMS FOR MVLE ADJUSTED ESTIMATORS;
        MDsq = 1/VRD**2;
        MRRsq = 1/VLRR**2;
        MORsq = 1/VLOR**2;     ** terms for random effects model;
AD = MD*RD;
ALRR = MRR*LRR;
ALOR = MOR*LOR;
run;
 
PROC MEANS SUM noprint;
   VAR MD AD MRR ALRR
       MOR ALOR
       morsq mrrsq mdsq;
   OUTPUT OUT=FOUR
   SUM=TMD TAD TMRR TALRR
       TMOR TALOR
       tmorsq tmrrsq tmdsq;
 
DATA five; SET FOUR
  (keep =TALOR TMOR TALRR TMRR TAD TMD
       tmorsq tmrrsq tmdsq);
 
*** ESTIMATORS AGGREGATED OVER STRATA;
 
   ALOR = TALOR/TMOR;  aor=exp(alor);
   VALOR = 1/TMOR;
   ALRR = TALRR/TMRR;  arr=exp(alrr);
   VALRR=1/TMRR;
   AD = TAD/TMD;
   VAD= 1/TMD;
   U95AD = AD + 1.96*SQRT(VAD);
   L95AD = AD - 1.96*SQRT(VAD);
   U95ALRR = ALRR + 1.96*SQRT(VALRR);
   L95ALRR = ALRR - 1.96*SQRT(VALRR);
   U95ARR = EXP(U95ALRR);
   L95ARR = EXP(L95ALRR);
   U95ALOR = ALOR + 1.96*SQRT(VALOR);
   L95ALOR = ALOR - 1.96*SQRT(VALOR);
   U95AOR = EXP(U95ALOR);
   L95AOR = EXP(L95ALOR);
run;
 
/*
PROC PRINT; VAR AD VAD L95AD U95AD; RUN;
 TITLE4 'RISK DIFFERENCE: ESTIMATE AND CONFIDENCE LIMITS';
PROC PRINT; VAR ARR ALRR VALRR
 L95ALRR U95ALRR L95ARR U95ARR;
 TITLE4 'RELALTIVE RISK: ESTIMATE AND CONFIDENCE LIMITS';
PROC PRINT; VAR AOR ALOR VALOR
 L95ALOR U95ALOR L95AOR U95AOR;
 TITLE4 'ODDS RATIO: ESTIMATE AND CONFIDENCE LIMITS';
run;
*/
 
%LET K=0;
%LET CONVERG=No;
 
data allest (rename=(ad=rad vad=tvad l95ad=l95rad u95ad=u95rad alor=ralor valor=tvalor
             l95alor=l95ralor u95alor=u95ralor aor=raor l95aor=l95raor u95aor=u95raor
             alrr=ralrr valrr=tvalrr l95alrr=l95ralrr u95alrr=u95ralrr arr=rarr
             l95arr=l95rarr u95arr=u95rarr));
  set five (keep=ad vad l95ad u95ad alor valor l95alor u95alor aor l95aor u95aor
                 alrr valrr l95alrr u95alrr arr l95arr u95arr);
  iter=&k;
  label iter='ITERATION';
run;
 
%DO %UNTIL ("&CONVERG"="Yes");
 
%LET K=%EVAL(&K+1);
 
data six; if _n_=1 then set five (keep = alor alrr ad tmor tmrr tmd
       tmorsq tmrrsq tmdsq);
 set three (keep = lor lrr rd mor mrr md);
run;
data six; set six end=eof;
retain chlor chlrr chrd dfh;
 
**** TESTS OF HOMOGENEITY (INTERACTION) USING MVLE;
 
if _n_ = 1 then do;
 chlor=0; chlrr=0; chrd=0; dfh=0;
end;
chlor = chlor + (mor*((lor-alor)**2));
chlrr = chlrr + (mrr*((lrr-alrr)**2));
chrd  = chrd  + (md*((rd-ad)**2));
dfh = dfh+1;
if eof then do;
  dfh = dfh -1;
  pchlor = 1 - probchi(chlor, dfh);
  pchlrr = 1 - probchi(chlrr, dfh);
  pchrd = 1 - probchi(chrd, dfh);
  vblor = tmor*(chlor - dfh) / ((tmor**2) - tmorsq);
  vblrr = tmrr*(chlrr - dfh) / ((tmrr**2) - tmrrsq);
  vbrd = tmd*(chrd - dfh) / ((tmd**2) - tmdsq);
  vblor = max(0, vblor);
  vblrr = max(0, vblrr);
  vbrd = max(0, vbrd);
  if vbrd lt 0.00000001 & vblor lt 0.00000001 & vblrr lt 0.00000001 then call symput('converg','Yes');
  output;
end;
run;
 
%IF "&CONVERG"="No" %THEN %DO;
 
proc print; var chlor pchlor chlrr pchlrr chrd pchrd dfh
  vblor vblrr vbrd;
TITLE3
 "TESTS OF HOMOGENEITY USING MVLE AND BETWEEN VARIANCE COMPONENTS for ITERATION=&k";
 
data seven; if _n_=1 then set six (keep = vblor vblrr vbrd);
 SET TWO;
TITLE3 'random effect MVLE ESTIMATORS';
 
tvrd = vrd + vbrd;
tvlrr = vlrr + vblrr;
tvlor = vlor + vblor;
rmD = 1/TVRD;
rmRR = 1/TVLRR;
rmOR = 1/TVLOR;     ** TErmS FOR rmVLE ADJUSTED ESTImATORS;
        rmDsq = 1/TVRD**2;
        rmRRsq = 1/TVLRR**2;
        rmORsq = 1/TVLOR**2;     ** TErmS FOR RANDOM EFFECTS MODEL;
rAD = rmD*RD;
rALRR = rmRR*LRR;
rALOR = rmOR*LOR;
run;
 
/*
proc print ;
  var stratum tvrd tvlrr tvlor rmd rmrr rmor rad ralrr ralor;
*/
 
PROC MEANS SUM noprint;
   VAR rmD rAD rmRR rALRR
       rmOR rALOR rmDsq rmRRsq rmORsq;
   OUTPUT OUT=FOUR
   SUM= TrmD TrAD TrmRR TrALRR
       TrmOR TrALOR TrmDsq TrmRRsq TrmORsq;
run;
 
 
DATA ESTNEW; SET FOUR
  (keep =TrALOR TrmOR TrALRR TrmRR TrAD TrmD TrmDsq TrmRRsq TrmORsq _FREQ_);
 
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
   iter=&k;
run;
 
data allest;
set allest estnew;
run;
 
proc sort data=estnew; by _freq_; run;
 
%IF &K GT 1 %THEN %DO;
 
  data est;
  merge estold (rename=(rad=radj tvad=tvadj ralrr=ralrrj tvalrr=tvalrrj
                       ralor=ralorj tvalor=tvalorj))
        estnew (rename=(rad=radk tvad=tvadk ralrr=ralrrk tvalrr=tvalrrk
                       ralor=ralork tvalor=tvalork));
        by _freq_;
                rad=abs(radj-radk); ralrr=abs(ralrrj-ralrrk); ralor=abs(ralorj-ralork);
                tvad=abs(tvadj-tvadk); tvalrr=abs(tvalrrj-tvalrrk);
                tvalor=abs(tvalorj-tvalork);
 
                If (RAD le 0.00001 & RALRR le 0.00001 & RALOR le 0.00001 &
                    TVAD le 0.00001 & TVALRR le 0.00001 & TVALOR le 0.00001) then
                        call symput('converg','Yes');
  run;
 
%END;
 
  data estold;
        set estnew;
  run;
 
data five ;
  set ESTold (rename=(rad=ad ralrr=alrr ralor=alor trmd=tmd trmrr=tmrr trmor=tmor
                   trmdsq=tmdsq trmrrsq=tmrrsq trmorsq=tmorsq)) ;
run;
data three ;
  set seven (rename=(rmd=md rmor=mor rmrr=mrr)) ;
run;
 
%END;   * END &CONVERG=NO LOOP;
 
%put CONVERGE = &converg;
 
%IF "&CONVERG"="Yes" %THEN %DO;
 
title3 "Over-dispersion parameter estimates for ITERATION=&k";
PROC PRINT DATA=SIX NOOBS LABEL; VAR VBRD VBLOR VBLRR; RUN;
 
title3 'Random effect aggregate estimators';
PROC PRINT DATA=allEST noobs label; VAR iter rAD tVAD L95rAD U95rAD;
 TITLE4 'RISK DIFFERENCE: ESTIMATE AND CONFIDENCE LIMITS';
run;
PROC PRINT DATA=allEST noobs label; VAR iter rARR rALRR tVALRR
 L95rALRR U95rALRR L95rARR U95rARR;
 TITLE4 'RELALTIVE RISK: ESTIMATE AND CONFIDENCE LIMITS';
PROC PRINT DATA=allEST noobs label; VAR iter rAOR rALOR tVALOR
 L95rALOR U95rALOR L95rAOR U95rAOR;
 TITLE4 'ODDS RATIO: ESTIMATE AND CONFIDENCE LIMITS';
RUN;
 
%END;           * CLOSE CONVERG=YES LOOP;
 
%END;   * END DO UNTIL LOOP;
 
%MEND RANDOM;
 
run;
