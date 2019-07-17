*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter5/pair2x2.sas                  *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: The macro %paired for the analysis of a matched or        *;
*          paired 2x2 table with McNemar's test, the conditional     *;
*          odds ratio and its confidence limits, and the population  *;
*          averaged relative risk and its limits. This is used for   *;
*          Examples 5.7 and 5.8                                      *;
*          This allows for multiple tables (indexed by k) but does   *;
*          not conduct a stratified-adjusted analysis                *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
%MACRO PAIRED;
 
TITLE1 'Paired or Matched 2x2 TABLES';
 
DATA TWO; SET ONE;
 
*** McNemar TEST;
 
CHISQ=(f-g)**2/(f+g);  P = 1 - PROBCHI(CHISQ,1);
 
*** ESTIMATORS;
 
N=e + f + g + h;
P11=e/N; P12=f/N; P21=g/N; P22=h/N;
Nc1=e+g; nr1=e+f;
P1=nr1/N;  P2=nc1/N;
M=f+g;
orc=f/g; lorc=log(orc); vlorc=M/(f*g);
za=probit(0.975);
U95LORC = LORC + 1.96*SQRT(VLORC);
L95LORC = LORC - 1.96*SQRT(VLORC);
U95ORC = EXP(U95LORC);
L95ORC = EXP(L95LORC);
 
RRA=P1/P2; LRRa=LOG(RRa); vlrra=M/(nr1*nc1);
U95lrra = lrra + 1.96*SQRT(Vlrra);
L95lrra = lrra - 1.96*SQRT(Vlrra);
U95rra = EXP(U95lrra);
L95rra = EXP(L95lrra);
 
PROC PRINT; VAR K e f g h Nr1 Nc1 p12 p21 P1 P2 chisq p;
 TITLE4
 'FREQUENCIES AND PROPORCTIONS WITHIN the table and McNemars test';
PROC PRINT; VAR ORC LORC VLORC
 L95LORC U95LORC L95ORC U95ORC;
 TITLE4 'conditional ODDS RATIOS, LOG ODDS AND CONFIDENCE LIMITS';
PROC PRINT; VAR rra Lrra VLrra
 L95Lrra U95Lrra L95rra U95rra;
 TITLE4 'Population averaged RR, log rr AND CONFIDENCE LIMITS';
 
%MEND;
 
DATA ONE;
INPUT K e f g h;
CARDS;
 1 10 31 14 10
*****;
TITLE2 'Example 5.7 -- DCCT pregnancy. *** e and h are dummy numbers, ignore RR';
%PAIRED;
 
run;
