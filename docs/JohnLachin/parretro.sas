*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter5/parretro.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: A macro %retropar to compute the Population Attributable  *;
*          Risk for a case-control 2x2 Table and its confidence      *;
*          limits.                                                   *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
*** PE is the prevalence of exposure to the risk factor in the population;
 
%MACRO retropar;
 
TITLE1 'Population Attributable Risk for a case-control 2 X 2 TABLE';
 
DATA ONE; SET ONE;
 
A=D1; B=D2; C=N1-D1; D=N2-D2;
N=N1+N2; DT=D1+D2;  P=DT/N;
P1=D1/N1;P2=D2/N2;
OR=A*D/(B*C); LOR=LOG(OR);
 
VLOR = (1/(N1*P1*(1-P1))) + (1/(N2*P2*(1-P2)));
U95LOR = LOR + 1.96*SQRT(VLOR);
L95LOR = LOR - 1.96*SQRT(VLOR);
U95OR = EXP(U95LOR);
L95OR = EXP(L95LOR);
 
par=pe*(or-1)/((pe*or)+(1-pe));
lgtpar=log(par/(1-par));
vlgtpar=((or/(or-1))**2)*vlor;
U95LgPAR = Lgtpar + 1.96*SQRT(VLgtpar);
L95LgPAR = Lgtpar - 1.96*SQRT(VLgtpar);
U95PAR = 1/(1+EXP(-U95LgPAR));
L95PAR = 1/(1+EXP(-L95LgPAR));
 
PROC PRINT; VAR A B C D N1 N2 P1 P2 P;
 TITLE4
 'FREQUENCIES AND PROPORTIONS ';
PROC PRINT; VAR
 OR LOR VLOR
 L95LOR U95LOR L95OR U95OR;
 TITLE4 'ODDS RATIOS, LOG ODDS AND CONFIDENCE LIMITS';
PROC PRINT; VAR PE PAR LgtPAR VLgtPAR
 L95LgPAR U95LgPAR L95Par U95PAR;
 TITLE4 'PE PAR, LOGit PAR AND CONFIDENCE LIMITS';
%MEND;
 
data one; input pe d1 n1 d2 n2;
cards;
0.30 154 200 80 200
;
title2 'Example 5.1';
 
%retropar;
 
 
RUN;
 
RUN;
