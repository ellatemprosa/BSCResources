*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter3/ppwrplot.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: Plot the power function for the test for two proportions  *;
*          for different sample sizes as a function of the risk      *;
*          difference, relative risk or odds ratios for different    *;
*          sample sizes.                                             *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
***
 **    POWER OF ANALYSES OF PROPORTIONS FOR 2 INDEPENDENT GROUPS
 **
 **    PARAMETERS ARE:
 **       TA = ALPHA CRITICAL VALUE, 1.645 FOR ALPHA 0,05 ONE-SIDED
 **       N  = TOTAL SAMPLE SIZE
 **       QC = SAMPLE FRACTION FOR THE C GROUP, NC = QC*N
 **       QE = SAMPLE FRACTION FOR THE E GROUP, NE = QE*N
 **       PC = PROBABILITY FOR GROUP 1 ( << 0.5)
 **       DD = DIFFERENCE DD = PE - PC
 **       RR = RELATIVE RISK RR = PE/PC
 **       OR = ODDS RATIO OR = PE(1-PC) / PC(1-PE)
 **
 ;
*****  specify output filename;
     filename out '/jml/biostatmethods/chapter3/plot.out';
 %MACRO POWR;
       NE=N*QE;
       NC=N*QC;
       PZTL = 0.0;
       S0  = PH*(1-PH) * ((1/NE)+(1/NC));
       S1E = PE*(1-PE) / NE;
       S1C = PC*(1-PC) / NC;
       NUMA= TA*SQRT(S0);
       NUM= SQRT( S1E + S1C);
       PZT =  (DD - NUMA)/NUM;
       PZTL=PROBNORM(PZT );
 %MEND POWR;
 
 
 DATA ONE;
 ***** specify parameters;
       PC=0.2; TA=1.96;
       QC=0.5; QE=0.5;
 TITLE1 'POWER CURVES FOR THE TEST OF A DIFFERENCE IN PROPORTIONS';
 TITLE2 ' '  ;
 TITLE3 'PC = 0.20,   QC = 0.5,  QE = 0.5,   ZA = 1.96';
 TITLE4 ' '  ;
 TITLE5 'A : N = 200,        B : N = 400,         C : N = 800';
 *****************************************************************;
 C=2.0;
 C=SQRT(C);
 MK=20;                        *** number of points;
   MK1=MK-1;
 
 DO K=1 TO MK1;
 DD=(1-PC) * 0.5 * K/MK;       *** dd ranges from 0+ to 0.5*(1-PC);
 PE=PC + DD;
 PH = (QE*PE) + (QC*PC);
 RR = PE/PC;
 OR = PE*(1-PC) / (PC*(1-PE));
 
 ***** specify parameters;
 N=200;
 %POWR;
 P1=PZTL;
 
 ***** specify parameters;
 N=400;
 %POWR;
 P2=PZTL;
 
 ***** specify parameters;
 N=800;
 %POWR;
 P3=PZTL;
 
 DIFF=DD;
 FILE OUT;
   PUT PC TA QE QC DD RR OR PE P1 P2 P3;
 OUTPUT;
 END;
 PROC PRINT;
 PROC PLOT NOLEGEND;
      PLOT P1*DIFF = 'A'
           P2*DIFF = 'B'
           P3*DIFF = 'C'
    / OVERLAY
      VAXIS = 0 TO 1.0 BY .10
      VZERO HZERO
   ;
TITLE6 'IN TERMS OF THE DIFFERENCE IN PROPORTIONS';
 PROC PLOT NOLEGEND;
      PLOT P1*RR = 'A'
           P2*RR = 'B'
           P3*RR = 'C'
    / OVERLAY
      VAXIS = 0 TO 1.0 BY .10
      VZERO HZERO
   ;
TITLE6 'IN TERMS OF THE RELATIVE RISK PE/PC';
 PROC PLOT NOLEGEND;
      PLOT P1*OR = 'A'
           P2*OR = 'B'
           P3*OR = 'C'
    / OVERLAY
      VAXIS = 0 TO 1.0 BY .10
      VZERO HZERO
   ;
TITLE6 'IN TERMS OF THE ODDS RATIO PE*(1-PC) / PC*(1-PE)';
 PROC PLOT NOLEGEND;
      PLOT OR*DIFF = 'O'
           RR*DIFF = 'R'
    / OVERLAY
      VZERO HZERO
   ;
TITLE1
  ' PLOTS OF THE RELATIVE RISK AND ODDS RATIO VERSUS THE DIFFERENCE';
 PROC PLOT NOLEGEND;
      PLOT OR*RR = 'O'
    /
      VZERO HZERO
   ;
TITLE1 'PLOT OF THE ODDS RATIO VERSUS THE RELATIVE RISK';
 
 run;
