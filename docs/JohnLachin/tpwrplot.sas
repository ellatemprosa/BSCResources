*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter3/tpwrplot.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: Plot the power function for the test for two means for    *;
*          different sample sizes as a function of the non-central   *;
*          factor for different sample sizes. A variation of this    *;
*          was used to generate Figure 3.2.                          *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
*
*     PARAMETERS ARE:
*        TA = ALPHA CRITICAL VALUE, 1.645 FOR ALPHA 0,05 ONE-SIDED
*        N  = TOTAL SAMPLE SIZE
*        QC = SAMPLE FRACTION FOR THE C GROUP, NC = QC*N
*        QE = SAMPLE FRACTION FOR THE E GROUP, NE = QE*N
*        K  = ¦MU(E) - MU(C)¦ / 2*SIGMA
*
 ;
*****  specify output filename;
     filename out '/jml/biostatmethods/chapter3/Tpwrplot.out';
 %MACRO POWR;
       PZT =  SQRT(N)*K - TA;
       PZTL=PROBNORM(PZT );
 %MEND POWR;
 
 
 DATA ONE;
 TITLE1 'POWER CURVES FOR THE TEST OF A DIFFERENCE IN MEANS';
 TITLE2 '   ';
 TITLE3 '      QC = 0.5,  QE = 0.5,   ZA = 1.96';
 TITLE4 '   ';
 TITLE5 'A : N = 100,       B : N = 200,         C : N = 400';
 *****************************************************************;
 TA = PROBIT(0.975);
 DO L=1 TO 30;
 K=L/100;
 
 ***** specify parameters;
 N=100;
 %POWR;
 P1=PZTL;
 
 ***** specify parameters;
 N=200;
 %POWR;
 P2=PZTL;
 
 ***** specify parameters;
 N=400;
 %POWR;
 P3=PZTL;
 
 FILE OUT;
   PUT K P1 P2 P3;
 OUTPUT;
 END;
 PROC PRINT;
 PROC PLOT NOLEGEND;
      PLOT P1*K = 'A'
           P2*K = 'B'
           P3*K = 'C'
    / OVERLAY
      VAXIS = 0 TO 1.0 BY .10
      VZERO HZERO
   ;
 run;
