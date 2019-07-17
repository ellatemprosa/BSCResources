*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter9/nepha1c.sas                  *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: fits the PH model for the effect of the current annual    *;
*          mean HbA1c on the risk of developing microalbuminuria in  *;
*          the conventional treatment group of the secondary cohort  *;
*          of the DCCT. This was used to generate the results        *;
*          presented in Example 9.9.                                 *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
filename dcctneph '/jml/biostatmethods/datasets/nephdata.dat';
 
***********************************************************;
 
data neph; infile dcctneph;
input mask_pat int primary etdpatb neur aer0 neph2flg neph2vis
duration female age adult bcval5 hbael bmi mhba1-mhba9;
 
FLAG=NEPH2FLG;
TIME=NEPH2VIS/4;
if aer0 < 40;
if primary=0;
if int=0;
laer0=log(aer0);
retbase = primary;
IF FLAG NE .;
 
title1 "Risk of microalbuminuria as a Function of Time-Dependent Mean HbA1c";
 
PROC PHREG;
  MODEL TIME*FLAG(0)= LMHBA
        / TIES=DISCRETE ALPHA=0.05 RL;
  ARRAY MHBA (9) MHBA1-MHBA9;
  DO J=1 TO 9;
   IF time eq (J) THEN DO;
     LMHBA=LOG(MHBA(J));
   END;
  END;
 
RUN;
