*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter5/matchtab.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: SAS PROC FREQ analysis of matched 2 x 2 tables as         *;
*          described in Table 5.1 and applied to the data from       *;
*          Examples 5.3 and 5.4.                                     *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
%macro doit;
data two; set one;
Emember=1; Nmember=1; x=e; output;
Emember=2; Nmember =1; x=f; output;
Emember =1; Nmember =2; x=g; output;
Emember =2; Nmember =2; x=h; output;
proc freq; tables Emember* Nmember / all nopercent nocol agree;
     weight x; exact mcnem;
title2 'SAS PROC FREQ Analysis of Matched 2 x 2 tables';
run;
%mend;
 
data one; input e f g h;
**** note that the values of e and h are irrelevant;
cards;
10 2 8 10
;
title1 'Chapter 5, Example 5.3';
%doit;
 
data one; input e f g h;
**** note that the values of e and h are irrelevant;
cards;
10 80 55 10
;
title1 'Chapter 5, Example 5.4';
%doit;
 
run;
