*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter2/table23.sas                  *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: SAS Proc FREQ analysis of basic 2x2 tables.               *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
title1 'SAS PROC FREQ Analysis of 2x2 tables';
%macro doit;
data two; set one;
group=1; response=1; x=a; output;
group=2; response=1; x=b; output;
group=1; response=2; x=c; output;
group=2; response=2; x=d; output;
proc freq; tables group*response / all nopercent nocol; weight x;
%mend;
 
data one; input a b c d; cards;
53 40 47 60
;
title1 'Neuropathy Clinical Trial, Example 2.2';
%doit;
 
run;
