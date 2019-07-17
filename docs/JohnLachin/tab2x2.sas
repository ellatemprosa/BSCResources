*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter2/tab2x2.sas                   *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: SAS Proc FREQ analysis of basic 2x2 tables.               *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
options nodate linesize=65 pagesize=67;
 
%macro doit;
data two; set one;
group=1; response=1; x=a; output;
group=2; response=1; x=b; output;
group=1; response=2; x=c; output;
group=2; response=2; x=d; output;
proc freq; tables group*response / all nopercent nocol;
     weight x;
title2 'SAS PROC FREQ Analysis of 2 x 2 tables';
%mend;
 
data one; input a b c d; cards;
53 40 47 60
;
title1 'Neuropathy Clinical Trial, Example 2.2';
%doit;
 
data one; input a b c d;
cards;
7 8 12 2
;
title1 'Chapter 2, Example 2.7';
%doit;
 
data one; input a n1 b n2;
keep a b c d;
c = n1 - a;
d = n2 - b;
cards;
60 113 36 90
;
title1 'Chapter 2, Example 2.10';
%doit;
 
data one; input a b c d;
cards;
72 20 684 553
;
title1 'Example 2.11, Cornfield CHD example';
%doit;
 
data one; input a n1 b n2;
keep a b c d;
c = n1 - a;
d = n2 - b;
cards;
154 200 80 200
;
title1 'Chapter 5, Example 5.1';
%doit;
 
run;
