*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter4/tab4-4.sas                   *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: The SAS program that appears in Table 4.4.                *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
data one;
input k a b c d;
cards;
1  16  20  26  27
2   9   4   3   5
3  28  16  18  28
;
Title1 'Example 4.1: Ulcer Clinical Trial';
data two; set one;
keep i j k f;
**K=Stratum, I=Group, J=Response, F=Frequency;
i = 1; j = 1; f =a; output;
i = 2; j = 1; f =b; output;
i = 1; j = 2; f =c; output;
i = 2; j = 2; f =d; output;
proc freq; table k*(i j) / chisq nocol nopercent; weight f;
Title2
'Association Of Stratum By Group (K*I) And By Response (K*J)';
proc freq; table k*i*j / cmh; weight f;
Title2 'SAS Mantel-Haenszel Analysis';
run;
