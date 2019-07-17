*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter4/ktab2x2.sas                  *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: This job invokes the macros for the analysis of the data  *;
*          in Examples 4.1, 4.6 and 4.24 using fixed and random      *;
*          effects models.                                           *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
*** The macros can be accessed either by "submitting" each job to sas before submitting this
    job, or by statements such as the following (remove the * first);
 
* filename freqanal '/jml/biostatmethods/chapter4/freqanal.sas';
* %include freqanal;
* filename freqrand '/jml/biostatmethods/chapter4/freqrand.sas';
* %include freqrand;
* filename randiter '/jml/biostatmethods/chapter4/randiter.sas';
* %include randiter;
* filename freqmh '/jml/biostatmethods/chapter4/freqmh.sas';
* %include freqmh;
 
*** To save paper, I recommend that you use the following;
 
options formdlim='-';
 
DATA ONE;
INPUT K N1 D1 N2 D2;
CARDS;
 1  42  16  47  20
 2  12   9   9   4
 3  46  28  44  16
*****;
TITLE1 'Ulcer Clinical trial Data (Example 4.1)';
 
** choose one of the following;
 
%freqanal; %freqrand;
**OR
**%freqmh;
 
DATA ONE;
INPUT K N1 D1 N2 D2;
CARDS;
 1 35  4  42  5
 2 21  4  31  13
 3 89  2  62  2
 4 73  8  45  9
*****;
TITLE1 'Religion and Mortality (Example 4.6)';
 
** choose one of the following;
 
%freqanal; %freqrand;
**OR
**%freqmh;
 
 
data one; input k d1 n1 d2 n2;
cards;
1 14 131 14 136
2 21 385 17 134
3 14 57 24 48
4 6 38 18 40
5 12 1011 35 760
6 138 1370 175 1336
7 15 506 20 524
8 6 108 2 103
9 65 153 40 102
*****;
TITLE1 'Meta-analysis of Diuretics and Pre-eclampsia (Example 4.24)';
 
*** fixed effects and one-step random effects analysis;
%freqanal;
%freqrand;
 
*** iterative random effects analysis;
%random;
 
data one; input k d1 n1 d2 n2;
cards;
1 24 25 22 26
2 35 37 35 39
3 31 36 38 41
4 46 53 42 57
5 60 73 51 79
6 39 53 32 52
*****;
TITLE1 'Example 5.2 - Frequency Matched Analysis';
 
*** fixed effects analysis;
%freqanal;
 
 
 
run;
