*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter5/macksetp.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: This uses the data set from Norm Breslow at U. Washington *;
*          This generates a data set for matched pairs within strata *;
*          defined by hypertensive versus not as described in        *;
*          Example 5.13. The data set is used by the job             *;
*          MackPair.sas to generate the overall (marginal) 2x2       *;
*          table, and the tables withiin hypertension strata         *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
filename mackdata '/jml/biostatmethods/datasets/bresmack.dat';
libname mack      '/jml/biostatmethods/datasets';
 
DATA ONE;
infile mackdata;
INPUT caseset CASE age GBDX HYPERT OBESE ESTROGEN DOSECODE
 DURcode NONEST;
 
DATA one; SET;
RETAIN CONTROLN 0;
 
IF CASE=1 THEN DO;
   CONTROLN=0;
END;
 
IF CASE=0 THEN CONTROLN=CONTROLN+1;
 
conjest=dosecode;
 if dosecode > 1 then conjest=1;
 
if obese = . then obese = 0;  **** combine missing with not obese;
 
title1 'Analysis of Breslow-Day endometrial cancer matched-pairs data';
 
data one; set one;
if case=1 or controln=1;
 
proc sort; by caseset;
/*;
proc print; var caseset controln conjest obese hypert;
*/;
 
data two; set one; by caseset;
retain e f g h casex controlx casez controlz;
 
array x (4) f e g h ;
**** e, f, g, h notation as in JL book, array index (cell) as in table below;;
 
if first.caseset then do j=1 to 4;
   x(j) = 0;                                   *         CELL                 ;
   casez=0;
   controlz=0;                                 *                     Control  ;
   casex=0;                                    *                   x=1     0  ;
   controlx=0;                                 *        x=1 (+)    e=2    f=1 ;
end;                                           *   case                       ;
if case=1 then do;                             *        0          g=4    h=3 ;
   casez=hypert;                               *                              ;
   casex=conjest;
end;
if case=0 then do;
   controlz=hypert;
   controlx=conjest;                           *         STRATUM              ;
end;
if last.caseset then do;                       *                    Control   ;
   if casez=0 then stratum = controlz + 1;     *                  z=1      0  ;
   if casez=1 then stratum = controlz + 3;     *       z=1 (+)  stratum=4  3  ;
   if casex=0 then cell = controlx + 3;        *  case                        ;
   if casex=1 then cell = controlx + 1;        *       0            2      1  ;
   if cell>. then x(cell)=1;                   *                              ;
   if cell=. then delete;                      *                              ;
   ***IF CELL=1 OR CELL=4 THEN;
       output;
end;
 
proc print; var caseset casex controlx casez controlz stratum cell e f g h;
proc sort; by stratum cell;
 
data mack.pairs; set two;
 
data three; set mack.pairs end=eof;
retain te tf tg th 0;
***if f=1 or g=1;
te = te + e;
tf = tf + f;
tg = tg + g;
th = th + h;
if eof then output;
 
proc print; var stratum te tf tg th;
 
data one; set three;
keep k e f g h;
k=stratum; e=te; f=tf; g=tg; h=th;
if f*g>0 then output;
 
%paired;
title1 'Mack unadjusted analysis';
 
data three; set mack.pairs; by stratum;
retain te tf tg th;
***if f=1 or g=1;
if first.stratum then do;
  te=0; tf=0;
  tg=0; th=0;
end;
te = te + e;
tf = tf + f;
tg = tg + g;
th = th + h;
if last.stratum then output;
 
proc print; var stratum casez controlz te tf tg th;
 
data one; set three;
keep k e f g h;
k=stratum; e=te; f=tf; g=tg; h=th;
if k = 1 or k=4;
if f*g>0 then output;
 
%Kpaired;
title1 'Mack stratified-adjusted analysis';
 
run;
 
 
 
run;
 
run;
