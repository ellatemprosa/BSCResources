*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter5/mackpair.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: generates the marginal summary 2x2 table shown in         *;
*          Example 5.8 and calls macro %paired, and also generates   *;
*          the tables within strata defined by the baseline          *;
*          hypertension status as described in Example 5.13. The     *;
*          program calls %Kpaired but the analysis fails due to a    *;
*          zero frequency within one stratum.                        *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
*** this uses the data set from Norm Breslow at U. Washington;
 
libname mack '/jml/biostatmethods/datasets';
 
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
 
title1 'Example 5.8. Mack et al. marginal unadjusted analysis';
%paired;
 
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
 
title1 'Example 5.13. Mack et al. stratified-adjusted analysis by hypertension status';
%Kpaired;
 
run;
