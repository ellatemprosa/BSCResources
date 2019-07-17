*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter8/fhcntan.sas                  *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: reads the file fhcnt.dat and defines additional variables *;
*          that can be used for analyses of numbers of events using  *;
*          rates and Poisson models as requested in Problem 8.10.    *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
title1 'CGD data from Fleming and Harrington (1991)';
 
filename fhcgd '/jml/biostatmethods/chapter8/fhcnt.dat';
 
data one;
infile fhcgd;
input id z1 z2 z3 z4 z5 z6 z7 z8 z9 nevents futime;
trt=0; if z1=1 then trt=1;
inherit=z2-1;      *** 0=x-linked, 1=autosomal recessive;
age=z3;
height=z4;
weight=z5;
steroids=(3-z6)-1; ***0=no corticosteroids on entry, 1=yes;
antibio= (3-z7)-1; ***0=no antiobiotics on entry, 1=yes;
female = z8-1;
NIH = (z9=1);
Amstdm = (z9=3);
EurOther = (z9=4);
 
logtime=log(futime);
 
*** other statements here;
 
run;
