*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter9/fh-cgdph.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: reads the data from Fleming and Harrington (1991) with    *;
*          the times of recurrent infections among children in a     *;
*          clinical trial of interferon versus placebo. This uses    *;
*          the data set <I>FH-cgd.dat</I> that is suitable for use   *;
*          with PHREG to fit multiplicative intensity proportonal    *;
*          intensity models.                                         *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
title1 'CGD data from Fleming and Harrington (1991)';
 
 Filename cgddata '/jml/biostatmethods/datasets/FH-CGD.dat';
 
data one; infile cgddata;
input id rdt idt z1 z2 z3 z4 z5 z6 z7 z8 z9 t1 t2 d s;
trt=0; if z1=1 then trt=1;
delta=0; if d=1 then delta=1;
time=t1;
 
 
*** additional statements here as necessary, such as PHreg;
 
 
run;
