*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter9/fh-cgdag.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: reads the data from Fleming and Harrington (1991) with    *;
*          the times of recurrent infections among children in a     *;
*          clinical trial of interferon versus placebo. This uses    *;
*          the data set <I>cgdtimes</I> that is suitable for use     *;
*          with the macro <I>AGtests</I> to compute Aalen-Gill test  *;
*          statistics for recurrent events.                          *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
title1 'CGD data from Fleming and Harrington (1991)';
 
 libname biosbook '/jml/biostatmethods/datasets';
 
data one; set biosbook.cgdtimes;
 
*** additional statements here as necessary, such as AGtests;
 
 
run;
