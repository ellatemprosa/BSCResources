*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter9/hypotest.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: This program uses the macro AGtests in AGtests.sas to     *;
*          compute  the cumulative intensities and the Aalen-Gill    *;
*          tests for the recurrent hypoglycemia counting process     *;
*          described in Example 9.11.                                *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*---------- ----------------------------------------------------------*;
 
libname hypo '/jml/biostatmethods/datasets/hypoglycemia';
 
filename AGtests '/jml/biostatmethods/chapter9/AGtests.sas';
 
%include AGtests;
 
data intimes;
 set hypo.hytimes;
 
%AGtests;
 
run;
