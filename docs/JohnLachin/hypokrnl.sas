*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter9/hypokrnl.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: This program uses the macro smooth in kernelsm.sas to     *;
*          computes the kernel smoothed estimate of the intensity    *;
*          intensity function for the recurrent hypoglycemia counting*;
*          process described in Example 9.11 and presented in Figure *;
*          9.4.                                                      *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*---------- ----------------------------------------------------------*;
 
libname hypo '/jml/biostatmethods/datasets/hypoglycemia';
 
filename kernel '/jml/biostatmethods/chapter9/kernelsm.sas';
 
%include kernel;
 
data intimes;
 set hypo.hytimes;
 
*** data in days, band width is 9 months (in days), rate is per 100 patient years,
    plottime is in months;
 
%smooth(data=intimes,
   intime=day,
   width=(9/12)*365.25,
   ratek=100,
   ratetime=year,
   plottime=month,
   axis = vaxis = 0 to 120 by 20 haxis = 0 to 108 by 12);
 
run;
