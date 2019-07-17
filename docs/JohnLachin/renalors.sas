*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter7/renalors.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: performed the computations of the odds ratio for a unit   *;
*          increase in HbA1c as a function of the level of systolic  *;
*          blood pressure presented in Figure 7.1 of Example 7.13.   *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
data one;
bh=3.7478;  *** hbael coefficient;
bi=-0.0272; *** interaction coefficient;
 
*** covariance matrix of estimates;
sh=2.5208241577;  shi=-0.021237075;
                   si=0.0001804729;
 
sbpl=95;
sbpu=135;
range=sbpu-sbpl;
do j = 1 to range;
  sbp=sbpl + (j)*(sbpu-sbpl)/(range);
  b = bh + sbp*bi;
  or = exp(b);
  vb= sh + (sbp**2)*si + 2*sbp*shi;
  seb=sqrt(vb);
  orl=exp(b - 1.96*seb);
  oru=exp(b + 1.96*seb);
  output;
end;
 
proc print; var sbp b vb seb or orl oru;
proc plot;
  plot or*sbp = '*'
       orl*sbp = '-'
       oru*sbp = '-'
   / overlay vaxis = 0 to 4 by 0.5;
 
run;
