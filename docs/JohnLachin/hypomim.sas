*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter9/hypomim.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: This program fits the multiplicative intensity model to   *;
*          the recurrent  hypoglycemia events in the intensive group *;
*          patients of the secondary cohort of the DCCT as described *;
*          in Example 9.12. The input data set hypomimi contains the *;
*          event data for the intensive group subjects.              *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*---------- ----------------------------------------------------------*;
 
libname hypobook '/jml/biostatmethods/datasets/hypoglycemia';
 
data events; set hypobook.hypomimi;
if retbase='SCND'; ** select the secondary cohort only;
 
lhba1c2=lhba1c**2;
 
proc phreg;
strata phase2;
model (starts,stops)*hypoflg(0) = lhba1c lhba1c2 nprior;
 
run;
