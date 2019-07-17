*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter8/frch.sas                     *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: reads the data from Frome and Checkoway, 1985, that are   *;
*          used in Problem 8.6.                                      *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
data one; input age $ d1 t1 d2 t2;
cards;
15   4 181343   1 172675
25  38 146207  16 123065
35 119 121374  30  96216
45 221 111353  71  92051
55 259  83004 102  72159
65 310  55932 130  54722
75 226  29007 133  32185
85  65   7538  40   8328
;
 
proc means sum; var d1 t1 d2 t2;
output out=two sum=d1 t1 d2 t2;
 
data two; set two;
age='all';
drop _type_ _freq_;
 
data three; set one two;
proc print;
 
run;
