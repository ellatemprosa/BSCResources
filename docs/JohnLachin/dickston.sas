*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter7/dickston.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: Program for the logistic regression analysis of the       *;
*          unmatched retrospective study of Dick and Stone (1973)    *;
*          presented in Example 7.5.                                 *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
options nodate ls=65;
 
data one; input hl sm ht ihd frequncy @@; cards;
0 0 0 1 15 0 0 0 0 82
0 0 1 1 10 0 0 1 0 37
0 1 0 1 39 0 1 0 0 81
0 1 1 1 23 0 1 1 0 28
1 0 0 1 18 1 0 0 0 16
1 0 1 1 7  1 0 1 0 12
1 1 0 1 19 1 1 0 0 19
1 1 1 1 15 1 1 1 0 8
;
proc logistic descending;
  model ihd = hl sm ht / rl;
  weight frequncy;
run;
 
run;
