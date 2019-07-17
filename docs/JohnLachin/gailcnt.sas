*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter8/gailcnt.sas                  *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: This reads the data from Gail, Santner and Brown, 1980    *;
*          that are used in Problem 8.4. In their paper, the actual  *;
*          times start at 60. The data here start with 60 as time    *;
*          zero. A11 animals were exposed for 122 days (182-60).     *;
*          group 1 = retinoid, 2= control                            *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
Data allevts; input group m;
rate=m/122;
lndays=log(122);
patient=_n_;
im=m; if im> 1 then im=1;
nm=m; futime=122;
Cards;
1   1
1   0
1   2
1   1
1   4
1   3
1   6
1   1
1   1
1   5
1   2
1   1
1   5
1   2
1   3
1   4
1   5
1   5
1   1
1   2
1   6
1   0
1   1
0   7
0  11
0   9
0   2
0   9
0   4
0   6
0   7
0   6
0   1
0  13
0   2
0   1
0  10
0   4
0   5
0  11
0  11
0   9
0  12
0   1
0   3
0   1
0   3
0   3
;
 
proc sort; by group;
proc univariate; var rate; by group;
 
*** add additional programming statements;
 
run;
