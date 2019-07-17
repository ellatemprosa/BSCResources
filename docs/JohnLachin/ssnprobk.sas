*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter3/ssnprobk.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: Used to generate Table 3.1 which presents the             *;
*          non-centrality parameter for the test for two             *;
*          proportions as a function of the probabilities in each    *;
*          group.                                                    *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
data one;
 
do j = 2 to 10;
 p1=j/10;
 mt=j-1;
do m = 1 to mt;
 p2=m/10;
 pb=(p1+p2)/2;
 k=abs(p1-p2)/sqrt(4*pb*(1-pb));
 rd=p1-p2;
 rr=p1/p2;
 or=p1*(1-p2)/(p2*(1-p1));
 output;
end;
end;
 
proc print;
 
run;
