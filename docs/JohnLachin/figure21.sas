*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter2/figure21.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: generates figures relating the odds ratio, relative risk  *;
*          and risk difference when one measure is held constant     *;
*          over a range of values of p2. A version of this was used  *;
*          to generate Figure 2.1 where the risk difference was      *;
*          held constant.                                            *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
Data one;
do i=10 to 100;
  p2=i/100;
  rd=-0.1;
  p1=p2+rd;
  rr=p1/p2;
  or=rr*(1-p2)/(1-p1);
  output;
end;
 
proc print;
proc plot;
  plot rd*p2 = 'D'
       rr*p2 = 'R'
       or*p2 = 'O'
  / overlay haxis = 0 to 1;** vaxis = 0 to 4;
title 'risk difference held to -0.1 over the range of p2';
run;
 
Data one;
do i=1 to 100;
  p2=i/100;
  rr=0.5;
  p1=p2*rr;
  rd=p1-p2;
  or=rr*(1-p2)/(1-p1);
  output;
end;
 
proc print;
proc plot;
  plot rd*p2 = 'D'
       rr*p2 = 'R'
       or*p2 = 'O'
  / overlay haxis = 0 to 1;** vaxis = 0 to 5;
title 'relative risk held to 0.5 over the range of p2';
 
run;
 
Data one;
do i=0 to 100;
  p2=i/100;
  or=0.5;
  p1=p2*or/((1-p2)+(or*p2));
  rd=p1-p2;
  rr=p1/p2;
  if p1<=1.0 then output;
end;
 
proc print;
proc plot;
  plot rd*p2 = 'D'
       rr*p2 = 'R'
       or*p2 = 'O'
  / overlay haxis = 0 to 1;** vaxis = 0 to 2;
title 'Odds ratio held to 0.5 over the range of p2';
 
run;
