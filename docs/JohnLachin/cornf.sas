*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter7/cornf.sas                    *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: The SAS Program for the Framingham CHD Data presented     *;
*          in Table 7.3. The output from this program appears in     *;
*          Table 7.4. This also contains the interaction model       *;
*          described in Example 7.10                                 *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
options nodate formdlim='-'; *** ls=65;
 
Title1 'Cornfield Framingham data on cholesterol and blood pressure';
data one; input hichol hisbp chd frequncy;
interact=hichol*hisbp;
 cards;
0 0 1 10
0 0 0 421
0 1 1 10
0 1 0 132
1 0 1 38
1 0 0 494
1 1 1 34
1 1 0 190
;
proc logistic descending;
  model chd = hichol hisbp / rl;
  weight frequncy;
title2 'main effects models';
 
proc logistic descending;
  model chd = hichol hisbp interact/ rl;
  weight frequncy;
  test hichol=interact=0;
  test hisbp=interact=0;
title2 'interaction effects model';
 
run;
 
run;
