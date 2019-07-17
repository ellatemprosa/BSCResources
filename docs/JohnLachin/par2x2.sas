*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter2/par2x2.sas                   *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: computes the population attributable risk and C.I. using  *;
*          Walter's limits, Leung-Kupper limits, and the logit       *;
*          limits derived in the text.                               *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
options ls=80 formdlim='-';;
 
%macro par;
data two; set one;
p1=d1/n1; q1=1-p1;
p2=d2/n2; q2=1-p2;
m=d1+d2; m1=m;
par=(d1-n1*p2)/m;
N=n1+n2;
 
*** Walter 1978 variance in his a,b,c,d notation (different from mine);
a=d1; b=n1-d1; c=d2; d=n2-d2;  *** Walter (1978) notation;
vpar=( c*N* ( (a*d*(N-c)) + (b*(c**2)) ) ) / ( ((a+c)**3) * ((c+d)**3) );
separ=sqrt(vpar);
za=probit(0.975);
parcl = par - za*separ;
parcu = par + za*separ;
 
*** JL variance;
lpar=log(par/(1-par));
vlp=( p1/ ((p1-p2)**2) )*( (n2*p2*q1 + n1*p1*q2) / (n1*n2*p2) );
selp=sqrt(vlp);
lpcl = lpar - za*selp;
lpcu = lpar + za*selp;
pcl=exp(lpcl)/(1+exp(lpcl));
pcu=exp(lpcu)/(1+exp(lpcu));
 
*** Leung-Kupper variance and limits in Walter notation;
vlp_LK=( (a+c)*(c+d) / (((a*d)-(b*c))**2) )*
       ( ( (a*d*(N-c)) + ((c**2)*b) )/(N*c) );  *** L-K notation = Walter;
selp_LK=sqrt(vlp_LK);
lpcl_LK = lpar - za*selp_LK;
lpcu_LK = lpar + za*selp_LK;
pcl_LK=exp(lpcl_LK)/(1+exp(lpcl_LK));
pcu_LK=exp(lpcu_LK)/(1+exp(lpcu_LK));
 
proc print; var par vpar separ parcl parcu;
title3 'Walter 1978 variance and confidence limits';
proc print; var lpar vlp selp lpcl lpcu pcl pcu;
title3 'JL logit variance and confidence limits';
proc print; var lpar vlp_LK selp_LK lpcl_LK lpcu_LK pcl_LK pcu_LK;
title3 'Leung-Kupper variance and confidence limits';
run;
%mend;
 
title1 'population attributable risk and confidence limits';
**** 1 = exposed, 2 = non-exposed;
data one; input d1 n1 d2 n2;
cards;
72 756 20 573
;
title2 'Cornfield cholesterol data';
%par;
 
data one; input d1 n1 d2 n2;
cards;
72 756 40 573
;
title2 'Hypothetical cholesterol data';
%par;
 
data one; input d1 n1 d2 n2;
cards;
72 756 10 573
;
title2 'Hypotheticalcholesterol data';
%par;
 
 
run;
