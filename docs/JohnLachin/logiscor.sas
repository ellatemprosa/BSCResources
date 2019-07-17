*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter7/logiscor.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: computes the score vector, and the estimated expected     *;
*          information under the null hyothesis and computes the     *;
*          model score test as illustrated in Example 7.6. First     *;
*          run Renal.sas to set up the renal data set.               *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
data one; set renal;
y= micro24;
 
data inall; set one;
if y ne .;
 
%let nvar=6;    *** including the intercept;
%let nvar2=36;
%let ivars= int hbael yearsdm sbp female;
 
PROC LOGISTIC data=inall DESCENDING COVOUT OUTEST=modfit;
   MODEL Y = &ivars / COVB;
   OUTPUT OUT=TWO XBETA=XBETA PRED=PRED;
TITLE2 'LOGISTIC REGRESSION MODEL';
 
DATA BETAS; set modfit;
ARRAY BETA {&NVAR} BETA1-BETA&NVAR;
ARRAY X {&NVAR} INTERCEP &IVARS;
if _TYPE_ = 'PARMS';
DO J = 1 TO &NVAR;
  BETA(J) = X(J);
END;
*** PROC PRINT;
 
proc means data=inall n sum noprint; var y;
output out=four n=yn sum=ysum;
 
DATA THREE;
 if _n_=1 then set four (keep = yn ysum);
 SET TWO (KEEP = Y &ivars xbeta PRED) end=eof;
array x (&nvar) intercep &ivars;
retain predo;
RETAIN Io1-Io&NVAR2 0;                *** outer covariate products;
RETAIN uo1-uo&NVAR 0;
ARRAY uo {&NVAR} uo1-uo&NVAR;
ARRAY Io {&NVAR, &NVAR} Io1-Io&NVAR2;
intercep=1;
if _n_=1 then do;
 predo=ysum/yn;
 do i = 1 to &nvar;
 do j = 1 to &nvar;
  Io(i,j)=0;
 end;
 end;
end;
vo = predo*(1 - predo);
resido = (y - predo);
do i = 1 to &nvar;
  uio=x(i)*resido;
  uo(i) = uo(i) + uio;
do j = 1 to &nvar;
  Io(i,j)=Io(i,j) + x(i)*x(j)*vo;
end;
end;
 
if eof then output;
***proc print;
***title2
'Elements of the estimated information matrix under H0';
 
data uos; set three (keep = uo1-uo&nvar);
data Ios; set three (keep = Io1-Io&nvar2);
 
proc iml;
reset print;
use Ios;
read all into Io;
use uos;
read all into uo;
uo = uo`;
Io = shape(Io,&nvar,&nvar);
uo›1,1®=0;
score = uo`*inv(Io)*uo;
scorep=1-probchi(score,(&nvar-1));
 
quit;
 
RUN;
 
RUN;
