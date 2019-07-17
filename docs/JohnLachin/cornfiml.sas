*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter7/cornfiml.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: applies the binomial logit model to the Cornfield data    *;
*          as shown in Example 7.1 using the SAS PROC IML sample     *;
*          library routine LOGIT. The output was used to generate    *;
*          the iterative solution presented in Table 7.2.            *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
 
 /****************************************************************/
 /*          S A S   S A M P L E   L I B R A R Y                 */
 /*                                                              */
 /*    NAME: LOGIT                                               */
 /*   TITLE: LOGISTIC AND PROBIT REGRESSION FOR BINARY           */
 /*          RESPONSE MODELS                                     */
 /* PRODUCT: IML                                                 */
 /*  SYSTEM: ALL                                                 */
 /*    KEYS: MATRIX  REGR    SUGI6                               */
 /*   PROCS: IML                                                 */
 /*    DATA:                                                     */
 /*                                                              */
 /* SUPPORT: RHD                         UPDATE:                 */
 /*     REF:                                                     */
 /*    MISC:                                                     */
 /*                                                              */
 /****************************************************************/
 
proc iml ;
 
 /*----------Routine for Estimating Binary Response Models---------------*/
 /* Y is the binary response, X are regressors, wgt are count weights,   */
 /* model is choice of LOGIT PROBIT, parm has the names of the parameters*/
start binest;
 
   b0=log((y`*wgt)/((1-y)`*wgt));  *** specify intercept starting value;
   b = repeat(0,ncol(x),1);
   b›1,®=b0;
   oldb=b+1; /* starting values */
 
   do iter=1 to 20 while(max(abs(b-oldb))>1e-8);  oldb=b;
      z = x*b;
      run f;
      loglik =sum( ((y=1)#log(p) + (y=0)#log(1-p))#wgt);
      btransp = b`;
      print iter loglik btransp;
      w = wgt/(p#(1-p));
      xx = f # x;
      xpxi = inv(xx`*(w#xx));
      u=(xx`*(w#(y-p)));
      print b u xpxi;
      b = b + xpxi*u;
      end;
 
   p0 = sum((y=1)#wgt)/sum(wgt);  /* average response */
   loglik0 =sum( ((y=1)#log(p0) + (y=0)#log(1-p0))#wgt);
   chisq =  ( 2 # (loglik-loglik0));
   df    = ncol(x)-1;
   prob  = 1-probchi(chisq,df);
   print ,'Likelihood Ratio with Intercept-only Model' chisq df prob,;
   stderr = sqrt(vecdiag(xpxi));
   tratio = b/stderr;
   print parm b stderr tratio,,;
 
finish;
 
 /*---routine to yield distribution function and density---*/
start f;
if model='LOGIT'  then do; p=1/(1+exp(-z)); f=p#p#exp(-z);                 end;
if model='PROBIT' then do; p=probnorm(z);   f=exp(-z#z/2)/sqrt(8*atan(1)); end;
finish;
 
 /* Data From Cornfield (1962)*/
data={
0 0 1 10      ,
0 0 0 421     ,
0 1 1 10      ,
0 1 0 132     ,
1 0 1 38      ,
1 0 0 494     ,
1 1 1 34      ,
1 1 0 190    };
 
y=data›,3®;
wgt=data›,4®;
n=nrow(data);
x = repeat(1,n,1)¶¶( data›,{1 2}®);   /* intercept, local, female */
parm = {intercept, local, female};    /* names of regressors */
 
model={logit};  run binest;           /* run logit  model */
/*
model={probit}; run binest;           /* run probit model */
*/;
 
run;
