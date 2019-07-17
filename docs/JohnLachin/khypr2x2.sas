*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter6/khypr2x2.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: computes MLE of conditional odds ratio for K 2x2 tables   *;
*          using the Ulcer clinical trial data as presented in       *;
*          Example 6.7. This program uses a modification of the      *;
*          macro %Newton that performs the Newton-Raphson iterative  *;
*          solution for a scalar function. In this version, %FOFX    *;
*          also evaluates the derivative of the function at x and    *;
*          returns the value using the variable name                 *;
*          "dyofx". A macro DFOFX is not used.                       *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
options formdlim='-';
 
%MACRO Newton;
**;
 
    xlast = x0;
    %fofx(xlast);
    ylast = yofx;
    dylast = dyofx;
    itn=0;
 
startit:
 
    itn=itn+1;
    xnext = xlast - (ylast / dylast);
    %fofx(xnext);
    ynext  = yofx;
    dynext = dyofx;
 
**;
     k=0;
***; output three;
 
    if abs(ynext) < &error then go to stopit;
    if itn=&maxitn then go to abort;
 
**;
      ylast = ynext;
      dylast= dynext;
      xlast = xnext;
      vlast=-1/dylast;
 
go to startit;
 
abort: put 'terminated after max # iterations reached';
 
stopit:
    xroot=xnext;
 
%mend Newton;
 
run;
 
*** define the macro for the hypergeometric probability used in this computation;
 
%macro tablep;
  cp = probhypr(n, n1, m1, ia, or);
  if ia>al then p = cp - probhypr(n, n1, m1, ia-1, or);
  else p=cp;  *** for lower extreme table;
%mend;
 
%macro fofx(newor);
yofx=0;
dyofx=0;
or=&newor;
do k = 1 to &nstrat;
  a=as(k); b=bs(k); c=cs(k); d=ds(k);
  m1=a+b; m2=c+d; n1=a+c; n2=b+d; N=n1+n2;
  al=max(0, m1-n2); au=min(m1,n1);
  ea=0; Easq=0;
  do ia = al to au;
     %tablep;
     ea = ea + ia*p;
     easq=easq+ (ia*ia*p);
  end;
  u=a-Ea;
  va = easq - (ea**2);
  du=-va;
  yofx= yofx + u;
  dyofx = dyofx + du;
  output three;
end;
%mend;
 
%macro doit;
data two; set one end=eof;
retain a1-a&nstrat b1-b&nstrat c1-c&nstrat d1-d&nstrat;
array as (&nstrat) a1-a&nstrat;
array bs (&nstrat) b1-b&nstrat;
array cs (&nstrat) c1-c&nstrat;
array ds (&nstrat) d1-d&nstrat;
as(k)=xa;
bs(k)=xb;
cs(k)=xc;
ds(k)=xd;
if eof then output;
 
proc print;
 
data three test; set two;
array as (&nstrat) a1-a&nstrat;
array bs (&nstrat) b1-b&nstrat;
array cs (&nstrat) c1-c&nstrat;
array ds (&nstrat) d1-d&nstrat;
 
x0=1.0;  *** initial value for OR;
%let maxitn=25;
%let error=0.00005;
%newton
vlast=-1/dylast;
 
proc print data=three; var itn xlast ylast dylast vlast;
proc sort; by itn k;
 
data four; set three end=eof;
keep lastitn;
lastitn=itn;
if eof then output;
 
data five; if _n_=1 then set four;  set three;
if itn=lastitn;
if k > 0;
 
data six; set five end=eof;
retain tv 0;
tv=tv + va;
if eof then do;
  mleor=xnext;
  mlelor=log(xnext);
  vlor=1/tv;
  output;
end;
 
proc print; var itn mleor mlelor vlor;
 
%mend;
 
DATA ONE;
INPUT K xN1 xa xN2 xb;
xc= xn1-xa; xd=xn2-xb;
CARDS;
 1  42  16  47  20
 2  12   9   9   4
 3  46  28  44  16
*****;
TITLE1 'Example 6.7, using Ulcer Clinical Trial Data';
 
%let nstrat=3;
%doit;
 
run;
