*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter6/newton.sas                   *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: The macro %Newton that provides the iterative solution    *;
*          of a scalar function using either Newton-Raphson or       *;
*          Fisher scaling. See Fisher Recomb.sas for the use of this *;
*          macro to solve for the recombination fraction in the      *;
*          Fisher maize data using either Newton-Raphson or Fisher   *;
*          scoring. You must also specify: 1) A macro %FOFX(x) which *;
*          evaluates the function of x to be solved and returns the  *;
*          value using the variable name "yofx". FOFX must be such   *;
*          that the value is zero at the solution. 2) A macro        *;
*          %DFOFX(x) which evaluates the derivative of the function  *;
*          at x and returns the value using the variable name        *;
*          "dyofx". For Newton-Raphson the hessian is used, for      *;
*          Fisher scoring the negative expected information is used. *;
*          3) A macro variable ERROR which gives the desired error   *;
*          of the result which when reached terminates the           *;
*          iteration. 4) A starting value in the variable x0 (not a  *;
*          macro variable). 5) A macro variable maxitn that fixes    *;
*          the maximum number of iterations at which the program     *;
*          aborts. Macro variables are defined using a %let          *;
*          statement, e.g. %let error=0.00000001.                    *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
%MACRO Newton;
**
**  Macro Newton for the iterative solution of a **Scalar** function
**
**  You must also specify:
      1.  A macro %FOFX(x) which evaluates the function of x to be solved
          and returns the value using the variable name "yofx".
          FOFX must be such that the value is zero at the solution.
      2.  A macro %DFOFX(x) which evaluates the derivative of the function at x
          and returns the value using the variable name "dyofx".
          For Newton-Raphson the hessian is used, for Fisher scoring the negative expected
          information is used.
      3.  A macro variable ERROR which gives the desired error of the result which
          when reached terminates the iteration.
      4.  A starting value in the variable x0 (not a macro variable).
      5.  A macro variable maxitn that fixes the maximum number of iterations at which
          the program aborts.
    Macro variables are defined using a %let statement, e.g. %let error=0.00000001;
**;
 
    xlast = x0;
    %fofx(xlast);
    ylast = yofx;
    %dfofx(xlast);
    dylast = dyofx;
    itn=0;
 
startit:
 
    itn=itn+1;
    xnext = xlast - (ylast / dylast);
    %fofx(xnext);
    ynext  = yofx;
    %dfofx(xnext);
    dynext = dyofx;
 
**;
 
output;
    if abs(ynext) < &error then go to stopit;
    if itn=&maxitn then go to abort;
 
**;
      ylast = ynext;
      dylast= dynext;
      xlast = xnext;
 
go to startit;
 
abort: put 'terminated after max # iterations reached';
 
stopit:
    xroot=xnext;
    output;
 
%mend Newton;
 
run;
