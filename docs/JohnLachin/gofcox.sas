*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter9/gofcox.sas                   *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: is a macro written by Dr. Oliver Bautista to evaluate     *;
*          the proportional hazards assumption for each covariate    *;
*          in the model using the method of Lin (1991).              *;
*          Instructions for use of the macro are provided in the     *;
*          program as comments.                                      *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
/*------------------------------------------------------------------*/
/* File:    GOFCOX MACRO on OMBACT1(191)                            */
/* Updated: 12 May 1998 at 14:00:21                                 */
/* Purpose: Goodness-of-Fit test for the Cox PH regression model    */
/*          based on the procedure of D.Y. Lin (1991).              */
/* Code by: O.M. Bautista                                           */
/*                                                                  */
/* Required Parameters:                                             */
/*          Data    = input data set                                */
/*          Time    = event or censoring time                       */
/*          EvtFlag = indicator variable, 1=with event, 0=censored  */
/*          Var     = list of baseline/time dependent covariates    */
/*                                                                  */
/* Optional Parameters:                                             */
/*          rho      = see D.Y. Lin (1991).  Default=1.0            */
/*          tau      = see D.Y. Lin (1991).  Default=0.0            */
/*          printopt = controls the prinout of the IML routine      */
/*                     that computes the weighted parameter         */
/*                     estimates. Recognized values are:            */
/*                        0 : No results printed (default).         */
/*                        1 : Prints iteration history, as well as  */
/*                            summary of optimization start and     */
/*                            termination.                          */
/*                        2 : Prints initial and final parameter    */
/*                            estimates in addition to those in     */
/*                            option 1.                             */
/*                        3 : Prints values of termination criteria */
/*                            and other control parameters in addi- */
/*                            tion to those in option 2.            */
/*                        4 : Prints values of the parameter vector */
/*                            after each iteration, in addition     */
/*                            to those in option 3.                 */
/*                        5 : Prints values of the gradient vector  */
/*                            after each iteration, in addition     */
/*                            to those in option 4.                 */
/*                                                                  */
/*------------------------------------------------------------------*/
%MACRO  GOFCOX(Data     =      ,
               Time     =      ,
               EvtFlag  =      ,
               Vars     =      ,
               rho      = 1.0  ,
               tau      = 0.0  ,
               printopt = 0    );
 
  *---------------- Preliminary Data Processing -----------------*;
  /***************************************************************;
  /*               Extract the variable names                   */;
  /*------------------------------------------------------------*/;
     %Let N=0;
     %Do %While(%Scan(&Vars,&N+1)¬=);
       %Let N=%Eval(&N+1);
       %Let Var&N=%Scan(&Vars,&N);
     %End;
 
     Data Varnames;
        length PARAM $ 14 ;
          %DO i=1 %to &n;
            PARAM=upcase("&&Var&i"); output;
          %end;
     Data CovNames;
        length PARAM $ 9  ;
          %DO i=1 %to &n;
            PARAM=upcase("&&Var&i"); output;
          %end;
     run;
 
  /***************************************************************;
  /*    Drop observations with missing covariates, event flag,  */;
  /*    or event/censoring times.                               */;
  /*------------------------------------------------------------*/;
     Data _Use_ _drop_ ;
       set &Data; drop=0;
       if
          %DO i=1 %to &n;
            (&&Var&i > .Z) and
          %end;
            (&EvtFlag > .Z) and (&time > .Z)
       then output _use_;
       else drop=1;
       output _drop_;
     Data _Use_ ; set _Use_(drop=drop); run;
 
  /*******************************************************************/;
  /* Count Records with Missing Covariates, Event Flag, or Event Time*/;
  /*-----------------------------------------------------------------*/;
       data _drop_; set _drop_; if drop=1;
       proc summary data=_drop_;
         var drop;
         output out=_miss_ n=nmiss;
       proc datasets; delete _drop_;
       run;
 
  /***************************************************************;
  /*          Create Analysis Data Set(s)                       */;
  /*------------------------------------------------------------*/;
       %Let Nvars=%Eval(&N+2);
       data _temp_;
          set _use_ ;
          array usevars {&Nvars}
                &time &EvtFlag %DO i=1 %to &n; &&Var&i %end; ;
          array matvars {&Nvars} _var1 - _var&nvars.        ;
          do i=1 to &nvars ;
             matvars{i}=usevars{i};
          end;
          rename
                _var1=&time
                _var2=&EvtFlag
                %do i=3 %to &nvars;
                     %let j = %Eval(&i-2);
                     _var&i = &&Var&j
                %end;
             ;
          keep _var1 - _var&nvars. ;
       data _cntrl_;
         set _temp_;
       proc datasets;
        delete _temp_;
       run;
  *-------------  End Preliminary Data Processing ---------------*;
 
  *-----------------------*;
  *--- Run Proc PH Reg ---*;
  *-----------------------*;
      proc phreg data=&data outest=parms covout;
        model &time*&EvtFlag(0)=&vars /alpha=0.05 rl;
      run;
      data _Beta Cov_Beta;
        set parms;
        if _TYPE_="PARMS" then output _Beta;
        if _TYPE_="COV"   then output Cov_Beta;
        keep &vars;
      run;
  *------------------------------------------------*;
  * Run Lifetest to get Kaplan-Meier Estimate      *;
  * of survival function to be used as weights.    *;
  *------------------------------------------------*;
    proc lifetest data=_cntrl_ outsurv=surv noprint;
      time &time*&Evtflag(0);
    run;
    data surv;
     set surv;
     if _CENSOR_=0;
    run;
  *---------------------------------*;
  * Assign weights to original data *;
  *---------------------------------*;
    proc sort data=_cntrl_;
      by &time;
    data surv;
      merge _cntrl_(in=cntrl) surv(keep=&time survival);
      if cntrl;
      by &time;
    data surv;
      set surv;
      retain surv;
      if _N_=1 and survival = . then survival=1;
      if survival > . then surv=survival;
      keep &time surv;
    run;
 *------------------------------------------*;
 * Run IML Routines                         *;
 * Var names used:                          *;
 *     n = number of observations           *;
 *     p = number of covariates             *;
 *     X = Failure/Censoring time           *;
 *     d = 1 if X=Failure time              *;
 *     Z = Matrix of Covariates             *;
 *     W = Vector of Kaplan-Meier Survival  *;
 *         probabilities (combined groups). *;
 *------------------------------------------*;
 
  proc iml worksize=9000 symsize=1000 ;
    use _cntrl_ ; read all into input ; close _cntrl_;
    use surv    ; read all into evtime; close surv;
    n=nrow(input)  ; r=ncol(input)  ;
    X=input (¦,1¦) ; d=input (¦,2¦) ; Z=input (¦,3:r¦); p=ncol(Z);
    FT=evtime(¦,2¦); _rho_=&rho     ; _tau_=&tau      ;
    if _rho_=0 & _tau_>0 then W=( (1-FT)##_tau_ );
    if _rho_>0 & _tau_=0 then W=(    FT ##_rho_ );
    if _rho_=0 & _tau_=0 then W=J(n,1,1);
 
   *-----------------------------------*;
   *-- Start of function Definitions --*;
   *-----------------------------------*;
      start S0_Cox(Beta) global(X,d,Z,n,p);
         ZB = Z*Beta`; *- nx1 vector of linear predictors;
         do j=1 to n;
             _s0_ = _s0_ // sum( ( X >= X(¦j¦) )#exp(ZB) );
         end;
         S0=_s0_;
         return(S0);
      finish S0_Cox;
 
      start S1_Cox(Beta) global(X,d,Z,n,p);
         ZB = Z*Beta`; *- nx1 vector of linear predictors;
         do j=1 to n;
           tmp  = ( X >= X(¦j¦) )#exp(ZB)#Z ;
           _s1_ = _s1_ // tmp(¦+,¦);
         end;
         S1=_s1_;
         return(S1);
      finish S1_Cox;
 
      start S2_Cox(Xj,Beta) global (X,Z,n,p);
        hj=J(p,p,0);
        do k=1 to n;
           if X(¦k¦) >= Xj then do;
              Zk=Z(¦k,¦);
              hj=hj+(exp(Zk*Beta`)#(Zk`*Zk));
           end;
        end;
        S2j=hj;
        return(S2j);
      finish S2_Cox;
 
   *******************************************;
   ** Un-Weighted Cox PH Parameter Estimates *;
   *******************************************;
      *--------------------------------------------------*;
      *--- Objective Function of the Cox PH model -------*;
      *---      (Log of Partial Likelihood)       -------*;
      *--------------------------------------------------*;
        start F_Cox(Beta) global(X,d,Z,n,p);
           ZB = Z*Beta`; *- nx1 vector of linear predictors;
           S0 = S0_cox(Beta);
           f  = sum( d#( ZB-log(S0) ) );
           return(f);
        finish F_Cox;
 
      *-------------------------------------------*;
      *--- Gradient Vector of Cox PH model -------*;
      *-------------------------------------------*;
        start G_Cox(Beta) global(X,d,Z,n,p);
           ZB  = Z*Beta`; *- nx1 vector of linear predictors;
           S0_ = S0_cox(Beta); S1_ = S1_cox(Beta);
           tmp = d#( Z-(S1_#(S0_##(-1))) );
           g   = (tmp(¦+,¦)) ;
           return(g);
        finish G_Cox;
 
      *------------------------------------------*;
      *--- Hessian Matrix of Cox PH model -------*;
      *------------------------------------------*;
        start H_Cox(Beta) global(X,d,Z,n,p);
           S0=S0_cox(Beta); S1=S1_cox(Beta);
           hj=J(p,p,0);
           do j=1 to n ;
            if d(¦j¦)=1 then do;
               S1T_S1=S1(¦j,¦)`*S1(¦j,¦);
               S0j=S0(¦j¦);
               S0j2=S0j*S0j;
               S2=S2_Cox(X(¦j¦),Beta);
               hj=hj+( d(¦j¦)#( (-S2/S0j) + (S1T_S1/S0j2) ) );
            end;
           end;
           h = hj;
           return(h);
        finish H_Cox;
 
   ****************************************;
   ** Weighted Cox PH Parameter Estimates *;
   ****************************************;
      *-------------------------------------------------------*;
      *--- Objective Function of the Weighted Cox PH model ---*;
      *-------------------------------------------------------*;
        start F_WTCox(Beta) global(W,X,d,Z,n,p);
           ZB = Z*Beta`; *- nx1 vector of linear predictors;
           S0 = S0_cox(Beta);
           f  = sum( d#W#( ZB-log(S0) ) );
           return(f);
        finish F_WTCox;
 
      *----------------------------------------------------*;
      *--- Gradient Vector of the Weighted Cox PH model ---*;
      *----------------------------------------------------*;
        start G_WTCox(Beta) global(W,X,d,Z,n,p);
           ZB  = Z*Beta`; *- nx1 vector of linear predictors;
           S0_ = S0_cox(Beta); S1_ = S1_cox(Beta);
           tmp = d#W#( Z-(S1_#(S0_##(-1))) );
           g   = (tmp(¦+,¦)) ;
           return(g);
        finish G_WTCox;
 
      *---------------------------------------------------*;
      *--- Hessian Matrix of the Weighted Cox PH model ---*;
      *---------------------------------------------------*;
        start H_WTCox(Beta) global(W,X,d,Z,n,p);
           S0=S0_cox(Beta); S1=S1_cox(Beta);
           hj=J(p,p,0);
           do j=1 to n ;
            if d(¦j¦)=1 then do;
               S1T_S1=S1(¦j,¦)`*S1(¦j,¦);
               S0j=S0(¦j¦);
               S0j2=S0j*S0j;
               S2=S2_Cox(X(¦j¦),Beta);
               hj=hj+( W(¦j¦)#( (-S2/S0j) + (S1T_S1/S0j2) ) );
            end;
           end;
           h = hj;
           return(h);
        finish H_WTCox;
 
      *--------------------------------------------------*;
      *--- For use in Computing the Covariance of the ---*;
      *--- Weighted Parameter Estimates               ---*;
      *--------------------------------------------------*;
        start H_WT2Cox(Beta) global(W,X,d,Z,n,p);
           S0=S0_cox(Beta); S1=S1_cox(Beta);
           hj=J(p,p,0);
           W2=W#W;
           do j=1 to n ;
            if d(¦j¦)=1 then do;
               S1T_S1=S1(¦j,¦)`*S1(¦j,¦);
               S0j=S0(¦j¦);
               S0j2=S0j*S0j;
               S2=S2_Cox(X(¦j¦),Beta);
               hj=hj+(W2(¦j¦)#((-S2/S0j)+(S1T_S1/S0j2)));
            end;
           end;
           h = hj;
           return(h);
        finish H_WT2Cox;
 
   *---------------------------------*;
   *-- End of function Definitions --*;
   *---------------------------------*;
 
   Beta0=J(1,p,.)     ; *-- Initial Parameter Estimates  --*;
   opt=J(1,11,.)      ; *-- Null Options vector          --*;
   opt(¦1¦)=1         ; *-- Maximize objective function  --*;
   opt(¦2¦)=&printopt ;
   *-----------------------------------------------*;
   *--- Get Maximum Partial Likelihood Estimate ---*;
   *-----------------------------------------------*;
   /* call NLPNRA(rc,beta,"F_cox",beta0,opt) grd="G_Cox" hes="H_Cox" */;
   /* Obj=F_cox(Beta); _2logl=-2*obj; Grad=G_cox(Beta);              */;
    * Hess=H_cox(Beta);
    *  Beta_Cov=inv(-Hess);
    *   SE_Beta=sqrt(vecdiag(Beta_Cov));
    use _Beta;
     read all into Beta;
      close _Beta;
    use Cov_Beta;
     read all into Beta_Cov;
      close Beta;
        SE_Beta=sqrt(vecdiag(Beta_Cov))`;
   *---------------------------------------*;
   *--- Get Weighted Parameter Estimate ---*;
   *---------------------------------------*;
   if opt(¦2¦) > 0 then print "----- Weighted Parameter Estimate -----";
   call NLPNRA(rc,betaw,"F_WTCox",beta0,opt) grd="G_WTCox" hes="H_WTCox";
   /* Obj=F_WTCox(Betaw); _2logl=-2*obj; Grad=G_WTCox(Betaw); */;
      Hess=H_WTCox(Betaw);
       A=inv(-Hess);
        B=-H_WT2Cox(Betaw);
         BetawCov=A*B*A;
          SE_Betaw=sqrt(vecdiag(BetawCov))`;
   *--------------------------------------*;
   *--- Chi-Square for Goodness of Fit ---*;
   *--------------------------------------*;
   Diff     = Betaw - Beta;
   Diff_Cov = BetawCov - Beta_Cov;
   Diff_SE  = sqrt(vecdiag(Diff_Cov))`;
   X2       = Diff*inv(Diff_Cov)*Diff`;
   prob     = 1-probchi(X2,p);
   * print betaw , se_betaw, beta, se_beta, diff, diff_se, p;
   temp = (betaw`) ¦¦ (se_betaw`) ¦¦ (Beta`)  ¦¦ (se_beta`) ¦¦ (diff`) ¦¦
          (diff_se`) ¦¦ J(p,1,x2) ¦¦ J(p,1,p) ¦¦ J(p,1,prob);
   create params from temp;
     append from temp;
       close params;
  quit;
  data params; merge varnames params;
   rename col1=betaw col2=se_betaw col3=beta col4=se_beta
          col5=diff  col6=se_diff col7=x2 col8=df col9=prob;
  run;
   *--------------------------*;
   *--- Print Test Results ---*;
   *--------------------------*;
   %macro prt_test;
    data _null_;
     set params; file print header=newpage n=ps;
     _rho=&rho; _tau=&tau;
     if _N_=1 then do;
        put //// @75 "Wald Chi-Square"
              // @33 "Test of Proportional Hazards Assumption:"
                 + 2  x2 Best7. +1 "with " df 1.0 +1 "DF (p=" prob 6.4 ")"
            //// @43 "Weighted and Maximum Partial Likelihood Estimates"
             /// @45 "Weighted"
                 @65 "Maximum Likelihood"
                 @93 "Difference"
       OVERPRINT @40 21*"_" @64 21*"_" @88 21*"_"
              // @27  'Parameter'
                 @39  ' Estimate '
                 @51  'Std. Error'
                 @63  ' Estimate '
                 @75  'Std. Error'
                 @87  ' Estimate '
                 @99  'Std. Error'
              //
        ;
     end;
        put @27  Param
            @40  betaw    Best8.
            @52  se_betaw Best8.
            @64  beta     Best8.
            @76  se_beta  Best8.
            @88  diff     Best8.
            @100 se_diff  Best8. / ;
     return;
      newpage:
        put /// @37 61* "*"
              / @37 "*  GOODNESS-OF-FIT ANALYSIS OF THE COX PH" +1
                    "REGRESSION MODEL"
                @97 "*"
              / @37 "*" @47 "BASED ON A CLASS OF PARAMETER ESTIMATORS"
                @97 "*"
              / @37 "*" @59 "D.Y. Lin (1991)"
                @97 "*"
              / @37 "*" @57 "(rho=" _rho 4.2 +1 "tau=" _tau 4.2 ")"
                @97 "*"
              / @37 "*"
                @97 "*"
              / @37 "*" @49 "Programmed in SAS by: O.M. Bautista"
                @97 "*"
              / @37 61* "*";
      return;
    run;
   %mend;
   %prt_test;
  run;
%MEND GOFCOX  ;
