*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter8/table81.sas                  *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: generates the data set displayed in Table 8.1 derived     *;
*          from the data set dccthypo.dat. Note that all             *;
*          computations are based on the latter data set, not        *;
*          Table81.dat.                                              *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
options ls=80 nodate;
 
filename table81 '/jml/biostatmethods/chapter8/table81.dat';
filename dccthypo '/jml/biostatmethods/datasets/hypoglycemia/dccthypo.dat';
 
data allevts;
infile dccthypo;
input grp nevents fuday iu duration
 female adult bcval5 hbael hxcoma obweight;
 
  ievents = 0; if nevents > 0 then ievents=1;
  fuyears = fuday/365.25;
  rate =  nevents/fuyears;
  lnyears = log(fuyears);
  insulin = iu/obweight;
  if grp=1 then group='Int ';
  if grp=0 then group='Conv';
  Label grp='tx group (1=Int, 0=Conv)'
        group='treatment group: Int or Conv'
        nevents='# severe hypoglycemia episodes'
        ievents='any hypoglycemia episodes (1=Y, 0=N)'
        fuday='total days of follow-up in the study'
        fuyears='total years of follow-up in the study'
        lnyears='log years of follow-up'
        rate='rate of episodes per year of follow-up'
        insulin='insulin units per kg weight'
        duration='months diabetes duration'
        female='female (1) or male (0)'
        adult='adult (1) or adolescent (0)'
        bcval5='C-peptide in pmol/mL'
        hbael='level of HbA1c at initial screening'
        hxcoma='Prior history of coma/seizure'
        obweight='body weight in kg.';
 
data test; set allevts;
If _n_ < 50;
proc print;
 
data two; set allevts;
 
file table81;
put (grp nevents fuyears rate insulin duration
 female adult bcval5 hbael hxcoma) (2 * 2. 3 * 9.5 5. 3. 4. 2 * 6.2 3.);
 
run;
