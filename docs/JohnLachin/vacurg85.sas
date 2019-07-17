*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter9/VACURG85.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: This program reads the data from the VA Cooperative       *;
*          Urology Research Group study of prostate cancer described *;
*          by Byar (1985) that are used in Problem 9.17. The         *;
*          variables described in Problem 9.17 are computed from the *;
*          raw data file VACURG85.dat that is provided in the        *;
*          datasets directory.                                       *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
filename vacurgin '/jml/biostatmethods/datasets/vacurg85.dat';
 
title1 'Prognostic variables for survival of prostatic cancer';
title2 'D.P. Byar, in Andrews, D.F., Herzberg, A.M. (1986), pp. 261-274';
 
data one; infile vacurgin;
input patid stage rx startm startd starty mosfu
status age wt pf hx sbp dbp ekg hg sz sg ap bm;
 
dead = status;
if status > 0 then dead=1;      * survival indicator variable;
 
if pf = -9999 then pf = .;      * Performance status;
if hx = -9999 then hx = .;      * History of cardiovascular disease;
if ekg = -9999 then ekg = .;    * EKG code;
if age = -9999 then age = .;    * age in years;
if sbp = -9999 then sbp = .;    * systolic blood pressure mm/hg;
if dbp = -9999 then dbp = .;    * diastolic blood pressure mm/hg;
if hg = -9999 then hg = .;      * Serum hemoglobin g/100 ml.;
if sz  = -9999 then sz  = .;    * tumor size in cm**2;
if sg  = -9999 then sg  = .;    * combined index of tumor stage and grade;
if wt = -9999 then wt = .;      * weight index kg - cm. height + 200;
if ap = 9999 then ap = .;       * alkaline phosphatase KA units, a measure of liver function;
if bm = 9 then bm = .;          * bone metastases 0=no, 1=yes;
 
ap = ap*(.1);                   * restore implied decimals;
hg = hg*(.1);                   * restore implied decimals;
sbp = sbp*10;                   * restore implied decimals;
dbp = dbp*10;                   * restore implied decimals;
 
lnap = log(ap);
if pf > 1 then pf = 1;            * cateforize pf =1 for bedridden, 0 for normal;
 
if 1<rx<4 then trt=0;
 If rx=4 then trt=1;              * high dose (1) versus low (0);
censor=(status=0);                * censoring indicator;
if (1<=age<75) then agegroup=1;   * age groups;
if (75<=age<80) then agegroup=2;
if age>=80 then agegroup=3;
if hx=0 then history=0;           * history of cardiovascular disease (1=yes,0-no);
if hx=1 then history=1;
if 0<=sz<30 then size=0;          * size of tumor, 0=small, 1=large;
if sz>=30 then size=1;
if 0<=sg<=10 then grade=0;        * combined index high (1) versus not (0);
if sg>10 then grade=1;
 
if ekg <=1 then do;               * EKG normal (0) versus abnormal (1);
ekgc = 0;
end;
if ekg >=2 then do;
ekgc = 1;
end;
 
rx1=0; rx2 = 0; rx3 = 0; rx4 = 0;     * indicator variables for each treatment group;
if rx = 1 then rx1 = 1;
if rx = 2 then rx2 = 1;
if rx = 3 then rx3 = 1;
if rx = 4 then rx4 = 1;
 
run;
 
proc means;
 
** add additional SAS statements as desired;
 
run;
