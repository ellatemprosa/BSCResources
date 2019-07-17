*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter7/prostate.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: includes the prostate cancer data set from Collett        *;
*          (1991) presented in Table 7.11. This data set was         *;
*          obtained (downloaded) from the SAS online data sets for   *;
*          the SAS book Logistic Regression Examples Using the SAS   *;
*          System, where it appears on p. 17-18.                     *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
data prostate;
title 'Prostate Data';
input case age acid xray size grade nodalinv @@;
lacd=log(acid);
datalines;
1  66  .48 0 0 0 0  2  68  .56 0 0 0 0  3  66  .50 0 0 0 0
4  56  .52 0 0 0 0  5  58  .50 0 0 0 0  6  60  .49 0 0 0 0
7  65  .46 1 0 0 0  8  60  .62 1 0 0 0  9  50  .56 0 0 1 1
10 49  .55 1 0 0 0  11 61  .62 0 0 0 0  12 58  .71 0 0 0 0
13 51  .65 0 0 0 0  14 67  .67 1 0 1 1  15 67  .47 0 0 1 0
16 51  .49 0 0 0 0  17 56  .50 0 0 1 0  18 60  .78 0 0 0 0
19 52  .83 0 0 0 0  20 56  .98 0 0 0 0  21 67  .52 0 0 0 0
22 63  .75 0 0 0 0  23 59  .99 0 0 1 1  24 64 1.87 0 0 0 0
25 61 1.36 1 0 0 1  26 56  .82 0 0 0 1  27 64  .40 0 1 1 0
28 61  .50 0 1 0 0  29 64  .50 0 1 1 0  30 63  .40 0 1 0 0
31 52  .55 0 1 1 0  32 66  .59 0 1 1 0  33 58  .48 1 1 0 1
34 57  .51 1 1 1 1  35 65  .49 0 1 0 1  36 65  .48 0 1 1 0
37 59  .63 1 1 1 0  38 61 1.02 0 1 0 0  39 53  .76 0 1 0 0
40 67  .95 0 1 0 0  41 53  .66 0 1 1 0  42 65  .84 1 1 1 1
43 50  .81 1 1 1 1  44 60  .76 1 1 1 1  45 45  .70 0 1 1 1
46 56  .78 1 1 1 1  47 46  .70 0 1 0 1  48 67  .67 0 1 0 1
49 63  .82 0 1 0 1  50 57  .67 0 1 1 1  51 51  .72 1 1 0 1
52 64  .89 1 1 0 1  53 68 1.26 1 1 1 1
;
 
* Fit the model of interest;
 
run;
