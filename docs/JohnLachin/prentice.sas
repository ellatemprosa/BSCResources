*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter9/Prentice.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: This program reads the data in Table 9.9 from Prentice    *;
*          (1973) in a format suitable for use with the other macros *;
*          provided or with SAS procedures.                          *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*---------------------------------------------------------------------*;
 
data in; input patid time delta z1 z2 z3 z4 z5;
cards;
16   999  1  90    12     54     1     1
21   991  1  70     7     50     1     1
24   587  1  60     3     58     0     1
29   457  1  90     2     64     0     1
 2   411  1  70     5     64     1     0
25   389  1  90     2     62     0     1
28   357  1  70    13     58     0     1
 9   314  1  50    18     43     0     0
34   283  1  90     2     51     0     1
20   242  1  50     1     70     0     1
19   231  0  50     8     52     1     1
 3   228  1  60     3     38     0     0
30   201  1  80    28     52     1     1
13   141  1  30     4     63     0     0
 4   126  1  60     9     63     1     0
 5   118  1  70    11     65     1     0
17   112  1  80     6     60     0     1
22   111  1  70     3     62     0     1
 8   110  1  80    29     68     0     0
10   100  0  70     6     70     0     0
18    87  0  80     3     48     0     1
 7    82  1  40    10     69     1     0
 1    72  1  60     7     69     0     0
33    44  1  60    13     70     1     1
11    42  1  60     4     81     0     0
26    33  1  30     6     64     0     1
32    30  1  70    11     63     0     1
27    25  1  20    36     63     0     1
14    25  0  80     9     52     1     0
35    15  1  50    13     40     1     1
15    11  1  70    11     48     1     0
 6    10  1  10     5     49     0     0
12     8  1  40    58     63     1     0
23     1  1  20    21     65     1     1
31     1  1  50     7     35     0     1
;
title1 'Prentice carcinoma data';
 
* add additional SAS statements here;
 
 
RUN;
