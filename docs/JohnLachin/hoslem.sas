*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/datasets/hoslem.sas                   *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: generates the data set used as the basis for the          *;
*          analyses in the text, Tables 7.8 and 7.9. The original    *;
*          source of the data was a table in the SAS Technical       *;
*          Report P-229, p. 465-6. The data were scanned and used    *;
*          in this program to generate the data set HLData.dat that  *;
*          is used in the program HLBwt.sas for the analyses shown   *;
*          in Chapter 7. Because the data were scanned from a        *;
*          secondary source, the data and analyses may differ from   *;
*          those shown by Hosmer and Lemeshow.                       *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
**** Hosmer-Lemeshow data from SAS Technical Report P-229, pgs. 465-67;
**** set up data set for used in Table 7.8;
 
filename HosLem '/jml/biostatmethods/DataSets/HLData.dat';
 
Data one; input id age case lwt smoke ht ui @@;
cards;
  25  16   1    130  0    0   0    143  16   0    110  0    0   0
 166  16   0    112  0    0   0    167  16   0    135  1    0   0
 189  16   0    135  1    0   0    206  16   0    170  0    0   0
 216  16   0     95  0    0   0     37  17   1    130  1    0   1
  45  17   1    110  1    0   0     68  17   1    120  1    0   0
  71  17   1    120  0    0   0     83  17   1    142  0    1   0
  93  17   0    103  0    0   0    113  17   0    122  1    0   0
 116  17   0    113  0    0   0    117  17   0    113  0    0   0
 147  17   0    119  0    0   0    148  17   0    119  0    0   0
 180  17   0    120  1    0   0     49  18   1    148  0    0   0
  50  18   1    110  1    0   0     89  18   0    107  1    0   1
 100  18   0    100  1    0   0    101  18   0    100  1    0   0
 132  18   0     90  1    0   1    133  18   0     90  1    0   1
 168  18   0    229  0    0   0    205  18   0    120  1    0   0
 208  18   0    120  0    0   0     23  19   1     91  1    0   1
  33  19   1    102  0    0   0     34  19   1    112  1    0   1
  85  19   0    182  0    0   1     96  19   0     95  0    0   0
  97  19   0    150  0    0   0    124  19   0    138  1    0   0
 129  19   0    189  0    0   0    135  19   0    132  0    0   0
 142  19   0    115  0    0   0    181  19   0    105  0    0   0
 187  19   0    235  1    1   0    192  19   0    147  1    0   0
 193  19   0    147  1    0   0    197  19   0    184  1    1   0
 224  19   0    120  1    0   0     27  20   1    150  1    0   0
  31  20   1    125  0    0   1     40  20   1    120  1    0   0
  44  20   1     80  1    0   1     47  20   1    109  0    0   0
  51  20   1    121  1    0   1     60  20   1    122  1    0   0
  76  20   1    105  0    0   0     87  20   0    105  1    0   0
 104  20   0    120  0    0   1    146  20   0    103  0    0   0
 155  20   0    169  0    0   1    160  20   0    141  0    0   1
 172  20   0    121  1    0   0    177  20   0    127  0    0   0
 201  20   0    120  0    0   0    211  20   0    170  1    0   0
 217  20   0    158  0    0   0     20  21   1    165  1    1   0
  28  21   1    200  0    0   1     30  21   1    103  0    0   0
  52  21   1    100  0    0   0     84  21   1    130  1    1   0
  88  21   0    108  1    0   1     91  21   0    124  0    0   0
 128  21   0    185  1    0   0    131  21   0    160  0    0   0
 144  21   0    110  1    0   1    186  21   0    134  0    0   0
 219  21   0    115  0    0   0     42  22   1    130  1    0   1
  67  22   1    130  1    0   0     92  22   0    118  0    0   0
  98  22   0     95  0    1   0    137  22   0     85  1    0   0
 138  22   0    120  0    1   0    140  22   0    130  1    0   0
 161  22   0    158  0    0   0    162  22   0    112  1    0   0
 174  22   0    131  0    0   0    184  22   0    125  0    0   0
 204  22   0    169  0    0   0    220  22   0    129  0    0   0
  17  23   1     97  0    0   1     59  23   1    187  1    0   0
  63  23   1    120  0    0   0     69  23   1    110  1    0   0
  82  23   1     94  1    0   0    130  23   0    130  0    0   0
 139  23   0    128  0    0   0    149  23   0    119  0    0   0
 164  23   0    115  1    0   0    173  23   0    190  0    0   0
 179  23   0    123  0    0   0    182  23   0    130  0    0   0
 200  23   0    110  0    0   0     18  24   1    128  0    0   0
  19  24   1    132  0    1   0     29  24   1    155  1    0   0
  36  24   1    138  0    0   0     61  24   1    105  1    0   0
 118  24   0     90  1    0   0    136  24   0    115  0    0   0
 150  24   0    110  0    0   0    156  24   0    115  0    0   0
 185  24   0    133  0    0   0    196  24   0    110  0    0   0
 199  24   0    110  0    0   0    225  24   0    116  0    0   0
  13  25   1    105  0    1   0     15  25   1     85  0    0   1
  24  25   1    115  0    0   0     26  25   1     92  1    0   0
  32  25   1     89  0    0   0     46  25   1    105  0    0   0
 103  25   0    118  1    0   0    111  25   0    120  0    0   1
 120  25   0    155  0    0   0    121  25   0    125  0    0   0
 169  25   0    140  0    0   0    188  25   0     95  1    0   1
 202  25   0    241  0    1   0    215  25   0    120  0    0   0
 221  25   0    130  0    0   0     35  26   1    117  1    0   0
  54  26   1     96  0    0   0     75  26   1    154  0    1   0
  77  26   1    190  1    0   0     95  26   0    113  1    0   0
 115  26   0    168  1    0   0    154  26   0   133   1    0   0
 218  26   0    160  0    0   0     16  27   1   150   0    0   0
  43  27   1    130  0    0   1    125  27   0   124   1    0   0
   4  28   1    120  1    0   1     79  28   1    95   1    0   0
 105  28   0    120  1    0   0    109  28   0   120   0    0   0
 112  28   0    167  0    0   0    151  28   0   140   0    0   0
 159  28   0    250  1    0   0    212  28   0   134   0    0   0
 214  28   0    130  0    0   0     10  29   1   130   0    0   1
  94  29   0    123  1    0   0    114  29   0   150   0    0   0
 123  29   0    140  1    0   0    190  29   0   135   0    0   0
 191  29   0    154  0    0   0    209  29   0   130   1    0   0
  65  30   1    142  1    0   0     99  30   0   107   0    0   1
 141  30   0     95  1    0   0    145  30   0   153   0    0   0
 176  30   0    110  0    0   0    195  30   0   137   0    0   0
 203  30   0    112  0    0   0     56  31   1   102   1    0   0
 107  31   0    100  0    0   1    126  31   0   215   1    0   0
 163  31   0    150  1    0   0    222  31   0   120   0    0   0
  22  32   1    105  1    0   0    106  32   0   121   0    0   0
 134  32   0    132  0    0   0    170  32   0   134   1    0   0
 175  32   0    170  0    0   0    207  32   0   186   0    0   0
;
 
proc sort; by age descending case;
 
DATA two; SET one; by age descending case;
RETAIN matchset casen CONTROLN 0;
 
if first.age then do;
   casen=0; controln=0; matchset=matchset+1;
   end;
 
IF CASE=1 THEN DO;
   casen=casen+1;
   CONTROLN=0;
   END;
 
IF CASE=0 THEN do;
   CONTROLN=CONTROLN+1;
   casen=0;
   end;
 
time = 2 - case;
 
proc print; var matchset casen controln lwt smoke ht ui case time;
 
data three; set two;
file hoslem;
put matchset casen controln lwt smoke ht ui case time;
 
 
run;