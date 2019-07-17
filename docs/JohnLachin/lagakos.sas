*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/datasets/lagakos.sas                  *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: reads the data from Lagakos (1978) and creates a SAS      *;
*          data set that is used for the analyses in Chapter 9.      *;
*          This job should be run on your platform to create the     *;
*          SAS data set. The data set was originally used by         *;
*          Lagakos to describe an approach to the analysis of        *;
*          competing risks, there being two modes or causes of       *;
*          failure (spread of disease) - metastatic versus not. For  *;
*          the analyses herein, however, a single outcome is         *;
*          employed - spread of disease of any cause.                *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
libname bioweb '/jml/biostatmethods/datasets';
 
** 194 patients with squamous cell carcinoma:
  perfstat:0=ambulatory, 1=non-ambulatory performance status
  treatmnt: 0=A, 1=B
  age in years
  time in weeks
  cause of failure: 0=censored, 1=local spread of disease, 2=metastatic spread,
  fail (defined below): spread of disease, either local or metastatic
;
 
data inone; input perfstat treatmnt age time cause @@;
cards;
1  0  61  11  1  0  1  69  34  2  1  1  60  13  2  0  1  74  25  1
1  1  52  11  2  0  1  60  22  1  1  0  76   1  0  1  1  76   7  0
1  0  49   5  2  0  1  52  39  1  0  1  54  18  0  0  0  52  12  0
0  1  58  24  2  0  1  55   3  2  1  1  59  12  0  0  1  39  18  1
0  1  43   7  1  1  1  68  12  2  0  0  50  30  1  0  1  71  35  1
0  0  65  17  2  1  1  64  10  0  0  1  51   6  2  1  1  61   2  1
0  0  68  37  0  0  1  66  61  0  1  0  65  88  1  1  0  42  12  1
0  1  47  59  1  1  1  48  11  1  0  1  69  14  0  1  1  66   3  2
1  1  67  22  0  1  0  57  30  0  1  1  63   4  0  1  1  66  11  0
0  0  44  13  1  0  0  48  10  2  1  1  57   2  0  0  1  68  31  0
0  1  73  41  2  0  1  67  24  0  1  1  51  21  0  1  1  71  55  2
1  1  49  11  1  0  1  63  10  0  0  1  43  19  1  0  1  61  14  1
1  1  46   8  1  0  1  52  12  1  0  0  69  12  2  0  1  74  13  2
1  1  60   1  0  0  0  70  26  0  0  1  70  14  2  0  1  75  17  1
1  1  70  21  1  0  0  60  12  2  0  1  71  11  2  0  1  71  63  2
0  1  60  17  2  1  1  54  18  0  1  1  55  10  1  0  1  50   6  0
0  1  67  40  1  0  0  60  31  0  0  1  62  29  1  1  0  63  13  0
0  1  57  56  1  0  1  58   4  0  0  1  49   8  1  1  1  57  11  2
1  1  68  13  0  0  0  68  15  2  0  1  62   7  1  1  1  39   1  0
1  0  67   4  0  1  0  58  13  1  0  0  64  15  0  0  1  61  13  1
0  1  68  27  0  1  1  53   5  1  0  0  48  27  1  0  1  53   7  1
1  0  50   4  0  0  1  58  34  2  0  1  48   7  1  1  0  49   8  1
0  1  57  19  0  0  0  66  39  1  1  0  74  11  1  0  1  48  53  0
0  1  60  20  1  1  1  64   6  1  0  1  54  19  1  1  0  38  19  2
1  0  54   5  0  1  1  74  20  2  0  1  39  27  1  0  1  55  15  1
0  1  65  14  0  1  1  62  19  1  0  1  59   7  1  0  1  62  38  2
1  1  70  17  0  1  1  63  12  0  0  1  60   4  1  1  1  57   2  0
1  1  62  30  1  0  1  49 101  0  0  0  52  21  1  1  1  61   3  0
1  0  55   1  0  0  1  57  11  2  0  1  59   9  2  1  1  63  26  1
1  1  62  21  2  1  1  52   4  1  0  1  52  40  1  0  1  59  15  1
1  1  48   8  1  1  1  64  17  2  0  1  66   6  0  1  1  46  11  2
0  1  53   5  1  0  1  79  36  1  0  1  59  27  0  0  1  64  10  1
1  0  60  14  0  0  1  58  25  2  1  1  65   8  1  1  0  67   4  0
0  1  65  11  2  0  0  57  34  0  1  1  41   8  1  1  1  71   5  1
0  1  44  11  1  0  0  48  12  0  0  0  71  34  1  0  1  57  11  0
1  1  64   8  2  1  1  73  21  1  1  0  62  12  1  0  1  52  20  0
1  1  43   1  0  1  1  60   5  1  0  1  54  78  1  1  1  75  18  0
0  1  65   9  2  1  1  43   9  1  0  1  57  11  1  0  1  44  10  0
0  1  60  14  0  0  1  64  62  1  0  1  72  41  0  0  0  49  11  0
0  1  55  24  2  0  1  59   9  1  0  0  63   6  0  0  0  38   7  0
0  1  42  33  0  1  1  70   1  1  0  1  74  14  1  0  1  69  15  1
1  0  46  84  2  1  1  59   2  0  1  1  58   2  2  0  0  71  32  2
1  1  67  16  2  1  1  60   7  1  1  1  58   3  1  0  1  64  34  0
0  1  59  28  0  1  0  71  21  1  0  1  49  21  1  1  1  68  22  1
1  1  62  29  2  1  1  44  28  0  0  1  58  45  0  0  1  64  18  1
1  1  66  22  1  1  0  48  14  0  0  1  59   9  1  1  1  60   6  0
1  1  50  15  0  0  1  59  41  2  0  0  73   5  1  1  1  70   5  2
1  1  56   9  1  1  0  76  34  1  1  0  61   4  0
0  1  39  22  2  0  0  44  38  1  0  0  60  20  1
;
 
***proc print;
 
data intwo; set inone;
fail = cause; if cause=2 then fail=1;
patient=_n_;
 
label perfstat = 'ambulatory (1=no, 0=yes)'
  treatmnt = 'treatment group (0=A, 1=B)'
  age= 'age in years'
  time= 'time in weeks'
  cause= 'failure type (0=censored, 1=local spread, 2=metastatic)'
  fail = 'spread of disease, either local or metastatic (0=censored, 1=spread)';
 
data bioweb.lagakos; set intwo;
proc freq; tables perfstat*treatmnt*fail;
title1 'Lagakos squamous cell carcinoma data';
 
run;
