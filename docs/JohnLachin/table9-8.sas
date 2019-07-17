*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter9/Table9-8.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: This program reads the data in Table 9.8 in a format      *;
*          suitable for use with the other macros provided or with   *;
*          SAS procedures.                                           *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*---------------------------------------------------------------------*;
 
*** note: "Not Healed" changed to "NotHeal" in the status variable;
*** group 1 = drug, 2 = placebo;
 
data table98; input group status $ month lastxray;
cards;
 1   Healed       2   .
 1   Healed       2   .
 1   Healed       2   .
 1   Healed       2   .
 1   Healed       2   .
 1   Healed       2   .
 1   Lost         1   0
 1   Lost         2   0
 1   Surgery      1   0
 1   Healed       4   0
 1   Healed       4   .
 1   Healed       4   .
 1   Healed       4   .
 1   Death        3   2
 1   Lost         3   2
 1   Lost         4   2
 1   Lost         4   2
 1   Surgery      3   2
 1   Surgery      4   2
 1   Healed       6   .
 1   Healed       6   .
 1   Lost         5   4
 1   Lost         5   4
 1   Surgery      6   4
 1   NotHeal      6   6
 1   NotHeal      6   6
 1   NotHeal      6   6
 1   NotHeal      6   6
 2   Healed       2   .
 2   Healed       2   .
 2   Lost         1   0
 2   Lost         1   0
 2   Lost         2   0
 2   Surgery      2   0
 2   Surgery      2   0
 2   Death        1   0
 2   Healed       4   .
 2   Healed       4   .
 2   Lost         3   2
 2   Lost         4   2
 2   Surgery      4   2
 2   Surgery      3   2
 2   Surgery      4   2
 2   Healed       6   .
 2   Healed       6   .
 2   Lost         5   4
 2   Lost         5   4
 2   Lost         5   4
 2   Surgery      6   4
 2   Surgery      5   4
 2   Death        5   4
 2   NotHeal      6   6
;
*  add additional statements;
