*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/chapter5/kpairtab.sas                 *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: Analysis of the data in Example 5.12 using the macro      *;
*          %Kpaired for the stratified-adjusted analysis of K        *;
*          matched or paired 2x2 tables.                             *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*--------------------------------------------------------------------*;
 
*** The macro can be accessed either by "submitting" the job to sas before submitting this
    job, or by statements such as the following (remove the * first);
 
* filename Kpaired '/jml/biostatmethods/chapter5\kpair2x2.sas';
* %include Kpaired;
 
*** To save paper, I recommend that you use the following;
 
options formdlim='-';
 
DATA ONE;
INPUT K e f g h;
CARDS;
 1 10  6  5 11
 2 10 10  5 11
 3 10  6  1 11
 4 10  9  3 11
*****;
TITLE1 'Example 5.12 -- DCCT pregnancy data, *** e and h are dummy numbers, ignore RR';
%KPAIRED;
 
run;
