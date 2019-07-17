*--------------------------------------------------------------------*;
* Program: /jml/biostatmethods/datasets/hypoglycemia/impthypo.sas    *;
* Source:  Biostatistical Methods: The assessment of Relative Risks  *;
*          John Wiley and Sons, 2000                                 *;
* Author:  John M. Lachin <jml@biostat.bsc.gwu.edu>                  *;
*          Copyright (c) 2000 by John M. Lachin                      *;
* Purpose: This program reads the SAS export data set dccthypo.xpt   *;
*          and generates four SAS data sets on your platform:        *;
*          hyevents, hytimes, hypomimi and hypomimc that contain     *;
*          the recurrent  hypoglycemia events from the DCCT          *;
* Web URL: http://www.bsc.gwu.edu/jml/biostatmethods                 *;
*---------- ----------------------------------------------------------*;
 
** destination library;
 
    libname dcctdata 'c:\sasjobs\bioweb\datasets\hypoglycemia';
 
** input export data file;
 
    libname dcct xport 'c:\sasjobs\bioweb\datasets\hypoglycemia\dccthypo.xpt';
 
proc copy in=dcct out=dcctdata; run;
 
proc contents data=dcctdata.hyevents;
proc contents data=dcctdata.hytimes;
proc contents data=dcctdata.hypomimi;
proc contents data=dcctdata.hypomimc;
 
 run;
 
run;
