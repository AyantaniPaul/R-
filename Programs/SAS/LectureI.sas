LIBNAME stat480 '/home/u18376423/Rakesh';  *Specifies the SAS data
                                      library (directory);
DATA stat480.temp2;
  input subj 1-4 gender 6 height 8-9 weight 11-13;
  DATALINES;
1024 1 65 125
1167 1 68 140
1168 2 68 190
1201 2 72 190
1302 1 63 115
  ;
RUN;

PROC PRINT data=stat480.temp2;
  title 'Output dataset: STAT480.TEMP2';
RUN;

options yearcutoff=1950;
data _null_;
   a='26oct02'd;
   put 'SAS date='a date.;
   put 'formatted date='a date9.;
run;

options yearcutoff=1920 nodate pageno=1 linesize=80 pagesize=60;

data schedule;
   input @1 jobid $ @6 projdate mmddyy10.;
   datalines;
A100 01/15/25
A110 03/15/2025
A200 01/30/96
B100 02/05/00
B200 06/15/2000
;

proc print data=schedule;
   format projdate mmddyy10.;
run;

data test;
   Time1=86399;
   format Time1 datetime.;
   Date1=86399;
   format Date1 date.;
   Time2=86399;
   format Time2 timeampm.;
run;
proc print data=test;
   title  'Same Number, Different SAS Values';
   footnote1 'Time1 is a SAS DATETIME value';
   footnote2 'Date1 is a SAS DATE value';
   footnote3 'Time2 is a SAS TIME value'.;
run;
footnote;


DATA rakesh.icecream;
INFILE '/home/u18376423/Rakesh/IceCreamSales.dat' FIRSTOBS = 3;
INPUT Flavor $ 1-9 Location BoxesSold;
RUN;

Proc print data=icecream;
run;

/*Export Procedure in SAS*/
PROC EXPORT DATA=icecream
            OUTFILE='/home/u18376423/Rakesh/IceCreamSales.csv'
            DBMS=DLM REPLACE;
     DELIMITER=',';
     PUTNAMES=NO;
RUN;

DATA _NULL_;
SET icecream;
FILE '/home/u18376423/Rakesh/Newfile.dat';
PUT Flavor $ 1-9 Location BoxesSold;
RUN;
