/* The simplest form of the INPUT statement uses list input. List input is used to read
data values that are separated by a delimiter character. With list input, SAS reads a data
value until it encounters a blank space. SAS assumes the value has ended and assigns the 
data to the appropriate variable in the program data vector.*/
data scores;
   input name $ score1 score2;
   datalines;
Riley 1132 1187
Henderson 1015 1102
;
PROC PRINT data = scores;
title 'Scores';
run;

/*This program also uses the above data but notice that here the data
is delimited by a comma instead of a blank space, the default delimiter.*/

data scores;
infile datalines dlm=',';
   input name $ score1 score2;
   datalines;
Riley, 1132, 1187
Henderson, 1015, 1102
;
PROC PRINT data = scores;
title 'Scores';
run;

/*Modified List Input,  Use of ampersand(&) modifier. It allows you to read character
 values that contain embdded blanks.*/

DATA citypops1;
   infile DATALINES FIRSTOBS = 2;
   length city $ 12;
   input city & pop2000;
   DATALINES;
City  Yr2000Popn
New York  8008278
Los Angeles  3694820
Chicago  2896016
Houston  1953631
Philadelphia  1517550
Phoenix  1321045
San Antonio  1144646
San Diego  1223400
Dallas  1188580
San Jose  894943
;
RUN;

PROC PRINT data = citypops1;
   title 'The citypops data set1';
   format pop2000 comma10.;
RUN;

/*Modified List Input,  Use of Colon(:) modifier. It allows you to read non standard datavalues 
and character values that are longer than eight character. 
FIRSTOBS=2 option in the infile statement tells SAS to begin reading 
data at line-2 instead of line-1*/

DATA citypops;
   infile DATALINES FIRSTOBS = 2;
   input city & $12. pop2000: comma.;
   DATALINES;
City  Yr2000Popn
New York  8,008,278
Los Angeles  3,694,820
Chicago  2,896,016
Houston  1,953,631
Philadelphia  1,517,550
Phoenix  1,321,045
San Antonio  1,144,646
San Diego  1,223,400
Dallas  1,188,580
San Jose  894,943
;
RUN;

/*A format statement above and below tells SAS to format of input data and format
 to printing procedure in SAS for variable pop2000. */

PROC PRINT data = citypops;
   title 'The citypops data set';
   format pop2000 comma10.;
RUN;

/* Column Input-This input method allows you to read variable value that 
occupied the same column within each record. TRUNCOVER option forces the 
input statement to stop reading when it gets to end of a short line 
This option will not skip the information*/

data scores;
   infile datalines truncover;
   input name $ 1-12 score2 17-20 score1 27-30; /*Column Input*/
   datalines;
Riley           1132       987 
Henderson       1015      1102
;

proc print data=scores;
run;


data hw;
   infile datalines truncover;
   input height 1-4 weight 5-8; /*Column Input*/
   datalines;
170 67 
173 68
160 57  
;

proc means data=hw;
run;

/*Formatted Input allows us to read both standard and nonstandard numeric data*/
/*+n relative pointer control tells SAS to move the input pointer forward n columns
 to a column number that is relative to the current position.
 comma w.d conatians commas and decimal places*/
data scores;
   input name $12. +4 score1 comma5. +6 score2 comma5.; /*Formatted Input */
   datalines;
Riley           1,132      1,187
Henderson       1,015      1,102
;

proc print data=scores;
run;

/*Named Input- Named input reads the input data records that contain a variable 
name followed by an equal sign and a value for the variable. The INPUT
 statement reads the input data record at the current location of the input pointer*/
data games;
   input name=$ score1= score2=; 
   datalines;
name=riley score1=1132 score2=1187  
;

proc print data=games;
run;


/*The MEANS and SUMMARY Procedures*/


data my_data;
    input Method $ Score @@;
    datalines;
A 78 A 81
A 82 A 82
A 85 A 88
A 88 A 90
B 81 B 83
B 83 B 85
B 86 B 88
B 90 B 91
C 84 C 88
C 88 C 89
C 90 C 93
C 95 C 98
;
run;

proc means data=my_data;
class Method;
var Score;
run;

proc means data=my_data median;
class Method;
var Score;
run;

proc summary data=my_data PRINT;
class Method;
variable Score;
run;
