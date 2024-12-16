/*#n line pointer Control*/
data club2;
input #2 Team $ 1-6 #1 Name $ 6-23 IdNumber 1-4
#3 StartWeight 1-3 EndWeight 5-7;
datalines;
1023 David Shaw
red
189 165
1049 Amelia Serrano
yellow
145 124
1219 Alan Nance
red
210 192
1246 Ravi Sinha
yellow
194 177
1078 Ashley McKnight
red
127 118
1221 Jim Brown
yellow
220 .
;
proc print data=club2;
run;

/* @ line pointer Control a data record is read by more than one INPUT statement (trailing @).
 Use a single trailing @ to allow the next INPUT statement to read from the same record.*/
data test;
   input a @1 b;
   datalines;
 4
 5
;

proc print data=test;
run;

DATA temp;
  input @1  subj 4. 
        @6  f_name $11. 
		@18 l_name $6.
		+3 height 2. 
        +5 wt_date mmddyy8. 
        +1 calorie comma5.;
  format wt_date mmddyy8. calorie comma5.;
  DATALINES;
1024 Alice       Smith  1 65 125 12/1/95  2,036
1167 Maryann     White  1 68 140 12/01/95 1,800
1168 Thomas      Jones  2    190 12/2/95  2,302
1201 Benedictine Arnold 2 68 190 11/30/95 2,432
1302 Felicia     Ho     1 63 115 1/1/96   1,972
  ;
RUN;

PROC PRINT data = temp;
  title 'Output dataset: TEMP';
  id subj;
RUN;

/*Dataset procedure in SAS is used to manipulate SAS datasets.
 It can be used to perform a variety of operations on datasets, 
 such as copying, renaming, deleting, and modifying datasets.*/

proc datasets library=work;
contents data=club2;
run;

proc datasets library=sashelp nolist;
contents data=air;
run;
/* It can also be used to copy a dataset from one library to another */

proc datasets;
copy in=work out=rakesh;
run;


/*It can also be used to copy a dataset from one library to another*/
proc datasets;
copy in=work out=rakesh;
select temp;
run;


/* Logistic Regression using SAS logistic procedure. To model 1s rather than 0s, we use the descending option.
 The default proc logistic models 0s rather than 1s*/

proc means data="/home/u18376423/Rakesh/SAS_Programs/binary.sas7bdat";
  var gre gpa;
run;

proc freq data="/home/u18376423/Rakesh/SAS_Programs/binary.sas7bdat";
  tables rank admit admit*rank;
run;

proc logistic data="/home/u18376423/Rakesh/SAS_Programs/binary.sas7bdat" descending;
  class rank / param=ref;
  model admit = gre gpa rank;
run;

/*The CLASS statement names the classification variables to be used in the analysis.*/

