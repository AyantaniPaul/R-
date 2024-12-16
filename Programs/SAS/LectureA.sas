/*SAS programs are constructed from two basic building blocks: DATA steps and PROC steps
/*A typical program starts with a DATA step to create a SAS data set and then passes the 
  data to a PROC step for processing.*/
DATA distance;
input Miles;
Kilometers = 1.61 * Miles;
datalines;
12
23
24
24
;
PROC PRINT DATA = distance;
title 'Miles in Kms';
RUN;

/*The $ sign indicates that the variable is a string*/
/*The print procedure is used to print a dataset*/
/*We can also create a new variable from the existing one*/
/*In the program given below, we have created the variable LOSS*/
/*In the SAS program, stands for missing value*/
data club1;
input IdNumber Name $ Team $ StartWeight EndWeight;
Loss=StartWeight-EndWeight;
datalines;
1023 David red 189 165 
1049 Amelia yellow 145 124
1219 Alan red 210 192
1246 Ravi yellow 194 177
1078 Ashley red 127 118
1221 Jim yellow 220 . 
; 
proc print data=club1;
title "Weight of Club Members";
run;


DATA Employee; 
  INPUT empid ename $ salary DEPT $ ; 
DATALINES; 
1 Rick 623.3 	IT 		 
2 Dan 515.2 	OPS	
3 Mike 611.5 	IT 	
4 Ryan 729.1    HR 
5 Gary 843.25   FIN 
6 Tusar 578.6   IT 
7 Pranab 632.8  OPS
8 Rasmi 722.5   FIN 
;
RUN;
Proc print data=Employee;
run;
/*To read data from an existing SAS data set, you can use SET statement
The KEEP statement creates a new data set that contains only the variables 
listed in the KEEP statement*/

DATA OnlyDept;
 SET Employee;
 KEEP DEPT;
  RUN;
 PROC PRINT DATA=OnlyDept; 
RUN;
/*The following example uses the DROP statement to read all of the variables
from the CITY data set and to exclude the variables that are listed in the DROP
statement from being written to the new data set.*/
DATA OnlyDept;
 SET Employee;
 DROP empid salary;
  RUN;
 PROC PRINT DATA=OnlyDept; 
RUN; 

/*Often programmers find that they want to use some of the observations in a
data set and exclude the rest. The most common way to do this is with a
subsetting IF statement in a DATA step.*/

DATA OnlyDept;
 SET Employee;
 IF salary < 700 THEN DELETE;
  RUN;
 PROC PRINT DATA=OnlyDept; 
RUN; 

