/*Combining a Grand Total with the Original Data*/
DATA shoes;
INFILE '/home/u18376423/Rakesh/Shoesales.dat';
INPUT Style $ 1-15 ExerciseType $ Sales;
RUN;
* Output grand total of sales to a data set and print;
PROC MEANS NOPRINT DATA = shoes;
VAR Sales;
OUTPUT OUT = summarydata SUM(Sales) = GrandTotal;
PROC PRINT DATA = summarydata;
TITLE 'Summary Data Set';
RUN;
* Combine the grand total with the original data;
DATA shoesummary;
IF _N_ = 1 THEN SET summarydata;
SET shoes;
Percent = Sales / GrandTotal * 100;
PROC PRINT DATA = shoesummary;
VAR Style ExerciseType Sales GrandTotal Percent;
TITLE 'Overall Sales Share';
RUN;

/*Concatenating SAS data sets*/
options pagesize=60 linesize=80 pageno=1 nodate;
data sales;
input EmployeeID $ 1-9 Name $ 11-29 @30 HireDate date9.
Salary HomePhone $;
format HireDate date9.;
datalines;
429685482 Martin, Virginia   09aug1990 34800 493-0824
244967839 Singleton, MaryAnn 24apr1995 27900 929-2623
996740216 Leighton, Maurice  16dec1993 32600 933-6908
675443925 Freuler, Carl      15feb1998 29900 493-3993
845729308 Cage, Merce        19oct1992 39800 286-0519
;
proc print data=sales;
title 'Sales Department Employees';
run;
data customer_support;
input EmployeeID $ 1-9 Name $ 11-29 @30 HireDate date9.
Salary HomePhone $;
format HireDate date9.;
datalines;
324987451 Sayre, Jay         15nov1994 44800 933-2998
596771321 Tolson, Andrew     18mar1998 41200 929-4800
477562122 Jensen, Helga      01feb1991 47400 286-2816
894724859 Kulenic, Marie     24jun1993 41400 493-1472
988427431 Zweerink, Anna     07jul1995 43700 929-3885
;
proc print data=customer_support;
title 'Customer Support Department Employees';
run;

data dept1_2;
set sales customer_support;
run;
proc print data=dept1_2;
title 'Employees in Sales and Customer Support Departments';
run;

/*Concatenating SAS data sets using append procedure*/
options pagesize=60 linesize=80 pageno=1 nodate;
proc append base=sales data=customer_support;
run;
proc print data=sales;
title 'Employees in Sales and Customer Support Departments';
run;

/*Example*/
data one;
input EmployeeID $ 1-9 Name $ 11-29 Age;
datalines;
12345     xyx                10
12346     xyz                10
;

data two;
input EmployeeID $ 1-9 Name $ 11-29 Age;
datalines;
12349     abc                11
12348     uvw                13
;

data three;
set one two;
run;

proc append base=one data=two;
run;

proc print data=one;
run;
