/*Creating a new transformed variable using if then statement*/
  DATA grades;
	input name $ 1-15 e1 e2 e3 e4 p1 f1;
	* if the first exam is less than 65 indicate failed;
	if (e1 < 65) then status = 'Failed';
	DATALINES;
Alexander Smith  78 82 86 69  97 80
John Simon       88 72 86  . 100 85
Patricia Jones   98 92 92 99  99 93
Jack Benedict    54 63 71 49  82 69
Rene Porter     100 62 88 74  98 92
;
RUN;

PROC PRINT data = grades;
	var name e1 status;
RUN;
/*=======================================================================*/
/*Creating a new transformed variable using if then else statement*/
DATA grades;
	input name $ 1-15 e1 e2 e3 e4 p1 f1;
	* if the first exam is less than 65 indicate failed;
	if (e1 < 65) then status = 'Failed';
	* otherwise indicate passed;
	else status = 'Passed';
	DATALINES;
Alexander Smith  78 82 86 69  97 80
John Simon       88 72 86  . 100 85
Patricia Jones   98 92 92 99  99 93
Jack Benedict    54 63 71 49  82 69
Rene Porter     100 62 88 74  98 92
;
RUN;

PROC PRINT data = grades;
	var name e1 status;
RUN;
/*=======================================================================*/
DATA grades;
	input name $ 1-15 e1 e2 e3 e4 p1 f1;
	if p1 in (98, 99, 100) then project = 'Excellent';
	DATALINES;
Alexander Smith  78 82 86 69  97 80
John Simon       88 72 86  . 100 85
Patricia Jones   98 92 92 99  99 93
Jack Benedict    54 63 71 49  82 69
Rene Porter     100 62 88 74  98 92
;
RUN;

PROC PRINT data = grades;
	var name p1 project;
RUN;
/*=======================================================================*/
DATA grades;
    length overall $ 10;
   	input name $ 1-15 e1 e2 e3 e4 p1 f1;
	avg = round((e1+e2+e3+e4)/4,0.1);   
		 if (avg = .)                   then overall = 'Incomplete';
	else if (avg >= 90)                 then overall = 'A';
	else if (avg >= 80) and (avg < 90)  then overall = 'B';
	else if (avg >= 70) and (avg < 80)  then overall = 'C';
	else if (avg >= 65) and (avg < 70)  then overall = 'D';
	else if (avg < 65)                  then overall = 'F';	
	DATALINES;
Alexander Smith  78 82 86 69  97 80
John Simon       88 72 86  . 100 85
Patricia Jones   98 92 92 99  99 93
Jack Benedict    54 63 71 49  82 69
Rene Porter     100 62 88 74  98 92
;
RUN;

PROC PRINT data = grades;
	var name avg overall;
RUN;
/*==================================================================*/
DATA grades;
    length overall $ 10;
   	input name $ 1-15 e1 e2 e3 e4 p1 f1;
	avg = round((e1+e2+e3+e4)/4,0.1);   
		 if (avg EQ .)         then overall = 'Incomplete';
	else if (90 LE avg LE 100) then overall = 'A';
	else if (80 LE avg LT  90) then overall = 'B';
	else if (70 LE avg LT  80) then overall = 'C';
	else if (65 LE avg LT  70) then overall = 'D';
	else if (0  LE avg LT  65) then overall = 'F';	
	DATALINES;
Alexander Smith  78 82 86 69  97 80
John Simon       88 72 86  . 100 85
Patricia Jones   98 92 92 99  99 93
Jack Benedict    54 63 71 49  82 69
Rene Porter     100 62 88 74  98 92
;
RUN;

PROC PRINT data = grades;
	var name avg overall;
RUN;
/*======================================================================*/
DATA grades;
    length action $ 7 
           action2 $ 7; 
    input name $ 1-15 e1 e2 e3 e4 p1 f1 status $;
	     if (status = 'passed') then action = 'none';
    else if (status = 'failed') then action = 'contact';
	else if (status = 'incomp') then action = 'contact';
	     if (upcase(status) = 'PASSED') then action2 = 'none';
	else if (upcase(status) = 'FAILED') then action2 = 'contact';
	else if (upcase(status) = 'INCOMP') then action2 = 'contact';
	DATALINES;
Alexander Smith  78 82 86 69  97 80 passed
John Simon       88 72 86  . 100 85 incomp
Patricia Jones   98 92 92 99  99 93 PAssed
Jack Benedict    54 63 71 49  82 69 FAILED
Rene Porter     100 62 88 74  98 92 PASSED
;
RUN;

PROC PRINT data = grades;
	var name status action action2;
RUN;
/*The IFC and IFN functions in SAS are used to create a new 
data column based on a conditional expression. They can be 
used instead of IF-THEN-ELSE statements. The difference 
between IFC and IFN is that the values returned from IFC
are character, while they are numeric for IFN*/


data _null_;
   input name $ grade;
   performance = ifc(grade>80, 'Pass             ', 'Needs Improvement');
   put name= performance=;
   datalines;
John 74
Kareem 89
Kati 100
Maria 92
;

run;

/*======================================*/
data _null_;
   input TotalSales;
   commission=ifn(TotalSales > 10000, TotalSales*.05, TotalSales*.02);
   put commission=;
   datalines;
25000
10000
500
10300
;
run;
