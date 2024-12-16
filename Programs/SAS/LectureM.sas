DATA grades;
	input name $ 1-15 e1 e2 e3 e4 p1 f1;
	* add up each students four exam scores
	  and store it in examtotal;
	examtotal = e1 + e2 + e3 + e4;
	DATALINES;
Alexander Smith  78 82 86 69  97 80
John Simon       88 72 86  . 100 85
Patricia Jones   98 92 92 99  99 93
Jack Benedict    54 63 71 49  82 69
Rene Porter     100 62 88 74  98 92
;
RUN;

PROC PRINT data = grades;
	var name e1 e2 e3 e4 examtotal;
RUN;

DATA grades;
	input name $ 1-15 e1 e2 e3 e4 p1 f1;
	e2 = e2 + 8;  * add 8 to each student's 
	                second exam score (e2);
	DATALINES;
Alexander Smith  78 82 86 69  97 80
John Simon       88 72 86  . 100 85
Patricia Jones   98 92 92 99  99 93
Jack Benedict    54 63 71 49  82 69
Rene Porter     100 62 88 74  98 92
;
RUN;

PROC PRINT data = grades;
	var name e1 e2 e3 e4 p1 f1;
RUN;

DATA grades;
	input name $ 1-15 e1 e2 e3 e4 p1 f1;
	final = 0.6*e1+e2+e3+e4/4 + 0.2*p1 + 0.2*f1;
	DATALINES;
Alexander Smith  78 82 86 69  97 80
John Simon       88 72 86  . 100 85
Patricia Jones   98 92 92 99  99 93
Jack Benedict    54 63 71 49  82 69
Rene Porter     100 62 88 74  98 92
;
RUN;

PROC PRINT data = grades;
	var name e1 e2 e3 e4 p1 f1 final;
RUN;

DATA grades;
	input name $ 1-15 e1 e2 e3 e4 p1 f1;
	final = 0.6*((e1+e2+e3+e4)/4) + 0.2*p1 + 0.2*f1;
	DATALINES;
Alexander Smith  78 82 86 69  97 80
John Simon       88 72 86  . 100 85
Patricia Jones   98 92 92 99  99 93
Jack Benedict    54 63 71 49  82 69
Rene Porter     100 62 88 74  98 92
;
RUN;

PROC PRINT data = grades;
	var name e1 e2 e3 e4 p1 f1 final;
RUN;

DATA grades;
	input name $ 1-15 e1 e2 e3 e4 p1 f1;
	* calculate the average by definition;
	avg1 = (e1+e2+e3+e4)/4;   
	* calculate the average using the mean function;
	avg2 = mean(e1,e2,e3,e4); 
	DATALINES;
Alexander Smith  78 82 86 69  97 80
John Simon       88 72 86  . 100 85
Patricia Jones   98 92 92 99  99 93
Jack Benedict    54 63 71 49  82 69
Rene Porter     100 62 88 74  98 92
;
RUN;

PROC PRINT data = grades;
	var name e1 e2 e3 e4 avg1 avg2;
RUN;

DATA grades;
	input name $ 1-15 phone e1 e2 e3 e4 p1 f1;
	areacode = int(phone/10000000);
	DATALINES;
Alexander Smith 8145551212  78 82 86 69  97 80
John Simon      8145562314  88 72 86  . 100 85
Patricia Jones  7175559999  98 92 92 99  99 93
Jack Benedict   5705551111  54 63 71 49  82 69
Rene Porter     8145542323 100 62 88 74  98 92
;
RUN;

PROC PRINT data = grades;
	var name phone areacode;
RUN;
