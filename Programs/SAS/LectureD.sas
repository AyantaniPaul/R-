/* Testing the condition before creating an observation and use of @ Trailing*/
data red_team;
input team $ 13-18 @;
if team='red';
input Id 1-4 StartWeight 20-22 EndWeight 24-26;
datalines;
1001        red    50  45
1002        blue   55  47
1003        red    51  48
;
run;

proc print data=red_team;
title 'Red Team';
run;

/*Creating multiple observations from single record using double trailing @@ 
line hold specifier*/

data body_fat;
input gender $ Percentfat @@;
datalines;
m 13 f 14
m 22 f 32
m 16 f 17
;
run;

proc print data=body fat;
run;

/* Reading Multiple record to create a single observation*/
/*Mutipule Input statement one for each recod  can read each record
 into a single observation*/

data club2;
input IdNumber 1-4; 
input;    /*Null input Statement  or you can use input Team $ 1-6; */
input StartWeight 1-3 EndWeight 5-7; 
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
title 'Weight Club Members';
run;

/* Using  (/) Line Pointer Control: 
The forward slash (/) line pointer control is used to read records in a 
raw data file sequentially */

data club2;
input IdNumber 1-4 / / StartWeight 1-3 EndWeight 5-7;
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

/* The same task can be performed in SAS using multiple input statement.
An example is given below*/

data club3;
input Id 1-4;
input;
input stw 1-3 enw 5-7;
datalines;
1023 David Shaw
red
189 165
1049 Amelia Serrano
yellow
145 124
;
proc print data=club3;
title 'Weight Club Members';
run;


/* Analysis of Variance (ANOVA) in SAS*/

/*Suppose a researcher recruits 24 students to participate in a study.
The students are randomly assigned to use one of three studying methods
to prepare for an exam. The exam results for each student are shown in data*/

/*create dataset*/
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

proc ANOVA data=my_data;
class Method;
model Score = Method;
run;

/*We can use the means function to specify that a Tukey post-hoc test. This test
should be performed if the overall p-value of the one-way ANOVA is statistically
significant.*/

proc ANOVA data=my_data;
class Method;
model Score = Method;
means Method / tukey cldiff;
run;

/*Two way anova. There is one continuous dependent varialble retention and two
 categorical variables Fe and Zn. The two way anova can performed as follows*/

data retention;
	input retention 1-6 Fe $ 7-11 Zn $ 16-20  ;
	datalines;
0.2089	Low	    High
0.1563	Low	    High
0.1339	Low	    High
0.2407	Low	    High
0.2027	Low	    High
0.2429	Low	    High
0.2940	Low     High
0.2394	Low	    High
0.7255	Low	    Low
0.7285	Low	    Low
0.7375	Low	    Low
0.6855	Low	    Low
0.7093	Low	    Low
0.6996	Low	    Low
0.6619	Low	    Low
0.7052	Low	    Low
1.0635	High    Low
1.0532	High    Low
0.953	High    Low
0.9355	High    Low
0.9389	High    Low
0.9852	High    Low
0.9923	High    Low
1.0267	High    Low
1.4334	High    High
1.4961	High    High
1.3809	High    High
1.3718	High    High
1.4063	High    High
1.375	High    High
1.4286	High    High
1.4004	High    High
;
run;

proc anova data=retention;
      class Fe Zn;
	  model retention = Fe Zn Fe*Zn;
run;
