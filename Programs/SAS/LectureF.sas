/*Using (@) Line hold Identifier*/
data red_team;
input Team $ 13-18 @; 
if Team='red'; 
input IdNumber 1-4 StartWeight 20-22 EndWeight 24-26; 
datalines;
1023 David  red    189 165
1049 Amelia yellow 145 124
1219 Alan   red    210 192
1246 Ravi   yellow 194 177
1078 Ashley red    127 118
1221 Jim    yellow 220 .
; 
proc print data=red_team;
title 'RED TEAM';
run;

/*Using (@@) Line hold Identifier*/
data body_fat;
input Gender $ PercentFat @@;
datalines;
m 13.3 f 22
m 22 f 23.2
m 16 m 12
;
proc print data=body_fat;
run;

/*Combine datasets using one to one match merging*/
DATA descriptions;
INFILE '/home/u18376423/Rakesh/chocolate.dat' TRUNCOVER;
INPUT CodeNum $ 1-4 Name $ 6-14 Description $ 15-60;
RUN;
DATA sales;
INFILE '/home/u18376423/Rakesh/chocsales.dat';
INPUT CodeNum $ 1-4 PiecesSold 6-7;

PROC PRINT DATA = descriptions;
RUN;
PROC PRINT DATA = sales;
RUN;

PROC SORT DATA = sales;
BY CodeNum;
RUN;
* Merge data sets by CodeNum;
DATA chocolates;
MERGE sales descriptions;
BY CodeNum;
PROC PRINT DATA = chocolates;
TITLE "Today's Chocolate Sales";
RUN;

/* Cluster Analysis in SAS*/
/*The PROC VARCLUS procedure in SAS/STAT performs clustering of variables, 
it divides a set of variables by hierarchical clustering. */

proc print data=SASHELP.IRIS;
run;

proc varclus data=SASHELP.IRIS MAXCLUSTERS=4;
            var PetalWidth SepalWidth;
run;

/*The FASTCLUS procedure cluster analysis procedure performs k-means clustering 
on the basis of distances computed from one or more variables*/

proc print data=sashelp.cars;
run;

proc fastclus data=sashelp.cars maxclusters=20;
            var EngineSize Cylinders;

