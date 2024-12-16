data demo;
 input patient 1-2 sex $ 3-4 age 5-7 ps 8-9;
 datalines;
 1 F 45 0
 4 M 63 2
 3 M 57 1
 5 F 72 3
 2 F 39 0
 3 M 57 1
 4 M 63 0
 ;
 run;
 proc sort data=demo out=demo1;
 by patient;
 run;

 proc sort data=demo out=demo2;
 by decending patient;
 run;

 proc sort data=demo out=demo3(keep=patient age);
 by patient;
 run; 

 proc sort data=demo(where=(age>50)) out=demo8;
 by patient;
 run; 

 proc sort data=demo out=demo9;
 where age>50 AND sex='M';
 by patient;
 run;


/*Format Procedure*/
/*Informats tell SAS how to read data, while formats
 tell SAS how to write (or print) data.*/
/*The FORMAT procedure enables you to define your own informats and formats for
 variables. In addition, you can perform these actions:
print the parts of a catalog that contain informats or formats
store descriptions of informats or formats in a SAS data set
use a SAS data set to create informats or formats*/

proc format;
value clsfmt 1='Freshman' 2='Sophmore' 3='Junior' 4='Senior';
run;
data class;
   format z clsfmt.;
   label x='ID NUMBER'
      y='AGE'
      z='CLASS STATUS';
   input x y z;
datalines;
1 20 4
2 18 1
;

proc print data=class;
run;

/*The CONTENTS procedure generates summary 
information about the contents of a dataset, including:
1.The variables' names, types, and attributes(including formats, informats, and labels)
2.How many observations are in the dataset
3.How many variables are in the dataset
4.When the dataset was created
This procedure is especially useful if you have imported your data from a file and want
to check that your variables have been read correctly, and have the appropriate variable
type and format.
*/
proc contents data=class;
run;
/*Using dataset procedure given below, all labels and formats are removed
 by using modify statement and attrib option. The contents statement in dataset procedure 
 displays the contents of the dataset*/
proc datasets lib=work memtype=data;
   modify class; 
     attrib _all_ label=' '; 
     attrib _all_ format=;
run;

/*The DATASETS procedure is a utility procedure that manages your SAS files. With PROC DATASETS, you can do the following:
1. copy SAS files from one SAS library to another
2. rename SAS files
3. repair SAS files
4. delete SAS files
5. list the SAS files that are contained in a SAS library
6. list the attributes of a SAS data set:
7. the date on which the data was last modified
8. whether the data is compressed
9. whether the data is indexed
10.manipulate passwords on SAS files
11.append SAS data sets
12.modify attributes of SAS data sets and variables within the data sets
13.create and delete indexes on SAS data sets
14.create and manage audit files for SAS data sets
15.create and delete integrity constraints on SAS data sets
16.create and manage extended attributes of data sets
*/
proc contents data=class;
run;
