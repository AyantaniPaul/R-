/*Controlling Input with options in the infile Statement*/
/*If your raw data is already stored in a file, then you do not have to
bring that file into data stream. You can use infile statement to specify the file 
containing raw data*/
DATA icecream;
INFILE '/home/u18376423/Rakesh/IceCreamSales.dat' FIRSTOBS = 3;
INPUT Flavor $ 1-9 Location BoxesSold;
RUN;

proc print data=icecream;
run;

DATA icecream1;
INFILE '/home/u18376423/Rakesh/IceCreamSales.dat' FIRSTOBS = 3 OBS=4;
INPUT Flavor $ 1-9 Location BoxesSold;
RUN;

proc print data=icecream1;
run;
/*You can also use the alternate method to referencing the file as given
 below in the program*/

filename rr '/home/u18376423/Rakesh/IceCreamSales.dat';
DATA icecream2;
INFILE  rr FIRSTOBS = 3 OBS=4;
INPUT Flavor $ 1-9 Location BoxesSold;
RUN;

proc print data=icecream2;
run;

/* FLOWOVER is default behaviour of the infile statement. Input pointer moves to
the next record to fulfil the number of variables used in 
the input statement if some of the variables are missing in the first line.*/

/*STOPOVER stops the data step processing if any of the 
variable mentioned in input statement is missing.
*/

data weather;
   infile datalines stopover;
   input temp1-temp5;
   datalines;
97.9 98.1 98.3
98.6 99.2 99.1 98.5 97.5
96.2 97.3 98.3 97.6 96.5
;
/*Because SAS does not find a TEMP4 value in the first data record, it sets _ERROR_ to 1, 
stops building the data set, and prints data line 1.*/

proc print data=weather;
run;


/* MISSOVER option prevents the DATA step from going to the next line if it
 does not find values in the current record for all of the variables in the 
INPUT statement. Instead, the DATA step assigns a missing value for all variables
 that do not have values.*/
DATA class102;
INFILE '/home/u18376423/Rakesh/AllScores.dat' MISSOVER;
INPUT Name $ Test1 Test2 Test3 Test4 Test5;
RUN;

proc print data=class102;
run;

/*TRUNCOVER option causes the DATA step to assign the raw data value to the variable 
even if the value is shorter than expected by the INPUT statement. If, when the DATA step
encounters the end of an input record, there are variables without values, the
variables are assigned missing values for that observation.*/

DATA homeaddress;
INFILE '/home/u18376423/Rakesh/Address.dat' TRUNCOVER;
INPUT Name $ 1-15 Number 16-19 Street $ 22-37;
RUN;

proc print data=homeaddress;
run;


/*Regression procedure in SAS can be used to perform regression analysis in SAS. 
An example is given below*/

 data regan;
  input  Height Weight Age;
      datalines;
  69.0 112.5 14    
  62.8 102.5 14    
  59.8  84.5 12  
  59.0  99.5 12     
  56.3  77.0 12   
  64.8 128.0 12   
  66.5 112.0 15
  56.5  84.0 13
  65.3  98.0 13
  57.3  83.0 12
  63.5 102.5 14
  62.5 112.5 15
  62.5  84.0 13
  51.3  50.5 11
  64.3  90.0 14
  66.5 112.0 15
  72.0 150.0 16
  67.0 133.0 15 
  57.5  85.0 11
   ;
   
 /*ODS Statistical Graphics (also known as ODS Graphics) is functionality for easily creating statistical graphics. 
  Over 100 statistical procedures can produce graphs as automatically as they do tables*/
 
  ods graphics on;
   
   proc reg;
      model Weight = Height;
   run;
   
  ods graphics off;
  
 
 /*An example is given below to perform regression for more than one regressor variable*/
  
   ods graphics on;
   
   proc reg;
      model Weight = Height Age;
   run;
   
  ods graphics off;
