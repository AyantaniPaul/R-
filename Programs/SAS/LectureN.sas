/* Modifying SAS Data Sets
The examples use the SAS data set INVENTORY as input. */

options pagesize=60 linesize=80 pageno=1 nodate;
data inventory;
input PartNumber $ Description $ InStock @17
Received Date date9. @27 Price;
format ReceivedDate date9.;
datalines;
K89R seal   34  27jul1998 245.00
M4J7 sander 98  20jun1998 45.88
LK43 filter 121 19may1999 10.99
MN21 brace 43   10aug1999 27.87
BC85 clamp 80   16aug1999 9.55
NCF3 valve 198  20mar1999 24.50
KJ66 cutter 6   18jun1999 19.77
UYN7 rod 211    09sep1999 11.55
JD03 switch 383 09jan2000 13.99
BV1E timer 26   03aug2000 34.50
;
proc print data=inventory;
title 'Tool Warehouse Inventory';
run;

/*You can use the MODIFY statement to replace all values for a specific variable or
variables in a data set. In the following program, the price of each part in the
inventory is increased by 15%. */

data inventory;
modify inventory;
price=price+(price*.15);
run;
proc print data=inventory;
title 'Tool Warehouse Inventory';
title2 '(Price reflects 15% increase)';
format price 8.2;
run;