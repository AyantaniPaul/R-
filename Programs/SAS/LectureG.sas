/*One To Many Match Merge
Sometimes we need to combine two data sets by matching 
one observation from one data set with more than one 
observation in another.Before you merge two data sets, 
they must be sorted by one or more common variables.
If your data sets are not already sorted in the proper order, 
then use PROC SORT to do the job. 
*/
DATA regular;
INFILE '/home/u18376423/Rakesh/Shoes.dat';
INPUT Style $ 1-15 ExerciseType $ RegularPrice;
PROC SORT DATA = regular;
BY ExerciseType;
RUN;
DATA discount;
INFILE '/home/u18376423/Rakesh/Disc.dat';
INPUT ExerciseType $ Adjustment;
RUN;
* Perform many-to-one match merge;
DATA prices;
MERGE regular discount;
BY ExerciseType;
NewPrice = ROUND(RegularPrice - (RegularPrice * Adjustment), .01);
PROC PRINT DATA = prices;
TITLE 'Price List for May';
RUN;

/*Merge Summary Statistics with original Data
Once in a while you need to combine summary statistics 
with your data, such as when you want to compare each 
observation to the group mean, or when you want to 
calculate a percentage using the group total. To do this, 
summarize your data using PROC MEANS, and put the 
results in a new data set. Then merge the summarized data 
back with the original data using a one-to-many match 
merge. */

DATA shoes;
INFILE '/home/u18376423/Rakesh/Shoesales.dat';
INPUT Style $ 1-15 ExerciseType $ Sales;
PROC SORT DATA = shoes;
BY ExerciseType;
RUN;
* Summarize sales by ExerciseType and print;
PROC MEANS NOPRINT DATA = shoes;
VAR Sales;
BY ExerciseType;
OUTPUT OUT = summarydata SUM(Sales) = Total;
PROC PRINT DATA = summarydata;
TITLE 'Summary Data Set';
RUN;
* Merge totals with the original data set;
DATA shoesummary;
MERGE shoes summarydata;
BY ExerciseType;
Percent = Sales / Total * 100;

PROC PRINT DATA = shoesummary;
BY ExerciseType;
ID ExerciseType;
VAR Style Sales Total Percent;
TITLE 'Sales Share by Type of Exercise';
RUN;
