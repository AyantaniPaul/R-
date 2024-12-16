/* This Lecture is about obtaining maximum likelihood estimate in SAS
In the program given below, we will obtain the MLE of Binomial distribution*/ 

data Binomial;
input N @@;
datalines;
   7 7 13 9 8 8 9 9 5 6 6 9 5 10 4 5 3 8 4
;

proc means data=Binomial;
run;
proc nlmixed data=Binomial;
   parms p = 0.5;             * initial value for parameter;
   bounds 0 < p < 1;
   NTrials = 20;
   model N ~ binomial(NTrials, p);
run;

/* Here is the program to fit Poisson distribution for the Given Dataset */

proc nlmixed data=Binomial;
   parms p = 0.5;             * initial value for parameter;
   bounds 0 < p < 100;
   NTrials = 20;
   model N ~ poisson(p);
run;

data gam;
input t @@;
datalines;
  3.0568005  3.4967798  5.0790670  0.6473754  5.7175502  3.6831836
  0.7183728  2.9193013  0.6285473  3.2015139  6.2094257 12.5801915
  0.7864980  7.6370416  1.3344010  4.2110403  7.3025858 12.1536682
  1.0991529 17.6610024
;

proc nlmixed data=gam;
   parms a=1.4 b=2.9;             * initial value for parameter;
   model t ~ gamma(a, b);
run;