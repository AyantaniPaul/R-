/*Root Finding by Newton Raphson Method*/
/*Code for finding square root of 5 using Newton Raphson Method*/
data root;
        x=5; 
        y0=1;       
        count=0;
        do until (w<1e-8);  
                count=count+1;  
                y=(y0+x/y0)/2; 
                w=abs(y-y0);
                y0=y;
                if y =sqrt(x) then is_eq_sqrt="YES";
                else is_eq_sqrt="NO";
                output;
        end;
run;

proc print data=root;
run;

/*Output from data step Apllication*/
data _null_;
set work.root;
file '/home/u18376423/Rakesh/rr1.dat';
put  x y0 count;
run;

/*Generation of Uniform Random Number*/

data A;
call streaminit(123); /* set random number seed */
do i = 1 to 10;
   u = rand("Uniform"); /* u ~ U[0,1] */
   output;
end;
run;

proc print data=A;
run;

/*Program for the Correlation Procedure*/

data bodyfat;
Input age pctfat;
datalines;
23 28
39 31
41 26
49 25
50 31
53 35
53 42
54 29
56 33
57 30
58 33
58 34
60 41
61 34
;
proc print data=bodyfat;
run;
proc corr data=bodyfat;
run;

proc means data=bodyfat;
run;

/* SAS Program for Performing ANOVA*/
title1 'Nitrogen Content of Red Clover Plants';
   data Clover;
      input Strain $ Nitrogen @@;
      datalines;
   3DOK1  19.4 3DOK1  32.6 3DOK1  27.0 3DOK1  32.1 3DOK1  33.0
   3DOK5  17.7 3DOK5  24.8 3DOK5  27.9 3DOK5  25.2 3DOK5  24.3
   3DOK4  17.0 3DOK4  19.4 3DOK4   9.1 3DOK4  11.9 3DOK4  15.8
   3DOK7  20.7 3DOK7  21.0 3DOK7  20.5 3DOK7  18.8 3DOK7  18.6
   3DOK13 14.3 3DOK13 14.4 3DOK13 11.8 3DOK13 11.6 3DOK13 14.2
   COMPOS 17.3 COMPOS 19.4 COMPOS 19.1 COMPOS 16.9 COMPOS 20.8
   ;
   proc anova data=clover;
    class strain;
      model Nitrogen = Strain;
   run;
   
/* SAS Program for Simple Random Sampling with and Without Replacement*/ 
data club;
input IdNumber 1-4 Name $ 6-24 Team $ StartWeight EndWeight;
Loss = StartWeight - EndWeight;
datalines;
1023 David Shaw         red 189 165
1049 Amelia Serrano     yellow 145 124
1219 Alan Nance         red 210 192
1246 Ravi Sinha         yellow 194 177
1078 Ashley McKnight    red 127 118
1221 Jim Brown          yellow 220 .
1095 Susan Stewart      blue 135 127
1157 Rosa Gomez         green 155 141
1331 Jason Schock       blue 187 172
1067 Kanoko Nagasaka    green 135 122
1251 Richard Rose       blue 181 166
1333 Li-Hwa Lee         green 141 129
1192 Charlene Armstrong yellow 152 139
1352 Bette Long         green 156 137
1262 Yao Chen           blue 196 180
1087 Kim Sikorski       red 148 135
1124 Adrienne Fink      green 156 142
1197 Lynne Overby       red 138 125
1133 John VanMeter      blue 180 167
1036 Becky Redding      green 135 123
1057 Margie Vanhoy      yellow 146 132
1328 Hisashi Ito        red 155 142
1243 Deanna Hicks       blue 134 122
1177 Holly Choate       red 141 130
1259 Raoul Sanchez      green 189 172
1017 Jennifer Brooks    blue 138 127
1099 Asha Garg          yellow 148 132
1329 Larry Goss         yellow 188 174
;
proc print data=club;
title'Rakesh Ranjan';
run;

/*Random sample consisting of roughly 10% of observation in a data set*/
data sample;
set club;
if ranuni(12345)<0.1;
run;

proc print data=sample;
run;
/* Random sampling to extract exactly n observation from a data set*/
data some;
retain k 15 n ;
drop k n;
set club nobs=nn;
if _n_ = 1 then n = nn;
if ranuni(0) < k / n then do;
output;
k = k - 1;
end;
if k = 0 then stop;
n = n - 1;
run;







