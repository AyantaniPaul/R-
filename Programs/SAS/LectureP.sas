/*The KEEP statement causes a DATA step to write only the variables that you specify to one or more SAS data sets*/
data rr;
set stpsamp.stpsale;
keep region citysize pop;
run;

proc print data=rr;
run; 

data rr (keep=region citysize pop);
set stpsamp.stpsale;
run;

proc print data=rr;
run; 

data rr1;
set stpsamp.stpsale;
drop region citysize pop; *Excludes variables from output SAS data sets
run;

proc print data=rr1;
run; 

data rr4 (keep=region citysize pop)
     rr3 (drop=region citysize pop);
set stpsamp.stpsale;
run;

data rr5;
set stpsamp.stpbgt;
run;

data rr6;
set rr5;
NEWBGT=BUDGET*1.10;
run;


data rr7;
set rr5;
if budget gt 20000 then comments='High';
else if budget le 20000 then comments ='Low';
run;

proc print data=rr7;
run;

data rr8;
set rr5;
budget=budget*1.1;
run;

proc print data=rr8;
run;

data rr9;
set rr5;
length comments $ 30;
if budget gt 20000 then comments='Budget of Department is High';
else if budget le 20000 then comments ='Budget of Department is Low';
run;

data rr10;
set rr5;
if budget gt 20000 then delete;
run;



