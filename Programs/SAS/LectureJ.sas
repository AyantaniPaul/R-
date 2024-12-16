OPTIONS PS=58 LS=80 NODATE NONUMBER;

DATA s210006 s310032 s410010;
   set rr.icdblog;
        if (subj = 210006) then output s210006;
   else if (subj = 310032) then output s310032;
   else if (subj = 410010) then output s410010;
 RUN;

 PROC PRINT data = s210006 NOOBS;
   title 'The s210006 data set';
 RUN;
 
 PROC PRINT data = s310032 NOOBS;
   title 'The s310032 data set';
 RUN;

PROC PRINT                NOOBS;
   title 'The s410010 data set';
 RUN;

LIBNAME stat481 'C:\Myrawdata';
DATA stat481.subj210006 stat481.subj310032;
   set stat481.icdblog;
   if (subj = 210006) then output stat481.subj210006;
 RUN;

 PROC PRINT data = subj210006 NOOBS;
   title 'The subj210006 data set';
 RUN;
 
 PROC PRINT data = subj310032 NOOBS;
   title 'The subj310032 data set';
 RUN;

DATA stat481.subj210006 stat481.subj310032;
   set stat481.icdblog;
   if (subj = 210006) then output stat481.subj210006;
   else output stat481.subj310032;
 RUN;

 PROC PRINT data = stat481.subj210006 NOOBS;
   title 'The subj210006 data set';
 RUN;
 
 PROC PRINT data = stat481.subj310032 NOOBS;
   title 'The subj310032 data set';
 RUN;


  DATA symptoms visitsix;
   set stat481.icdblog;
   if form = 'sympts' then output symptoms;
   if v_type = 6 then output visitsix;
 RUN;

 PROC PRINT data = symptoms NOOBS;
   title 'The symptoms data set';
 RUN;

 PROC PRINT data = visitsix NOOBS;
   title 'The visitsix data set';
 RUN;
 
DATA grades;
    input idno 1-2 l_name $ 5-9 gtype $ 12-13 grade 15-17;
    cards;
10  Smith  E1  78
10  Smith  E2  82
10  Smith  E3  86
10  Smith  E4  69
10  Smith  P1  97
10  Smith  F1 160
11  Simon  E1  88
11  Simon  E2  72
11  Simon  E3  86
11  Simon  E4  99
11  Simon  P1 100
11  Simon  F1 170
12  Jones  E1  98
12  Jones  E2  92
12  Jones  E3  92
12  Jones  E4  99
12  Jones  P1  99
12  Jones  F1 185
;
RUN;

PROC PRINT data = grades NOOBS;
   title 'The grades data set';
RUN;

DATA exams;
  set grades (where = (gtype in ('E1', 'E2', 'E3', 'E4')));
RUN;

/*One of the most powerful uses of a RETAIN statement is to compare values 
across observations. The following program uses the RETAIN statement to 
compare values across observations, and in doing so determines each 
student's lowest grade of the four semester exams*/
DATA lowest (rename = (lowtype = gtype)); 
   set exams;
   by idno;
   retain lowgrade lowtype;
   if first.idno then lowgrade = grade;
   lowgrade = min(lowgrade, grade);
   if grade = lowgrade then lowtype = gtype;
   if last.idno then output;
   drop gtype;
RUN;
/*Rename Specifies new names for variables in output SAS data sets*/
/*The BY statement is used in SAS to instruct the DATA step 
 to process dataset observations in groups, rather than singly*/
/*The RETAIN statement specifies variables whose values are not set 
to missing at the beginning of each iteration of the DATA step*/
PROC PRINT data=exams;
RUN;
PROC PRINT data=lowest;
   title 'Output Dataset: LOWEST';
RUN;







