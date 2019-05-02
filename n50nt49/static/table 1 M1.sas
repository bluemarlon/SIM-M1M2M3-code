proc import datafile="C:\Users\rec14003\Desktop\staticont2iid.csv" out=mydata1 dbms=csv replace;
    getnames=yes;
run;

proc mixed data=mydata1;
	class PID Day;
	model yjt=X1 X2/
				 s;
	random X2/subject=pid;
run;

proc glimmix data=mydata1;
	class PID Day;
	model yjt=X1 X2/
				 s dist=normal;
	random X2/subject=pid;
	output out=igausout1 pred(blup ilink)=prediction resid(blup ilink)=r;
run;


   proc nlmixed data=mydata1;
      parms b0=2 b2=1 b1=1 s2u=1 s2e=1;
      den = b0 + b1*X1 + (b2+u1)*X2;
      model yjt ~ normal(den,s2e);
      random u1 ~ normal(0,s2u) subject=PID;
   run;
