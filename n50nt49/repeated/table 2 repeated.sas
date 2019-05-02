proc import datafile="C:\Users\rec14003\Desktop\ofercontsim.csv" out=mydata1 dbms=csv replace;
    getnames=yes;
run;

proc mixed data=mydata1;
	class PID Day;
	model yjt=X1 X2/
				 s;
	random X2/subject=pid;
	repeated day / subject= pid type=ar(1);
run;

proc glimmix data=mydata1;
	class PID Day;
	model yjt=X1 X2/
				 s dist=normal;
	random X2/subject=pid;
	random day / subject= pid type=ar(1) residual;
	output out=igausout1 pred(blup ilink)=prediction resid(blup ilink)=r;
run;
