clear

[header,data]=hdrload('correlogramsOut.dat');
auto1=data(:,1);
xcor1=data(:,2);
auto2=data(:,3);
xcor2=data(:,4);

plot(xcor2)