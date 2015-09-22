WS,6;
per=read_ascii("vper.dat");
den=read_ascii("vden.dat");
w=where((per>0)&(den>0));
per=per(w);den=den(w);
pair = [per,den];
x=spanl(min([per]),max([per]),100);
y=spanl(min([den]),max([den]),100);
h=histo2d(pair,x,y);
//logxy,1,1;
pli,smooth(h),x(min),y(min),x(max),y(max);
logxy,1,1;
xytitles,"persistence","density";

