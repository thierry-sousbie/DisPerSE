N=40
for (i=0,skl=[];i<N;i++) grow,skl,NDskel_read(fname+swrite(format="%d",i));
WS,0;
for (i=0;i<N;i++) plNDskel,skl(20-i),color=char([30+i*8,20,220-i*10]),width=2;
for (i=0,len=[];i<N;i++) grow,len,sqrt(((double(*skl(N-i).segpos))(,dif,)(,1,)^2)(sum,))(sum);
WS,1;
x=pow(10.,span(-8,-4+0.1,40));
plg,len,x,width=2,color=__red;
