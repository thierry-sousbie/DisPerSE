#! /usr/lib/yorick/bin/yorick -batch
#include "$HOME/Yorick/Chris/utls.i"

args=get_argv()(2:); // waits for arg

if (numberof(args) ==0) error,"Expecting arguments: fname [dmin dmax]";
fname=args(1);
if (numberof(args) == 2) error,"Expecting arguments: fname [dmin dmax]";
dmax=-1;
if (numberof(args) == 3) {
  dmin=str2double(args(2));
  dmax=str2double(args(3));
 }
  
file = xopen(fname, compress="auto");
header = rdline(file);
header = streplace(header,strfind("#",header),"");
header=strtrim(header);
header = strtok(header," ",100);
header=header(where(strlen(header)));
file=[];

ra=where(header == "ra");
dec=where(header == "dec");
dist=where(header == "dist");
if (is_void(dist)) dist=where(header == "z");
data=read_ascii(fname);

write,"statistics on distance :";
stat,data(dist,);

w=where((data(dec,)<280)&(data(dec,)>100)& ((data(dec,)<247)|(data(ra,)<52)));
data=data(,w);
if (dmax>0) {
  w=where((data(dist,) > dmin) & (data(dist,) < dmax));
  data=data(,w);
 }

h=header(1);
for (i=2;i<=numberof(header);i++) h+=" "+header(i);
write,format="Saving %d galaxies ...",numberof(data(1,));
smwrite("main_"+fname,transpose(data),head=h);
