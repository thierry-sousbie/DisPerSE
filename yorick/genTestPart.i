//#include "~/prog/morse/yorick/persistence.i"
#include "netPart.i"
#include "~/prog/disperse/yorick/network.i"

WS_DELTA=0 * 3;
RECOMP=1;
N=[512,512]*2;
npts=4000;
u=span(-1.1,1.1,N(1))(,-::N(2)-1);
v=span(-1.1,1.1,N(2))(-::N(1)-1,);
uv=[u,v];
cutval=0;
cutval=float(cutval);
X0=-0.2;

func val(x,y)
{
  t=abs(pow(x,1)*pow(y,1));
  a=10+(pi+atan(y,x));
  t*=a;
  w=where((sign(x)>0) & (sign(y)<0));
  //w=where((sign(x)<0) | (sign(y)>0));
  t(w)*=-1;
  t+=(y-x)*2;

  d=sqrt(x*x+y*y);
  t-=10*exp(-d*5);
  
  //t+=random(dimsof(t))*3.;
  stat,t;
  //if (sign(x)>0) t*=sign(y);
  return t;
}

func genCoords(np,frac=)
{
  a=array(double,[2,np,2]);
  if (is_void(frac)) frac=0.1;
  npc = long(np*frac);
  t=span(0.,2*pi,npc);
  t+=2.*pi/npc/2*(random(npc)-0.5);
  a(:npc,1) = cos(t);
  a(:npc,2) = sin(t);
  t=random(np-npc)*2*pi;
  d=sqrt(random(np-npc))
  a(npc+1:np,) = d*[cos(t),sin(t)];
  return a;
}

if (RECOMP) c=genCoords(npts); else c=transpose(read_ascii("pc.dat"));

v=double(val(c(,1),c(,2)));

if (RECOMP) {
  smwrite("pc.dat",c,head="px py");
  exec,"delaunay_2D pc.dat -btype void";
  NDfield_write(v,"pv.ND");
  NDfield_write(val(uv(,,1),uv(,,2)),"pc.ND");
 }

net=NDnetwork_read("pc.dat.NDnet");

w=where((*net.data).name=="field_value")(1);
(*net.data)(w).data=&v;
netP=netPart(net,x0=X0);

NDnetwork_write(*netP(1),"pc.dat.left");
NDnetwork_write(*netP(2),"pc.dat.right");

//net_sep=*netP(1);

if (cutval<=0) {
    write,"mse pc.dat.NDnet  -field pv.ND  -upSkl -downSkl -ppairs -compactify natural -debug";    
    exec,"mse pc.dat.left.NDnet   -upSkl -downSkl -ppairs -compactify natural  -debug";
    exec,"mv dg.NDnet dg.left.NDnet";
    exec,"mv tempNET.NDnet pc.dat.left.B.NDnet";
    exec,"mse pc.dat.right.NDnet  -upSkl -downSkl -ppairs -compactify natural  -debug";
    exec,"mv dg.NDnet dg.right.NDnet";
    exec,"mv tempNET.NDnet pc.dat.right.B.NDnet";
    exec,"mse pc.dat.NDnet  -field pv.ND  -upSkl -downSkl -ppairs -compactify natural  -debug";
  } else {
  write,"mse pc.dat.NDnet  -field pv.ND -cut "+swrite(format="%.3g",cutval)+" -upSkl -downSkl -ppairs -compactify natural -debug";  
  exec,"mse pc.dat.left.NDnet  -cut "+swrite(format="%.3g",cutval)+" -upSkl -downSkl -ppairs -compactify natural  -debug";
  exec,"mv dg.NDnet dg.left.NDnet";
  exec,"mv tempNET.NDnet pc.dat.left.B.NDnet";
  exec,"mse pc.dat.right.NDnet   -cut "+swrite(format="%.3g",cutval)+" -upSkl -downSkl -ppairs -compactify natural  -debug";
  exec,"mv dg.NDnet dg.right.NDnet";
  exec,"mv tempNET.NDnet pc.dat.right.B.NDnet";
  exec,"mse pc.dat.NDnet  -field pv.ND -cut "+swrite(format="%.3g",cutval)+" -upSkl -downSkl -ppairs -compactify natural  -debug";
 }

net_left=NDnetwork_read("pc.dat.left.B.NDnet");
net_right=NDnetwork_read("pc.dat.right.B.NDnet");


if (cutval<=0) strctval=""; else strctval = swrite(format="_c%.3g",cutval);

sklu=NDskel_read("pc.dat.NDnet"+strctval+".up.NDskl");
skld=NDskel_read("pc.dat.NDnet"+strctval+".down.NDskl");
sklu_left=NDskel_read("pc.dat.left.NDnet"+strctval+".up.NDskl");
skld_left=NDskel_read("pc.dat.left.NDnet"+strctval+".down.NDskl");
sklu_right=NDskel_read("pc.dat.right.NDnet"+strctval+".up.NDskl");
skld_right=NDskel_read("pc.dat.right.NDnet"+strctval+".down.NDskl");
pairs=NDnetwork_read("pc.dat.NDnet"+strctval+".ppairs.NDnet");
pairs_left=NDnetwork_read("pc.dat.left.NDnet"+strctval+".ppairs.NDnet");
pairs_right=NDnetwork_read("pc.dat.right.NDnet"+strctval+".ppairs.NDnet");
dg=NDnetwork_read("dg.NDnet");
dg_left=NDnetwork_read("dg.left.NDnet");
dg_right=NDnetwork_read("dg.right.NDnet");
/*
sklu=NDskel_read("pc.ND.up.NDskl");
skld=NDskel_read("pc.ND.down.NDskl");
*/

WS,WS_DELTA+1,dpi=175;
pli,val(uv(,,1),uv(,,2)),-1.1,-1.1,1.1,1.1;rmtick;
c_c;
//c_bb;
//plNDnet,net
plNDskel,sklu,color=__red,withnodes=1,width=2;
plNDskel,skld,color=__blue,width=2;
plNDnet,pairs,color=__orange,width=2,nonper=1;
plNDnet,pairs_right,color=__green,width=2;
plNDnet,pairs_left,color=__pink,width=2,nonper=1;

WS,WS_DELTA+2,dpi=175;
pli,val(uv(,,1),uv(,,2)),-1.1,-1.1,1.1,1.1;rmtick;
c_bb;
plNDnet,net;
plNDnet,dg,color=__red,width=3;

func NDNetDataIndex(net,type,txt)
{
    w=where((*net.data).type == type);
    if (numberof(w)==0) return -1;
    
    w=w(where(((*net.data).name)(w) == txt));
    if (numberof(w)==1) return w(*)(1);

    return -1;
}
//eps,"2d_withboundary.ps";
WS,WS_DELTA+3,dpi=150;

pairT_id=NDNetDataIndex(pairs,1,"type");
if (pairT_id<0) error,"";
type=*(*pairs.data)(pairT_id).data;
den_id=NDNetDataIndex(pairs,0,"field_value");
if (den_id<0) error,"";
den=*(*pairs.data)(den_id).data;

vid=*pairs.f_vertexIndex(1);
per=abs(den(vid)(dif,));
Xval=den(vid)(min,);
color=[__blue,__green,__red];

for (i=0;i<=1;i++)
  {
    w=where(type==i+1);
    plmk2,per(w),Xval(w),incolor=color(,i+1),width=11,msize=0.4,marker=i+1;
  }

WS,WS_DELTA+4,dpi=175;
plNDnetData1,net_right,"field_value";
plNDnet,net_right;
plNDskel,sklu_right,color=__red,withnodes=1,width=2;
plNDskel,skld_right,color=__blue,width=2;
plNDnet,pairs_right,color=__green,width=2
plNDskel,sklu,color=__black;
plNDskel,skld,color=__orange;

WS,WS_DELTA+5,dpi=175;
plNDnetData1,net_left,"field_value";
plNDnet,net_left;
plNDskel,sklu_left,color=__red,withnodes=1,width=2;
plNDskel,skld_left,color=__blue,width=2;
plNDnet,pairs_left,color=__green,width=2
plNDskel,sklu,color=__black;
plNDskel,skld,color=__orange;

  
WS,WS_DELTA+6,dpi=175;
pli,val(uv(,,1),uv(,,2)),-1.1,-1.1,1.1,1.1;rmtick;
c_bb;
plNDnet,net_right;
plNDnet,net;
//plNDnet,net_sep,color=__white,width=5;
plNDnet,dg,color=__red,width=5;
plNDnet,dg_right,color=__green,width=3;
plNDnet,dg_left,color=__purple,width=2;

WS,WS_DELTA+7,dpi=175;
//plNDnet,net_sep,width=5,color=__grey;
plNDnet,net_right;
plNDnet,net;
plNDnet,dg,color=__red,width=5;
plNDnet,dg_right,color=__green,width=2,nonper=1;
plNDnet,pairs_right,color=__yellow,width=2
plNDnet,dg_left,color=__orange,width=2;
plNDnet,pairs_left,color=__purple,width=2,nonper=1;
