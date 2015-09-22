#include "~/prog/morse/yorick/persistence.i"

windelta=1;

fname="ic5146_col_dens";
use_mask=1;
mask_fname="ic5146_cur_skel_mask";
RECOMP=0;
filename=fname+".fits";
f=fits_read(fname+".fits");

have_mask=0;
w=where(f!=f);
if (numberof(w)) {
  //w2=where((f==f)&(f!=0));
  //v=max(f(*)(w2))*1.E+4;
  //w2=[];
  f(*)(w)=0;//(random(numberof(w))-0.5)*v;
  mymask=array(char(0),dimsof(f));
  mymask(*)(w)=1;
  have_mask=1;
 } else if (use_mask) mymask=array(char(0),dimsof(f));

if (use_mask) {
  w=where(fits_read(mask_fname+".fits") == 0);
  mymask(*)(w)=1;
  have_mask=1;
 }
//d=fft_smooth(f,1);
d=double(f);

NDfield_write(d,"field.ND");
if (RECOMP) {
  if (have_mask) {
    NDfield_write(mymask,"mask.ND");
    cmd="mse "+filename+" -mask mask.ND -nsig 0 -cut 0 -no_sklDump -no_manifolds -ppairs";
  } else {
    cmd="mse "+filename+" -nsig 0 -cut 0 -no_sklDump -no_manifolds -ppairs";
  }
  cmd;exec,cmd;
 }

WS,3*windelta+0,dpi=150;
ppairs=NDnetwork_read(""+filename+".ppairs.NDnet");
pairs=getNDnetPairsId(ppairs);
crit=getCritDist(pairs);
NDnetPPlot(crit,pldelta=1);
logxy,1,1;

cutat="0.7E21";

if (have_mask) {
  cmd="mse "+filename+" -mask mask.ND -cut "+cutat+" -nsig 0 -no_downSkl -no_manifolds -loadMSC "+filename+".MSC -ppairs -outName "+filename+"_cut -noTags";
 } else {
  cmd="mse "+filename+" -nsig 0 -cut "+cutat+" -no_downSkl -no_manifolds -loadMSC "+filename+".MSC -outName "+filename+"_cut -noTags";
 }
cmd;exec,cmd;

skl=NDskel_read(""+filename+"_cut.up.NDskl");
x0=*skl.x0;delta=*skl.delta;
WS,3*windelta+1,dpi=150;
pli,d,x0(1),x0(2),x0(1)+delta(1),x0(2)+delta(2);
c_c;
plNDskel,skl,width=1,withnodes=1;

threshold="5E20";
angle="45";

cmd="skelconv "+filename+"_cut.up.NDskl -smooth 3 -breakdown -assemble "+threshold+" "+angle;
cmd;exec,cmd;
cmd="skelconv "+filename+"_cut.up.NDskl.S003.BRK.ASMB -smooth 3";
cmd;exec,cmd;
skl=NDskel_read(""+filename+"_cut.up.NDskl.S003.BRK.ASMB.S003");
x0=*skl.x0;delta=*skl.delta;
WS,3*windelta+2,dpi=150;
e=d;
w=where(e<str2double(threshold));
e(w)=str2double(threshold);w=[];
pli,d+0*log(e),x0(1),x0(2),x0(1)+delta(1),x0(2)+delta(2);
c_c;
plNDskel,skl,width=1,withnodes=1;
