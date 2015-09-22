#include "/home/thierry/prog/morse/yorick/ASCIIskel.i"
#include "/home/thierry/prog/morse/yorick/ASCIIpairs.i"



fname="ic5146_col_dens.fits";
mask_fname="ic5146_cur_skel_mask.fits~";
// the '~' means the '0' correspond to masked regions
// by default, the reverse convention is used by mse ...



cmd="mse "+fname;
if (strlen(mask_fname)) cmd+=" -mask "+mask_fname;
cmd+= " -ppairs_ASCII -no_sklDump -no_downSkl -no_manifolds -withBoundary";
write,"executing :"+cmd;
system,cmd;

pairs=ASCIIpairs_read(fname+".ASCIIpairs");
drawPDiagram(pairs,win=0,dpi=150);

level=float(0);
write,"Give me a persistence threashold :";
read(level);

cmd="mse "+fname;
if (strlen(mask_fname)) cmd+=" -mask "+mask_fname;
cmd+=" -loadMSC "+fname+".MSC";
cmd+=" -cut "+swrite(format="%e",level);
cmd+=" -no_downSkl -no_manifolds -withBoundary";
write,"executing :"+cmd;
system,cmd;

cmd="skelconv "+fname+".up.NDskl -smooth 10 ";
write,"executing :"+cmd;
system,cmd;

skl=ASCIIskel_read(fname+".up.NDskl.ASCIIskl");
d=fits_read(fname);
w=where(d!=d);
if (numberof(w)) d(w)=0;

winkill,1;
window,1,dpi=150;
pli,d;
plASCIIskel,skl,filcol=__blue;
winkill,2;
window,2,dpi=150;
g=skl2fits(skl,fname+"_skl.fits",fname);
pli,int(g);


