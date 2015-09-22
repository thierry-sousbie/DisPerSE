
struct ASCIIpairsStruct {
  int ndims;
  pointer npairs;
  pointer pairs;

}

func ASCIIpairs_read(filename)
{
  file = open(filename,"r");
  line=rdline(file);
  if (line != "PERSISTENCE PAIRS")
    error,"wrong format.";

  res=ASCIIpairsStruct();
  npairs=[];
  pairs=[];
  line=rdline(file);
  tmp=array(int,1);sread,line,tmp;res.ndims=tmp(1);
  for (i=0;i<res.ndims;i++)
    {
      line=rdline(file);
      tmp=array(int,2);sread,line,tmp;grow,npairs,tmp(1);
      v=array(double,[2,2,npairs(0)]);info,v;
      read,file,v;
      grow,pairs,&v;
    }
  res.npairs = &npairs;
  res.pairs = &pairs;
  return res;
}

func drawPDiagram(d,win=,dpi=)
{
  if (is_void(dpi)) dpi=100;
  if (!is_void(win))
    {
      winkill,win;
      window,win,dpi=dpi;
    }
  
  //d=ASCIIpairs_read(filename);
  p=*(*d.pairs)(1+1);
  w=where(p(1,)>0);
  plmk,p(dif,w)(*),(p(1,w)),msize=0.25,marker=4;//,width=10,incolor=__red;
  logxy,1,1;
}
