struct ASCIIskel_crit_Struct {
  int type;
  pointer pos;
  double val;
  int pair;
  int nfil;
  int boundary;
  
  pointer destID;
  pointer filID;
}

struct ASCIIskel_fil_Struct {
  int srcID;
  int destID;
  pointer pos;
}

struct ASCIIskelStruct {
  int ndims;
  int ncrit;
  int nfil;
  pointer crit;
  pointer fil;
}

func ASCIIskel2Segments(skl)
{
  /*
    DOCUMENT ASCIIskel2Segments(skl)
    converts an ASCIIskl to a simple [ndims,2,N] array of segments.
    save with :
    input="filename.ASCIIskl";
    skl=ASCIIskel_read(input);
    sklSegs=ASCIIskel2Segments(skl);
    smwrite(input+".SEGS",transpose(sklSegs(*,)),head="ra1 dec1 z1 ra2 dec2 z2");
  */
  
  p=[];
  q=[];
  for (i=1;i<=skl.nfil;i++)
    {
      f=*(*skl.fil)(i).pos;
      n=dimsof(f)(0);
      grow,p,f(,:n-1);
      grow,q,f(,2:);
    }
  //return [p,q];
  return transpose([transpose(p),transpose(q)],[3,2,1]);
}

func ASCIIskel_read(filename)
{
  file = open(filename,"r");
  line=rdline(file);
  if (line != "NDSKEL")
    error,"wrong format.";

  res=ASCIIskelStruct();
  line=rdline(file);
  tmp=array(int,1);sread,line,tmp;res.ndims=tmp(1);
  line=rdline(file);
  if (line != "[CRITICAL POINTS]")
    error,"wrong format.";
  ndims=res.ndims;
  
  line=rdline(file);
  tmp=array(int,1);sread,line,tmp;res.ncrit=tmp(1);
  crit=array(ASCIIskel_crit_Struct,res.ncrit);
  
  tmp=array(double,4+ndims);
  tmp1=array(int,1);
  for (i=1;i<=res.ncrit;i++)
    {
      line=rdline(file);
      sread,line,tmp;
      crit(i).type=int(tmp(1));
      crit(i).pos=&tmp(2:2+ndims-1);
      crit(i).val=tmp(ndims+2);
      crit(i).pair=int(tmp(ndims+3));
      crit(i).boundary=int(tmp(ndims+4));

      line=rdline(file);
      sread,line,tmp1;res.nfil=tmp1(1);
      if (res.nfil)
        {
          tmp2=array(int,2*res.nfil+1);
          sread,line,tmp2;
          crit(i).destID = &(tmp2(2::2));
          crit(i).filID = &(tmp2(3::2));
        }
    }
  
  line=rdline(file);
  if (line != "[FILAMENTS]")
    error,"wrong format.";
  
  line=rdline(file);
  tmp=array(int,1);sread,line,tmp;res.nfil=tmp(1);
  fil=array(ASCIIskel_fil_Struct,res.nfil);
  tmp=array(int,3);

  for (i=1;i<=res.nfil;i++)
    {
      line=rdline(file);
      sread,line,tmp;
      fil(i).srcID = tmp(1);
      fil(i).destID = tmp(2);
      np=tmp(3);
      pos=array(float,ndims*np);
      cur=0;nr=0;
      while (np>0) {
        if (np>=256) nr=256; else nr=np;
        tmp1=array(float,ndims*nr);
        np-=nr;
        
        line=rdline(file);
        sread,line,tmp1;
        pos(1+ndims*cur:ndims*(cur+nr)) = tmp1;
        //fil(i).pos=&reform(tmp1(4:),[2,ndims,tmp(3)]);
        cur+=nr;
      }
      //write,format="%d=%d\n",cur,tmp(3);
      fil(i).pos=&reform(pos,[2,ndims,tmp(3)]);
    }

  res.fil=&fil;
  res.crit=&crit;
  
  return res;
}

func plASCIIskel(skl,filcol=)
{
  fil=(*skl.fil);
  crit=(*skl.crit);
  for (i=1;i<=skl.nfil;i++)
    {
      //if (i>2) continue;
      p=*fil(i).pos;
      m=sqrt(pow(p(,dif),2)(sum,))(max); 
      if (m>2) {
      plg,p(2,),p(1,),color=filcol;
      //return;
      }
    }
  //return;
  col=char([[255,0,0],[0,255,0],[0,0,255]]);
  for (i=1;i<=skl.ncrit;i++)
    {
      //plmk2,(*crit(i).pos)(2),(*crit(i).pos)(1),msize=0.25,marker=crit(i).type,incolor=col(crit(i).type,),width=10;
      plmk,(*crit(i).pos)(2),(*crit(i).pos)(1),msize=0.25,marker=crit(i).type,color=col(crit(i).type,);//,width=10;
    }
}

func skl2fits(asciiSkl,output_fname,template_fits, res=, tagID=)
{
  /*
    DOCUMENT skl2fits(asciiSkl,output_fname,template_fits, res=)
    res : the resolution of the output image (array(int,2)) (optionnal)
    template_fit : the source image, used to get the resolution (optionnal)
    if 'template_fits' can be ommited only if 'res' is given.
   */
  
  if (is_void(res))
    {
      if (is_void(template_fits))
        error,"I need resolution or a fits template.";
      else {
        fh=fits_open(template_fits);
        res=fits_get_dims(fh)(2:);
        fits_close,fh;
      }
    }
  res=int(res);
  res2=array(int(2),1);
  grow,res2,res;
  g=array(int(0),res2);
  xres=res(1);
  fil=(*asciiSkl.fil);

  if (is_void(tagID)) dotag=0; else dotag=1;
  
  for (i=1;i<=asciiSkl.nfil;i++)
    {
      if (dotag) tagval=i; else tagval=1;
      p=*fil(i).pos;
      w=where((p(1,)<0)|(p(1,)>=res(1))|(p(2,)<0)|(p(2,)>=res(2)));
      if (numberof(w)) continue;
      w=1+ int(p(1,))+int(p(2,))*xres;
      g(w)=tagval;
      w=1+ int(p(1,zcen))+int(p(2,zcen))*xres;
      g(w)=tagval;
    }
  
  system,"rm "+output_fname;
  fits_write(output_fname,g);
  return g;
}
