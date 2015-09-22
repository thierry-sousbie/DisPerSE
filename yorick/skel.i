#include "Chris/randfield.i"
#include "Thierry/NDfield.i"
#include "Thierry/NDskel.i"


func skel2gad(fname,redshift=)
/* DOCUMENT 
     converts a skeleton into a gadget file
   SEE ALSO:
 */
{
  if(is_void(redshift)) redshift=0.;
  if(is_string(fname)) skl=skel_read(fname,arr=1); else skl=fname;
  q=GadgetSnapshot();
  q.pos=&skl;
  q.redshift=redshift;
  q.npart=[0,dimsof(skl)(0),0,0,0,0];
  q.nparttotal=[0,dimsof(skl)(0),0,0,0,0];
  q.boxsize=max(skl)-min(skl);
  return q;
}


struct SkelStruct {
  /* */
  double pos(3); // position 
  double dens; // local density
  double evec(3,3); //eigenvectors
  double eval(3); // eigenvalues
  double rad(2); // local curvature radii
}


func skel_read(fname,arr=,width=)
/* DOCUMENT 
     reads a skeleton file
     contains
     the list dimsof each segments
   SEE ALSO:
 */
{
  nd=1;
  if(width) arr=1;
  file=open(fname,"r");
  read,file,nd;
  nd;
  if(nd<0) { nd=-nd;
    d =array(0.,6*nd);
    read,file,d;
    return (reform(d,[2,3,2*nd]));
  }
  lst=[];
  d =array(0,nd);
  read,file,d;
  if(arr)
  {
    lst=array(0.,[2,3,d(sum)]);
  read,file,lst;
  }
  else
      for(i=1;i<=numberof(d);i++)
        {
          a=array(0.,[2,3,d(i)]);
          read,file,a;
          lst= _cat(lst,a);
        }
if(!width)      return lst;

 q=array(SkelStruct,d(sum));
 q.pos=lst;
 lst=array(0.,[1,d(sum)]); //read density
 read,file,lst;
 q.dens=lst;
 lst=array(0.,[3,3,3,d(sum)]); //read eigenvectors
 read,file,lst; 
 q.evec=lst;
 lst=array(0.,[2,3,d(sum)]); //read eigenvalues
 read,file,lst; 
 q.eval=lst;
 lst=array(0.,[2,2,d(sum)]); //read radii
 read,file,lst; 
 q.rad=lst;
  return q;
}


func skel_write(skl,fname)
/* DOCUMENT 
     writes a skeleton file
   SEE ALSO:
 */
{
  error,"not finished yet !";
  file=open(fname,"w");
  if(is_list(skl))
    {
      write,file,format="%d \n",_len(skl);
      nn=lst2arr(_map(numberof,skl))/3;
      write,file,format="%d ",nn;
      for(i=1;i<=_len(skl);i++)
        {
          x=_nxt(skl);
          write,file,format="%1.6g %1.6g %1.6g \n",x,linesize=25;
        }
    }
      else error,"skl needs to be a list";
  close,file;
  return fname;
}





func skelsad_read(fname,&sdens,dens=)
/* DOCUMENT 
     reads the saddle points of skeleton file
     contains the list of maxima, saddle panckake, saddle filaments minima
     SEE ALSO:
 */
{
  nmax=nsadp=nsads=nmin=0;
  file=open(fname,"r");
  read,file,nmax,nsadp,nsads,nmin;
  a=array(0.,[2,3,nmax+nsadp+nsads+nmin]);
      read,file,a;
      lst= _cat(_lst(a(,1:nmax)),
                _lst(a(,nmax+1:nmax+nsadp)),
                _lst(a(,nmax+1+nsadp:nmax+nsadp+nsads)),
                _lst(a(,nmax+nsadp+nsads+1:0)));
if(dens)   a=array(0.,[2,13,nmax+nsadp+nsads+nmin]);
      read,file,a;
      sdens=_cat(_lst(a(1,1:nmax)),
                _lst(a(1,nmax+1:nmax+nsadp)),
                _lst(a(1,nmax+1+nsadp:nmax+nsadp+nsads)),
                _lst(a(1,nmax+nsadp+nsads+1:nmax+nsadp+nsads+nmin)));
  return lst;
}


func rsex_make(fname,verb=,peak=,nosmooth=,smooth=,tol=,voidpeak=,keep=,period=,debug=,withmax=,mask=)
/* DOCUMENT 
   Global Skekelon
u=array(0.,[3,128,128,128]); for(i1=1;i1<=4;i1++)for(i2=1;i2<=4;i2++)for(i3=1;i3<=4;i3++) u(32*(i1)-16,32*(i2)-16,32*(i3)-16)=1;
u1=fft_smooth(u,25);
skl=rsex_make(u1,verb=1);
slice3d,u1;
plNDskel,skl,offset=1,noerase=1;
   SEE ALSO:
 */
{

  if(!is_string(fname))
    {
      fname1=mktemp("temp.ND.",dir=".");
      NDfield_write,fname,fname1;
      //     if(numberof(mask)) NDfield_write,mask,fname1+".mask.ND";
    }
  else
    {
      fname1=fname;
      keep=1;
    }                                                                                 
  if(peak) post=" -peak "; else post="";
  if(voidpeak) post+=" -voidandpeak ";
  if(withmax) post+=" -withextrema ";
  if((mask)) post+=" -mask ";
  if(is_string(period)) post+=" -periodicity  "+period +" ";
  if(period==0) post+=" -periodicity 000";
  if(tol) post+=" -tol "+pr1(tol)+" ";
  if(nosmooth) post+=""; else if(is_void(smooth)) post+=" -smooth"; else post+=" -smooth"+pr1(smooth);
  if(verb) _system=system;  else  if(debug) _system=write; else _system=exec;
  
  er=_system("~/soft/skeleton/bin/rsex.exe -save_bin "+fname1+post);
  
  if( debug) return;
    skl=NDskel_read(fname1+".real.NDskl");
    if(!keep) er=exec("rm "+fname1+".real.NDskl");
  if(!keep) er=exec("rm "+fname1);
  return skl;
}


func GPP_seg(u,reso=,peakpatch=,tol=,skl=)
/* DOCUMENT 
     does segmentation of field u in the GPP sense;
   SEE ALSO:
 */
{
  require,"msort.i";
  if(is_void(reso)) reso=4;
  us=fft_smooth(u,reso);
  dd=dimsof(u);
  if(dd(1)==2) ud=abs(us(dif,)(pcen,))+abs(us(,dif)(,pcen));
  if(dd(1)==3) ud=abs(us(dif,,)(pcen,,))+abs(us(,dif,)(,pcen,))+
                 abs(us(,,dif)(,,pcen));
  if(skl) return rsex_make(NDfield_write(fft_smooth(ud,reso),"tmp.ND"));
  pp=GPP_make(NDfield_write(fft_smooth(ud,reso),"tmp.ND"),tol=tol);
  if(peakpatch)  return pp;
  pp2=pp*0.+max(pp);
  s=heapsort(pp);p1=pp(s);w=where(p1(dif)!=0);w=grow(0,w);w=grow(w,numberof(pp));
  for(i=1;i<numberof(w);i++) {w1=s(w(i)+1:w(i+1)); pp2(w1)=avg(u(w1));}
  return pp2;
}


func GPP_make(fname,&mm,&xyzc,&ee,withpos=,withellip=,verb=,proba=,wght=,peak=,tol=)
/* DOCUMENT 
   Global Peak Patch
   u=GPP_make(NDfield_write(fft_smooth(u0,2),"crap.ND"),mm);

for(i=1;i<4;i+=0.25){ u=fft_smooth(fft(genrandfield2D(1024,fun=PS,param=i),-1).re,5); up=GPP_make(NDfield_write(u,"crap.ND")); window,long(i*4); c_gg; pli,up; }


or by segment;
require,"Chris/split.i";
un=split(u,2);
ppn=[];for(i=1;i<=dimsof(un)(0);i++){ grow,ppn,[GPP_make(un(..,i))];}
pp=splitb(ppn);


   if(wght) is uses the density field to wght the center of mass;
   SEE ALSO:
 */
{
if(!is_string(fname))
    {
      fname1=mktemp("temp.ND.",dir=".");
      NDfield_write,fname,fname1;
    }
  else
    {
      fname1=fname;
       keep=1;
    }
  
  require,"msort.i";
  if(withellip) withpos=1;
  if(verb) _system=system;  else _system=exec;
  if(proba) post=" -save_proba "; else post="";
  if(peak) post+=" -peak ";
  if(tol) post+=" -tol "+pr1(tol)+" ";
  er=_system("~/soft/skeleton/bin/rsex.exe  -save_patch "+post+fname1);
  if(proba) return NDfield_read(fname+".proba");
if(peak)  u=NDfield_read(fname1+".ppatch"); else   u=NDfield_read(fname1+".vpatch");
  s=heapsort(u);
  u1=u(s); w=where(u1(dif)!=0);
  mm=w(dif); w=grow(1,w); ///w=grow(w,numberof(u));//FIXME: !!!!!!!!!
  if(withpos)
    {// FIXME: deal with periodicity
      d=dimsof(u);
      x=indgen(d(2))(,-,-);
      y=indgen(d(3))(-,,-);
      z=indgen(d(4))(-,-,);
      xyz=[x,y,z]; xyz=xyz(*,)(s,);
      if(wght) us=u(*,)(s,);
      xyzc=array(0.,[2,3,numberof(mm)]);
      for(i=1;i<=numberof(mm);i++){
      if(wght) xyzc(,i)=(xyz*us)(w(i):w(i+1),)(avg,)/(us(w(i):w(i+1))(avg));
        else  if(peak)  
          xyzc(,i)=xyz(w(i):w(i+1),)(us(w(i):w(i+1))(mxx),);
       else
          xyzc(,i)=xyz(w(i):w(i+1),)(,); 
      }
  }
  if(withellip)
    {
      ee=array(0.,[2,3,numberof(mm)]);
      for(i=1;i<=numberof(mm);i++){
        tt=xyz(w(i):w(i+1),)-xyzc(,i)(-,);
        x=tt(,1);y=tt(,2);z=tt(,3);
        x=abs(min(abs(x),d(2)-abs(x))); // FIXME: check periodicity
        y=abs(min(abs(y),d(3)-abs(y)));
        z=abs(min(abs(z),d(4)-abs(z)));
        II=[[x*x,x*y,x*z],[x*y,y*y,y*z],[x*z,z*y,z*z]];
        if(wght) II *= us(w(i):w(i+1))/avg(us(w(i):w(i+1)));
        inert=II(avg,,);
        ee(,i)=eigs3D(inert);
      }

    }
  return u;
}



func genus(sd,rcrit)
/* DOCUMENT 
EXAMPLE
skl=skel_make("simu.gad",sad,sd);
or
sad=skelsad_read("simu.gad.ext.ascii",sd,dens=1);
tt=genus(sd);
plg,tt(2,),tt(1,); logxy,1,0;
SEE ALSO:
 */
{
  M=_car(sd,1);
  s1=_car(sd,2);
  s2=_car(sd,3);
  m=_car(sd,4);
if(is_void(rcrit))  rcrit=spanl(1e-3,5,100);
  tt=array(0.,numberof(rcrit));
  for(i=1;i<=100;i+=1)
    {
      dcrit=rcrit(i);
      tt(i)=(M>dcrit)(sum)-(s1>dcrit)(sum)+(s2>dcrit)(sum)-(m>dcrit)(sum);
    }
  return transpose([rcrit,tt]);
}


func NDskel_curv(skl)
/* DOCUMENT 
     returns skel curvature
   SEE ALSO:
 */
{
    dx=(*skl.segpos)(,dif,)(,1,);
    d2x=dx(,dif)(,pcen);
    dx=transpose(dx);
    d2x=transpose(d2x);
    ndx=(dx*dx)(,sum);
    ncrs=cross(dx,d2x);if(skl.ndims>2) ncrs=(ncrs*ncrs)(,sum); else ncrs=(ncrs*ncrs);
    w=where(ncrs!=0);
    curv=ndx(w)^2/abs(ncrs(w))^0.5; // includes ds
    cv=curv(sum)/((abs(ndx(w))^0.5)(sum));
  return cv;
}



func NDskel_tors(skl)
/* DOCUMENT
     returns skel torsion
   SEE ALSO:
 */
{

  if(skl.ndims==2) return 0; // That's why it's called "ND"...

  dx=(*skl.segpos)(,dif,)(,1,);
  d2x=dx(,dif)(,pcen);
  d3x=d2x(,dif)(,pcen);

  dx=transpose(dx);
  d2x=transpose(d2x);
  d3x=transpose(d3x);

  ndx=(dx*dx)(,sum);

  ncrs=cross(dx,d2x);

  mixt=array(0.,skl.nsegs);
  for(i=1;i<=skl.nsegs;i++) {
    mixt(i)=ncrs(i,+)*d3x(i,+);
  }

  ncrs=(ncrs*ncrs)(,sum);

  w=where(ncrs!=0);
  tors=mixt(w)*ndx(w)^0.5/ncrs(w); // includes ds
  tr=tors(sum)/((ndx(w)^0.5)(sum));

  return tr;
}





func skel_curv(skl)
/* DOCUMENT 
     returns skel curvature
   SEE ALSO:
 */
{
  cv=array(0.,_len(skl));
  for(i=1;i<=_len(skl);i++) {
    x=transpose(_nxt(skl));
    dx=x(,dif)(,pcen); 
    d2x=dx(,dif)(,pcen);
    //    d3x=d3x(,dif)(,pcen);
    //    mat=[dx,d2x,d3x];
    //    d=map(zLUdet,mat);
    ndx=(dx*dx)(,sum);
    ncrs=cross(dx,d2x); ncrs=(ncrs*ncrs)(,sum);
    w=where(ncrs!=0);
    curv=ndx(w)^2/ncrs(w)^0.5; // includes ds
    cv(i)=curv(sum)/((ndx(w)^0.5)(sum));
  }
  return cv;
}








func skel_tors(skl)
/* DOCUMENT 
     returns skel torsion
   SEE ALSO:
 */
{
  tr=array(0.,_len(skl));
  for(i=1;i<=_len(skl);i++) {
    x=transpose(_nxt(skl));
    dx=x(,dif)(,pcen); 
    d2x=dx(,dif)(,pcen);
    d3x=d2x(,dif)(,pcen);
    mat=[dx,d2x,d3x];
    d=dimsof(x)(2); dd=array(0.,d);
    for(j=1;j<=d;j++) dd(j)=zLUdet(mat(j,,));
    ndx=(dx*dx)(,sum);
    ncrs=cross(dx,d2x); ncrs=(ncrs*ncrs)(,sum);
    w=where(dd!=0);
    tors=ncrs(w)*ndx(w)^0.5/dd(w); // includes ds
    tr(i)=tors(sum)/((ndx(w)^0.5)(sum));
  }
  return tr;
}



func skel_smooth(fname,fwhm=,arr=,boxsize=)
/* DOCUMENT 
     smooths a skeleton 
   SEE ALSO:
 */
{
  if(is_void(boxsize)) boxsize=1;
  if(is_void(fwhm)) fwhm=5;
  if(is_string(fname)) skl=skel_read(fname); else skl=fname;
  skl3=[];
  n=_len(skl);
  for(i=1;i<=n;i++) {
    x=_nxt(skl); d=dimsof(x)(0);
     w=(where(abs(x(,dif))(max,)>boxsize/2));
     np=numberof(w); 
     if(np){ w=grow([0],w);w=grow(w,d);} else {np=1; w=[0,d];}
     for(j=1;j<=np;j++)
      {
        x0=x(,w(j)+1:w(j+1));
        if(dimsof(x0)(0)>2)
          {
             x1=transpose(fsmooth(transpose(x0),min(fwhm,d)));
             x1(,1)=x0(,1);
             x1(,0)=x0(,0);
           } else x1=x0;
        if(arr) grow,skl3,x1; else skl3=_cat(skl3,x1);
      }
  }
  return skl3;
}


func skel_make(fname,&sad,&len,&sdens,smooth=,npix=,keep=,verb=,densfile=,execpath=,length=,debug=,savebin=,arr=,width=,post=,true=)
/* DOCUMENT 
     makes the skeleton of the corresponding
     gadget file
     sad contains the list of saddle points
     len contains the length / numberof max
     if keep is set it won't erase the skeleton files etc..
     if verb is not set it won't display the progress of the skeleton
     EXAMPLE
     cd,"/data1/pichon/SIMUS/";
     ll=exec("ls G*\/seed*_256_*_012");
     skl=skel_make(ll(1));     
   SEE ALSO:
 */
{
  if(is_void(smooth)) smooth=6;
  if(is_void(arr)) arr=1;
  if(is_void(execpath)) execpath="/home/pichon/soft/skeleton/bin/";
  if(is_void(npix)) npix=256;
  if(length)  len=" -length "; else len="";
  if(width)  wdth=" -width "; else wdth="";
  if(savebin)  savbin=" -save_bin "; else savbin="";
  if(is_void(densfile)) densfile=0;
  if(densfile) flag=" -d "; else flag="";
  if(debug) {_exec=write; _system=write; } else {_exec=exec; _system=system; }
  if(smooth)
          {
            er= _exec(execpath+"/smooth.exe  "+flag+fname+" -n "+pr1(npix)+" -s "+pr1(smooth)+" ");
            if(densfile)            sfname=splittok(fname,tok="/")(0)+".s"+pr1(smooth)+".00"; else
              sfname=splittok(fname,tok="/")(0)+".CIC.s"+pr1(smooth)+".00"; 
            if(is_void(post)) post=3;  
            if(true)
              {
                er=_system(execpath+"/rsex.exe  "+sfname+"  -save_ascii > /dev/null");
                skelname=sfname+".real.skl.ascii";
              }
            else
              {
                er=_system(execpath+"/skelex.exe -cut "+sfname+"  -nogl -post "+pr1(post)+" -ascii "+wdth+len+savbin +"> /dev/null");
                skelname=sfname+".skl.ascii";
              }
          }
        else
          {
            sfname=fname;
            keep=2;
            if(true)
              {
                er=_system(execpath+"/rsex.exe  "+sfname+"  -save_ascii > /dev/null");
                skelname=sfname+".real.skl.ascii";
              }
            else
              {
                er=_system(execpath+"/skelex.exe -cut "+sfname+"  -nogl -post "+pr1(post)+" -ascii "+wdth+len+savbin +"> /dev/null");
                skelname=sfname+".skl.ascii";
              }
            if(debug) return skelname;
          }
  skl=skel_read(skelname,arr=arr,width=width);
if(real)  sad=skelsad_read(sfname+".real.ext.ascii",sdens,dens=1);
else  sad=skelsad_read(sfname+".ext.ascii",sdens,dens=1); 
  if(length)
    {
      length=load(sfname+".skl.length");
      len=length/numberof(_car(sad,1));
      write,"skel length/halo",len;
      len=length;
    }
  if(is_void(keep)){
    if(!densfile)    er=system("rm "+sfname);
    er=system("rm "+sfname+".skl.ascii");
    er=system("rm "+sfname+".ext.ascii");
  }
  else
    {
      if(keep==1) if(!densfile)   er=system("rm "+sfname);
    }
  return skl;
}





func skel_dist(sk1,sk2,med=)
  /* DOCUMENT 
     returns the distance between two skeletons
     defined as the mean minimal euclidian distance between
     two segments. 
     EXAMPLE
sk2=skel_read("seed005555444_256_050_018.CIC.s10.00.skl.ascii");
sk1=skel_read("seed545457982_256_050_018.CIC.s10.00.skl.ascii");
dd=skel_dist(sk1,sk2,med=1);
SEE ALSO:
  */
{
  tr1=BuildOctreeFromPos(sk1);
  tr2=BuildOctreeFromPos(sk2);
  n1=dimsof(sk1)(0);
      n2=dimsof(sk2)(0);
      aa1=array(0.,[1,n1]);
      aa2=array(0.,[1,n2]);
      for(j=1;j<=n1;j++)
          {
            s1=sk1(,j); // for speed
             d=0.;id=FindPosClosest(tr2,s1,dist=&d);
            aa1(j)=d;
          }
      for(j=1;j<=n2;j++)
          {
            s2=sk2(,j); // for speed
             d=0.;id=FindPosClosest(tr1,s2,dist=&d);
            aa2(j)=d;
          }
      er=FreeOctree(tr1);
      er=FreeOctree(tr2);
      return [&aa1,&aa2];
}


func skel_pos_dist(sk1,pos) 
  /*     returns the distance between a skeleton and pos
     two segments. 
     EXAMPLE
sk2=skel_read("seed005555444_256_050_018.CIC.s10.00.skl.ascii");
dd=skel_dist(sk1,pos);
SEE ALSO:
  */
{
  tr1=BuildOctreeFromPos(sk1);
  if(dimsof(pos)(1)==2)
    {
  n1=dimsof(pos)(0);
      aa1=array(0.,[1,n1]);
      for(j=1;j<=n1;j++)
          {
            s1=pos(,j); 
             d=0.;id=FindPosClosest(tr1,s1,dist=&d);
            aa1(j)=d;
          }
    } 
  else    {
    tr1=BuildOctreeFromPos(sk1);
    n0=numberof(pos);
    dd=dimsof(pos);
    x=span(0,1,dd(2));
    y=span(0,1,dd(3));
    z=span(0,1,dd(4));
    pos=[x(,-,-),y(-,,-),z(-,,-)](*,); pos=transpose(pos);
    aa1=pos*0;
    for(j=1;j<=n0;j++)
      {
        s1=pos(,j); 
        d=0.;id=FindPosClosest(tr1,s1,dist=&d);
        aa1(j)=d;
      }
   aa1= reform(aa1,dimsof(cube));
  }
  er=FreeOctree(tr1);
  return aa1;
}


func skel_dist_old(sk1,sk2,med=)
  /* DOCUMENT 
     returns the distance between two skeletons
     defined as the mean minimal euclidian distance between
     two segments. 
     EXAMPLE
sk2=skel_read("seed005555444_256_050_018.CIC.s10.00.skl.ascii");
sk1=skel_read("seed545457982_256_050_018.CIC.s10.00.skl.ascii");
dd=skel_dist(sk1,sk2,med=1);
SEE ALSO:
  */
{
  if(is_list(sk1))
    {
      aa=array(0.,[2,_len(sk1),_len(sk2)]);
      /*      for(i=1;i<=_len(sk1);i++)
        {s1=_nxt(sk1)(,::5); i;
        sk=sk2;
        for(j=1;j<=_len(sk2);j++)
          { // this stuff should be faster ... doesn't seem to be !
            s2=_nxt(sk)(,::5);
            tt=s1(,,-)-s2(,-,); tt*=tt; tt=sqrt(tt(sum,,));
            aa(i,j)=min(tt(*));
          }
        }
      */
        for(i=1;i<=_len(sk1);i++)
          {
            s1=_car(sk1,i)(,::5); // for speed
            for(j=1;j<=_len(sk2);j++)
              {
                s2=_car(sk2,j)(,::5);
                tt=s1(,,-)-s2(,-,); tt*=tt; tt=sqrt(tt(sum,,));
                aa(i,j)=min(tt(*));
              }
          }
      if(med)   return (median(aa(min,))+median(aa(,min)))/2.;
      return (aa(min,)(avg)+aa(,min)(avg))/2.;
    }
  else
    {
      n1=dimsof(sk1)(0);
      n2=dimsof(sk2)(0);
      aa1=array(0.,[1,n1]);
      aa2=array(0.,[1,n2]);
      for(j=1;j<=n1;j++)
          {
            s1=sk1(,j); // for speed
            s2=sk2(,::5);
            tt=s1(,,-)-s2(,-,); tt*=tt; tt=sqrt(tt(sum,,));
            aa1(j)=min(tt(*));
          }
      for(j=1;j<=n2;j++)
          {
            s1=sk2(,j); // for speed
            s2=sk1(,::5);
            tt=s1(,,-)-s2(,-,); tt*=tt; tt=sqrt(tt(sum,,));
            aa2(j)=min(tt(*));
          }
      
      if(med)   return (median(aa1)+median(aa2))/2.;
      return (aa1(avg)+aa2(avg))/2.;
    }
}


func skelPDF_make(fname,&vv,&mu,&len,smooth=,npix=,keep=,verb=,densfile=,execpath=)
/* DOCUMENT 
     makes the skeleton of the corresponding
     gadget file
     sad contains the list of saddle points
     if keep is set it won't erase the skeleton files etc..
     if verb is not set it won't display the progress of the skeleton
     EXAMPLE
     cd,"/data1/pichon/SIMUS/";
     ll=exec("ls G*\/seed*_256_*_012");
     tt=[];for(i=1;i<=numberof(ll);i++) grow,tt,[skelPDF_make(ll(i),vv,mu)];

//   
   SEE ALSO:
 */
{
  if(is_void(smooth)) smooth=10;
  if(is_void(execpath)) execpath="/home/pichon/soft/skeleton/bin/";
  if(is_void(npix)) npix=256;
  dir="";
  if(is_void(densfile))  skelname=splittok(fname,tok="/")(0)+".CIC.s"+pr1(smooth)+".00.skl"; else skelname=fname+".skl";
  if(numberof(exec("ls "+skelname+dir))!=1)
    {
      if(is_void(densfile))
        {
          er= exec(execpath+"/smooth.exe "+fname+" -n "+pr1(npix)+" -s "+pr1(smooth)+" -center -overd "+dir);
          er=system(execpath+"/skelex.exe -cut "+splittok(fname,tok="/")(0)+".CIC.s"+pr1(smooth)+".00  -nogl -post 3 -save -length"+dir);
        } else
    er=system(execpath+"/skelex.exe -cut "+splittok(fname,tok="/")(0)+"  -nogl -post 3 -length -save"+dir);
      er=system(execpath+"/skelprop.exe -skl "+skelname+" -snap "+fname+" -moments");
    }
      er=system(execpath+"/skelprop.exe -skl "+skelname+" -snap "+fname+" -moments");
      moments=load(skelname+".PDF_V.dat");
  if(is_void(keep))  er=system("rm "+splittok(fname,tok="/")(0)+".CIC.s"+pr1(smooth)+".00");
  if(is_void(keep)) er=system("rm "+splittok(fname,tok="/")(0)+".CIC.s"+pr1(smooth)+".00.skl.ascii");
  if(is_void(keep)) er=system("rm "+splittok(fname,tok="/")(0)+".CIC.s"+pr1(smooth)+".00.ext.ascii");
  vv=reform(moments(,1),[2,100,100]);
  mu=reform(moments(,2),[2,100,100]);
  return reform(moments(,3),[2,100,100]);
}

func NDskelDist(skl,pos,&wght,node=)
/* DOCUMENT 
     returns the distance of pos to skl and the corresponding wghts;
     if node is set it will look at the distance to nodes.
   SEE ALSO:
 */
{
  require,"octree.i";
  if(!node) xyz=(*skl.segpos)(,avg,); else xyz=(*skl.nodepos)(avg,);
  xyz+=RMS(xyz)*random(dimsof(xyz))/1e5;
  tr=BuildOctreeFromPos(xyz(,*));
  d=array(0.,dimsof(pos)(0)); id=FindNeighbours(tr,pos,1,dist=&d);
  FreeOctree(tr);
  wght=(*skl.segdata)(1,id)(1,);
  return d;
}



func NDcountSaddle(skl,&list,debug=)
/* DOCUMENT 
     returns the number of saddle points at the border of PeakPatch


hh=[]; for(i=1;i<=50;i++) {i; u=fft_smooth(fft(genrandfield2D(64,fun=PS,param=0.5),-1),8).re; skl=rsex_make(u,voidpeak=1); grow,hh,[histo1d(NDcountSaddle(skl),span(2,30,29))];}
pler,hh/float(hh(sum,)(-,)),-0.5+span(2,30,29)(zcen),msize=0.5,incolor=-1,mean=1;
xytitles,"N","P(N)";pltitle,"PDF of saddle points";
     SEE ALSO:
 */
{
  p1=where(*skl.segdatainfo=="patch"+" p1")(1);
  p2=where(*skl.segdatainfo=="patch"+" p2")(1);
  seg_dif = (*skl.segdata)(p2,)-(*skl.segdata)(p1,);
  w = where (seg_dif != 0); 
  if(debug) pl2,(*skl.segpos)(,1,w),width=10;
  list = (*skl.segdata)(p2,w); 
  grow,list,(*skl.segdata)(p1,)(w);
  h = histogram(list+1) ;
  return h(sort(h));
}


func NDskelTag(sklfn,fieldfn,execpath=,verb=,debug=,smth=,tagname=,overwrite=)
/* DOCUMENT 

skl=rsex_make(NDfield_write(u1,"test.ND"));
  skl2=NDskelTag("test.ND.real.NDskl","test.ND",verb=1);    
   SEE ALSO:
 */
{
  if(is_void(tagname)) tagname="newtag";
  if(is_void(execpath)) execpath="/home/pichon/soft/skeleton/bin/";
    if(verb) _system=system;  else  if(debug) _system=write; else _system=exec;
    if(smth) er=_system(execpath+"/skelconv.exe "+sklfn+" -smooth"+pr1(smth)+"  -save_ND");
      else er=_system(execpath+"/skelconv.exe "+sklfn+" -tag  "+fieldfn+" "+tagname+ "  -save_ND");
  skl=NDskel_read(sklfn+".NDskl");
  if(overwrite){er=exec("mv "+sklfn+".NDskl "+sklfn+""); }
  return skl;
}
  
  
  func NDgetMinMax(u,pp)
/* DOCUMENT 
     returns the indx of the max of field u on each PP
     u=fft_smooth(fft(genrandfield2D(64),-1),10).re;
     pp=GPP_make(u,peak=1);
     ii=NDgetMinMax(u,pp);
     x=indgen(64)(,-:1:64)/64.; y=transpose(x,[1,2]);
     pli,u,0,0,1,1;
     PL,y(*)(ii(,2)),x(*)(ii(,2)),msize=1,incolor=__white

     or to get max along filaments
      u1=(*skl.segdata)(1,)
      p1=(*skl.segdata)(3,)
      jj=NDgetMinMax(u1,p1)
      PL,(*skl.segpos)(2,2,jj(,2)),(*skl.segpos)(1,2,jj(,2)),msize=1

     SEE ALSO:
 */
{

  require,"msort.i";
  
  idxmax=array(0,1+max(long(pp)));
  idxmin=array(0,1+max(long(pp)));
 s=heapsort(pp);pp1=pp(s);
  un = 1+int(unique(pp1));
  w=where(pp1(dif)!=0);
  idxmax=array(0,1+int(max(pp)));
  w=grow(0,w);w=grow(w,numberof(pp));
  for(i=1;i<numberof(w);i++) {w1=s(w(i)+1:w(i+1));
    idxmax(un(i))=w1(u(w1)(*)(mxx));
    idxmin(un(i))=w1(u(w1)(*)(mnx));
  }
  return [idxmin,idxmax];
}



func NDskel_pad(skl)
/* DOCUMENT
   SEE ALSO:
 */
{
  skl1=skl;
  skl2=skl;
  skl3=skl;
  skl3.comment="Merged skeleton";
  
}


func NDskelonedge(skl)
/* DOCUMENT 
     returns patch index of skeleton on edge of cube;

      tt=NDskelonedge(skl);
      hh=[];for(i=1;i<=numberof(tt);i++) hh=grow(hh,where(tt(i)==pp(*)));
      pp2=pp*0;pp2(hh)=1; pli,pp2;

      SEE ALSO:
 */
{
  p1=where(*skl.segdatainfo=="patch"+" p1")(1);
  p2=where(*skl.segdatainfo=="patch"+" p2")(1);
  pp=(*skl.segdata)(p1,);
  xmin=*skl.x0+*skl.delta/(*skl.dims);
  xmax=*skl.x0+*skl.delta*(1.-1./(*skl.dims));
  w=  where((((*skl.segpos)(,avg,)<=xmin)+((*skl.segpos)(,avg,)>=xmax))(sum,)!=0);
  return unique((*skl.segdata)([p1,p2],w)(*));
}


func NDskel2rskel(skl,debug=,cutedge=)
/* DOCUMENT 
     compute the rigid skeleton of a given skeleton:
     ie the set of segments connecting max to saddle points;
     EXAMPLE
     require,"split.i";
     u=fft_smooth(fft(genrandfield2D(128*2),-1),15).re;
     skl=rsex_make(NDfield_write(u,"test.ND"),voidpeak=1);
     skl2=NDskel2rskel(skl);
     limits,0.25,0.75,0.25,0.75;
   SEE ALSO:
 */
{
  p1=where(*skl.segdatainfo=="patch"+" p1")(1);
  p2=where(*skl.segdatainfo=="patch"+" p2")(1);

  d1=where(*skl.segdatainfo=="density"+" p1")(1);

  u=(*skl.segdata)(d1,);
  pp=(*skl.segdata)(p1,);
  
  seg_dif = (*skl.segdata)(p2,)-(*skl.segdata)(p1,);
  w2 = where (seg_dif != 0); 
  
  s=heapsort(pp);pp1=pp(s);
  un = 1+int(unique(pp1));
  w=where(pp1(dif)!=0);
  idxmax=array(0,1+int(max(pp)));
  w=grow(0,w);w=grow(w,numberof(pp));
  for(i=1;i<numberof(w);i++) {w1=s(w(i)+1:w(i+1));
    idxmax(un(i))=w1(u(w1)(*)(mxx));
  }
  x = (*skl.segpos)(,avg,);
  dat=(*skl.segdata)(4,idxmax(where(idxmax!=0)));
  segcol=1+long(199*normalise(dat));
  x1=x(,idxmax(where(idxmax!=0)));
  plcolor,segcol,x1(2,),x1(1,),npol=4,size=0.1;
  pl2,x1,width=3;
  dat=(*skl.segdata)(4,w2);
  segcol=1+long(199*normalise(dat));
  x1=x(,w2);
  plcolor,segcol,x1(2,),x1(1,),npol=10,size=0.1;
  XM=X1=X2=x(,w2);
  drp=array(1,dimsof(XM)(0));
  edgeid=NDskelonedge(skl);
  for(i=1;i<=numberof(w2);i++)
    {
      patchid=int((*skl.segdata)(p2,w2(i)));
      x1=x(,idxmax(1+patchid));x2=x(,w2(i));
      X1(,i)=x1;
      if (((patchid==edgeid)(sum))&&
            ((abs(x2-x1)>(*skl.delta)/2)(sum)))
        {
          drp(i)=0;
        }
      else
        pldj,x1(1),x1(2),x2(1),x2(2),width=2;
      
      patchid=int((*skl.segdata)(p1,w2(i)));
      x1=x(,idxmax(1+patchid));x2=x(,w2(i));
      X2(,i)=x1;
      if (((patchid==edgeid)(sum))&&
            ((abs(x2-x1)>(*skl.delta)/2)(sum)))
        {
          drp(i)=0;
        }
      else
        pldj,x1(1),x1(2),x2(1),x2(2),width=2;
      
    }
  drp=where(drp);          
  skl2=skl;
  skl2.comment="Rigid skeleton";
  skl2.nodepos=&(x(,idxmax(where(idxmax!=0))));
  skl2.nodedata=&((*skl.segdata)(,idxmax(where(idxmax!=0))));
  skl2.nnodes=numberof(where(idxmax!=0));
  XM=XM(,drp);  X1=X1(,drp);  X2=X2(,drp);
  skl2.segpos=&transpose([grow(XM,XM),grow(X1,X2)],[2,3]);
  skl2.nsegs=dimsof(*skl2.segpos)(0);
  return skl2;
}
  

  func plNDskel(skl,color=,noerase=,back=,last=,twoD=,threeD=,
                width=,xrange=,yrange=,zrange=,trange=,offset=,alt=,az=,
                slice=,delay=,fun=,win=,cut=,pal=,thresh=,debug=,withmax=,bar=,withnodes=)
/* DOCUMENT
   plots NDskel
   color= can be "density"  or "patch"
   if fun= e.g. log it will plot the log  e.g.density as a color table with
   color="density"
   SEE ALSO:
 */
  
{
  if(is_void(delay)) delay=300;
  if(is_void(bar)) bar=0;
  if(numberof(slice))
    { if(!is_void(twoD)) limits, (*skl.x0)(1,),
        (*skl.x0)(1,)+(*skl.delta)(1,), (*skl.x0)(2,), (*skl.x0)(2,)+(*skl.delta)(2,);
      if(threeD) limits, (*skl.x0)(1,),
        (*skl.x0)(1,)+(*skl.delta)(1,), (*skl.x0)(2,), (*skl.x0)(2,)+(*skl.delta)(2,);
      
      k=1;
      for(i=1;i<=slice;i++)
        {
          zrange=[(*skl.x0)(3)+(i-1)*1./slice*(*skl.delta)(3),(*skl.x0)(3)+i*1./slice*(*skl.delta)(3)];
          if(win) window,k++; else if(!noerase) fma;
          er=plNDskel(skl,color=color,noerase=noerase,back=back,last=last,
                   twoD=twoD,threeD=threeD,width=width,xrange=xrange,
                   yrange=yrange,zrange=zrange,trange=trange,offset=offset,
                   alt=alt,az=az,fun=fun,pal=pal);
          pause,delay;
        }
    }
  w=array(1,skl.nsegs);
  if(numberof(xrange))
  w*=((*skl.segpos)(1,1,)>xrange(1))*((*skl.segpos)(1,1,)<xrange(2));
  if(numberof(yrange))
  w*=((*skl.segpos)(2,1,)>yrange(1))*((*skl.segpos)(2,1,)<yrange(2));
  if(numberof(zrange))
  w*=((*skl.segpos)(3,1,)>zrange(1))*((*skl.segpos)(3,1,)<zrange(2));
  if(numberof(trange))
  w*=((*skl.segpos)(4,1,)>trange(1))*((*skl.segpos)(4,1,)<trange(2));
  if(thresh)
    {
      if(thresh>0) {
        cc=centile((*skl.segdata)(1,),lcut=thresh/2);
        w*=((*skl.segdata)(1,)> cc(2));
      }
      if(thresh<0) {
        cc=centile((*skl.segdata)(1,),lcut=-thresh/2);
        w*=((*skl.segdata)(1,)< cc(1));
      }
    }
  w=where(w); 
  if(!numberof(w)) return;
    
  W=array(1,skl.nnodes);
  if(numberof(xrange))
  W*=((*skl.nodepos)(1,)>xrange(1))*((*skl.nodepos)(1,)<xrange(2));
  if(numberof(yrange))
  W*=((*skl.nodepos)(2,)>yrange(1))*((*skl.nodepos)(2,)<yrange(2));
  if(numberof(zrange))
  W*=((*skl.nodepos)(3,)>zrange(1))*((*skl.nodepos)(3,)<zrange(2));
  if(numberof(trange))
  W*=((*skl.nodepos)(4,)>trange(1))*((*skl.nodepos)(4,)<trange(2));
  if(thresh)
    {
      if(thresh>0) {
        cc=centile((*skl.nodedata)(1,),lcut=thresh/2);
              W*=((*skl.nodedata)(1,)> cc(2));
      }
      if(thresh<0)  {
        cc=centile((*skl.nodedata)(1,),lcut=-thresh/2);
        W*=((*skl.nodedata)(1,)< cc(1));
      }
    }
  W=where(W);
  //if(!numberof(W)) return;

  //  if(is_void(width)) width=(*skl.delta)(avg)/2e2;

  width0=width;
  /*
  if(is_string(color))  width0=  2*median(abs((*skl.segpos)(1,1,)(dif)(2:-1)));
  else width0=width;
  
  if(is_void(width))  width=width0; else width=width0*width;
  */
  
  if(is_void(alt))   alt=pi/4;
  if(is_void(az)) az=pi/5;

  if(offset) off=-1; else off=0;
  if(offset) fac=2; else fac=1;
  if(skl.ndims==2) twoD=[1];
  if((threeD))
    {
      if(skl.ndims<3) error,"can't plot 2D skeleton in 3D you fool !";
      x=(*skl.segpos*fac+off)(1,1,w);
      y=(*skl.segpos*fac+off)(2,1,w);
      z=(*skl.segpos*fac+off)(3,1,w);
      yp1= x * (c= cos(az-pi/2)) - y * (s= sin(az-pi/2));
      xp1= z * cos(alt-pi/2) + (x * s + y * c) * sin(alt-pi/2);
      x=(*skl.segpos*fac+off)(1,2,w);
      y=(*skl.segpos*fac+off)(2,2,w);
      z=(*skl.segpos*fac+off)(3,2,w);
      yp2= x * (c= cos(az-pi/2)) - y * (s= sin(az-pi/2));
      xp2= z * cos(alt-pi/2) + (x * s + y * c) * sin(alt-pi/2);
      if(is_string(color))
        {
          w1=where(*skl.segdatainfo==color+" p1")(1);
          w2=where(*skl.segdatainfo==color+" p2")(1);
          if(is_func(fun)) dat=fun((*skl.segdata)([w1,w2],w)(1,)); else dat=(*skl.segdata)([w1,w2],w)(1,);
          segcol=1+long(239*normalise(dat));
          pldjc,Colors(240)(segcol),xp1,yp1,xp2,yp2,width=width;
        } else
        pldj,xp1,yp1,xp2,yp2,width=width,color=color,type=type;
    }
  else if(!is_void(twoD))
    {
      
      if(is_string(color))
        {
          w1=where(*skl.segdatainfo==color+" p1");
          w2=where(*skl.segdatainfo==color+" p2");
          if ((numberof(w1)==0)||(numberof(w2)==0))
          {
              w1=where(*skl.segdatainfo==color)(1);
              w2=w1;
          }
          else {
              w1=w1(1);
              w2=w2(1);
          }
          w1;w2;
          if(is_func(fun)) dat=fun((*skl.segdata)([w1,w2],w)(2,)); else
            dat=(*skl.segdata)([w1,w2],w)(2,);
          mx = max(dat);mn=min(dat);
          segcol=1+long(199*(dat-mn)/(mx-mn));
          
               if(withmax)
                 {
                   u1=(*skl.segdata)(1,);
                   p1=(*skl.segdata)(3,);
                   jj=NDgetMinMax(u1,p1);
                   if(is_func(fun)) dat=fun(u1(jj(,2))); else  dat=u1(jj(,2));
                   if(color=="density")   nodecol=1+long(199*(dat-mn)/(mx-mn))(::1); else nodecol=1+long(199*normalise(dat))(::1);
                   if(is_string(pal)) col1=Colors3(200,pal=pal)(,nodecol); else
                     col1=Colors(200)(::-1)(nodecol);
        plcolor,col1,(*skl.segpos)(2,2,jj(,2)),(*skl.segpos)(1,2,jj(,2)),size=width*35,npol=4,bar=bar;
                   plp,(*skl.segpos)(2,2,jj(,2)),(*skl.segpos)(1,2,jj(,2)),symbol=8,fill=1;

                 } else if (numberof(W))
                 {
                   w1=where(*skl.nodedatainfo==color);
                   if (numberof(w1))
                   {
                       w1=w1(1);
                       if(is_func(fun))
                           dat=fun((*skl.nodedata)(w1,W));
                       else dat=(*skl.nodedata)(w1,W);
                   } else
                   {
                       dat=indgen(numberof(W));
                   }
                   mn=min(dat);
                   mx=max(dat);
                   nodecol=1+long(199*(dat-mn)/(mx-mn))(::1);
                   if(is_string(pal)) col1=Colors3(200,pal=pal)(,nodecol); else
                       col1=Colors(200)(::-1)(nodecol);
                   
                   //plcolor,col1,(*skl.nodepos)(2,W),(*skl.nodepos)(1,W),size=width*25,npol=8,bar=bar;
                   plcolor,col1,(*skl.nodepos)(2,W),(*skl.nodepos)(1,W),size=0.01,npol=8,bar=bar;
                   
                 }
               // error,"";
           if(is_string(pal)) col1=Colors3(200,pal=pal)(,segcol); else
             col1=Colors(200)(::-1)(segcol);
/*
          pldjc,col1,(*skl.segpos*fac+off)(1,1,w),
            (*skl.segpos*fac+off)(2,1,w),(*skl.segpos*fac+off)(1,2,w),
            (*skl.segpos*fac+off)(2,2,w),width=width;*/
           
           y=(*skl.segpos*fac+off)(2,,w);
           x=(*skl.segpos*fac+off)(1,,w);
           c=col1(w);
           for (i=1;i<=numberof(c);i++)
               plg,y(,i),x(,i),color=c(i),width=width;
                   
           
          if(debug&&(numberof(W))) plt1,pr01(long(dat)),(*skl.nodepos)(1,W),(*skl.nodepos)(2,W),tosys=1;
          
        }
          else
          {
              
              //c=col1(w);
              if ((withnodes)&&(numberof(W))) {
                t=((*skl.node).type)(W);incolor=[__blue,__green,__red,__orange];stat,t;
                  for (i=min(t);i<=max(t);i++)
                  {
                      wt = where(t==i);
                      if (numberof(wt))
                        {
                        if (skl.ndims!=2)
                           {
                             y=(*skl.nodepos)(3,W(wt));
                             x=(*skl.nodepos)(2,W(wt));
                           }
                         else
                           {
                             y=(*skl.nodepos)(2,W(wt));
                             x=(*skl.nodepos)(1,W(wt));
                           }
                        
                        plmk2,y,x,incolor=incolor(,i+1),width=11,msize=0.3,marker=i+1;
                        }
                  }
              }

              if (skl.ndims!=2)
                {
                  if (numberof(twoD)!=2)
                    {
                      y=(*skl.segpos*fac+off)(3,,w);
                      x=(*skl.segpos*fac+off)(2,,w);
                    }
                  else
                    {
                      y=(*skl.segpos*fac+off)(twoD(2),,w);
                      x=(*skl.segpos*fac+off)(twoD(1),,w);
                    }
                }
              else {
              y=(*skl.segpos*fac+off)(2,,w);
              x=(*skl.segpos*fac+off)(1,,w);
              }
             
              for (i=1;i<=numberof(w);i++)
                  plg,y(,i),x(,i),color=color,width=width;

              
              
              /*
              if (numberof(W) && withnodes)
              {
                  col1=array(char(color/2),skl.nnodes);
                  mywidth=2*median(abs((*skl.segpos)(1,1,)(dif)(2:-1)));
                  plcolor,col1,(*skl.nodepos)(2,W),(*skl.nodepos)(1,W),size=mywidth/1000,npol=8,bar=bar;
              }
              pldj,(*skl.segpos*fac+off)(1,1,w),
                   (*skl.segpos*fac+off)(2,1,w),(*skl.segpos*fac+off)(1,2,w),
                   (*skl.segpos*fac+off)(2,2,w),color=color,width=width,type=type;
              */
          }
      /*
      if(skl.ndims>=3)
          if(!is_string(color))
            {
              segcol=1+long(199*normalise((*skl.segpos*fac+off)(skl.ndims,2,w)));
          if(is_string(pal)) col1=Colors3(200,pal=pal)(,segcol); else col1=Colors(200)(segcol)(::-1); 
              pldjc,col1,(*skl.segpos*fac+off)(1,1,w),(*skl.segpos*fac+off)(2,1,w),(*skl.segpos*fac+off)(1,2,w),(*skl.segpos*fac+off)(2,2,w),width=width;
            }
      */
    }
  else if(skl.ndims==3)
    {
      if(is_void(color)) color=__red/256.;
      if(!noerase) clear3d;
      if(!numberof(back))  back=[0.0,0.0,0.0]; back_rgb3d,back;
      draw3d;
      if(is_string(color))
        {
          w1=where(*skl.segdatainfo==color+" p1")(1);
          w2=where(*skl.segdatainfo==color+" p2")(1);
          segcol=1+long(199*normalise((*skl.segdata)([w1,w2],w)(1,)));
          if(is_string(pal)) fma;
          for(i=1;i<=numberof(w);i++) lines3d,2,((*skl.segpos*fac+off)(:3,,w(i))),Colors3(200,pal=pal)(,segcol(i))/200.;
          inc_seq3d;
          draw3d_trigger;
        } else
          { 
          for(i=1;i<=numberof(w);i++) lines3d,2,((*skl.segpos*fac+off)(:3,,w(i))),color;
          inc_seq3d;
          draw3d_trigger;
        }
    }
  else if(skl.ndims==4)
    {
      if(!noerase) clear3d;
      if(is_void(last)) last=4;
      rg=where(indgen(4)!=last);
      if(!numberof(back))  back=[0.0,0.0,0.0]; back_rgb3d,back;
          if(is_string(pal)) fma;
      for(i=1;i<=numberof(w);i++) lines3d,2,((*skl.segpos*fac+off)(rg,,w(i))),Colors3(200,pal=pal)(,1+long(199*(*skl.segpos*fac+off)(last,1,w(i))))/200.;
      inc_seq3d;
      draw3d_trigger;
    }
  if(cut&&(numberof(W))) {
    skl1=skl;
    skl1.nsegs=numberof(w);
    skl1.nnodes=numberof(W);
    skl1.segpos=&((*skl.segpos*fac+off)(,,w));
    skl1.segdata=&((*skl.segdata)(,w));
    skl1.nodedata=&((*skl.nodedata)(,W));
    skl1.nodepos=&((*skl.nodepos*fac+off)(,W));
    return skl1;
  } else return fname1;
}


func NDskel_info(void)
/* DOCUMENT 
     numberof of filaments connected to node i
info,(*skl.node)(i).nnext

numberof of segments in the jth filament of node i
(*(*skl.node)(i).nsegs)(j)

index of the node to which node i is connected by filament j
(*(*skl.node)(i).next)(j)
i.e. if I start at node i and follow filament j I arrive at node
(*(*skl.node)(i).next)(j)

index of the first segment of the jth filament of node i
(*(*skl.node)(i).seg)(j)

index of the two nodes which are connected by the filament containing the segment i
 (*skl.seg)(i).nodes 

*skl.segdatainfo contains
["density p1","density p2","patch p1","patch p2","size","dist_to_node",
"dist_from_node","dist_to_saddle"]
where, for segment i, dist_to_node contains the curvilinear distance to node (*skl.seg)(i).nodes(2) 

if this condition is satisfied
 (*(*skl.node)(i).seg)(j).nodes(1)   == i 
then it means I am moving away from my starting point.

if a segment next or prev is zero then you are at the end or the beginning of
the filament. For instance:
numberof(where((*skl.seg).next==0))
 is equal to  (*skl.node).nnext(sum)/2

to find saddle points
w = where((*skl.seg).flags==1)
plg2,(*skl.segpos)(2,avg,w),(*skl.segpos)(1,avg,w),msize=8,color=_black

   SEE ALSO:
 */
{}


func NDskel_extract(skl,id,recenter=)
/* DOCUMENT 
     extract filament which is centered on node id;
     !!!!! FIXME:      needs to fix  nodes 
     SEE ALSO:
 */
{
  idx=NDfilaments(skl,id,0);
  sdata=ndata=spos=npos=[];
    grow,ndata,(*skl.nodedata)(,id);
    grow,npos, (*skl.nodepos)(,id);
  for(i=1;i<=numberof(idx);i++) {
    ii=*(idx)(i);
    grow,sdata,(*skl.segdata)(,ii);
    grow,spos, (*skl.segpos)(,,ii);
  }
  skl1=NDSkelStruct();
  if(numberof(recenter)){center=spos(,,avg); spos=(spos+center) % 1.;}
  skl1.segdata=&sdata;
  skl1.segpos=&spos;
  skl1.nodedata=&ndata;
  skl1.nodepos=&npos;
  skl1.comment="Extracted "+pr01(id,ndigit=5);
  skl1.nsegs=dimsof(sdata)(0);
  skl1.nnodes=dimsof(ndata)(0);
  skl1.x0=skl.x0;
  skl1.ndims=skl.ndims;
  skl1.dims=skl.dims;
  skl1.nodedatainfo=skl.nodedatainfo;
  skl1.segdatainfo=skl.segdatainfo;
  skl1.delta=skl.delta;
  skl1.nsegdata=skl.nsegdata;
  skl1.nnodes=1;
  skl1.nnodedata=skl.nnodedata;
  return skl1;
}

//dir can be 0 for any direction, 1/-1 for postive/negative direction
func NDfilaments(skl,curnode,dir)
/* DOCUMENT 

*NDfilaments(skl,1,0)(1)
index of all segments which start at node 1 and belong to filament 1
whatever the orientation
*NDfilaments(skl,1,1)(1)
same but only downhill orientation (w.r.t height of maximum)
*NDfilaments(skl,1,-1)(1)
same but only uphill orientation;  (not totally functinonal now)

e.g
idx=NDfilaments(skl,1,0);
w1=where(*skl.segdatainfo=="density p1")(1);
for(i=1;i<=numberof(idx);i++) {ii=*(idx)(i); pldjc,(*skl.segdata)(w1,ii),(*skl.segpos)(1,1,ii),(*skl.segpos)(2,1,ii),(*skl.segpos)(1,2,ii),(*skl.segpos)(2,2,ii);}

// bug: 1 seg long filaments are always positive direction

SEE ALSO:
 */
{
  node=(*skl.node)(curnode);
  seg =*node.seg;
  nsegs = *node.nsegs;
  if (node.nnext==0) return [];
  ptr = array(pointer,node.nnext);
  
  index=1;

  for (i=1;i<=node.nnext;i++)
    {
      if (dir>=0)
        {
          myseg = seg(i);
          
          if (((*skl.seg)(myseg).next!=0)||(nsegs(i)==1))
            {
              arr = array(int,nsegs(i));
              arr(1)=myseg;
              for (j=2;j<=nsegs(i);j++)
                {
                  myseg = (*skl.seg)(myseg).next;
                  arr(j) = myseg;
                }
              ptr(index) = &arr;
              index++;
            }
          
        }
      if (dir<=0)
        {
          myseg = seg(i);
          if (((*skl.seg)(myseg).prev!=0))
            {
              arr = array(int,nsegs(i));
              arr(1)=myseg;
              for (j=2;j<=nsegs(i);j++)
                {
                  myseg = (*skl.seg)(myseg).prev;
                  arr(j) = myseg;
                }
              ptr(index) = &arr;
              index++;
            }
        }
    }
  
  return ptr(1:index-1);
}



func skel_orient(skl,boxsize=,x0=,npix=)
/* DOCUMENT 
    skl2=skel_smooth(skl,boxsize=1e4);
    oo=skel_orient(skl2,boxsize=2e4);
    plii,(oo),span(-1,1,2),span(-1,1,2);
    skl2=skel_smooth(skl,boxsize=1e4);n=_len(skl2); for(i=1;i<=n;i++){x=_nxt(skl2);if(dimsof(x)(0)>1) pl2,x/1e4-1,color=-5,width=0.3;}
     
   SEE ALSO:
 */
{
  if(is_void(npix)) npix=64;
  if(is_void(x0)) x0=0;
  if(is_void(boxsize)) boxsize=1.;
  
  n=_len(skl);n;
  u=array(0.,[2,npix,npix]);
  for(i=1;i<=n;i++)
    {
      x=_nxt(skl);
     if(dimsof(x)(0)>1)u+=histo2d(transpose(x(1:2,)),span(x0,x0+boxsize,npix+1),span(x0,x0+boxsize,npix+1),wght=abs(x(3,dif)(pcen))*dimsof(x)(0));
    }
  return u;
}



//if(is_void(u)) {u=fitsRead("/raid/pichon/SIMU-512-stephane/40Mpc/Redshift2/LambdaCDM_512x512_40Mpc_h0.70_Ob0v04.nHIN04.y001.fits") u=log(u); u=smooth(u,15);u=u(:256,:256);}

//uz=genrandfield3D(128); u=fft(uz,[-1,-1,-1]).re;u=smooth(u,5);





func hess3D(field,&s,ft=)
/* DOCUMENT hessian in 3D
     
   SEE ALSO:
 */
{
  local dxx,dyy,dzz,dxy,dxz,dyz,hess;
  if(ft)
    {
      d=dimsof(field);
      nx=d(2);
      ny=d(3);
      nz=d(4);
      kx=2*pi*fft_indgenNyquist(nx)(,-:1:ny,-:1:nz)/nx;
      ky=2*pi*fft_indgenNyquist(ny)(-:1:nx,,-:1:nz)/ny;
      kz=2*pi*fft_indgenNyquist(nz)(-:1:nx,-:1:ny,)/nz;
      uk=fft(field,1);
      H= fft(-[[uk*kx*kx,uk*ky*kx,uk*kz*kx], // Hessian by fft
           [uk*kx*ky,uk*ky*ky,uk*kz*ky],
           [uk*kx*kz,uk*ky*kz,uk*kz*kz]],[-1,-1,-1,0,0]).re/numberof(field);
      return H;
    }
  dxx=  field(dif,,)(pcen,,)(dif,,)(pcen,,);
 dyy=  field(,dif,)(,pcen,)(,dif,)(,pcen,);
 dzz=  field(,,dif)(,,pcen)(,,dif)(,,pcen);
 dxy=  0.5*field(dif,,)(pcen,,)(,dif,)(,pcen,)+
   0.5*field(,dif,)(,pcen,)(dif,,)(pcen,,);
 dxz=  0.5*field(dif,,)(pcen,,)(,,dif)(,,pcen)+
   0.5*field(,,dif)(,,pcen)(dif,,)(pcen,,);
 dyz=  0.5*field(,dif,)(,pcen,)(,,dif)(,,pcen)+
   0.5*field(,,dif)(,,pcen)(,dif,)(,pcen,);
 hess=[[dxx,dxy,dxz],[dxy,dyy,dyz],[dxz,dyz,dzz]];
 //  hess=(hess+transpose(hess,[4,5]))/2.;
  // s=SVdec(hess,U);
 return hess;
}





func hess2D(field,&s,ft=)
/* DOCUMENT hessian in 2D
     
   SEE ALSO:
 */
{
  if(ft)
    {
      d=dimsof(field);
      nx=d(2);
      ny=d(3);
      kx=2*pi/(nx+1)*fft_indgenNyquist(nx)(,-:1:ny);
      ky=2*pi/(ny+1)*fft_indgenNyquist(ny)(-:1:nx,);
      uk=fft(field,1);
      H= fft(-[[uk*kx*kx,uk*ky*kx], // Hessian by fft
           [uk*kx*ky,uk*ky*ky]],[-1,-1,0,0]).re/numberof(field);
      return H;
    }
  local hess,dxx,dyy,dxy;
 dxx=  field(dif,)(pcen,)(dif,)(pcen,);
 dyy=  field(,dif)(,pcen)(,dif)(,pcen);
 dxy=  0.5*field(dif,)(pcen,)(,dif)(,pcen)+
   0.5*field(,dif)(,pcen)(dif,)(pcen,);
 hess=[[dxx,dxy],[dxy,dyy]];
 hess=(hess+transpose(hess,[3,4]))/2.;
 // s=SVdec(hess,U);
 return hess;
}


func grad3D(field,ft=)
/* DOCUMENT gradiant
     
   SEE ALSO:
 */
{
  local dx,dy,dz,grad;
  if(ft)

    {
      d=dimsof(field);
      nx=d(2);
      ny=d(3);
      nz=d(4);
      kx=2*pi*fft_indgenNyquist(nx)(,-:1:ny,-:1:nz)/(nx+1);
      ky=2*pi*fft_indgenNyquist(ny)(-:1:nx,,-:1:nz)/(ny+1);
      kz=2*pi*fft_indgenNyquist(nz)(-:1:nx,-:1:ny,)/(nz+1);
      return fft(1.i*[kx,ky,kz]*fft(field,1)(,,,-),[-1,-1,-1,0]).re/numberof(field);
    }
  dx=  field(dif,,)(pcen,,);
  dy=  field(,dif,)(,pcen,);
  dz=  field(,,dif)(,,pcen);
 grad=[dx,dy,dz];
 return grad;
}


func grad2D(field,ft=)
/* DOCUMENT gradiant in 2D
     
   SEE ALSO:
 */
{
  local dx,dy,grad;
  if(ft)

    {
      d=dimsof(field);
      nx=d(2);
      ny=d(3);
      kx=2*pi/(nx+1)*fft_indgenNyquist(nx)(,-:1:ny);
      ky=2*pi/(ny+1)*fft_indgenNyquist(ny)(-:1:nx,);
      return fft(1.i*[kx,ky]*fft(field,1)(,,-),[-1,-1,0]).re/numberof(field);
    }

  dx=  field(dif,)(pcen,);
  dy=  field(,dif)(,pcen);
 grad=[dx,dy];
 return grad;
}


func cross(x,y)
/* DOCUMENT 
     crossproduct of 2 vectors in 2/3D
   SEE ALSO:
 */
{
  if(dimsof(x)(0)==3)    return [-x(..,3)*y(..,2)+x(..,2)*y(..,3),
    x(..,3)*y(..,1)-x(..,1)*y(..,3),
    x(..,1)*y(..,2)-x(..,2)*y(..,1)];
  else return  x(..,1)*y(..,2)-x(..,2)*y(..,1);
}

func dot(A,x)
/* DOCUMENT 
     dotproduct of 2 vectors in 2/3D
   SEE ALSO:
 */
{

  if(dimsof(x)(0)==3)
     return [A(..,1,1)*x(..,1)+A(..,1,2)*x(..,2)+A(..,1,3)*x(..,3),
           A(..,2,1)*x(..,1)+A(..,2,2)*x(..,2)+A(..,2,3)*x(..,3),
           A(..,3,1)*x(..,1)+A(..,3,2)*x(..,2)+A(..,3,3)*x(..,3)];
  else return [A(..,1,1)*x(..,1)+A(..,1,2)*x(..,2),
           A(..,2,1)*x(..,1)+A(..,2,2)*x(..,2)];
    
}



func crit2D(u,ft=)
/* DOCUMENT 
returns cross[Hessian.grad,grad]
   
   SEE ALSO:
 */
{
  gd=grad2D(u,ft=ft);
  hs=hess2D(u,ft=ft);
  tt1=dot(hs,gd);
  tt=cross(gd,tt1);
  //  tt /= max(1e-15,sqrt((gd^2)(..,sum)));
  //  tt /= max(1e-15,sqrt((tt1^2)(..,sum)));
  return tt;
}

func crit3D(u)
/* DOCUMENT 
returns cross[Hessian.grad,grad]
   
   SEE ALSO:
 */
{
  gd=grad3D(u);
  hs=hess3D(u);
  tt1=dot(hs,gd);
  tt=cross(gd,tt1);
  //  tt /= max(1e-15,sqrt((gd^2)(..,sum)));
  //  tt /= max(1e-15,sqrt((tt1^2)(..,sum)));
  return tt;
}
