
struct NDnetworkData
{
    int type; //0 for vertex, t for faces of type (e.g. dimension) t.
   string name;
    pointer data;
};

NDNETWORKMAXDIMS=6;
    
struct NDnetwork
{
  string comment;
  int periodicity;
  int ndims;
  int ndims_net;
  int f_areSimplex;
  pointer x0;
  pointer delta;
  int indexSize;
  int cumIndexSize;
  int floatSize;

  int nvertex;
  pointer v_coord;
  int nfaces(NDNETWORKMAXDIMS);
    
  int haveVertexFromFace(NDNETWORKMAXDIMS);
  pointer f_numVertexIndexCum(NDNETWORKMAXDIMS);
  pointer f_vertexIndex(NDNETWORKMAXDIMS);

  int haveFaceFromVertex(NDNETWORKMAXDIMS);
  pointer v_numFaceIndexCum(NDNETWORKMAXDIMS);
  pointer v_faceIndex(NDNETWORKMAXDIMS);
    
  int haveFaceFromFace(NDNETWORKMAXDIMS,NDNETWORKMAXDIMS);
  pointer f_numFaceIndexCum(NDNETWORKMAXDIMS,NDNETWORKMAXDIMS);
  pointer f_faceIndex(NDNETWORKMAXDIMS,NDNETWORKMAXDIMS);

  int haveVFlags;
  int haveFFlags(NDNETWORKMAXDIMS);
  pointer v_flag;
  pointer f_flag(NDNETWORKMAXDIMS);

  int ndata;
  pointer data;

  char dummy(160-3*4);
} ;



func NDnetwork_ascii_write(net,fname)
{
  ff=open(fname,"w");
  write,ff,format="%s\n","ANDNET";
  write,ff,format="%d\n",net.ndims;
  //BBOX
  write,ff,format="%ld\n",net.nvertex;
  p=*net.v_coord;
  fmt_="%g";
  fmt=fmt_;
  for (i=2;i<=net.ndims;i++) fmt+=" "+fmt_;
  fmt+="\n";//fmt;
  
  if (net.ndims==2) write,ff,format=fmt,p(1,),p(2,);
  if (net.ndims==3) write,ff,format=fmt,p(1,),p(2,),p(3,);
  
  for (i=1;i<=net.ndims;i++)
    {
      if (net.f_vertexIndex(i)==pointer()) continue;
      id=*net.f_vertexIndex(i)-1;
      write,ff,format="%ld %ld\n",i,net.nfaces(i);
      if (i==1) write,ff,format="%ld %ld\n",id(1,),id(2,);
      if (i==2) write,ff,format="%ld %ld %ld\n",id(1,),id(2,),id(3,);
      if (i==3) write,ff,format="%ld %ld %ld %ld\n",id(1,),id(2,),id(3,),id(4,);
    }
  
  write,ff,format="%s\n","[ADDITIONAL_DATA]";
  data=*net.data;
  for (i=1;i<=numberof(data);i++)
    {
      write,ff,format="%s\n",data(i).name;
      write,ff,format="%ld\n",data(i).type;
      write,ff,format="%g\n",*data(i).data;
    }
  
  close,ff;
  return fname;
}

func NDnetwork_write(net,fname)
{
  fname2=fname+".ascii_TMP";
  NDnetwork_ascii_write(net,fname2);
  cmd="netconv "+fname2+" -to NDnet -outName "+fname;
  exec,cmd;
}

func NDnetwork_read(fname,headeronly=)
/* DOCUMENT 


SEE ALSO:
 */
{
  if(!is_void(last)) error,"not implemented yet !";
  if(!is_void(first)) error,"not implemented yet !";
  if (is_void(fname)) error," NDfield_read needs an argument";

  stream = _grid_open(fname);
  address = 0;
  dummy = _grid_read(int);
  if (dummy!=16)
    {
      close, stream;
      stream = _grid_open(fname,encoding = "sun" );
      address=0;
      dummy = _grid_read(int);
    } 
  tag =_grid_read(char, 16);
  if(strchar(tag)(1)!="NDNETWORK") error,"Sorry, this is not an NDnetwork file ...";
  dummy = _grid_read(int);

  net=NDnetwork();

  dummy = _grid_read(int);
  net.ndims = _grid_read(int);
  net.ndims_net = _grid_read(int);
  dummy = _grid_read(int);
  
  dummy = _grid_read(int);
  comment =_grid_read(char, 80);
  net.comment=strchar(comment)(1);
  net.periodicity=_grid_read(int);
  net.f_areSimplex=_grid_read(int);
  net.x0 = &_grid_read(double,net.ndims);
  net.delta =&_grid_read(double,net.ndims);
  net.indexSize =_grid_read(int);
  net.cumIndexSize =_grid_read(int);
  net.floatSize =_grid_read(int);
  net.dummy =_grid_read(char,160-3*4);
  net.nvertex = _grid_read(int);
  dummy = _grid_read(int);
  
  if(headeronly) return net;
  dummy = _grid_read(int);
  if (net.floatSize==8)
    net.v_coord = &_grid_read(double,[2,net.ndims,net.nvertex]);
  else
    net.v_coord = &_grid_read(float,[2,net.ndims,net.nvertex]);
  dummy = _grid_read(int);

  dummy = _grid_read(int);
  net.nfaces(1:net.ndims) = _grid_read(int,1+net.ndims)(2:);
  dummy = _grid_read(int);

  dummy = _grid_read(int);
  net.haveVertexFromFace(1:net.ndims) = _grid_read(int,1+net.ndims)(2:);
  dummy = _grid_read(int);
  
  for (i=1;i<=net.ndims;i++)
  {
      if (net.haveVertexFromFace(i))
      {
          if (!net.f_areSimplex)
          {
              dummy = _grid_read(int);
              if (net.cumIndexSize==8)
                net.f_numVertexIndexCum(i)=&(1+_grid_read(long,1+net->nfaces(i)));
              else
                 net.f_numVertexIndexCum(i)=&(1+_grid_read(int,1+net->nfaces(i)));
              dummy = _grid_read(int);
              
              dummy = _grid_read(int);
              if (net.indexSize==8)
                net.f_vertexIndex(i)=&(1+_grid_read(long,(*net.f_numVertexIndexCum(i))(1+net.nfaces(i))));
              else
                net.f_vertexIndex(i)=&(1+_grid_read(long,(*net.f_numVertexIndexCum(i))(1+net.nfaces(i))));
              dummy = _grid_read(int);
          }
          else
          {
              dummy = _grid_read(int);
              if (net.indexSize==8)
                net.f_vertexIndex(i)=&(1+_grid_read(long,[2,i+1,net.nfaces(i)]));
              else
                net.f_vertexIndex(i)=&(1+_grid_read(int,[2,i+1,net.nfaces(i)]));
              dummy = _grid_read(int);
          }
      }
  }

  dummy = _grid_read(int);
  net.haveFaceFromVertex(1:net.ndims) = _grid_read(int,1+net.ndims)(2:);
  dummy = _grid_read(int);
  
  for (i=1;i<=net.ndims;i++)
  {
      if (net.haveFaceFromVertex(i))
      {
          dummy = _grid_read(int);
          if (net.cumIndexSize==8)
            net.v_numFaceIndexCum(i)=&(1+_grid_read(long,net.nvertex+1));
          else
            net.v_numFaceIndexCum(i)=&(1+_grid_read(int,net.nvertex+1));
          dummy = _grid_read(int);
          
          dummy = _grid_read(int);
          if (net.cumIndexSize==8)
            net.v_faceIndex(i)=&(1+_grid_read(long,(*net.v_numFaceIndexCum(i))(1+net.nvertex)-1));
          else
            net.v_faceIndex(i)=&(1+_grid_read(int,(*net.v_numFaceIndexCum(i))(1+net.nvertex)-1));
          dummy = _grid_read(int);
      }
  }

  dummy = _grid_read(int);
  net.haveFaceFromFace(1:net.ndims,1:net.ndims) = _grid_read(int,[2,1+net.ndims,1+net.ndims])(2:,2:);
  dummy = _grid_read(int);
  
  for (i=1;i<=net.ndims;i++)
  {
     for (j=1;j<=net.ndims;j++)
     { 
         if (net.haveFaceFromFace(i,j))
         {
             dummy = _grid_read(int);
             if (net.cumIndexSize==8)
               net.f_numFaceIndexCum(i,j)=&(1+_grid_read(long,net.nfaces(i)));
             else
               net.f_numFaceIndexCum(i,j)=&(1+_grid_read(int,net.nfaces(i)));
             dummy = _grid_read(int);
             
             dummy = _grid_read(int);
             if (net.cumIndexSize==8)
               net.f_faceIndex(i,j)=&(1+_grid_read(long,(*net.f_numFaceIndexCum(i,j))(1+net.nfaces(i))-1));
             else
               net.f_faceIndex(i,j)=&(1+_grid_read(int,(*net.f_numFaceIndexCum(i,j))(1+net.nfaces(i))-1));
             dummy = _grid_read(int);
         }
     }
  }

  dummy = _grid_read(int);
  net.haveVFlags = _grid_read(int,1);
  dummy = _grid_read(int);

  if (net.haveVFlags)
  {
      dummy = _grid_read(int);
      net.v_flag=&(_grid_read(char,net.nvertex));
      dummy = _grid_read(int);
  }

  dummy = _grid_read(int);
  net.haveFFlags(1:net.ndims) = _grid_read(int,net.ndims+1)(2:);
  dummy = _grid_read(int);

  for (i=1;i<=net.ndims;i++)
      if (net.haveFFlags(i))
          {
              dummy = _grid_read(int);
              net.f_flag(i)=&(_grid_read(char,net.nfaces(i)));
              dummy = _grid_read(int);
          }

  dummy = _grid_read(int);
  net.ndata = _grid_read(int,1);
  dummy = _grid_read(int);

  if (net.ndata) net.data = &array(NDnetworkData,net.ndata);
  
  for (i=1;i<=net.ndata;i++)
  {
      dummy = _grid_read(int);
      (*net.data)(i).type = _grid_read(int,1);
      (*net.data)(i).name = strchar(_grid_read(char,255))(1);
      dummy = _grid_read(int);

      if ((*net.data)(i).type==0) N=net.nvertex;
      else N=net.nfaces((*net.data)(i).type);
      
      dummy = _grid_read(int);
      (*net.data)(i).data = &_grid_read(double,N);
      dummy = _grid_read(int,1);
  }
  
  close,stream;

  
  
  return net;
}
/*
  WS,0;
  v=(*net.f_vertexIndex(2))(,1+12819)
  nei=(*net.f_vertexIndex(2))(,1+[39040,12818,80913])
  

  //v=(*net.f_vertexIndex(2))(,12819)
  //index=(*net.v_numFaceIndexCum(2));
  //f=[];
  //for (i=1;i<=numberof(v);i++)
  //  grow,f,(*net.v_faceIndex(2))(index(v(i)):index(v(i)+1)-1);
  //nei=(*net.f_vertexIndex(2))(,f)

  
  pos=*net.v_coord;
  PL,pos(2,v),pos(1,v),msize=1,marker=1;
  PL,pos(2,nei(*)),pos(1,nei(*)),msize=1,marker=2,incolor=__blue;
 */

func NDNetDataIndex(net,type,txt)
{
    w=where((*net.data).type == type);
    if (numberof(w)==0) return -1;
    
    w=w(where(((*net.data).name)(w) == txt));
    if (numberof(w)==1) return w(*)(1);

    return -1;
}

func plNDnetPatches(net,edges=,peak=)
{
    c_c;
    if (peak)
        pi=NDNetDataIndex(net,net.ndims,"peak patch index");
    else
        pi=NDNetDataIndex(net,net.ndims,"void patch index");
    
    patch = *((*net.data)(pi).data);
    patch = double(sort(random(numberof(patch)))) (1+int(patch));
    c=bytscl(patch)+char(1);
    palette, red, green, blue, query=1;
    c=[red(c),green(c),blue(c)];
    //p=(patch-min(patch))/(max(patch)-min(patch));
    //c=char(255.*[p,abs(p-0.5),1.-p]);
    c=array(c,[2,2,2]);
    flist=*(net.f_vertexIndex(net.ndims));
    pos=*net.v_coord;
    x=pos(1,);y=pos(2,);
    p1=flist(1,);
    p2=flist(2,);
    p3=flist(3,);
    
    for (i=1;i<=net.nfaces(net.ndims);i++)
        plf,c(i,,,),[[y(p1(i)),y(p2(i))],[y(p3(i)),y(p3(i))]],[[x(p1(i)),x(p2(i))],[x(p3(i)),x(p3(i))]],edges=edges;
}

func plotFOF(net,color=,incolor=,randcol=,msize=,peak=)
{
    if (!is_void(peak))
        id=NDNetDataIndex(net,0,"Peak group ID");
    else
        id=NDNetDataIndex(net,0,"Void group ID");
    
    val = int(*((*net.data)(id).data));
    w=where(val);numberof(w);
    pos=*net.v_coord;
    if (!is_void(randcol))
    {
        col=char(random(val(max)+1,3)*255);
        for (i=1;i<=numberof(w);i++)
            PL,pos(2,w(i)),pos(1,w(i)),color=col(val(w(i)),),msize=msize,incolor=col(val(w(i)),);
    }
    else PL,pos(2,w),pos(1,w),color=color,msize=msize,incolor=incolor;
}


func remap(in,&result,&ct,wei=,srt=)
{
  // assign an integer group index to each value in 'in'.
  // index ranges from 1 to numberof of different values in 'in'.
  // ct is the number of identical values in each group
  // result is the weighted number of identical values in each group
  
  if (is_void(srt)) srt=sort(in);
  
  tmp=in(srt)(dif); 
  w=where(tmp);
  //simply count how many of each group
  ct=array(long(0),numberof(w)+1);
  ct(1)=w(1);
  ct(2:numberof(w))=w(dif);
  ct(0)=numberof(in)-ct(sum);
   
  out=long(in);
  out(srt)=long((tmp!=0)(cum)+1);
  
  if (!is_void(wei))
    {
      // weighted count of each group
      dm=dimsof(wei);
      dm(0)=numberof(w)+1;
      result=array(double(0),dm);
      for (i=1;i<=dimsof(wei)(0);i++) result(..,out(i))+=wei(..,i);
    }
  
  return out;
}


func myplcol(val,y,x,pal=,marker=,width=,type=,msize=)
/* DOCUMENT plcol(c,y,x,pal=,marker=,width=,type=)
    displays a cluster of points (x,y) colored by vector c with marker
    SEE ALSO: plcolor,plmk,plp;
*/
{  // see bytscl for palette;
  local c1,ll,i;
  if(is_void(pal)) pal="idl-03.gp";
  if(is_void(marker)) marker=2;  
  if(is_void(width)) width=5;
  if(is_void(type)) type=0;
  /*
    palette,pal;
    c1 =long(100+100*(c-c(min))/(c(max)-c(min)));
    ll=dimsof(c)(2);
  */
  c1=bytscl(val)+char(1);
  palette, red, green, blue, query=1;
  c1=[red(c1),green(c1),blue(c1)];
  info,c1;
  for(i=1;i<=dimsof(c1)(2);i++){ plg,y(i),x(i),type=type,marker=marker,width=width,msize=msize,color=c1(i,);}
}



func plNDnetFlags(net,mask=,equal=,flagtype=,msize=,color=,edges=)
{
  if (is_void(color)) color=__black;
  if (is_void(flagtype)) flagtype=net.ndims_net;
  if (flagtype==0)
    {
      pos=*net.v_coord;
      val=*net.v_flag;
      if (!is_void(mask)) val=val&mask;
      if (!is_void(equal))
        {
          w=where(val==equal);stat,w;
          plg2,pos(2,w),pos(1,w),color=color,msize=msize;
        }
      else
        {            
          myplcol,(remap(val)),pos(2,),pos(1,),width=1;
        }
    }
  if (flagtype==net.ndims_net)
    {
      pos=*net.v_coord;
      flist=*(net.f_vertexIndex(net.ndims_net));
      x=pos(1,);y=pos(2,);
      p1=flist(1,);
      p2=flist(2,);
      p3=flist(3,);

      if (net.periodicity)
      {
        w1 = abs([x(p1),x(p2),x(p3)] - [x(p1),x(p2),x(p3)](,avg))(,max) < (*net.delta)(1)/3;
        w2 = abs([y(p1),y(p2),y(p3)] - [y(p1),y(p2),y(p3)](,avg))(,max) < (*net.delta)(2)/3;
        w1=where(w1&w2);w2=[];
      }
      else w1=indgen(net.nfaces(net.ndims_net));
        
      
      val=*net.f_flag(net.ndims_net);
      if (!is_void(mask)) val=val&mask;
      if (!is_void(equal))
        {
          w=where(val==equal);
          for (j=1;j<=numberof(w);j++)
            {
              i=w1(w(j));
              plf,color,[[y(p1(i)),y(p2(i))],[y(p3(i)),y(p3(i))]],[[x(p1(i)),x(p2(i))],[x(p3(i)),x(p3(i))]],edges=edges;
            }
        }
      else
        {
          val=remap(val);
          c=bytscl(val)+char(1);
          palette, red, green, blue, query=1;
          c=[red(c),green(c),blue(c)];
          c=array(c,[2,2,2]);
          for (j=1;j<=net.nfaces(net.ndims_net);j++)
            {
              i=w1(j);
              plf,c(i,,,),[[y(p1(i)),y(p2(i))],[y(p3(i)),y(p3(i))]],[[x(p1(i)),x(p2(i))],[x(p3(i)),x(p3(i))]],edges=edges;
            }
        }
    }
}

func plNDnetData(net,name,edges=,randomize=,logscale=,valmin=,valmax=,fromvert=)
{
    c_c;
    if (is_void(fromvert)) fromvert=0;
    pi=NDNetDataIndex(net,net.ndims,name);
    if ((pi==-1)||(fromvert)) {
      pi=NDNetDataIndex(net,0,name);
      fromvert=1;
      if (pi==-1) return;
    }
    flist=*(net.f_vertexIndex(net.ndims));
    val = *((*net.data)(pi).data);
    /*
    id=sort(val);
    v=val(id)(dif);
    v=int(v);v(where(v))=1;
    v=v(cum);val=v;
    */
    if (fromvert) val=val(flist)(avg,);
    
    if (randomize)
        val = double(sort(random(numberof(val)))) (1+int(val));
    if (logscale)
        val=log(val);
    stat,val;
    if ((!is_void(valmin))||(!is_void(valmax)))
    {
        if (is_void(valmin)) valmin=min(val);
        if (is_void(valmax)) valmax=max(val);
    }
    
    c=bytscl(val)+char(1);
    palette, red, green, blue, query=1;
    c=[red(c),green(c),blue(c)];
    c=array(c,[2,2,2]);
    
    pos=*net.v_coord;
    x=pos(1,);y=pos(2,);
    p1=flist(1,);
    p2=flist(2,);
    p3=flist(3,);
    valmin;valmax;
    if (net.periodicity)
      {
        w1 = abs([x(p1),x(p2),x(p3)] - [x(p1),x(p2),x(p3)](,avg))(,max) < (*net.delta)(1)/3;info,w1;
        w2 = abs([y(p1),y(p2),y(p3)] - [y(p1),y(p2),y(p3)](,avg))(,max) < (*net.delta)(2)/3;
        w1=where(w1&w2);w2=[];
        if (is_void(valmin)&&is_void(valmax)) {
          for (i=1;i<=numberof(w1);i++)
            plf,c(w1(i),,,),[[y(p1(w1(i))),y(p2(w1(i)))],[y(p3(w1(i))),y(p3(w1(i)))]],[[x(p1(w1(i))),x(p2(w1(i)))],[x(p3(w1(i))),x(p3(w1(i)))]],edges=edges;
        } else {
          for (i=1;i<=numberof(w1);i++)
            if ((val(w1(i))<=valmax)&&(val(w1(i))>=valmin))
              plf,c(w1(i),,,),[[y(p1(w1(i))),y(p2(w1(i)))],[y(p3(w1(i))),y(p3(w1(i)))]],[[x(p1(w1(i))),x(p2(w1(i)))],[x(p3(w1(i))),x(p3(w1(i)))]],edges=edges;
        }
        
      } else
      {
        if (is_void(valmin)&&is_void(valmax)) {
          for (i=1;i<=net.nfaces(net.ndims);i++)
            plf,c(i,,,),[[y(p1(i)),y(p2(i))],[y(p3(i)),y(p3(i))]],[[x(p1(i)),x(p2(i))],[x(p3(i)),x(p3(i))]],edges=edges;
        } else {
          for (i=1;i<=net.nfaces(net.ndims);i++)
            if ((val(i)<=valmax)&&(val(i)>=valmin))
              plf,c(i,,,),[[y(p1(i)),y(p2(i))],[y(p3(i)),y(p3(i))]],[[x(p1(i)),x(p2(i))],[x(p3(i)),x(p3(i))]],edges=edges;
        }
      }
}

func plNDnetData1(net,name,edges=,randomize=,logscale=)
{
    c_c;
    pi=NDNetDataIndex(net,0,name);
    if (pi==-1) return;    
    val = *((*net.data)(pi).data);
    
    if (randomize)
        val = double(sort(random(numberof(val)))) (1+int(val));
    if (logscale)
        val=log10(val);
    
    c=bytscl(val)+char(1);
    palette, red, green, blue, query=1;
    c=[red(c),green(c),blue(c)];
    c=array(c,[2,2,2]);
    flist=*(net.f_vertexIndex(net.ndims));
    pos=*net.v_coord;
    x=pos(1,);y=pos(2,);
    p1=flist(1,);
    p2=flist(2,);
    p3=flist(3,);
    
    for (i=1;i<=net.nfaces(net.ndims);i++)
        plf,char(c([p1(i),p2(i),p3(i)],,,)(avg,,,)),[[y(p1(i)),y(p2(i))],[y(p3(i)),y(p3(i))]],[[x(p1(i)),x(p2(i))],[x(p3(i)),x(p3(i))]],edges=edges;
}

func periodize(&x0,&x1,&y0,&y1,delta,frac=, remove=)
{
  if (is_void(frac)) frac=0.5;
  
  if (!is_void(remove))
    {
      w=where((abs(x0-x1)<=delta(1)*frac)&(abs(y0-y1)<=delta(2)*frac));
      x0=x0(w);x1=x1(w);y0=y0(w);y1=y1(w);
    }
  else
    {
  
      w=where(abs(x0-x1)>delta(1)*frac);
      if (numberof(w)) {
        wp=where(x0(w)>x1(w));
        wn=where(x0(w)<x1(w));  
        if (numberof(wp)) x0(w(wp)) -= delta(1);
        if (numberof(wn)) x1(w(wn)) -= delta(1);
      }
  
  
      w=where(abs(y0-y1)>delta(2)*frac);
    
      if (numberof(w)) {
        wp=where(y0(w)>y1(w));
        wn=where(y0(w)<y1(w));
        if (numberof(wp)) y0(w(wp)) -= delta(2);
        if (numberof(wn)) y1(w(wn)) -= delta(2);
      }
    }
  
}

func plNDnet(net,withnodes=,posonly=, color=, size=, width=, noedges=, noparts=,fromtype=,msize=)
{
    if (is_void(withnodes)) withnodes=0;

    delta = *(net.delta);
    periodic = net.periodicity;
    pos=*net.v_coord;    
    if (is_void(noparts)) plg2,pos(2,),pos(1,),color=color,msize=msize;
    if (posonly) return;

    if (is_void(fromtype))
        fromtype=1;
    
    //if (net.ndata) vpatch = *(*net.data)(4).data;

    if (!net.haveVertexFromFace(fromtype))
      {
        error,"I don't have the right kind of faces ...";
      }

    if (is_void(noedges))
    {
        face=*net.f_vertexIndex(fromtype);
        if (fromtype==1)
        {
          if (periodic)
            {
              x0=pos(1,face(1,));y0=pos(2,face(1,));
              x1=pos(1,face(2,));y1=pos(2,face(2,));
              periodize,x0,x1,y0,y1,delta;
              
              pldj,x0,y0,x1,y1,width=width,color=color;
            }
          else
            {
              pldj,pos(1,face(1,)),pos(2,face(1,)),pos(1,face(2,)),pos(2,face(2,)),width=width,color=color;
            }
        }
        else
        {
          if (periodic) {
            x0=pos(1,face(1,));y0=pos(2,face(1,));
            x1=pos(1,face(2,));y1=pos(2,face(2,));
            periodize,x0,x1,y0,y1,delta;
            pldj,x0,y0,x1,y1,width=width,color=color;
            x0=pos(1,face(1,));y0=pos(2,face(1,));
            x1=pos(1,face(3,));y1=pos(2,face(3,));
            periodize,x0,x1,y0,y1,delta;
            pldj,x0,y0,x1,y1,width=width,color=color;
            x0=pos(1,face(3,));y0=pos(2,face(3,));
            x1=pos(1,face(2,));y1=pos(2,face(2,));
            periodize,x0,x1,y0,y1,delta;
            pldj,x0,y0,x1,y1,width=width,color=color;
          }
          else {
            pldj,pos(1,face(1,)),pos(2,face(1,)),pos(1,face(2,)),pos(2,face(2,)),color=color,width=width;
            pldj,pos(1,face(1,)),pos(2,face(1,)),pos(1,face(3,)),pos(2,face(3,)),color=color,width=width;
            pldj,pos(1,face(3,)),pos(2,face(3,)),pos(1,face(2,)),pos(2,face(2,)),color=color,width=width;
          }
        }
    }

        
    if (net.f_areSimplex)
    {
        if (net.haveVertexFromFace(1))
        {
            face=*net.f_vertexIndex(1);

            if (net.haveFFlags(1))
            {
                //"hello";
                /*
                fl=*(net.f_flag(1));
                w=where(fl&(1<<2));numberof(w);
                if (numberof(w))
                    pldj,pos(1,face(1,w)),pos(2,face(1,w)),pos(1,face(2,w)),pos(2,face(2,w)),width=4,color=__blue;
                */

              
                
                fl=*(net.f_flag(1));
                w=where(fl&(1<<7));numberof(w);
                if (numberof(w))
                    pldj,pos(1,face(1,w)),pos(2,face(1,w)),pos(1,face(2,w)),pos(2,face(2,w)),width=5,color=__blue;
                
                fl=*(net.f_flag(1));
                w=where(fl&(1<<0));numberof(w);
                
                if (numberof(w))
                    pldj,pos(1,face(1,w)),pos(2,face(1,w)),pos(1,face(2,w)),pos(2,face(2,w)),width=5,color=__red;
                
                
                
            }
            //else pldj,pos(1,face(1,)),pos(2,face(1,)),pos(1,face(2,)),pos(2,face(2,)),width=width,color=color;
                 
        }

    }
    else
    {
        error,"plNDnet not implemented for non simplicial networks so far ..."
    }
    
    if (withnodes)
    {
        
        fl=*(net.v_flag);
        
        w=where(fl&(1<<3));numberof(w);
        if (numberof(w)) PL,pos(2,w),pos(1,w),incolor=__red,msize=1;

         w=where(fl&(1<<4));numberof(w);
        if (numberof(w)) PL,pos(2,w),pos(1,w),incolor=__green,msize=1;
        
        w=where(fl&(1<<2));numberof(w);
        if (numberof(w)) PL,pos(2,w),pos(1,w),incolor=__blue,msize=1;
        
    }
    
}

func DelaunayConnStat(net,&ybin,&xbin,dataname=)
{
    local result;

    if (!is_void(dataname))
    {
        val = *(*net.data)(NDNetDataIndex(net,0,"Density")).data;
    }
    
    result=[];
    for (i=1;i<=numberof(net);i++)
    {
        id=*net(i).v_numFaceIndexCum(1);
        flag=((1<<0)|(1<<1));
        w=where(!((*net(i).v_flag)&(flag)));
        n = id(dif)(w);
        x=indgen(max(n)+1)-0.5;
        if (!is_void(dataname))
        {
            v=val(w);
            y=spanl(min(v),max(v),10);
            h=double(histo2d(([double(n),v]),x,y));
            for (j=1;j<=numberof(y)-1;j++)
            {
                h(,j)=double(h(,j))/(max(h(,j)));
            }
            tmp = array(0.,[])
            tmp=h;
            x=x(zcen);
            y=y(zcen);
            xbin=x(,-::numberof(y)-1);
            ybin=y(-::numberof(x)-1,);
        }
        else
        {
            h=double(histo1d(n,x));
            h=double(h)/h(sum);
            tmp=transpose([x(zcen),h]);
        }
        grow,result,&tmp;
        //error,"";
    }
    if (numberof(net)==1) return *result(1);
    
    return result;
}

func plNDrecPatches(name,num=,noplot=,peak=,minima=)
{
    if (is_void(num)) num=2;
    fname=name+".delaunay.NDnet";
    net=[];
    for (i=1;i<=num;i++)
    {
        grow,net,NDnetwork_read(fname+".skeleton");
    net=NDnetwork_read("grad_net.NDnet");    //WS,num+1;
        if (!noplot)
        {
            WS,i,dpi=150;
            plNDnetPatches(net(0),edges=1,peak=peak);
            
            h=DelaunayConnStat(net(0),y,x,dataname="Density");
            
            WS,num+i;
            plcont,h,y,x;
            logxy,0,1;
        }
        if (!minima)
            fname+=".max.delaunay.NDnet";
        else
            fname+=".min.delaunay.NDnet";
    }
    return net;
}

func getCoordsFromNDfield(grid,index,type)
{
    /*
     N=[1024,1024];x0=[0.,0.];delta=[1.,1.];
     grid = [(indgen(N(1))-1)(,-::N(2)-1),(indgen(N(2))-1)(-::N(1)-1,)];
     //grid(,,1) = x0(1)+grid(,,1)/N(1)*delta(1);
     //grid(,,2) = x0(2)+grid(,,2)/N(2)*delta(2);
     grid=grid(*,);
    */
    if (numberof(index)>1)
    {
        result = array(0.,[2,numberof(index),2]);
        w=where(type==0);
        if (numberof(w)) result(w,) = grid(1+index(w),);
        w=where(type==1);
        if (numberof(w))
        {
            result(w,) = grid(1+index(w)/2,);
            j=index(w)%2;
            w2=where(j==0);
            result(w(w2),2)+=0.5;
            w2=where(j==1);
            result(w(w2),1)+=0.5;
        }
        w=where(type==2);
        if (numberof(w))
        {
            result(w,) = grid(1+index(w),)+0.5;
        }
        return result;
    }
    
    
    if (type==0) return double(grid(1+index,));
    if (type==1)
    {
        result = double(grid(1+index/2,));
        if (index%2) result(2)+=0.5;
        else result(1)+=0.5;
        return result;
    }
    if (type==2)
    {
        return double(grid(1+index,))+0.5;
    }
}

func plotGradNDfield(grid,grad)
{
    wok = where(grad(3,)<(1<<30));
    g=grad(,wok);
    w=where((g(1,)==0)&(g(3,)>0));
    wm=where((g(1,)==0)&(g(3,)<0));
    
    p1=getCoordsFromNDfield(grid,g(2,w),array(char(0),numberof(w)));
    p2=getCoordsFromNDfield(grid,g(3,w),array(char(1),numberof(w)));
    tmp = abs(p2-p1);
    w=where((tmp(,1) < 2)&(tmp(,2) < 2))
    pldj,p1(w,1),p1(w,2),p2(w,1),p2(w,2);

    w=where((g(1,)==1)&(g(3,)>0));
    wm=where((g(1,)==1)&(g(3,)<0));
    
    p1=getCoordsFromNDfield(grid,g(2,w),array(char(1),numberof(w)));
    p2=getCoordsFromNDfield(grid,g(3,w),array(char(2),numberof(w)));
    tmp = abs(p2-p1);
    w=where((tmp(,1) < 2)&(tmp(,2) < 2))
    pldj,p1(w,1),p1(w,2),p2(w,1),p2(w,2),color=__green;
}

func getCoordsFromNet(net,index,type,utype=)
{  
    vc=*net.v_coord;
    
    if (numberof(type)>1)
    {
        result = array(float,[2,net.ndims,numberof(index)]);
        if (is_void(utype))
          tlist=unique(type);
        else tlist=utype;
        for (i=1;i<=numberof(tlist);i++)
        {
            cur_type=tlist(i);
            wc = where(type == cur_type);
            if (cur_type==0) result(,wc) = vc(,index(wc));
            else {
                vi = *net.f_vertexIndex(cur_type);
                result(,wc) = vc(,vi(,index(wc)))(,avg,);
            }
        }
        return result;
    }
    
    if (type(1)==0) return vc(,index(1));
    
    vi = *net.f_vertexIndex(type(1));
    return vc(,vi(,index(1)))(,avg,);
}

func getValFromNet(net,index,type, fieldname, fieldval=)
{
   
    if (is_void(fieldval)) {
        pi=NDNetDataIndex(net,0,fieldname);
        if (pi==-1) error, "No field "+fieldname+" available.";
        fieldval = *((*net.data)(pi).data);
    }
    
    if (numberof(type)>1)
    {
        result = array(double,numberof(index));
        tlist=unique(type);
        for (i=1;i<=numberof(tlist);i++)
        {
            cur_type=tlist(i);
            wc = where(type == cur_type);
            if (cur_type==0) result(wc) = fieldval(index(wc));
            else {
                vi = *net.f_vertexIndex(cur_type);
                result(wc) = fieldval(vi(,index(wc)))(avg,);
            }
        }
        return result;
    }
    type=type(1);
    if (type==0) return fieldval(index);
    
    vi = *net.f_vertexIndex(type);
    return fieldval(vi(,index))(avg,);
}


func plNDnetGrad(net,grad,skl=,skl_only=,plmax=, plcrit=)
{
     /* DOCUMENT
        
        //netd=NDnetwork_read("~/prog/dsex/simu_2D.ND.NDnet.dsex");
        skl=int(read_ascii("skl.dat"));
        net=NDnetwork_read("grad_net.NDnet");
        grad=int(read_ascii("grad.dat"));
        WS,1,dpi=150;
        plNDnetData1,net,"Field value",edges=1,logscale=1;
        plNDnetGrad(net,grad,skl=skl);
     */

    w0=where(grad(1,) == 0);
    w1=where(grad(1,) == 1);
    w2=where(grad(1,) == 2);

   
    
    w0r=w0(where((grad(3,w0)<net.nfaces(1)) & (grad(3,w0)>=0)));
    w1r=w1(where((grad(3,w1)<net.nfaces(2)) & (grad(3,w1)>=0)));
    //w2r=w2(where((grad(3,w2)<net.nfaces(2)) & (grad(3,w2)>=0)));
   
    flist=*(net.f_vertexIndex(1));
    flist2=*(net.f_vertexIndex(2));
   
    id=grad(2,w0r)+1;/*w1=where(id<=net.nvertex);*/
    p0 = (*net.v_coord)(,id);
    id=grad(3,w0r)+1;/*w2=where(id<=net.nfaces(1));*/id=flist(,id);
    p0d = (*net.v_coord)(,id)(,avg,);
    //w1=where(w1&w2);
    
    id=grad(2,w1r)+1;/*wa=where(id<=net.nfaces(1));*/id=flist(,id);    
    p1 = (*net.v_coord)(,id)(,avg,);
    id=grad(3,w1r)+1;/*w2=where(id<=net.nfaces(2));*/id=flist2(,id);    
    p1d = (*net.v_coord)(,id)(,avg,);
    //w2=where(wa&w2);

    if (!is_void(plcrit))
      {
        wcrit0 = w0(where(grad(3,w0)==2147483647));
        wcrit1 = w1(where(grad(3,w1)==2147483647));
        wcrit2 = w2(where(grad(3,w2)==2147483647));
    
        id = grad(2,wcrit0)+1;
        pc0 = (*net.v_coord)(,id);

        id = grad(2,wcrit1)+1;
        id=flist(,id);    
        pc1 = (*net.v_coord)(,id);
        
        id = grad(2,wcrit2)+1;
        id = flist2(,id);    
        pc2 = (*net.v_coord)(,id);

        
        
        
        N=dimsof(pc2)(0);
        z=array(__red,[2,2,2]);
        x=array(double,3);
        y=array(double,3);
        for (i=1;i<=N;i++)
          {
            x=pc2(1,,i);
            y=pc2(2,,i);
            x=[[x(1),x(2)],[x(3),x(3)]];
            y=[[y(1),y(2)],[y(3),y(3)]];
            plf,z,y,x
          }

        x0=pc1(1,1,);y0=pc1(2,1,);x1=pc1(1,2,);y1=pc1(2,2,);
        if (net.periodicity) periodize,x0,x1,y0,y1,*net.delta,frac=0.1,remove=1;
        pldj,x0,y0,x1,y1,color=__green,width=4;

        PL,pc0(2,),pc0(1,),incolor=__blue,color=__black,msize=1;
        
        
      }
    
    
    //p2 = (*net.v_coord)(,);
    if (is_void(skl_only))
    {
      
      //plNDnet,net;
        if (net.periodicity) {
          x0=p0(1,);y0=p0(2,);x1=p0d(1,);y1=p0d(2,);
          periodize,x0,x1,y0,y1,*net.delta,frac=0.1,remove=1;
          pldj,x0,y0,x1,y1,color=__green,width=2;
          
          x0=p1(1,);y0=p1(2,);x1=p1d(1,);y1=p1d(2,);
          periodize,x0,x1,y0,y1,*net.delta,frac=0.1,remove=1;
          pldj,x0,y0,x1,y1,color=__red,width=2;
        } else
          {
            pldj,p0(1,),p0(2,),p0d(1,),p0d(2,),color=__green,width=2;
            pldj,p1(1,),p1(2,),p1d(1,),p1d(2,),color=__red,width=2;
          }
        
    }
    
    if (!is_void(skl))
    {
        segp = (*net.v_coord)(,flist)(,avg,);
        nodep = (*net.v_coord);
        facep = (*net.v_coord)(,flist2)(,avg,);

        w0=where(skl(1,)==0);
        w1=where(skl(1,)==1);
        w2=where(skl(1,)==2);

        if (numberof(w0)) w01=w0(where(skl(3,w0)==1));
        if (numberof(w1)) w10=w1(where(skl(3,w1)==0));
        if (numberof(w1)) w12=w1(where(skl(3,w1)==2));
        if (numberof(w2)) w21=w2(where(skl(3,w2)==1));
        
        if (numberof(w01)) pldj,segp(1,skl(4,w01)+1),segp(2,skl(4,w01)+1),nodep(1,skl(2,w01)+1),nodep(2,skl(2,w01)+1),color=__pink,width=3;
        if (numberof(w10)) pldj,segp(1,skl(2,w10)+1),segp(2,skl(2,w10)+1),nodep(1,skl(4,w10)+1),nodep(2,skl(4,w10)+1),color=__pink,width=3;

        if (numberof(w12)) pldj,segp(1,skl(2,w12)+1),segp(2,skl(2,w12)+1),facep(1,skl(4,w12)+1),facep(2,skl(4,w12)+1),color=__orange,width=3;
        if (numberof(w21)) pldj,facep(1,skl(2,w21)+1),facep(2,skl(2,w21)+1),segp(1,skl(4,w21)+1),segp(2,skl(4,w21)+1),color=__orange,width=3;
        if (!is_void(plmax))PL,facep(2,skl(2,w21)+1),facep(1,skl(2,w21)+1),color=__red,msize=0.1;
    }
}

func readCycles(fname,nmin=)
{
    if (is_void(fname)) fname="cycles.dat";
    if (is_void(nmin)) nmin=0;
     
    file=open(fname,"r");
    N = 0;
    read,file,N;
    result=[];
    dataresult=[];
    for (i=0;i<N;i++)
    {
        nseg=([0,0,0,0,0]);
        read,file,nseg;
        //nseg(1)+=1;
        //nseg(3)+=1;
        tmp=[0,0,0,0];
        //nseg;
        if (nseg(5)>=nmin) {
            q=array(int,[2,4,nseg(5)+1]);
            q(,1) = nseg(1:4);
            for (j=2;j<=nseg(5)+1;j++) {
                read,file,tmp;                
                q(,j)=tmp;
            }
            //q(1,)+=1;
            //q(3,)+=1;
            grow,result,&q;
        }
    }
    
    close,file;
    return result;
}

func cycleVal(net,cycle,fieldname=,fieldindex=)
{
    if (is_void(fieldindex))
    {
        if (is_void(fieldname))
            fieldname="Field value";
    

        pi=NDNetDataIndex(net,0,fieldname);
        if (pi==-1) error, "No field "+fieldname+" available.";
        fieldval = *((*net.data)(pi).data);
    } else fieldval = *((*net.data)(fieldindex).data);

    
    result=[];
    /*
    c=array(int,[2,4,numberof(cycle)]);
    for (i=1;i<=numberof(cycles);i++)
    {
        c(,i) = [(*(cycles(i)))(,1)];
    }
    */
    
    //stat,c(1,);
    d1=getValFromNet(net,cycle(1,)+1,cycle(2,),fieldval=fieldval);
    d2=getValFromNet(net,cycle(3,)+1,cycle(4,),fieldval=fieldval);    

    result = [d1,d2,abs(d2-d1)];
    return result;
    /*
    for (i=1;i<numberof(cycles);i++)
    {
        i;
        
        c=(*(cycles(i)))(,1);
        
        d1=getValFromNet(net,c(1),c(2),fieldval=denval);
        d2=getValFromNet(net,c(3),c(4),fieldval=denval);
        ld1=getValFromNet(net,c(1),c(2),fieldval=ldenval);
        ld2=getValFromNet(net,c(3),c(4),fieldval=ldenval);

        grow,result,[d1,d2,ld1,ld2,abs(d2-d1)]; 
    }
*/
}

func plotCycles(net,cycles,randcol=,color=,nmin=,nmax=,noise=,width=,start=,stop=,profile=)
{
    if (is_void(nmin)) nmin=1;
    if (is_void(nmax)) nmax=100000000;
    if (is_void(width)) width=1;
    if (is_void(color)) col=char(random(numberof(cycles)+1,3)*255);

    if (is_void(start)) start=1;
    if (is_void(stop)) stop=numberof(cycles);
    
    cur_window=window();
    
    /*
    if (!is_void(profile)) {
        pi=NDNetDataIndex(net,net.ndims,name);
        if (pi==-1) {
            pi=NDNetDataIndex(net,0,name);
            if (pi==-1) return;
        }
        density = 
    }
    */
    
    if (!is_void(profile)) {
        if (!is_string(profile)) profile = "Field value";
        fieldindex = NDNetDataIndex(net,0,profile);
        window,cur_window+1;
        logxy,1,1;
        fma;
    }
    
    for (i=start;i<=stop;i++)
    {
        
        
        c=*(cycles(i));
        if (dimsof(c)(3)<=nmin) continue;
        if (dimsof(c)(3)>=nmax) continue;
        
        p1=getCoordsFromNet(net,c(1,)+1,c(2,));
        p2=getCoordsFromNet(net,c(3,)+1,c(4,));
        //info,p1;
        y=[p1(2,2:),p2(2,2:)];
        x=[p1(1,2:),p2(1,2:)];
        if (!is_void(noise))
        {
            //y(,1)+=0.75*(random(dimsof(y)(2)/2)(-::1,)(*)-0.5)*(y(,dif)(*)(rms));
            //x(,1)+=0.75*(random(dimsof(x)(2)/2)(-::1,)(*)-0.5)*(x(,dif)(*)(rms));
            dy=(abs(y(::2,2)-y(2::2,2)))(-::1,)(*);
            dx=(abs(x(::2,2)-x(2::2,2)))(-::1,)(*);
            //dy=(y(::2,dif)(-::1,)(*));
            //dx=(x(::2,dif)(-::1,)(*));
            y(,1)+=noise*0.75*(random(dimsof(y)(2)/2)(-::1,)(*)-0.5)*dy;
            x(,1)+=noise*0.75*(random(dimsof(x)(2)/2)(-::1,)(*)-0.5)*dx;
        }

        if (!is_void(profile)) window,cur_window;
        
        //dimsof(x)(2)
        if (is_void(color)) {
            for (j=1;j<=dimsof(x)(2);j++) plg,y(j,),x(j,),color=col(i,),width=width;
            pldj,p1(1,1),p1(2,1),p2(1,1),p2(2,1),color=col(i,),width=width+1;
        }
        else {
            for (j=1;j<=dimsof(x)(2);j++) plg,y(j,),x(j,),color=color,width=width;
            pldj,p1(1,1),p1(2,1),p2(1,1),p2(2,1),color=color,width=width+1;
        }

        if (!is_void(profile))
        {
            window,cur_window+1;
            val = cycleVal(net,c,fieldindex=fieldindex);
            if (is_void(color)) {
                PL,val(1,3),val(1,1),color=char(0.75*col(i,)),incolor=col(i,),msize=1.5;
                plg2,val(2:,3),val(2:,1),msize=2,color=char(col(i,));
            } else {
                plg2,val(1,3),val(1,1),msize=10,color=color;
                plg2,val(2:,3),val(2:,1),msize=2,color=char(color*0.75);
            }
            
            info,val;
            }
        //if (start-stop < 20)
            i;
        //pldj,p1(1,2:),p1(2,2:),p2(1,2:),p2(2,2:),color=col(i,);
        
    }
    if (!is_void(profile)) window,cur_window;
}

func persistencePlot(ppairs,delta=,xlim=,nbins=,noerase=,linetype=)
{
    /* DOCUMENT
       TO plot persistence pairs :

       ppairs=read_ascii("pairs.dat");
       p1 = getCoordsFromNet(net,int(1+ppairs(1,)),int(ppairs(2,)));
       p2 = getCoordsFromNet(net,int(1+ppairs(3,)),int(ppairs(4,)));
       t1=ppairs(2,);
       t2=ppairs(4,);
       wu=where((t1==net.ndims)|(t2==net.ndims));
       wd=where((t1==0)|(t2==0));
       //WS,dpi=175;
       //plNDnet,net,posonly=1;
       pldj,p1(1,),p1(2,),p2(1,),p2(2,),color=__green;
       pldj,p1(1,wu),p1(2,wu),p2(1,wu),p2(2,wu),color=__red;
       pldj,p1(1,wd),p1(2,wd),p2(1,wd),p2(2,wd),color=__blue;
     */

    if (is_void(delta)) delta = 0;
    if (is_void(nbins)) nbins=100;
/*
    x=spanl(min(ppairs(5,)),max(ppairs(5,)),nbins);
    y=spanl(xlim(1),xlim(2),nbins);
    */
    if (is_void(xlim)) x=spanl(min(ppairs(5,)),max(ppairs(5,)),nbins);
    else x=spanl(xlim(1),xlim(2),nbins);

    if (is_void(xlim)) y=spanl(min(ppairs(6,)),max(ppairs(6,)),nbins);
    else y=spanl(xlim(1),xlim(2),nbins);
    
    
    t1=ppairs(2,);
    t2=ppairs(4,);
    ndims=max([max(t1),max(t2)]);
    w=where((ppairs(5,)>0));
    wu=where((ppairs(5,)>0)&((t1==ndims)|(t2==ndims)));
    wd=where((ppairs(5,)>0)&((t1==0)|(t2==0)));
    wi=where((ppairs(5,)>0)&((t1!=ndims)&(t2!=ndims))&((t1!=0)&(t2!=0)));

    //moy=ppairs(5:6,)(*)(avg);
    moy=1;
    //per = abs(log(ppairs(5:6,)/moy)(dif,))(*);
    
    per = abs((ppairs(5:6,)/moy)(dif,))(*);
    
    wp=where(per >0);
    wpu=where((per>0)&((t1==ndims)|(t2==ndims)));
    wpd=where((per>0)&((t1==0)|(t2==0)));
    wpi=where((per>0)&((t1!=ndims)&(t2!=ndims))&((t1!=0)&(t2!=0)));

    h=histo2d(transpose(ppairs(5:6,)),x,y);
    hu=histo2d(transpose(ppairs(5:6,wu)),x,y);
    hd=histo2d(transpose(ppairs(5:6,wd)),x,y);
    if (numberof(wi)) hi=histo2d(transpose(ppairs(5:6,wi)),x,y);

    //moy = abs(ppairs(5:6,)(min,));
    /*
    moy=ppairs(5:6,)(*)(avg);
    if (numberof(where(moy==0))) moy(where(moy==0))=1;
    per = abs(ppairs(5:6,)(dif,))(*)/moy;
    */

    if (!is_void(xlim)) xp=spanl(xlim(1),xlim(2),nbins);
    else xp=spanl(min(per(wp)),max(per(wp)),nbins);

    
    hp=histo1d(per(wp),xp);
    hpu=histo1d(per(wpu),xp);
    hpd=histo1d(per(wpd),xp);
    if (numberof(wpi)) hpi=histo1d(per(wpi),xp);

    //if (is_void(noerase)) WS,5*delta+0,dpi=150;
    //else window,5*delta+0;
    WS,5*delta,dpi=150;
    plcont,log(1+hd),y(zcen)(-::98,),x(zcen)(,-::98);
    c_c;
    logxy,1,1;

    //if (is_void(noerase)) WS,5*delta+1,dpi=150;
    //else window,5*delta+1;
    WS,5*delta+1,dpi=150;
    plcont,log(1+hu),y(zcen)(-::98,),x(zcen)(,-::98);
    c_c;
    logxy,1,1;

    WS,5*delta+2,dpi=150;
    plcont,log(1+hu + transpose(hd)),y(zcen)(-::98,),x(zcen)(,-::98);
    c_c;
    logxy,1,1;

    if (numberof(wi)) {
        //if (is_void(noerase)) WS,5*delta+2,dpi=150;
        //else window,5*delta+2;
        WS,5*delta+2,dpi=150;
        plcont,log(1+hi),y(zcen)(-::98,),x(zcen)(,-::98);
        c_c;
        logxy,1,1;
    }

    if (is_void(noerase)) WS,5*delta+3,dpi=150;
    else window,3;
    plg,1+hp,xp(zcen),width=2,type=linetype;
    plg,1+hpu,xp(zcen),color=__red,width=2,type=linetype;
    plg,1+hpd,xp(zcen),color=__blue,width=2,type=linetype;
    if (numberof(wpi)) plg,1+hpi,xp(zcen),color=__green,width=2;
    logxy,1,1;
    
    if (is_void(noerase)) WS,5*delta+4,dpi=150;
    else window,4;
    plg,1+hp(cum),xp,width=2,type=linetype;
    plg,1+hpu(cum),xp,color=__red,width=2,type=linetype;
    plg,1+hpd(cum),xp,color=__blue,width=2,type=linetype;
    if (numberof(wpi)) plg,1+hpi(cum),xp,color=__green,width=2;
    logxy,1,1;
/*
    WS,5*delta+5,dpi=150;
    w=where(hpd>0);
    plg,(1.+hpu(w))/(1.+hpd(w)),xp(zcen)(w),width=2;


    WS,5*delta+5,dpi=150;
    plg,(double(hp)(dif)),xp(zcen)(zcen),width=2;
    plg,(double(hpu)(dif)),xp(zcen)(zcen),color=__red,width=2;
    plg,(double(hpd)(dif)),xp(zcen)(zcen),color=__blue,width=2;
    if (numberof(wpi)) plg,hpi(dif),xp(zcen)(zcen),color=__green,width=2;
    logxy,1,1;
    hpd(dif);
    xp(dif)(zcen);
*/
}

func EuclidianSort(tab, curindex=, index=, perm=)
{
  
  if (is_void(index)) index = dimsof(tab)(mxx)-1;
  if (is_void(curindex)) curindex = 1;
  
  
  local ndims;
  ndims= dimsof(tab)(1);
    
  t=tab;
  if (ndims!=index) {ndims;index;t = transpose(t,[index,ndims]);}
  init_dims = dimsof(t);
  
  if (ndims!=2) t=t(*,);
  if (perm)
    {
      inv_rec = array(int,dimsof(t));
      ind = indgen(dimsof(t)(2));
      for (i=1;i<dimsof(t)(0);i++)
        {
          h = heapsort(t(,i));
          t(,i) = t(h,i);
          inv_rec(h,i) = ind;
        }
    }
  
  t=t(,sort(t(curindex,)));
  
  if (dimsof(t)(2)==curindex) return t;
                                         
  w=where(t(curindex,)(dif)==0);
  if (numberof(w)) {
    if (numberof(w)==1) wi=[];
    else wi=where(w(dif)>1);
  start=array(int,numberof(wi)+1);
  stop=array(int,numberof(wi)+1);
  if (numberof(wi)) {
    start(2:) = w(wi()+1);
    stop(:-1) = w(wi)+1;
  }
  start(1)=w(1);
  stop(0)=w(0)+1;
  w=[];wi=[];
  nloop = numberof(start);
  newid=curindex+1;
  for (i=1;i<=nloop;i++)
    {
      t(,start(i):stop(i)) = EuclidianSort(t(,start(i):stop(i)),curindex=newid,index=index);
    }
  
  }

  if (perm)
    {
      for (i=1;i<dimsof(t)(0);i++)
        {
          t(,i) = t(inv_rec(,i),i);
        }
    }
  
  if (ndims!=2) t=reform(t,init_dims);
                  
  if (ndims!=index) t = transpose(t,[index,ndims]);

  return t;
}

func SetDifference(tab1, tab2, &set1,&set2,index=,perm=)
{
  
  if (is_void(index)) index = dimsof(tab1)(mxx)-1;

  w=where(dimsof(tab1)!=dimsof(tab2));
  if (numberof(w) != 0)
    {
      w-=1;
      if ((numberof(w)>1)||(w(1)!=index))
        error,"tab1 to tab2 ha zenzen chigau !!!";
    }
  
  local ndims;
  ndims= dimsof(tab1)(1);
  
  t1=EuclidianSort(tab1,index=index,perm=perm);
  t2=EuclidianSort(tab2,index=index,perm=perm);

  if (ndims!=index) {
    t1 = transpose(t1,[index,ndims]);
    t2 = transpose(t2,[index,ndims]);
  }
  
  i1=1;i2=1;
  n1=dimsof(t1)(0);
  n2=dimsof(t2)(0);

  set1=[];
  set2=[];
  
  do {
    
    ta=t1(*,i1);
    tb=t2(*,i2);
    if (perm)
      {
        heapsort,ta;
        heapsort,tb;
      }
    
    w=where(ta!=tb);
    
    if (numberof(w))
      {
        if (ta(w(1))<tb(w(1))) {
          grow,set1,t1(*,i1);
          i1++;
        }
        else {
          grow,set2,t2(*,i2);
          i2++;
        }
      }
    else {
      i1++;i2++;
    }
        
  } while ((n1>=i1)&&(n2>=i2));
  
  
  if (n1>i1) {
    grow,set1,t1(,i1:)(*);
  }
  if (n2>i2) {
    grow,set2,t2(,i2:)(*);
  }
  
  if (numberof(set1)) {
    tmp=t1;
    newdims = dimsof(tmp);
    newdims(0)=numberof(set1)/numberof(tmp(*,1));
    set1=reform(set1,newdims);
  }

  if (numberof(set2)) {
    tmp=t2;
    newdims = dimsof(tmp);
    newdims(0)=numberof(set2)/numberof(tmp(*,1));
    set2=reform(set2,newdims);
  }
  
  if (ndims!=index) {
    tab1 = transpose(tab1,[index,ndims]);
    tab2 = transpose(tab2,[index,ndims]);
    if (numberof(set1)) set1 = transpose(set1,[index,ndims]);
    if (numberof(set2)) set2 = transpose(set2,[index,ndims]);
  }
  
}


/*
Some codes ...


plNDnetData(net,"Log Field value")

sklt=skl_fof;
seg=*sklt.segpos;
npos=*sklt.nodepos;
node=(*sklt.node);
npos=(*sklt.nodepos);
nd = 0;

for (nd=1;nd<sklt.nnodes;nd++) {
//nd++;
N=0;
node(nd).nnext;
for (i=1;i<=node(nd).nnext;i++) {
if ((i%10)==0) i;
N=N+1;
col=char(100+random(3)*155);
sl=seg(,,*NDfilaments(sklt,nd,0)(N));pldj,sl(1,1,),sl(2,1,),sl(1,2,),sl(2,2,),color=col,width=2;
nnode = (*node(nd).next)(N);
PL,npos(2,nnode),npos(1,nnode),incolor=col;
}
}

 */

/*
  
cd,"/home/thierry/prog/morse_tmp";
#include "../morse/yorick/network.i"
#include "../morse/yorick/skel.i"

RELOAD=1;

nmin=["20","50","100"];
netname="net_id."+nmin+".NDnet";
cyclename="cycles."+nmin+".dat";
ppairsname="remaining_pairs."+nmin+".dat";
sklname="skl_up."+nmin+".NDskl";
if (RELOAD) {
ppairs_p=[];
cycle_p=[];
net_p=[];
skl_p=[];
}
for (i=1;i<=numberof(nmin);i++)
{
if (RELOAD) {
    grow,net_p,&NDnetwork_read(netname(i));
    grow,cycle_p,&readCycles(cyclename(i));
    grow,ppairs_p,&read_ascii(ppairsname(i));
    grow,skl_p,&NDskel_read(sklname(i));
    }

net=*net_p(i);cycle=*cycle_p(i);ppairs=*ppairs_p(i);skl=*skl_p(i);

WS,i,dpi=150;
plNDnetData,net,"Log Field value";
plotFOF(net,randcol=1,peak=1);
plNDskel,skl;

p1 = getCoordsFromNet(net,int(1+ppairs(1,)),int(ppairs(2,)));
       p2 = getCoordsFromNet(net,int(1+ppairs(3,)),int(ppairs(4,)));
       t1=ppairs(2,);
       t2=ppairs(4,);
       wu=where((t1==net.ndims)|(t2==net.ndims));
       wd=where((t1==0)|(t2==0));
       //WS,dpi=175;
       //plNDnet,net,posonly=1;
       pldj,p1(1,),p1(2,),p2(1,),p2(2,),color=__green;
       pldj,p1(1,wu),p1(2,wu),p2(1,wu),p2(2,wu),color=__red;

}
*/

/*
fma;
plotFOF(net,randcol=1,peak=1);
plNDnet,net,posonly=1;
//i=106;
for (i=116;i<126;i++)
{
j=i;
//i+=10;j=i+10;
col=char(127+random(3)*128);
plotCycles(net,cycles,start=i,stop=j,width=5,nmin=10,color=char(col*0.75));
plotCycles(net,ccycles,start=i,stop=j,width=3,nmin=10,color=col,noise=0.25);
}
plNDskel,sklc,withnodes=1,color=__red;
pldj,p1(1,),p1(2,),p2(1,),p2(2,),color=__green;

plNDskel,skl,withnodes=1,color=__red;

w=where(nodes.type==1)
i=9;
p=(*skl.segpos)(,,*NDfilaments(skl,w(i),0)(1))(,avg,);
q=(*skl.segpos)(,,*NDfilaments(skl,w(i),0)(2))(,avg,);
plg,p(2,),p(1,),width=5,color=__blue;
plg,q(2,),q(1,),width=5,color=__green;

i=1109;
p=(*skl.segpos)(,,*NDfilaments(skl,(i),0)(1))(,avg,);
q=(*skl.segpos)(,,*NDfilaments(skl,(i),0)(2))(,avg,);
plg,p(2,),p(1,),width=5,color=__blue;
plg,q(2,),q(1,),width=5,color=__green;

pldj,p1(1,),p1(2,),p2(1,),p2(2,),color=__orange;
pldj,p1(1,wu),p1(2,wu),p2(1,wu),p2(2,wu),color=__purple;

 */


/* LOADS stuff

#include "../morse/yorick/network.i"
#include "../morse/yorick/skel.i"

skld=NDskel_read("skl_down.NDskl");
sklc=NDskel_read("skl_up.NDskl");

cycles=readCycles();
ccycles=readCycles("canonic_cycles.dat");
ppairs=read_ascii("remaining_pairs.dat");
ppairs_all=read_ascii("pairs.dat");

net=NDnetwork_read("net_id.NDnet");
p1 = getCoordsFromNet(net,int(1+ppairs(1,)),int(ppairs(2,)));
p2 = getCoordsFromNet(net,int(1+ppairs(3,)),int(ppairs(4,)));
t1=ppairs(2,);
t2=ppairs(4,);
wu=where((t1==net.ndims)|(t2==net.ndims));
wd=where((t1==0)|(t2==0));
   
 */


/*
//WS,5,dpi=150;
fma;
plNDnetData,net,"Log Field value";
plotFOF(net,randcol=1,peak=1);
plNDnet,net,posonly=1;

//cplim,4,5;

i0=113;delta=10;
for (i=i0;i<i0+delta;i++) {
j=i;
//i+=10;j=i+10;
col=char(127+random(3)*128);
plotCycles(net,cycles,start=i,stop=j,width=5,nmin=10,color=char(col*0.75));
plotCycles(net,ccycles,start=i,stop=j,width=3,nmin=10,color=col,noise=0.25);
}
plNDskel,skld,withnodes=1,color=__green;
plNDskel,sklc,withnodes=1,color=__red;
pldj,p1(1,),p1(2,),p2(1,),p2(2,),color=__white;
pldj,p1(1,wu),p1(2,wu),p2(1,wu),p2(2,wu),color=__orange;
 */

/*
  cycle A detruire: 167

  fma
plNDnetData,net,"Log Field value";
plotFOF(net,randcol=1,peak=1);
plotCycles(net,ccycles,start=167,stop=167,width=3,nmin=800,nmax=10000,profile=1,color=__pink);
plNDskel,sklc,withnodes=1,color=__red;
pldj,p1(1,),p1(2,),p2(1,),p2(2,),color=__orange;
pldj,p1(1,wu),p1(2,wu),p2(1,wu),p2(2,wu),color=__purple;
plotCycles(net,cycles,start=167,stop=167,width=1,nmin=1,nmax=10000,color=__green);

*/



/*
//skel and grad

#include "~/prog/morse/yorick/network.i"
#include "~/prog/morse/yorick/skel.i"

net=NDnetwork_read("simu_2D.ND.NDnet");
grad=int(read_ascii("grad.dat"));
sklu=NDskel_read("simu_2D.ND.NDnet.up.NDskl");
skld=NDskel_read("simu_2D.ND.NDnet.down.NDskl");

WS,0,dpi=150;
plNDnetData,net,"Log Field value"
plNDnetGrad(net,grad);
index=[923952,34135,34208,43892,2699,772359,126861,126594,335466,335451]+1
type=[1,0,0,2,2,1,0,0,2,2]
for (i=1;i<=numberof(type);i++) PL,coords(2,i),coords(1,i),incolor=[__blue,__green,__red](1+type(i),),msize=1.5;
plNDskel,sklu,color=__yellow
plNDskel,skld,color=__pink

 */
