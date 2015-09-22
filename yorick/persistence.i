struct persistencePairsList {
  pointer pairs;
  pointer data;
  pointer data_name;
  pointer pos;
}

struct critPointsList {
  int minType;
  int maxType;
  int minPairType;
  int maxPairType;
  pointer density;
  pointer persistence;
  pointer persistenceR;
  
  pointer sign;
  pointer pairType;
  pointer type;
  pointer index;
}

struct persistenceHisto {
  int pairType;
  pointer x;
  pointer y;
  pointer h;
}

func getCutVal(nsig)
{
  r=exec("mse -dispSLevel -nsig "+swrite(format="%f",double(nsig)))(2:);
  return str2double(r);
}

func getNDnetPairsId(pairs)//,&x,&dat)
{
  result = persistencePairsList();
  id=*(pairs.f_vertexIndex)(1);
  result.pos = &((*pairs.v_coord)(,id));

  w=where((*pairs.data).type==0);
  data = (*pairs.data)(w);

  v_data=array(double,[2,numberof(data),pairs.nvertex]);
  for (i=1;i<=numberof(data);i++) v_data(i,) = *(data(i).data);
  result.pairs= &(long(v_data([1,2],id)));
  if (dimsof(v_data)(2)==4) result.data = &(v_data([3,4,4],id));
  else result.data = &(v_data([3,4,5],id));
  result.data_name= &((data.name)(3:));
  return result;
}

func getPairType(pairs)
{
  return pairs([1,3],..)(min,);
}

func getCritDist(pairs, ravg=)//,data,&id)
{
  if (is_void(ravg)) ravg=1.;
  data = *pairs.data;
  p=(*pairs.pairs);
  list=critPointsList();

  list.density = &(data(1,*)/ravg);
  list.persistence = &(data(2,*));
  list.persistenceR = &(data(3,*));

  list.type = &(char(p(1,*)));  
  list.index = &(long(p(2,*)));
  sgn=(p(1,..)*0+1);
  w=where(p(1,1,)<p(1,2,));
  if (numberof(w)) sgn(1,w)=-1;
  w=where(p(1,2,)<p(1,1,));
  if (numberof(w)) sgn(2,w)=-1;
  w=where(p(1,2,)==p(1,1,));
  if (numberof(w)) sgn(,w)=0;
  list.sign=&(sgn(*));

  list.pairType = &(*list.type);
  w=where(*list.sign == 1);
  if (numberof(w)) (*list.pairType)(w)-=1;
  
  list.minType = (*list.type)(*)(min);
  list.maxType = (*list.type)(*)(max);
  list.minPairType = (*list.pairType)(*)(min);
  list.maxPairType = (*list.pairType)(*)(max);

  return list;
}

func cutCritDist(crit,cut)
{
  w=[];
  for (i=crit.minPairType;i<=crit.maxPairType;i++)
    {
      t=where((*crit.persistenceR>cut(i+1))&(*crit.pairType==i));
      if (numberof(t)) grow,w,t;
    }
  
  newcrit=crit;
  if (numberof(w))
    {
      newcrit.density = &(*newcrit.density)(w);
      newcrit.persistence = &(*newcrit.persistence)(w);
      newcrit.persistenceR = &(*newcrit.persistenceR)(w);
      newcrit.sign = &(*newcrit.sign)(w);
      newcrit.type = &(*newcrit.type)(w);
      newcrit.pairType = &(*newcrit.pairType)(w);
      newcrit.index = &(*newcrit.index)(w);
    }
  else {
    newcrit.density = &[];
    newcrit.persistence = &[];
    newcrit.persistenceR = &[];
    newcrit.sign = &[];
    newcrit.type = &[];
    newcrit.pairType = &[];
    newcrit.index = &[];
    newcrit.minType = -1;newcrit.maxType = -1;
    newcrit.minPairType = -1;newcrit.maxPairType = -1;
    
    return newcrit;
  }
  
  newcrit.minType = (*newcrit.type)(*)(min);
  newcrit.maxType = (*newcrit.type)(*)(max);
  newcrit.minPairType = (*newcrit.pairType)(*)(min);
  newcrit.maxPairType = (*newcrit.pairType)(*)(max);
  
  return newcrit;
}

func getBetti(crit)
{
  
  color=[__blue,__green,__yellow,__red];
  B=array(pointer,crit.maxPairType+1);
  for (i=crit.minPairType;i<=crit.maxPairType;i++)
    {
      B(i+1)=&[];
      w=where(*crit.pairType == i);
      //numberof(w);i;
      if (numberof(w)) {
        bv=double((*crit.sign)(w));
        bd=double(*crit.density)(w);
        s=sort(bd);
        bd=bd(s);
        bv=-bv(s);
        bc = bv(cum)(2:);
        B(i+1)=&[bd,bc,bv];
      }      
    }
  return B;
}

func euler(betti,draw=, color=, width=, type=, ravg=, rescale=, nbins=)
{
  if (is_void(ravg)) ravg=1;
  d=[];
  val=[];
  
  for (i=1;i<=numberof(betti);i++)
    {
      b=*betti(i);
      if (numberof(b))
        {
          grow,d,b(,1);
          if (i%2)
            grow,val,b(,3);
          else
            grow,val,-b(,3);
        }
    }
 
  if (!numberof(d)) return [];
  s = sort(d);
  //total = val(sum);
  //val=total-val(s)(cum)(2:);
  val = val(s)(cum)(2:);
  d=d(s);
  if (!is_void(rescale))
    {
      x=1+(d-ravg)/ravg;
      x0=rescale(1);
      delta=rescale(2);
      x1=rescale(3);
      x=exp(x1+(log(x)-x0)*delta);
    }
  else x=1+(d-ravg)/ravg;
  if (!is_void(draw)) {
    if (is_void(nbins)) delta=1;
    else delta = int(numberof(val)/nbins);
    if (delta<1) delta=1;
    plg,val(::delta),x(::delta),width=width,color=color,type=type;
  }
  return transpose([d,val]);
}

func getppairPDF(crit,&grid,&lims,forcelims=,samegrid=, for_pli=, nbins=)
{
  if (is_void(nbins)) nbins=150;
  
  w=where((*crit.sign) == -1);
  local index,d,pr;
  local xdmin,xdmax,ydmin,ydmax;
  
  for (i=0;i<=crit.maxPairType;i++)
    {    
      w0=where((*crit.pairType)(w) == i);
      if (numberof(w0)) {
        grow, index,&w(w0);
        grow,d,&((*crit.density)(*index(0)));stat,*d(0);
        grow,pr,&((*crit.persistenceR)(*index(0)));
      }
      else {
        grow,index,&[];
        grow,d,&[];
        grow,pr,&[];
        grow,xdmin,0;
        grow,xdmax,0;
        grow,ydmin,0;
        grow,ydmax,0;
      }
    }
  result=[];
  lims=[];
  h=[];
  grid=[];
  for (i=0;i<=crit.maxPairType;i++)
  {
    x=*d(i+1);
    if (!numberof(x)) continue;
    y=x*(*pr(i+1));

    if (numberof(forcelims))
      grow,lims,[forcelims(,i+1)];
    else
      grow,lims,[[x(min),x(max),y(min),y(max)]];
    lims(,0);
    
    xg=spanl(lims(1,0),lims(2,0),nbins+1);
    yg=spanl(lims(3,0),lims(4,0),nbins+1);

    grow,h,&histo2d([x,y],xg,yg);

    if (!is_void(for_pli)) grow,grid,&[xg(zcen),yg(zcen)];
    else grow,grid,&([xg(zcen)(,-),yg(zcen)(-,)]); 
  }

  return h;
}

func NDnetPPHisto(crit,type=,color=,delta=, useRatio=)
{
  if (structof(crit)==structof(NDnetwork())) {
    crit = getCritDist(getNDnetPairsId(crit));
  }
  res=[];
  if (is_void(type)) type=indgen(crit.maxPairType+1)-1;
  for (i=1;i<=numberof(type);i++)
    {
      t=type(i);
      w=where((*crit.pairType) == t);
      if (numberof(w)==0) write,format="WARNING: No pair of type %d. SKIPPED\n",t;
      wd=w(where((*crit.type)(w)==t));
      if (is_void(useRatio)) wp=wd(where((*crit.persistence)(wd)>0));
      else wp=wd(where((*crit.persistenceR)(wd)>1));
      x=(*crit.density)(wp);
      if (is_void(useRatio)) {y=(*crit.persistence)(wp);if (t==0) x+=y;}
      else y=x*(*crit.persistenceR)(wp);
      stat,y;
      if (!is_void(delta)) y-=x;
      if (is_void(useRatio))
        {
          xx=span(x(min),x(max),101);
          yy=spanl(y(min),y(max),101);
        }
      else
        {
          xx=spanl(x(min),x(max),101);
          yy=spanl(y(min),y(max),101);
        }
      h=histo2d(([x,y]),xx,yy);
      tmp=persistenceHisto();
      tmp.x=&xx(zcen);
      tmp.y=&yy(zcen);
      tmp.h=&double(h);
      tmp.pairType=t;
      grow,res,&tmp;
      /*
      plg2,y,x,color=color;
      
      return h;
      stat,x;
      stat,y;
      plg2,y,x,color=color;
      */
      /*
      x=spanl(x(min),x(max),200);
      y=spanl(y(min),y(max),200);
      h=histo2D([x,y],x,y);
      grow,res,persistenceHisto(&x,&y,&h);*/
    }

  if (numberof(res)==1) return res(1);
  return res;
}


func NDnetPPairsPlot(ppairs,bbox=,width=,cl=)
{
  wper = where(((*ppairs.data).name=="persistence Ratio")&((*ppairs.data).type==1))(1);  
  pr=(*(*ppairs.data)(wper).data);
  wper=where(pr<1);
  if (numberof(wper)) {
  }
  
  w=where(((*ppairs.data).name=="type")&((*ppairs.data).type==1))(1);
  t=int(*(*ppairs.data)(w).data);
  x=*ppairs.v_coord;
  s=*(ppairs.f_vertexIndex(1));
  p=x(,s);
  //t=int(*crit.pairType);
  if (is_void(cl)) cl=[__red,__green,__blue];
  //cl=[__purple,__purple,__purple];
  x0=p(1,1,);x1=p(1,2,);
  y0=p(2,1,);y1=p(2,2,);
  if (!is_void(bbox)) periodize(x0,x1,y0,y1,bbox);
  /*
  for (j=0;j<=numberof(wper);j++) {
    i=wper(j);
    plg,[y0(i),y1(i)],[x0(i),x1(i)],color=cl(,t(i)),width=width;
    }
  */
  
  for (i=1;i<=dimsof(p)(0);i++) {  
    plg,[y0(i),y1(i)],[x0(i),x1(i)],color=cl(,t(i)),width=width;
  }
  
}

func NDnetPPlot(crit, pldelta=,msize=,marker=, useRatio=, level=, xyscale=)
{
  /* DOCUMENT
     From a NDnet file:
     
     ppairs=NDnetwork_read("FNAME.ppairs.NDnet");
     pairs=getNDnetPairsId(ppairs);
     crit=getCritDist(pairs);
     NDnetPersistencePlot(crit,pldelta=1);

     to plot the pairs:
     ppairs=NDnetwork_read("FNAME.ppairs.NDnet");
     pairs=getNDnetPairsId(ppairs);
     crit=getCritDist(pairs);
     
     x=*ppairs.v_coord;
     s=*(ppairs.f_vertexIndex(1));
     p=x(,s);
     t=int(*crit.pairType)+1;
     cl=[__red,__green,__blue];
     x0=p(1,1,);x1=p(1,2,);
     y0=p(2,1,);y1=p(2,2,);
     periodize(x0,x1,y0,y1,array(50000.,2));
     for (i=1;i<=dimsof(p)(0);i++) {  
     plg,[y0(i),y1(i)],[x0(i),x1(i)],color=cl(,t(i));
     }
   */

  if (structof(crit)==structof(NDnetwork())) {
    crit = getCritDist(getNDnetPairsId(crit));
  }
  
  if (is_void(xyscale)) xyscale=[1.,1.];
  if (numberof(xyscale)==1) xyscale=[xyscale,xyscale];
  if (is_void(msize)) msize=3;
  w=where((*crit.sign) == -1);
  local index,d,p,pr;
  //index=[];
  //d=[];pr=[];
  for (i=0;i<=crit.maxPairType;i++)
    {    
      w0=where((*crit.pairType)(w) == i);
      if (numberof(w0)) {
        grow, index,&w(w0);
        grow,d,&((*crit.density)(*index(0)));
        grow,p,&((*crit.persistence)(*index(0)));        
        grow,pr,&((*crit.persistenceR)(*index(0)));        
      }
      else {
        grow,index,&[];
        grow,d,&[];
        grow,p,&[];
        grow,pr,&[];
      }
    }
  
  color=[__blue,__green,__orange,__red];
  
  crit.maxPairType; 
  for (i=0;i<=crit.maxPairType;i++)
    {
      if (!is_void(level)) {
        w=where(*pr(i+1) >= level(i+1));
        x=*d(i+1);x=x(w);
        if (!numberof(x)) continue;
        y=x+(*p(i+1))(w);
        //y=x*((*pr(i+1))(w));
      }
      else {
        x=*d(i+1);
        if (!numberof(x)) continue;
        y=x+(*p(i+1));
        //y=x*(*pr(i+1));
      }
      y*=xyscale(2);
      x*=xyscale(1);

      //plmk2,y,x,incolor=incolor(i+1,),width=11,msize=0.3,marker=i+1;
      
      //if (pldelta) plg2,abs(y-x),x,marker,msize=msize,color=color(,i+1+2);
      //else plg2,y,x,marker,msize=msize,color=color(,i+1+2);

      if (pldelta) plmk2,abs(y-x),x,incolor=color(,i+1),width=11,msize=0.4,marker=i+1;
      else plmk2,y,x,incolor=color(,i+1),width=11,msize=0.4,marker=i+1;
    }
  
}
