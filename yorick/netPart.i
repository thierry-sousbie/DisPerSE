func netPartID(net,x0=,dim=)
{
  if (is_void(dim)) dim=1;
  if (is_void(x0)) x0=0;

  p=*net.v_coord;
  id=array(long(-1),net.nvertex);
  w1=where(p(dim,)<x0);
  if (numberof(w1)==0) error,"";
  id(w1)=1;

  w2=where(p(dim,)>=x0);
  if (numberof(w2)==0) error,"";
  id(w2)=2;

  seg=*net.f_vertexIndex(1);

  w=where((p(dim,seg)(1,)<x0)!=(p(dim,seg)(2,)<x0));
  segw=seg(,w);
  vid=p(dim,segw)(mxx,);
  for (i=1;i<=numberof(vid);i++) id(segw(vid(i),i))=0;

  return id;
}

func netPart(net,x0=,dim=)
{
  if (is_void(dim)) dim=1;
  if (is_void(x0)) x0=0;

  id=netPartID(net,x0=x0,dim=dim);

  //w1=where((id==0)|(id==1));
  //w2=where((id==0)|(id==2));
  
  new_net=[];
  for (i=1;i<=2;i++) {
    w=where((id==0)|(id==i));
    nnet=net;
    p=(*net.v_coord)(,w);
    newID=array(int(0),net.nvertex);
    nnet.nvertex=numberof(w);
    newID(w)=indgen(numberof(w));
    nnet.v_coord=&p;
    newdata=*nnet.data;
    
    for (k=1;k<=nnet.ndata;k++)
      {
        data=newdata(k);
        if (data.type==0)
          {
            d=*data.data;
            data.data=&d(w);
            newdata(k)=data;
          }
      } 
    
    for (j=1;j<=net.ndims;j++)
      //for (j=net.ndims;j<=net.ndims;j++)
      {
        if (j!=net.ndims) {
          nnet.nfaces(j)=0;
          nnet.f_vertexIndex(j)=&[];
          continue;
        }
        fid=(*net.f_vertexIndex(j));
        w1=where((id(fid)(1,)==i)|(id(fid)(1,)==0));      
        for (k=2;k<=j+1;k++)
          {
            w2=where((id(fid)(k,w1)==i)|(id(fid)(k,w1)==0));
            w1=w1(w2);            
          }
        if ((j==net.ndims)&&(i==1))
          {
            wt=where(id(fid)(max,w1)!=0);
            if (numberof(wt)) w1=w1(wt);
          }     
        nnet.f_vertexIndex(j)=&newID(fid(,w1));

        for (k=1;k<=net.ndata;k++)
          {
            data=newdata(k);
            if (data.type==k)
              {
                d=*data.data;
                data.data=&d(w1);
                newdata(k)=data;
              }
          }

        nnet.nfaces(j)=numberof(w1);
      }
    nnet.data=&newdata;
    
    grow,new_net,&nnet;
  }
  return new_net;
}
