
struct NDfield_info
{
  string comment;
  int fdims_index; //index where dims are the fields dims
  pointer x0;
  pointer delta;
} 

  func NDfield_read(fname,&header,last=,first=,headeronly=)
/* DOCUMENT 
   header.fdims_index is set to 1 if the the field represent coordinates of points, 0 if it reprersents values in an array.
   for instance, a dims=[3,4] array represent the 3D coordinates of 4 points if header.fdims_index=1, and the values of a 3x4
   if header.fdims_index=0.

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
  if(strchar(tag)(1)!="NDFIELD") error,"Unknown format";
  dummy = _grid_read(int);

  header=NDfield_info();
  
  dummy = _grid_read(int);
  comment =_grid_read(char, 80);
  header.comment=strchar(comment)(1);
  ndims =_grid_read(int);
  dims =_grid_read(int,ndims);dummy = _grid_read(int,20-ndims);
  header.fdims_index=_grid_read(int);
  if (header.fdims_index)
    {
      real_ndims=dims(1);
      dims=dims([1,2]);
    }
  else real_ndims=ndims;
  
  datatype=_grid_read(int);
  header.x0 = &_grid_read(double,real_ndims);dummy = _grid_read(double,20-real_ndims);
  header.delta =&_grid_read(double,real_ndims);dummy = _grid_read(double,20-real_ndims);
  dummy = _grid_read(char,160);
  dummy = _grid_read(int);
  if(headeronly) return header;
  if(datatype==(1<<0))  ndtype=char(0);
  if(datatype==(1<<1))  ndtype=char(0);
  if(datatype==(1<<2))  ndtype=short(0);
  if(datatype==(1<<3))  ndtype=short(0);
  if(datatype==(1<<4))  ndtype=int(0);
  if(datatype==(1<<5))  ndtype=int(0);
  if(datatype==(1<<6))  ndtype=long(0);
  if(datatype==(1<<7))  ndtype=long(0);
  if(datatype==(1<<8))  ndtype=float(0);
  if(datatype==(1<<9))  ndtype=double(0);

  if((sizeof(long)==4)*(ndtype=="long")) error,"long on 32 bit machine";
  dummy = _grid_read(int);
  tab=_grid_read(ndtype,grow(numberof(dims),dims));
  close,stream;
  return tab;
}


func NDfield_write(tab,fname,header,last=,first=, x0=, delta=, comment=, isCoords=, scaleBBox=)
/* DOCUMENT 
   set isCoords to 1 if the field store coordinates, set to 0 if it stores values in an array (easy way to set fdims_index).

   header.fdims_index is set to 1 if the the field represent coordinates of points, 0 if it reprersents values in an array.
   for instance, a dims=[3,4] array represent the 3D coordinates of 4 points if header.fdims_index=1, and the values of a 3x4
   if header.fdims_index=0.

   parameters x0, delta, comment and isCoord replace the values in header even if header is non void. 
SEE ALSO:
 */
{
  if (is_void(fname)) error," NDfield_write needs 2 arguments";
  dd=int(dimsof(tab));
  
  if(is_void(header))
    {
      header=NDfield_info();
      header.comment="Generated automatically";
      header.x0=&array(double,20);
      if (is_void(isCoords)) header.delta=&array(1.,20);
      else header.delta=&array(0.,20);
      header.fdims_index=0; // C indexing 
    }
  if (!is_void(isCoords))  header.fdims_index=1;
  if (!is_void(x0)) {
    tmp=*header.x0;
    tmp(:numberof(x0)) = x0;
    header.x0=&tmp ;
  }
  if (!is_void(delta)) {
    tmp=*header.delta;
    tmp(:numberof(delta)) = delta;
    header.delta=&tmp ;
  }

  if (!is_void(scaleBBox)) {
    tx0=*header.x0;
    tdelta=*header.delta;
    c=tx0+tdelta/2;
    tdelta*=scaleBBox;
    tx0=c-tdelta/2;
    header.x0=&tx0;
    header.delta=&tdelta;
  }
  
  if (!is_void(comment)) header.comment=comment;
  
  if (typeof(tab)=="complex") {write,"casting input to double"; tab=double(tab);}

  if (typeof(tab)=="char") datatype=(1<<1)
  if (typeof(tab)=="short") datatype=(1<<3);
  if (typeof(tab)=="int") datatype=(1<<5);
  if (typeof(tab)=="long") datatype=(1<<7);
  if (typeof(tab)=="float") datatype=(1<<8);
  if (typeof(tab)=="double") datatype=(1<<9);
  if ((sizeof(long)==4)&&(typeof(tab)=="long"))
    datatype=(1<<5);
  
  
  stream = _grid_new(fname, overwrite=1);
  address = 0;
  
  dummy =int(16);
  _grid_write,dummy;

  dummy=array(char,16);
  tag="NDFIELD";
  dummy(1:strlen(tag)+1)=strchar(tag);
  _grid_write,dummy;
  
  dummy =int(16);
  _grid_write,dummy;

  dummy=int(sizeof(int)*(20+3) + sizeof(double)*(2*20) + 160*sizeof(char));
  _grid_write,dummy;
  
  dummy=array(char,80);
  dummy(1:strlen(header.comment)+1)=strchar(header.comment);
  _grid_write,dummy;
  _grid_write,dd(1);
  dims=array(int,20);
  dims(1:dd(1)) = dd(2:);
  _grid_write,dims;
  _grid_write,header.fdims_index;
  _grid_write,int(datatype);
  x0=array(double,20);
  x0(:numberof(*header.x0))=*header.x0;
  _grid_write,x0;
  x0(:numberof(*header.delta))=*header.delta;
  _grid_write,x0;
  dummy=array(char,160);
  _grid_write,dummy;
  
  dummy=int(sizeof(int)*(20+3) + sizeof(double)*(2*20) + 160*sizeof(char));
  _grid_write,dummy;

  
  dummy=int(sizeof(tab));
  _grid_write,dummy;
  
  _grid_write,tab;
  
  _grid_write,dummy;
  close,stream;
  return fname;
}
