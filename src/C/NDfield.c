/*
 # Copyright (C) 2009-2013 Thierry Sousbie
 # University of Tokyo / CNRS, 2009
 #
 # This file is part of porject DisPerSE
 # 
 #  Author          : Thierry Sousbie
 #  Contact         : tsousbie@gmail.com	
 #
 #  Licenses        : This file is 'dual-licensed', you have to choose one
 #                    of the two licenses below to apply.
 #
 #                    CeCILL-C
 #                    The CeCILL-C license is close to the GNU LGPL.
 #                    ( http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html )
 #
 #                or  CeCILL v2.0
 #                    The CeCILL license is compatible with the GNU GPL.
 #                    ( http://www.cecill.info/licences/Licence_CeCILL_V2-en.html )
 #
 #  This software is governed either by the CeCILL or the CeCILL-C license
 #  under French law and abiding by the rules of distribution of free software.
 #  You can  use, modify and or redistribute the software under the terms of
 #  the CeCILL or CeCILL-C licenses as circulated by CEA, CNRS and INRIA
 #  at the following URL : "http://www.cecill.info".
 #
 #  As a counterpart to the access to the source code and  rights to copy,
 #  modify and redistribute granted by the license, users are provided only
 #  with a limited warranty  and the software's author,  the holder of the
 #  economic rights,  and the successive licensors  have only  limited
 #  liability.
 #
 #  In this respect, the user's attention is drawn to the risks associated
 #  with loading,  using,  modifying and/or developing or reproducing the
 #  software by the user in light of its specific status of free software,
 #  that may mean  that it is complicated to manipulate,  and  that  also
 #  therefore means  that it is reserved for developers  and  experienced
 #  professionals having in-depth computer knowledge. Users are therefore
 #  encouraged to load and test the software's suitability as regards their
 #  requirements in conditions enabling the security of their systems and/or
 #  data to be ensured and,  more generally, to use and operate it in the
 #  same conditions as regards security.
 #
 #  The fact that you are presently reading this means that you have had
 #  knowledge of the CeCILL and CeCILL-C licenses and that you accept its terms.
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "myendian.h"
#include "NDfield.h"

#ifdef HAVE_CFITS_IO
#include "fitsio.h"
#endif


char *print_dataType(char *dest, int type)
{
  switch (type)
    {
    case ND_UCHAR: return strcpy(dest,"UCHAR"); break;
    case ND_CHAR: return strcpy(dest,"CHAR"); break;
    case ND_USHORT: return strcpy(dest,"USHORT"); break;
    case ND_SHORT: return strcpy(dest,"SHORT"); break;
    case ND_UINT: return strcpy(dest,"UINT"); break;
    case ND_INT: return strcpy(dest,"INT"); break;
    case ND_ULONG: return strcpy(dest,"ULONG"); break;
    case ND_LONG: return strcpy(dest,"LONG"); break;
    case ND_FLOAT: return strcpy(dest,"FLOAT"); break;
    case ND_DOUBLE: return strcpy(dest,"DOUBLE"); break;
    }

  return dest;
}

int sizeof_NDfield(int type)
{
  switch (type)
    {
    case ND_UCHAR: return sizeof(unsigned char); break;
    case ND_CHAR: return sizeof(char); break;
    case ND_USHORT: return sizeof(unsigned short); break;
    case ND_SHORT: return sizeof(short); break;
    case ND_UINT: return sizeof(unsigned int); break;
    case ND_INT: return sizeof(int); break;
    case ND_ULONG: return sizeof(unsigned long); break;
    case ND_LONG: return sizeof(long); break;
    case ND_FLOAT: return sizeof(float); break;
    case ND_DOUBLE: return sizeof(double); break;
    }

  return 0;
}

int Get_NDtype(int size,int is_integer, int not_signed)
{
  if (size==sizeof(double))
    {
      if (is_integer)
	return ND_LONG;
      else
	return ND_DOUBLE;
    }
  else
    {
      if (is_integer)
	{
	  if (not_signed)
	    {
	      if (size==sizeof(char))
		return ND_UCHAR;
	      else if (size==sizeof(short))
		return ND_USHORT;
	      else if (size==sizeof(int))
		return ND_UINT;
	    }
	  else
	    {
	      if (size==sizeof(char))
		return ND_CHAR;
	      else if (size==sizeof(short))
		return ND_SHORT;
	      else if (size==sizeof(int))
		return ND_INT;
	    }
	}
      else
	return ND_FLOAT;
    }

  return -1;
}

NDfield *Create_NDfield(int *dims,int ndims,int fdims_index,int type,double *x0, double *delta, void *data,const char *comment)
{
  NDfield *field;
  int i;

  if (dims==NULL)
    {
      fprintf (stderr,"Error in Create_NDfield: dimensions needed.");
      return NULL;
    }
  if (ndims==0)
    {
      fprintf (stderr,"Error in Create_NDfield: numberof of dimensions must be >0.");
      return NULL;
    }
 
  field=calloc(1,sizeof(NDfield));
  if (fdims_index) memcpy(field->dims,dims,(ndims-fdims_index)*sizeof(int));
  else memcpy(field->dims,dims,ndims*sizeof(int));
  memcpy(field->dims,dims,ndims*sizeof(int));
  field->ndims = ndims;

  field->nval=1;
  if (fdims_index) field->nval=(long)field->dims[0]*(long)field->dims[1];
  else
    for (i=0;i<ndims;i++)
      field->nval *= field->dims[i];

  field->fdims_index = fdims_index;
  if (field->fdims_index) field->n_dims=2;
  else field->n_dims=field->ndims;
  /*
  if (fdims_index>=0)
    field->fdims_index = fdims_index;
  else
    field->fdims_index = ndims-1;
  */
  
  field->datatype=type;
  if ((type==ND_LONG)&&(sizeof(long)==4)) field->datatype=ND_INT;
  if ((type==ND_ULONG)&&(sizeof(long)==4)) field->datatype=ND_UINT;

  field->datasize = sizeof_NDfield(field->datatype);

  if (x0!=NULL)
    memcpy(field->x0,x0,sizeof(double)*ndims);
  if (delta!=NULL)
    memcpy(field->delta,delta,sizeof(double)*ndims);

  if (data==NULL)
    field->val = calloc(field->nval,field->datasize);
  else
    field->val=data;

  if (comment!=NULL)
    strcpy(field->comment,comment);

  return field;
}

int Init_NDfield(NDfield *field,int *dims,int ndims,int fdims_index,int type, double *x0,double *delta, void *data,char *comment)
{
  int i;

  memset(field,0,sizeof(NDfield));

  if (dims==NULL)
    {
      fprintf (stderr,"Error in Create_NDfield: dimensions needed.");
      return -1;
    }
  if (ndims==0)
    {
      fprintf (stderr,"Error in Create_NDfield: numberof of dimensions must be >0.");
      return -1;
    }
 
  memcpy(field->dims,dims,ndims*sizeof(int));
  field->ndims = ndims;

  field->nval=1;
  for (i=0;i<ndims;i++)
    field->nval *= field->dims[i];

  field->fdims_index = fdims_index;
  /*
  if (fdims_index>=0)
    field->fdims_index = fdims_index;
  else
    field->fdims_index = ndims-1;
  */
  
  field->datatype=type;
  if ((type==ND_LONG)&&(sizeof(long)==4)) type=ND_INT;
  if ((type==ND_ULONG)&&(sizeof(long)==4)) type=ND_UINT;

  field->datasize = sizeof_NDfield(field->datatype);

  if (x0!=NULL)
    memcpy(field->x0,x0,sizeof(double)*ndims);
  if (delta!=NULL)
    memcpy(field->delta,delta,sizeof(double)*ndims);

  if (data==NULL)
    field->val = calloc(field->nval,field->datasize);
  else
    field->val=data;

  if (comment!=NULL)
    strcpy(field->comment,comment);

  return 0;
}

int Free_NDfield(NDfield **field)
{
  if (*field == NULL) return -1;

  free((*field)->val);
  free(*field);
  *field=NULL;

  return 0;
}

int Save_NDfield(NDfield *field,const char *filename)
{
  int i;
  char tag[16];
  char dummy[160];
  FILE *f;

  char dims_msg[255];
  sprintf(dims_msg,"[%d",field->dims[0]);
  for (i=1;i<field->n_dims;i++) sprintf(dims_msg,"%s,%d",dims_msg,field->dims[i]);
  sprintf(dims_msg,"%s]",dims_msg);
  
  char msg[255];
  print_dataType(msg,field->datatype);
  if (field->fdims_index!=0) sprintf(msg,"%s coords",msg);
  else sprintf(msg,"%s array",msg);
  printf ("Saving %s %s to file %s ...",dims_msg,msg,filename);
  fflush(0);

  //printf ("Saving %dD field to file %s ...",field->ndims,filename);fflush(0);

  memset(tag,0,16*sizeof(char));
  memset(dummy,0,160*sizeof(char));
  strcpy(tag,NDFIELD_TAG);
  i=16;

  f=fopen(filename,"w");

  fwrite(&i,sizeof(int),1,f);
  fwrite(tag,sizeof(char),16,f);
  fwrite(&i,sizeof(int),1,f);

  i=sizeof(int)*(NDFIELD_MAX_DIMS+3) + sizeof(double)*(2*NDFIELD_MAX_DIMS) + (160+80)*sizeof(char);
  fwrite(&i,sizeof(int),1,f);
  fwrite(field->comment,sizeof(char),80,f);
  fwrite(&field->ndims,sizeof(int),1,f);
  fwrite(field->dims,sizeof(int),field->ndims,f);
  if (field->ndims<NDFIELD_MAX_DIMS) fwrite(dummy,sizeof(int),NDFIELD_MAX_DIMS-field->ndims,f);
  fwrite(&field->fdims_index,sizeof(int),1,f);
  fwrite(&field->datatype,sizeof(int),1,f);
  fwrite(field->x0,sizeof(double),field->ndims,f);
  if (field->ndims<NDFIELD_MAX_DIMS) fwrite(dummy,sizeof(double),NDFIELD_MAX_DIMS-field->ndims,f);
  fwrite(field->delta,sizeof(double),field->ndims,f);
  if (field->ndims<NDFIELD_MAX_DIMS) fwrite(dummy,sizeof(double),NDFIELD_MAX_DIMS-field->ndims,f);
  fwrite(field->dummy,sizeof(char),160,f);
  fwrite(&i,sizeof(int),1,f);
  
  i=field->nval*sizeof_NDfield(field->datatype);
  fwrite(&i,sizeof(int),1,f);
  fwrite(field->val,sizeof_NDfield(field->datatype),field->nval,f);
  fwrite(&i,sizeof(int),1,f);

  fclose(f);
  printf (" done.\n");
  return 0;
}

int LoadDensity_CIC(const char *fname,OLD_density_grid *density)
{
    FILE *f;
    unsigned int i;
    char test[30];
    int swap = 0;
    long offs=0;
    
    if(!(f = fopen(fname,"r")))
    {
	fprintf(stderr,"File %s does not exist.\n",fname);
	return 1;
    }
    int ret;

    ret=fread(&i,sizeof(int),1,f);
    ret=fread(test,sizeof(char)*30,1,f);
    ret=fread(&i,sizeof(int),1,f);

    if (i!=30*sizeof(char)) swap=1;
 
    if (strcmp(test,"Density grid file")) 
	{
	    printf("-->");
	    return -1;
	}
    else printf ("Loading density grid from %s ... ",fname);fflush(0);

    if (swap) printf ("(Swapping)");fflush(0);
    
    ret=fread(&i,sizeof(int),1,f);
    fread_sw(&(density->Nx),sizeof(int),1,f,swap);
    fread_sw(&(density->Ny),sizeof(int),1,f,swap);
    fread_sw(&(density->Nz),sizeof(int),1,f,swap);
    fread_sw(&(density->NNodes),sizeof(int),1,f,swap);

    fread_sw(&(density->x0),sizeof(float),1,f,swap);
    fread_sw(&(density->y0),sizeof(float),1,f,swap);
    fread_sw(&(density->z0),sizeof(float),1,f,swap);
    
    fread_sw(&(density->dx),sizeof(float),1,f,swap);
    fread_sw(&(density->dy),sizeof(float),1,f,swap);
    fread_sw(&(density->dz),sizeof(float),1,f,swap);
    ret=fread(&i,sizeof(int),1,f);

    if (density->NNodes != (INT)density->Nx*(INT)density->Ny*(INT)density->Nz)
	density->NNodes = (INT)density->Nx*(INT)density->Ny*(INT)density->Nz;

    if (!(density->grid=(FLOAT *)malloc((size_t)density->NNodes*sizeof(FLOAT))))
    {
	fprintf(stderr,"Not enough memory for density->grid while loading\n");
	exit(0);
    }
    
    ret=fread(&i,sizeof(int),1,f);
    for (i=0;i<density->Nx;i++)
      {
	fread_sw(&(density->grid[offs]),sizeof(FLOAT),density->NNodes/density->Nx,f,swap);
	offs += density->NNodes/density->Nx;
      }
    ret=fread(&i,sizeof(int),1,f);
    
    i=sizeof(float);
    ret= fread(&i,sizeof(int),1,f);
    fread_sw(&density->redshift,sizeof(float),1,f,swap);
    fread_sw(&density->smoothing,sizeof(float),1,f,swap);
    ret=fread(&i,sizeof(int),1,f);
    
    density->HasGradiant=0;
    printf ("done.\n");
    return 0;
}



inline int C2I(INT *coord, INT *index,int *dims, int ndims,int periodic)
{
  INT val;
  INT dec;
  int i;
  int out=0;

  *index=0;
  dec=1;
  for (i=0;i<ndims;i++)
    {   
      val=coord[i];
      if (periodic&(1<<i))
	{
	  if (val==dims[i]) val=0;
	  if (val>=dims[i]) val = val%(dims[i]);
	  //else if (val<0) val = (100*(dims[i]) + val)%(dims[i]);
	  else if (val<0) val = dims[i]-((-val)%dims[i]);
	}
      else
	{
	    if (val>=dims[i]) {val = dims[i]-1;out=1;}
	    if (val<0) {val=0;out=1;}
	}
      *index += val*dec;
      dec*=dims[i];
    }
  
  return out;
}

int Save_NDfieldPartial(char *filename, NDfield *header, NDfield *field, int periodic)
{
  int i,k;
  long n;
  long j;
  char tag[16];
  char dummy[160];
  FILE *f;
  long pos;
  long long tmp=0;
  INT index,old_index;
  double coord[NDFIELD_MAX_DIMS];
  INT coordi[NDFIELD_MAX_DIMS];
  INT x0[NDFIELD_MAX_DIMS];
  int elSize;
  int newfile=0;

  printf ("Saving partial %dD field to file %s ...",field->ndims,filename);fflush(0);

  memset(tag,0,16*sizeof(char));
  memset(dummy,0,160*sizeof(char));
  strcpy(tag,NDFIELD_TAG);
  i=16;

  if ((f=fopen(filename,"r"))==NULL) 
    newfile=1;
  else fclose(f);

  if (!newfile) 
    f=fopen(filename,"r+");
  else
    f=fopen(filename,"w");

  fseek(f,0,SEEK_SET);

  fwrite(&i,sizeof(int),1,f);
  fwrite(tag,sizeof(char),16,f);
  fwrite(&i,sizeof(int),1,f);

  i=sizeof(int)*(NDFIELD_MAX_DIMS+3) + sizeof(double)*(2*NDFIELD_MAX_DIMS) + (160+80)*sizeof(char);
  fwrite(&i,sizeof(int),1,f);
  fwrite(header->comment,sizeof(char),80,f);
  fwrite(&header->ndims,sizeof(int),1,f);
  fwrite(header->dims,sizeof(int),header->ndims,f);
  if (header->ndims<NDFIELD_MAX_DIMS) fwrite(dummy,sizeof(int),NDFIELD_MAX_DIMS-header->ndims,f);
  fwrite(&header->fdims_index,sizeof(int),1,f);
  fwrite(&header->datatype,sizeof(int),1,f);
  fwrite(header->x0,sizeof(double),header->ndims,f);
  if (header->ndims<NDFIELD_MAX_DIMS) fwrite(dummy,sizeof(double),NDFIELD_MAX_DIMS-header->ndims,f);
  fwrite(header->delta,sizeof(double),header->ndims,f);
  if (header->ndims<NDFIELD_MAX_DIMS) fwrite(dummy,sizeof(double),NDFIELD_MAX_DIMS-header->ndims,f);
  fwrite(header->dummy,sizeof(char),160,f);
  fwrite(&i,sizeof(int),1,f);

  elSize = sizeof_NDfield(field->datatype);
  for (i=0,n=1;i<header->ndims;i++) n*=(long)header->dims[i];

  i=header->nval*sizeof_NDfield(header->datatype);
  fwrite(&i,sizeof(int),1,f);
  
  if (newfile)
    {
      printf(" (new file)");fflush(0);
      pos = ftell(f);
      for (j=0;j<n;j++) 
	fwrite(&tmp,elSize,1,f);
      i=header->nval*sizeof_NDfield(header->datatype);
      fwrite(&i,sizeof(int),1,f);
      fseek(f,pos,SEEK_SET);
    }

  //do the partial write HERE ...
  for (i=0;i<header->ndims;i++) coord[i]=field->x0[i];
  for (i=0;i<header->ndims;i++) x0[i]=coordi[i] = header->dims[i]*(coord[i]-header->x0[i])/header->delta[i];
  for (i=0,n=1;i<field->ndims;i++) n*=field->dims[i];

  j=0;
  old_index=0;
  
  for (j=0;j<n*elSize;j+=elSize)
    {
	k=C2I(coordi,&index,header->dims,header->ndims,periodic);
	
	if (index-old_index) fseek(f,elSize*(index-old_index),SEEK_CUR);
	if (!k) fwrite(field->val + j,elSize,1,f);
	
	i=0;coordi[i]++;
	if (*coordi>=x0[i]+field->dims[i])
	    do {coordi[i]=x0[i];coordi[++i]++;} while(coordi[i]>=x0[i]+header->dims[i]);
	old_index=index+1;
    } 
 
  fclose(f);

  printf (" done.\n");

  return 1;
}

NDfield *Load_NDfieldHeader(const char *filename)
{
  int i;
  char tag[16];
  char dummy[160];
  FILE *f;
  int swap=0;
  NDfield *field;

  char comment[80];
  int dims[NDFIELD_MAX_DIMS];
  int ndims;
  int fdims_index; //index where dims are the fields dims
  int datatype;
  double x0[NDFIELD_MAX_DIMS];
  double delta[NDFIELD_MAX_DIMS];
  //int isgrafic=0;
 

  memset(tag,0,16*sizeof(char));
  memset(dummy,0,160*sizeof(char));
  strcpy(tag,NDFIELD_TAG);
  i=16;

  if(!(f = fopen(filename,"r")))
    {
	fprintf(stderr,"File %s does not exist.\n",filename);
	return NULL;
    }

  fread_sw(&i,sizeof(int),1,f,swap);
  fclose(f);
  
  if (i!=16) swap=1-swap;

  f=fopen(filename,"r");
  fread_sw(&i,sizeof(int),1,f,swap);
  fread_sw(tag,sizeof(char),16,f,swap);
  fread_sw(&i,sizeof(int),1,f,swap);
  tag[15]='\0';
  
  if (strcmp(tag,NDFIELD_TAG))
    {
      float tmpf[8];
      fclose(f);
      f=fopen(filename,"r");
      fread_sw(&i,sizeof(int),1,f,swap);
      if (i!=44) {
	swap=1-swap;
	fclose(f);f=fopen(filename,"r");
	fread_sw(&i,sizeof(int),1,f,swap);
      }

      if (i==44)
	{
	  fread_sw(dims,sizeof(int),3,f,swap);
	  fread_sw(tmpf,sizeof(float),8,f,swap);
	  fread_sw(&i,sizeof(int),1,f,swap);
	  strcpy(comment,"GRAFIC file");
	  ndims=3;fdims_index=0;
	  datatype=ND_FLOAT;
	  for (i=0;i<3;i++) delta[i] = tmpf[0]*dims[i]*tmpf[7]/100.;
	  for (i=0;i<3;i++) x0[i]=0;
	  //isgrafic=1;
	}
      else
	{
	  fclose(f);
	  fprintf (stderr,"File %s has an unknown format.\n",filename);
	  return NULL;
	}
    }
  else
    {
      fread_sw(&i,sizeof(int),1,f,swap);
      fread_sw(comment,sizeof(char),80,f,swap);
      fread_sw(&ndims,sizeof(int),1,f,swap);
      fread_sw(dims,sizeof(int),NDFIELD_MAX_DIMS,f,swap);
      fread_sw(&fdims_index,sizeof(int),1,f,swap);
      fread_sw(&datatype,sizeof(int),1,f,swap);
      fread_sw(x0,sizeof(double),NDFIELD_MAX_DIMS,f,swap);
      fread_sw(delta,sizeof(double),NDFIELD_MAX_DIMS,f,swap);
      fread_sw(dummy,sizeof(char),160,f,swap);
      fread_sw(&i,sizeof(int),1,f,swap);
    }

  field = Create_NDfield(dims,ndims,fdims_index,datatype,x0,delta,(void*)(&i),comment);
  field->val=NULL;
  fclose(f);

  memcpy(field->dummy,dummy,sizeof(char)*160);

  return field;
}

NDfield *Load_NDfieldChunk(char *filename, double *x0, double *delta, int periodic)
{
  char tag[16];
  char dummy[160];
  FILE *f;
  int swap=0;
  char comment[80];
  int ldims[NDFIELD_MAX_DIMS];
  int lndims;
  int fdims_index; //index where dims are the fields dims
  int datatype;
  double lx0[NDFIELD_MAX_DIMS];
  double ldelta[NDFIELD_MAX_DIMS];

  NDfield *field;
  int i;
  INT index;
  INT oldindex;
  int imin[NDFIELD_MAX_DIMS];
  int imax[NDFIELD_MAX_DIMS];
  int deltai[NDFIELD_MAX_DIMS];
  double deltax[NDFIELD_MAX_DIMS];
  double xmin[NDFIELD_MAX_DIMS];
  double xmax[NDFIELD_MAX_DIMS];
  INT coord[NDFIELD_MAX_DIMS];
  double tmp;
  INT tmpi;
  INT n;
  int elSize;
  int isgrafic=0;
  
  memset(tag,0,16*sizeof(char));
  memset(dummy,0,160*sizeof(char));
  strcpy(tag,NDFIELD_TAG);
  i=16;

  if(!(f = fopen(filename,"r")))
    {
	fprintf(stderr,"File %s does not exist.\n",filename);
	return NULL;
    }
 
  fread_sw(&i,sizeof(int),1,f,swap);
  fclose(f);
  
  if (i!=16) swap=1-swap;

  f=fopen(filename,"r");
  fread_sw(&i,sizeof(int),1,f,swap);
  fread_sw(tag,sizeof(char),16,f,swap);
  fread_sw(&i,sizeof(int),1,f,swap);
  tag[15]='\0';
  
  if (strcmp(tag,NDFIELD_TAG))
    {
      fclose(f);
      if ((field=Load_NDfieldHeader(filename))==NULL) return NULL;
      strcpy(comment,field->comment);
      lndims=field->ndims;
      memcpy(ldims,field->dims,sizeof(int)*NDFIELD_MAX_DIMS);
      memcpy(lx0,field->x0,sizeof(float)*NDFIELD_MAX_DIMS);
      memcpy(ldelta,field->delta,sizeof(float)*NDFIELD_MAX_DIMS);
      memcpy(dummy,field->dummy,sizeof(char)*160);
      fdims_index=field->fdims_index;
      datatype=field->datatype;
      f=fopen(filename,"r");
      fseek(f,44+3*sizeof(int),SEEK_SET);  
      isgrafic=1;
    }
  else
    {
      fread_sw(&i,sizeof(int),1,f,swap);
      fread_sw(comment,sizeof(char),80,f,swap);
      fread_sw(&lndims,sizeof(int),1,f,swap);
      fread_sw(ldims,sizeof(int),NDFIELD_MAX_DIMS,f,swap);
      fread_sw(&fdims_index,sizeof(int),1,f,swap);
      fread_sw(&datatype,sizeof(int),1,f,swap);
      fread_sw(lx0,sizeof(double),NDFIELD_MAX_DIMS,f,swap);
      fread_sw(ldelta,sizeof(double),NDFIELD_MAX_DIMS,f,swap);
      fread_sw(dummy,sizeof(char),160,f,swap);
      fread_sw(&i,sizeof(int),1,f,swap);
      fread_sw(&i,sizeof(int),1,f,swap);
    }
  /*
  for (i=0;i<lndims;i++)
    if (ldelta[i] <= delta[i])
      {
	fprintf (stderr,"ERROR in Load_NDfieldChunk, chunk larger than array.\n");
	return NULL;
      }
  */
  for (i=0;i<lndims;i++)
    {
      tmpi = (int)(1.E-8 + (x0[i]-lx0[i])/(ldelta[i])*ldims[i]) - 1;
      //printf("%d: [%e %e %e-> %d]",i,(x0[i]-lx0[i]),(ldelta[i]),((x0[i]-lx0[i])/(ldelta[i])*ldims[i]),tmpi);
      
      imin[i] = tmpi;
      xmin[i] = lx0[i] +  tmpi * (ldelta[i]/ldims[i]);

      tmpi = (int)(1.E-8 + (x0[i]+delta[i]-lx0[i])/(ldelta[i])*ldims[i]) + 1;
      //printf(" [%e %e %e-> %d]\n",(x0[i]+delta[i]-lx0[i]),(ldelta[i]),((x0[i]+delta[i]-lx0[i])/(ldelta[i])*ldims[i]),tmpi);      
      
      imax[i] = tmpi;
      xmax[i] = lx0[i] +  tmpi * (ldelta[i]/ldims[i]);
         
      if (periodic&(1<<i))
	{
	  tmp = (xmin[i]-lx0[i])/ ldelta[i];
	  tmpi = - ((int)tmp);
	  if (tmpi)
	    {
	      imax[i] += tmpi * ldims[i];
	      imin[i] += tmpi * ldims[i];
	      xmax[i] += tmpi * ldelta[i];
	      xmin[i] += tmpi * ldelta[i];
	    }
	}
      else
	{
	  if (imax[i] > ldims[i]) {imax[i]=ldims[i];xmax[i]=lx0[i]+ldelta[i];}
	  if (imin[i] < 0) {imin[i]=0;xmin[i]=lx0[i];}
	}
      
      deltai[i] = imax[i]-imin[i];
      deltax[i] = xmax[i]-xmin[i];
    }


  
  tmpi=1;
  for (i=0;i<lndims;i++) coord[i]=imin[i];
  for (i=0;i<lndims;i++) tmpi*=deltai[i];

  /*
  for (i=0;i<lndims;i++) printf ("[%f,%f]",x0[i],x0[i]+delta[i]);
  printf("\n");
  for (i=0;i<lndims;i++) printf ("[%f,%f]",xmin[i],xmin[i]+deltax[i]);
  printf("\n");
  for (i=0;i<lndims;i++) printf ("[%d,%d]",imin[i],imax[i]);
  printf("\n");
  */

  field = Create_NDfield(deltai,lndims,fdims_index,datatype,xmin,deltax,NULL,comment);
  memcpy(field->dummy,dummy,160*sizeof(char));
  //for (i=0;i<2;i++) printf ("final %f %f\n",xmin[i],deltax[i]);
  //for (i=0;i<2;i++) printf ("cut %f %f\n",x0[i],delta[i]);
  //for (i=0;i<2;i++) printf ("init %f %f\n",lx0[i],ldelta[i]);
  

  oldindex=0;
  elSize = sizeof_NDfield(field->datatype);
  printf ("Loading %dD chunk from file %s ...",lndims,filename);fflush(0);
  for(n=0;n<tmpi*elSize;n+=elSize)
    {
      C2I(coord,&index,ldims,lndims,periodic);
      
      if (index-oldindex) fseek(f,elSize*(index-oldindex),SEEK_CUR);

      if (isgrafic)
	{
	  int a,b,c;
	  c=ldims[1]*ldims[0];
	  a=(int)(index/c);
	  b=(int)(oldindex/c);
	  
	  if (a!=b) 
	    fseek(f,sizeof(int)*2*(a-b),SEEK_CUR);
	  
	}

      fread_sw(field->val + n,elSize,1,f,swap);
      i=0;coord[i]++;
      if (*coord>=*imax)
	do {coord[i]=imin[i];coord[++i]++;} while(coord[i]>=imax[i]);
      oldindex=index+1;
    }
 
  fclose(f);
  printf(" done.\n");
  return field;
}

NDfield *Load_NDfieldChunkHeader(char *filename, double *x0, double *delta, int periodic)
{
  char tag[16];
  char dummy[160];
  FILE *f;
  int swap=0;
  char comment[80];
  int ldims[NDFIELD_MAX_DIMS];
  int lndims;
  int fdims_index; //index where dims are the fields dims
  int datatype;
  double lx0[NDFIELD_MAX_DIMS];
  double ldelta[NDFIELD_MAX_DIMS];

  NDfield *field;
  int i;
  int imin[NDFIELD_MAX_DIMS];
  int imax[NDFIELD_MAX_DIMS];
  int deltai[NDFIELD_MAX_DIMS];
  double deltax[NDFIELD_MAX_DIMS];
  double xmin[NDFIELD_MAX_DIMS];
  double xmax[NDFIELD_MAX_DIMS];
  double tmp;
  INT tmpi;
  //int isgrafic=0;
    
  memset(tag,0,16*sizeof(char));
  memset(dummy,0,160*sizeof(char));
  strcpy(tag,NDFIELD_TAG);
  i=16;

  if(!(f = fopen(filename,"r")))
    {
	fprintf(stderr,"File %s does not exist.\n",filename);
	return NULL;
    }
 
  fread_sw(&i,sizeof(int),1,f,swap);
  fclose(f);
  
  if (i!=16) swap=1-swap;

  f=fopen(filename,"r");
  fread_sw(&i,sizeof(int),1,f,swap);
  fread_sw(tag,sizeof(char),16,f,swap);
  fread_sw(&i,sizeof(int),1,f,swap);
  tag[15]='\0';
  
  if (strcmp(tag,NDFIELD_TAG))
    {
      fclose(f);
      if ((field=Load_NDfieldHeader(filename))==NULL) return NULL;
      strcpy(comment,field->comment);
      lndims=field->ndims;
      memcpy(ldims,field->dims,sizeof(int)*NDFIELD_MAX_DIMS);
      memcpy(lx0,field->x0,sizeof(float)*NDFIELD_MAX_DIMS);
      memcpy(ldelta,field->delta,sizeof(float)*NDFIELD_MAX_DIMS);
      memcpy(dummy,field->dummy,sizeof(char)*160);
      fdims_index=field->fdims_index;
      datatype=field->datatype;
      f=fopen(filename,"r");
      fseek(f,44+3*sizeof(int),SEEK_SET);  
      //isgrafic=1;
    }
  else
    {
      fread_sw(&i,sizeof(int),1,f,swap);
      fread_sw(comment,sizeof(char),80,f,swap);
      fread_sw(&lndims,sizeof(int),1,f,swap);
      fread_sw(ldims,sizeof(int),NDFIELD_MAX_DIMS,f,swap);
      fread_sw(&fdims_index,sizeof(int),1,f,swap);
      fread_sw(&datatype,sizeof(int),1,f,swap);
      
      fread_sw(lx0,sizeof(double),NDFIELD_MAX_DIMS,f,swap);
      fread_sw(ldelta,sizeof(double),NDFIELD_MAX_DIMS,f,swap);
      fread_sw(dummy,sizeof(char),160,f,swap);
      fread_sw(&i,sizeof(int),1,f,swap);
      fread_sw(&i,sizeof(int),1,f,swap);
    }

  for (i=0;i<lndims;i++)
    {
      tmpi = (int)(1.E-8 + (x0[i]-lx0[i])/(ldelta[i])*ldims[i]) - 1;
      imin[i] = tmpi;
      xmin[i] = lx0[i] +  tmpi * (ldelta[i]/ldims[i]);
      
      tmpi = (int)(1.E-8 + (x0[i]+delta[i]-lx0[i])/(ldelta[i])*ldims[i]) + 1;
      imax[i] = tmpi;
      xmax[i] = lx0[i] +  tmpi * (ldelta[i]/ldims[i]);
      
    
      if (periodic&(1<<i))
	{
	  tmp = (xmin[i]-lx0[i])/ ldelta[i];
	  tmpi = - ((int)tmp);
	  if (tmpi)
	    {
	      imax[i] += tmpi * ldims[i];
	      imin[i] += tmpi * ldims[i];
	      xmax[i] += tmpi * ldelta[i];
	      xmin[i] += tmpi * ldelta[i];
	    }
	}
      else
	{
	  if (imax[i] > ldims[i]) {imax[i]=ldims[i];xmax[i]=lx0[i]+ldelta[i];}
	  if (imin[i] < 0) {imin[i]=0;xmin[i]=lx0[i];}
	}
      
      deltai[i] = imax[i]-imin[i];
      deltax[i] = xmax[i]-xmin[i];
    }

  field = Create_NDfield(deltai,lndims,fdims_index,datatype,xmin,deltax,(void *)(&tmpi),comment);
  memcpy(field->dummy,dummy,160*sizeof(char));
  field->val=NULL;
  
  fclose(f);
  
  return field;
}

int IsNDfield(const char *filename)
{
  int i;
  char tag[16];
  char comment[80];
  int ndims;
  FILE *f;
  int swap=0;
  
  memset(tag,0,16*sizeof(char));
   
  if(!(f = fopen(filename,"r")))
    {
	fprintf(stderr,"File %s does not exist.\n",filename);
	return 0;
    }

  fread_sw(&i,sizeof(int),1,f,swap);
  if (i!=16) swap=1-swap;
  fread_sw(tag,sizeof(char),16,f,swap);
  fread_sw(&i,sizeof(int),1,f,swap);
  
  
  tag[15]='\0';
  if (strcmp(tag,NDFIELD_TAG)) {fclose(f); return 0;}

  fread_sw(&i,sizeof(int),1,f,swap);
  fread_sw(comment,sizeof(char),80,f,swap);
  fread_sw(&ndims,sizeof(int),1,f,swap);
  fclose(f);
  return ndims;
}
#ifdef HAVE_CFITS_IO
void NDFIELD_FITS_printerror( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/

    if (status)
    {
       fits_report_error(stderr, status); /* print error report */

       exit( status );    /* terminate the program, returning error status */
    }
    return;
}
#endif

int IsNDfieldFromFITS(const char *filename)
{
#ifdef HAVE_CFITS_IO
  int status = 0;
  fitsfile *fptr;
  int nbpp;
  int ndims;
  long dimsl[20];

  if (!fits_open_file(&fptr, filename, READONLY, &status))
    {
      int res = fits_get_img_param(fptr,20,&nbpp,&ndims,dimsl,&status);
      fits_close_file(fptr, &status);
      if (res) return 0;
      else return ndims;
    }
  return 0;
#endif
  return 0;
}

NDfield *Load_NDfieldFromFITS(const char *filename)
{
#ifdef HAVE_CFITS_IO
  int status = 0;
  fitsfile *fptr;
  int ndims;
  int nbpp;
  int dims[20];
  long dimsl[20];
  double x0[20];
  double delta[20];
  long fp[20];
  long lp[20];
  long inc[20];
  long i;
  double *data;
  long npix=1;
  NDfield *field;
  
  if (!fits_open_file(&fptr, filename, READONLY, &status))
    {		
      if (fits_get_img_param(fptr,20,&nbpp,&ndims,dimsl,&status)) {
	NDFIELD_FITS_printerror( status );
	exit(-1);
      }
      //printf("ndims=%d\n",ndims);
      for (i=0;i<ndims;i++) {
	x0[i]=0;
	delta[i]=dimsl[i];
	dims[i]=dimsl[i];
	fp[i]=1;lp[i]=dimsl[i];inc[i]=1;
	npix*=dimsl[i];
      }
      if (dimsl[ndims-1]==1) ndims--;

      char dims_msg[255];
      sprintf(dims_msg,"[%d",dims[0]);
      for (i=1;i<ndims;i++) sprintf(dims_msg,"%s,%d",dims_msg,dims[i]);
      sprintf(dims_msg,"%s]",dims_msg);
      
      printf ("Loading %s array from file %s ... ",dims_msg,filename);fflush(0);
      
      int anynul;
      data=(double*)malloc(sizeof(double)*npix);
      if (fits_read_subset(fptr,TDOUBLE,fp,lp,inc,NULL,data,&anynul,&status))
	{
	  NDFIELD_FITS_printerror( status );
	  exit(-1);
	}	
      
      field = Create_NDfield(dims,ndims,0,ND_DOUBLE,x0,delta,(void *)data,"Loaded from FITS file");
      printf("done.\n");
    }
  else
    {
      fprintf (stderr,"Error loading file %s.\n",filename);
      NDFIELD_FITS_printerror( status );
      exit(-1);
    }

  return field;
#endif
  fprintf(stderr,"ERROR: CfitsIO library is needed to load '%s'.\n",filename);
  fprintf (stderr," install CFitsIO and recompile ...\n");
  exit(0);
  return NULL;
}

NDfield *Load_NDfield(const char *filename)
{
  int i;
  char tag[16];
  char dummy[160];
  FILE *f;
  int swap=0;
  NDfield *field;

  char comment[80];
  int dims[NDFIELD_MAX_DIMS];
  int ndims;
  int fdims_index; //index where dims are the fields dims
  int datatype;
  double x0[NDFIELD_MAX_DIMS];
  double delta[NDFIELD_MAX_DIMS];
  OLD_density_grid density;

  memset(tag,0,16*sizeof(char));
  memset(dummy,0,160*sizeof(char));
  strcpy(tag,NDFIELD_TAG);
  i=16;

  if(!(f = fopen(filename,"r")))
    {
	fprintf(stderr,"File %s does not exist.\n",filename);
	return NULL;
    }

  fread_sw(&i,sizeof(int),1,f,swap);
  fclose(f);
  
  if (i!=16) swap=1-swap;

  f=fopen(filename,"r");
  fread_sw(&i,sizeof(int),1,f,swap);
  fread_sw(tag,sizeof(char),16,f,swap);
  fread_sw(&i,sizeof(int),1,f,swap);
  tag[15]='\0';
  
  if (strcmp(tag,NDFIELD_TAG))
    {
      fclose(f);      
      if (LoadDensity_CIC(filename,&density))
	{
	  fprintf (stderr,"File %s has an unknown format.\n",filename);
	  return NULL;
	}
      if (sizeof(FLOAT)==8) datatype=ND_DOUBLE;
      else datatype=ND_FLOAT;
      strcpy(comment,"CIC format");
      ndims=3;dims[0]=density.Nx;dims[1]=density.Ny;dims[2]=density.Nz;
      fdims_index=0;
      x0[0]=density.x0;x0[1]=density.y0;x0[2]=density.z0;
      delta[0]=density.dx*density.Nx;delta[1]=density.dy*density.Ny;delta[2]=density.dz*density.Nz;
      field = Create_NDfield(dims,ndims,fdims_index,datatype,x0,delta,density.grid,comment);
      //printf ("d=%e %e\n",x0[0],delta[1]);
      return field;
    }

  fread_sw(&i,sizeof(int),1,f,swap);
  fread_sw(comment,sizeof(char),80,f,swap);
  fread_sw(&ndims,sizeof(int),1,f,swap);
  fread_sw(dims,sizeof(int),NDFIELD_MAX_DIMS,f,swap);
  fread_sw(&fdims_index,sizeof(int),1,f,swap);
  fread_sw(&datatype,sizeof(int),1,f,swap);
  fread_sw(x0,sizeof(double),NDFIELD_MAX_DIMS,f,swap);
  fread_sw(delta,sizeof(double),NDFIELD_MAX_DIMS,f,swap);
  fread_sw(dummy,sizeof(char),160,f,swap);
  fread_sw(&i,sizeof(int),1,f,swap);
  
  field = Create_NDfield(dims,ndims,fdims_index,datatype,x0,delta,NULL,comment);
  
  char dims_msg[255];
  sprintf(dims_msg,"[%d",dims[0]);
  for (i=1;i<field->n_dims;i++) sprintf(dims_msg,"%s,%d",dims_msg,dims[i]);
  sprintf(dims_msg,"%s]",dims_msg);
  
  char msg[255];
  print_dataType(msg,datatype);
  if (fdims_index!=0) sprintf(msg,"%s coords",msg);
  else sprintf(msg,"%s array",msg);
  printf ("Loading %s %s from file %s ...",dims_msg,msg,filename);
  fflush(0);

  if (swap) printf ("(swapping)\n");
  
  memcpy(field->dummy,dummy,160*sizeof(char));

  fread_sw(&i,sizeof(int),1,f,swap);
  fread_sw(field->val,sizeof_NDfield(field->datatype),field->nval,f,swap);
  fread_sw(&i,sizeof(int),1,f,swap);
  
  fclose(f);
  printf (" done.\n");
  return field;
}

int Convert_NDfield(NDfield *field,int type)
{
  void *new_tab;
  long i;

  DEFINE_NDVARS(ptr);
  DEFINE_NDVARS(new);
  /*
  unsigned char *ptruc;
  char *ptrc;
  unsigned short *ptrus;
  short *ptrs;
  unsigned int*ptrui;
  int *ptri;
  unsigned long *ptrul;
  long *ptrl;
  float *ptrf;
  double *ptrd;

  unsigned char *newuc;
  char *newc;
  unsigned short *newus;
  short *news;
  unsigned int*newui;
  int *newi;
  unsigned long *newul;
  long *newl;
  float *newf;
  double *newd;
  */

  if (type==field->datatype)
    return 0;

  if (sizeof_NDfield(type) > sizeof_NDfield(field->datatype))
    new_tab = calloc(field->nval,sizeof_NDfield(type));
  else
    new_tab = field->val;

  switch (field->datatype)
    {
    case ND_UCHAR: ptr_uc = (unsigned char*)field->val;break;
    case ND_CHAR: ptr_c = (char*)field->val;break;
    case ND_USHORT: ptr_us = (unsigned short*)field->val;break;
    case ND_SHORT: ptr_s = (short*)field->val;break;
    case ND_UINT: ptr_ui = (unsigned int*)field->val;break;
    case ND_INT: ptr_i = (int*)field->val;break;
    case ND_ULONG:ptr_ul = (unsigned long*)field->val;break;
    case ND_LONG: ptr_l = (long*)field->val;break;
    case ND_FLOAT: ptr_f = (float*)field->val;break;
    case ND_DOUBLE: ptr_d = (double*)field->val;break;
    }

  switch (type)
    {
    case ND_UCHAR: new_uc = (unsigned char*)new_tab;break;
    case ND_CHAR: new_c = (char*)new_tab;break;
    case ND_USHORT: new_us = (unsigned short*)new_tab;break;
    case ND_SHORT: new_s = (short*)new_tab;break;
    case ND_UINT: new_ui = (unsigned int*)new_tab;break;
    case ND_INT: new_i = (int*)new_tab;break;
    case ND_ULONG:new_ul = (unsigned long*)new_tab;break;
    case ND_LONG: new_l = (long*)new_tab;break;
    case ND_FLOAT: new_f = (float *)new_tab;break;
    case ND_DOUBLE: new_d = (double*)new_tab;break;
    }

  switch (field->datatype)
    {
    case ND_UCHAR:
	switch (type)
	  {
	  case ND_UCHAR: for (i=0;i<field->nval;i++) new_uc[i]=(unsigned char)ptr_uc[i]; break; 
	  case ND_CHAR: for (i=0;i<field->nval;i++) new_c[i]=(char)ptr_uc[i]; break; 
	  case ND_USHORT: for (i=0;i<field->nval;i++) new_us[i]=(unsigned short)ptr_uc[i]; break; 
	  case ND_SHORT: for (i=0;i<field->nval;i++) new_s[i]=(short)ptr_uc[i]; break; 
	  case ND_UINT: for (i=0;i<field->nval;i++) new_ui[i]=(unsigned int)ptr_uc[i]; break; 
	  case ND_INT: for (i=0;i<field->nval;i++) new_i[i]=(int)ptr_uc[i]; break; 
	  case ND_ULONG: for (i=0;i<field->nval;i++) new_ul[i]=(unsigned long)ptr_uc[i]; break; 
	  case ND_LONG: for (i=0;i<field->nval;i++) new_l[i]=(long)ptr_uc[i]; break; 
	  case ND_FLOAT: for (i=0;i<field->nval;i++) new_f[i]=(float)ptr_uc[i]; break; 
	  case ND_DOUBLE: for (i=0;i<field->nval;i++) new_d[i]=(double)ptr_uc[i]; break; 
	  }
      break;
    case ND_CHAR:
	switch (type)
	  {
	  case ND_UCHAR: for (i=0;i<field->nval;i++) new_uc[i]=(unsigned char)ptr_c[i]; break; 
	  case ND_CHAR: for (i=0;i<field->nval;i++) new_c[i]=(char)ptr_c[i]; break; 
	  case ND_USHORT: for (i=0;i<field->nval;i++) new_us[i]=(unsigned short)ptr_c[i]; break; 
	  case ND_SHORT: for (i=0;i<field->nval;i++) new_s[i]=(short)ptr_c[i]; break; 
	  case ND_UINT: for (i=0;i<field->nval;i++) new_ui[i]=(unsigned int)ptr_c[i]; break; 
	  case ND_INT: for (i=0;i<field->nval;i++) new_i[i]=(int)ptr_c[i]; break; 
	  case ND_ULONG: for (i=0;i<field->nval;i++) new_ul[i]=(unsigned long)ptr_c[i]; break; 
	  case ND_LONG: for (i=0;i<field->nval;i++) new_l[i]=(long)ptr_c[i]; break; 
	  case ND_FLOAT: for (i=0;i<field->nval;i++) new_f[i]=(float)ptr_c[i]; break; 
	  case ND_DOUBLE: for (i=0;i<field->nval;i++) new_d[i]=(double)ptr_c[i]; break; 
	  }
      break;
    case ND_USHORT:
	switch (type)
	  {
	  case ND_UCHAR: for (i=0;i<field->nval;i++) new_uc[i]=(unsigned char)ptr_us[i]; break; 
	  case ND_CHAR: for (i=0;i<field->nval;i++) new_c[i]=(char)ptr_us[i]; break; 
	  case ND_USHORT: for (i=0;i<field->nval;i++) new_us[i]=(unsigned short)ptr_us[i]; break; 
	  case ND_SHORT: for (i=0;i<field->nval;i++) new_s[i]=(short)ptr_us[i]; break; 
	  case ND_UINT: for (i=0;i<field->nval;i++) new_ui[i]=(unsigned int)ptr_us[i]; break; 
	  case ND_INT: for (i=0;i<field->nval;i++) new_i[i]=(int)ptr_us[i]; break; 
	  case ND_ULONG: for (i=0;i<field->nval;i++) new_ul[i]=(unsigned long)ptr_us[i]; break; 
	  case ND_LONG: for (i=0;i<field->nval;i++) new_l[i]=(long)ptr_us[i]; break; 
	  case ND_FLOAT: for (i=0;i<field->nval;i++) new_f[i]=(float)ptr_us[i]; break; 
	  case ND_DOUBLE: for (i=0;i<field->nval;i++) new_d[i]=(double)ptr_us[i]; break; 
	  }
      break;
    case ND_SHORT:
	switch (type)
	  {
	  case ND_UCHAR: for (i=0;i<field->nval;i++) new_uc[i]=(unsigned char)ptr_s[i]; break; 
	  case ND_CHAR: for (i=0;i<field->nval;i++) new_c[i]=(char)ptr_s[i]; break; 
	  case ND_USHORT: for (i=0;i<field->nval;i++) new_us[i]=(unsigned short)ptr_s[i]; break; 
	  case ND_SHORT: for (i=0;i<field->nval;i++) new_s[i]=(short)ptr_s[i]; break; 
	  case ND_UINT: for (i=0;i<field->nval;i++) new_ui[i]=(unsigned int)ptr_s[i]; break; 
	  case ND_INT: for (i=0;i<field->nval;i++) new_i[i]=(int)ptr_s[i]; break; 
	  case ND_ULONG: for (i=0;i<field->nval;i++) new_ul[i]=(unsigned long)ptr_s[i]; break; 
	  case ND_LONG: for (i=0;i<field->nval;i++) new_l[i]=(long)ptr_s[i]; break; 
	  case ND_FLOAT: for (i=0;i<field->nval;i++) new_f[i]=(float)ptr_s[i]; break; 
	  case ND_DOUBLE: for (i=0;i<field->nval;i++) new_d[i]=(double)ptr_s[i]; break; 
	  }
      break;
    case ND_UINT:
	switch (type)
	  {
	  case ND_UCHAR: for (i=0;i<field->nval;i++) new_uc[i]=(unsigned char)ptr_ui[i]; break; 
	  case ND_CHAR: for (i=0;i<field->nval;i++) new_c[i]=(char)ptr_ui[i]; break; 
	  case ND_USHORT: for (i=0;i<field->nval;i++) new_us[i]=(unsigned short)ptr_ui[i]; break; 
	  case ND_SHORT: for (i=0;i<field->nval;i++) new_s[i]=(short)ptr_ui[i]; break; 
	  case ND_UINT: for (i=0;i<field->nval;i++) new_ui[i]=(unsigned int)ptr_ui[i]; break; 
	  case ND_INT: for (i=0;i<field->nval;i++) new_i[i]=(int)ptr_ui[i]; break; 
	  case ND_ULONG: for (i=0;i<field->nval;i++) new_ul[i]=(unsigned long)ptr_ui[i]; break; 
	  case ND_LONG: for (i=0;i<field->nval;i++) new_l[i]=(long)ptr_ui[i]; break; 
	  case ND_FLOAT: for (i=0;i<field->nval;i++) new_f[i]=(float)ptr_ui[i]; break; 
	  case ND_DOUBLE: for (i=0;i<field->nval;i++) new_d[i]=(double)ptr_ui[i]; break; 
	  }
      break;
    case ND_INT:
	switch (type)
	  {
	  case ND_UCHAR: for (i=0;i<field->nval;i++) new_uc[i]=(unsigned char)ptr_i[i]; break; 
	  case ND_CHAR: for (i=0;i<field->nval;i++) new_c[i]=(char)ptr_i[i]; break; 
	  case ND_USHORT: for (i=0;i<field->nval;i++) new_us[i]=(unsigned short)ptr_i[i]; break; 
	  case ND_SHORT: for (i=0;i<field->nval;i++) new_s[i]=(short)ptr_i[i]; break; 
	  case ND_UINT: for (i=0;i<field->nval;i++) new_ui[i]=(unsigned int)ptr_i[i]; break; 
	  case ND_INT: for (i=0;i<field->nval;i++) new_i[i]=(int)ptr_i[i]; break; 
	  case ND_ULONG: for (i=0;i<field->nval;i++) new_ul[i]=(unsigned long)ptr_i[i]; break; 
	  case ND_LONG: for (i=0;i<field->nval;i++) new_l[i]=(long)ptr_i[i]; break; 
	  case ND_FLOAT: for (i=0;i<field->nval;i++) new_f[i]=(float)ptr_i[i]; break; 
	  case ND_DOUBLE: for (i=0;i<field->nval;i++) new_d[i]=(double)ptr_i[i]; break; 
	  }
      break;
    case ND_ULONG:
	switch (type)
	  {
	  case ND_UCHAR: for (i=0;i<field->nval;i++) new_uc[i]=(unsigned char)ptr_ul[i]; break; 
	  case ND_CHAR: for (i=0;i<field->nval;i++) new_c[i]=(char)ptr_ul[i]; break; 
	  case ND_USHORT: for (i=0;i<field->nval;i++) new_us[i]=(unsigned short)ptr_ul[i]; break; 
	  case ND_SHORT: for (i=0;i<field->nval;i++) new_s[i]=(short)ptr_ul[i]; break; 
	  case ND_UINT: for (i=0;i<field->nval;i++) new_ui[i]=(unsigned int)ptr_ul[i]; break; 
	  case ND_INT: for (i=0;i<field->nval;i++) new_i[i]=(int)ptr_ul[i]; break; 
	  case ND_ULONG: for (i=0;i<field->nval;i++) new_ul[i]=(unsigned long)ptr_ul[i]; break; 
	  case ND_LONG: for (i=0;i<field->nval;i++) new_l[i]=(long)ptr_ul[i]; break; 
	  case ND_FLOAT: for (i=0;i<field->nval;i++) new_f[i]=(float)ptr_ul[i]; break; 
	  case ND_DOUBLE: for (i=0;i<field->nval;i++) new_d[i]=(double)ptr_ul[i]; break; 
	  }
      break;
    case ND_LONG:
	switch (type)
	  {
	  case ND_UCHAR: for (i=0;i<field->nval;i++) new_uc[i]=(unsigned char)ptr_l[i]; break; 
	  case ND_CHAR: for (i=0;i<field->nval;i++) new_c[i]=(char)ptr_l[i]; break; 
	  case ND_USHORT: for (i=0;i<field->nval;i++) new_us[i]=(unsigned short)ptr_l[i]; break; 
	  case ND_SHORT: for (i=0;i<field->nval;i++) new_s[i]=(short)ptr_l[i]; break; 
	  case ND_UINT: for (i=0;i<field->nval;i++) new_ui[i]=(unsigned int)ptr_l[i]; break; 
	  case ND_INT: for (i=0;i<field->nval;i++) new_i[i]=(int)ptr_l[i]; break; 
	  case ND_ULONG: for (i=0;i<field->nval;i++) new_ul[i]=(unsigned long)ptr_l[i]; break; 
	  case ND_LONG: for (i=0;i<field->nval;i++) new_l[i]=(long)ptr_l[i]; break; 
	  case ND_FLOAT: for (i=0;i<field->nval;i++) new_f[i]=(float)ptr_l[i]; break; 
	  case ND_DOUBLE: for (i=0;i<field->nval;i++) new_d[i]=(double)ptr_l[i]; break; 
	  }
      break;
    case ND_FLOAT:
	switch (type)
	  {
	  case ND_UCHAR: for (i=0;i<field->nval;i++) new_uc[i]=(unsigned char)ptr_f[i]; break; 
	  case ND_CHAR: for (i=0;i<field->nval;i++) new_c[i]=(char)ptr_f[i]; break; 
	  case ND_USHORT: for (i=0;i<field->nval;i++) new_us[i]=(unsigned short)ptr_f[i]; break; 
	  case ND_SHORT: for (i=0;i<field->nval;i++) new_s[i]=(short)ptr_f[i]; break; 
	  case ND_UINT: for (i=0;i<field->nval;i++) new_ui[i]=(unsigned int)ptr_f[i]; break; 
	  case ND_INT: for (i=0;i<field->nval;i++) new_i[i]=(int)ptr_f[i]; break; 
	  case ND_ULONG: for (i=0;i<field->nval;i++) new_ul[i]=(unsigned long)ptr_f[i]; break; 
	  case ND_LONG: for (i=0;i<field->nval;i++) new_l[i]=(long)ptr_f[i]; break; 
	  case ND_FLOAT: for (i=0;i<field->nval;i++) new_f[i]=(float)ptr_f[i]; break; 
	  case ND_DOUBLE: for (i=0;i<field->nval;i++) new_d[i]=(double)ptr_f[i]; break; 
	  }
	break;
    case ND_DOUBLE:
	switch (type)
	  {
	  case ND_UCHAR: for (i=0;i<field->nval;i++) new_uc[i]=(unsigned char)ptr_d[i]; break; 
	  case ND_CHAR: for (i=0;i<field->nval;i++) new_c[i]=(char)ptr_d[i]; break; 
	  case ND_USHORT: for (i=0;i<field->nval;i++) new_us[i]=(unsigned short)ptr_d[i]; break; 
	  case ND_SHORT: for (i=0;i<field->nval;i++) new_s[i]=(short)ptr_d[i]; break; 
	  case ND_UINT: for (i=0;i<field->nval;i++) new_ui[i]=(unsigned int)ptr_d[i]; break; 
	  case ND_INT: for (i=0;i<field->nval;i++) new_i[i]=(int)ptr_d[i]; break; 
	  case ND_ULONG: for (i=0;i<field->nval;i++) new_ul[i]=(unsigned long)ptr_d[i]; break; 
	  case ND_LONG: for (i=0;i<field->nval;i++) new_l[i]=(long)ptr_d[i]; break; 
	  case ND_FLOAT: for (i=0;i<field->nval;i++) new_f[i]=(float)ptr_d[i]; break; 
	  case ND_DOUBLE: for (i=0;i<field->nval;i++) new_d[i]=(double)ptr_d[i]; break; 
	  }
      break;
    }
  
  

  if (sizeof_NDfield(type) > sizeof_NDfield(field->datatype))
    {
      free(field->val);
      field->val = new_tab;
    }
  else if (sizeof_NDfield(type) < sizeof_NDfield(field->datatype))
    field->val = realloc(field->val,field->nval*sizeof_NDfield(type));

  field->datatype=type;

  return 0;
}

// delta_bbox and x0_bbox are the coordinates of the GLOBAL bounding box
// This can be NULL if periodic=0
// result stored in x0_cut and x0_delta
int NDIntersection(double *x0_a,double *delta_a,double *x0_b,double *delta_b, double *x0_bbox,double *delta_bbox, double *x0_cut,double *delta_cut, int ndims, int periodic)
{
  int i,j;
  double xmin[ndims];
  double xmax[ndims];
  double x[ndims];
  int tmpa[ndims];
  int tmpb[ndims];

  for (j=0;j<ndims;j++) { xmin[j]=1.E40; xmax[j]=-1.E40;}

  //printf ("a:[%f,%f][%f,%f]\n",x0_a[0],x0_a[1],x0_a[0]+delta_a[0],x0_a[1]+delta_a[1]);
  //printf ("b:[%f,%f][%f,%f]\n",x0_b[0],x0_b[1],x0_b[0]+delta_b[0],x0_b[1]+delta_b[1]);
  
  /*
  if (periodic)
    for (j=0;j<ndims;j++)
      {
	if (periodic&(1<<j))
	  {
	    tmpa[j]=0;tmpb[j]=0;
	    if (x0_a[j]<x0_bbox[j])
	      tmpa[j] = (int)(((x0_bbox[j]+delta_bbox[j])-x0_a[j])/delta_bbox[j]);
	    else if (x0_a[j]>x0_bbox[j]+delta_bbox[j])
	      tmpa[j] = -(int)((x0_a[j]-x0_bbox[j])/delta_bbox[j]);

	    if (x0_b[j]<x0_bbox[j])
	      tmpb[j] = (int)(((x0_bbox[j]+delta_bbox[j])-x0_b[j])/delta_bbox[j]);
	    else if (x0_b[j]>x0_bbox[j]+delta_bbox[j])
	      tmpb[j] = -(int)((x0_b[j]-x0_bbox[j])/delta_bbox[j]);

	    if (tmpa[j]) x0_a[j]+=tmpa[j]*delta_bbox[j];
	    if (tmpb[j]) x0_b[j]+=tmpb[j]*delta_bbox[j];
	  }
      }
  */
  if (periodic)
    for (j=0;j<ndims;j++)
      {
	if (periodic&(1<<j))
	  {
	    tmpa[j]=0;tmpb[j]=0;
	    if (x0_a[j]<x0_bbox[j])
	      tmpa[j] = (int)(((x0_bbox[j]+delta_bbox[j])-x0_a[j])/delta_bbox[j]);
	    else if (x0_a[j]>x0_bbox[j]+delta_bbox[j])
	      tmpa[j] = -(int)((x0_a[j]-x0_bbox[j])/delta_bbox[j]);

	    if (tmpa[j]) x0_a[j]+=tmpa[j]*delta_bbox[j];
	    //printf ("%f %f, %f %d\n",x0_a[j],x0_b[j],(x0_a[j]-x0_b[j]),(int)(x0_a[j]-x0_b[j]));

	    if (x0_b[j]+delta_b[j]<x0_a[j]) 
	      tmpb[j] = (int)((x0_a[j]-x0_b[j]-delta_bbox[j])/delta_bbox[j]) +1;
	      
	    else if (x0_a[j]+delta_a[j]<x0_b[j])
	      tmpb[j] = (int)((x0_a[j]-x0_b[j]+delta_bbox[j])/delta_bbox[j]) -1;
	      
	  
	     if (tmpb[j]) x0_b[j]+=tmpb[j]*delta_bbox[j];
	  }
      }
  //printf ("a2:[%f,%f][%f,%f]\n",x0_a[0],x0_a[1],x0_a[0]+delta_a[0],x0_a[1]+delta_a[1]);
  //printf ("b2:[%f,%f][%f,%f]\n",x0_b[0],x0_b[1],x0_b[0]+delta_b[0],x0_b[1]+delta_b[1]);
  
  for (i=0;i<(1<<ndims);i++)
    for (j=0;j<ndims;j++)
      {
	x[j]=x0_a[j];
	if (i&(1<<j)) x[j] += delta_a[j]; 
	
	if ((x[j]>=x0_b[j])&&(x[j]<=x0_b[j]+delta_b[j]))
	  {
	    if (x[j]<xmin[j]) xmin[j]=x[j];
	    if (x[j]>xmax[j]) xmax[j]=x[j];
	  }
	
	x[j]=x0_b[j];
	if (i&(1<<j)) x[j] += delta_b[j]; 
	
	if ((x[j]>=x0_a[j])&&(x[j]<=x0_a[j]+delta_a[j]))
	  {
	    if (x[j]<xmin[j]) xmin[j]=x[j];
	    if (x[j]>xmax[j]) xmax[j]=x[j];
	  }
      }

  if (periodic)
    for (j=0;j<ndims;j++)
      if (periodic&(1<<j))
	{
	  if (tmpa[j]) x0_a[j]-=tmpa[j]*delta_bbox[j];
	  if (tmpb[j]) x0_b[j]-=tmpb[j]*delta_bbox[j];
	  //printf ("%d %d\n",tmpa[j],tmpb[j]);
	}
 
  for (j=0;j<ndims;j++) if (xmin[j]>xmax[j]) return 0;
  for (j=0;j<ndims;j++) {x0_cut[j]=xmin[j];delta_cut[j]=xmax[j]-xmin[j];}
  //for (j=0;j<ndims;j++) if ((tmpa[j])||(tmpb[j])) printf ("PERIODICITY USEDDDD !!!!!\n");
  return 1;
}

/*
NDfield *ComputeNodeField(NDfield *field,int periodic)
{
  //int *dims=field->dims;
  //int ndims=field->ndims;

  INT coord[field->ndims];
  INT dm[field->ndims];
  INT dp[field->ndims];
  INT dx[field->ndims+1];
  INT delta;
  int i,j,k;
  //INT tmp[3];
  
  FLOAT *val=field->val;
  FLOAT *new_val;

  NDfield *new_field;

  new_field=Create_NDfield(field->dims,field->ndims,field->fdims_index,field->datatype,field->x0, field->delta, NULL, "Node value");
  new_val=(FLOAT *)new_field->val;

  memset(coord,0,sizeof(INT)*field->ndims);
  memset(new_val,0,field->nval*sizeof(FLOAT));
  coord[0]=-1;

  //dims=field->dims;
  //ndims=field->ndims;

  for (i=0;i<field->nval;i++)
    {
      for (j=0,coord[0]++;coord[j]>=field->dims[j];coord[j]=0,coord[j+1]++,j++) {}
      
      dx[0]=1;
      
      for (j=0;j<field->ndims;j++)
	{
	  dx[j+1]=dx[j]*field->dims[j];
	  
	  if (coord[j]+1>=field->dims[j]) dp[j]=dx[j]-dx[j+1];
	  else dp[j]=dx[j];
	  
	  if (coord[j]-1<0) dm[j]=dx[j+1]-dx[j];
	  else dm[j]=-dx[j];
	}
      
      for (j=0;j<(1<<field->ndims);j++)
	{
	  delta=i;
	  for (k=0;k<field->ndims;k++)
	    {
	      if (j&(1<<k)) delta += dm[k];
	    }
	  new_val[i]+=val[delta];
	}
      new_val[i]/=(1<<field->ndims);
    }

  return new_field;
}
*/



// interpolates NDfield field at pos pos
// result is returned in result
// pos is the position in grid units (!!! not x0 and delta units !!!) (C EST FAUX)
// keepdim<0 for regular interpolation
// else interpolate along every dimension except keepdim.
// in this case result should be allocated for field->dims[keepdim] elements.
// field_p should be an NDfield
int InterpolateND(void *field_p,FLOAT *result,FLOAT *pos_p, int keepdim, int periodic)
{
  INT i,j,k;
  NDfield *field=(NDfield *)field_p;

  int *dims = field->dims;
  int ndims = field->ndims;
  FLOAT du[ndims];
  FLOAT p_ddu[ndims];
  FLOAT pos[ndims];

  long dx[ndims+1];
  long dy[ndims+1];
  INT a[ndims];
  //FLOAT *val;

  long index;
  long oldindex;
  FLOAT fac;
  FLOAT tmp=0;
  int nval;

  int retval=0;

  //delta = 0.5 * field->delta[i]/dims[i];

  for (i=0;i<ndims;i++)
  {
      if (periodic&(1<<i))
	  pos[i]=(FLOAT)(pos_p[i]-field->x0[i])/field->delta[i]*(dims[i]) - 0.5;
      else
	  pos[i]=(FLOAT)(pos_p[i]-field->x0[i])/field->delta[i]*(dims[i]-1);
  }
  DEFINE_NDVARS(val);

  SETNDPOINTER(val,field->val,field->datatype);

  memcpy(p_ddu,pos,sizeof(FLOAT)*ndims);

  dx[0]=dy[0]=1;
  index=0;

  for (i=0;i<ndims;i++)
  {
      dx[i+1] = dx[i]*dims[i];
      dy[i+1] = dx[i+1];
      if (i==keepdim)
      {
	  a[i]=dx[i]=du[i]=0;
	  continue;
      }

      if (periodic&(1<<i))
      {
	  if (p_ddu[i]>=dims[i]) p_ddu[i]-= ((long)(p_ddu[i] / dims[i])) * dims[i];
	  if (p_ddu[i]<0) p_ddu[i]+=((long)(1. - p_ddu[i] / dims[i])) * dims[i];
	  if (p_ddu[i]==dims[i]) p_ddu[i]=0;

      }
      else
      {
	  if (p_ddu[i]>=(double)dims[i]-1)
	  {
	      if (p_ddu[i]>(double)dims[i]-1) retval|=(1<<i);
	      p_ddu[i]= dims[i]-1;
	  }
	  if (p_ddu[i]<0) {p_ddu[i]=0;retval|=(1<<i);}
      }

      a[i] = (INT) p_ddu[i];
      du[i] = p_ddu[i]-a[i];

      index += a[i]*dx[i];

      if (a[i]==(dims[i]-1))
      {
	  if (periodic&(1<<i)) dx[i]-=dx[i+1];
	  else dx[i]=0;

      }
  }
  //printf ("%d %d %d - %f %f %f, %d %d %d\n",a[0],a[1],a[2],du[0],du[1],du[2],dx[0],dx[1],dx[2]);
  //val = p_val;


  if (keepdim>=field->ndims)
  {
      fprintf (stderr,"Error in InterpolateND, keepdim must be lower than number of dims.\n");
      exit(0);
  }

  if (keepdim>=0) nval = field->dims[keepdim];
  else nval=1;
  oldindex=index;

  for (k=0;k<nval;k++)
  {
      result[k]=0;
      if ((k!=0)&&(keepdim>=0)) oldindex += dy[keepdim];

      for (j=0;j<(1<<ndims);j++)
      {
	if ((keepdim>=0)&&!(j&(1<<keepdim))) j+=(1<<keepdim);
	  for (i=0,index=0,fac=1.;i<ndims;i++)
	  {
	    if (i==keepdim) continue;

	    if (j&(1<<i))
	      {
		fac*=(1.-du[i]);
	      }
	    else
	      {
		fac*=du[i];
		index+=dx[i];
	      }
	  }
	  index+=oldindex;
	  /* debug here */
	  SETNDFIELD_VAL(tmp,fac*val,index,field->datatype);
	  result[k]+=tmp;
      }

  }

  return retval;

}


// inverse of Index2Coords
// if periodic!=0, assume periodic boundary conditions
// else truncate coords.
void Coords2IndexfND(float *coord_p, INT *index,double *x0, double *delta,int *dims, int ndims,int periodic)
{
  float coord[ndims];
  INT val;
  INT dec;
  int i;

  for (i=0;i<ndims;i++)
      //if (periodic&(1<<i))
	  coord[i]=(coord_p[i]-x0[i])/delta[i]*dims[i] - 0.5;
  //else
  //coord[i]=(coord_p[i]-x0[i])/delta[i]*(dims[i]-1);

  *index=0;
  dec=1;
  for (i=0;i<ndims;i++)
    {
      val=(INT)coord[i];

      if (periodic&(1<<i))
	{
	  if (val>=dims[i]) val = val%(dims[i]);
	  if (val==dims[i]) val=0;
	  //else if (val<0) val = (100*(dims[i]) + val)%(dims[i]);
	  else if (val<0) val = dims[i]-((-val)%dims[i]);
	}
      else
	{
	  if (val>=dims[i]) val = dims[i]-1;
	  if (val<0) val=0;
	}

      *index += val*dec;
      dec*=dims[i];
    }

  return;
}

int printNDfieldStat(NDfield *field, int dec)
{
  char dataT[50];
  long i;

  print_dataType(dataT,field->datatype);
  for (i=0;i<dec;i++) printf (" ");
  for (i=0;i<80;i++) if ((i=='\0')||(i=='\n')) break;
  if (i!=80) printf("comment: '%s'\n",field->comment);
  for (i=0;i<dec;i++) printf (" ");
  printf ("Field is a [%ld",(long)field->dims[0]);
  for (i=1;i<field->n_dims;i++) printf(",%ld",(long)field->dims[i]);
  printf ("] array of type '%s'.\n",(field->fdims_index)?"coords":"grid");
  for (i=0;i<dec;i++) printf (" ");
  printf("Datatype is '%s'.\n",dataT);

  for (i=0;i<dec;i++) printf (" ");
  printf("Bounding box: x0=[%g",field->x0[0]);
  for (i=1;i<field->ndims;i++) printf(",%g",field->x0[i]);
  printf("],\n");
  for (i=0;i<dec;i++) printf (" ");
  printf("              delta=[%g",field->delta[0]);
  for (i=1;i<field->ndims;i++) printf(",%g",field->delta[i]);
  printf("].\n");
  
  return 0;
}
