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
#include <string.h>
#include <math.h>
#include <assert.h>

#include "NDnetwork.h"
#include "myendian.h"
#include "sortID.h"
#include "simplex.h"
#include "mystring.h"

#include "global.h"

size_t fwrite_dummy(size_t size, size_t nmemb,FILE *stream)
{
  char buffer[FWRITE_DUMMY_BFR_SIZE];
  unsigned long bsize=FWRITE_DUMMY_BFR_SIZE/size;
  unsigned long nleft=nmemb;
  size_t ret=0;
    
  while(nleft>0)
    {
      unsigned long nw=(nleft<bsize)?nleft:bsize;
      ret+=fwrite(buffer,size,nw,stream);
      nleft-=nw;
    };

  return ret;
}

int fread_sw_ului(void *data,size_t sizeOut,size_t nb,size_t sizeIn,FILE *f,int swap)
{
  int ret;
  if (sizeIn==sizeOut)
    return fread_sw(data,sizeIn,nb,f,swap);
  else if (sizeIn<sizeOut)
    {
      long i;
      ret=fread_sw(data,sizeIn,nb,f,swap);
      unsigned long int *uo = (unsigned long int *)(data);
      unsigned int *ui = (unsigned int *)(data);
      
      for (i=nb;i>=0;i--)
	uo[i]=(unsigned long int) ui[i];
    }
  else
    {      
      if (nb==1)
	{
	  unsigned long int in;
	  unsigned int *uo = (unsigned int *)(data);
	  ret=fread_sw(&in,sizeIn,nb,f,swap);
	  (*uo)= (unsigned int) in;
	  return ret;
	}
      unsigned long i;
      data=realloc(data,sizeIn*nb);
      ret=fread_sw(data,sizeIn,nb,f,swap);
      unsigned long int *ui = (unsigned long int *)(data);
      unsigned int *uo = (unsigned int *)(data);
      //assert(0);
      for (i=0;i<nb;i++)
	uo[i]=(unsigned int) ui[i];
	
      data=realloc(data,sizeOut*nb);
      //assert(0);
    }
  return ret;
}

void* fread_sw_fd(void *data,size_t sizeOut,size_t nb,size_t sizeIn,FILE *f,int swap)
{
  //int ret;  
  if (sizeIn==sizeOut)
    {
      fread_sw(data,sizeIn,nb,f,swap);
      return data;
    }
  else if (sizeIn<sizeOut)
    {
      long i;
      fread_sw(data,sizeIn,nb,f,swap);
      double *uo = (double *)(data);
      float *ui = (float *)(data);
      
      for (i=nb;i>=0;i--)
	uo[i]=(double) ui[i];
    }
  else
    {      
      if (nb==1)
	{
	  double in;
	  float *uo = (float *)(data);
	  fread_sw(&in,sizeIn,nb,f,swap);
	  (*uo)= (float) in;
	  return data;
	}
      unsigned long i;
      data=realloc(data,sizeIn*nb);
      
      fread_sw(data,sizeIn,nb,f,swap);
      double *ui = (double *)(data);
      float *uo = (float *)(data);
      //assert(0);
      for (i=0;i<nb;i++)
	uo[i]=(float) ui[i];
	
      data=realloc(data,sizeOut*nb);
      
      //assert(0);
    } 
  return data;
}

int printNDnetStat(NDnetwork *net, int dec)
{
  long i;
  long m_dims=0;

  for (i=0;i<=net->ndims;i++) 
    if (net->haveVertexFromFace[i]) m_dims=i;

  for (i=0;i<dec;i++) printf (" ");
  for (i=0;i<80;i++) if ((i=='\0')||(i=='\n')) break;
  if (i!=80) printf("comment: '%s'\n",net->comment);
  for (i=0;i<dec;i++) printf (" ");
  printf ("%dD-Network has %ld vertices.\n",net->ndims,(long)net->nvertex);

  for (i=0;i<dec;i++) printf (" ");
  if (!net->periodicity)
    printf ("periodicity: non-periodic.\n");
  else
    {
      printf ("periodicity: ");
      for (i=0;i<net->ndims;i++)
	if (net->periodicity&(1<<i)) printf ("1");
	else printf ("0");
      printf (".\n");
    }

  for (i=0;i<dec;i++) printf (" ");  
  printf ("Available faces: ");
  for (i=0;i<=net->ndims;i++)
    {
      if (i!=m_dims)
	{
	  if (net->nfaces[i]) 
	    printf ("%ld %ld-F, ",(long)net->nfaces[i],i);
	}
      else
	{
	  if (net->nfaces[i]) 
	    printf ("%ld %ld-F.\n",(long)net->nfaces[i],i);
	}
    }

  for (i=0;i<dec;i++) printf (" ");
  printf("Bounding box: x0=[%g",net->x0[0]);
  for (i=1;i<net->ndims;i++) printf(",%g",net->x0[i]);
  printf("],\n");
  for (i=0;i<dec;i++) printf (" ");
  printf("              delta=[%g",net->delta[0]);
  for (i=1;i<net->ndims;i++) printf(",%g",net->delta[i]);
  printf("].\n");

  
  for (i=0;i<dec;i++) printf (" ");
  printf("Available fields: ");
  if (net->ndata) printf("'%s'(%d)",net->data[0].name,(int)net->data[0].type);
  for (i=1;i<net->ndata;i++)
    {
      if ((i%3) == 0)
	{
	  long j;
	  printf("\n");
	  for (j=0;j<dec+strlen("Available fields: ");j++) printf (" ");
	  printf("'%s'(%d)",net->data[i].name,(int)net->data[i].type);
	}
      else
	printf(", '%s'(%d)",net->data[i].name,(int)net->data[i].type);
    }
  if (net->ndata) printf("\n");

  return 0;
}

int IsNDnetwork(const char *filename)
{
  int i;
  char tag[NDNETWORK_DATA_STR_SIZE+2];
  FILE *f;
  int swap=0;
  
  memset(tag,0,16*sizeof(char));
   
  if(!(f = fopen(filename,"r")))
    {
	fprintf(stderr,"File %s does not exist.\n",filename);
	return 0;
    }

  fread_sw(&i,sizeof(int),1,f,swap);
  if (i!=NDNETWORK_DATA_STR_SIZE) swap=1-swap;
  fread_sw(tag,sizeof(char),NDNETWORK_DATA_STR_SIZE,f,swap);
  fread_sw(&i,sizeof(int),1,f,swap);
  
  fclose(f); 
  tag[NDNETWORK_DATA_STR_SIZE+1]='\0';

  if (strcmp(tag,NDNETWORK_TAG)) 
    {

    //if (isAscii(filename)) {
      f = fopen(filename,"r");
      char *line=NULL;int n;
      Mygetline (&line,&n,f);      
      fclose(f);
      
      line[strlen(line)-1]='\0';
      //printf("line(%d) = '%s'\n",n,line);
      if (!strcmp(line,NDNETWORK_ASCII_TAG)) 
	{
	  free(line);
	  return 2;
	}
      /*
      if (sscanf(line,"%d %d",&i,&n)==1) {
	free(line);
	
	return 2;
      }
      */
      free(line);
      return 0;
    }

  return 1;
}


void FreeNDnetwork(NDnetwork **net_p)
{
    NDnetwork *net=*net_p;
    int i,j;
    
    if (*net_p==NULL) return;
    free(net->x0);
    free(net->delta);
    free(net->nfaces);
    free(net->v_coord);
    for (i=0;i<=net->ndims;i++)
    {
	if (net->haveVertexFromFace[i])
	{
	  free(net->f_vertexIndex[i]);
	    if (!net->isSimpComplex)
	      free(net->f_numVertexIndexCum[i]);
	}
    }
    free(net->haveVertexFromFace);
    free(net->f_vertexIndex);
    free(net->f_numVertexIndexCum);
    
    for (i=0;i<=net->ndims;i++)
    {
	if (net->haveFaceFromVertex[i])
	{
	    free(net->v_faceIndex[i]);
	    free(net->v_numFaceIndexCum[i]);
	}
    }
    free(net->haveFaceFromVertex);
    free(net->v_faceIndex);
    free(net->v_numFaceIndexCum);

    for (i=0;i<=net->ndims;i++)
    {
	for (j=0;j<=net->ndims;j++)
	{
	    if (net->haveFaceFromFace[i][j])
	    {
		free(net->f_faceIndex[i][j]);
		free(net->f_numFaceIndexCum[i][j]);
	    }
	}
	free(net->haveFaceFromFace[i]);
	free(net->f_faceIndex[i]);
	free(net->f_numFaceIndexCum[i]);
    }
    free(net->haveFaceFromFace);
    free(net->f_faceIndex);
    free(net->f_numFaceIndexCum);

    if (net->haveVFlags) free(net->v_flag);
    for (i=0;i<=net->ndims;i++)
	if (net->haveFFlags[i]) free(net->f_flag[i]);
    free(net->f_flag);
    free(net->haveFFlags);

    for (i=0;i<net->ndata;i++)
    {
	free(net->data[i].data);
    }
    free(net->data);

    for (i=0;i<net->nsupData;i++)
    {
	free(net->supData[i].data);
    }
    free(net->supData);
    free(*net_p);
    *net_p=NULL;
}


NDnetwork *CreateNetwork(int ndims, int nvertex, int vflags)
{
    NDnetwork *net;
    long k;

    net = calloc(1,sizeof(NDnetwork));
    net->isSimpComplex=1;
    net->ndims=ndims;
    net->ndims_net=ndims;
    net->x0=malloc(net->ndims*sizeof(double));
    net->delta=malloc(net->ndims*sizeof(double));

    strcpy(net->comment,"");

    net->nvertex=nvertex;
    // FIX ME : should be NDNET_FLOAT not float
    if (nvertex!=0) 
   	net->v_coord=malloc(sizeof(float)*(size_t)net->ndims*(size_t)net->nvertex);
    else net->v_coord=NULL;

    net->nfaces=calloc(((size_t)net->ndims+1),sizeof(NDNET_UINT));
    
    net->haveVertexFromFace=calloc(net->ndims+1,sizeof(int));
    net->f_vertexIndex = calloc((1+net->ndims),sizeof(NDNET_UINT*));
    net->f_numVertexIndexCum = calloc((1+net->ndims),sizeof(NDNET_IDCUMT*));

    net->haveFaceFromVertex=calloc(net->ndims+1,sizeof(int));
    net->v_faceIndex = calloc((1+net->ndims),sizeof(NDNET_UINT*));
    net->v_numFaceIndexCum = calloc((1+net->ndims),sizeof(NDNET_IDCUMT*));

    net->haveFaceFromFace =calloc(net->ndims+1,sizeof(int*));
    net->f_faceIndex = calloc(net->ndims+1,sizeof(NDNET_UINT**));
    net->f_numFaceIndexCum = calloc(net->ndims+1,sizeof(NDNET_IDCUMT**));

    for (k=0;k<net->ndims+1;k++)
    {
      net->haveFaceFromFace[k]=calloc(((size_t)net->ndims+1),sizeof(int));
      net->f_faceIndex[k] = calloc((1+net->ndims),sizeof(NDNET_UINT*));
      net->f_numFaceIndexCum[k] = calloc(1+net->ndims,sizeof(NDNET_IDCUMT*));
    }
    
    if (nvertex) 
	net->haveVFlags=vflags;
    else
	net->haveVFlags=0;

    if (vflags)
    {
	if (nvertex)
	    net->v_flag=calloc(net->nvertex,sizeof(unsigned char));
	else
	    net->v_flag=NULL;
    }
    
    net->haveFFlags=calloc(net->ndims+1,sizeof(int));
    net->f_flag=calloc(net->ndims+1,sizeof(unsigned char *));

    net->ndata=net->nsupData=0;
    net->data=NULL;
    net->supData=NULL;

    net->indexSize=sizeof(NDNET_UINT);
    net->cumIndexSize=sizeof(NDNET_IDCUMT);
    net->floatSize=sizeof(NDNET_FLOAT);
    //printf("net : %d %d\n",net->indexSize,net->cumIndexSize);
    return net;
}

int Save_NDnetwork(NDnetwork *net,const char *filename)
{
  int j;
  long i,k;
  char tag[NDNETWORK_DATA_STR_SIZE];
  FILE *f;
  
  if (verbose>1) printf ("Saving %dD network to file %s ...",net->ndims,filename);fflush(0);

  memset(tag,0,NDNETWORK_DATA_STR_SIZE*sizeof(char));
  strcpy(tag,NDNETWORK_TAG);
  i=NDNETWORK_DATA_STR_SIZE;
  
  f=fopen(filename,"w");
  fwrite(&i,sizeof(int),1,f);
  fwrite(tag,sizeof(char),NDNETWORK_DATA_STR_SIZE,f);
  fwrite(&i,sizeof(int),1,f);

  j=sizeof(int)*2;
  fwrite(&j,sizeof(int),1,f);
  fwrite(&net->ndims,sizeof(int),1,f);
  fwrite(&net->ndims_net,sizeof(int),1,f);
  fwrite(&j,sizeof(int),1,f);  
  j=sizeof(NDNET_UINT)+sizeof(int)*4+80*sizeof(char)+net->ndims*sizeof(double)*2+152*sizeof(char);
  fwrite(&j,sizeof(int),1,f);
  fwrite(net->comment,sizeof(char),80,f);
  fwrite(&net->periodicity,sizeof(int),1,f);
  fwrite(&net->isSimpComplex,sizeof(int),1,f);
  fwrite(net->x0,sizeof(double),net->ndims,f);
  fwrite(net->delta,sizeof(double),net->ndims,f);
  fwrite(&net->indexSize,sizeof(int),1,f);
  fwrite(&net->cumIndexSize,sizeof(int),1,f);
  fwrite(net->dummy,sizeof(char),160-8,f);
  fwrite(&net->nvertex,sizeof(NDNET_UINT),1,f);
  fwrite(&j,sizeof(int),1,f);
  
  j=net->ndims*net->nvertex;
  fwrite(&j,sizeof(int),1,f);
  fwrite(net->v_coord,sizeof(float),(size_t)net->ndims*(size_t)net->nvertex,f);
  fwrite(&j,sizeof(int),1,f);
  
  j=(1+net->ndims)*sizeof(NDNET_UINT);
  fwrite(&j,sizeof(int),1,f);
  fwrite(net->nfaces,sizeof(NDNET_UINT),((size_t)net->ndims+1),f);
  fwrite(&j,sizeof(int),1,f);

  j=(1+net->ndims)*sizeof(int);
  fwrite(&j,sizeof(int),1,f);
  fwrite(net->haveVertexFromFace,sizeof(int),((size_t)net->ndims+1),f);
  fwrite(&j,sizeof(int),1,f);

  
  for (i=0;i<1+net->ndims;i++)
  {
      
      if (net->haveVertexFromFace[i])
      {
	  if (!net->isSimpComplex)
	  {
	      j=sizeof(NDNET_IDCUMT)*((size_t)net->nfaces[i]+1);
	      fwrite(&j,sizeof(int),1,f);
	      fwrite(net->f_numVertexIndexCum[i],sizeof(NDNET_IDCUMT),((size_t)net->nfaces[i]+1),f);
	      fwrite(&j,sizeof(int),1,f);

	      j=sizeof(NDNET_UINT)*((size_t)net->f_numVertexIndexCum[i][net->nfaces[i]]);
	      fwrite(&j,sizeof(int),1,f);
	      fwrite(net->f_vertexIndex[i],sizeof(NDNET_UINT),((size_t)net->f_numVertexIndexCum[i][net->nfaces[i]]),f);
	      fwrite(&j,sizeof(int),1,f);
	  }
	  else
	  {
	      j=sizeof(NDNET_UINT)*((size_t)(i+1)*net->nfaces[i]);
	      fwrite(&j,sizeof(int),1,f);
	      fwrite(net->f_vertexIndex[i],sizeof(NDNET_UINT),((size_t)(i+1)*net->nfaces[i]),f);
	      fwrite(&j,sizeof(int),1,f);
	  }
	  
      }
  }
  
  j=(1+net->ndims)*sizeof(int);
  fwrite(&j,sizeof(int),1,f);
  fwrite(net->haveFaceFromVertex,sizeof(int),((size_t)net->ndims+1),f);
  fwrite(&j,sizeof(int),1,f);

  for (i=0;i<1+net->ndims;i++)
  {
      
      if (net->haveFaceFromVertex[i])
      {
	  j=sizeof(NDNET_IDCUMT)*((size_t)net->nvertex+1);
	  fwrite(&j,sizeof(int),1,f);
	  fwrite(net->v_numFaceIndexCum[i],sizeof(NDNET_IDCUMT),((size_t)net->nvertex+1),f);
	  fwrite(&j,sizeof(int),1,f);
	  
	  j=sizeof(NDNET_UINT)*((size_t)net->v_numFaceIndexCum[i][net->nvertex]);
	  fwrite(&j,sizeof(int),1,f);
	  fwrite(net->v_faceIndex[i],sizeof(NDNET_UINT),((size_t)net->v_numFaceIndexCum[i][net->nvertex]),f);
	  fwrite(&j,sizeof(int),1,f);	  
      }
  }
  
  j=(1+net->ndims)*(1+net->ndims)*sizeof(int);
  fwrite(&j,sizeof(int),1,f);
  for (k=0;k<net->ndims+1;k++)
    fwrite(net->haveFaceFromFace[k],sizeof(int),((size_t)net->ndims+1),f);
  fwrite(&j,sizeof(int),1,f);

  for (i=0;i<1+net->ndims;i++)
  {
      for (k=0;k<1+net->ndims;k++)
      {

	  if (net->haveFaceFromFace[i][k])
	  {
	      j=sizeof(NDNET_IDCUMT)*(1+(size_t)net->nfaces[i]);
	      fwrite(&j,sizeof(int),1,f);
	      fwrite(net->f_numFaceIndexCum[i][k],sizeof(NDNET_IDCUMT),((size_t)net->nfaces[i]+1),f);
	      fwrite(&j,sizeof(int),1,f);
	      
	      j=sizeof(NDNET_UINT)*((size_t)net->f_numFaceIndexCum[i][k][net->nfaces[i]]);
	      fwrite(&j,sizeof(int),1,f);
	      fwrite(net->f_faceIndex[i][k],sizeof(NDNET_UINT),((size_t)net->f_numFaceIndexCum[i][k][net->nfaces[i]]),f);
	      fwrite(&j,sizeof(int),1,f);	  
	  }
      }
  }

  j=sizeof(int);
  fwrite(&j,sizeof(int),1,f);
  fwrite(&net->haveVFlags,sizeof(int),1,f);
  fwrite(&j,sizeof(int),1,f);

  if (net->haveVFlags)
  {
      j=sizeof(unsigned char)*net->nvertex;
      fwrite(&j,sizeof(int),1,f);
      fwrite(net->v_flag,sizeof(unsigned char),(size_t)net->nvertex,f);
      fwrite(&j,sizeof(int),1,f);
  }

  
  j=sizeof(int)*(net->ndims+1);
  fwrite(&j,sizeof(int),1,f);
  fwrite(net->haveFFlags,sizeof(int),(net->ndims+1),f);
  fwrite(&j,sizeof(int),1,f);
  
  for (i=0;i<1+net->ndims;i++)
      if (net->haveFFlags[i])
      {
	  j=sizeof(unsigned char)*net->nfaces[i];
	  fwrite(&j,sizeof(int),1,f);
	  if (j) fwrite(net->f_flag[i],sizeof(unsigned char),(size_t)net->nfaces[i],f);
	  fwrite(&j,sizeof(int),1,f);
      }

  j=sizeof(int);
  fwrite(&j,sizeof(int),1,f);
  fwrite(&net->ndata,sizeof(int),1,f);
  fwrite(&j,sizeof(int),1,f);
  
  for (i=0;i<net->ndata;i++)
  {
      
      j=sizeof(int)+255*sizeof(char);
      fwrite(&j,sizeof(int),1,f);
      fwrite(&(net->data[i].type),sizeof(int),1,f);
      fwrite(net->data[i].name,sizeof(char)*255,1,f);
      fwrite(&j,sizeof(int),1,f);

      if (net->data[i].type==0)
      {
	j=sizeof(double)*net->nvertex;
	fwrite(&j,sizeof(int),1,f);
	fwrite(net->data[i].data,sizeof(double),(size_t)net->nvertex,f);
	fwrite(&j,sizeof(int),1,f);
      }
      else
      {
	j=sizeof(double)*net->nfaces[net->data[i].type];
	fwrite(&j,sizeof(int),1,f);
	fwrite(net->data[i].data,sizeof(double),(size_t)net->nfaces[net->data[i].type],f);
	fwrite(&j,sizeof(int),1,f);
      }  
  }

  j=sizeof(int);
  fwrite(&j,sizeof(int),1,f);
  fwrite(&net->nsupData,sizeof(int),1,f);
  fwrite(&j,sizeof(int),1,f);
  
  for (i=0;i<net->nsupData;i++)
  {
      
      j=2*sizeof(int)+2*255*sizeof(char);
      fwrite(&j,sizeof(int),1,f);
      fwrite(&(net->supData[i].type),sizeof(int),1,f);
      fwrite(net->supData[i].name,sizeof(char)*255,1,f);
      fwrite(&(net->supData[i].datasize),sizeof(int),1,f);
      fwrite(net->supData[i].datatype,sizeof(char)*255,1,f);
      fwrite(&j,sizeof(int),1,f);

      if (net->supData[i].type==0)
      {
	  j=(size_t)net->supData[i].datasize*net->nvertex;
	  fwrite(&j,sizeof(int),1,f);
	  fwrite(net->supData[i].data,(size_t)net->supData[i].datasize,(size_t)net->nvertex,f);
	  fwrite(&j,sizeof(int),1,f);
      }
      else
      {
	  j=(size_t)net->supData[i].datasize*net->nfaces[net->supData[i].type];
	  fwrite(&j,sizeof(int),1,f);
	  fwrite(net->supData[i].data,(size_t)net->supData[i].datasize,(size_t)net->nfaces[net->supData[i].type],f);
	  fwrite(&j,sizeof(int),1,f);
      }  
  }

 

  fclose(f);
  if (verbose>1) printf (" done.\n");
  return 0;
}

int Save_NDnetwork_ASCII(NDnetwork *net, const char *filename)
{
  FILE *f;  
  long i,j,k;
 
  if (verbose>1) printf ("Saving %dD network to ASCII file %s ...",net->ndims,filename);fflush(0);
  f=fopen(filename,"w");
  if (f==NULL) return -1;

  fprintf(f,"%s\n",NDNETWORK_ASCII_TAG);
  fprintf(f,"%d\n",net->ndims);
  fprintf(f,"#%s\n",net->comment);
  fprintf(f,"BBOX [%g",net->x0[0]);
  for (i=1;i<net->ndims;i++) fprintf(f,",%g",net->x0[i]);
  fprintf(f,"] [%g",net->delta[0]);
  for (i=1;i<net->ndims;i++) fprintf(f,",%g",net->delta[i]);
  fprintf(f,"]\n");
  
  fprintf(f,"%ld\n",(long)net->nvertex);
  for (i=0;i<net->nvertex;i++)
    {
      fprintf(f,"%g",net->v_coord[i*net->ndims]);
      for (j=1;j<net->ndims;j++)
	fprintf(f," %g",net->v_coord[i*net->ndims+j]);
      fprintf(f,"\n");
    }

  for (i=0;i<=net->ndims;i++)
    {
      if (net->haveVertexFromFace[i])
	{
	  fprintf(f,"%ld %ld\n",i,(long)net->nfaces[i]);
	  for (j=0;j<net->nfaces[i];j++)
	    {
	      fprintf(f,"%ld",(long)net->f_vertexIndex[i][(i+1)*j]);
	      for (k=1;k<(i+1);k++)
		fprintf(f," %ld",(long)net->f_vertexIndex[i][(i+1)*j+k]);
	      fprintf(f,"\n");
	    }
	}
    }
  
  int dataPart=0;

  if (net->ndata)
    {
      fprintf(f,"%s\n",NDNETWORK_ASCII_ADDTAG);
      dataPart=1;
      for (i=0;i<net->ndata;i++)
	{
	  fprintf(f,"%s\n",net->data[i].name);
	  fprintf(f,"%d\n",(int)net->data[i].type);
	  long N=(net->data[i].type==0)?net->nvertex:net->nfaces[net->data[i].type];
	  for (j=0;j<N;j++)
	    fprintf(f,"%g\n",net->data[i].data[j]);
	}
    }

  if (net->haveVFlags)
    {
      if (!dataPart)
	{
	  fprintf(f,"%s\n",NDNETWORK_ASCII_ADDTAG);
	  dataPart=1;
	}
      fprintf(f,"%s\n",NDNETWORK_V_FLAG_TAG);
      fprintf(f,"%d\n",(int)0);
      for (j=0;j<net->nvertex;j++)
	fprintf(f,"%u\n",(unsigned int)net->v_flag[j]);
    }
  for (i=0;i<=net->ndims;i++)
    {
      if (net->haveFFlags[i])
	{
	  if (!dataPart)
	    {
	      fprintf(f,"%s\n",NDNETWORK_ASCII_ADDTAG);
	      dataPart=1;
	    }
	  char name[255];
	  sprintf(name,"%s_%d",NDNETWORK_F_FLAG_TAG,(int)i);
	  fprintf(f,"%s\n",name);
	  fprintf(f,"%d\n",(int)i);
	  for (j=0;j<net->nfaces[i];j++)
	    fprintf(f,"%d\n",net->f_flag[i][j]);
	}
    }
  fclose(f);

  if (verbose>1) printf("done.\n");

  return 0;
}

NDnetwork *Load_NDnetwork_ASCII(const char *filename)
{
  FILE *f;
  int ndims;
  long nv;
  long nt;
  float p[4];
  long i,j;
  double x0[10];
  double delta[10];
  //long n;
  int in;
  int dummy;
  char *line=NULL;
  char comment[80];
  f=fopen(filename,"r");
  Mygetline (&line,&in,f);
  line[strlen(line)-1]='\0';
  if (strcmp(line,NDNETWORK_ASCII_TAG))
    {
      fprintf(stderr,"ERROR trying to load file '%s'.\n",filename);
      fprintf(stderr,"      not a valid '%s' file.\n",NDNETWORK_ASCII_TAG);
      exit(0);
    }
  dummy=fscanf(f,"%d\n",&ndims); 
  
  dummy=Mygetline (&line,&in,f);
  if (line[0]=='#')
    {
      if (dummy>80) line[79]='\0';
      strcpy(comment,line);
      Mygetline (&line,&in,f);
    }
  //else strcpy(comment,"loaded from ASCII file");
  //printf("comment : '%s'. %d\n",comment,ndims);

  for (i=0;i<10;i++) x0[i]=delta[i]=-1;
  char *c=NULL;
  if ((c=strstr(line,"BBOX"))!=NULL)
    {
      char *str[1000];
      int ntok=str2tok(c,"[ ],;|",6,0,0,str);
      if (ntok != ndims*2+1)
	{
	  fprintf(stderr,"ERROR trying to load file '%s'.\n",filename);
	  fprintf(stderr,"  Bad definition of the bounding box, should have ndims=%d\n",ndims);
	  exit(-1);
	}
      for (i=1;i<ntok;i++)
	{
	  if (i<=ndims) x0[i-1]=strtod(str[i],NULL);
	  else delta[i-1-ndims]=strtod(str[i],NULL);
	}
      Mygetline(&line,&in,f);
    }
  
  dummy=sscanf(line,"%ld\n",&nv);  
  Mygetline(&line,&in,f);
  dummy=sscanf(line,"%g %g %g %g",&p[0],&p[1],&p[2],&p[3]); 
 
  if (dummy!=ndims)
    {
      fprintf (stderr,"ERROR in Load_NDnetwork_ASCII: I expect %d coordinates (%d) for each vertex.\n",ndims,dummy);
    exit(-1);
    }

  if ((ndims>3)||(ndims<1)) {
    free(line);
    fprintf (stderr,"ERROR in Load_NDnetwork_ASCII: Only 1D/2D/3D networks supported from ASCII files.\n");
    exit(-1);
  }
  NDnetwork *net = CreateNetwork(ndims,nv,0);
  strcpy(net->comment,comment);
  net->ndata=0;

  if (verbose>1) printf ("Loading %dD network from ASCII file \"%s\" ...",net->ndims,filename);fflush(0);
  
  for (i=0;i<ndims;i++) net->v_coord[i]=p[i];

  if (ndims ==1)
    for (i=1;i<nv;i++) 
      dummy=fscanf(f,"%g\n",&net->v_coord[i]);
  else if (ndims==2)
    for (i=1;i<nv;i++) 
      dummy=fscanf(f,"%g %g\n",&net->v_coord[i*2],&net->v_coord[i*2+1]);
  else if (ndims==3)
    for (i=1;i<nv;i++) 
      dummy=fscanf(f,"%g %g %g\n",&net->v_coord[i*3],&net->v_coord[i*3+1],&net->v_coord[i*3+2]);
  else if (ndims==4)
    for (i=1;i<nv;i++) 
      dummy=fscanf(f,"%g %g %g %g\n",&net->v_coord[i*4],&net->v_coord[i*4+1],&net->v_coord[i*4+2],&net->v_coord[i*4+3]);

  int addData=0;
  int ret=0;
  while ((ret=Mygetline (&line,&in,f))>=0) 
    {      
      if (ret==0) continue;
      //printf("%s : IN = %d\n",line,ret);

      if (!addData)
	{
	  char tmpL[1024];
	  strcpy(tmpL,line);
	  tmpL[strlen(tmpL)-1]='\0';
	  if (strstr(tmpL,NDNETWORK_ASCII_ADDTAG)!=NULL) 
	    {
	      addData=1;
	      continue;
	    }	  
	}

      if (!addData)
	{
	  int type;
	  dummy=sscanf(line,"%d %ld\n",&type,&nt);
	  net->f_vertexIndex[type] = (NDNET_UINT *) malloc(nt*(type+1)*sizeof(NDNET_UINT));	  
	  net->haveVertexFromFace[type] = 1;
	  net->nfaces[type]=nt;
	  if (sizeof(NDNET_UINT)==4)
	    {
	      unsigned int *vi = (unsigned int *)net->f_vertexIndex[type];
	      if (type ==1)
		for (i=0;i<nt;i++) 
		  dummy=fscanf(f,"%d %d\n",&vi[i*2],&vi[i*2+1]);
	      else if (type==2)
		for (i=0;i<nt;i++) 
		  dummy=fscanf(f,"%d %d %d\n",&vi[i*3],&vi[i*3+1],&vi[i*3+2]);
	      else if (type==3)
		for (i=0;i<nt;i++) 
		  dummy=fscanf(f,"%d %d %d %d\n",&vi[i*4],&vi[i*4+1],&vi[i*4+2],&vi[i*4+3]);
	      else if (type==4)
		for (i=0;i<nt;i++) 
		  dummy=fscanf(f,"%d %d %d %d %d\n",&vi[i*5],&vi[i*5+1],&vi[i*5+2],&vi[i*5+3],&vi[i*5+4]);
	    }
	  else
	    {
	      unsigned long int *vi = (unsigned long int *)net->f_vertexIndex[type];
	      if (type ==1)
		for (i=0;i<nt;i++) 
		  dummy=fscanf(f,"%ld %ld\n",&vi[i*2],&vi[i*2+1]);
	      else if (type==2)
		for (i=0;i<nt;i++) 
		  dummy=fscanf(f,"%ld %ld %ld\n",&vi[i*3],&vi[i*3+1],&vi[i*3+2]);
	      else if (type==3)
		for (i=0;i<nt;i++) 
		  dummy=fscanf(f,"%ld %ld %ld %ld\n",&vi[i*4],&vi[i*4+1],&vi[i*4+2],&vi[i*4+3]);
	      else if (type==4)
		for (i=0;i<nt;i++) 
		  dummy=fscanf(f,"%ld %ld %ld %ld %ld\n",&vi[i*5],&vi[i*5+1],&vi[i*5+2],&vi[i*5+3],&vi[i*5+4]);
	    }
	}
      else
	{
	  if (strstr(line,NDNETWORK_V_FLAG_TAG)!=NULL)
	    {
	      int dt=-1;
	      dummy=fscanf(f,"%d\n",&dt); //printf("t=%d\n",data->type);
	      if (dt!=0) 
		{
		  fprintf(stderr,"ERROR: field '%s' has type %d, should be 0.\n",NDNETWORK_V_FLAG_TAG,dt);
		  fprintf(stderr,"  info: this field is use to store internal flags.\n");
		  exit(0);
		}
	      net->v_flag=(unsigned char*)calloc(net->nvertex,sizeof(unsigned char));
	      net->haveVFlags=1;
	      for (i=0;i<net->nvertex;i++) 
		{
		  dummy=fscanf(f,"%d\n",&dt);
		  net->v_flag[i]=(unsigned char)dt;
		}
	      
	      continue;
	    }

	   if (strstr(line,NDNETWORK_F_FLAG_TAG)!=NULL)
	    {
	      int dt=-2;
	      int tp=-1;
	      tp = strstr(line,NDNETWORK_F_FLAG_TAG)[strlen(NDNETWORK_F_FLAG_TAG)+1] - '0';
	      line[strlen(line)-1]='\0';
	      dummy=fscanf(f,"%d\n",&dt); //printf("t=%d\n",data->type);
	      if (dt!=tp) 
		{
		  fprintf(stderr,"ERROR: field '%s' has type %d, should probably be %d.\n",line,dt,tp);
		  fprintf(stderr,"  info: this field is use to store internal flags.\n");
		  exit(0);
		}
	      net->f_flag[tp]=(unsigned char*)calloc(net->nfaces[tp],sizeof(unsigned char));
	      net->haveFFlags[tp]=1;;
	      for (i=0;i<net->nfaces[tp];i++) 
		{
		  dummy=fscanf(f,"%d\n",&dt);
		  net->f_flag[tp][i]=(unsigned char)dt;
		}
	      continue;
	    }

	  net->ndata++;
	  net->data = (NDnetwork_Data*) realloc(net->data,net->ndata*sizeof(NDnetwork_Data));
	  NDnetwork_Data *data=&net->data[net->ndata-1];
	  //Mygetline (&line,&in,f);
	  line[strlen(line)-1]='\0';
	  strcpy(data->name,line);
	  //printf("ADD: %s\n",data->name);
	  //dummy=fscanf(f,"%s\n",data->name);
	  dummy=fscanf(f,"%d\n",&data->type); //printf("t=%d\n",data->type);
	  //data->type=0;
	  long N=(data->type)?net->nfaces[data->type]:nv;
	  data->data=(double *)malloc(sizeof(double)*N);
	  for (i=0;i<N;i++) dummy=fscanf(f,"%lg\n",&data->data[i]);
	}
    }

  if (verbose>1) 
    {
      printf (" done. (%ld v",nv);
      for (i=1;i<=net->ndims;i++)
	{
	  if (!net->nfaces[i]) continue;
	  printf (", %ld %ld-F ",(long)net->nfaces[i],i);
	}
      printf ("and %d fields)\n",net->ndata);
    }

  fclose(f);
  free(line);

  if (delta[0]<0)
    {
      for (j=0;j<ndims;j++) net->delta[j]=net->x0[j]=net->v_coord[j];
      for (i=0;i<nv;i++)
	for (j=0;j<ndims;j++)
	  {
	    if (net->v_coord[i*ndims+j]<net->x0[j]) 
	      net->x0[j]=net->v_coord[i*ndims+j];
	    if (net->v_coord[i*ndims+j]>net->delta[j]) 
	      net->delta[j]=net->v_coord[i*ndims+j];
	  }
      for (j=0;j<ndims;j++) net->delta[j]-=net->x0[j];
    }
  else
    {
      for (i=0;i<net->ndims;i++) 
	{
	  net->x0[i]=x0[i];
	  net->delta[i]=delta[i];
	}
    }
  //for (i=0;i<=net->ndims;i++) printf (" %d ",net->nfaces[i]);
  //printf("\n");
  
  return net;
}

NDnetwork *Load_NDnetwork(const char *filename)
{
  long i,j,k;
  int dummy;
  char tag[NDNETWORK_DATA_STR_SIZE];
  FILE *f;
  int swap=0;
  NDnetwork *net;
  
  if (IsNDnetwork(filename) == 2) return Load_NDnetwork_ASCII(filename);
  
  net = calloc(1,sizeof(NDnetwork));
  memset(tag,0,NDNETWORK_DATA_STR_SIZE*sizeof(char));
  
  f=fopen(filename,"r");
  
  fread_sw(&dummy,sizeof(int),1,f,swap);
  if (dummy!=NDNETWORK_DATA_STR_SIZE) swap = 1-swap;
  fread_sw(tag,sizeof(char),NDNETWORK_DATA_STR_SIZE,f,swap);
  fread_sw(&dummy,sizeof(int),1,f,swap);
  
  if (strcmp(tag,NDNETWORK_TAG))
    {
      fclose(f);
      fprintf (stderr,"File %s has an unknown format.\n",filename);
      return NULL;
    }
  
  j=sizeof(int);
  fread_sw(&dummy,sizeof(int),1,f,swap);
  fread_sw(&net->ndims,sizeof(int),1,f,swap);
  fread_sw(&net->ndims_net,sizeof(int),1,f,swap);
  fread_sw(&dummy,sizeof(int),1,f,swap);

  if (verbose>1) printf ("Loading %dD network from file \"%s\" ...",net->ndims,filename);fflush(0);

  net->x0=malloc(net->ndims*sizeof(double));
  net->delta=malloc(net->ndims*sizeof(double));

  fread_sw(&dummy,sizeof(int),1,f,swap);
  fread_sw(net->comment,sizeof(char),80,f,swap);
  fread_sw(&net->periodicity,sizeof(int),1,f,swap);
  fread_sw(&net->isSimpComplex,sizeof(int),1,f,swap);
  fread_sw(net->x0,sizeof(double),net->ndims,f,swap);
  fread_sw(net->delta,sizeof(double),net->ndims,f,swap);
  fread_sw(&net->indexSize,sizeof(int),1,f,swap);
  if (net->indexSize != 8) net->indexSize=4;
  fread_sw(&net->cumIndexSize,sizeof(int),1,f,swap);
  if (net->cumIndexSize != 8) net->cumIndexSize=4;
  fread_sw(&net->floatSize,sizeof(int),1,f,swap);
  if (net->floatSize != 8) net->floatSize=4;
  fread_sw(net->dummy,sizeof(char),160-3*sizeof(int),f,swap);
  //fread_sw(&net->nvertex,sizeof(NDNET_UINT),1,f,swap);
  fread_sw_ului(&net->nvertex,sizeof(NDNET_UINT),1,net->indexSize,f,swap);
  fread_sw(&dummy,sizeof(int),1,f,swap);
  
  // FIXME : NO SUPPORT FOR DOUBLE COORDS IMPLEMENTED, JUST CONVERTING NOW
  // should be:
  // fread_sw_fd(net->v_coord,sizeof(NDNET_FLOAT),(size_t)net->nvertex*(size_t)net->ndims,net->floatSize,f,swap);
  //j=net->ndims*net->nvertex;
  fread_sw(&dummy,sizeof(int),1,f,swap);
  net->v_coord=malloc(sizeof(float)*(size_t)net->ndims*(size_t)net->nvertex);
  //printf("->%ld\n",(unsigned long)net->v_coord);
  net->v_coord=fread_sw_fd(net->v_coord,sizeof(float),(size_t)net->nvertex*(size_t)net->ndims,net->floatSize,f,swap);
  // whatever the input format, it is always converted single precision !
  net->floatSize = sizeof(float);
  //printf("=->%ld\n",(unsigned long)net->v_coord);
      //fread_sw(net->v_coord,sizeof(float)*(size_t)net->ndims*(size_t)net->nvertex,1,f,swap);   
  fread_sw(&dummy,sizeof(int),1,f,swap);

  //j=(1+net->ndims)*sizeof(uint);
  net->nfaces=malloc(sizeof(NDNET_UINT)*((size_t)net->ndims+1));
  fread_sw(&dummy,sizeof(int),1,f,swap);
  fread_sw_ului(net->nfaces,sizeof(NDNET_UINT),((size_t)net->ndims+1),net->indexSize,f,swap);
  fread_sw(&dummy,sizeof(int),1,f,swap);

  //j=(1+net->ndims)*sizeof(int);
  net->haveVertexFromFace=malloc(sizeof(int)*((size_t)net->ndims+1));
  fread_sw(&dummy,sizeof(int),1,f,swap);
  fread_sw(net->haveVertexFromFace,sizeof(int)*((size_t)net->ndims+1),1,f,swap);
  fread_sw(&dummy,sizeof(int),1,f,swap);

  net->f_vertexIndex = calloc(sizeof(NDNET_UINT*),(1+net->ndims));
  net->f_numVertexIndexCum = calloc(sizeof(NDNET_IDCUMT*),(1+net->ndims));

  for (i=0;i<1+net->ndims;i++)
  {    
      
      if (net->haveVertexFromFace[i])
      {
	if (!net->isSimpComplex)
	  {
	    
	    //j=sizeof(uint)*((size_t)net->nfaces[i]+1);
	    net->f_numVertexIndexCum[i]=malloc(sizeof(NDNET_IDCUMT)*((size_t)net->nfaces[i]+1));
	    fread_sw(&dummy,sizeof(int),1,f,swap);
	    fread_sw_ului(net->f_numVertexIndexCum[i],sizeof(NDNET_IDCUMT),((size_t)net->nfaces[i]+1),net->cumIndexSize,f,swap);
	    fread_sw(&dummy,sizeof(int),1,f,swap);
	    
	    //j=sizeof(uint)*((size_t)net->f_numVertexIndexCum[i][net->nfaces[i]]);
	    net->f_vertexIndex[i]=malloc(sizeof(NDNET_UINT)*((size_t)net->f_numVertexIndexCum[i][net->nfaces[i]]));
	    fread_sw(&dummy,sizeof(int),1,f,swap);
	    fread_sw_ului(net->f_vertexIndex[i],sizeof(NDNET_UINT),((size_t)net->f_numVertexIndexCum[i][net->nfaces[i]]),net->indexSize,f,swap);
	    fread_sw(&dummy,sizeof(int),1,f,swap); 
	  }
	else
	  {
	    //j=sizeof(uint)*((size_t)(i+1)*net->nfaces[i]);
	    net->f_vertexIndex[i]=malloc(sizeof(NDNET_UINT)*((size_t)(i+1)*net->nfaces[i]));
	    fread_sw(&dummy,sizeof(int),1,f,swap);
	    fread_sw_ului(net->f_vertexIndex[i],sizeof(NDNET_UINT),((size_t)(i+1)*net->nfaces[i]),net->indexSize,f,swap);
	    fread_sw(&dummy,sizeof(int),1,f,swap);
	  }
	  
      }
  }
  
  j=(1+net->ndims)*sizeof(int);
  net->haveFaceFromVertex=malloc(sizeof(int)*((size_t)net->ndims+1));
  fread_sw(&dummy,sizeof(int),1,f,swap);
  fread_sw(net->haveFaceFromVertex,sizeof(int)*((size_t)net->ndims+1),1,f,swap);
  fread_sw(&dummy,sizeof(int),1,f,swap);

  net->v_faceIndex = calloc(sizeof(NDNET_UINT*),(1+net->ndims));
  net->v_numFaceIndexCum = calloc(sizeof(NDNET_IDCUMT*),(1+net->ndims));
  //printf("%ld == %d\n",sizeof(NDNET_IDCUMT),net->cumIndexSize);
  for (i=0;i<1+net->ndims;i++)
  {     
    //printf("HELLO %ld\n",i);
      if (net->haveFaceFromVertex[i])
      {
	//j=sizeof(uint)*((size_t)net->nvertex+1);
	
	  net->v_numFaceIndexCum[i]=malloc(sizeof(NDNET_IDCUMT)*((size_t)net->nvertex+1));
	  fread_sw(&dummy,sizeof(int),1,f,swap);
	  fread_sw_ului(net->v_numFaceIndexCum[i],sizeof(NDNET_IDCUMT),((size_t)net->nvertex+1),net->cumIndexSize,f,swap);
	  fread_sw(&dummy,sizeof(int),1,f,swap);
	  
	  //j=sizeof(uint)*((size_t)net->v_numFaceIndexCum[i][net->nvertex]);
	  net->v_faceIndex[i]=malloc(sizeof(NDNET_UINT)*((size_t)net->v_numFaceIndexCum[i][net->nvertex]));
	  fread_sw(&dummy,sizeof(int),1,f,swap);
	  fread_sw_ului(net->v_faceIndex[i],sizeof(NDNET_UINT),((size_t)net->v_numFaceIndexCum[i][net->nvertex]),net->indexSize,f,swap);
	  fread_sw(&dummy,sizeof(int),1,f,swap);	  
      }
  }
  
  net->haveFaceFromFace =calloc(sizeof(int *),(1+net->ndims));
  net->f_faceIndex = calloc(sizeof(NDNET_UINT**),(1+net->ndims));
  net->f_numFaceIndexCum = calloc(sizeof(NDNET_IDCUMT**),(1+net->ndims));

  fread_sw(&dummy,sizeof(int),1,f,swap);
  for (k=0;k<net->ndims+1;k++)
  {
      net->f_faceIndex[k] = calloc(sizeof(NDNET_UINT*),(1+net->ndims));
      net->f_numFaceIndexCum[k] = calloc(sizeof(NDNET_IDCUMT*),(1+net->ndims));

      net->haveFaceFromFace[k]=malloc(sizeof(int)*((size_t)net->ndims+1));
      fread_sw(net->haveFaceFromFace[k],sizeof(int)*((size_t)net->ndims+1),1,f,swap);
  }   
  fread_sw(&dummy,sizeof(int),1,f,swap);

  for (i=0;i<1+net->ndims;i++)
  {
      for (k=0;k<1+net->ndims;k++)
      {
	  
	  if (net->haveFaceFromFace[i][k])
	  {
	    //j=sizeof(uint)*(1+(size_t)net->nfaces[i]);
	      net->f_numFaceIndexCum[i][k]=malloc(sizeof(NDNET_IDCUMT)*((size_t)net->nfaces[i]+1));
	      fread_sw(&dummy,sizeof(int),1,f,swap);
	      fread_sw_ului(net->f_numFaceIndexCum[i][k],sizeof(NDNET_IDCUMT),((size_t)net->nfaces[i]+1),net->cumIndexSize,f,swap);
	      fread_sw(&dummy,sizeof(int),1,f,swap);
	      
	      //j=sizeof(uint)*((size_t)net->f_numFaceIndexCum[i][k][net->nfaces[i]]);
	      net->f_faceIndex[i][k]=malloc(sizeof(NDNET_UINT)*((size_t)net->f_numFaceIndexCum[i][k][net->nfaces[i]]));
	      fread_sw(&dummy,sizeof(int),1,f,swap);
	      fread_sw_ului(net->f_faceIndex[i][k],sizeof(NDNET_UINT),((size_t)net->f_numFaceIndexCum[i][k][net->nfaces[i]]),net->indexSize,f,swap);
	      fread_sw(&dummy,sizeof(int),1,f,swap);	  
	  }
      }
  }

  j=sizeof(int);
  fread_sw(&dummy,sizeof(int),1,f,swap);
  fread_sw(&net->haveVFlags,sizeof(int),1,f,swap);
  fread_sw(&dummy,sizeof(int),1,f,swap);

  if (net->haveVFlags)
  {
      j=sizeof(unsigned char)*net->nvertex;
      net->v_flag=malloc(sizeof(unsigned char)*(size_t)net->nvertex);
      fread_sw(&dummy,sizeof(int),1,f,swap);
      fread_sw(net->v_flag,sizeof(unsigned char),(size_t)net->nvertex,f,swap);
      fread_sw(&dummy,sizeof(int),1,f,swap);
  }

  net->haveFFlags=calloc(net->ndims+1,sizeof(int));
  fread_sw(&dummy,sizeof(int),1,f,swap);
  fread_sw(net->haveFFlags,sizeof(int),(net->ndims+1),f,swap);
  fread_sw(&dummy,sizeof(int),1,f,swap);
  
  net->f_flag=calloc(net->ndims+1,sizeof(unsigned char *));
  
  for (i=0;i<1+net->ndims;i++)
      if (net->haveFFlags[i])
      {	  
	  j=sizeof(unsigned char)*net->nfaces[i];
	 
	  if (j) net->f_flag[i]=malloc(sizeof(unsigned char)*(size_t)net->nfaces[i]);
	  fread_sw(&dummy,sizeof(int),1,f,swap);
	  if (j) fread_sw(net->f_flag[i],sizeof(unsigned char),(size_t)net->nfaces[i],f,swap);
	  fread_sw(&dummy,sizeof(int),1,f,swap);
      }

  
  fread_sw(&dummy,sizeof(int),1,f,swap);
  fread_sw(&net->ndata,sizeof(int),1,f,swap);
  fread_sw(&dummy,sizeof(int),1,f,swap);
  

  net->data=malloc(sizeof(NDnetwork_Data)*net->ndata);
  for (i=0;i<net->ndata;i++)
  {     
    //j=sizeof(int)+255*sizeof(char);
      fread_sw(&dummy,sizeof(int),1,f,swap);
      fread_sw(&(net->data[i].type),sizeof(int),1,f,swap);
      fread_sw(net->data[i].name,sizeof(char)*255,1,f,swap);
      fread_sw(&dummy,sizeof(int),1,f,swap);

      if (net->data[i].type==0)
      {
	//j=sizeof(double)*net->nvertex;
	  net->data[i].data=malloc(sizeof(double)*(size_t)net->nvertex);
	  fread_sw(&dummy,sizeof(int),1,f,swap);
	  fread_sw(net->data[i].data,sizeof(double),(size_t)net->nvertex,f,swap);
	  fread_sw(&dummy,sizeof(int),1,f,swap);
      }
      else
      {
	//j=sizeof(double)*net->nfaces[net->data[i].type];
	  net->data[i].data=malloc(sizeof(double)*(size_t)net->nfaces[net->data[i].type]);
	  fread_sw(&dummy,sizeof(int),1,f,swap);
	  fread_sw(net->data[i].data,sizeof(double),(size_t)net->nfaces[net->data[i].type],f,swap);
	  fread_sw(&dummy,sizeof(int),1,f,swap);
      }  
  }

  fread_sw(&dummy,sizeof(int),1,f,swap);
  fread_sw(&net->nsupData,sizeof(int),1,f,swap);
  fread_sw(&dummy,sizeof(int),1,f,swap);
  
  net->supData=malloc(sizeof(NDnetwork_SupData)*net->nsupData);
  for (i=0;i<net->nsupData;i++)
  {     
      j=sizeof(int)+255*sizeof(char);
      fread_sw(&dummy,sizeof(int),1,f,swap);
      fread_sw(&(net->supData[i].type),sizeof(int),1,f,swap);
      fread_sw(net->supData[i].name,sizeof(char)*255,1,f,swap);
      fread_sw(&(net->supData[i].datasize),sizeof(int),1,f,swap);
      fread_sw(net->supData[i].datatype,sizeof(char)*255,1,f,swap);
      fread_sw(&dummy,sizeof(int),1,f,swap);

      if (net->supData[i].type==0)
      {
	//j=(size_t)net->supData[i].datasize*net->nvertex;
	  net->supData[i].data=malloc((size_t)net->supData[i].datasize*(size_t)net->nvertex);
	  fread_sw(&dummy,sizeof(int),1,f,swap);
	  fread_sw(net->supData[i].data,(size_t)net->supData[i].datasize,(size_t)net->nvertex,f,swap);
	  fread_sw(&dummy,sizeof(int),1,f,swap);
      }
      else
      {
	//j=(size_t)net->supData[i].datasize*net->nfaces[net->supData[i].type];
	  net->supData[i].data=malloc((size_t)net->supData[i].datasize*(size_t)net->nfaces[net->supData[i].type]);
	  fread_sw(&dummy,sizeof(int),1,f,swap);
	  fread_sw(net->supData[i].data,(size_t)net->supData[i].datasize,(size_t)net->nfaces[net->supData[i].type],f,swap);
	  fread_sw(&dummy,sizeof(int),1,f,swap);
      }  
  }

  fclose(f);
  if (verbose>1) printf (" done.\n");
  return net;
}


// Computes the lists of type 'type' faces the vertice
// belong to given we have the list of vertice that compose
// all type 'type' faces. 
int ComputeFaceFromVertex(NDnetwork *net, int type)
{
    NDNET_UINT nvert;
    NDNET_UINT *vlist;
    NDNET_IDCUMT *nfaces;
    NDNET_UINT nalloc=((1<<2)-1);
    NDNET_UINT **vface;
    long i,j,k,v;

    if ((type>net->ndims)||(type<1)) return -1;

    if (net->haveFaceFromVertex[type]) return 1;
    
    if (ComputeVertexFromFace(net,type)<0)
    {
	fprintf(stderr,"I cannot compute type %d faces from vertice.\n",type);
	fprintf(stderr,"I need vertice of type %d faces to do that ...\n",type);
	return -1;
    }

    if (verbose>1) printf ("Computing %d-faces lists for vertice ... ",type);fflush(0);

    //calloc here v_numFaceIndexCum[type] and v_faceIndex[type]
    vface=calloc(net->nvertex,sizeof(NDNET_UINT*));
    nfaces = net->v_numFaceIndexCum[type] = calloc(net->nvertex+1,sizeof(NDNET_IDCUMT));

    for (i=0;i<net->nfaces[type];i++)
    {
	nvert=NUM_VERTEX_IN_FACE(net,type,i);
	vlist=VERTEX_IN_FACE(net,type,i);
	
	for (j=0;j<nvert;j++)
	{
	    v=vlist[j];

	    if (!(nfaces[v]&nalloc)) 
	      vface[v]=realloc(vface[v],sizeof(NDNET_UINT)*((size_t)nfaces[v]+nalloc+1));

	    vface[v][nfaces[v]++] = i;
	}
    }
    
    //reallocate the proper amount of memory (not to waste too much for high nalloc values)
    if (nalloc!=1)
	for (i=0;i<net->nvertex;i++)
	    vface[i]=realloc(vface[i],sizeof(NDNET_UINT)*(size_t)nfaces[i]);
    //if (type==1) exit(0);
    //Cumulative count ...
    for (i=0,j=0;i<=net->nvertex;i++)
    {
	k=nfaces[i];
	nfaces[i]=j;
	j+=k;
    }
    //if (type==2) exit(0);
    //reallocate v_faceIndex properly ...
    net->v_faceIndex[type] = malloc((size_t)sizeof(NDNET_UINT) * (size_t)nfaces[net->nvertex]);
    for (i=0;i<net->nvertex;i++)
    {
      //for (j=0;j<nfaces[i+1]-nfaces[i];j++)
      //net->v_faceIndex[type][nfaces[i]+j] = (NDNET_UINT) vface[i][j];
      memcpy(&(net->v_faceIndex[type][nfaces[i]]),vface[i],sizeof(NDNET_UINT)*(nfaces[i+1]-nfaces[i]));
      free(vface[i]);
    }

    free(vface);
    
    /*
    if (net->haveVFlags)
	
    {
	NDNetFlags_enable(net,type);
	for (i=0;i<net->nfaces[type];i++)
	{
	    nvert=NUM_VERTEX_IN_FACE(net,type,i);
	    vlist=VERTEX_IN_FACE(net,type,i);
	    for (j=0;j<nvert;j++)
	    {
		if (net->v_flag[vlist[j]]&NDNETFLAG_OUT)
		    net->f_flag[type][i]|=NDNETFLAG_OUT;
		else
		    net->f_flag[type][i]&=~NDNETFLAG_OUT;
	    }
	}
    }
*/
    if (verbose>1) printf ("done.\n");fflush(0);
    net->haveFaceFromVertex[type]=1;
    return 0;
}

inline int CNP(int n,int p)
{
  unsigned long int i;
  unsigned long int up=1;
  unsigned long int down=1;

  for (i=n;i>p;i--)
    up*=i;
  for (i=n-p;i>0;i--)
    down*=i;

  return up/down;
}


int ComputeVertexFromFace(NDnetwork *net,int type)
{
    NDNET_UINT nfaces;
    NDNET_UINT *vlist;
    NDNET_UINT *vlist_ref;
    NDNET_UINT i,j,k,l,n;
    long ref_type;
    long nv_min;	
    long nv_max;
    long *sorted_id;
    
    if ((type>net->ndims)||(type<1)) return -1;

    if (net->haveVertexFromFace[type]) return 1;
    
    if (!net->isSimpComplex)
    {
	fprintf(stderr,"I cannot compute vertice list from %d type faces.\n",type);
	fprintf(stderr,"This was only implemented for simplices so far ...\n");
	return -1;
    }
    
    for (ref_type=type+1;ref_type<=net->ndims;ref_type++)
    	if (ComputeFaceFromVertex(net,ref_type)>=0) break;
    
    if (ref_type>net->ndims)
    {
	fprintf(stderr,"I cannot compute vertice list from %d type faces.\n",type);
	fprintf(stderr,"I need vertice in type_ref>%d type faces to do that ...\n",type);
	return -1;
    }
   
    //if (ref_type!=type+1) ComputeVertexFromFace(net,type+1);

    if (verbose>1) printf ("Computing vertice list for %d-faces ... ",type);fflush(0);

    
    // look for all possible combinations of NUM_VERTEX_IN_FACE(net,type, ) vertice
    // among NUM_VERTEX_IN_FACE(net,ref_type==type+1, ).
    nv_min=NUM_VERTEX_IN_FACE(net,type,0);
    nv_max=NUM_VERTEX_IN_FACE(net,ref_type,0);

    nfaces = CNP(nv_max,nv_min) * net->nfaces[ref_type];
    vlist = net->f_vertexIndex[type] = malloc ((size_t)nfaces * nv_min  *sizeof(NDNET_UINT));
    
    if (verbose>1) printf ("\rComputing vertice list for %d-faces ... (full list)",type);fflush(0);
    for (i=(1<<NUM_VERTEX_IN_FACE(net,type,0))-1;i<(1<<NUM_VERTEX_IN_FACE(net,ref_type,0))-1;i++) 
    {
	int id[nv_max];
	
	for (k=0,n=0;k<nv_max;k++)
	{
	    if ((1<<k)&i) 
		id[n++]=k;
	    
	}
	//That's one possible combination :)
	if (n==nv_min)
	{
	    
	    // so we compute for all ref_type faces the type faces
	    // corresponding to this combination 
	    for (j=0;j<net->nfaces[ref_type];j++)
	    {
		vlist_ref=VERTEX_IN_FACE(net,ref_type,j);
	
		for (k=0,l=net->nvertex;k<nv_min;k++)
		{
		    vlist[k]=vlist_ref[id[k]];
		    if (vlist[k]<l) {l=vlist[k];n=k;}
		}

		if (n!=0) {l=vlist[0];vlist[0]=vlist[n];vlist[n]=l;}

		vlist += nv_min;
	    }	    
	}    
    }
    
    vlist = net->f_vertexIndex[type];

    if (verbose>1) printf ("\rComputing vertice list for %d-faces ... (sorting)     ",type);fflush(0);
    //that's cool but we have too many faces now :(
    sorted_id = LongSortArrayUI(vlist,nfaces,nv_min,0);
    
    if (verbose>1) printf ("\rComputing vertice list for %d-faces ... (removing dup.)     ",type);fflush(0);
    //so lets remove the duplicates
    for (i=0;i<nfaces-1;i++)
    {
	if (vlist[sorted_id[i]*nv_min]==net->nvertex+1) continue;
	j=i+1;
	while (vlist[sorted_id[j]*nv_min]==net->nvertex+1) if (++j >=nfaces) break;
	
	if (j<nfaces)
	    do {
		for (k=0;k<nv_min;k++)
		{
		    for (l=0;l<nv_min;l++)
			if (vlist[sorted_id[j]*nv_min+l]==vlist[sorted_id[i]*nv_min+k])
			    break;
		    if (l==nv_min) break;
		}

		if (k==nv_min)
		    vlist[sorted_id[j]*nv_min]=net->nvertex+1; //sorted_id[j] is a duplicate

		if (++j >=nfaces) break;
		while (vlist[sorted_id[j]*nv_min]==net->nvertex+1) if ((++j) >=nfaces) break;
		// en fait faut en trouver 2 AU PLUS si type==ref_type-1 ...
	    } while ((j<nfaces)&&(vlist[sorted_id[j]*nv_min]==vlist[sorted_id[i]*nv_min]));
    }

    free(sorted_id);

    // and reallocate the list properly ...
    for (i=0,j=0;i<nfaces;i++)
    {
	if (vlist[i*nv_min]!=net->nvertex+1)
	{
	    if (i!=j) 
		memcpy(&vlist[j*nv_min],&vlist[i*nv_min],(size_t)nv_min*sizeof(NDNET_UINT));
	    j++;
	}
    }
  
    nfaces=j;
    vlist = realloc(vlist,(size_t)nv_min*nfaces*sizeof(NDNET_UINT));
    net->nfaces[type]=nfaces;
    net->f_vertexIndex[type] = vlist;
    
    net->haveVertexFromFace[type]=1;
    if (verbose>1) printf ("\rComputing vertice list for %d-faces ... done.          \n",type);fflush(0);

    return 0;
}


// (*list) is set to the list of faces of type 'type' that belong to (or encompass)
// the face 'ref_face' of type 'ref_type'. Depending on the data available in 'net', 
// the list can be obtained by a simple memcpy or can necessitate some computations.
//
// The function returns the number of such faces (or -1 in case of failure).
// (*list) should be NULL or already allocated, this function will reallocate it apropriately. 
int FacesInFace(NDnetwork *net, NDNET_UINT ref_face, int ref_type,int type, NDNET_UINT **list)
{
    long i,m;
    long j,n;
    long k,o;
    long l;
    long nfaces,nfaces_alloc;
    NDNET_UINT *vlist;
    NDNET_UINT *wlist;
    NDNET_UINT *xlist;
    
    if (type==ref_type) {*list=realloc(*list,sizeof(NDNET_UINT));**list=ref_face;return 1;}

    if (net->haveFaceFromFace[ref_type][type])
    {
	n=NUM_FACE_IN_FACE(net,ref_type,ref_face,type);
	*list=realloc(*list,sizeof(NDNET_UINT)*n);
	memcpy(*list,FACE_IN_FACE(net,ref_type,ref_face,type),sizeof(NDNET_UINT)*n);
	return n;
    }
    
    if (((ComputeFaceFromVertex(net,type)<0)&&(type!=0))||((ComputeVertexFromFace(net,ref_type)<0)&&(ref_type!=0)))
    {
	fprintf (stderr,"ERROR: in FacesInFace, not enough information in network to retrieve\n");
	fprintf (stderr,"       faces of type %d that belong to faces of type %d. \n",type,ref_type);
	exit(0);
	return -1;
    }

    if (ref_type==0)
    {
	n=NUM_FACE_IN_VERTEX(net,type,ref_face);
	//printf ("n=%d\n",n);
	*list=realloc(*list,sizeof(NDNET_UINT)*n);
	memcpy(*list,FACE_IN_VERTEX(net,type,ref_face),sizeof(NDNET_UINT)*n);
	return n;
    }

    if (type==0)
    {
	n=NUM_VERTEX_IN_FACE(net,ref_type,ref_face);
	//printf ("n=%d\n",n);
	*list=realloc(*list,sizeof(NDNET_UINT)*n);
	memcpy(*list,VERTEX_IN_FACE(net,ref_type,ref_face),sizeof(NDNET_UINT)*n);
	return n;
    }

    m = NUM_VERTEX_IN_FACE(net,ref_type,ref_face);
    vlist = VERTEX_IN_FACE(net,ref_type,ref_face);
    nfaces=0;
    nfaces_alloc=100;
    *list=realloc(*list,sizeof(NDNET_UINT)*nfaces_alloc);

    if (type<ref_type)
    {
	for (i=0;i<m-1;i++)
	{
	    n=NUM_FACE_IN_VERTEX(net,type,vlist[i]);
	    wlist=FACE_IN_VERTEX(net,type,vlist[i]);

	    for (j=0;j<n;j++)
	    {
		o=NUM_VERTEX_IN_FACE(net,type,wlist[j]);
		xlist=VERTEX_IN_FACE(net,type,wlist[j]);
		
		for (k=0;k<o;k++)
		{
		    for (l=0;l<m;l++)
			if (xlist[k]==vlist[l]) break;
		    if (l==m) break;
		}
		// we have a face that belongs to ref_face ...
		if (k==o) 
		{
		    for (k=0;k<nfaces;k++)
			if ((*list)[k] == wlist[j]) break;
		    if (k==nfaces)
		    {
			nfaces++;
			if (nfaces>nfaces_alloc)
			{
			    nfaces_alloc+=100;
			    *list=realloc(*list,sizeof(NDNET_UINT)*nfaces_alloc);
			}
			
			(*list)[nfaces-1]=wlist[j];
		    }
		}
	    }
	}
    }
    else 
    {
	for (i=0;i<m;i++)
	{
	    n=NUM_FACE_IN_VERTEX(net,type,vlist[i]);
	    wlist=FACE_IN_VERTEX(net,type,vlist[i]);
	    
	    for (j=0;j<n;j++)
	    {
		o=NUM_VERTEX_IN_FACE(net,type,wlist[j]);
		xlist=VERTEX_IN_FACE(net,type,wlist[j]);
		
		for (k=0;k<m;k++)
		{
		    for (l=0;l<o;l++)
			if (xlist[l]==vlist[k]) break;
		    if (l==o) break;
		}
		// we have a face that encompasses ref_face ...
		if (k==m) 
		{
		    for (k=0;k<nfaces;k++)
			if ((*list)[k] == wlist[j]) break;
		    if (k==nfaces)
		    {
			nfaces++;
			if (nfaces>nfaces_alloc)
			{
			    nfaces_alloc+=100;
			    *list=realloc(*list,sizeof(NDNET_UINT)*nfaces_alloc);
			}
			(*list)[nfaces-1]=wlist[j];
		    }
		}
	    }
	}
    }
    
    return nfaces;
}

// vertice if (type==0)
// interface returns the surface of the interface with the different neighbours ... 
int NDFindNeighbours(NDnetwork *net, int type, NDNET_UINT index, NDNET_UINT **list, double **interface)
{
    long n=0;
    NDNET_UINT *flist=NULL;
    NDNET_UINT *nlist=NULL;
    NDNET_UINT  i,j;
    NDNET_UINT ntot;
    long cur;

    if (type<0) 
    {
	fprintf (stderr,"In NDFindNeighbours, 'type=%d', should be be 0 or positive. \n",type);
	return -1;
    }


    if (type==0)
    {
	float *coord[2];

	ComputeFaceFromVertex(net,1);

	n=NUM_FACE_IN_VERTEX(net,1,index);
	flist=FACE_IN_VERTEX(net,1,index);

	*list = realloc(*list,n*sizeof(NDNET_UINT));
	if (interface!=NULL) 
	    *interface = realloc(*interface,n*sizeof(double));
	
	coord[0]=&net->v_coord[index*net->ndims];
	
	for (i=0;i<n;i++)
	{
	  //int nin=0;
	    
	    nlist=VERTEX_IN_FACE(net,1,flist[i]);
	    
	    if (nlist[0]==index)
		(*list)[i]=nlist[1];
	    else
		(*list)[i]=nlist[0];

	    coord[1]=&net->v_coord[(*list)[i]*net->ndims];
	    
	    if (interface!=NULL) 
		(*interface)[i] =  SimplexVolumef(coord,2,net->ndims);
	   
	}

	return n;
    }
    else if (type==1) 
    {
	fprintf (stderr,"NDFindNeighbours not implemented for type=1 sor far. \n");
	return -1;
    }
   
    n=FacesInFace(net, index, type, type-1, &flist);
    ntot=n;
    *list = realloc(*list,ntot*sizeof(NDNET_UINT));
    if (interface!=NULL) 
	*interface = realloc(*interface,n*sizeof(double));
    
    cur=0;
    // actually find the neighbours (which share a (type-1) face with face index)
    for (i=0;i<n;i++) 
    {
	int nnei;
	
	nnei = FacesInFace(net, flist[i], type-1, type, &nlist);
	if (nnei<2) continue;
	else if (nnei==2)
	{
	    if (ntot==cur) 
	    {
		ntot++;*list = realloc(*list,ntot*sizeof(NDNET_UINT));
	    }
	    if (nlist[0]==index)
		(*list)[cur]=nlist[1];
	    else
		(*list)[cur]=nlist[0];
	    cur++;
	}
	else
	{
	    for (j=0;j<nnei;j++)
	    {
		if (nlist[j]==index) continue;
		if (ntot==cur) 
		{
		    ntot++;*list = realloc(*list,ntot*sizeof(NDNET_UINT));
		}
		(*list)[cur++]=nlist[j];
	    }
	}

	if (interface!=NULL) 
	{
	    
	    long nvert = NUM_VERTEX_IN_FACE(net,type-1,(*list)[i]);
	    NDNET_UINT *vlist = VERTEX_IN_FACE(net,type-1,(*list)[i]);
	    float *coord[nvert];
	 
	    for (j=0;j<nvert;j++)
		coord[j] = &net->v_coord[vlist[j]*net->ndims];
    
	    if (interface!=NULL) 
		(*interface)[i] =  SimplexVolumef(coord,nvert,net->ndims);
	}
    }

    free(flist);
    free(nlist);
    
    return cur;
}

void NDnet_smooth(NDnetwork *net, int n)
{
    
    double *tmp = (double*)malloc(net->nvertex*net->ndims*sizeof(double));
    int *nn =  (int*)malloc(net->nvertex*sizeof(int));
    long i,j,k,l,m,d;
    
    //double x[net->ndims+1];

    //memcpy(tmp,net->v_coord,net->nvertex*net->ndims*sizeof(float));
    for (k=0;k<n;k++)
    {
	memset(nn,0,net->nvertex*sizeof(int));
	memset(tmp,0,net->nvertex*net->ndims*sizeof(double));
	long ntot=0;

	for (d=1;d<=net->ndims;d++)
	  {
	    for (i=0;i<net->nfaces[d];i++)
	      {    
		NDNET_UINT *v=VERTEX_IN_FACE(net,d,i);
		for (l=0;l<d;l++)
		  for (m=l+1;m<d+1;m++)
		    {
		      for (j=0;j<net->ndims;j++)
			{
			  double x0 = net->v_coord[v[l]*(net->ndims)+j];
			  double x1 = net->v_coord[v[m]*(net->ndims)+j];
			  
			  if (x1-x0>net->delta[j]*0.5)
			    {
			      if (nn[v[l]]>=0) tmp[v[l]*(net->ndims)+j] += x1-net->delta[j];
			      if (nn[v[m]]>=0) tmp[v[m]*(net->ndims)+j] += x0+net->delta[j];
			    }
			  else if (x0-x1>net->delta[j]*0.5)
			    {
			      if (nn[v[l]]>=0) tmp[v[l]*(net->ndims)+j] += x1+net->delta[j];
			      if (nn[v[m]]>=0) tmp[v[m]*(net->ndims)+j] += x0-net->delta[j];
			    }
			  else
			    {
			      if (nn[v[l]]>=0) tmp[v[l]*(net->ndims)+j] +=x1;
			      if (nn[v[m]]>=0) tmp[v[m]*(net->ndims)+j] +=x0;
			    }
			}
		      if (nn[v[l]]>=0) nn[v[l]]++;
		      if (nn[v[m]]>=0) nn[v[m]]++;  
		    }
	      }

	    for (i=0;i<net->nvertex;i++)
	      {
		if (nn[i]>0)
		  {
		    for (j=0;j<net->ndims;j++)
		      net->v_coord[i*(net->ndims)+j] = tmp[i*(net->ndims)+j]/nn[i];
		    nn[i]=-1;
		    ntot++;
		  }
	      }

	    if (ntot == net->nvertex) break;
	  }
    }

    free(nn);
    free(tmp);    
}

int NDNetFlags_enable(NDnetwork *net,int type)
{
    if (type==0)
    {
	if (!net->haveVFlags)
	{
	    net->haveVFlags=1;
	    net->v_flag=(unsigned char*)calloc((size_t)net->nvertex,sizeof(unsigned char));
	    return 0;
	}
    }
    else
    {
	if (!net->haveFFlags[type])
	{
	    net->haveFFlags[type]=1;
	    net->f_flag[type] = (unsigned char*)calloc(net->nfaces[type],sizeof(unsigned char));
	    return 0;
	}
    }

    return 1;
}

int Simplex_CS(NDnetwork *net,int type, long index, float *result)
{
  long i,j;

  int n = NUM_VERTEX_IN_FACE(net,type,index);
  NDNET_UINT *vert= VERTEX_IN_FACE(net,type,index);
  double v[net->ndims+1][net->ndims];
  double *pv[net->ndims+1];
  double center[net->ndims];
   
  if ((n!=net->ndims+1)&&(n!=2)) {
     fprintf(stderr,"ERROR: cannot compute circumcircle, ndims=%d, n=%d\n",net->ndims,n);
     fprintf(stderr,"NOT IMPLEMENTED YET !!! (using COM instead)\n");
     return Simplex_COM(net,type,index,result);
     //exit(0);
   }
   
   for (i=0;i<n;i++)
     {
       pv[i]=&v[i][0];
       for (j=0;j<net->ndims;j++) v[i][j] = net->v_coord[net->ndims * vert[i] + j];
     }

   if (net->periodicity) {
     int region=0;
     for (i=0;i<net->ndims+1;i++)
       {
	 for (j=0;j<net->ndims;j++)
	   {
	     if (fabs(v[i][j]-net->x0[j])>2*net->delta[j]/3) region|=(1<<(2*j+1));
	     if (fabs(v[i][j]-net->x0[j])<net->delta[j]/3) region|=(1<<(2*j));
	   }
       }
     for (j=0;j<net->ndims;j++)
       {
	 if ((region&(1<<(2*j+1)))&&(region&(1<<(2*j))))
	   {
	     for (i=0;i<net->ndims+1;i++)
	       if (fabs(v[i][j]-net->x0[j])>net->delta[j]/2)
		 v[i][j]-=net->delta[j];
	   }
       }
   }

   if (n==2)
    {
      for (i=0;i<net->ndims;i++)
	result[i] = 0.5*(v[0][i]+v[1][i]);
      //net->v_coord[net->ndims * vert[0]+i]+
      //net->v_coord[net->ndims * vert[1]+i]);
      return 0;
    }


   SimplexSphere(pv,n-1,center);
   for (i=0;i<net->ndims;i++) result[i]=center[i];
   
   return 0;
}

int Simplex_CS_fc(NDnetwork *net,int type, double **pv, float *result)
{
  int i,j;

  int n = type+1;
  //int n = NUM_VERTEX_IN_FACE(net,type,index);
  //uint *vert= VERTEX_IN_FACE(net,type,index);
  //double v[net->ndims+1][net->ndims];
  //double *pv[net->ndims+1];
  double center[net->ndims];
  /* 
  if ((n!=net->ndims+1)&&(n!=2)) {
     fprintf(stderr,"ERROR: cannot compute circumcircle, ndims=%d, n=%d\n",net->ndims,n);
     fprintf(stderr,"NOT IMPLEMENTED YET !!! (using COM instead)\n");
     return -1;
     //return Simplex_COM(net,type,index,result);
     //exit(0);
   }
*/

 

  /*
   for (i=0;i<n;i++)
     {
       pv[i]=&v[i][0];
       //for (j=0;j<net->ndims;j++) v[i][j] = net->v_coord[net->ndims * vert[i] + j];
     }
  */

   if (net->periodicity) {
     int region=0;
     for (i=0;i<n;i++)
       {
	 for (j=0;j<net->ndims;j++)
	   {
	     if (fabs(pv[i][j]-net->x0[j])>2*net->delta[j]/3) region|=(1<<(2*j+1));
	     if (fabs(pv[i][j]-net->x0[j])<net->delta[j]/3) region|=(1<<(2*j));
	   }
       }
     for (j=0;j<net->ndims;j++)
       {
	 if ((region&(1<<(2*j+1)))&&(region&(1<<(2*j))))
	   {
	     for (i=0;i<n;i++)
	       if (fabs(pv[i][j]-net->x0[j])>net->delta[j]/2)
		 pv[i][j]-=net->delta[j];
	   }
       }
   }
    if (n==1)
    {
      for (j=0;j<net->ndims;j++) result[j]=pv[0][j];
    }
    else if (n==2)
    {
      for (i=0;i<net->ndims;i++)
	result[i] = 0.5*(pv[0][i]+pv[1][i]);
      return 0;
    }


   SimplexSphere(pv,n-1,center);
   for (i=0;i<net->ndims;i++) result[i]=center[i];
   
   return 0;
}


// WARNING : periodic boundary conditions not properly implemented yet
// See Simplex_CS() ...
int Simplex_COM(NDnetwork *net,int type, long index, float *result)
{
    long i,j,n;
    NDNET_UINT *vlist;
    float *coord;
    float tmp;
    //int region=0;
    //double val[net->ndims];
    
    n=NUM_VERTEX_IN_FACE(net,type,index);
    vlist=VERTEX_IN_FACE(net,type,index);
    
    for (i=0;i<net->ndims;i++) result[i]=0;
    
    for (j=0;j<n;j++)
      {
	coord = &(net->v_coord[net->ndims*vlist[j]]);
	for (i=0;i<net->ndims;i++)	
	  result[i] += *(coord++);
      }
    
    tmp=1./n;
    for (i=0;i<net->ndims;i++) result[i]*=tmp;

    return 0;
}

int tagBoundaries(NDnetwork *net)
{
  long i,j,t;
  NDNET_UINT n;
  int a,b;
  float *vpos=net->v_coord;
  printf ("* Tagging network boundary ... ");fflush(0);
  
  NDNetFlags_enable(net,0);
  n=0;
  for (i=0;i<net->nvertex;i++)
    {
      net->v_flag[i] &= (~NDNETFLAG_BOUNDARY);
      if (net->v_flag[i]&NDNETFLAG_OUT) n++;
      
      for (j=0;j<net->ndims;j++)
	{
	  if (vpos[j+net->ndims*i]< net->x0[j]) 
	    net->v_flag[i] |= (NDNETFLAG_BOUNDARY);
	  if (vpos[j+net->ndims*i]> net->x0[j]+net->delta[j]) 
	    net->v_flag[i] |= (NDNETFLAG_BOUNDARY);
	}
    }
  //printf("n=%d\n",n);
  for (t=1;t<net->ndims+1;t++)
    {
      NDNetFlags_enable(net,t);
      
      for (i=0;i<net->nfaces[t];i++)
	net->f_flag[t][i] &= (~NDNETFLAG_BOUNDARY);
    }
  
  for (i=0,n=0;i<net->nfaces[1];i++)
    {
      NDNET_UINT *v = VERTEX_IN_FACE(net,1,i);
      
      a=(net->v_flag[v[0]]&NDNETFLAG_OUT);
      b=(net->v_flag[v[1]]&NDNETFLAG_OUT);
      if (a!=b)
	{
	  if (b) net->v_flag[v[0]]|=NDNETFLAG_BOUNDARY;
	  else net->v_flag[v[1]]|=NDNETFLAG_BOUNDARY;
	}				
    }
  for (i=0,n=0;i<net->nvertex;i++)
    if (net->v_flag[i]&NDNETFLAG_BOUNDARY)
      n++;
  
  printf("(%ld V)",(long)n);fflush(0);
  for (t=1;t<net->ndims+1;t++)
    {
      n=0;
      for (i=0;i<net->nfaces[t];i++)
	{
	  NDNET_UINT *v = VERTEX_IN_FACE(net,t,i);
	  
	  for (j=0;j<t+1;j++)
	    {
	      if ( (net->v_flag[v[j]]&NDNETFLAG_BOUNDARY) )
		break;
	    }
	  
	  if (j!=t+1) {net->f_flag[t][i]|=NDNETFLAG_BOUNDARY;n++;}
	}
      printf("(%ld %ld-F)",(long)n,(long)t);fflush(0);
    }
  
  for (t=1;t<net->ndims+1;t++)
    {
      n=0;
      for (i=0;i<net->nfaces[t];i++)
	{
	  NDNET_UINT *v = VERTEX_IN_FACE(net,t,i);
	  
	  for (j=0;j<t+1;j++)
	    {
	      if (! (net->v_flag[v[j]]&NDNETFLAG_OUT) )
		break;
	    }
	  
	  if (j==t+1) {net->f_flag[t][i]|=NDNETFLAG_OUT;n++;}
	}
      //printf("(%d %d-faces)",n,t);fflush(0);
    }
  printf(" done.\n");
  return 1;
}


int getPeriodicityCut(NDnetwork *net,int type, float *coords, int unfold)
{
  int flag=0;
  int i,j;
  int n=type+1;

  if (type==1)
    {
      for (j=0;j<net->ndims;j++)
	{
	  if (fabs(coords[j]-coords[j+net->ndims])>0.5*net->delta[j]) 
	    flag|=(1<<j);
	}
    }
  else {
    for (i=0;i<n;i++)
      for (j=i+1;j<n;j++)
	{
	  float new_coord[2*net->ndims];
	  memcpy(new_coord,&coords[net->ndims*i],net->ndims*sizeof(float));
	  memcpy(&new_coord[net->ndims],&coords[net->ndims*j],net->ndims*sizeof(float));
	  flag|=getPeriodicityCut(net,1,new_coord,0);
	}
  }

  flag &= net->periodicity;

  if (unfold) 
    {
      for (j=0;j<net->ndims;j++)
	{
	  if (flag&(1<<j)) {
	    for (i=0;i<n;i++)
	      if (coords[j+i*net->ndims]-net->x0[j]>0.5*net->delta[j]) coords[j+i*net->ndims]-=net->delta[j];
	  }
	}
    }

  return flag;
}

int getSimplexCoords(NDnetwork *net,int type, NDNET_UINT index, float *coords, int unfold)
{
  if (!net->haveVertexFromFace[type])
    ComputeVertexFromFace(net,type);
    
  int n=NUM_VERTEX_IN_FACE(net,type,index);
  NDNET_UINT *vlist=VERTEX_IN_FACE(net,type,index);
  long i;

  for (i=0;i<n;i++) 
    //for (j=0;j<net->ndims;j++) 
      memcpy(&coords[net->ndims*i],&net->v_coord[net->ndims*vlist[i]],net->ndims*sizeof(float));
  
  if (net->f_flag[type]!=NULL)
    if (unfold&&(net->f_flag[type][index]&NDNETFLAG_PERIODIC_CUT))
      {
	return getPeriodicityCut(net,type,coords,1);
      }
 
  return 0;
}

int setPeriodicCutFlag(NDnetwork *net)
{
  long i,j;
  //double d;
  //NDNET_UINT *list=NULL;
  float coords[(net->ndims+1)*net->ndims];
  //int flags[net->nfaces[1]];

  printf ("Tagging %d\n",net->periodicity);
  if (net->periodicity==0)
    return 0;
  
  for (j=1;j<=net->ndims;j++)
    {
      if (!net->haveVertexFromFace[j]) continue;
      NDNetFlags_enable(net,j);
      for (i=0;i<net->nfaces[j];i++)
	{
	  net->f_flag[j][i]&= ~NDNETFLAG_PERIODIC_CUT;
	  getSimplexCoords(net,j,i,coords,0);
	  if (getPeriodicityCut(net,j,coords,0))
	    {
	      net->f_flag[j][i]|= NDNETFLAG_PERIODIC_CUT;
	      //printf ("tagged %d %d\n",i,j);
	    }
	}
    }
  
  return 1;
}

NDnetwork *Load_NDnetwork_header(const char *filename, NDnet_fstruct_info *info)
{
  long i,j,k;
  int dummy;
  char tag[NDNETWORK_DATA_STR_SIZE];
  FILE *f;
  int swap=0;
  NDnetwork *net;
  
  //if (IsNDnetwork(filename) == 2) return Load_NDnetwork_ASCII(filename);
  
  net =(NDnetwork *) calloc(1,sizeof(NDnetwork));
  memset(tag,0,NDNETWORK_DATA_STR_SIZE*sizeof(char));
  
  f=fopen(filename,"r");
  
  fread_sw(&dummy,sizeof(int),1,f,swap);
  if (dummy!=NDNETWORK_DATA_STR_SIZE) swap = 1-swap;
  fread_sw(tag,sizeof(char),NDNETWORK_DATA_STR_SIZE,f,swap);
  fread_sw(&dummy,sizeof(int),1,f,swap);

  info->swap=swap;
  
  if (strcmp(tag,NDNETWORK_TAG))
    {
      fclose(f);
      fprintf (stderr,"ERROR in Load_NDnetwork_header: File %s has an unknown format.\n",filename);
      return NULL;
    }
  
  j=sizeof(int);
  fread_sw(&dummy,sizeof(int),1,f,swap);
  fread_sw(&net->ndims,sizeof(int),1,f,swap);
  fread_sw(&net->ndims_net,sizeof(int),1,f,swap);
  fread_sw(&dummy,sizeof(int),1,f,swap);
 
  if (verbose>1) printf ("Loading %dD network from file \"%s\" ...",net->ndims,filename);fflush(0);

  net->x0=(double*)calloc(net->ndims,sizeof(double));
  net->delta=(double*)calloc(net->ndims,sizeof(double));

  fread_sw(&dummy,sizeof(int),1,f,swap);
  fread_sw(net->comment,sizeof(char),80,f,swap);
  fread_sw(&net->periodicity,sizeof(int),1,f,swap);
  fread_sw(&net->isSimpComplex,sizeof(int),1,f,swap);
  fread_sw(net->x0,sizeof(double),net->ndims,f,swap);
  fread_sw(net->delta,sizeof(double),net->ndims,f,swap);
  fread_sw(&net->indexSize,sizeof(int),1,f,swap);
  if (net->indexSize != 8) net->indexSize=4;
  fread_sw(&net->cumIndexSize,sizeof(int),1,f,swap);
  if (net->cumIndexSize != 8) net->cumIndexSize=4;
  fread_sw(&net->floatSize,sizeof(int),1,f,swap);
  if (net->floatSize != 8) net->floatSize=4;
  //fread_sw(net->dummy,sizeof(char),152,f,swap);
  fread_sw(net->dummy,sizeof(char),160-3*sizeof(int),f,swap);
  fread_sw_ului(&net->nvertex,sizeof(NDNET_UINT),1,net->indexSize,f,swap);
  fread_sw(&dummy,sizeof(int),1,f,swap);

  fread_sw(&dummy,sizeof(int),1,f,swap);
  info->s_v_coord=sizeof(float)*(size_t)net->ndims*(size_t)net->nvertex;
  info->rss_v_coord=sizeof(float);
  info->rn_v_coord=(size_t)net->ndims*(size_t)net->nvertex;
  info->rs_v_coord=sizeof(float);
  fgetpos(f,&info->v_coord);
  fseek(f,info->rn_v_coord*info->rs_v_coord,SEEK_CUR);
  net->v_coord=NULL;
    // here
  fread_sw(&dummy,sizeof(int),1,f,swap);
  
  net->nfaces=(NDNET_UINT *)calloc((size_t)net->ndims+1,sizeof(NDNET_UINT));
  fread_sw(&dummy,sizeof(int),1,f,swap);
  fread_sw_ului(net->nfaces,sizeof(NDNET_UINT),((size_t)net->ndims+1),net->indexSize,f,swap);
  fread_sw(&dummy,sizeof(int),1,f,swap);

  net->haveVertexFromFace=(int*)calloc((size_t)net->ndims+1,sizeof(int));
  fread_sw(&dummy,sizeof(int),1,f,swap);
  fread_sw(net->haveVertexFromFace,sizeof(int)*((size_t)net->ndims+1),1,f,swap);
  fread_sw(&dummy,sizeof(int),1,f,swap);

  net->f_vertexIndex = (NDNET_UINT**)calloc(1+net->ndims,sizeof(NDNET_UINT*));
  net->f_numVertexIndexCum = (NDNET_IDCUMT**)calloc(1+net->ndims,sizeof(NDNET_IDCUMT*));

  for (i=0;i<1+net->ndims;i++)
    {    
	
      if (net->haveVertexFromFace[i])
	{
	  if (!net->isSimpComplex)
	    {
	      fread_sw(&dummy,sizeof(int),1,f,swap);		
	      info->s_f_numVertexIndexCum[i]=sizeof(NDNET_IDCUMT)*((size_t)net->nfaces[i]+1);
	      info->rss_f_numVertexIndexCum[i]=sizeof(NDNET_IDCUMT);
	      info->rn_f_numVertexIndexCum[i]=((size_t)net->nfaces[i]+1);
	      info->rs_f_numVertexIndexCum[i]=net->cumIndexSize;
	      fgetpos(f,&info->f_numVertexIndexCum[i]);
	      fseek(f,info->rn_f_numVertexIndexCum[i]*info->rs_f_numVertexIndexCum[i],SEEK_CUR);
	      net->f_numVertexIndexCum[i]=NULL;
	      // here
	      fread_sw(&dummy,sizeof(int),1,f,swap);

	      fread_sw(&dummy,sizeof(int),1,f,swap);
	      info->s_f_vertexIndex[i]=sizeof(NDNET_UINT)*((size_t)net->f_numVertexIndexCum[i][net->nfaces[i]]);
	      info->rss_f_vertexIndex[i]=sizeof(NDNET_UINT);
	      info->rn_f_vertexIndex[i]=((size_t)(i+1)*net->nfaces[i]);
	      info->rs_f_vertexIndex[i]=net->indexSize;
	      fgetpos(f,&info->f_vertexIndex[i]);
	      fseek(f,info->rn_f_vertexIndex[i]*info->rs_f_vertexIndex[i],SEEK_CUR);
	      net->f_vertexIndex[i]=NULL;
	      // here
	      fread_sw(&dummy,sizeof(int),1,f,swap);
	    }
	  else
	    {
	      fread_sw(&dummy,sizeof(int),1,f,swap);
		
	      info->s_f_vertexIndex[i]=sizeof(NDNET_UINT)*((size_t)(i+1)*net->nfaces[i]);
	      info->rss_f_vertexIndex[i]=sizeof(NDNET_UINT);
	      info->rn_f_vertexIndex[i]=((size_t)(i+1)*net->nfaces[i]);
	      info->rs_f_vertexIndex[i]=net->indexSize;
	      fgetpos(f,&info->f_vertexIndex[i]);
	      
	      fseek(f,info->rn_f_vertexIndex[i]*info->rs_f_vertexIndex[i],SEEK_CUR);
	      net->f_vertexIndex[i]=NULL;
	      // here
	      fread_sw(&dummy,sizeof(int),1,f,swap);
	    }
	}
    }
 
  j=(1+net->ndims)*sizeof(int);
  net->haveFaceFromVertex=(int*)calloc((size_t)net->ndims+1,sizeof(int));
  fread_sw(&dummy,sizeof(int),1,f,swap);
  fread_sw(net->haveFaceFromVertex,sizeof(int)*((size_t)net->ndims+1),1,f,swap);
  fread_sw(&dummy,sizeof(int),1,f,swap);

  net->v_faceIndex = calloc(1+net->ndims,sizeof(NDNET_UINT*));
  net->v_numFaceIndexCum = calloc(1+net->ndims,sizeof(NDNET_IDCUMT*));

  for (i=0;i<1+net->ndims;i++)
    {     
      if (net->haveFaceFromVertex[i])
	{
	  fread_sw(&dummy,sizeof(int),1,f,swap);
	  info->s_v_numFaceIndexCum[i]=sizeof(NDNET_IDCUMT)*((size_t)net->nvertex+1);
	  info->rss_v_numFaceIndexCum[i]=sizeof(NDNET_IDCUMT);
	  info->rn_v_numFaceIndexCum[i]=((size_t)net->nvertex+1);
	  info->rs_v_numFaceIndexCum[i]=net->cumIndexSize;
	  fgetpos(f,&info->v_numFaceIndexCum[i]);
	  fseek(f,info->rn_v_numFaceIndexCum[i]*info->rs_v_numFaceIndexCum[i],SEEK_CUR);
	  net->v_numFaceIndexCum[i]=NULL;
	  //here
	  fread_sw(&dummy,sizeof(int),1,f,swap);

	  fread_sw(&dummy,sizeof(int),1,f,swap);
	  info->s_v_faceIndex[i]=sizeof(NDNET_UINT)*((size_t)net->v_numFaceIndexCum[i][net->nvertex]);
	  info->rss_v_faceIndex[i]=sizeof(NDNET_UINT);
	  info->rn_v_faceIndex[i]=((size_t)net->v_numFaceIndexCum[i][net->nvertex]);
	  info->rs_v_faceIndex[i]=net->indexSize;
	  fgetpos(f,&info->v_faceIndex[i]);	  
	  fseek(f,info->rn_v_faceIndex[i]*info->rs_v_faceIndex[i],SEEK_CUR);
	  net->v_faceIndex[i]=NULL;
	  //here
	  fread_sw(&dummy,sizeof(int),1,f,swap);
	}
    }

  net->haveFaceFromFace =calloc(1+net->ndims,sizeof(int *));
  net->f_faceIndex = calloc(1+net->ndims,sizeof(NDNET_UINT**));
  net->f_numFaceIndexCum = calloc((1+net->ndims),sizeof(NDNET_IDCUMT**));

  fread_sw(&dummy,sizeof(int),1,f,swap);
  for (k=0;k<net->ndims+1;k++)
    {
      net->f_faceIndex[k] = (NDNET_UINT**)calloc((1+net->ndims),sizeof(NDNET_UINT*));
      net->f_numFaceIndexCum[k] = (NDNET_IDCUMT**)calloc((1+net->ndims),sizeof(NDNET_IDCUMT*));
	
      net->haveFaceFromFace[k]=(int*)calloc(((size_t)net->ndims+1),sizeof(int));
      fread_sw(net->haveFaceFromFace[k],sizeof(int),((size_t)net->ndims+1),f,swap);
    }   
  fread_sw(&dummy,sizeof(int),1,f,swap);

  for (i=0;i<1+net->ndims;i++)
    {
      for (k=0;k<1+net->ndims;k++)
	{
	  
	  if (net->haveFaceFromFace[i][k])
	    {
	      fread_sw(&dummy,sizeof(int),1,f,swap);
	      info->s_f_numFaceIndexCum[i][k]=sizeof(NDNET_IDCUMT)*((size_t)net->nfaces[i]+1);
	      info->rss_f_numFaceIndexCum[i][k]=sizeof(NDNET_IDCUMT);
	      info->rn_f_numFaceIndexCum[i][k]=((size_t)net->nfaces[i]+1);
	      info->rs_f_numFaceIndexCum[i][k]=net->cumIndexSize;
	      fgetpos(f,&info->f_numFaceIndexCum[i][k]);
	      fseek(f,info->rn_f_numFaceIndexCum[i][k]*info->rs_f_numFaceIndexCum[i][k],SEEK_CUR);
	      net->f_numFaceIndexCum[i][k]=NULL;
	      //here
	      fread_sw(&dummy,sizeof(int),1,f,swap);
	      
	      fread_sw(&dummy,sizeof(int),1,f,swap);
	      info->s_f_faceIndex[i][k]=sizeof(NDNET_UINT)*((size_t)net->f_numFaceIndexCum[i][k][net->nfaces[i]]);
	      info->rss_f_faceIndex[i][k]=sizeof(NDNET_UINT);
	      info->rn_f_faceIndex[i][k]=((size_t)net->f_numFaceIndexCum[i][k][net->nfaces[i]]);
	      info->rs_f_faceIndex[i][k]=net->indexSize;
	      fgetpos(f,&info->f_faceIndex[i][k]);
	      fseek(f,info->rn_f_faceIndex[i][k]*info->rs_f_faceIndex[i][k],SEEK_CUR);
	      net->f_faceIndex[i][k]=NULL;
	      //here
	      fread_sw(&dummy,sizeof(int),1,f,swap);	  
	    }
	}
    }
   
  j=sizeof(int);
  fread_sw(&dummy,sizeof(int),1,f,swap);
  fread_sw(&net->haveVFlags,sizeof(int),1,f,swap);
  fread_sw(&dummy,sizeof(int),1,f,swap);

  if (net->haveVFlags)
    {
      fread_sw(&dummy,sizeof(int),1,f,swap);
      info->s_v_flag=sizeof(unsigned char)*(size_t)net->nvertex;
      info->rs_v_flag=sizeof(unsigned char);
      info->rn_v_flag=(size_t)net->nvertex;
      info->rss_v_flag=sizeof(unsigned char);
      fgetpos(f,&info->v_flag);
      fseek(f,info->rn_v_flag*info->rs_v_flag,SEEK_CUR);
      net->v_flag=NULL;
      //here
      fread_sw(&dummy,sizeof(int),1,f,swap);
    }

  net->haveFFlags=(int*)calloc(net->ndims+1,sizeof(int));
  fread_sw(&dummy,sizeof(int),1,f,swap);
  fread_sw(net->haveFFlags,sizeof(int)*(net->ndims+1),1,f,swap);
  fread_sw(&dummy,sizeof(int),1,f,swap);
  net->f_flag=calloc(net->ndims+1,sizeof(unsigned char *));

  for (i=0;i<=net->ndims;i++)
    if (net->haveFFlags[i])
      {	  
	j=sizeof(unsigned char)*net->nfaces[i];
	fread_sw(&dummy,sizeof(int),1,f,swap);
	
	if (j) 
	  {
	    info->s_f_flag[i]=sizeof(unsigned char)*(size_t)net->nfaces[i];
	    info->rs_f_flag[i]=sizeof(unsigned char);
	    info->rn_f_flag[i]=(size_t)net->nfaces[i];
	    info->rss_f_flag[i]=sizeof(unsigned char);
	    fgetpos(f,&info->f_flag[i]);
	    fseek(f,info->rn_f_flag[i]*info->rs_f_flag[i],SEEK_CUR);
	    net->f_flag[i]=NULL;
	    //here
	  }
	fread_sw(&dummy,sizeof(int),1,f,swap);
      }

  fread_sw(&dummy,sizeof(int),1,f,swap);
  fread_sw(&net->ndata,sizeof(int),1,f,swap);
  fread_sw(&dummy,sizeof(int),1,f,swap);
  //printf("found ndata = %d\n",net->ndata);
  net->data=(NDnetwork_Data*)calloc(net->ndata,sizeof(NDnetwork_Data));
  for (i=0;i<net->ndata;i++)
    {     
      fread_sw(&dummy,sizeof(int),1,f,swap);
      fread_sw(&(net->data[i].type),sizeof(int),1,f,swap);
      fread_sw(net->data[i].name,sizeof(char)*255,1,f,swap);
      fread_sw(&dummy,sizeof(int),1,f,swap);

      if (net->data[i].type==0)
	{
	  fread_sw(&dummy,sizeof(int),1,f,swap);
	  info->s_data[i]=sizeof(double)*(size_t)net->nvertex;
	  info->rs_data[i]=sizeof(double);
	  info->rn_data[i]=(size_t)net->nvertex;
	  info->rss_data[i]=sizeof(double);
	  fgetpos(f,&info->data[i]);
	  fseek(f,info->rn_data[i]*info->rs_data[i],SEEK_CUR);
	  net->data[i].data=NULL;
	  // here
	  fread_sw(&dummy,sizeof(int),1,f,swap);
	}
      else
	{
	  fread_sw(&dummy,sizeof(int),1,f,swap);
	  info->s_data[i]=sizeof(double)*(size_t)net->nfaces[net->data[i].type];
	  info->rs_data[i]=sizeof(double);
	  info->rn_data[i]=(size_t)net->nfaces[net->data[i].type];
	  info->rss_data[i]=sizeof(double);
	  fgetpos(f,&info->data[i]);
	  fseek(f,info->rn_data[i]*info->rs_data[i],SEEK_CUR);
	  net->data[i].data=NULL;
	  //here
	  fread_sw(&dummy,sizeof(int),1,f,swap);
	}  
    }

  net->supData=(NDnetwork_SupData*)calloc(net->nsupData,sizeof(NDnetwork_SupData));
  for (i=0;i<net->nsupData;i++)
    {     
      j=sizeof(int)+255*sizeof(char);
      fread_sw(&dummy,sizeof(int),1,f,swap);
      fread_sw(&(net->supData[i].type),sizeof(int),1,f,swap);
      fread_sw(net->supData[i].name,sizeof(char)*255,1,f,swap);
      fread_sw(&(net->supData[i].datasize),sizeof(int),1,f,swap);
      fread_sw(net->supData[i].datatype,sizeof(char)*255,1,f,swap);
      fread_sw(&dummy,sizeof(int),1,f,swap);

      if (net->supData[i].type==0)
	{
	  fread_sw(&dummy,sizeof(int),1,f,swap);
	  info->s_supData[i]=(size_t)net->supData[i].datasize*(size_t)net->nvertex;
	  info->rs_supData[i]=(size_t)net->supData[i].datasize;
	  info->rn_supData[i]=(size_t)net->nvertex;
	  info->rss_supData[i]=(size_t)net->supData[i].datasize;
	  fgetpos(f,&info->supData[i]);
	  fseek(f,info->rn_supData[i]*info->rs_supData[i],SEEK_CUR);
	  net->supData[i].data=NULL;
	  // here
	  fread_sw(&dummy,sizeof(int),1,f,swap);
	}
      else
	{
	  fread_sw(&dummy,sizeof(int),1,f,swap);
	  info->s_supData[i]=(size_t)net->supData[i].datasize*(size_t)net->nfaces[net->supData[i].type];
	  info->rs_supData[i]=(size_t)net->supData[i].datasize;
	  info->rn_supData[i]=(size_t)net->nfaces[net->supData[i].type];
	  info->rss_supData[i]=(size_t)net->supData[i].datasize;
	  fgetpos(f,&info->supData[i]);
	  fseek(f,info->rn_supData[i]*info->rs_supData[i],SEEK_CUR);
	  net->supData[i].data=NULL;
	  // here
	  fread_sw(&dummy,sizeof(int),1,f,swap);
	}  
    }
  fclose(f);
  return net;
}

void NDnetReadData(void *res,FILE *f,fpos_t *where,size_t rs, size_t rn, size_t rss, int swap)
{
  fsetpos(f,where);
  fread_sw_ului(res,rs,rn,rss,f,swap);
}

void NDnetWriteData(void *data,FILE *f,fpos_t *where,size_t rss, size_t rn,long count)
  {
    long n=count;
    
    if (n<0) n=rn;
        
    fsetpos(f,where);
    fwrite(data,rss,n,f);
    fgetpos(f,where);
  }

NDnet_fstruct_w_info *Save_NDnetwork_header(NDnetwork *net,const char *filename, FILE **fout)
{
  NDnet_fstruct_w_info *info=(NDnet_fstruct_w_info *)malloc(sizeof(NDnet_fstruct_w_info));
  int j;
  long i,k;
  char tag[NDNETWORK_DATA_STR_SIZE];
  FILE *f;
   
  memset(tag,0,NDNETWORK_DATA_STR_SIZE*sizeof(char));
  strcpy(tag,NDNETWORK_TAG);
  i=NDNETWORK_DATA_STR_SIZE;
  
  f=(*fout)=fopen(filename,"w");
  fwrite(&i,sizeof(int),1,f);
  fwrite(tag,sizeof(char),NDNETWORK_DATA_STR_SIZE,f);
  fwrite(&i,sizeof(int),1,f);

  j=sizeof(int)*2;
  fwrite(&j,sizeof(int),1,f);
  fwrite(&net->ndims,sizeof(int),1,f);
  fwrite(&net->ndims_net,sizeof(int),1,f);
  fwrite(&j,sizeof(int),1,f);  
  j=sizeof(NDNET_UINT)+sizeof(int)*4+80*sizeof(char)+net->ndims*sizeof(double)*2+152*sizeof(char);
  fwrite(&j,sizeof(int),1,f);
  fwrite(net->comment,sizeof(char),80,f);
  fwrite(&net->periodicity,sizeof(int),1,f);
  fwrite(&net->isSimpComplex,sizeof(int),1,f);
  fwrite(net->x0,sizeof(double),net->ndims,f);
  fwrite(net->delta,sizeof(double),net->ndims,f);
  fwrite(&net->indexSize,sizeof(int),1,f);
  fwrite(&net->cumIndexSize,sizeof(int),1,f);
  fwrite(net->dummy,sizeof(char),160-8,f);
  fwrite(&net->nvertex,sizeof(NDNET_UINT),1,f);
  fwrite(&j,sizeof(int),1,f);
  
  j=net->ndims*net->nvertex;
  fwrite(&j,sizeof(int),1,f);
  fgetpos(f,&info->v_coord);
  info->rn_v_coord=(size_t)net->ndims*(size_t)net->nvertex;
  info->rss_v_coord=sizeof(float);
  fwrite_dummy(sizeof(float),(size_t)net->ndims*(size_t)net->nvertex,f);
  fwrite(&j,sizeof(int),1,f);
  
  j=(1+net->ndims)*sizeof(NDNET_UINT);
  fwrite(&j,sizeof(int),1,f);
  fwrite(net->nfaces,sizeof(NDNET_UINT)*((size_t)net->ndims+1),1,f);
  fwrite(&j,sizeof(int),1,f);

  j=(1+net->ndims)*sizeof(int);
  fwrite(&j,sizeof(int),1,f);
  fwrite(net->haveVertexFromFace,sizeof(int)*((size_t)net->ndims+1),1,f);
  fwrite(&j,sizeof(int),1,f);

  
  for (i=0;i<1+net->ndims;i++)
  {
      
      if (net->haveVertexFromFace[i])
      {
	  if (!net->isSimpComplex)
	  {
	      j=sizeof(NDNET_IDCUMT)*((size_t)net->nfaces[i]+1);
	      fwrite(&j,sizeof(int),1,f);
	      fgetpos(f,&info->f_numVertexIndexCum[i]);
	      info->rn_f_numVertexIndexCum[i]=((size_t)net->nfaces[i]+1);
	      info->rss_f_numVertexIndexCum[i]=sizeof(NDNET_IDCUMT);
	      fwrite_dummy(sizeof(NDNET_IDCUMT),((size_t)net->nfaces[i]+1),f);
	      fwrite(&j,sizeof(int),1,f);

	      j=sizeof(NDNET_UINT)*((size_t)net->f_numVertexIndexCum[i][net->nfaces[i]]);
	      fwrite(&j,sizeof(int),1,f);
	      fgetpos(f,&info->f_vertexIndex[i]);
	      info->rss_f_vertexIndex[i]=sizeof(NDNET_UINT);
	      info->rn_f_vertexIndex[i]=((size_t)(i+1)*net->nfaces[i]);
	      fwrite_dummy(sizeof(NDNET_UINT),((size_t)net->f_numVertexIndexCum[i][net->nfaces[i]]),f);
	      fwrite(&j,sizeof(int),1,f);
	  }
	  else
	  {
	      j=sizeof(NDNET_UINT)*((size_t)(i+1)*net->nfaces[i]);
	      fwrite(&j,sizeof(int),1,f);
	      fgetpos(f,&info->f_vertexIndex[i]);
	      info->rn_f_vertexIndex[i]=((size_t)(i+1)*net->nfaces[i]);
	      info->rss_f_vertexIndex[i]=sizeof(NDNET_UINT);
	      fwrite_dummy(sizeof(NDNET_UINT),((size_t)(i+1)*net->nfaces[i]),f);
	      fwrite(&j,sizeof(int),1,f);
	  }
	  
      }
  }
  
  j=(1+net->ndims)*sizeof(int);
  fwrite(&j,sizeof(int),1,f);
  fwrite(net->haveFaceFromVertex,sizeof(int)*((size_t)net->ndims+1),1,f);
  fwrite(&j,sizeof(int),1,f);

  for (i=0;i<1+net->ndims;i++)
  {
      
      if (net->haveFaceFromVertex[i])
      {
	  j=sizeof(NDNET_IDCUMT)*((size_t)net->nvertex+1);
	  fwrite(&j,sizeof(int),1,f);
	  fgetpos(f,&info->v_numFaceIndexCum[i]);
	  info->rn_v_numFaceIndexCum[i]=((size_t)net->nvertex+1);
	  info->rss_v_numFaceIndexCum[i]=sizeof(NDNET_IDCUMT);
	  fwrite_dummy(sizeof(NDNET_IDCUMT),((size_t)net->nvertex+1),f);
	  fwrite(&j,sizeof(int),1,f);
	  
	  j=sizeof(NDNET_UINT)*((size_t)net->v_numFaceIndexCum[i][net->nvertex]);
	  fwrite(&j,sizeof(int),1,f);
	  fgetpos(f,&info->v_faceIndex[i]);
	  info->rss_v_faceIndex[i]=sizeof(NDNET_UINT);
	  info->rn_v_faceIndex[i]=((size_t)net->v_numFaceIndexCum[i][net->nvertex]);
	  fwrite_dummy(sizeof(NDNET_UINT),((size_t)net->v_numFaceIndexCum[i][net->nvertex]),f);
	  fwrite(&j,sizeof(int),1,f);	  
      }
  }
  
  j=(1+net->ndims)*(1+net->ndims)*sizeof(int);
  fwrite(&j,sizeof(int),1,f);
  for (k=0;k<net->ndims+1;k++)
      fwrite(net->haveFaceFromFace[k],sizeof(int)*((size_t)net->ndims+1),1,f);
  fwrite(&j,sizeof(int),1,f);

  for (i=0;i<1+net->ndims;i++)
  {
      for (k=0;k<1+net->ndims;k++)
      {

	  if (net->haveFaceFromFace[i][k])
	  {
	      j=sizeof(NDNET_IDCUMT)*(1+(size_t)net->nfaces[i]);
	      fwrite(&j,sizeof(int),1,f);
	      info->rss_f_numFaceIndexCum[i][k]=sizeof(NDNET_IDCUMT);
	      info->rn_f_numFaceIndexCum[i][k]=((size_t)net->nfaces[i]+1);
	      fgetpos(f,&info->f_numFaceIndexCum[i][k]);
	      fwrite_dummy(sizeof(NDNET_IDCUMT),((size_t)net->nfaces[i]+1),f);
	      fwrite(&j,sizeof(int),1,f);
	      
	      j=sizeof(NDNET_UINT)*((size_t)net->f_numFaceIndexCum[i][k][net->nfaces[i]]);
	      fwrite(&j,sizeof(int),1,f);
	      info->rss_f_faceIndex[i][k]=sizeof(NDNET_UINT);
	      info->rn_f_faceIndex[i][k]=((size_t)net->f_numFaceIndexCum[i][k][net->nfaces[i]]);
	      fgetpos(f,&info->f_faceIndex[i][k]);
	      fwrite_dummy(sizeof(NDNET_UINT),((size_t)net->f_numFaceIndexCum[i][k][net->nfaces[i]]),f);
	      fwrite(&j,sizeof(int),1,f);	  
	  }
      }
  }

  j=sizeof(int);
  fwrite(&j,sizeof(int),1,f);
  fwrite(&net->haveVFlags,sizeof(int),1,f);
  fwrite(&j,sizeof(int),1,f);

  if (net->haveVFlags)
  {
      j=sizeof(unsigned char)*net->nvertex;
      fwrite(&j,sizeof(int),1,f);
      fgetpos(f,&info->v_flag);
      info->rn_v_flag=(size_t)net->nvertex;
      info->rss_v_flag=sizeof(unsigned char);
      fwrite_dummy(sizeof(unsigned char),(size_t)net->nvertex,f);
      fwrite(&j,sizeof(int),1,f);
  }

  
  j=sizeof(int)*(net->ndims+1);
  fwrite(&j,sizeof(int),1,f);
  fwrite(net->haveFFlags,sizeof(int)*(net->ndims+1),1,f);
  fwrite(&j,sizeof(int),1,f);
  
  for (i=0;i<1+net->ndims;i++)
      if (net->haveFFlags[i])
      {
	  j=sizeof(unsigned char)*net->nfaces[i];
	  fwrite(&j,sizeof(int),1,f);
	  if (j) 
	    {
	      fgetpos(f,&info->f_flag[i]);
	      info->rn_f_flag[i]=(size_t)net->nfaces[i];
	      info->rss_f_flag[i]=sizeof(unsigned char);
	      fwrite_dummy(sizeof(unsigned char),(size_t)net->nfaces[i],f);
	    }
	  fwrite(&j,sizeof(int),1,f);
      }

  j=sizeof(int);
  fwrite(&j,sizeof(int),1,f);
  fwrite(&net->ndata,sizeof(int),1,f);
  fwrite(&j,sizeof(int),1,f);
  
  for (i=0;i<net->ndata;i++)
  {
      
      j=sizeof(int)+255*sizeof(char);
      fwrite(&j,sizeof(int),1,f);
      fwrite(&(net->data[i].type),sizeof(int),1,f);
      fwrite(net->data[i].name,sizeof(char)*255,1,f);
      fwrite(&j,sizeof(int),1,f);

      if (net->data[i].type==0)
      {
	j=sizeof(double)*net->nvertex;
	fwrite(&j,sizeof(int),1,f);
	fgetpos(f,&info->data[i]);
	info->rn_data[i]=(size_t)net->nvertex;
	info->rss_data[i]=sizeof(double);
	fwrite_dummy(sizeof(double),(size_t)net->nvertex,f);
	fwrite(&j,sizeof(int),1,f);
      }
      else
      {
	j=sizeof(double)*net->nfaces[net->data[i].type];
	fwrite(&j,sizeof(int),1,f);
	fgetpos(f,&info->data[i]);
	info->rn_data[i]=(size_t)net->nfaces[net->data[i].type];
	info->rss_data[i]=sizeof(double);
	fwrite_dummy(sizeof(double),(size_t)net->nfaces[net->data[i].type],f);
	fwrite(&j,sizeof(int),1,f);
      }  
  }

  j=sizeof(int);
  fwrite(&j,sizeof(int),1,f);
  fwrite(&net->nsupData,sizeof(int),1,f);
  fwrite(&j,sizeof(int),1,f);
  
  for (i=0;i<net->nsupData;i++)
  {
      
      j=2*sizeof(int)+2*255*sizeof(char);
      fwrite(&j,sizeof(int),1,f);
      fwrite(&(net->supData[i].type),sizeof(int),1,f);
      fwrite(net->supData[i].name,sizeof(char)*255,1,f);
      fwrite(&(net->supData[i].datasize),sizeof(int),1,f);
      fwrite(net->supData[i].datatype,sizeof(char)*255,1,f);
      fwrite(&j,sizeof(int),1,f);

      if (net->supData[i].type==0)
      {
	  j=(size_t)net->supData[i].datasize*net->nvertex;
	  fwrite(&j,sizeof(int),1,f);
	  fgetpos(f,&info->supData[i]);
	  info->rn_supData[i]=(size_t)net->nvertex;
	  info->rss_supData[i]=(size_t)net->supData[i].datasize;
	  fwrite_dummy((size_t)net->supData[i].datasize,(size_t)net->nvertex,f);
	  fwrite(&j,sizeof(int),1,f);
      }
      else
      {
	  j=(size_t)net->supData[i].datasize*net->nfaces[net->supData[i].type];
	  fwrite(&j,sizeof(int),1,f);
	  fgetpos(f,&info->supData[i]);
	  info->rn_supData[i]=(size_t)net->nfaces[net->supData[i].type];
	  info->rss_supData[i]=(size_t)net->supData[i].datasize;
	  fwrite_dummy((size_t)net->supData[i].datasize,(size_t)net->nfaces[net->supData[i].type],f);
	  fwrite(&j,sizeof(int),1,f);
      }  
  }

  //fclose(f);
  
  return info;
}
