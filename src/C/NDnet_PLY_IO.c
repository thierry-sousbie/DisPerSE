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

#include "NDnet_PLY_IO.h"
#include "NDnetwork.h"
#include "NDnetwork_tags.h"
#include "rply/rply.h"
#include "global.h"

typedef struct cb_bbox_data_str {
  NDnetwork **net;
  long nx0;
  long ndelta;
} cb_bbox_data;

static int cb_read_bbox(p_ply_argument argument) {
  p_ply_property property=NULL;
  long value_index;  
  cb_bbox_data *data=NULL;   
  const char *name;
  ply_get_argument_user_data(argument, (void**)&data, NULL); 
  NDnetwork *net = *(data->net);
  
  ply_get_argument_property(argument, &property, NULL,&value_index);  
  ply_get_property_info(property, &name, NULL,NULL,NULL);
  
  if (value_index<0)
    {
      if (net->ndims!=(double)ply_get_argument_value(argument))
	fprintf(stderr,"WARNING: incorrect numberof of arguments for bounding box!\n");
      return 1;
    }
  
  if ((!strcmp(name,"x0"))||(!strcmp(name,"X0")))
    net->x0[data->nx0++]=(double)ply_get_argument_value(argument);
  else if ((!strcmp(name,"delta"))||(!strcmp(name,"DELTA")))
    net->delta[data->ndelta++]=(double)ply_get_argument_value(argument);
  //else fprintf(stderr,"");
  
  return 1;
}

typedef struct cb_coord_data_str {
  NDnetwork **net;
  NDnetwork_Data *data;
  long cur;
  int flags;
} cb_coord_data;

static int cb_read_coord(p_ply_argument argument) {
  cb_coord_data *data=NULL;  
  long dec;

  ply_get_argument_user_data(argument, (void**)&data, &dec); 
  NDnetwork *net = *(data->net);
  if (dec>0)
    {
      net->v_coord[data->cur*net->ndims + dec-1]=(float)ply_get_argument_value(argument);
      data->cur++;
    }
  else 
    {
      if (data->flags)
	net->v_flag[data->cur++]=(unsigned char)ply_get_argument_value(argument);
      else
	data->data->data[data->cur++]=(double)ply_get_argument_value(argument);
    }
  
  return 1;
}

typedef struct cb_face_data_str {
  NDnetwork **net;
  NDnetwork_Data *data;
  int type;
  long cur;
  long ntot;
  //int flags;
} cb_face_data;

static int cb_read_face(p_ply_argument argument) {
  long length, value_index;
  long dec;
 
  cb_face_data *data=NULL;// = (cb_face_data *)dummy;
  ply_get_argument_user_data(argument, (void**)&data, &dec);
  
  NDnetwork *net = *(data->net);

  ply_get_argument_property(argument, NULL, &length, &value_index);
  
  if (dec<0)
    {
      if (value_index<0)
	{
	  if (data->cur==0) 
	    {
	      data->type = (int)ply_get_argument_value(argument)-1;
	      net->f_vertexIndex[data->type] = (NDNET_UINT*)calloc(data->ntot*(data->type+1),sizeof(NDNET_UINT));   
	      net->haveVertexFromFace[data->type]=1;
	    }
	  net->nfaces[data->type]++;
	}
      else 
	{
	  net->f_vertexIndex[data->type][data->cur++]=
	    (NDNET_UINT)ply_get_argument_value(argument);
	}
    }
  else
    {
      data->data->data[data->cur++]=(double)ply_get_argument_value(argument);
    }

  return 1;
}

void dummy_error_cb(p_ply ply, const char *message)
{
  return;
}

int IsPLYNetwork(const char *fname)
{
  p_ply ply;
  ply = ply_open(fname, dummy_error_cb, 0, NULL);
  if (!ply) return 0; 
  if (!ply_read_header(ply)) return 0; 
  if (!ply_close(ply)) return 0;
  return 1;
}

NDnetwork *Load_NDnetworkFromPLY(const char *fname)
{
  p_ply iply;
  NDnetwork *net=NULL;
  long nvertex=-1;
  long ndims=0;
  long ndata=0;
  NDnetwork_Data data[100];
  int *actual_type[100];

  long nf=0;
  long nc=0;
  cb_face_data face_data[100];
  cb_coord_data coord_data[100];
  cb_bbox_data bbox_data;
  int bbox_read=0;
  long i;
  int haveVFlags=0;

  bbox_data.net=&net;
  bbox_data.nx0=0;
  bbox_data.ndelta=0;
  for (i=0;i<100;i++)
    {      
      coord_data[i].net=&net;
      coord_data[i].cur=0;
      coord_data[i].flags=0;
      coord_data[i].data=NULL;

      face_data[i].net=&net;
      face_data[i].cur=0;
      //face_data[i].flags=0;
      face_data[i].ntot=-1;

      actual_type[i]=NULL;
    }
  
  iply = ply_open(fname, NULL, 0, NULL);
  if (!iply) return NULL; 
  if (!ply_read_header(iply)) return NULL;

  p_ply_element element = NULL;

  while ((element = ply_get_next_element(iply, element))) {
    p_ply_property property = NULL;
    const char *element_name;
    long ninstances = 0;
    ply_get_element_info(element, &element_name, &ninstances);

    if (!strcmp(element_name,"vertex"))
      {
	while ((property = ply_get_next_property(element, property))) 
	  {
	    const char *property_name;
	    e_ply_type type, length_type, value_type;
	    long index=0;
	    ply_get_property_info(property, &property_name, &type,&length_type, &value_type);
	    if ((!strcmp(property_name,"x"))||(!strcmp(property_name,"X"))) index=1;
	    if ((!strcmp(property_name,"y"))||(!strcmp(property_name,"Y"))) index=2;
	    if ((!strcmp(property_name,"z"))||(!strcmp(property_name,"Z"))) index=3;
	    if ((strlen(property_name)==strlen("x0"))&&
		((property_name[0]=='X')||(property_name[0]=='x')))
	      {
		index = (long)property_name[1]-(long)'0';
		if ((index>=0)&&(index<10)) index++;
		else index=0;
	      }

	    if (index>0)
	      {
		long nread=ply_set_read_cb(iply,element_name,property_name,cb_read_coord,&coord_data[nc++],index);
		ndims++;
		if (nvertex<0) nvertex=nread;
		else if (nread!=nvertex) {fprintf(stderr,"ERROR reading PLY file, wrong format!\n");exit(-1);}
	      }
	    else
	      {
		long nread=ply_set_read_cb(iply,element_name,property_name,cb_read_coord,&coord_data[nc],-ndata);

		if (!strcmp(property_name,"flags")) 
		  {
		    haveVFlags=1;
		    coord_data[nc].flags=1;
		  }
		else data[ndata].data=(double *)calloc(nread,sizeof(double));

		coord_data[nc++].data=&data[ndata];
		data[ndata].type=0;
		strcpy(data[ndata].name,property_name);
		
		ndata++;
	      }
	  }		
      }

    else if ((!strcmp(element_name,"face"))
	     ||(!strcmp(element_name,"segment"))
	     ||(!strcmp(element_name,"tetrahedron"))
	     ||((strstr(element_name,"-face")!=NULL)&&(strlen(element_name)==strlen("-face")+1))
	     ||((strstr(element_name,"-simplex")!=NULL)&&(strlen(element_name)==strlen("-simplex")+1)))	      
      {	
	int ok=0;	
	long nstart=ndata;
	int *dataTypePtr=NULL;
	while ((property = ply_get_next_property(element, property))) 
	  {
	    const char *property_name;
	    e_ply_type type, length_type, value_type;
	    
	    ply_get_property_info(property, &property_name, &type, &length_type, &value_type);
	    if (type != PLY_LIST)
	      {
		long nread=ply_set_read_cb(iply,element_name,property_name,cb_read_face,&face_data[nf],ndata);
		face_data[nf].ntot=nread;
		face_data[nf].data=&data[ndata];
		data[ndata].data=(double *)calloc(nread,sizeof(double));
		strcpy(data[ndata].name,property_name);
		
		nf++;ndata++;				
	      }
	    else
	      {
		if (ok)
		  {
		    fprintf(stderr,"WARNING: reading PLY file '%s':\n",fname);
		    fprintf(stderr,"         '%s' is supposed to contain exactly one LIST element.\n",element_name);
		  }
		long nread=ply_set_read_cb(iply,element_name,property_name,cb_read_face,&face_data[nf],-1);
		dataTypePtr=&face_data[nf].type;
		face_data[nf++].ntot=nread;
		ok=1;
	      }
	  }
	
	for (i=nstart;i<ndata;i++) actual_type[i]=dataTypePtr;
      }
    else if ((!strcmp(element_name,"bbox"))||(!strcmp(element_name,"BBOX")))
      {		
	bbox_read=1;
	while ((property = ply_get_next_property(element, property))) 
	  {
	    const char *property_name;
	    ply_get_property_info(property, &property_name, NULL,NULL,NULL);
	    ply_set_read_cb(iply,element_name,property_name,cb_read_bbox,&bbox_data,0);
	  }
      }
  }

  net=CreateNetwork(ndims,nvertex,haveVFlags); 
  net->data=(NDnetwork_Data*) calloc(ndata,sizeof(NDnetwork_Data));
  long j=0;
  for (i=0;i<ndata;i++) 
    {      
      if (strcmp(data[i].name,"flags"))
	memcpy(&net->data[j++],&data[i],sizeof(NDnetwork_Data));
    }
  net->ndata=j;
  net->data=(NDnetwork_Data*) realloc(net->data,net->ndata*sizeof(NDnetwork_Data));

  for (i=0;i<net->ndims;i++) net->x0[i]=net->delta[i]=0;
  if (!ply_read(iply)) return NULL;

  j=0;
  for (i=0;i<ndata;i++) 
    {      
     if (actual_type[i]!=NULL) data[i].type=*actual_type[i];
     if (strcmp(data[i].name,"flags"))
       {
	 if (actual_type[i]!=NULL) net->data[j].type=*actual_type[i];
	 j++;
       }
    }
  for (i=0;i<ndata;i++) 
    {   
      if (!strcmp(data[i].name,"flags"))
	{
	  if (data[i].type==0) continue;
	  
	  net->haveFFlags[data[i].type]=1;
	  net->f_flag[data[i].type]=calloc(net->nfaces[data[i].type],sizeof(unsigned char));
	  long k;
	  for (k=0;k<net->nfaces[data[i].type];k++)
	    {
	      net->f_flag[data[i].type][k]=(unsigned char)data[i].data[k];
	      //printf("%ld: %ld\n",k,(long)net->f_flag[data[i].type][k]);
	    }

	  free(data[i].data);
	  data[i].data=NULL;
	}
    }
  //for (i=0;i<net->ndata;i++) 
  //if (actual_type[i]!=NULL) net->data[i].type=*actual_type[i];

  if (!ply_close(iply)) fprintf(stderr,"ERROR closing PLY file '%s'.\n",fname);
  if (!bbox_read)
    {
      long j;
      for (i=0;i<net->ndims;i++) 
	net->x0[i]=net->delta[i]=net->v_coord[i];
	
      for (i=0;i<net->nvertex;i++)
	for (j=0;j<net->ndims;j++)
	  {
	    if (net->v_coord[i*net->ndims+j]<net->x0[j]) net->x0[j]=net->v_coord[i*net->ndims+j];
	    if (net->v_coord[i*net->ndims+j]>net->delta[j]) net->delta[j]=net->v_coord[i*net->ndims+j];
	  }
      
      for (i=0;i<net->ndims;i++)
	{
	  net->delta[i]-=net->x0[i];
	  net->x0[i]-=net->delta[i]*0.05;
	  net->delta[i]*=1.1;
	}
    }
  return net;
}

int Save_NDnetworkToPLY(NDnetwork *net,const char *fname,int type)
{
  e_ply_storage_mode storage_mode = PLY_LITTLE_ENDIAN;
  p_ply oply = NULL;
  long i,j,k;
  char *coord_prop[] = {"x","y","z","x3","x4","x5","x6","x7","x8","x9"};
  int maxFT=-1;

  if (type==NDNET_PLY_ASCII) storage_mode=PLY_ASCII;
    
  for (i=0;i<=net->ndims;i++) if (net->nfaces[i]) maxFT=i;

  oply = ply_create(fname, storage_mode, NULL, 0, NULL);
  if (!oply) {fprintf(stderr,"Unable to create file '%s'", fname);return -1;}  
  
  ply_add_element(oply,"vertex", net->nvertex);
  for (j=0;j<net->ndims;j++)
    {
      ply_add_property(oply, coord_prop[j],PLY_FLOAT,0,0);
    }
  
  for (j=0;j<net->ndata;j++)
    {
      if (net->data[j].type==0)
	{
	  ply_add_property(oply, net->data[j].name, PLY_DOUBLE,0,0);
	}
    }
  
  if (net->haveVFlags) ply_add_property(oply,"flags",PLY_UCHAR,0,0);
  
  if (maxFT>=0)
    {
      //for (k=0;k<=maxFT;k++)
      for (k=maxFT;k>=0;k--)
	{
	  if (net->nfaces[k]==0) continue;
	  char element_name[255];
	  if (k==maxFT) strcpy(element_name,"face");
	  else sprintf(element_name,"%ld-face",k);
	  ply_add_element(oply,element_name, net->nfaces[k]);
	  ply_add_property(oply,"vertex_indices", PLY_LIST,PLY_UCHAR,PLY_UINT);
	  for (j=0;j<net->ndata;j++)
	    {
	      if (net->data[j].type==k)
		{
		  ply_add_property(oply, net->data[j].name, PLY_DOUBLE,0,0);
		}
	    }
	  if (net->haveFFlags[k]) ply_add_property(oply,"flags",PLY_UCHAR,0,0);
	}
    }

  ply_add_element(oply,"bbox",1);
  ply_add_property(oply,"x0",PLY_LIST,PLY_UCHAR,PLY_DOUBLE);
  ply_add_property(oply,"delta",PLY_LIST,PLY_UCHAR,PLY_DOUBLE);
  
  if (!ply_write_header(oply)) 
    {
      fprintf(stderr,"Failed writing '%s' header", fname);
      return -1;
    }

  int ndata=0;
  double *data[50];

  
  
  for (j=0;j<net->ndata;j++)
    if (net->data[j].type==0) 
	data[ndata++]=net->data[j].data;
  
  for (i=0;i<net->nvertex;i++)
    {
      for (j=0;j<net->ndims;j++) ply_write(oply,net->v_coord[i*net->ndims+j]);
      for (j=0;j<ndata;j++) ply_write(oply,data[j][i]);
      if (net->haveVFlags) ply_write(oply,net->v_flag[i]);
    }
    
  if (maxFT>=0)
    {
      //for (k=0;k<=maxFT;k++)
      for (k=maxFT;k>=0;k--)
	{
	  if (net->nfaces[k]==0) continue;
	  ndata=0;

	  for (j=0;j<net->ndata;j++)
	    if (net->data[j].type==k) 
	      data[ndata++]=net->data[j].data;
	  
	  for (j=0;j<net->nfaces[k];j++)
	    {
	      ply_write(oply,k+1);
	      for (i=0;i<=k;i++) ply_write(oply,net->f_vertexIndex[k][j*(k+1)+i]);
	      for (i=0;i<ndata;i++) ply_write(oply,data[i][j]);
	      if (net->haveFFlags[k]) ply_write(oply,net->f_flag[k][j]);
	    }
	}
    }

  ply_write(oply,net->ndims);
  for (i=0;i<net->ndims;i++) ply_write(oply,net->x0[i]);
  ply_write(oply,net->ndims);
  for (i=0;i<net->ndims;i++) ply_write(oply,net->delta[i]);
   
  if (!ply_close(oply)) {fprintf(stderr,"Error closing file '%s'", fname);return -1;}

  return 0;
}
