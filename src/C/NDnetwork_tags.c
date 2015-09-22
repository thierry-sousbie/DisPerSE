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
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "simplex.h"
#include "NDnetwork_tags.h"
#include "global.h"
#include "NDfield.h"


int NDDataIndex(NDnetwork *net, int type, const char *name)
{
    int i;

    for (i=0;i<net->ndata;i++)
      {
	//printf("'%d'/'%s' <=> '%d'/'%s'\n",type,name,net->data[i].type,net->data[i].name);
	if ((net->data[i].type==type)&&(!strcmp(name,net->data[i].name)))
	  {
	    //printf("FOUND '%d'/'%s' @%d.\n",type,name,i);
	    return i;
	  }
      }
	
    //printf("COULD NOT FIND '%d'/'%s'.\n",type,name);
    return -1;
    /*
    for (i=0;i<net->ndata;i++)
    {
      printf("'%d'/'%s' <=> '%d'/'%s'\n",type,name,net->data[i].type,net->data[i].name);
      if ((net->data[i].type==type)&&(net->data[i].name[0]==name[0]))
	if (!strcmp(name,net->data[i].name))
	  {
	    printf("FOUND\n");
	    return i;
	  }
      printf("FAILED.\n");
    }

    return -1;
    */
}

int NDSupDataIndex(NDnetwork *net, int type, const char *name)
{
    int i;
    
    for (i=0;i<net->nsupData;i++)
    {
	if ((net->supData[i].type==type)&&(net->supData[i].name[0]==name[0]))
	    if (!strcmp(name,net->supData[i].name))
		return i;
    }

    return -1;
}

NDnetwork_Data *addNDData(NDnetwork *net, int type, const char *name, 
			     double(*compute_val)(NDnetwork* net,int type,int index,void* info),
			     void *info)
{
    long i;
    int index=-1;
    double *data;
    long N;
    char tmp[255];

    if (type==0) 
	N=net->nvertex;
    else
	N=net->nfaces[type];

    if (N==0)
    {
	fprintf(stderr,"In addNDData: Cannot add data for %d-faces, as they are not defined.\n",type);
	return NULL;
    }

    index=NDDataIndex(net,type,name);
    
    if (type==0) sprintf(tmp,"vertices");
    else if (type==1) sprintf(tmp,"edges");
    else sprintf(tmp,"%d faces",type);

    if (index>=0)
    {if (verbose>1) printf ("Replacing '%s' data field for %s ... ",name,tmp);fflush(0);}
    else
    {if (verbose>1) printf ("Adding '%s' data field for %s ... ",name,tmp);fflush(0);}
    
    if (index<0)
    {
	net->ndata++;
	net->data=realloc(net->data,sizeof(NDnetwork_Data)*net->ndata);
	index=net->ndata-1;
	net->data[index].type=type;
	strcpy(net->data[index].name,name);
	net->data[index].data=malloc(N*sizeof(double));
    }
 
    data=net->data[index].data;
    
    for (i=0;i<N;i++)
	data[i] = compute_val(net,type,i,info);

    if (verbose>1) printf("done.\n");

    return &net->data[index];
}

// be carefull, Data_p is not COPIED and *supData_p is set to NULL
// use ptr =  addNDDataArr(...,&ptr,...) to keep it.
double *addNDDataArr(NDnetwork *net, int type, const char *name, double **data_p)
{
    int index=-1;
    long N;
    char tmp[255];
    void *data;

    if (data_p!=NULL)
      data=*data_p;
    else data=NULL;

    if (type==0) 
	N=net->nvertex;
    else
	N=net->nfaces[type];

    if (N==0)
    {
	fprintf(stderr,"In addNDDataArr: Cannot add data for %d-faces, as they are not defined.\n",type);
	return NULL;
    }

    index=NDDataIndex(net,type,name);
    
    if (type==0) sprintf(tmp,"vertice");
    else if (type==1) sprintf(tmp,"edges");
    else sprintf(tmp,"%d faces",type);

    if (index>=0)
    {if (verbose>1) printf ("Replacing '%s' data field for %s ... ",name,tmp);fflush(0);}
    else
    {if (verbose>1) printf ("Adding '%s' data field for %s ... ",name,tmp);fflush(0);}
    
    if (index<0)
    {
	net->ndata++;
	net->data=realloc(net->data,sizeof(NDnetwork_Data)*net->ndata);	
	index=net->ndata-1;
    }
    else
    {
	free(net->data[index].data);
    }
    
    net->data[index].type=type;
    memset(net->data[index].name,0,255);
    strcpy(net->data[index].name,name);
    if (data!=NULL)
      net->data[index].data=data;
    else net->data[index].data=calloc(N,sizeof(double));
   
    if (data_p != NULL) *data_p=NULL;

    if (verbose>1) printf("done.\n");

    return net->data[index].data;
}

// be carefull, supData_p is not COPIED and *supData_p is set to NULL
// use ptr =  addNDSupDataArr(...,&ptr,...) to keep it.
void *addNDSupDataArr(NDnetwork *net, int type, const char *name, int datasize, const char *datatype, void **supData_p)
{
    int index=-1;
    long N;
    char tmp[255];
    void *supData=*supData_p;

    if (type==0) 
	N=net->nvertex;
    else
	N=net->nfaces[type];

    if (N==0)
    {
	fprintf(stderr,"In addNDSupDataArr: Cannot add supData for %d-faces, as they are not defined.\n",type);
	return NULL;
    }

    index=NDSupDataIndex(net,type,name);
    
    if (type==0) sprintf(tmp,"vertice");
    else if (type==1) sprintf(tmp,"edges");
    else sprintf(tmp,"%d faces",type);

    if (index>=0)
    {if (verbose>1) printf ("Replacing '%s' supData field for %s ... ",name,tmp);fflush(0);}
    else
    {if (verbose>1) printf ("Adding '%s' supData field for %s ... ",name,tmp);fflush(0);}
    
    if (index<0)
    {
	net->nsupData++;
	net->supData=realloc(net->supData,sizeof(NDnetwork_SupData)*net->nsupData);
	index=net->nsupData-1;
    }
    else
    {
	free(net->supData[index].data);
    }
    
    net->supData[index].type=type;
    strcpy(net->supData[index].name,name);
    net->supData[index].datasize=datasize;
    strcpy(net->supData[index].datatype,datatype);
    net->supData[index].data=supData;
    
    *supData_p=NULL;

    if (verbose>1) printf("done.\n");

    return net->supData[index].data;
}
/*
double NDNetNPatchesTag(NDnetwork *net, int type, int index, void *info)
{
    NetPatch_type *proba = &(((NetPatch_type *)info)[index]);

    return proba->n;
}
*/

 /*
double NDNetPatchProbaTag(NDnetwork *net, int type, int index, void *info)
{
    int i,j;
    int max;
   
    if (type==0)
    {
	NetPatch_type *proba = &(((NetPatch_type *)info)[index]);

	max=0;
	for (i=0;i<proba->n;i++)
	{
	    if (proba->p[i]>proba->p[max]) max=i;
	}
	
	return (double) proba->p[max];
    }
    else
    {
	int n=NUM_VERTEX_IN_FACE(net,type,index);
	int patchProba[n];
	int patchIndex[n];
	int npatches=0;
	NetPatch_type *proba;
	uint *vlist = VERTEX_IN_FACE(net,type,index);
	    
	for (i=0;i<n;i++)
	{
	    proba = &(((NetPatch_type *)info)[vlist[i]]);
	    max=0;
	    for (j=1;j<proba->n;j++)
		if (proba->p[j]>proba->p[max]) max=j;
	    
	    for (j=0;j<npatches;j++)
		if (proba->index[max]==patchIndex[j])
		{
		    patchProba[j]+=proba->p[max];
		    break;
		}

	    if (j==npatches)
	    {
		patchProba[npatches]=proba->p[max];
		patchIndex[npatches]=proba->index[max];
		npatches++;
	    }
	}
	
	max=0;
	for (i=1;i<npatches;i++)
	    if (patchProba[i]>patchProba[max]) 
		max=i;
	
	double tmpd;

	for (i=0,tmpd=0;i<npatches;i++) tmpd+=patchProba[i];
	for (i=0,tmpd=0;i<npatches;i++) patchProba[i]/=tmpd;
	return (double)patchProba[max];
    }
}
 */

/*
// info should be set to patches proba as returned by ComputeNetPatches
double NDNetPatchTag(NDnetwork *net, int type, int index, void *info)
{
    int i,j;
    int max;
   
    if (type==0)
    {
	NetPatch_type *proba = &(((NetPatch_type *)info)[index]);

	max=0;
	for (i=0;i<proba->n;i++)
	{
	    if (proba->p[i]>proba->p[max]) max=i;
	}
	
	return (double) proba->index[max];
    }
    else
    {
	int n=NUM_VERTEX_IN_FACE(net,type,index);
	int patchProba[n];
	int patchIndex[n];
	int npatches=0;
	NetPatch_type *proba;
	uint *vlist = VERTEX_IN_FACE(net,type,index);
	    
	for (i=0;i<n;i++)
	{
	    proba = &(((NetPatch_type *)info)[vlist[i]]);
	    max=0;
	    for (j=1;j<proba->n;j++)
		if (proba->p[j]>proba->p[max]) max=j;
	    
	    for (j=0;j<npatches;j++)
		if (proba->index[max]==patchIndex[j])
		{
		    patchProba[j]+=proba->p[max];
		    break;
		}

	    if (j==npatches)
	    {
		patchProba[npatches]=proba->p[max];
		patchIndex[npatches]=proba->index[max];
		npatches++;
	    }
	}
	
	max=0;
	for (i=1;i<npatches;i++)
	    if (patchProba[i]>patchProba[max]) 
		max=i;
	
	return (double)patchIndex[max];
    }
}
*/

 /*
int TagNetworkPatches(NDnetwork *net,int type, NetPatch_type **vpatches,NetPatch_type **ppatches,int storeProba, int storePatches)
{
    if (vpatches != NULL)
    {
	addNDData(net,type,VOID_PATCH_INDEX_TAG,NDNetPatchTag,*vpatches);
	if (storeProba) addNDData(net,type,VOID_PATCH_PROBA_TAG,NDNetPatchProbaTag,*vpatches);
	//addNDData(net,net->ndims,VOID_PATCH_INDEX_TAG,NDNetPatchTag,*vpatches);
    }

    if (ppatches != NULL)
    {
	addNDData(net,type,PEAK_PATCH_INDEX_TAG,NDNetPatchTag,*ppatches);
	if (storeProba) addNDData(net,type,PEAK_PATCH_PROBA_TAG,NDNetPatchProbaTag,*ppatches);
	//addNDData(net,net->ndims,PEAK_PATCH_INDEX_TAG,NDNetPatchTag,*ppatches);
    }
    
    //NDNetTagCritical(net);

    if ((storePatches)&&(type==0))
    {
	if (vpatches != NULL) *vpatches = (NetPatch_type *)addNDSupDataArr(
	    net,type,"void patches proba", sizeof(NetPatch_type),"NetPatch_type", (void **)vpatches);

	if (ppatches != NULL) *ppatches = (NetPatch_type *)addNDSupDataArr(
	    net,type,"peak patches proba", sizeof(NetPatch_type),"NetPatch_type", (void **)ppatches);
    }

    return 1;
}


// call something like TagNetworkPatches(...) or
// addNDData(net,net->ndims,VOID_PATCH_INDEX_TAG,NDNetPatchTag,vpatches);
// addNDData(net,net->ndims,PEAK_PATCH_INDEX_TAG,NDNetPatchTag,ppatches);
// before calling this function ...
int NDNetTagCritical(NDnetwork *net)
{
    int type;
    int vpindex;
    int ppindex;
    int nfaces;
    uint *faces=NULL;
    int i,j,k,l;
    double *vp;
    double *pp;
    int *pplist=NULL;
    int *vplist=NULL;
    int nppatch=0;
    int nvpatch=0;
    int nppatch_max=0;
    int nvpatch_max=0;
    int ntagp[net->ndims_net];
    int ntagv[net->ndims_net];

    //vpindex=NDDataIndex(net, net->ndims, VOID_PATCH_INDEX_TAG);
    //ppindex=NDDataIndex(net, net->ndims, PEAK_PATCH_INDEX_TAG);
    vpindex=NDDataIndex(net, 0, VOID_PATCH_INDEX_TAG);
    ppindex=NDDataIndex(net, 0, PEAK_PATCH_INDEX_TAG);
    //printf ("\n%d %d \n",ppindex,vpindex);
    if ((vpindex<0)&&(ppindex<0))
    {
	fprintf (stderr,"ERROR in NDNetTagCritical(NDnetwork *net).\n");
	fprintf (stderr,"%d-faces should be tagged with '%s' or '%s'.\n",0,VOID_PATCH_INDEX_TAG,PEAK_PATCH_INDEX_TAG);
	return -1;
    }

    if (ppindex>=0) pp = net->data[ppindex].data; else pp=NULL;
    if (vpindex>=0) vp = net->data[vpindex].data; else vp=NULL;
    //printf ("\n%d %d \n",ppindex,vpindex);
    if (verbose) printf ("* Identifying critical lines and surfaces ... ");fflush(0);
    //unsigned char *ninter=NULL;
    NDNetFlags_enable(net,0);
    for (type=1;type<=net->ndims_net;type++)
    {
	
	ntagv[type]=ntagp[type]=0;

	if (net->nfaces[type]==0)
	{
	    if (verbose) printf ("(%d-faces skipped) ",type);
	    fflush(0);
	    continue;
	}
	//else if (verbose) printf ("(%d-faces)",type);
	fflush(0);

	NDNetFlags_enable(net,type);
	
	//ninter = (unsigned char *)realloc (net->nfaces[type],sizeof(unsigned char));
	for (i=0;i<net->nfaces[type];i++)
	{
	    uint *v_list;
	    int nvert;


	    NDNETFLAG_SET_VOID_INTERFACE(net->f_flag[type][i],0);
	    NDNETFLAG_SET_PEAK_INTERFACE(net->f_flag[type][i],0);


	    if (net->f_flag[type][i]&NDNETFLAG_OUT) continue;

	    
	    v_list = VERTEX_IN_FACE(net,type,i);
	    nvert = NUM_VERTEX_IN_FACE(net,type,i);
	    nppatch=nvpatch=0;
	    for (j=0;j<nvert;j++)
	    {
		if (net->v_flag[v_list[j]]&NDNETFLAG_OUT) {nppatch=nvpatch=0;break;}
		
		if (pp!=NULL)
		{
		    for (k=0;k<nppatch;k++)
			if (pp[v_list[j]] == pplist[k]) break;
		    if (k==nppatch)
		    {
			if (nppatch == nppatch_max)
			{
			    nppatch_max+=10;
			    pplist=(int*)realloc(pplist,sizeof(int)*nppatch_max);
			}
			pplist[nppatch] = pp[v_list[j]];
			nppatch++;
		    }
		}

		if (vp!=NULL)
		{
		    for (k=0;k<nvpatch;k++)
			if (vp[v_list[j]] == vplist[k]) break;
		    if (k==nvpatch)
		    {
			if (nvpatch == nvpatch_max)
			{
			    nvpatch_max+=10;
			    vplist=(int*)realloc(vplist,sizeof(int)*nvpatch_max);
			}
			vplist[nvpatch] = vp[v_list[j]];
			nvpatch++;
		    }
		}
	    }
	    
	    if (nppatch>1)
	    {
		NDNETFLAG_SET_PEAK_INTERFACE(net->f_flag[type][i],nppatch-1);
		ntagp[type]++;
	    }
	    
	    if (nvpatch>1)
	    {
		NDNETFLAG_SET_VOID_INTERFACE(net->f_flag[type][i],nvpatch-1);
		ntagv[type]++;
	    }
	}
	//free(ninter);
	if (verbose) printf ("(%dp/%dv %d-faces) ",ntagp[type],ntagv[type],type);
    }

    free (pplist);
    free (vplist);

    if (verbose) printf ("done.\n");

    return 0;
}



double NDLogOf(NDnetwork *net, int type, NDNET_UINT index, void *info)
{
    double val;

    val = net->data[*((int*) info)].data[index];
    if (val<=0) 
	return 0;
    else 
	return log(val);

}
*/

double NDVolume(NDnetwork *net, int type, NDNET_UINT index, void *info)
{
    double vol;
    float **coord=NULL;
    NDNET_UINT *vert;
    NDNET_UINT nvert;
    long i;

    if (type<=0)
    {
	fprintf(stderr,"In NDVolume: cannot compute volume of a point (type==0).\n");
	return -1;
    }
    
    nvert = NUM_VERTEX_IN_FACE(net,type,index);
    vert=VERTEX_IN_FACE(net,type,index);
    coord=realloc(coord,nvert*sizeof(float*));
    for (i=0;i<nvert;i++) 
    {
	coord[i]=&net->v_coord[vert[i]*net->ndims];
    }
    vol=SimplexVolumef(coord,nvert,net->ndims);
 
    free(coord);
    return vol;
}

double NDDensity(NDnetwork *net, int type, NDNET_UINT index, void *info)
{
    double *volume;

    if (info != NULL)
    {
	volume = (double *)info;
    }
    else volume=NULL;

    if (type==0)
	return NDVertexDensity(net,index,volume);
    else
	return NDFaceDensity(net,type,index);
}

double NDFaceDensity(NDnetwork *net, int type, NDNET_UINT index)
{
    return 1./NDVolume(net,type,index,NULL);
}

double NDVertexDensity(NDnetwork *net, NDNET_UINT index, double *volume)
{
    long i,j;
    double vol=0;
    NDNET_UINT nvert;
    NDNET_UINT *face;
    NDNET_UINT *vert;
    float **coord=NULL;

    face = FACE_IN_VERTEX(net,net->ndims,index);
    if (volume != NULL)
    {
	for (i=0;i<NUM_FACE_IN_VERTEX(net,net->ndims,index);i++)
	    vol+=volume[face[i]];
    }
    else
    {
	for (i=0;i<NUM_FACE_IN_VERTEX(net,net->ndims,index);i++)
	{
	    nvert = NUM_VERTEX_IN_FACE(net,net->ndims,face[i]);
	    vert=VERTEX_IN_FACE(net,net->ndims,face[i]);
	    coord=realloc(coord,nvert*sizeof(float*));
	    for (j=0;j<nvert;j++) 
		coord[j]=&net->v_coord[vert[j]*net->ndims];
	    vol+=SimplexVolumef(coord,nvert,net->ndims);
	}
	free(coord);
    }
    
    return 1./vol;
}

NDnetwork_Data *SmoothNDData(NDnetwork *net, int type,const char *name, uint ntimes)
{
    int index=-1;
    //int val;
    NDNET_UINT *nei=NULL;
    double *data;
    
    double *tmp_pd;
    double d;
    NDNET_UINT N;
    NDNET_UINT i,j,k;
    int nnei;

    index=NDDataIndex(net,type,name);
    if (index<0)
    {
	fprintf(stderr,"ERROR in SmoothNDData: I need a '%s' tag field.\n",name);
	exit(0);
    }

    if (verbose>1) printf ("Smoothing '%s' data field %d times ... ",name,ntimes);fflush(0);

    data = net->data[index].data;
    if (type==0) N=net->nvertex;
    else N=net->nfaces[type];

    double *buffer1 = data;
    double *buffer2 = malloc(N*sizeof(double));

    double *source=buffer1;
    double *dest=buffer2;
      
    for (k=0;k<ntimes;k++)
      {     
	for (i=0;i<N;i++)
	  {	  
	    nnei = NDFindNeighbours(net, type, i, &nei, NULL);
	    for (j=0,d=0;j<nnei;j++)
	      d+=source[nei[j]];

	    dest[i] = (d+source[i])/(1+nnei);
	  }

	tmp_pd=source;
	source=dest;
	dest=tmp_pd;
      }

    //if (ntimes&1)
      {
	tmp_pd=source;
	source=dest;
	dest=tmp_pd;
      }
    
    net->data[index].data=dest;
    free(source);
    
    free(nei);

    if (verbose>1) printf("done.\n");

    return &net->data[index];
}


int AddNDDataTagFromGrid(NDnetwork *net,NDfield *field,const char *name)
{
  long i,j,k;
  FLOAT tmp;
  INT index;
  FLOAT pos[net->ndims];
  float fpos[net->ndims];
  double *tag;
  DEFINE_NDVARS(d);
  int periodic=0;

  if (field->datatype&(ND_CHAR|ND_UCHAR|ND_SHORT|ND_USHORT|ND_INT|ND_UINT|ND_LONG))
    {
      if (verbose>1) printf ("Interpolating '%s' integer data field to vertices ...",name);fflush(0);
     
      SETNDPOINTER(d,field->val,field->datatype);
      tag=calloc(net->nvertex,sizeof(double));
      for (i=0,j=0;i<net->nvertex;i++,j++)
	{
	  for (k=0;k<net->ndims;k++) fpos[k] = net->v_coord[net->ndims*i+k];
	  Coords2IndexfND(fpos, &index,field->x0,field->delta,field->dims, field->ndims,periodic);
	  SETNDFIELD_VAL(tag[j],d,index,field->datatype);
	}
      if (verbose>1) printf (" done.\n");
      addNDDataArr(net,0,name,&tag);
    }
  else
    {
      if (verbose>1) printf ("Interpolating '%s' data field to vertices ...",name);fflush(0);
      tag=calloc(net->nvertex,sizeof(double));
      for (i=0,j=0;i<net->nvertex;i++,j++)
	{
	  for (k=0;k<net->ndims;k++) pos[k] = net->v_coord[net->ndims*i+k];
	  InterpolateND(field,&tmp,pos,-1,periodic);tag[j]=tmp;
	}
      if (verbose>1) printf (" done.\n");
      addNDDataArr(net,0,name,&tag);
    }

  //if (verbose>1) printf (" done.\n");
  return 0;
}
