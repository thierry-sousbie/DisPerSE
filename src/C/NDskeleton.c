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
#include "NDskel_tags.h"
#include "NDskeleton.h"

#include "mystring.h"
#include "global.h"

#ifdef HAVE_CFITS_IO
#include "fitsio.h"
#endif


int Free_NDskeleton(NDskel **skl_p)
{
  long i;
  NDskel *skl=*skl_p;


  for (i=0;i<skl->nnodes;i++)
    {
      free(skl->Node[i].Next);
      free(skl->Node[i].Seg);
      free(skl->Node[i].nsegs);
    }

  for (i=0;i<skl->nsegdata;i++)
    {
      free(skl->segdata_info[i]);
    }
  free(skl->segdata_info);

  for (i=0;i<skl->nnodedata;i++)
    {
      free(skl->nodedata_info[i]);    
    }
  free(skl->nodedata_info);

  free(skl->segpos);free(skl->nodepos);
  if (skl->nsegdata) free(skl->segdata);
  if (skl->nnodedata) free(skl->nodedata);
  free(skl->Seg);
  free(skl->Node);
  free(skl->dims);
  free(skl->x0);
  free(skl->delta);

  free(*skl_p);
  *skl_p=NULL;
  
  return 0;
}

int NDskel_SegDataIndex(NDskel *skl, const char *name)
{
    int i;
    
    for (i=0;i<skl->nsegdata;i++)
    {
      if (!strcmp(name,skl->segdata_info[i]))
	return i;
    }

    return -1;
}

int NDskel_NodeDataIndex(NDskel *skl, const char *name)
{
    int i;
    
    for (i=0;i<skl->nnodedata;i++)
    {
      if (!strcmp(name,skl->nodedata_info[i]))
	return i;
    }
    
    return -1;
}

double ComputeSegLen(NDskel *skl, NDskl_seg* seg)
{
  int i;
  double d=0;
  
  for (i=0;i<skl->ndims;i++)
    d+= pow(seg->pos[i+skl->ndims]-seg->pos[i],2);

  return sqrt(d);
}

double ComputeDistToNext(NDskel *skl, NDskl_seg* seg_p)
{
  int i;
  double d=0;
  double dt=0;
  NDskl_seg* seg = seg_p;

  do 
    {
      //printf ("index:%d, seg=%d\n",seg->index,seg);
      for (i=0,d=0;i<skl->ndims;i++)
	d+= pow(seg->pos[i+skl->ndims]-seg->pos[i],2);
      dt+=sqrt(d);
    } while ((seg=seg->Next)!=NULL);

  return dt;
}

double ComputeDistFromPrev(NDskel *skl, NDskl_seg* seg_p)
{
  int i;
  double d=0;
  double dt=0;
  NDskl_seg* seg = seg_p;

  while ((seg=seg->Prev)!=NULL)
  {
      for (i=0,d=0;i<skl->ndims;i++)
	d+= pow(seg->pos[i+skl->ndims]-seg->pos[i],2);
      dt+=sqrt(d);
  }

  return dt;
}
/*
double ComputeDistToSaddle(NDskel *skl, NDskl_seg* seg_p)
{
  int i;
  double d=0;
  double dt=0;
  NDskl_seg* seg = seg_p;

  while ((seg!=NULL)&&(!(seg->flags&FLAG_NDSKL_SAD)))
    {
	for (i=0,d=0;i<skl->ndims;i++)
	    d+= pow(seg->pos[i+skl->ndims]-seg->pos[i],2);
	dt+=sqrt(d);	
	seg = seg->Next;
    }
  if (seg!=NULL) 
  {
      for (i=0,d=0;i<skl->ndims;i++)
	    d+= pow(seg->pos[i+skl->ndims]-seg->pos[i],2);
      dt +=sqrt(d)/2;
  }
  else
  {
      seg = seg_p->Prev;
      dt=0;
      while ((seg!=NULL)&&(!(seg->flags&FLAG_NDSKL_SAD)))
      {
	  for (i=0,d=0;i<skl->ndims;i++)
	      d+= pow(seg->pos[i+skl->ndims]-seg->pos[i],2);
	  dt-=sqrt(d);	
	  seg = seg->Prev;
      }
      if (seg!=NULL)
      {
	  for (i=0,d=0;i<skl->ndims;i++)
	      d+= pow(seg->pos[i+skl->ndims]-seg->pos[i],2);
	  dt -=sqrt(d)/2;
      }
  }

  //if (seg==NULL) printf ("ERROOORRRR\n");

  //{
  //   int n=0;
  //   for (i=0;i<skl->nsegs;i++)
  //	  if (skl->Seg[i].flags&FLAG_NDSKL_SAD) n++;
  //
  //      printf ("n=%d\n",n);exit(0);
  //}

  return dt;
}
*/
int NDskelCheckSanity(NDskel *skl, int periodicity)
{
    int i,j,k;
    int nbugs;
    int ndisc;
    int npass=1;

    if (skl==NULL) return -1;

    printf("* Checking skeleton sanity ... ");

    do {
	nbugs=0;
	ndisc=0;
	for (i=0;i<skl->nnodes;i++)
	{
	    //printf ("node %d\n",i);
	    NDskl_node *node;
	    node = &skl->Node[i];
	    
	    for (j=0;j<node->nnext;j++)
	    {
		NDskl_seg *seg;
		NDskl_node *nextnode;
		NDskl_seg *curseg;
		NDskl_seg *oldseg;
		double cmpval[skl->ndims];
		
		for (k=0;k<skl->ndims;k++)
		    cmpval[k] = 0.01*skl->delta[k]/skl->dims[k];

		seg = node->Seg[j];
		oldseg=seg;
		nextnode = node->Next[j];

		for (k=0;k<skl->ndims;k++)
		    if ((fabs(node->pos[k]-seg->pos[k])>cmpval[k])&&
			(fabs(node->pos[k]-seg->pos[k+skl->ndims])>cmpval[k])&&
			((!(periodicity & (1<<k)))||
			 ((fabs(fabs(node->pos[k]-seg->pos[k])-skl->delta[k])>cmpval[k])&&
			  (fabs(fabs(node->pos[k]-seg->pos[k])-skl->delta[k])<0.5*skl->delta[k])&&
			  (fabs(fabs(node->pos[k]-seg->pos[k+skl->ndims])-skl->delta[k])>cmpval[k])&&
			  (fabs(fabs(node->pos[k]-seg->pos[k+skl->ndims])-skl->delta[k])<0.5*skl->delta[k]))))
			break;

		if (k!=skl->ndims)
		{
		    if (!ndisc)
		    {
			printf ("Node discontinuity :\n");
			printf ("nodepos=[%f %f %f]\n",
				node->pos[0],node->pos[1],node->pos[2]);
			printf ("nextpos=[%f %f %f],[%f %f %f]\n",
			    seg->pos[0],seg->pos[1],seg->pos[2],
			    seg->pos[3],seg->pos[4],seg->pos[5]);
		    }
/*
		    for (k=0;k<skl->nsegs;k++)
		    {
			printf ("LOOKING FOR LOST NODE ;)");
		    }
		    */
		    nbugs++;
		    ndisc++;
		}

		//else printf ("OK\n");
		
		if ((seg->Next!=NULL)&&(seg->Prev!=NULL))
		{
		    //printf("PROBLEM 1\n");
		    nbugs++;
		}
		
		if ((seg->Next!=NULL)||(seg->Prev!=NULL))
		{
		    if (seg->Next!=NULL)
		    {
			
			for (oldseg=seg,curseg=seg->Next;curseg!=NULL;oldseg=curseg,curseg=curseg->Next)
			{

			    for (k=0;k<skl->ndims;k++)
				if (curseg->pos[k]!=oldseg->pos[k]) break;
			    if (k!=skl->ndims)
				for (k=0;k<skl->ndims;k++)
				    if (curseg->pos[k]!=oldseg->pos[k+skl->ndims]) break;
			    if (k!=skl->ndims)
				for (k=0;k<skl->ndims;k++)
				    if (curseg->pos[k+skl->ndims]!=oldseg->pos[k+skl->ndims]) break;
			    if (k!=skl->ndims)
				for (k=0;k<skl->ndims;k++)
				    if (curseg->pos[k+skl->ndims]!=oldseg->pos[k]) break;
			    if (k!=skl->ndims) 
			    {
				for (k=0;k<skl->ndims;k++)
				    if (fabs(curseg->pos[k]-oldseg->pos[k])>skl->delta[k]*0.5)
					break;

				if (k==skl->ndims)
				{
				    printf ("segment discontinuity ERROR (Next dir)\n");
				    printf ("curpos=[%f %f %f],[%f %f %f]\n",
					    curseg->pos[0],curseg->pos[1],curseg->pos[2],
					    curseg->pos[3],curseg->pos[4],curseg->pos[5]);
				    printf ("prepos=[%f %f %f],[%f %f %f]\n",
					    oldseg->pos[0],oldseg->pos[1],oldseg->pos[2],
					    oldseg->pos[3],oldseg->pos[4],oldseg->pos[5]);
				}
			    }
			    //else printf ("OK\n");
			    //printf ("seg->index = %d\n",seg->index);
			    if (curseg->Prev != oldseg)
			    {
				
				if (curseg->Prev == NULL)
				{
/*
				    //printf("Fixed NULL value Prev.\n");
				    printf ("\ncurpos=[%f %f %f],[%f %f %f]\n",
					    curseg->pos[0],curseg->pos[1],curseg->pos[2],
					    curseg->pos[3],curseg->pos[4],curseg->pos[5]);
				    printf ("\nprepos=[%f %f %f],[%f %f %f]\n",
					    oldseg->pos[0],oldseg->pos[1],oldseg->pos[2],
					    oldseg->pos[3],oldseg->pos[4],oldseg->pos[5]);
*/
				    curseg->Prev = oldseg;
				    nbugs++;
				}
				else
				{
				    //printf("Fixed non-NULL value Prev.\n");
				    curseg->Prev = oldseg;
				    nbugs++;
				}
				
			    }
			}
			
		    }
		    if (seg->Prev!=NULL)
		    {
			for (oldseg=seg,curseg=seg->Prev;curseg!=NULL;oldseg=curseg,curseg=curseg->Prev)
			{
			    for (k=0;k<skl->ndims;k++)
				if (curseg->pos[k]!=oldseg->pos[k]) break;
			    if (k!=skl->ndims)
				for (k=0;k<skl->ndims;k++)
				    if (curseg->pos[k]!=oldseg->pos[k+skl->ndims]) break;
			    if (k!=skl->ndims)
				for (k=0;k<skl->ndims;k++)
				    if (curseg->pos[k+skl->ndims]!=oldseg->pos[k+skl->ndims]) break;
			    if (k!=skl->ndims)
				for (k=0;k<skl->ndims;k++)
				    if (curseg->pos[k+skl->ndims]!=oldseg->pos[k]) break;
			    
			    if (k!=skl->ndims) 
			    {
			
				for (k=0;k<skl->ndims;k++)
				    if (fabs(curseg->pos[k]-oldseg->pos[k])>skl->delta[k]*0.5)
					break;

				if (k==skl->ndims)
				{
				    printf ("segment discontinuity ERROR (Prev dir)\n");
				    printf ("curpos=[%f %f %f],[%f %f %f]\n",
					    curseg->pos[0],curseg->pos[1],curseg->pos[2],
					    curseg->pos[3],curseg->pos[4],curseg->pos[5]);
				    printf ("prepos=[%f %f %f],[%f %f %f]\n",
					    oldseg->pos[0],oldseg->pos[1],oldseg->pos[2],
					    oldseg->pos[3],oldseg->pos[4],oldseg->pos[5]);
				    
				}
			    }

			//else printf ("OK\n");
			    //printf ("seg->index = %d\n",seg->index);
			    if (curseg->Next != oldseg)
			    {
				
				if (curseg->Next == NULL)
				{
				    //printf("Fixed NULL value Next.\n");
				    curseg->Next = oldseg;
				    nbugs++;
				}
				else
				{
				    //printf("Fixed non-NULL value Next.\n");
				    curseg->Next = oldseg;
				    nbugs++;
				}
				
			    }
			}
		    }
		}
		else {}
		
		seg=oldseg;
		for (k=0;k<nextnode->nnext;k++)
		{
		    if (nextnode->Seg[k] == seg)
			break;
		}
		if (k==nextnode->nnext)
		{
		    
		    nbugs++;
		    for (k=0;k<nextnode->nnext;k++)
		    {
			if (nextnode->Next[k] == node)
			break;
		    }

		    printf ("node %d does not lead to node %d\n",node->index,nextnode->index);
		    seg=node->Seg[j];
		    printf ("seg = %d ",seg->index);
		    for (oldseg=seg,curseg=seg->Prev;curseg!=NULL;oldseg=curseg,curseg=curseg->Prev)
			printf ("(%d->%d)",oldseg->index,curseg->index);
		    for (oldseg=seg,curseg=seg->Next;curseg!=NULL;oldseg=curseg,curseg=curseg->Next)
			printf ("(%d->%d)",oldseg->index,curseg->index);
		    printf ("= %d.\n",nextnode->Seg[k]->index);
		    
		    seg=nextnode->Seg[k];
		    printf ("seg = %d ",seg->index);
		    for (oldseg=seg,curseg=seg->Prev;curseg!=NULL;oldseg=curseg,curseg=curseg->Prev)
			printf ("(%d->%d)",oldseg->index,curseg->index);
		    for (oldseg=seg,curseg=seg->Next;curseg!=NULL;oldseg=curseg,curseg=curseg->Next)
			printf ("(%d->%d)",oldseg->index,curseg->index);
		    printf ("= %d.\n",node->Seg[j]->index);

		    exit(0);
		    
		}
		
	    }
	}
	printf ("(pass %d: %d bugs (%d disc.)) ",npass++,nbugs,ndisc);
	//exit(0);
    } while ((nbugs)&&(npass<3));

    if ((npass>1) && (nbugs))
    {
	printf (" could not fix %d bugs (%d disc).\n",nbugs,ndisc);
	if (ndisc==nbugs)
	{
	    printf (" This may be because the skeleton has been computed from an old version of rsex.\n");
	    printf (" Don't worry too much, this should be OK :)\n");
	}
    }
    else
	printf ("done.\n");


    return nbugs;
}
// the segments will go frm a node of higher nodefieldname value down to a lower one
// fieldname can be "" in which case no sort is done
// After calling this function, the two points within each segment are
// sorted so that the first point of next segment is the same as the second point
// of current segment 
int SortNDskelSegments(NDskel *skl,char *nodefieldname,int periodic)
{
  int index = -1;
  long i,j,k;
  int reverse=0;
  NDskl_seg *tmpseg;
  NDskl_seg *seg;
  double d1,d2,d;
  float pos[skl->ndims];

  int dataexch[skl->nsegdata];

  for (i=0;i<skl->nsegdata;i++) dataexch[i]=-1;
  for (i=0;i<skl->nsegdata;i++)
    {
      if (!strcmp(&(skl->segdata_info[i][strlen(skl->segdata_info[i])-3]),"_p1"))
	{
	  for (j=0;j<skl->nsegdata;j++)
	    if (!strncmp(skl->segdata_info[i],skl->segdata_info[j],strlen(skl->segdata_info[i])-1))
	      {
		if (i!=j) 
		  dataexch[i]=j;
	      }
	}
    }
  
  for (i=0;i<skl->nnodedata;i++)
    {
      if (!strcmp(skl->nodedata_info[i],nodefieldname)) index = i;
    }
  //printf ("node %s : %d\n",nodefieldname,index);

  //printf ("nnodes : %d\n",skl->nnodes);
  for (i=0;i<skl->nnodes;i++)
    {
      //printf ("nnext : %d\n",skl->Node[i].nnext);
      for (j=0;j<skl->Node[i].nnext;j++)
	{

	  if (skl->Node[i].Next[j]->index!=skl->Node[i].Seg[j]->nodes[1]) continue;

	  if (index!=-1)
	    {
	      if (skl->Node[i].Next[j]->data[index]<skl->Node[i].data[index])
		reverse=0;
	      else
		reverse=1;
	    }

	  seg = skl->Node[i].Seg[j];

	  do
	    {
	      if (seg->Prev==NULL)
		{
		  for (k=0,d1=d2=0;k<skl->ndims;k++)
		    {
		      d=fabs(seg->pos[k]-skl->Node[i].pos[k]);
		      if (d>skl->delta[k]/2) d-=skl->delta[k];
		      d1+= d*d;
		      d=fabs(seg->pos[k+skl->ndims]-skl->Node[i].pos[k]);
		      if (d>skl->delta[k]/2) d-=skl->delta[k];
		      d2+= d*d;
		    }
		  //printf ("Next: %e %e \n",d1,d2);
		}
	      else if (seg->Next==NULL)
		{
		  for (k=0,d1=d2=0;k<skl->ndims;k++)
		    {
		      d=fabs(seg->pos[k]-skl->Node[i].Next[j]->pos[k]);
		      if (d>skl->delta[k]/2) d-=skl->delta[k];
		      d2+= d*d;
		      d=fabs(seg->pos[k+skl->ndims]-skl->Node[i].Next[j]->pos[k]);
		      if (d>skl->delta[k]/2) d-=skl->delta[k];
		      d1+= d*d;
		    }
		  //printf ("Next: %e %e \n",d1,d2);
		}
	      else 
		for (k=0,d1=d2=0;k<skl->ndims;k++)
		  {
		    d=fabs(seg->Prev->pos[k+skl->ndims]-seg->pos[k]);
		    if (d>skl->delta[k]/2) d-=skl->delta[k];
		    d1+= d*d;
		    
		    d=fabs(seg->Prev->pos[k+skl->ndims]-seg->pos[k+skl->ndims]);
		    if (d>skl->delta[k]/2) d-=skl->delta[k];
		    d2+= d*d;
		  }
	      
	      if ((d2<d1))//||(reverse&(d1<=d2)))
	      {
		  //wrong direction exchange pos
		  memcpy(pos,seg->pos,skl->ndims*sizeof(float));
		  memcpy(seg->pos,&(seg->pos[skl->ndims]),skl->ndims*sizeof(float));
		  memcpy(&(seg->pos[skl->ndims]),pos,skl->ndims*sizeof(float));		  
		  //data p1/p2
		  for (k=0;k<skl->nsegdata;k++)
		    {
		      if (dataexch[k]!=-1)
			{
			  //data = seg->data[k];
			  seg->data[k]=seg->data[dataexch[k]];
			  seg->data[dataexch[k]] = seg->data[k];
			}
		    }		  
		}

	      seg=seg->Next;
	
	    } while (seg!=NULL);

	  if (reverse)
	    {
	      seg = skl->Node[i].Seg[j];
	      do {
		memcpy(pos,seg->pos,skl->ndims*sizeof(float));
		memcpy(seg->pos,&(seg->pos[skl->ndims]),skl->ndims*sizeof(float));
		memcpy(&(seg->pos[skl->ndims]),pos,skl->ndims*sizeof(float));		  
		//data p1/p2
		for (k=0;k<skl->nsegdata;k++)
		  {
		    if (dataexch[k]!=-1)
		      {
			//data = seg->data[k];
			seg->data[k]=seg->data[dataexch[k]];
			seg->data[dataexch[k]] = seg->data[k];
		      }
		  }	
		
		k=seg->nodes[0];
		seg->nodes[0]=seg->nodes[1];
		seg->nodes[1]=k;
		tmpseg=seg->Next;
		seg->Next=seg->Prev;
		seg->Prev=tmpseg;
		seg=seg->Prev;
	      } while (seg!=NULL);
	    }
	}
    }
  
  return 0;
}

NDskel *Create_NDskel(int *dims,int ndims,double *x0, double *delta,const char *comment, int nnodes, int nsegs)
{
  NDskel *skl;
  
  skl=calloc(1,sizeof(NDskel));
  skl->ndims=ndims;
  skl->dims = malloc(sizeof(int)*ndims);
  skl->x0= malloc(sizeof(double)*ndims);
  skl->delta= malloc(sizeof(double)*ndims);
  if (dims!=NULL)
    memcpy(skl->dims,dims,sizeof(int)*ndims);
  if (x0!=NULL)
    memcpy(skl->x0,x0,sizeof(double)*ndims);
  if (delta!=NULL)
    memcpy(skl->delta,delta,sizeof(double)*ndims);
  if (comment!=NULL)
    strcpy(skl->comment,comment);
  else
    strcpy(skl->comment,"No comment");

  if (nnodes>0)
    {
      long i;
      skl->nnodes=nnodes;
      skl->Node = (NDskl_node *) calloc(skl->nnodes,sizeof(NDskl_node));
      skl->nodepos = (float *) calloc(skl->nnodes,sizeof(float)*skl->ndims);
      for (i=0;i<skl->nnodes;i++)
	{
	  skl->Node[i].pos=&skl->nodepos[i*skl->ndims];      
	  skl->Node[i].data=NULL;
	}
    }

  if (nnodes>0)
    {
      long i;
      skl->nsegs=nsegs;
      skl->Seg = (NDskl_seg *) calloc(skl->nsegs,sizeof(NDskl_seg));
      skl->segpos = (float *) calloc(skl->nsegs,sizeof(float)*skl->ndims*2);
      for (i=0;i<skl->nsegs;i++) 
	{
	  skl->Seg[i].pos=&skl->segpos[i*skl->ndims*2];
	  skl->Seg[i].data=NULL;
	}
    }

  return skl;
}

NDskel *Copy_NDskel(NDskel *oskl)
{
  NDskel *skl;
  long i,j;

  skl=calloc(1,sizeof(NDskel));
  memcpy(skl,oskl,sizeof(NDskel));
  skl->dims=calloc(oskl->ndims,sizeof(int));
  skl->x0=calloc(oskl->ndims,sizeof(double));
  skl->delta=calloc(oskl->ndims,sizeof(double));
  
  memcpy(skl->dims,oskl->dims,oskl->ndims*sizeof(int));
  memcpy(skl->x0,oskl->x0,oskl->ndims*sizeof(double));
  memcpy(skl->delta,oskl->delta,oskl->ndims*sizeof(double));
  
  skl->segdata_info = calloc(skl->nsegdata,sizeof(char *));
  skl->nodedata_info = calloc(skl->nnodedata,sizeof(char *));
 
  for (i=0;i<skl->nsegdata;i++)
    {
      skl->segdata_info[i]=calloc(NDSKEL_DATA_STR_SIZE,sizeof(char));
      strcpy(skl->segdata_info[i],oskl->segdata_info[i]);
      //memcpy(&skl->segdata_info[i],&oskl->segdata_info[i],sizeof(char)*NDSKEL_DATA_STR_SIZE);
    }
  for (i=0;i<skl->nnodedata;i++)
    {
      skl->nodedata_info[i]=calloc(NDSKEL_DATA_STR_SIZE,sizeof(char));
      strcpy(skl->nodedata_info[i],oskl->nodedata_info[i]);
      //memcpy(&skl->nodedata_info[i],&oskl->nodedata_info[i],sizeof(char)*NDSKEL_DATA_STR_SIZE);
    }
  
  skl->segpos = malloc((long)sizeof(float)*skl->nsegs*2*skl->ndims);
  memcpy(skl->segpos,oskl->segpos,sizeof(float)*skl->nsegs*2*skl->ndims);
 
  skl->nodepos = malloc((long)sizeof(float)*skl->nnodes*skl->ndims);
  memcpy(skl->nodepos,oskl->nodepos,sizeof(float)*skl->nnodes*skl->ndims);
  
  skl->segdata = malloc((long)sizeof(double)*skl->nsegs*skl->nsegdata);
  memcpy(skl->segdata,oskl->segdata,sizeof(double)*skl->nsegs*skl->nsegdata);
  
  skl->nodedata = malloc((long)sizeof(double)*skl->nnodes*skl->nnodedata);
  memcpy(skl->nodedata,oskl->nodedata,sizeof(double)*skl->nnodes*skl->nnodedata);
  
  skl->Node = malloc((long)sizeof(NDskl_node)*skl->nnodes);
  skl->Seg = malloc((long)sizeof(NDskl_seg)*skl->nsegs);
  
  memcpy(skl->Node,oskl->Node,skl->nnodes*sizeof(NDskl_node));
  memcpy(skl->Seg,oskl->Seg,skl->nsegs*sizeof(NDskl_seg));

  for (i=0;i<skl->nnodes;i++)
    {
      memcpy(&skl->Node[i],&oskl->Node[i],sizeof(NDskl_node));
      skl->Node[i].Next = malloc(skl->Node[i].nnext*sizeof(NDskl_node *));
      skl->Node[i].Seg = malloc(skl->Node[i].nnext*sizeof(NDskl_seg *));
      skl->Node[i].nsegs = malloc(skl->Node[i].nnext*sizeof(int));
      skl->Node[i].data = oskl->Node[i].data+(skl->nodedata-oskl->nodedata);
      skl->Node[i].pos = oskl->Node[i].pos + (skl->nodepos-oskl->nodepos);

      for (j=0;j<skl->Node[i].nnext;j++)
	{
	  skl->Node[i].Next[j] = oskl->Node[i].Next[j]+(skl->Node-oskl->Node);
	  skl->Node[i].Seg[j] = oskl->Node[i].Seg[j]+(skl->Seg-oskl->Seg); 
	  skl->Node[i].nsegs[j] = oskl->Node[i].nsegs[j];
	}
    }
  
  for (i=0;i<skl->nsegs;i++)
    {
      memcpy(&skl->Seg[i],&oskl->Seg[i],sizeof(NDskl_seg));
      //memcpy(skl->Seg[i].nodes,oskl->Seg[i].nodes,sizeof(int)*2);
      skl->Seg[i].pos = oskl->Seg[i].pos + (skl->segpos-oskl->segpos);
      skl->Seg[i].data = oskl->Seg[i].data + (skl->segdata-oskl->segdata);
      if (oskl->Seg[i].Next!=NULL) skl->Seg[i].Next = oskl->Seg[i].Next + (skl->Seg-oskl->Seg);
      if (oskl->Seg[i].Prev!=NULL) skl->Seg[i].Prev = oskl->Seg[i].Prev + (skl->Seg-oskl->Seg);
    }

  return skl;
}

int NDskel_SetNNextInNode(NDskl_node *node,int nnodes)
{
  if (nnodes==0) return 1;
  
  node->Next = realloc (node->Next,nnodes * sizeof(NDskl_node*));
  node->Seg = realloc (node->Seg,nnodes * sizeof(NDskl_seg*));
  node->nsegs = realloc(node->nsegs,nnodes*sizeof(int));
  
  return 0;
}

int NDskel_realloc(NDskel *skl,int nsegs,int nnodes,int nsegalloc,int nnodealloc, int delta)
{
  long i,j;
  int result = 0;
  NDskl_seg *oldseg=skl->Seg;
  NDskl_node *oldnode=skl->Node;
  float *oldsegpos=skl->segpos;
  float *oldnodepos=skl->nodepos;
  double *oldnodedata=skl->nodedata;
  double *oldsegdata=skl->segdata;

  if ((nsegalloc==0)||(!(nsegs%nsegalloc))) 
    {
      skl->Seg = realloc(skl->Seg,(nsegs+nsegalloc+delta)*sizeof(NDskl_seg));
      memset(&skl->Seg[skl->nsegs],0,(nsegs+nsegalloc+delta-skl->nsegs)*sizeof(NDskl_seg));
      skl->segpos = realloc(skl->segpos,(nsegs+nsegalloc+delta)*skl->ndims*2*sizeof(float));
      if (skl->nsegdata!=0)
	skl->segdata = realloc(skl->segdata,(nsegs+nsegalloc+delta)*skl->nsegdata*sizeof(double));

      if (oldsegpos != skl->segpos)
	{
	  for (i=0;i<(nsegs+nsegalloc+delta);i++)
	    skl->Seg[i].pos = &skl->segpos[skl->ndims*2*i];
	}
      else
	{
	  for (i=skl->nsegs;i<(nsegs+nsegalloc+delta);i++)
	    skl->Seg[i].pos = &skl->segpos[skl->ndims*2*i];
	}
      if (oldsegdata != skl->segdata)
	{
	  for (i=0;i<(nsegs+nsegalloc+delta);i++)
	    skl->Seg[i].data = &skl->segdata[i*skl->nsegdata];
	}
      else
	{
	  for (i=skl->nsegs;i<(nsegs+nsegalloc+delta);i++)
	    skl->Seg[i].data = &skl->segdata[i*skl->nsegdata];
	}
      result|=(1<<0);
    }

  if ((nnodealloc==0)||(!(nnodes%nnodealloc))) 
    {
      //printf ("hello (%d %d %d %d)\n",nsegs,nnodes,nnodealloc,nnodes%nnodealloc);
      skl->Node = realloc(skl->Node,(nnodes+nnodealloc+delta)*sizeof(NDskl_node));
      memset(&skl->Node[skl->nnodes],0,(nnodes+nnodealloc+delta-skl->nnodes)*sizeof(NDskl_node));
      skl->nodepos = realloc(skl->nodepos,(nnodes+nnodealloc+delta)*skl->ndims*sizeof(float));
      if (skl->nnodedata!=0)
	skl->nodedata = realloc(skl->nodedata,(nnodes+nnodealloc+delta)*skl->nnodedata*sizeof(double));

      if (oldnodepos != skl->nodepos)
	{
	  for (i=0;i<(nnodes+nnodealloc+delta);i++)
	    skl->Node[i].pos = &skl->nodepos[skl->ndims*i];
	}
      else
	{
	  for (i=skl->nnodes;i<(nnodes+nnodealloc+delta);i++)
	    skl->Node[i].pos = &skl->nodepos[skl->ndims*i];
	}
      if (oldnodedata != skl->nodedata)
	{
	  for (i=0;i<(nnodes+nnodealloc+delta);i++)
	    skl->Node[i].data = &skl->nodedata[i*skl->nnodedata];
	}
      else
	{
	  for (i=skl->nnodes;i<(nnodes+nnodealloc+delta);i++)
	    skl->Node[i].data = &skl->nodedata[i*skl->nnodedata];
	}
      result|=(1<<1);
    }

  if (oldseg!=skl->Seg)
    {
      //printf ("Hello ALLLLLLLLLL\n\n\n");
      for (i=0;i<skl->nsegs;i++)
	{
	  if (skl->Seg[i].Next!=NULL) skl->Seg[i].Next += skl->Seg-oldseg;
	  if (skl->Seg[i].Prev!=NULL) skl->Seg[i].Prev += skl->Seg-oldseg;
	}
      for (i=0;i<skl->nnodes;i++)
	{
	  for (j=0;j<skl->Node[i].nnext;j++)
	    if (skl->Node[i].Seg[j]!=NULL) skl->Node[i].Seg[j] += skl->Seg-oldseg;
	}
    }
  
  if (oldnode!=skl->Node)
    {
      for (i=0;i<skl->nnodes;i++)
	{
	  for (j=0;j<skl->Node[i].nnext;j++)
	    if (skl->Node[i].Next[j]!=NULL) skl->Node[i].Next[j] += skl->Node-oldnode;
	}
    }
  

  return result;
}

int IsNDskeleton(const char *filename)
{
    int i;
    char tag[NDSKEL_DATA_STR_SIZE];
    FILE *f;
    int swap=0;
        
    memset(tag,0,16*sizeof(char));
    
    f=fopen(filename,"r");
    fread_sw(&i,sizeof(int),1,f,swap);
    if (i!=16) swap = 1-swap;
    fread_sw(tag,sizeof(char),16,f,swap);
    fread_sw(&i,sizeof(int),1,f,swap);
    fclose(f);
    tag[15]='\0';
    if (strcmp(tag,NDSKEL_TAG)) return 0;
    
    return 1;

}

NDskel *Load_NDskelHeader(const char *filename)
{
    int i,j;
  char tag[NDSKEL_DATA_STR_SIZE];
  char dummy[160];
  FILE *f;
  int swap=0;
  NDskel *skl;
  //int index;

  //NDskl_node *node;
  //NDskl_seg *seg;

  skl = calloc(1,sizeof(NDskel));
  memset(tag,0,16*sizeof(char));
  
  f=fopen(filename,"r");
  fread_sw(&i,sizeof(int),1,f,swap);
  if (i!=16) swap = 1-swap;
  fread_sw(tag,sizeof(char),16,f,swap);
  fread_sw(&i,sizeof(int),1,f,swap);

  if (strcmp(tag,NDSKEL_TAG))
    {
      fclose(f);
      fprintf (stderr,"File %s has an unknown format.\n",filename);
      return NULL;
    }
  
  fread_sw(&j,sizeof(int),1,f,swap);
  fread_sw(skl->comment,sizeof(char),80,f,swap);
  fread_sw(&skl->ndims,sizeof(int),1,f,swap);

  skl->dims = malloc(skl->ndims*sizeof(int));
  skl->x0 = malloc(skl->ndims*sizeof(double));
  skl->delta = malloc(skl->ndims*sizeof(double));

  fread_sw(skl->dims,sizeof(int),skl->ndims,f,swap);
  if (skl->ndims<NDSKEL_MAX_DIMS) fread_sw(dummy,sizeof(int),NDSKEL_MAX_DIMS-skl->ndims,f,swap);
  fread_sw(skl->x0,sizeof(double),skl->ndims,f,swap);
  if (skl->ndims<NDSKEL_MAX_DIMS) fread_sw(dummy,sizeof(double),NDSKEL_MAX_DIMS-skl->ndims,f,swap);
  fread_sw(skl->delta,sizeof(double),skl->ndims,f,swap);
  if (skl->ndims<NDSKEL_MAX_DIMS) fread_sw(dummy,sizeof(double),NDSKEL_MAX_DIMS-skl->ndims,f,swap);
  
  fread_sw(&skl->nsegs,sizeof(int),1,f,swap);
  fread_sw(&skl->nnodes,sizeof(int),1,f,swap);
  fread_sw(&skl->nsegdata,sizeof(int),1,f,swap);
  fread_sw(&skl->nnodedata,sizeof(int),1,f,swap);
  fread_sw(&j,sizeof(int),1,f,swap);

  skl->segdata_info = calloc(skl->nsegdata,sizeof(char *));
  skl->nodedata_info = calloc(skl->nnodedata,sizeof(char *));
 
  if (skl->nsegdata) fread_sw(&j,sizeof(int),1,f,swap);
  for (i=0;i<skl->nsegdata;i++)
    {
      skl->segdata_info[i]=calloc(NDSKEL_DATA_STR_SIZE,sizeof(char));
      fread_sw(skl->segdata_info[i],sizeof(char),NDSKEL_DATA_STR_SIZE,f,swap);
    }
  if (skl->nsegdata) fread_sw(&j,sizeof(int),1,f,swap);

  if (skl->nnodedata) fread_sw(&j,sizeof(int),1,f,swap);
  for (i=0;i<skl->nnodedata;i++)
    {
      skl->nodedata_info[i]=calloc(NDSKEL_DATA_STR_SIZE,sizeof(char));
      fread_sw(skl->nodedata_info[i],sizeof(char),NDSKEL_DATA_STR_SIZE,f,swap);
    }
  if (skl->nnodedata) fread_sw(&j,sizeof(int),1,f,swap);

  fclose(f);
  
  return skl;
}

NDskel *Load_NDskel(const char *filename)
{
  int i,j;
  char tag[NDSKEL_DATA_STR_SIZE];
  char dummy[160];
  FILE *f;
  int swap=0;
  NDskel *skl;
  int index;

  NDskl_node *node;
  NDskl_seg *seg;

  skl = calloc(1,sizeof(NDskel));
  memset(tag,0,16*sizeof(char));
  
  f=fopen(filename,"r");
  fread_sw(&i,sizeof(int),1,f,swap);
  if (i!=16) swap = 1-swap;
  fread_sw(tag,sizeof(char),16,f,swap);
  fread_sw(&i,sizeof(int),1,f,swap);

  if (strcmp(tag,NDSKEL_TAG))
    {
      fclose(f);
      //fprintf (stderr,"File %s has an unknown format.\n",filename);
      return NULL;
    }
  
  fread_sw(&j,sizeof(int),1,f,swap);
  fread_sw(skl->comment,sizeof(char),80,f,swap);
  fread_sw(&skl->ndims,sizeof(int),1,f,swap);

  printf ("Loading %dD skeleton from file %s ...",skl->ndims,filename);fflush(0);

  skl->dims = malloc(skl->ndims*sizeof(int));
  skl->x0 = malloc(skl->ndims*sizeof(double));
  skl->delta = malloc(skl->ndims*sizeof(double));

  fread_sw(skl->dims,sizeof(int),skl->ndims,f,swap);
  if (skl->ndims<NDSKEL_MAX_DIMS) fread_sw(dummy,sizeof(int),NDSKEL_MAX_DIMS-skl->ndims,f,swap);
  fread_sw(skl->x0,sizeof(double),skl->ndims,f,swap);
  if (skl->ndims<NDSKEL_MAX_DIMS) fread_sw(dummy,sizeof(double),NDSKEL_MAX_DIMS-skl->ndims,f,swap);
  fread_sw(skl->delta,sizeof(double),skl->ndims,f,swap);
  if (skl->ndims<NDSKEL_MAX_DIMS) fread_sw(dummy,sizeof(double),NDSKEL_MAX_DIMS-skl->ndims,f,swap);
  
  fread_sw(&skl->nsegs,sizeof(int),1,f,swap);
  fread_sw(&skl->nnodes,sizeof(int),1,f,swap);
  fread_sw(&skl->nsegdata,sizeof(int),1,f,swap);
  fread_sw(&skl->nnodedata,sizeof(int),1,f,swap);
  fread_sw(&j,sizeof(int),1,f,swap);

  skl->segdata_info = calloc(skl->nsegdata,sizeof(char *));
  skl->nodedata_info = calloc(skl->nnodedata,sizeof(char *));
 
  if (skl->nsegdata) fread_sw(&j,sizeof(int),1,f,swap);
  for (i=0;i<skl->nsegdata;i++)
    {
      skl->segdata_info[i]=calloc(NDSKEL_DATA_STR_SIZE,sizeof(char));
      fread_sw(skl->segdata_info[i],sizeof(char),NDSKEL_DATA_STR_SIZE,f,swap);
    }
  if (skl->nsegdata) fread_sw(&j,sizeof(int),1,f,swap);

  if (skl->nnodedata) fread_sw(&j,sizeof(int),1,f,swap);
  for (i=0;i<skl->nnodedata;i++)
    {
      skl->nodedata_info[i]=calloc(NDSKEL_DATA_STR_SIZE,sizeof(char));
      fread_sw(skl->nodedata_info[i],sizeof(char),NDSKEL_DATA_STR_SIZE,f,swap);
    }
  if (skl->nnodedata) fread_sw(&j,sizeof(int),1,f,swap);
  
  skl->segpos = malloc((long)sizeof(float)*skl->nsegs*2*skl->ndims);
  fread_sw(&j,sizeof(int),1,f,swap);
  fread_sw(skl->segpos,sizeof(float),skl->nsegs*2*skl->ndims,f,swap);
  fread_sw(&j,sizeof(int),1,f,swap);
  
  skl->nodepos = malloc((long)sizeof(float)*skl->nnodes*skl->ndims);
  fread_sw(&j,sizeof(int),1,f,swap);
  fread_sw(skl->nodepos,sizeof(float),skl->nnodes*skl->ndims,f,swap);
  fread_sw(&j,sizeof(int),1,f,swap);

  skl->segdata = malloc((long)sizeof(double)*skl->nsegs*skl->nsegdata);
  fread_sw(&j,sizeof(int),1,f,swap);
  fread_sw(skl->segdata,sizeof(double),skl->nsegs*skl->nsegdata,f,swap);
  fread_sw(&j,sizeof(int),1,f,swap);

  skl->nodedata = malloc((long)sizeof(double)*skl->nnodes*skl->nnodedata);
  fread_sw(&j,sizeof(int),1,f,swap);
  fread_sw(skl->nodedata,sizeof(double),skl->nnodes*skl->nnodedata,f,swap);
  fread_sw(&j,sizeof(int),1,f,swap);

  skl->Node = malloc((long)sizeof(NDskl_node)*skl->nnodes);
  skl->Seg = malloc((long)sizeof(NDskl_seg)*skl->nsegs);

  fread_sw(&j,sizeof(int),1,f,swap);
  for(i=0;i<skl->nnodes;i++)
    {
      node = &skl->Node[i];
      fread_sw(&index,sizeof(int),1,f,swap);
      node->pos = &skl->nodepos[(long)index*skl->ndims];
      if (skl->nnodedata) node->data = &skl->nodedata[(long)index*skl->nnodedata];
      else node->data=NULL;
      fread_sw(&node->flags,sizeof(int),1,f,swap);
      fread_sw(&node->nnext,sizeof(int),1,f,swap);
      fread_sw(&node->type,sizeof(int),1,f,swap);
      fread_sw(&node->index,sizeof(int),1,f,swap);

      node->nsegs=malloc(sizeof(int)*node->nnext);
      node->Next=malloc(sizeof(NDskl_node *)*node->nnext);
      node->Seg=malloc(sizeof(NDskl_seg *)*node->nnext);
      fread_sw(node->nsegs,sizeof(int),node->nnext,f,swap);
      for (j=0;j<node->nnext;j++)
	{
	  fread_sw(&index,sizeof(int),1,f,swap);
	  if (index>=0) node->Next[j] = &skl->Node[index];
	  else node->Next[j] =NULL;

	  fread_sw(&index,sizeof(int),1,f,swap);

	   if (index>=0) node->Seg[j] = &skl->Seg[index];
	   else node->Seg[j] =NULL;
	}
    }
  fread_sw(&j,sizeof(int),1,f,swap);

  fread_sw(&j,sizeof(int),1,f,swap);
  for(i=0;i<skl->nsegs;i++)
    {
      seg = &skl->Seg[i];
      fread_sw(&index,sizeof(int),1,f,swap);
      seg->pos = &skl->segpos[(long)index*2*skl->ndims];
      if (skl->nsegdata) seg->data = &skl->segdata[(long)index*skl->nsegdata];
      else seg->data = NULL;
      fread_sw(seg->nodes,sizeof(int),2,f,swap);
      fread_sw(&seg->flags,sizeof(int),1,f,swap);
      fread_sw(&seg->index,sizeof(int),1,f,swap);
      
      fread_sw(&index,sizeof(int),1,f,swap);
      if (index<0) seg->Next=NULL;
      else seg->Next = &skl->Seg[index];
      
      fread_sw(&index,sizeof(int),1,f,swap);
      if (index<0) seg->Prev=NULL;
      else seg->Prev = &skl->Seg[index];
    
    }
  fread_sw(&j,sizeof(int),1,f,swap);

  fclose(f);
  printf (" done.\n");
  return skl;

}

int Save_ASCIIskel(NDskel *skl,const char *filename)
{
  long i,j,l;//,k
  //char tag[NDSKEL_DATA_STR_SIZE];
  FILE *f;
  long **filTab=NULL;
  int *filSegCount=NULL;
  NDskl_seg **filSegTab=NULL;
  int nfil;

  if (verbose>1) printf ("Saving %dD skeleton to ASCII file %s ...",skl->ndims,filename);fflush(0);

  int nodeDenId=getDataFieldID(skl,0,VALUE_TAG);
  int nodePairId=getDataFieldID(skl,0,PERSISTENCE_PAIR_TAG);
  
  if (nodeDenId==-1) fprintf(stderr,"\nWARNING : could not find '%s' information for nodes.\n",VALUE_TAG);
  if (nodePairId==-1) fprintf(stderr,"\nWARNING : could not find '%s' information for nodes.\n",PERSISTENCE_PAIR_TAG);

  f=fopen(filename,"w");
  //memset(tag,0,NDSKEL_DATA_STR_SIZE*sizeof(char));
  //strcpy(tag,NDSKEL_TAG);
  fprintf(f,"%s\n",NDSKEL_ASCII_TAG);
  fprintf(f,"%d\n",skl->ndims);
  fprintf(f,"#%s\n",skl->comment);
  fprintf(f,"BBOX [%g",skl->x0[0]);
  for (i=1;i<skl->ndims;i++) fprintf(f,",%g",skl->x0[i]);
  fprintf(f,"] [%g",skl->delta[0]);
  for (i=1;i<skl->ndims;i++) fprintf(f,",%g",skl->delta[i]);
  fprintf(f,"]\n");
  
  nfil=getNDskelFilTab(skl,&filSegTab,&filSegCount);

  filTab=(long**)calloc(skl->nnodes,sizeof(long*));
  for (i=0;i<skl->nnodes;i++) 
    if (skl->Node[i].nnext)
      filTab[i]=(long*)calloc(skl->Node[i].nnext,sizeof(long));
    else filTab[i]=NULL;
  

  NDskl_node *oldNode=NULL;
  NDskl_node *oldNext=NULL;
  int jNext=-1;
  int jNode=-1;
  for (i=0;i<nfil;i++)
    {
      NDskl_node *node=&skl->Node[filSegTab[i]->nodes[0]];
      NDskl_node *next=&skl->Node[filSegTab[i]->nodes[1]];

      if ((node == oldNode)&&(next==oldNext))
	{
	  jNode++;
	  jNext++;
	}
      else 
	{
	  jNode=0;
	  jNext=0;	  
	}

      for (;jNode<node->nnext;jNode++) 
	if (node->Next[jNode]==next) break;
      if (jNode==node->nnext) fprintf(stderr,"ERROR in Save_ASCIIskel: invalid skeleton file!\n");
      filTab[node->index][jNode]=i;

      for (;jNext<next->nnext;jNext++)
	if (next->Next[jNext]==node) break;
      if (jNext==next->nnext) fprintf(stderr,"ERROR in Save_ASCIIskel: invalid skeleton file!\n");
      filTab[next->index][jNext]=i;
      oldNode=node;
      oldNext=next;
      /*
      for (j=;j<node->nnext;j++) 
	if (node->Next[j]==next) break;
      if (j==node->nnext) fprintf(stderr,"ERROR in Save_ASCIIskel: invalid skeleton file!\n");
      filTab[node->index][j]=i;
      for (j=0;j<next->nnext;j++)
	if (next->Next[j]==node) break;
      if (j==next->nnext) fprintf(stderr,"ERROR in Save_ASCIIskel: invalid skeleton file!\n");
      filTab[next->index][j]=i;
*/
    }

  fprintf(f,"[CRITICAL POINTS]\n");
 
  fprintf(f,"%d\n",skl->nnodes);
  for (i=0;i<skl->nnodes;i++)
    {
      NDskl_node *node=&(skl->Node[i]);
      
      fprintf(f,"%d",node->type);
      for (j=0;j<skl->ndims;j++) fprintf(f," %g",node->pos[j]);

      if (nodeDenId>=0) fprintf(f," %g",node->data[nodeDenId]);
      else fprintf(f," 0");
      if (nodePairId>=0) fprintf(f," %d",(int)(node->data[nodePairId]));
      else fprintf(f," %d",(int)i);
      int boundary=0;
      if (node->flags&FLAG_NDNODE_INFINITE) boundary=3;
      else if (node->flags&FLAG_NDNODE_OUT) boundary=2;
      else if (node->flags&FLAG_NDNODE_BOUNDARY) boundary=1;
      fprintf(f," %d\n",boundary);
      fprintf(f," %d\n",node->nnext);
      
      for (j=0;j<node->nnext;j++)
	fprintf(f," %d %ld\n",node->Next[j]->index,(long)filTab[i][j]);
      
      //fprintf(f,"\n");
    }

  for (i=0;i<skl->nnodes;i++) free(filTab[i]);
  free(filTab);
 
  fprintf(f,"[FILAMENTS]\n");
  fprintf(f,"%d\n",nfil);
 
  for (i=0;i<nfil;i++)
    {
      NDskl_seg *seg=filSegTab[i];
      NDskl_node *node=&skl->Node[seg->nodes[0]];
      NDskl_node *next=&skl->Node[seg->nodes[1]];
     
      fprintf(f,"%ld %ld %ld\n",(long)node->index,(long)next->index,(long)filSegCount[i]+1);
      for (j=0;j<skl->ndims;j++) fprintf(f," %g",node->pos[j]);fprintf(f,"\n");
      if (seg->Next!=NULL)
	{
	  seg=seg->Next;
	  do {	    
	    for (j=0;j<skl->ndims;j++) 
	      fprintf(f," %g",seg->pos[j]);fprintf(f,"\n");
	    seg=seg->Next;
	  } while(seg!=NULL);
	}
      for (j=0;j<skl->ndims;j++) fprintf(f," %g",next->pos[j]);fprintf(f,"\n");
    }
 
  /*
  int reverse;
  for (i=0,k=0;i<skl->nnodes;i++)
    {
      NDskl_node *node=&(skl->Node[i]);
      for (j=0;j<node->nnext;j++) 
	{
	  NDskl_node *nextnode=node->Next[j];
	  if (node->Seg[j]->nodes[0]!=i) continue; 
	  //if (node->Next[j]->index>i) continue;

	  NDskl_seg *seg=node->Seg[j];
	  if (seg->nodes[0]==i) reverse=0; else reverse=1;
	  NDskl_seg *next = (reverse)?seg->Prev:seg->Next;
	  long delta1=(reverse)?skl->ndims:0;
	  long delta2=(reverse)?0:skl->ndims;

	  fprintf(f,"%d %d %d\n",node->index,nextnode->index,node->nsegs[j]+1); 	  
	  if (next==NULL)
	    {
	      for (l=0;l<skl->ndims;l++) fprintf(f," %g",seg->pos[l+delta1]);
	      for (l=0;l<skl->ndims;l++) fprintf(f," %g",seg->pos[l+delta2]);
	    }
	  else while (next!=NULL)
	    {
	      //printf("start (%ld) ... ",(long)seg);fflush(0);
	      int count=1;
	      for (l=0;l<skl->ndims;l++) fprintf(f," %g",seg->pos[l+delta1]);
	      while (next!=NULL)
		{
		  seg=next;
		  next = (reverse)?seg->Prev:seg->Next;
		  if (next!=NULL)
		    {
		      for (l=0;l<skl->ndims;l++) fprintf(f," %g",seg->pos[l+delta1]);
		      if (++count == 256) {fprintf(f,"\n");count=0;}
		    }
		}
	      NDskl_seg *prev=(reverse)?seg->Next:seg->Prev;
	      if (prev!=NULL)
		{
		  for (l=0;l<skl->ndims;l++) fprintf(f," %g",seg->pos[l+delta1]);
		  if (++count == 256) {fprintf(f,"\n");count=0;}
		  for (l=0;l<skl->ndims;l++) fprintf(f," %g",seg->pos[l+delta2]);
		}
	    }
	  fprintf(f,"\n");
	}
    }
  */
  fprintf(f,"[CRITICAL POINTS DATA]\n");
  fprintf(f,"%d\n",skl->nnodedata);
  for (i=0;i<skl->nnodedata;i++)
    fprintf(f,"%s\n",skl->nodedata_info[i]);
  for (i=0;i<skl->nnodes;i++)
    {	  
      for (j=0;j<skl->nnodedata-1;j++)
	fprintf(f,"%.7g ",skl->Node[i].data[j]);
      fprintf(f,"%.7g\n",skl->Node[i].data[j]);
    }

  fprintf(f,"[FILAMENTS DATA]\n");  
  int pid[skl->nsegdata*2];
  char tmp[skl->nsegdata][255];
  for (i=0,j=0;i<skl->nsegdata;i++)
    {
      char *w1=strstr(skl->segdata_info[i],"_p1");
      char *w2=strstr(skl->segdata_info[i],"_p2");
      
      if ((w1==NULL)&&(w2==NULL))
	{
	  pid[2*j]=pid[2*j+1]=i;	  
	  strcpy(tmp[j],skl->segdata_info[i]);
	  j++;
	}
      else if (w1!=NULL)
	{	  
	  strReplace(skl->segdata_info[i],tmp[j],"_p1","_p2");
	  int id=NDskel_SegDataIndex(skl,tmp[j]);
	  pid[2*j]=i;
	  pid[2*j+1]=id;	  
	  strcpy(tmp[j],skl->segdata_info[i]);
	  w1=strstr(tmp[j],"_p1");
	  (*w1)='\0';
	  j++;
	}
      else continue;       
    }
  int nsegdata=j;
  fprintf(f,"%d\n",nsegdata);
  for (i=0;i<nsegdata;i++)
    fprintf(f,"%s\n",tmp[i]);

  for (i=0;i<nfil;i++)
    {
      NDskl_seg *seg=filSegTab[i];
            
      for (l=0;l<nsegdata-1;l++)
	fprintf(f,"%.7g ",seg->data[pid[2*l]]);
      fprintf(f,"%.7g\n",seg->data[pid[2*l]]);

      if (seg->Next!=NULL)
	{	  
	  do {	    
	    seg=seg->Next;
	    for (l=0;l<nsegdata-1;l++)
	      fprintf(f,"%.7g ",seg->data[pid[2*l]]);
	    fprintf(f,"%.7g\n",seg->data[pid[2*l]]);
	  } while(seg->Next!=NULL);
	}
      
      for (l=0;l<nsegdata-1;l++)
	fprintf(f,"%.7g ",seg->data[pid[2*l+1]]);
      fprintf(f,"%.7g\n",seg->data[pid[2*l+1]]);
    }
  
  /*
  for (i=0,k=0;i<skl->nnodes;i++)
    {
      NDskl_node *node=&(skl->Node[i]);
      for (j=0;j<node->nnext;j++) 
	{
	  //NDskl_node *nextnode=node->Next[j];
	  //if (node->Next[j]->index>i) continue;
	  if (node->Seg[j]->nodes[0]!=i) continue;

	  NDskl_seg *seg=node->Seg[j];
	  if (seg->nodes[0]==i) reverse=0; else reverse=1;
	  NDskl_seg *next = (reverse)?seg->Prev:seg->Next;
	  long delta1=(reverse)?1:0;
	  long delta2=(reverse)?0:1;

	  if (next==NULL)
	    {
	      for (l=0;l<nsegdata-1;l++)
		fprintf(f,"%.7g ",seg->data[pid[2*l+delta1]]);
	      fprintf(f,"%.7g\n",seg->data[pid[2*l+delta1]]);
	      for (l=0;l<nsegdata-1;l++)
		fprintf(f,"%.7g ",seg->data[pid[2*l+delta2]]);
	      fprintf(f,"%.7g\n",seg->data[pid[2*l+delta2]]);
	    }
	  else while (next!=NULL)
	    {
	      for (l=0;l<nsegdata-1;l++)
		fprintf(f,"%.7g ",seg->data[pid[2*l+delta1]]);
	      fprintf(f,"%.7g\n",seg->data[pid[2*l+delta1]]);
	      while (next!=NULL)
		{
		  seg=next;
		  next = (reverse)?seg->Prev:seg->Next;
		  if (next!=NULL)
		    {
		      for (l=0;l<nsegdata-1;l++)
			fprintf(f,"%.7g ",seg->data[pid[2*l+delta1]]);
		      fprintf(f,"%.7g\n",seg->data[pid[2*l+delta1]]);
		    }
		}
	      NDskl_seg *prev=(reverse)?seg->Next:seg->Prev;
	      if (prev!=NULL)
		{
		  for (l=0;l<nsegdata-1;l++)
		    fprintf(f,"%.7g ",seg->data[pid[2*l+delta1]]);
		  fprintf(f,"%.7g\n",seg->data[pid[2*l+delta1]]);
		  for (l=0;l<nsegdata-1;l++)
		    fprintf(f,"%.7g ",seg->data[pid[2*l+delta2]]);
		  fprintf(f,"%.7g\n",seg->data[pid[2*l+delta2]]);
		}
	    }
	}
    }
  */
  fclose(f);
  freeNDskelFilTab(&filSegTab,&filSegCount);

  if (verbose>1) printf(" done.\n");
  return 0;
}

int Save_NDskel(NDskel *skl,const char *filename)
{
  int i,j,k;
  char tag[NDSKEL_DATA_STR_SIZE];
  char dummy[160];
  FILE *f;
  int index;

  NDskl_node *node;
  NDskl_seg *seg;
  
  if (verbose>1) printf ("Saving %dD skeleton to file %s ...",skl->ndims,filename);fflush(0);
  memset(dummy,0,160);
  memset(tag,0,NDSKEL_DATA_STR_SIZE*sizeof(char));
  strcpy(tag,NDSKEL_TAG);
  i=16;
  
  f=fopen(filename,"w");
  fwrite(&i,sizeof(int),1,f);
  fwrite(tag,sizeof(char),16,f);
  fwrite(&i,sizeof(int),1,f);

  j=sizeof(char)*(80+NDSKEL_DATA_STR_SIZE*(skl->nsegdata+skl->nnodedata))+
    sizeof(int)*(4+1+NDSKEL_MAX_DIMS)+
    sizeof(double)*(2*NDSKEL_MAX_DIMS)+
    sizeof(char)*NDSKEL_DATA_STR_SIZE*(skl->nnodedata+skl->nsegdata);

  fwrite(&j,sizeof(int),1,f);
  fwrite(skl->comment,sizeof(char),80,f);
  fwrite(&skl->ndims,sizeof(int),1,f);
  fwrite(skl->dims,sizeof(int),skl->ndims,f);
  if (skl->ndims<NDSKEL_MAX_DIMS) fwrite(dummy,sizeof(int),NDSKEL_MAX_DIMS-skl->ndims,f);
  fwrite(skl->x0,sizeof(double),skl->ndims,f);
  if (skl->ndims<NDSKEL_MAX_DIMS) fwrite(dummy,sizeof(double),NDSKEL_MAX_DIMS-skl->ndims,f);
  fwrite(skl->delta,sizeof(double),skl->ndims,f);
  if (skl->ndims<NDSKEL_MAX_DIMS) fwrite(dummy,sizeof(double),NDSKEL_MAX_DIMS-skl->ndims,f);
  
  fwrite(&skl->nsegs,sizeof(int),1,f);
  fwrite(&skl->nnodes,sizeof(int),1,f);
  fwrite(&skl->nsegdata,sizeof(int),1,f);
  fwrite(&skl->nnodedata,sizeof(int),1,f);
  fwrite(&j,sizeof(int),1,f);

  j=skl->nsegdata*NDSKEL_DATA_STR_SIZE*sizeof(char);
  if (skl->nsegdata) fwrite(&j,sizeof(int),1,f);
  for (i=0;i<skl->nsegdata;i++)
    {
      memset(tag,0,NDSKEL_DATA_STR_SIZE*sizeof(char));
      strcpy(tag,skl->segdata_info[i]);
      fwrite(tag,sizeof(char),NDSKEL_DATA_STR_SIZE,f);
    }
  if (skl->nsegdata) fwrite(&j,sizeof(int),1,f);

  j=skl->nnodedata*NDSKEL_DATA_STR_SIZE*sizeof(char);
  if (skl->nnodedata) fwrite(&j,sizeof(int),1,f);
  for (i=0;i<skl->nnodedata;i++)
    {
      memset(tag,0,NDSKEL_DATA_STR_SIZE*sizeof(char));
      strcpy(tag,skl->nodedata_info[i]);
      fwrite(tag,sizeof(char),NDSKEL_DATA_STR_SIZE,f);
    }
  if (skl->nnodedata) fwrite(&j,sizeof(int),1,f);
    
  j=sizeof(float)*skl->nsegs*2*skl->ndims;
  fwrite(&j,sizeof(int),1,f);
  fwrite(skl->segpos,sizeof(float),skl->nsegs*2*skl->ndims,f);
  fwrite(&j,sizeof(int),1,f);
  
  j=sizeof(float)*skl->nnodes*2*skl->ndims;
  fwrite(&j,sizeof(int),1,f);
  fwrite(skl->nodepos,sizeof(float),skl->nnodes*skl->ndims,f);
  fwrite(&j,sizeof(int),1,f);
  /*
  for (i=0;i<4*skl->ndims*2;i++)
    {
      if ((i%(2*skl->ndims))==0) printf ("\n");
      printf ("%f ",skl->segpos[i]);
      
    }
  printf ("\n");
  */
  /*
  for (i=0;i<skl->nnodes;i++)
    if (skl->Node[i].index != i)
      printf("Node: %d != %d\n",skl->Node[i].index,i);

  for (i=0;i<skl->nsegs;i++)
    if (skl->Seg[i].index != i)
      printf("SEG: %d != %d\n",skl->Node[i].index,i);
  */
  j=sizeof(double)*skl->nsegs*skl->nsegdata;
  fwrite(&j,sizeof(int),1,f);
  fwrite(skl->segdata,sizeof(double),skl->nsegs*skl->nsegdata,f);
  fwrite(&j,sizeof(int),1,f);
  
  j=sizeof(double)*skl->nnodes*skl->nnodedata;
  fwrite(&j,sizeof(int),1,f);
  fwrite(skl->nodedata,sizeof(double),skl->nnodes*skl->nnodedata,f);
  fwrite(&j,sizeof(int),1,f);
  j=0;
  fwrite(&j,sizeof(int),1,f);
  for(i=0;i<skl->nnodes;i++)
    {
      node = &skl->Node[i];
      index = (long)(node->pos-skl->nodepos)/(skl->ndims);
      fwrite(&index,sizeof(int),1,f);
      fwrite(&node->flags,sizeof(int),1,f);
      fwrite(&node->nnext,sizeof(int),1,f);
      fwrite(&node->type,sizeof(int),1,f);
      fwrite(&node->index,sizeof(int),1,f);
      fwrite(node->nsegs,sizeof(int),node->nnext,f);
      for (j=0;j<node->nnext;j++)
	{
	    if (node->Next[j] != NULL)
		fwrite(&node->Next[j]->index,sizeof(int),1,f);
	    else {k=-1;fwrite(&k,sizeof(int),1,f);}
	    
	    if (node->Seg[j] != NULL)
		fwrite(&node->Seg[j]->index,sizeof(int),1,f);
	    else {k=-1;fwrite(&k,sizeof(int),1,f);}
	}
    }
  fwrite(&j,sizeof(int),1,f);

  k=-1;
  fwrite(&j,sizeof(int),1,f);
  for(i=0;i<skl->nsegs;i++)
    {
      seg = &skl->Seg[i];
      index = (long)(seg->pos-skl->segpos)/(skl->ndims*2);
      fwrite(&index,sizeof(int),1,f);
      fwrite(seg->nodes,sizeof(int),2,f);
      fwrite(&seg->flags,sizeof(int),1,f);
      fwrite(&seg->index,sizeof(int),1,f);
      if (seg->Next==NULL) fwrite(&k,sizeof(int),1,f);
      else fwrite(&seg->Next->index,sizeof(int),1,f);
      if (seg->Prev==NULL) fwrite(&k,sizeof(int),1,f);
      else fwrite(&seg->Prev->index,sizeof(int),1,f);
    }
  fwrite(&j,sizeof(int),1,f);

  fclose(f);
  if (verbose>1) printf (" done.\n");
  return 0;

}

NDskel *NDskel_Subregion(NDskel *skl_p, double margin, int **seglist, long *nsegs, int inplace)
{
    int j,k;
    NDskel *skl=skl_p;
    int *segindex;
    int *nodeindex;
    int curindex;
    const long Nalloc = 1000;

    if (!inplace) skl=Copy_NDskel(skl_p);
  
    segindex = malloc(sizeof(int)*skl->nsegs);
    nodeindex = malloc(sizeof(int)*skl->nnodes);
	  
    curindex=0;
    for (j=0;j<skl->nnodes;j++)
    {
	for (k=0;k<skl->ndims;k++)
	{
	    if ((skl->Node[j].pos[k]<skl->x0[k]+margin)||(skl->Node[j].pos[k]>skl->x0[k]+skl->delta[k]-margin))
		break;
	}
	if (k==skl->ndims)
	{
	    //this node is inside
	    nodeindex[j]=curindex++;
	}
	else
	{
	    //this node is outside
	    nodeindex[j]=-1;
	}
    }
    //printf ("nnodes=%d\n",curindex);
    
    curindex=0;(*nsegs)=0;
    for (j=0;j<skl->nsegs;j++)
    {
	int ok;
	
	ok=0;
	for (k=0;k<skl->ndims;k++)
	{
	    if (((skl->Seg[j].pos[k]<skl->x0[k]+margin)&&(skl->Seg[j].pos[k+skl->ndims]<skl->x0[k]+margin))||
		((skl->Seg[j].pos[k]>skl->x0[k]+skl->delta[k]-margin)&&(skl->Seg[j].pos[k+skl->ndims]>skl->x0[k]+skl->delta[k]-margin)))
		break;
	    
	    if (((skl->Seg[j].pos[k]<skl->x0[k]+margin)&&(skl->Seg[j].pos[k+skl->ndims]>skl->x0[k]+margin))||
		((skl->Seg[j].pos[k+skl->ndims]<skl->x0[k]+margin)&&(skl->Seg[j].pos[k]>skl->x0[k]+margin)))
		ok=1;
	    
	    if (((skl->Seg[j].pos[k]<skl->x0[k]+skl->delta[k]-margin)&&(skl->Seg[j].pos[k+skl->ndims]>skl->x0[k]+skl->delta[k]-margin))||
		((skl->Seg[j].pos[k+skl->ndims]<skl->x0[k]+skl->delta[k]-margin)&&(skl->Seg[j].pos[k]>skl->x0[k]+skl->delta[k]-margin)))
		ok=1;
	}
	
	if (k==skl->ndims)
	{
	    // this segment is inside
	    if (ok)
	    {
		// this segment is on the edge
		if ((*nsegs)%Nalloc == 0) 
		    (*seglist) = (int *) realloc((*seglist),sizeof(int)*((*nsegs)+Nalloc));

		segindex[j]=curindex++;
		(*seglist)[(*nsegs)++] = segindex[j];
		//nsegs++;
	    }
	    else
	    {
		// this segment is fine 
		segindex[j]=curindex++;
	    }		  
	}
	else
	{
	    // this egment is outside, should be removed
	    segindex[j]=-1;
	}
    }
    //printf ("nsegs=%d\n",curindex);
    // fix the nodes ...
    for (j=0;j<skl->nnodes;j++)
    {
	NDskl_node *curnode;
	long tmp;

	if (nodeindex[j]<0) 
	{
	    curnode = &skl->Node[j];
	    free(curnode->Next);
	    free(curnode->Seg);
	    free(curnode->nsegs);
	    continue;
	}
	
	for (k=0;k<skl->Node[j].nnext;k++)
	{
	    tmp = nodeindex[skl->Node[j].Next[k]->index];
	    if (tmp>=0) skl->Node[j].Next[k] = &skl->Node[tmp];
	    else skl->Node[j].Next[k]=NULL;

	    tmp = segindex[skl->Node[j].Seg[k]->index];
	    if (tmp>=0)  skl->Node[j].Seg[k] = &skl->Seg[tmp];
	    else printf ("No one should ever read this ...\n");
	}  
	
	curnode = &skl->Node[nodeindex[j]];
		
	if (nodeindex[j]!=j)
	{
/*
	    free(curnode->Next);
	    free(curnode->Seg);
	    free(curnode->nsegs);
*/	
	    memcpy(curnode,&skl->Node[j],sizeof(NDskl_node));
	    memcpy(&skl->nodepos[nodeindex[j]*skl->ndims],curnode->pos,skl->ndims*sizeof(float));
	    memcpy(&skl->nodedata[nodeindex[j]*skl->nnodedata],curnode->data,skl->nnodedata*sizeof(double));
	    curnode->data = &skl->nodedata[nodeindex[j]*skl->nnodedata];
	    curnode->pos = &skl->nodepos[nodeindex[j]*skl->ndims];
	}
    }
    j=skl->nnodes-1;do {skl->nnodes = 1+nodeindex[j--];} while(skl->nnodes<=0);
    

    for (j=0;j<skl->nsegs;j++)
    {
	NDskl_seg *curseg = &skl->Seg[segindex[j]];
	long tmp;
	if (segindex[j]<0) continue;
	
	curseg = &skl->Seg[segindex[j]];

	if (skl->Seg[j].Next!=NULL)
	{
	    tmp = segindex[skl->Seg[j].Next->index];
	    if (tmp>=0) skl->Seg[j].Next = &skl->Seg[tmp];
	    else skl->Seg[j].Next=NULL;
	}
	if (skl->Seg[j].Prev!=NULL)
	{
	    tmp = segindex[skl->Seg[j].Prev->index];
	    if (tmp>=0) skl->Seg[j].Prev = &skl->Seg[tmp];
	    else skl->Seg[j].Prev=NULL;
	}

	if (segindex[j]!=j)
	    memcpy(curseg,&skl->Seg[j],sizeof(NDskl_seg));
	
	curseg->nodes[0]=nodeindex[curseg->nodes[0]];
	curseg->nodes[1]=nodeindex[curseg->nodes[1]];
	
	memcpy(&skl->segpos[segindex[j]*2*skl->ndims],curseg->pos,skl->ndims*2*sizeof(float));
	memcpy(&skl->segdata[segindex[j]*skl->nsegdata],curseg->data,skl->nsegdata*sizeof(double));
	curseg->data = &skl->segdata[segindex[j]*skl->nsegdata];
	curseg->pos = &skl->segpos[segindex[j]*skl->ndims*2];
    }
    j=skl->nsegs-1;do {skl->nsegs = 1+segindex[j--];} while(skl->nsegs<=0);
    
    free(segindex);
    free(nodeindex);
    
    for (j=0;j<skl->nsegs;j++) skl->Seg[j].index = j;
    for (j=0;j<skl->nnodes;j++) skl->Node[j].index = j;
    
    skl->segpos = realloc(skl->segpos,skl->nsegs*2*skl->ndims*sizeof(float));
    skl->segdata = realloc(skl->segdata,skl->nsegs*skl->nsegdata*sizeof(double));
    skl->nodepos = realloc(skl->nodepos,skl->nnodes*skl->ndims*sizeof(float));
    skl->nodedata = realloc(skl->nodedata,skl->nnodes*skl->nnodedata*sizeof(double));
    
    return skl;
}
/*
long getNDskelFilTabWithDuplicates(NDskel *skl, NDskl_seg ***filTab, int **filSize)
{
  long nfil=0;  
  long i,j;

  for (i=0,nfil=0;i<skl->nnodes;i++)
    {
      NDskl_node *node=&(skl->Node[i]);
      for (j=0;j<node->nnext;j++) 
	{
	  //if (node->Seg[j]->nodes[0]==i) 
	    {
	      nfil++;
	      //(*filTab)=(NDskl_seg**) realloc((*filTab),sizeof(NDskl_seg*)*(nfil+1));
	      //(*filTab)[nfil++]=node->Seg[j];
	    }
	}
    }

  (*filTab)=(NDskl_seg**) realloc((*filTab),sizeof(NDskl_seg*)*(nfil));

  for (i=0,nfil=0;i<skl->nnodes;i++)
    {
      NDskl_node *node=&(skl->Node[i]);
      for (j=0;j<node->nnext;j++) 
	{
	  //if (node->Seg[j]->nodes[0]==i) 
	    {
	      //(*filTab)=(NDskl_seg**) realloc((*filTab),sizeof(NDskl_seg*)*(nfil+1));
	      (*filTab)[nfil++]=node->Seg[j];
	    }
	}
    }

  if (filSize!=NULL)
    {
      (*filSize)=realloc((*filSize),nfil*sizeof(int));
      for (i=0;i<nfil;i++)
	{
	  NDskl_seg *seg=(*filTab)[i];
	  if (seg==NULL) (*filSize)[i]=0;
	  else
	    {
	      j=1;
	      while (seg->Next!=NULL)
		{
		  seg=seg->Next;
		  j++;
		}
	      (*filSize)[i]=j;
	    }
	}
    }

  return nfil;
}
*/
long getNDskelFilTab(NDskel *skl, NDskl_seg ***filTab, int **filSize)
{
  long nfil=0;  
  long i,j;

  for (i=0,nfil=0;i<skl->nnodes;i++)
    {
      NDskl_node *node=&(skl->Node[i]);
      for (j=0;j<node->nnext;j++) 
	{	  
	  if (node->Seg[j]->nodes[0]==i) 
	    {
	      nfil++;
	      //(*filTab)=(NDskl_seg**) realloc((*filTab),sizeof(NDskl_seg*)*(nfil+1));
	      //(*filTab)[nfil++]=node->Seg[j];
	    }
	}
    }

  (*filTab)=(NDskl_seg**) realloc((*filTab),sizeof(NDskl_seg*)*(nfil));

  for (i=0,nfil=0;i<skl->nnodes;i++)
    {
      NDskl_node *node=&(skl->Node[i]);
      for (j=0;j<node->nnext;j++) 
	{
	  if (node->Seg[j]->nodes[0]==i) 
	    {
	      //(*filTab)=(NDskl_seg**) realloc((*filTab),sizeof(NDskl_seg*)*(nfil+1));
	      (*filTab)[nfil++]=node->Seg[j];
	    }
	}
    }

  if (filSize!=NULL)
    {
      (*filSize)=realloc((*filSize),nfil*sizeof(int));
      for (i=0;i<nfil;i++)
	{
	  NDskl_seg *seg=(*filTab)[i];
	  j=1;
	  while(seg->Next!=NULL)
	    {
	      seg=seg->Next;
	      j++;
	    }
	  (*filSize)[i]=j;
	}
    }

  return nfil;
}

long getNDskelFilTabInfo(NDskel *skl, NDskl_seg **filTab, int nFil, char ***dataName, double ***data)
{
  long i,j;
  double len;
  int ndims=skl->ndims;
  int filTID=getDataFieldID(skl,1,TYPE_TAG);
  int ndata=0;

  (*dataName)=(char**)realloc(*dataName,5*sizeof(char*));
  for (i=0;i<5;i++) (*dataName)[i]=(char*)malloc(50*sizeof(char));
  (*data)=(double**)realloc(*data,nFil*sizeof(double*));
  for (i=0;i<nFil;i++) (*data)[i]=(double*)malloc(5*sizeof(double));

  strcpy((*dataName)[ndata++],ARC_ID_TAG);
  strcpy((*dataName)[ndata++],TYPE_TAG);
  strcpy((*dataName)[ndata++],DOWN_TAG(INDEX_TAG));
  strcpy((*dataName)[ndata++],UP_TAG(INDEX_TAG));
  strcpy((*dataName)[ndata++],LENGTH_TAG);
  
  for (i=0;i<nFil;i++)
    {
      NDskl_seg *seg=filTab[i];
      NDskl_node *n1=&skl->Node[seg->nodes[0]];
      NDskl_node *n2=&skl->Node[seg->nodes[1]];
      
      (*data)[i][0]=i; 
      //(*data)[i][1]=(n1->type<n2->type)?n1->type:n2->type;
      if (filTID<0) (*data)[i][1]=(n1->type<n2->type)?n1->type:n2->type;	
      else (*data)[i][1]=seg->data[filTID];
      (*data)[i][2]=(n1->type<n2->type)?n1->index:n2->index;
      (*data)[i][3]=(n1->type>n2->type)?n1->index:n2->index;
      
      do {
	for (j=0,len=0;j<ndims;j++) 
	  {
	    double tmp=(seg->pos[ndims+j]-seg->pos[j]);
	    len+=tmp*tmp;
	  }
	(*data)[i][4]+=sqrt(len);	
	seg=seg->Next;
      } while (seg!=NULL);
    }

  return ndata;
}

void freeNDskelFilTab(NDskl_seg ***filTab, int **filSize)
{
  if (filTab!=NULL)
    {
      free(*filTab);
      (*filTab)=NULL;
    }
  
  if (filSize!=NULL)
    {      
      free(*filSize);      
      (*filSize)=NULL;
    } 
}

void freeNDskelFilTabInfo(char ***dataName, double ***data, int nData, int nFil)
{
  long i;

  if (dataName!=NULL)
    {
      for (i=0;i<nData;i++)
	free((*dataName)[i]);
      free(*dataName);
      (*dataName)=NULL;
    }

  if (data!=NULL)
    {
      for (i=0;i<nFil;i++)
	 free((*data)[i]);
      free(*data);
      (*data)=NULL;
    }
}

int printNDskelStat(NDskel *skl, int dec)
{
  long i;
  NDskl_seg **filTab=NULL;
  long nfil=getNDskelFilTab(skl,&filTab,NULL);
  freeNDskelFilTab(&filTab,NULL);

  for (i=0;i<dec;i++) printf (" ");
  for (i=0;i<80;i++) if ((i=='\0')||(i=='\n')) break;
  if (i!=80) printf("comment: '%s'\n",skl->comment);
  for (i=0;i<dec;i++) printf (" ");
  printf ("Skeleton has %ld nodes and %ld segs.\n",(long)skl->nnodes, (long)skl->nsegs);
  for (i=0;i<dec;i++) printf (" ");
  printf ("Number of filaments: %ld.\n",nfil);
  for (i=0;i<dec;i++) printf (" ");
  printf("Bounding box: x0=[%g",skl->x0[0]);
  for (i=1;i<skl->ndims;i++) printf(",%g",skl->x0[i]);
  printf("],\n");
  for (i=0;i<dec;i++) printf (" ");
  printf("              delta=[%g",skl->delta[0]);
  for (i=1;i<skl->ndims;i++) printf(",%g",skl->delta[i]);
  printf("].\n");
  for (i=0;i<dec;i++) printf (" ");
  printf("Available node fields: ");
  if (skl->nnodedata) printf("'%s'",skl->nodedata_info[0]);
  for (i=1;i<skl->nnodedata;i++) 
    {
      if ((i%3) == 0)
	{
	  long j;
	  printf("\n");
	  for (j=0;j<dec+strlen("Available node fields: ");j++) printf (" ");
	  printf("'%s'",skl->nodedata_info[i]);
	}
      else printf(", '%s'",skl->nodedata_info[i]);     
    }
  printf(".\n");

  for (i=0;i<dec;i++) printf (" ");
  printf("Available seg fields: ");
  if (skl->nsegdata) printf("'%s'",skl->segdata_info[0]);
  for (i=1;i<skl->nsegdata;i++) 
    {
      if ((i%3) == 0)
	{
	  long j;
	  printf("\n");
	  for (j=0;j<dec+strlen("Available seg fields: ");j++) printf (" ");
	  printf("'%s'",skl->segdata_info[i]);
	}
      else printf(", '%s'",skl->segdata_info[i]);     
    }
  printf(".\n");  
  
  return 0;
}

#ifdef HAVE_CFITS_IO
void NDSKEL_FITS_printerror( int status)
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

long NDskel_pos2pix(NDskel *skl,float *pos, long *dims)
{
  long i;
  double p;
  long pix=0;
  long fac=1;
  for (i=0;i<skl->ndims;i++)
    {
      p=(pos[i]-skl->x0[i])/(skl->delta[i]) * dims[i];
      if (p<0) p=0;
      else if (p==dims[i]) p=dims[i]-1;
      else if (p>dims[i]) p-=dims[i];
      pix += ((long)p)*fac;
      fac *= dims[i];
    }
  return pix;
}

int Save_FITSskel(NDskel *skl,const char *filename, long *naxes_p)
{
  fitsfile *fptr;       
  int status;
  long i,j;
  long nelements;
  int *array;
  int bitpix=32;         
  long naxis=skl->ndims;        
  long naxes[skl->ndims]; 
  status = 0;       

  if (naxes_p==NULL)
    for (i=0;i<skl->ndims;i++) naxes[i] = (long)skl->delta[i];
  else
    for (i=0;i<skl->ndims;i++) naxes[i] = naxes_p[i];

  if (verbose>1) 
    {
      printf("Sampling skeleton to a [%ld",naxes[0]);
      for (i=1;i<skl->ndims;i++) printf(",%ld",naxes[i]);
      printf("] FITS image ...");fflush(0);
    }

  if (fits_create_file(&fptr, filename, &status)) 
    NDSKEL_FITS_printerror( status );

  if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
    NDSKEL_FITS_printerror( status );

  for (i=0,nelements=1;i<skl->ndims;i++) nelements*=naxes[i];
  array=(int*)calloc(nelements,sizeof(int));

  
  NDskl_seg **filTab=NULL;
  long nfil=getNDskelFilTab(skl,&filTab,NULL);  
  float pos[skl->ndims];
  for (i=0;i<nfil;i++)
    {
      NDskl_seg *seg=filTab[i];
      if (seg==NULL) continue;
      do {
	for (j=0;j<skl->ndims;j++) pos[j]=0.5*(seg->pos[j]+seg->pos[j+skl->ndims]);
	array[NDskel_pos2pix(skl,seg->pos,naxes)] = i+1;
	array[NDskel_pos2pix(skl,seg->pos+skl->ndims,naxes)] = i+1;
	array[NDskel_pos2pix(skl,pos,naxes)] = i+1;
	seg=seg->Next;
      } while (seg!=NULL);
      }
  freeNDskelFilTab(&filTab,NULL);
  
  /*
  float pos[skl->ndims];
  for (j=0;j<skl->nsegs;++j)
    {
      NDskl_seg *seg=&skl->Seg[j];
      long i=seg->nodes[0]+1;
      //if (seg->nodes[0]<0) printf("node=%d\n",seg->nodes[0]);
      int k;
      for (k=0;k<skl->ndims;k++) pos[k]=0.5*(seg->pos[k]+seg->pos[k+skl->ndims]);

      long pix=NDskel_pos2pix(skl,seg->pos,naxes);
      if ((array[pix]==0)||(array[pix]>i))
	array[pix] = i;
      
      pix=NDskel_pos2pix(skl,seg->pos+skl->ndims,naxes);
      if ((array[pix]==0)||(array[pix]>i))
	array[pix] = i;

      pix=NDskel_pos2pix(skl,pos,naxes);
      if ((array[pix]==0)||(array[pix]>i))
	array[pix] = i;
    }
  */
  if (verbose>1) printf(" done.\nSaving file %s ...",filename);
  if ( fits_write_img(fptr, TINT, 1, nelements, array, &status) )
    NDSKEL_FITS_printerror( status );

  if ( fits_write_key(fptr, TDOUBLE, "X0", &skl->x0[0],NULL, &status) )
    NDSKEL_FITS_printerror( status );
  if (skl->ndims>=2)
    if ( fits_write_key(fptr, TDOUBLE, "Y0", &skl->x0[1],NULL, &status) )
      NDSKEL_FITS_printerror( status );
  if (skl->ndims>=3)
    if ( fits_write_key(fptr, TDOUBLE, "Z0", &skl->x0[2],NULL, &status) )
      NDSKEL_FITS_printerror( status );

  if ( fits_write_key(fptr, TDOUBLE, "deltaX", &skl->delta[0],NULL, &status) )
    NDSKEL_FITS_printerror( status );
  if (skl->ndims>=2)
    if ( fits_write_key(fptr, TDOUBLE, "deltaY", &skl->delta[1],NULL, &status) )
      NDSKEL_FITS_printerror( status );	
  if (skl->ndims>=3)
    if ( fits_write_key(fptr, TDOUBLE, "deltaZ", &skl->delta[2],NULL, &status) )
      NDSKEL_FITS_printerror( status );

  fits_close_file(fptr, &status);           
  NDSKEL_FITS_printerror( status );

  if (verbose>1) printf(" done.\n");

  return 0;

}
/*
int Save_FITSskel_layered(NDskel *skl,const char *filename, long *naxes_p)
{
  fitsfile *fptr;       
  int status;
  long i,j;
  long nelements;
  int *array;
  int bitpix=32;         
  long naxis=skl->ndims+1;        
  long naxes[skl->ndims+1]; 
  status = 0;       

  if (naxes_p==NULL)
    for (i=0;i<skl->ndims;i++) naxes[i] = (long)skl->delta[i];
  else
    for (i=0;i<skl->ndims;i++) naxes[i] = naxes_p[i];  

  NDskl_seg **filTab=NULL;
  long nfil=getNDskelFilTab(skl,&filTab,NULL);  
  naxes_p[skl->ndims]=nfil;

  if (verbose>1) 
    {
      printf("Sampling skeleton to a [%ld",naxes[0]);
      for (i=1;i<skl->ndims+1;i++) printf(",%ld",naxes[i]);
      printf("] FITS image ...");fflush(0);
    }

  if (fits_create_file(&fptr, filename, &status)) 
    NDSKEL_FITS_printerror( status );

  if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
    NDSKEL_FITS_printerror( status );

  for (i=0,nelements=1;i<skl->ndims;i++) nelements*=naxes[i];
  array=(int*)calloc(nelements,sizeof(int));
  long filPadding=nelements/nfil;

  float pos[skl->ndims];
  for (i=0;i<nfil;i++)
    {
      NDskl_seg *seg=filTab[i];
      do {
	for (j=0;j<skl->ndims;j++) pos[j]=0.5*(seg->pos[j]+seg->pos[j+skl->ndims]);
	array[filPadding*i+NDskel_pos2pix(skl,seg->pos,naxes)] = i;
	array[filPadding*i+NDskel_pos2pix(skl,seg->pos+skl->ndims,naxes)] = i;
	array[filPadding*i+NDskel_pos2pix(skl,pos,naxes)] = i;
	seg=seg->Next;
      } while (seg!=NULL);
    }
  freeNDskelFilTab(&filTab,NULL);
  if (verbose>1) printf(" done.\nSaving file %s ...",filename);
  if ( fits_write_img(fptr, TINT, 1, nelements, array, &status) )
    NDSKEL_FITS_printerror( status );

  if ( fits_write_key(fptr, TDOUBLE, "X0", &skl->x0[0],NULL, &status) )
    NDSKEL_FITS_printerror( status );
  if (skl->ndims>=2)
    if ( fits_write_key(fptr, TDOUBLE, "Y0", &skl->x0[1],NULL, &status) )
      NDSKEL_FITS_printerror( status );
  if (skl->ndims>=3)
    if ( fits_write_key(fptr, TDOUBLE, "Z0", &skl->x0[2],NULL, &status) )
      NDSKEL_FITS_printerror( status );

  if ( fits_write_key(fptr, TDOUBLE, "deltaX", &skl->delta[0],NULL, &status) )
    NDSKEL_FITS_printerror( status );
  if (skl->ndims>=2)
    if ( fits_write_key(fptr, TDOUBLE, "deltaY", &skl->delta[1],NULL, &status) )
      NDSKEL_FITS_printerror( status );	
  if (skl->ndims>=3)
    if ( fits_write_key(fptr, TDOUBLE, "deltaZ", &skl->delta[2],NULL, &status) )
      NDSKEL_FITS_printerror( status );

  fits_close_file(fptr, &status);           
  NDSKEL_FITS_printerror( status );

  if (verbose>1) printf(" done.\n");

  return 0;

}
*/
#endif
