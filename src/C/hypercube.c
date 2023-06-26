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
#include "mytypes.h"
#include "hypercube.h"

#include <stdlib.h>
#include <stdio.h>

// each face is represented by a tag, composed of a mask M and a displacement D 
// tag = M+(D<<ndims);
// in M and D, the value of the Nth bit represents a displacement along the Nth direction

// if (M&(1<<N)) == true then no displacement along the Nth direction is allowed
// for instance, in 3D, if M=001b, one is allowed to move along Z and Y, which defines the face
// in the (ZY) plane, starting at origin. M=011b represents the segments along Z starting at origin ...

// if (D&(1<<N)) == true then the face defined by M is displaced along the Nth direction.
// for instance, in 3D,the tag=001001b=9 (M=001b and D=001b) is the face in (ZY) plane with origin moved along X.
// tag=010011b=19 (M=011b and D=010b) is the segment along Z, with origin displaced along Y.

#define COUNTBITS(val,result,tmp) {tmp=val;result=0;while (tmp) {result++;tmp=(((tmp)!=0)?((tmp)&((tmp)-1)):0);}}


inline int hCNP(int n, int p)
{
  unsigned long int i;
  unsigned long int up=1;
  unsigned long int down=1;

  if (n==p) return 1;

  for (i=n;i>p;i--)
    up*=i;
  for (i=n-p;i>0;i--)
    down*=i;

  return up/down;
}

hypercube_type *GenerateHypercube(int ndims)
{
  hypercube_type *h=NULL;
  nface_type *nface=NULL;
  nface_type *nface2=NULL;
  nfacelist_type *nfacelist=NULL;

  INT i,j,k,n;
  INT tmp;
  INT mask;
  INT mask2;
  int nbits;
  
  h=calloc (1,sizeof(hypercube_type));
  h->ndims=ndims;
  h->NFaceList = calloc(ndims+1,sizeof(nfacelist_type));
  h->faceFromTag = calloc(1<<(2*ndims),sizeof(nface_type*));
  // change < to <= to also consider 1-faces
  for (i=0;i<=ndims;i++)
    {
      h->NFaceList[i].ndims = ndims-i;
      h->NFaceList[i].nFacesPerCube = hCNP(ndims,ndims-i);
      h->NFaceList[i].nfaces = hCNP(ndims,ndims-i)*(1<<i);
      h->NFaceList[i].NFace = calloc(h->NFaceList[i].nfaces,sizeof(nface_type));
      

      for (j=0;j<h->NFaceList[i].nfaces;j++)
	{
	  nface = &(h->NFaceList[i].NFace[j]);
	  nface->ndims = ndims-i;
	  nface->nvertice = 0; //nvertice = 1<<(ndims-i)
	  nface->vertex = calloc(1<<(ndims-i),sizeof(int));
	  nface->tag=-1;
	  nface->nmothers = 0;
	  nface->ndaughters = 0;
	  nface->Mothers = calloc(i,sizeof(nface_type *));
	  //remove condition to consider 1-faces
	  //if (i!=ndims-1)
	  nface->Daughters = calloc(2*(ndims-i),sizeof(nface_type *));	  	 
	}
    }

  h->NFaceList[0].mother=NULL;
  h->NFaceList[0].daughter=&(h->NFaceList[1]);
  for (j=1;j<i;j++)
    {
      h->NFaceList[j].mother=&(h->NFaceList[j-1]);
      h->NFaceList[j].daughter=&(h->NFaceList[j+1]);
    }
  h->NFaceList[i-1].mother=&(h->NFaceList[i-2]);
  h->NFaceList[i-1].daughter=NULL;
  
  // i is a given vertex
  for (i=0;i<(1<<ndims);i++)
    {
      // j is the mask
      for (j=0;j<(1<<ndims);j++)
	{
	  COUNTBITS(j,nbits,tmp);
	  //we are considering a (ndims-nbits)-face
	  //remove this condition to consider 1-faces
	  //if (nbits!=ndims)
	    {
	      nfacelist = &(h->NFaceList[nbits]);
	      // n is our tag, each k-face has a unique tag
	      n = ((i&j)<<ndims) + j;
	      // this should be improved for speed
	      for (k=0;k<nfacelist->nfaces;k++)
		{
		  nface = (&nfacelist->NFace[k]);
		  if (nface->tag<0) 
		  {
		      nface->tag=n;
		  }
		  if (nface->tag==n) break;
		}
	      h->faceFromTag[n]=nface;
	      //if ((i&j)==0) 
	      nface->vertex[nface->nvertice++] = i;
	    }
	}
    }
  

  //change < to <= here to add the 1-faces
  mask = (1<<h->ndims)-1;
  for (n=0;n<=ndims-1;n++)
    {
      for (i=0;i<h->NFaceList[n].nfaces;i++)
	{
	  for (j=0;j<h->NFaceList[n+1].nfaces;j++)
	    {
	      nface = &(h->NFaceList[n].NFace[i]);
	      nface2 = &(h->NFaceList[n+1].NFace[j]);
	      INT MF=nface->tag & mask;
	      INT DF=(nface->tag & (mask<<h->ndims))>>h->ndims;

	      INT Mf=nface2->tag & mask;
	      INT Df=(nface2->tag & (mask<<h->ndims))>>h->ndims;

	      if (((Mf&MF)==MF)&&((DF&MF)==(Df&MF)))
		{
		    nface->Daughters[nface->ndaughters++] = nface2;
		    nface2->Mothers[nface2->nmothers++] = nface;
		}
	    }
	}
    }

  return h;
}

unsigned int HC_getCubeIndex(hypercube_type *h, char type, unsigned int id,int periodic)
{
    return id/h->NFaceList[h->ndims-type].nFacesPerCube;
}

unsigned int HC_getIndex(hypercube_type *h,const nface_type *face,const unsigned int *cube_coords, const int *dims, int periodic)
{
    nfacelist_type *flist = &(h->NFaceList[h->ndims-face->ndims]);
    int mask = (1<<h->ndims)-1;
    int M=face->tag & mask;
    int D=(face->tag & (mask<<h->ndims))>>h->ndims;

    long id=0;
    // now find the index of the NON-DISPLACED equivalent faces
    // this is slow and should be (easily) improved
    id = (long) (h->faceFromTag[M] - flist->NFace);
    //while (flist->NFace[id].tag != M) {id++;}
    
    unsigned int result=0;
    long dx=1;
    int i;
    for (i=0;i<h->ndims;i++)
    {
	if (D & (1<<i)) 
	{
	    if (periodic&(1<<i))
	    {
		if (cube_coords[i]+1>=dims[i])
		    result -= dx*(dims[i]-1);
		else
		    result += dx;
	    }
	    else
	    {
		printf("ERROR: non periodic boundary NOT implemented\n"); // find something to return ...
		exit(0);
		
	    }
	}
	result += cube_coords[i]*dx;
	dx*=dims[i];
    }

    return result*flist->nFacesPerCube + id;
}

nface_type *HC_getFace(hypercube_type *h, char type, unsigned int id, const int *dims,int periodic)
{
    nfacelist_type *flist = &(h->NFaceList[h->ndims-type]);
    int cube_index = id/flist->nFacesPerCube;

    if (periodic == ~0)
	return &(flist->NFace[id - cube_index*flist->nFacesPerCube]);
    else
    {
	return &(flist->NFace[id - cube_index*flist->nFacesPerCube]);
    }
}

void HC_getSymetricConfig(hypercube_type *h,const nface_type *face,const unsigned int *cube_coords,nface_type **new_face,unsigned int *new_cube_index,const int *dims,int periodic)
{
    int mask = (1<<h->ndims)-1;
    int M=face->tag & mask;
    int D=(face->tag & (mask<<h->ndims))>>h->ndims;
    long dx=1;
    unsigned int result=0;
    int i;
    
    for (i=0;i<h->ndims;i++)
    {
	if (M&(1<<i))
	{
	    if (D&(1<<i))
	    {
		// move in the positive direction
		if (periodic&(1<<i))
		{
		    if (cube_coords[i]+1>=dims[i])
			result -= dx*(dims[i]-1);
		    else
			result += dx;
		}
		else
		{
		    printf("ERROR: non periodic boundary NOT implemented\n"); // find something to return ...
		    exit(0);
		}
	    }
	    else
	    {
		// move in the negative direction
		if (periodic&(1<<i))
		{
		    if (cube_coords[i]==0)
			result += dx*(dims[i]-1);
		    else
			result -= dx;
		}
		else
		{
		    printf("ERROR: non periodic boundary NOT implemented\n"); // find something to return ...
		    exit(0);
		}
	    }
	    
	}
	result += cube_coords[i]*dx;
	dx*=dims[i];
    }

    *new_cube_index = result;
    *new_face = h->faceFromTag[(((~D)&M)<<h->ndims) + M];
}

void FreeHypercube(hypercube_type *h)
{
  int i,j;
  nface_type *nface;
  
  for (i=0;i<h->ndims-1;i++)
  {
      for (j=0;j<h->NFaceList[i].nfaces;j++)
      {
	  nface = &(h->NFaceList[i].NFace[j]);
	  free(nface->Mothers);
	  free(nface->Daughters);
	  free(nface->vertex);
      }
      free(h->NFaceList[i].NFace);
  }
  free(h->NFaceList);
  free(h);
  
}
