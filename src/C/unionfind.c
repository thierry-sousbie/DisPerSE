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
#include <time.h>
#include <assert.h>

#include "unionfind.h"
 
UF_type *uf_create(int N,int preserveID)
{
    UF_type *uf;
    long i;
    
    uf = calloc(1,sizeof(UF_type));
    uf->NNodes = N;
    uf->parent=calloc(N+1,sizeof(int));

    if (!preserveID) {
      uf->size=calloc(N+1,sizeof(int));
      for (i=0;i<N+1;i++) uf->size[i]=1;
    }
    else uf->size=NULL;
    
    uf->preserveID = preserveID;

    return uf;
}

int uf_free(UF_type **uf)
{
    
  free((*uf)->parent);
  (*uf)->parent=NULL;
  free((*uf)->size);
  (*uf)->size=NULL;
  free(*uf);
  *uf=NULL;
  return 0;
}

int uf_isOwnGroup(UF_type *uf, long node)
{
  return (uf->parent[node] == 0);
}

int uf_find(UF_type *uf, int node)
{
  if (uf->parent[node] == 0)
    return node;
 
  uf->parent[node] = uf_find(uf,uf->parent[node]-1)+1;
  return uf->parent[node]-1;
}
 
void uf_union(UF_type *uf, int node1, int node2)
{  
  int set1 = uf_find(uf,node1), set2 = uf_find(uf,node2);
  if (!uf->preserveID) 
    {
      if (set1 != set2) 
	{
	  int parentSet = (uf->size[set1] >= uf->size[set2]) ? set1 : set2;
	  int sonSet = (uf->size[set1] < uf->size[set2]) ? set1 : set2;
	  uf->parent[sonSet] = parentSet+1;
	  uf->size[parentSet] += uf->size[sonSet];
	}
    } 
  else 
    {
      if (set1 != set2) 
	{
	  uf->parent[set2] = set1+1;
	  //uf->size[set1] += uf->size[set2];
	}
    }
}

void uf_gr_union(UF_type *uf,int set1, int set2)
{
  //int set1=s1+1,set2=s2+1;
  if (!uf->preserveID) 
    {
      if (set1 != set2) 
	{
	  int parentSet = (uf->size[set1] >= uf->size[set2]) ? set1 : set2;
	  int sonSet = (uf->size[set1] < uf->size[set2]) ? set1 : set2;     
	  uf->parent[sonSet] = parentSet+1;
	  uf->size[parentSet] += uf->size[sonSet];
	}
    } 
  else 
    {
      if (set1 != set2) 
	{
	  uf->parent[set2] = set1+1;
	  
	  //assert(set2!=3798);
	  //uf->size[set1] += uf->size[set2];
	}
    }
}

/*
void uf_gr_union_preserve(UF_type *uf,int s1, int s2)
{    
  int set1=s1+1,set2=s2+1;
  if (set1 != set2) {
    uf->parent[set2] = set1;
    uf->size[set1] += uf->size[set2];
  }    
}
*/

UFL_type *ufl_create(long N,int preserveID)
{
    UFL_type *uf;
    long i;
    
    uf = calloc(1,sizeof(UFL_type));
    uf->NNodes = N;
    uf->parent=calloc(N+1,sizeof(long));

    if (!preserveID) {
      uf->size=calloc(N+1,sizeof(long));
      for (i=0;i<N+1;i++) uf->size[i]=1;
    } 
    else uf->size=NULL;
    
    uf->preserveID = preserveID;

    return uf;
}

int ufl_free(UFL_type *uf)
{
    
    free(uf->parent);
    uf->parent=NULL;
    free(uf->size);
    uf->size=NULL;
	    
    return 0;
}

int ufl_isOwnGroup(UFL_type *uf, long node)
{
  return (uf->parent[node] == 0);
}

long ufl_find(UFL_type *uf, long node)
{
    if (uf->parent[node] == 0)
        return node;
 
    uf->parent[node] = ufl_find(uf,uf->parent[node]-1)+1;
    return uf->parent[node]-1;
}

void ufl_union(UFL_type *uf, long node1, long node2)
{  
  long set1 = ufl_find(uf,node1), set2 = ufl_find(uf,node2);
  if (!uf->preserveID) 
    {
      if (set1 != set2) 
	{
	  long parentSet = (uf->size[set1] >= uf->size[set2]) ? set1 : set2;
	  long sonSet = (uf->size[set1] < uf->size[set2]) ? set1 : set2;
	  uf->parent[sonSet] = parentSet+1;
	  uf->size[parentSet] += uf->size[sonSet];
	}
    } 
  else 
    {
      if (set1 != set2) 
	{
	  uf->parent[set2] = set1+1;
	  //uf->size[set1] += uf->size[set2];
	}
    }
}


void ufl_gr_union(UFL_type *uf, long set1, long set2)
{
  //long set1=s1+1,set2=s2+1;
  if (!uf->preserveID) 
    {
      if (set1 != set2) 
	{
	  long parentSet = (uf->size[set1] >= uf->size[set2]) ? set1 : set2;
	  long sonSet = (uf->size[set1] < uf->size[set2]) ? set1 : set2;
	  uf->parent[sonSet] = parentSet+1;
	  uf->size[parentSet] += uf->size[sonSet];
	}
    } 
  else 
    {
      if (set1 != set2) 
	{
	  uf->parent[set2] = set1+1;
	  //uf->size[set1] += uf->size[set2];
	}
    }
}

/*
void ufl_gr_union_preserve(UFL_type *uf, long s1, long s2)
{
  long set1=s1+1,set2=s2+1;
  if (set1 != set2) {
      uf->parent[set2] = set1;
      uf->size[set1] += uf->size[set2];
    }    
}
*/
