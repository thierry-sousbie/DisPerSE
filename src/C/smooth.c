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

#include "smooth.h"

// does the job in MatrixSmoothSkeleton
// val is a list of vectors of size ndims
// val[1..ndims]: pos of first point in first segment
// val[ndims..2*ndims]: pos of second point in first segment
// val[2*ndims..3*ndims]: pos of first point in second segment,
// this should be the same as val[ndims..2*ndims] 
float *MatrixSmooth(float *val,int nval, int ndims)
{
  long i,d;
  float *newval=NULL;
  float tmp;



	
  newval=malloc (nval*sizeof(float)*ndims);
	
  memcpy(newval,val,nval*sizeof(float)*ndims);
	
  for (i=1;i<nval-2;i+=2)
    {
      for (d=0;d<ndims;d++)
	{
	  tmp = val[i*ndims+d]-val[(i+1)*ndims+d];
	  if (tmp==0)
	    {
	      newval[i*ndims+d] = newval[(i+1)*ndims+d] = 0.25*(val[i*ndims+d] + val[(i-1)*ndims+d] + val[(i+1)*ndims+d]+ val[(i+2)*ndims+d]);
	    }
	  else
	    {
	      newval[i*ndims+d]  = 0.25*(val[i*ndims+d] + val[(i-1)*ndims+d] + val[(i+1)*ndims+d]+ val[(i+2)*ndims+d]+2*tmp);
	      newval[(i+1)*ndims+d] = 0.25*(val[i*ndims+d] + val[(i-1)*ndims+d] + val[(i+1)*ndims+d]+ val[(i+2)*ndims+d]-2*tmp);
	    }
	}      
    }
	
  memcpy(val,newval,nval*sizeof(float)*ndims);
  free(newval);
  return val;
}

// Smooth a skeleton using a simple neighbour position averaging.
// Call n times to smooth over n segments ...
NDskel *MatrixSmoothSkeleton(NDskel *skl, int periodic)
{
  long i,j,n;
  float *pos=NULL;
  long nalloc=0;
  NDskl_seg *seg;
	
  for (i=0;i<skl->nnodes;i++)
    {
      for (j=0;j<skl->Node[i].nnext;j++)
	{
	  seg = skl->Node[i].Seg[j];
	  if (seg->Next!=NULL)
	    {
	      if (skl->Node[i].nsegs[j]>nalloc)
		{
		  nalloc = skl->Node[i].nsegs[j];
		  pos = realloc(pos,sizeof(float)*2*skl->ndims*nalloc);
		}
	      n=0;
	      do {
		memcpy(&(pos[n]),seg->pos,2*skl->ndims*sizeof(float));
		n+=2*skl->ndims;
		seg=seg->Next;
	      } while (seg!=NULL);
				
	      MatrixSmooth(pos,2*skl->Node[i].nsegs[j],skl->ndims);
	      seg = skl->Node[i].Seg[j];
	      n=0;
	      do {
		memcpy(seg->pos,&(pos[n]),2*skl->ndims*sizeof(float));
		n+=2*skl->ndims;
		seg=seg->Next;
	      } while (seg!=NULL);
	    }
	}
    }
	
  if (pos!=NULL) free(pos);
	
  return skl;
}

NDskel *smoothSkeleton(NDskel *skl, int N, int periodic)
{
	NDskel *myskl=skl;
	int i;

	for (i=0;i<N;i++)
		myskl = MatrixSmoothSkeleton(myskl,periodic);

	return myskl;
}

