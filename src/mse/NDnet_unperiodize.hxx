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
#ifndef __NDNET_UNPERIODIZE_HXX__
#define __NDNET_UNPERIODIZE_HXX__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <vector>
#include <set>

#include "NDnetwork.h"

class NetworkUnperiodizer
{
public:
  NetworkUnperiodizer(NDnetwork *network_)
  {
    net=network_;
    NDIM=net->ndims;
    x0.assign(net->x0,net->x0+NDIM);
    delta.assign(net->delta,net->delta+NDIM);
    halfDelta=delta;
    xmax=x0;
    for (int i=0;i<NDIM;++i)
      {
	halfDelta[i]/=2;
	xmax[i] += delta[i];
      }
  } 

  NDnetwork *get()
  {
    std::vector<NewVertex> newVertices;
    std::set<NewVertex> nvSet;
    typedef std::set<NewVertex>::iterator Iterator;

    //long nVert = net->ndims_net+1;
    long type  = net->ndims_net;

    for (type=1;type<=net->ndims_net;++type)
      {
	if (!net->haveVertexFromFace[type]) continue;
	for (unsigned long i=0;i<net->nfaces[type];++i)
	  {
	    NDNET_UINT *vertexId = VERTEX_IN_FACE(net,type,i);
	    float refCoord[NDIM];
	    std::copy(&xmax[0],&xmax[0]+NDIM,refCoord);
	    for (int j=0;j<(type+1);++j)
	      {
		float *coord = &net->v_coord[vertexId[j]*NDIM];
		for (int k=0;k<NDIM;++k)
		  if (coord[k]<refCoord[k]) refCoord[k]=coord[k];
	      }
	
	    //float *refCoord = &net->v_coord[vertexId[0]*NDIM];
	
	    for (int j=0;j<(type+1);++j)
	      {
		float *coord = &net->v_coord[vertexId[j]*NDIM];
		NewVertex nv;
		int needCopy=0;
		for (int k=0;k<NDIM;++k)
		  {
		    int changed=0;
		    nv.newCoord[k]=checkCoordConsistency(coord[k],refCoord[k],k,changed);
		    if (changed) needCopy=1;
		  }
		if (needCopy)
		  {		
		    nv.index = vertexId[j];
		    nv.newIndex = newVertices.size();
		    std::pair<Iterator,bool> result=nvSet.insert(nv);
		    if (result.second)
		      {
			newVertices.push_back(nv);
			vertexId[j]=net->nvertex + nv.newIndex;
		      }
		    else 
		      {
			vertexId[j]=net->nvertex + (*result.first).newIndex;
		      }
		  }
	      }
	  }
      }
    
    net->v_coord = (float*)realloc(net->v_coord,sizeof(float)*(net->nvertex+newVertices.size())*NDIM);    
    for (unsigned long i=0;i<newVertices.size();++i)
      {
	std::copy(newVertices[i].newCoord,newVertices[i].newCoord+NDIM,
		  net->v_coord+(net->nvertex+i)*NDIM);
      }

    if (net->haveVFlags)
      {
	net->v_flag = 
	  (unsigned char*)realloc(net->v_flag,
				  sizeof(unsigned char)*(net->nvertex+newVertices.size()));
	for (unsigned long i=0;i<newVertices.size();++i)
	  net->v_flag[net->nvertex+i]=net->v_flag[newVertices[i].index];
      }

    for (int j=0;j<net->ndata;++j)
      {
	if (net->data[j].type != 0) continue;
	net->data[j].data=
	  (double*)realloc(net->data[j].data,
			   sizeof(double)*(net->nvertex+newVertices.size()));
	double *data=net->data[j].data;
	for (unsigned long i=0;i<newVertices.size();++i)
	  data[net->nvertex+i]=data[newVertices[i].index];	  
      }

    net->nvertex+=newVertices.size();

    return net;
  }
  
private:
  struct NewVertex
  {
    float newCoord[3];
    NDNET_UINT index;
    unsigned long newIndex;

    bool operator<(NewVertex const& rhs) const
    {
      //NewVertex const& lhs = *this;
      if (newCoord[0]<rhs.newCoord[0]) return true;
      else if (newCoord[0]==rhs.newCoord[0])
	{
	  if (newCoord[1]<rhs.newCoord[1]) return true;
	  else if (newCoord[1]==rhs.newCoord[1])
	    {
	      if (newCoord[2]<rhs.newCoord[2]) return true;
	    }
	}
      return false;
    }
  };
  
  template <class T>
  T correctCoordsDiff(T len, int which, int &changed) const
  {
    T result=len;
    changed=0;
    if (which<NDIM)
      {
	if (len>=halfDelta[which]) 
	  {
	    result= (len-delta[which]);
	    changed=1;
	  }
	else if (len<-halfDelta[which])
	  {
	    result= delta[which]+len;
	    changed=1;
	  }
      }
    return result;    
  }

  template <class T, class T2>
  T checkCoordConsistency(T vCoord, const T2 refPoint, int dim, int &changed) const
  {         
    double dif=correctCoordsDiff<T>(T2(vCoord)-refPoint,dim,changed);
    if (changed) return refPoint+dif;
    return vCoord;
    //return refPoint+correctCoordsDiff<T>(T2(vCoord)-refPoint,dim,changed);  
  }


  NDnetwork *net;
  std::vector<double> x0;
  std::vector<double> xmax;
  std::vector<double> delta;
  std::vector<double> halfDelta;
  int NDIM;
};

#endif
