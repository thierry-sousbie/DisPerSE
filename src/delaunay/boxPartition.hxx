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
#ifndef __BOX_PARTITION_HXX__
#define __BOX_PARTITION_HXX__

#include <stdio.h>
#include <math.h>
#include <vector>

class boxPartition
{
private:
  int ndims;
  long npart_t;
  long npart;
  std::vector<double> x0;
  std::vector<double> delta;
  std::vector<long> div;

  //std::vector< std::vector<double> > px0;
  //std::vector< std::vector<double> > pdelta;
  
  bool index2coord(int index, std::vector<long> &v) const
  {
    long i;
    long cur=index;
    long cur_fac=1;

    if (index<0) return false;
    if (index>=npart) return false;

    for (i=0;i<ndims;i++)
      {
	v[i]=(long)(cur/cur_fac) % div[i];
	cur-=cur_fac*v[i];
	cur_fac*=div[i];
      }

    if (cur!=0) return false;
    
    return true;
      
  }

  void computePartitions()
  {    
    double mindelta;
    double minid;
    long i;
   
    div.assign(ndims,1);
 
    do {
      minid=0;
      mindelta=delta[0]/div[0];
      for (i=1;i<ndims;i++) 
	{
	  if (delta[i]/div[i] > mindelta)
	    {
	      mindelta=delta[i]/div[i];
	      minid=i;
	    }
	}
      div[minid]++;

      for (i=0,npart=1;i<ndims;i++) npart*=div[i];
    } while (npart<npart_t);

  
  }

public:

  boxPartition():ndims(-1)
  {
  }

  ~boxPartition()
  {
  }

  long build(long _npart, int _ndims)
  {
    return build(_npart,std::vector<double>(_ndims,0),std::vector<double>(_ndims,1));
  }

  template <class T>
  long build(long _npart,const std::vector<T> &_x0,const std::vector<T> &_delta)
  {
     npart_t=_npart;

     assert(_x0.size()==_delta.size());
     assert(_x0.size()>0);
     
     x0.assign(_x0.begin(),_x0.end());
     delta.assign(_delta.begin(),_delta.end());  
     ndims=x0.size();

     computePartitions();

     return npart;
  }
  
  template <class OutputIterator>
  bool getX0(int index, OutputIterator out) const
  {
    if ((index<0)||(index>=npart)) return false;
    std::vector<long> v(ndims);
    if (!index2coord(index,v)) return false;
    for (long i=0;i<ndims;i++) *out=x0[i]+((double)v[i]/(double)div[i])*delta[i];
    return true;
  }
  
  template <class OutputIterator>
  bool getDelta(int index, OutputIterator out) const
  {
    if ((index<0)||(index>=npart)) return false;
    for (long i=0;i<ndims;i++) (*out)=delta[i]/div[i];
    return true;
  }

  long getNPart() const
  {
    return npart;
  }

};

#endif
