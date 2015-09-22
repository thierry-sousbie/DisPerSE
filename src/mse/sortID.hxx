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
#ifndef __SORT_BY_ID_HXX__
#define __SORT_BY_ID_HXX__

#include <iostream>
#include <algorithm>
#include <vector>

namespace srtID {

  template<class RT, class T> struct compareID_less {
    compareID_less(const T arr,const int stride) : arr(arr),stride(stride) {}
    
    bool operator()(const RT a, const RT b) const
    { 
      if (arr[a]==arr[b])
	{
	  if (stride==1) return a<b;
	  //for (int i=1;i<stride;i++) if (arr[a+i]<arr[b+i]) return true;
	  return a<b;
	}
      return arr[a]<arr[b];
    }
    const T arr;
    const int stride;
  };

  template<class RT, class T> struct compareID_more {
    compareID_more(const T arr,const int stride) : arr(arr),stride(stride) {}

    bool operator()(const RT a, const RT b) const
    { 
      if (arr[a]==arr[b])
	{
	  if (stride==1) return a>b;
	  //for (int i=1;i<stride;i++) if (arr[a+i]>arr[b+i]) return true;
	  return a>b;
	}
      return arr[a]>arr[b];
    }
    const T arr;
    const int stride;
  };

  template<class RT, class T>
  RT *sortByID_ptr(T *arr, long size, int stride, int reverse)
  {
    RT* id=(RT*)malloc(sizeof(RT)*size);
    long i;
    for (i=0;i<size;i++) id[i]=i*stride;

    if (!reverse) std::sort(id, id+size, compareID_more<RT,T*>(arr,stride));
    else std::sort(id, id+size, compareID_less<RT,T*>(arr,stride));

    if (stride)
      {
	for (i=0;i<size;i++) id[i]/=stride;
      }
    
    return id;
  }
};

#endif
