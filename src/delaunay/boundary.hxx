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
#ifndef __BOUNDARY_CONDITIONS_H__
#define __BOUNDARY_CONDITIONS_H__

#include <vector>

class bConditions {
public: 
  enum BoundaryType {BType_Mirror, BType_Periodic, BType_None};

 private:
  std::vector<double> x0;
  std::vector<double> delta;
  BoundaryType bType;
  int ndims;

public: 

  std::vector<double> getX0() const
  {
    return x0;
  }

  std::vector<double> getDelta() const
  {
    return delta;
  }
 
  template <class T>
  bConditions(BoundaryType bTypep, std::vector<T> &x0p, std::vector<T> &deltap)
  {
    long i;
    x0.clear();delta.clear();
    for (i=0;i<=x0p.size();i++) x0.push_back(x0p[i]);
    for (i=0;i<=deltap.size();i++) delta.push_back(deltap[i]);
    ndims=x0p.size();
    bType = bTypep;
  }

  template <class T>
  bConditions(BoundaryType bTypep,T *x0p, T *deltap, int N)
  {
    long i;
    x0.clear();delta.clear();
    for (i=0;i<N;i++) x0.push_back(x0p[i]);
    for (i=0;i<N;i++) delta.push_back(deltap[i]);
    ndims=N;
    bType = bTypep;
  }

  ~bConditions()
  {
    
  }
  
  template <class T, class T2, class T3>
  static bool toSubBox_P(std::vector<T> &x, std::vector<T> &out, std::vector<T2> &x0s, std::vector<T2> &deltas, std::vector<T3> &x0p, std::vector<T3> &deltap)
  {
    int flag=0;
    long j,k,l;
    bool retval=false;
    //printf("P\n");
    for (j=0;j<x.size();j++)
      if ((x[j] < x0s[j])||(x[j] >= x0s[j]+deltas[j])) 
	  break;
      
    if (j==x.size()) {
      out.assign(x.begin(),x.end());
      retval=true;
    } else out.clear();

    flag=0;
    for (j=0;j<x.size();j++)
      {
	if (x0p[j]>x0s[j]) 
	  if (x[j] >= x0p[j]+deltap[j] - (x0p[j]-x0s[j])) flag|=(1<<(2*j));

	if (x0p[j]+deltap[j]<x0s[j]+deltas[j]) 
	  if (x[j] <= x0p[j]+ (x0s[j]+deltas[j]) - (x0p[j]+deltap[j])) flag|=(1<<(2*j+1));
      }
    
    if (flag) {
      int ncp=0;
      for (j=1;j<=flag;j++)
	{
	  int res=j&flag;
	  int res2=j&(~flag);
	  if ((res)&&(!res2)) {
	    T new_x[x.size()];
	    ncp++;
	    
	    for (k=0;k<x.size();k++)
	      {
		if (res&(1<<(2*k))) new_x[k]= x[k]-deltap[k];
		else if (res&(1<<(2*k+1))) new_x[k]= x[k]+deltap[k];
		else new_x[k]=x[k];
		//new_c.push_back(new_x[k]);
	      }
	    for (k=0;k<x.size();k++)
	      if ((new_x[k] < x0s[k])||(new_x[k] > x0s[k]+deltas[k])) 
		break;
	    if (k==x.size()) 
	      for (l=0;l<x.size();l++) out.push_back(new_x[l]);
	    
	  }
	}
    }
    return retval;
  }


  template <class T, class T2, class T3>
  static bool toSubBox_M(std::vector<T> &x, std::vector<T> &out, std::vector<T2> &x0s, std::vector<T2> &deltas, std::vector<T3> &x0p, std::vector<T3> &deltap)
  {
    int flag=0;
    long j,k,l;
    bool retval=false;
    //printf("M\n");
    for (j=0;j<x.size();j++)
      if ((x[j] < x0s[j])||(x[j] > x0s[j]+deltas[j])) 
	  break;
      
    if (j==x.size()) {
      out.assign(x.begin(),x.end());
      retval=true;
    } else out.clear();

    flag=0;
    for (j=0;j<x.size();j++)
      {
	if ((x0p[j]>x0s[j])&&(x[j]!=x0p[j])) 
	  if (x[j] <= x0p[j]+ (x0p[j]-x0s[j])) flag|=(1<<(2*j));

	if ((x0p[j]+deltap[j]<x0s[j]+deltas[j])&&(x[j]!=deltap[j]+x0p[j])) 
	  if (x[j] >= x0p[j]+deltap[j] - ((x0s[j]+deltas[j])-(x0p[j]+deltap[j]))) flag|=(1<<(2*j+1));
	  
      }

    if (flag) {
      int ncp=0;
      for (j=1;j<=flag;j++)
	{
	  int res=j&flag;
	  int res2=j&(~flag);
	  if ((res)&&(!res2)) {
	    T new_x[x.size()];
	    ncp++;
	    
	    for (k=0;k<x.size();k++)
	      {
		if (res&(1<<(2*k))) new_x[k]= x[k]-2*(x[k]-x0p[k]);
		else if (res&(1<<(2*k+1))) new_x[k]= x[k]+2*(x0p[k]+deltap[k]-x[k]);
		else new_x[k]=x[k];
		//new_c.push_back(new_x[k]);
	      }
	    
	    for (k=0;k<x.size();k++)
	      if ((new_x[k] < x0s[k])||(new_x[k] > x0s[k]+deltas[k])) 
		break;
	    if (k==x.size()) 
	      for (l=0;l<x.size();l++) out.push_back(new_x[l]);
	    
	  }
	}
    }
    return retval;
  }

  
  template <class T, class T2, class T3>
  static bool toSubBox_N(std::vector<T> &x, std::vector<T> &out, std::vector<T2> &x0s, std::vector<T2> &deltas, std::vector<T3> &x0p, std::vector<T3> &deltap)
  {
    int flag=0;
    long j,k,l;
    bool retval=false;

    for (j=0;j<x.size();j++)
      if ((x[j] < x0s[j])||(x[j] > x0s[j]+deltas[j])) 
	  break;
      
    if (j==x.size()) {
      out.assign(x.begin(),x.end());
      retval=true;
    } else out.clear();

    return retval;
  }
  
  template <class T, class T2>
  bool inBBox(const std::vector<T> &x, const std::vector<T2> &x0s, const std::vector<T2> &deltas)
  {
    long j;
    
    for (j=0;j<x.size();j++)
      if ((x[j] < x0s[j])||(x[j] >= x0s[j]+deltas[j])) 
	break;

    return (j==x.size());
  }
  
  template <class T, class T2>
  bool inBBox(const T *x, const std::vector<T2> &x0s, const std::vector<T2> &deltas)
  {
    long j;
    
    for (j=0;j<x0s.size();j++)
      if ((x[j] < x0s[j])||(x[j] >= x0s[j]+deltas[j])) 
	break;

    return (j==x0s.size());
  }

  template <class T>
  bool inBBox(const T *x)
  {
    long j;
    
    for (j=0;j<x0.size();j++)
      if ((x[j] < x0[j])||(x[j] >= x0[j]+delta[j])) 
	break;

    return (j==x0.size());
  }

  template <class T>
  bool inBBox(const std::vector<T> &x)
  {
    long j;
    
    for (j=0;j<x.size();j++)
      if ((x[j] < x0[j])||(x[j] >= x0[j]+delta[j])) 
	break;

    return (j==x.size());
  }

  template <class T, class T2>
  bool toSubBox(std::vector<T> &x, std::vector<T> &out, std::vector<T2> &x0s, std::vector<T2> &deltas)
  {
    switch (bType)
      {
      case BType_Periodic: 
	return toSubBox_P(x,out,x0s,deltas,x0,delta);
	break;
	
      case BType_Mirror:
	return toSubBox_M(x,out,x0s,deltas,x0,delta);
	break;
	
      case BType_None:
	return toSubBox_N(x,out,x0s,deltas,x0,delta);
	break;
      }
  return false;
  }  

  template <class T, class T2>
  bool toSubBox(const T *px, std::vector<T> &out, std::vector<T2> &x0s, std::vector<T2> &deltas)
  {
    std::vector<T> x(px,px+ndims);
    switch (bType)
      {
      case BType_Periodic: 
	return toSubBox_P(x,out,x0s,deltas,x0,delta);
	break;
	
      case BType_Mirror:
	return toSubBox_M(x,out,x0s,deltas,x0,delta);
	break;
	
      case BType_None:
	return toSubBox_N(x,out,x0s,deltas,x0,delta);
	break;
      }
  }

  

  
  
};

#endif
