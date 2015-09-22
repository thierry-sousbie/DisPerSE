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
#ifndef __NDSKEL_BREAKDOWN_TOOLS_HXX__
#define __NDSKEL_BREAKDOWN_TOOLS_HXX__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <list>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include "NDsubSkeleton.hxx"
#include "NDskeleton.h"
#include "global.h"

class compareSegPosEqual //: public std::binary_function <NDskl_seg *,NDskl_seg *,bool>
{
private:
  int ndims;
  bool strict;

public:
  compareSegPosEqual(int ndims_, bool strict_=true)
  {
    ndims=ndims_;
    strict=strict_;
  }

  bool operator()(NDskl_seg const &a,NDskl_seg const &b)
  {
    int j=0;
    for (int i=0;i<ndims;i++)
      {
	if (a.pos[i]!=b.pos[i]) j|=(1<<0);
	if (a.pos[ndims+i]!=b.pos[ndims+i]) j|=(1<<1);
	if (a.pos[i]!=b.pos[ndims+i]) j|=(1<<2);
	if (a.pos[ndims+i]!=b.pos[i]) j|=(1<<3);
      }
   
    j=~j;   
    if (( (j&(1<<0)) && (j&(1<<1)) )||
	( (j&(1<<2)) && (j&(1<<3)) ) )
      {
	return true;
      }  

    return false;
  }

  bool operator()(NDskl_seg * const &a,NDskl_seg * const &b)
  {
    return operator()(*a,*b);
  }
};

  
class compareSegPosLess //: public std::binary_function <NDskl_seg *,NDskl_seg *,bool>
{
private:
  int ndims;
  compareSegPosEqual cmpEq;

public:
  compareSegPosLess(int ndims_, bool strict_=true):cmpEq(ndims_,strict_)
  {
    ndims=ndims_;
  }

  bool operator()(NDskl_seg* const &a,NDskl_seg* const &b)
  {
    // arbitrary but definite ordering of equal values ...
    if (cmpEq(a,b)) 
      return a<b;
      
    for (int i=0;i<ndims;i++)
      {
	if (a->pos[i]!=b->pos[i]) return a->pos[i]<b->pos[i];
	if (a->pos[ndims+i]!=b->pos[ndims+i]) return a->pos[ndims+i]<b->pos[ndims+i];
      }   
    return false;
  }
};


class compareNodePosEqual //: public std::binary_function <NDskl_node *,NDskl_node *,bool>
{
private:
  int ndims;
  bool strict;

public:
  compareNodePosEqual(int ndims_, bool strict_=true)
  {
    ndims=ndims_;
    strict=strict_;
  }

  bool operator()(NDskl_node_str const &a,NDskl_node_str const &b)
  {
    if (&a == &b) return true;

    for (int i=0;i<ndims;i++)
      {
	if (a.pos[i]!=b.pos[i]) return false;
      }

    return true;
  }

  bool operator()(NDskl_node_str* const &a,NDskl_node_str* const &b)
  {
    return operator()(*a,*b);
  }
};

class compareNodePosLess //: public std::binary_function <NDskl_node *,NDskl_node *,bool>
{
private:
  int ndims;
  compareNodePosEqual cmpEq;

public:
  compareNodePosLess(int ndims_, bool strict_=true):cmpEq(ndims_,strict_)
  {
    ndims=ndims_;
  }

  bool operator()(NDskl_node const &a,NDskl_node const &b)
  {
    for (int i=0;i<ndims;i++)
      {
	if (a.pos[i]!=b.pos[i]) return a.pos[i]<b.pos[i];
      }   

    return false;
  }

  bool operator()(NDskl_node* const &a,NDskl_node* const &b)
  {
    return operator()(*a,*b);
  }
};


struct compareSignatureEq
{
  bool operator()(std::vector<long> const &a,std::vector<long> const &b)
  {
    if (a.size()!=b.size()) return false;

    for (unsigned long i=0;i<a.size();i++)
      if (a[i]!=b[i]) return false;

    return true;
  }

  bool operator()(std::vector<long>* const &a,std::vector<long>* const &b)
  {
    return operator()(*a,*b);
  }
};

struct compareSignatureLess
{
  bool operator()(std::vector<long> const &a,std::vector<long> const &b)
  {
    if (a.size()!=b.size()) return (a.size()<b.size());

    for (unsigned long i=0;i<a.size();i++)
      if (a[i]!=b[i]) return a[i]<b[i];

    return false;
  }

  bool operator()(std::vector<long>* const &a,std::vector<long>* const &b)
  {
    return operator()(*a,*b);
  }
};

#endif
