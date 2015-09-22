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
/* This is NOT finished */

#ifndef __NDSKEL_TRIM_HXX__
#define __NDSKEL_TRIM_HXX__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include "float.h"
#include <list>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <limits>

#include <assert.h>

#include "NDsubSkeleton.hxx"

#include "mytypes.h"
#include "NDskeleton.h"
#include "global.h"

class NDskelTrimT
{
private:

  template <class compareLess = std::less<double> >
  class sklSelectThresholdT
  {
  private:
    double thresh;
    int valID;
    int valID_p1;
    int valID_p2;
    
  public:
    
    sklSelectThresholdT(NDskel *skl, double threshold, std::string valueField = std::string(VALUE_TAG))
    {
      valID=NDskel_NodeDataIndex(skl,valueField.c_str());
      valID_p1=NDskel_SegDataIndex(skl,(valueField+std::string("_p1")).c_str());
      valID_p2=NDskel_SegDataIndex(skl,(valueField+std::string("_p2")).c_str());  
      thresh=threshold;

      if ((valID_p1<0)||(valID_p2<0))
      {
	valID_p1=NDskel_SegDataIndex(skl,valueField.c_str());
	valID_p2=valID_p1;
      }
      
      if ((valID<0)||(valID_p1<0)||(valID_p2<0))
	{
	  fprintf(stderr,"ERROR in NDskelAssembleT, data field '%s' is absent for nodes and/or segments.\n",valueField.c_str());
	  exit(0);
	}
    }
  
    bool test(const NDskl_node& node) const 
    {
      return compareLess()(node.data[valID],thresh);
    }    

    std::pair<bool,bool> test(const NDskl_seg& seg) const 
    {
      bool res1=compareLess()(seg.data[valID_p1],thresh);
      bool res2=compareLess()(seg.data[valID_p2],thresh);
      return std::make_pair(res1,res2);
    }    

  };

  class sklSelectBoundaryT
  {
  private:
    bool rmAll;
    NDskel *skl;

  public:
    
    sklSelectBoundaryT(NDskel *skl_, bool all=true)
    {
      rmAll=all;
      skl=skl_;
    }
  
    bool test(const NDskl_node& node) const 
    {
      if (rmAll) 
	return (node.flags&(FLAG_NDNODE_OUT|FLAG_NDNODE_BOUNDARY|FLAG_NDNODE_INFINITE));
      else
	return (node.flags&(FLAG_NDNODE_OUT|FLAG_NDNODE_INFINITE));
    }    

    std::pair<bool,bool> test(const NDskl_seg& seg) const 
    {
      std::pair<bool,bool> res;
      res.first=test(skl->Node[seg.nodes[0]])&&test(skl->Node[seg.nodes[1]]);
      res.second=res.first;
      return res;
    }    

  };

private:

  NDskel *skl;

public:

  template <class CondT> 
  NDskel *trim(CondT cond)
  {
    
    typedef std::list<NDskl_seg *> branchesT;
    typedef branchesT::iterator branches_it;
    long total=0;
    long i,j;
 
    branchesT branch;
    //bool AllowFilCrossingTh=true;
    branch.clear();
    printf ("Trimming skeleton ...");fflush(0);
    
    for (i=0;i<skl->nnodes;i++)
      {
	for (j=0;j<skl->Node[i].nnext;j++)
	  {
	    if (skl->Node[i].Next[j]->index<skl->Node[i].index) continue;
	    branch.push_back(skl->Node[i].Seg[j]);
	  }
      }

    skeletonBuildT newSkl;
    //std::list<NDskl_seg *> fil;
    //std::vector<int> fil_dir;
    
    for (branches_it it=branch.begin();it!=branch.end();it++)
      {
	skeletonBuildT::newNodeIt na;
	skeletonBuildT::newNodeIt nb;
	std::list<NDskl_seg *>::iterator fil_beg;
	std::list<NDskl_seg *>::iterator fil_it;

	std::list<NDskl_seg *> fil;
	int segdir = skeletonBuildT::getFilament(*it,std::back_inserter(fil));
	if (!fil.size()) continue;
	NDskl_node *node1,*node2;

	if (segdir>0) 
	  {
	    node1=&skl->Node[(*it)->nodes[0]];
	    node2=&skl->Node[(*it)->nodes[1]];
	  }
	else
	  {
	    node1=&skl->Node[(*it)->nodes[1]];
	    node2=&skl->Node[(*it)->nodes[0]];
	  }
	
	bool isAbove=!cond.test(*node1);

	if (isAbove) 
	  {
	    
	    na=newSkl.addNode(node1,skl->ndims);
	    fil_beg= fil.begin();
	  }
	
	i=0;
	for (fil_it = fil.begin() ;fil_it!=fil.end();fil_it++,i++)
	  {
	    bool crossing;
	    std::pair<bool,bool> res=cond.test(**fil_it);

	    if (isAbove)
	      crossing=res.first||res.second;
	    else
	      crossing=(!res.first)||(!res.second);
	      
	    if (crossing)
	      {	
		if (isAbove)
		  {
		    if (segdir>0) nb=newSkl.addNode(*fil_it,skl->ndims,1);
		    else nb=newSkl.addNode(*fil_it,skl->ndims,0);
		      
		    if (na!=nb) 
		      {
			newSkl.addFilament(na,nb);
			newSkl.AddSegToFilament(fil_beg,fil_it);
			newSkl.AddSegToFilament(*fil_it);
		      }
		  }
		else
		  {
		    if (segdir>0) na=newSkl.addNode(*fil_it,skl->ndims,0);
		    else na=newSkl.addNode(*fil_it,skl->ndims,1);

		    fil_beg=fil_it;
		  }
		isAbove = !isAbove;
	      }
	  }
	
	if (!cond.test(*node2))
	  nb=newSkl.addNode(node2,skl->ndims);
	
	if (isAbove)
	  {
	    if (cond.test(*node2))
	      {
		if (segdir>0)
		  nb=newSkl.addNode(fil.back(),skl->ndims,1);
		else
		  nb=newSkl.addNode(fil.back(),skl->ndims,0);
	      }	    
	    
	    if (na!=nb) 
	      {
		newSkl.addFilament(na,nb);
		newSkl.AddSegToFilament(fil_beg,fil.end());
	      }
	    	      
	  }
	
	
      }
    //printf("Total = %ld\n",total);
    NDskel* ret=newSkl.toNDskel(skl);	
    printf("done.\n");
    return ret;
  }

  NDskel *trimBelow(double threshold, std::string valueField = std::string(VALUE_TAG))
  {
    sklSelectThresholdT< std::less<double> > cond(skl,threshold,valueField);
    return trim(cond);
  }

  NDskel *trimAbove(double threshold, std::string valueField = std::string(VALUE_TAG))
  {
    sklSelectThresholdT< std::greater<double> > cond(skl,threshold,valueField);
    return trim(cond);
  }

  NDskel *trimBoundary(bool all=true)
  {
    sklSelectBoundaryT cond(skl,all);
    return trim(cond);
  }
  
  NDskelTrimT(NDskel *skl_)
  {
    skl=skl_;
  }

};

#endif
