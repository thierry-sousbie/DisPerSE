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
#include <vector>
#include <set>
#include <map>
#include <list>
#include <cmath>
#include <functional>
#include <stack>
#include <iterator>
#include <algorithm>

#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "arcs_nodes.hxx"
#include "Z2set.hxx"

namespace persistenceCycles
{
  typedef NDcomplex_node node;
  typedef node::node_it node_it;
  typedef NDcomplex_arc arc;
  typedef node::arc_it arc_it;

  typedef Z2set<node_it,node::compareItLess> cycle;
  typedef cycle::iterator cycle_it;
  typedef std::map<node_it,cycle,node::compareItLess> cycle_map;
  typedef cycle_map::value_type cycle_map_value;
  typedef cycle_map::iterator cycle_map_it;
  
  typedef std::map<node_it,std::list<node_it>,node::compareItLess> list_map;
  typedef list_map::value_type list_map_value;
  typedef list_map::iterator list_map_it;

  cycle_map_it getCycle(node_it nd, arc_it arcs_end, cycle_map &cycleMap)
  {
    cycle_map_it cmit=cycleMap.find(nd);
    if (cmit!=cycleMap.end()) return cmit;

    cycle c;    
    node_it refNode=nd->getPNode();

    if (nd->getType() == refNode->getType()) return cycleMap.end();
    
    std::vector<node_it> nei;
    std::vector<node_it>::iterator nei_it;
    refNode->getNeighbors(refNode,arcs_end,back_inserter(nei),nd->getType());
    c.insert(nei.begin(),nei.end());

    cycle_it cit=c.begin();
    std::set<node_it,node::compareItLess> inside;
    inside.insert(refNode);

    while (cit!=c.end())
      {
	node_it cur_refNode=(*cit)->getPNode();
	
	if (cur_refNode->getType() != refNode->getType())
	  {
	    cit++;
	    continue;
	  }

	if (*cit == nd) 
	  {
	    cit++;
	    continue;
	  }

	if (!inside.insert(cur_refNode).second)
	  {
	    cit++;
	    continue;
	  }
	
	node_it cur=*cit;
	cycle_map_it mit=cycleMap.find(cur);	
	if (mit != cycleMap.end())
	  {
	    c.insert(mit->second.begin(),mit->second.end());
	    cit=c.begin();
	    continue;
	  }
	
	nei.clear();	
	cur_refNode->getNeighbors(cur_refNode,arcs_end,back_inserter(nei),nd->getType());
	/*
	printf("new: %s\n",cur_refNode->getInfo(true).c_str());
	for (nei_it = nei.begin();nei_it!=nei.end();nei_it++)
	  {
	    printf("    ++ %s\n",(*nei_it)->getInfo(true).c_str());
	    printf("     >> %s\n",(*nei_it)->getPNode()->getInfo(true).c_str());
	  }
	*/
	c.insert(nei.begin(),nei.end());
	cit=c.begin();
      };
    
    cmit=cycleMap.insert(std::make_pair(nd,c)).first;
    
    return cmit;
  }

}
