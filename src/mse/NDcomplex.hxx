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
#ifndef __ND_COMPLEX_HEADER__
#define __ND_COMPLEX_HEADER__

#include <sys/time.h>

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
#include "unionfind.h"

#include "global.h"

#include "Z2set.hxx"
#include "persistenceCycles.hxx"
#include "boundaryMatrix.hxx"
#include "arcs_nodes.hxx"
#include "updtQueue.hxx"
#include "cells.hxx"
#include "ompPSort.hxx"

#include <fstream>
//#include "NDsubNetwork.hxx"

#ifdef USE_OPENMP
#include <parallel/algorithm>
#endif

#ifdef USE_THREADS
#include "pthreadBarrier.hxx"
#endif

class NDcomplex {
public:
    
  //static const double infinity = DBL_MAX;
  //arcs and nodes typedefs
  typedef NDcomplex_cellType cellType;

  typedef NDcomplex_arc arc;
  typedef NDcomplex_node node;
  typedef arc::arcs_list arcs_list;
  typedef node::nodes_list nodes_list;
  typedef arcs_list::iterator arcs_list_it;
  typedef nodes_list::iterator nodes_list_it;
  typedef arc::arc_it arc_it;
  typedef node::node_it node_it;
  typedef arc::arc_rit arc_rit;
  typedef node::node_rit node_rit;

  // arcs and nodes geometry typedefs
  
  typedef arc::geom_it arcGeom_it;
  typedef node::geom_it nodeGeom_it;
  
  typedef persistenceCycles::cycle cycle;
  typedef persistenceCycles::cycle_map cycle_map;
  typedef persistenceCycles::cycle_it cycle_it;
  typedef persistenceCycles::cycle_map_it cycle_map_it;
  
  
  enum cancellationMode {noBoundary=0, preserveBoundary=1};

private:
  cancellationMode currentMode;
  bool havePPairs;
  std::vector<bool> haveCycles;
  node::geom nodesGeom;
  arc::geom arcsGeom;
  
  arcs_list arcs;
  nodes_list nodes;

  arcs_list deleted_arcs;
  nodes_list deleted_nodes;

  std::vector<arcs_list> arcPool ;
  std::vector<nodes_list> nodePool ;

  int ndims;
  int ndims_net;
  
  void updatePersistencePair(const node_it n);
  //void computePersistencePairs_pairs();
  sparseZ2Matrix *buildBoundaryMatrix(std::vector<node_it> &nodeTab,char curType, bool reverse, bool skipAlreadyPaired=true);
  long setPairsFromBoundaryMatrix(const sparseZ2Matrix &M, const std::vector<node_it> &nodeTab, bool reverse);
  bool setCyclesFromBoundaryMatrix(int cycleType, sparseZ2Matrix &M, const std::vector<node_it> &nodeTab, bool reverse);
  void propagateReferenceFrom1Cycles();

  struct matrix_parms {
    int crit_index;
    bool up;
    bool pairs;
    bool cycles;
    matrix_parms (int a, bool b, bool c, bool d)
    {crit_index=a;up=b;pairs=c;cycles=d;}
  };

  void computePersistencePairs_matrix(const std::vector<int> compCycles=std::vector<int>());

  //void computePersistencePairs_incremental();    

  arc_it cancelPair(arc_it it, bool &canceled,   std::list<node_it> &p1, std::list<node_it> &p2, std::list<node_it> &p3, arc_it cur_it, bool force = false);

  //inline bool computeReferenceNodes_rec();

public:
  //node_it getInfiniteNode() {return infinite_node;}
  node_it getParentNode(node_it node, bool useRatio);
  std::pair<node_it,node_it> getReferenceNodes(arc_it arc);
  std::pair<node_it,node_it> getReferenceNodes_cycles(arc_it arc);
  //std::pair<node_it,node_it> getReferenceNodes_neighbors(arc_it arc);
  void setCancellationMode(cancellationMode mode) {currentMode=mode;}
  cancellationMode getCancellationMode() {return currentMode;}

  int getNDims(bool network=false) const {return (network)?ndims_net:ndims;}
  bool getPPairsState() {return havePPairs;}
  void resetPPairsState() {havePPairs=false;}
  std::vector<bool> getCyclesState() {std::vector<bool> res(haveCycles);res.resize(ndims+1);return res;}
  void resetCyclesState() {haveCycles.clear();}

  void init(int ndims_p, int ndims_net_p=-1, cancellationMode mode = noBoundary)
  {
    if (ndims!=-1) fprintf (stderr,"WARNING: in NDcomplex : multiple initializations.\n");		
    if (ndims_p<1) return;
    havePPairs=false;
    setCancellationMode(mode);
    ndims=ndims_p;
    if (ndims_net_p<0) ndims_net=ndims;
    else ndims_net=ndims_net_p;
    arcsGeom.init();		
    nodesGeom.init();
    //infinite_node = addNode(ndims,0,0,0,NODE_FLAG_INFINITE);
  }

    NDcomplex(int ndims_p=-1, int ndims_net_p=-1, cancellationMode mode = noBoundary)
	{
	    setCancellationMode(mode);
	    ndims=-1;
	    
	    init(ndims_p,ndims_net_p,mode);	    
	}

    ~NDcomplex()
	{}

  void clearNodeTags() {for (nodes_list_it it = nodes.begin();it!=nodes.end();it++) it->clearTag();}
  //void clearNodeTags2() {for (nodes_list_it it = nodes.begin();it!=nodes.end();it++) it->clearTag2();}

  int num_arcs() const {return arcs.size();}
  int num_nodes() const {return nodes.size();}

  arc_it arcs_begin() {return arcs.begin();}
  arc_it arcs_end() {return arcs.end();}
  arc_it rm_arcs_begin() {return deleted_arcs.begin();}
  arc_it rm_arcs_end() {return deleted_arcs.end();}

  node_it nodes_begin() {return nodes.begin();}
  node_it nodes_end() {return nodes.end();}
  node_it rm_nodes_begin() {return deleted_nodes.begin();}
  node_it rm_nodes_end() {return deleted_nodes.end();}

  arcGeom_it arcsGeom_end() {return arcsGeom.end();}
  nodeGeom_it nodesGeom_end() {return nodesGeom.end();}

  inline void sortNodes(bool (*cmp)(const node &,const node &)=node::comparePersistenceLess);
  
  inline arc_it addArc(node_it n1, node_it n2, int pool=0);
  inline arc_it removeArc(arc_it &it, bool storeDeleted=true);
  inline node_it addNode(cellType cell, double val,double val2=0, char flags=0, char type=-1, int pool=0);
  inline node_it removeNode(node_it node, bool storeDeleted=true);
  inline node_it removeNode(node_it node, arc_it &arc, bool storeDeleted=true);

  //inline bool computeReferenceNodes();
  inline void restoreRemovedArcsLinks(double robustness);
  inline void unrestoreRemovedArcsLinks();
  
  inline void createPools(int N);
  inline void commitPools();

  template <class InputIterator>
  inline arc::geom_it insertArcGeometry(InputIterator begin, InputIterator end, int pool=0);
  template <class InputIterator>
  inline node::geom_it insertNodeGeometry(InputIterator begin, InputIterator end, char type, int pool=0);  

  //inline arc::geom_it insertArcNullGeometry(bool store_geometry);
  inline arc::geom_it setArcNullGeometry(arc_it a, bool store_geometry);  
  inline node::geom_it insertNodeNullGeometry(char type, bool store_manifolds);  

  inline arc::geom_it duplicateArcGeometry(arcGeom_it ag);

  //bool defragArcsGeometry();
  //bool defragNodesGeometry();

  void toAscii(const char *fname="skl.dat");
  //void cyclesToAscii(const char *fname="cycles.dat");

  void simplify(std::vector<double> level, bool useRatio=false);
  int cancelTaggedNodes_Raw(bool force=false,bool solve=false);
  int cancelTaggedNodes_Smart(bool force=false, bool solve=false);

  template <class OutputIterator>
  int cancelTaggedNodes_Smart(OutputIterator imp_list, OutputIterator con_list, bool force, bool solve);
  
  int cancelTaggedNodes(bool force=false, bool solve=false) {return cancelTaggedNodes_Smart(force,solve);}
    

  inline arc_it cancelPair(arc_it it, bool &canceled, arc_it cur_it, bool force = false)
  {
    std::list<node_it> p1;
    std::list<node_it> p2;
    std::list<node_it> p3;

    return cancelPair(it,canceled,p1,p2,p3,cur_it,force);
  }
    
  inline arc_it cancelPair(arc_it it, bool &canceled, bool force = false) 
  {
    std::list<node_it> p1;
    std::list<node_it> p2;
    std::list<node_it> p3;

    return cancelPair(it,canceled,p1,p2,p3,it,force);
  }

  template<class OutputIterator>
  OutputIterator cancelPair(arc_it it, bool &canceled ,OutputIterator oit, bool force = false) 
  {
    std::list<node_it> p1;
    std::list<node_it> p2;
    std::list<node_it> p3;

    cancelPair(it,canceled,p1,p2,p3,it,force);
    for (std::list<node_it>::iterator lit=p1.begin();lit!=p1.end();lit++) *oit = *lit;
    for (std::list<node_it>::iterator lit=p2.begin();lit!=p2.end();lit++) *oit = *lit;
    for (std::list<node_it>::iterator lit=p3.begin();lit!=p3.end();lit++) *oit = *lit;

    return oit;
  }

  bool computePersistencePairs(const std::vector<int> compCycles=std::vector<int>());     
    
  void buildReflectiveBoundary();
  //bool sanitize(bool cancelNullPairs=true);
  void sanitizeBoundary(bool vertexAsMaxima, int store_manifolds, int store_arcsGeom);
  //void sanitizeBoundary_old();
  void removeBoundaries(bool removeAll=true);
  void removeOutNodes();
  bool haveBoundaries();

    class pPairExclude: std::binary_function <node_it,node_it,bool> {
    public :
	static char unMark(node_it node) {return node->flagsAnd(~NODE_FLAG_PAIRED);}
	static char mark(node_it node) {return node->flagsOr(NODE_FLAG_PAIRED);}
	static bool isMarked(node_it node) {return node->getFlags()&NODE_FLAG_PAIRED;}
	bool operator()(const node_it& n1, const node_it& n2)
	    {
		if ((n1->getFlags()&NODE_FLAG_PAIRED)||(n2->getFlags()&NODE_FLAG_PAIRED)) return true;
		if ((n1->getFlags()&(NODE_FLAG_POSITIVE|NODE_FLAG_NEGATIVE))==
		    (n2->getFlags()&(NODE_FLAG_POSITIVE|NODE_FLAG_NEGATIVE))) return true;
		
		return false;
	    }
    };

  void write(std::ofstream &str);
  void read(std::ifstream &str);

  std::string nodeType2Str(int type)
  {
    std::ostringstream stm;

    if (type==0) stm<<"minima";
    else if (type==ndims_net) stm<<"maxima";
    else stm<<type<<"-saddle";

    return stm.str();
  }

  std::string getInfo()
  {
    std::ostringstream stm;
    std::vector<long> node_count(ndims+1,0);
    std::vector<long> arc_count(ndims,0);
    int i;

    for (node_it n =nodes.begin(); n!=nodes.end(); n++)
      {
	node_count[n->getType()]++;
      }

     for (arc_it a =arcs.begin(); a!=arcs.end(); a++)
      {
	std::pair<node_it,node_it> p =a->endpoints();
	int id = std::min(p.first->getType(),p.second->getType());
	arc_count[id]++;
      }
     
     
     stm<<"    ";
     for (i=0;i<=ndims_net;i++)
       {
	 stm<<node_count[i]<<" ";
	 stm<<nodeType2Str(i);
	 
	 if (i!=ndims_net) stm<<", ";
       }
     stm<<"\n";

     for (i=0;i<ndims_net;i++)
       {
	 stm<<"    "<<arc_count[i]<<" ";
	 if (i==0) stm<<"(minimum / 1-saddle) arcs\n";
	 else if (i==ndims_net-1) stm<<"("<<i<<"-saddle / maximum) arcs\n";
	 else stm<<"("<<i<<"-saddle /"<<i+1<<"-saddle) arcs\n";
	 
       }
     return stm.str();
  }
    
};

inline void NDcomplex::
sortNodes(bool (*cmp)(const node &,const node &)) 
{
 
  std::vector<node_it> dummy_list;
  dummy_list.reserve(nodes.size());
  for (node_it n = nodes.begin();n!=nodes.end();n++) dummy_list.push_back(n);
  disperse::ompPSort(dummy_list.begin(),dummy_list.end(),
		     glob_num_omp_threads, node::compareItLess());
  node_it n=nodes.begin();
  for (std::vector<node_it>::iterator it=dummy_list.begin();it!=dummy_list.end();it++)
    {
      nodes.splice(n,nodes,*it);
      n=*it;
      n++;
    }
  dummy_list.clear();

}

inline NDcomplex::arc_it NDcomplex::
addArc(NDcomplex::node_it n1, NDcomplex::node_it n2, int pool)
{
  //printf("add %d\n",pool);
  NDcomplex::arc_it it;
  if (n1->getType()<n2->getType())
    {
      if (pool)
	it=arcPool[pool-1].insert(arcPool[pool-1].end(),arc(n1,n2,arcs.end(),arcsGeom.end()));
      else
	it=arcs.insert(arcs.end(),arc(n1,n2,arcs.end(),arcsGeom.end()));

#pragma omp critical 
      {
	n1->insertArc(it,n1,arcs.end());
	n2->insertArc(it,n2,arcs.end());
      }
    }	  
  else
    {
      if (pool)
	it=arcPool[pool-1].insert(arcPool[pool-1].end(),arc(n2,n1,arcs.end(),arcsGeom.end()));
      else
	it=arcs.insert(arcs.end(),arc(n2,n1,arcs.end(),arcsGeom.end()));
	    
#pragma omp critical 
      {
	n1->insertArc(it,n1,arcs.end());
	n2->insertArc(it,n2,arcs.end());
      }
    }
	  
  return it;
}

inline void NDcomplex::
restoreRemovedArcsLinks(double robustnesT)
{
  nodes.splice(nodes.end(),deleted_nodes,deleted_nodes.begin(),deleted_nodes.end());

  for (arc_it arc=rm_arcs_begin();arc!=rm_arcs_end();arc++)
    {
      node_it n1 = arc->endpoints().first;
      node_it n2 = arc->endpoints().second;
      n1->insertArc(arc,n1,arcs.end());
      n2->insertArc(arc,n2,arcs.end());
    }
}

inline void NDcomplex::
unrestoreRemovedArcsLinks()
{ 
  for (arc_it arc=rm_arcs_begin();arc!=rm_arcs_end();arc++)
    arc->unLink(arc,arcs.end()); 
  
  for (node_it node=nodes_begin();node!=nodes_end();)
    if (node->isDeleted()) 
      node=removeNode(node,true);
    else node++;
}

inline NDcomplex::arc_it NDcomplex::
removeArc(NDcomplex::arc_it &it, bool storeDeleted)
{
  bool store=storeDeleted;
  it->unLink(it,arcs.end());

  if (store)
    {
      arcGeom_it g=it->getGeometry();
      if (g==arcsGeom_end()) store=false;
      else if (g->getNRefs()>1) store=false;
    }

  arcGeom_it g=it->deleteGeometry(arcsGeom_end());
  if (g!=arcsGeom_end()) arcsGeom.erase(g);
  return arcs.erase(it);
  
}

inline NDcomplex::node_it NDcomplex::
addNode(cellType cell, double val,double val2, char flags,char type, int pool)
{       
  char nodeType;

  if (type<0) nodeType=cell.type();
  else nodeType=type;

  if (pool)
    return nodePool[pool-1].insert(nodePool[pool-1].end(),node(cell,nodeType,val,arcs.end(),nodes.end(),nodesGeom.end(),val2,flags));
  else
    return nodes.insert(nodes.end(),node(cell,nodeType,val,arcs.end(),nodes.end(),nodesGeom.end(),val2,flags));
}

inline NDcomplex::node_it NDcomplex::
removeNode(NDcomplex::node_it node, NDcomplex::arc_it &arc, bool storeDeleted)
{
  NDcomplex::arc_it it=node->getArc();

  while (it!=arcs.end())
    {
      if (it==arc)
	arc=removeArc(it,storeDeleted);
      else
	removeArc(it,storeDeleted);
      it=node->getArc();
    }

  nodeGeom_it g;
  g=node->deleteUpGeometry(nodesGeom_end());
  if (g!=nodesGeom_end()) nodesGeom.erase(g);
  g=node->deleteDownGeometry(nodesGeom_end());
  if (g!=nodesGeom_end()) nodesGeom.erase(g);

  if (storeDeleted)
    {
      node_it next=node;
      next++;
      node->setDeleted();
      deleted_nodes.splice(deleted_nodes.end(),nodes,node);
      return next;
    }
  else return nodes.erase(node);
}

inline NDcomplex::node_it NDcomplex::
removeNode(NDcomplex::node_it node, bool storeDeleted)
{
  NDcomplex::arc_it arc=arcs.end();
  return removeNode(node,arc,storeDeleted);
}

inline void NDcomplex::
createPools(int N)
{
  nodesGeom.createPools(N);
  arcsGeom.createPools(N);

  //printf("Creating %d pools\n",N);
  if (arcPool.size()) {
    fprintf(stderr,"ERROR in NDcomplex::createPools : pools are not empty\n");
    exit(0);
  }
  if (N>1) arcPool.resize(N-1);
  
  if (nodePool.size())  {
    fprintf(stderr,"ERROR in NDcomplex::createPools : pools are not empty\n");
    exit(0);
  }
  if (N>1) nodePool.resize(N-1);
  
}

inline void NDcomplex::
commitPools()
{
  nodesGeom.commitPools();
  arcsGeom.commitPools();

  int i;
 
  for (i=0;i<(int)arcPool.size();i++)
    arcs.splice(arcs.end(),arcPool[i]);
  arcPool.clear();

  for (i=0;i<(int)nodePool.size();i++)
    nodes.splice(nodes.end(),nodePool[i]);
  nodePool.clear();
}

template <class InputIterator>
inline NDcomplex::node::geom_it NDcomplex::
insertNodeGeometry(InputIterator begin, InputIterator end, char type,int pool)
{
  return nodesGeom.insertInPool(begin,end,type,pool);
}

template <class InputIterator>
inline NDcomplex::arc::geom_it NDcomplex::
insertArcGeometry(InputIterator begin, InputIterator end,int pool)
{
  return arcsGeom.insertInPool(begin,end,0,pool);
}

inline NDcomplex::node::geom_it NDcomplex::
insertNodeNullGeometry( char type, bool store_manifolds)
{
  if (!store_manifolds) return nodesGeom_end();

  std::vector<cellType> dummy; 
  return nodesGeom.insert(dummy.begin(),dummy.end(),type);
}

inline NDcomplex::arc::geom_it NDcomplex::
setArcNullGeometry(arc_it a, bool store_geometry)
{
  if (!store_geometry) return arcsGeom_end();

  //std::vector<cellType> dummy(500,a->endpoints().second->getCell());
  std::vector<cellType> dummy(1,a->endpoints().first->getCell());
  arcGeom_it ag=arcsGeom.insert(dummy.begin(),dummy.end(),0);
  a->setManifoldGeometry(ag);
  return ag;
}

inline NDcomplex::arc::geom_it NDcomplex::
duplicateArcGeometry(arcGeom_it ag)
{
  return arcsGeom.duplicate(ag);
}

inline void NDcomplex::buildReflectiveBoundary()
{
  node_it mnode;
  arc_it marc;
  long i;
  std::map<node_it,node_it,node::compareItLess> node_map;

  nodeGeom_it ng;
  arcGeom_it ag;
  std::vector<cellType> dummyCells;
  node_it nit;

  for (nit=nodes.begin();nit!=nodes.end();nit++)
    if (nit->isBoundary()) break;

  if (nit==nodes.end())
    {
      printf("Manifold has no boundaries, skipping boundary conditions.\n");
      return;
    }
  //removeOutNodes();
  printf("Building reflective boundaries ... ");fflush(0);

  
  
  printf("(nodes) ");fflush(0);
  long nnodes=nodes.size();
  i=0;
  for (node_it it=nodes.begin();i<nnodes;it++,i++)
    {
      if (it->isBoundary()) continue;
      
      mnode=nodes.insert(nodes.end(),node(*it,arcs.end(),nodes.end(),nodesGeom.end()));
     
      ng=it->getUpGeometry();
      if (ng != nodesGeom_end())
	{
	  ng = insertNodeGeometry(dummyCells.begin(), dummyCells.end(),it->getType());
	  mnode->setManifoldGeometry(ng,true);
	}
      ng=it->getDownGeometry();
      if (ng != nodesGeom_end())
	{
	  ng = insertNodeGeometry(dummyCells.begin(), dummyCells.end(),it->getType());
	  mnode->setManifoldGeometry(ng,false);
	}

      mnode->flagsOr(NODE_FLAG_OUT);
      node_map.insert(std::make_pair(it,mnode));
    }

  printf("(arcs) ");fflush(0);
  long narcs=arcs.size();
  i=0;
  for (arc_it it=arcs.begin();i<narcs;it++,i++)
    {
      std::pair<node_it,node_it> endpts = it->endpoints();
      node_it n1,n2;
      int nb=0;
      if (endpts.first->isBoundary()) {n1=endpts.first;nb++;}
      else n1=node_map[endpts.first];

      if (endpts.second->isBoundary()) {n2=endpts.second;nb++;}
      else n2=node_map[endpts.second];

      if (nb==2) continue;

      marc = addArc(n1,n2);

      ag = it->getGeometry();
      if (ag != arcsGeom_end())
	{
	  ag = insertArcGeometry(dummyCells.begin(),dummyCells.end());
	  marc->setManifoldGeometry(ag);
	}
    }
  
  printf("done. (+%ld n./+%ld a.)\n",nodes.size()-nnodes,arcs.size()-narcs);
}


inline void NDcomplex::
sanitizeBoundary(bool vertexAsMaxima, int store_manifolds, int store_arcsGeom)
{
  typedef std::map<node_it,node_it,node::compareItLess> nodeMap;
  typedef nodeMap::iterator nodeMap_it;

  node_it cur;
  long i;
  long count=0;
  long deltaIndex=vertexAsMaxima?-1:1;  
  std::vector<long> nadd(ndims+1,0);
  long narcs=0;
  nodeMap boundary;

  printf("Sanityzing complex boundary ... ");fflush(0);

  for (cur=nodes.begin();cur!=nodes.end();cur++)
    {
      if (!cur->isBoundary()) continue;
      
      nodeGeom_it ng;
      arcGeom_it ag;
      node_it newNode = addNode(cur->getCell(),cur->getVal(),cur->getVal2(),
				NODE_FLAG_OUT,cur->getType()+deltaIndex);
   
      ng = insertNodeNullGeometry(newNode->getType(),store_manifolds);
      newNode->setManifoldGeometry(ng,true);
      assert((!store_manifolds)||(ng->getType()>=0));
      assert((!store_manifolds)||(ng->getType()<=ndims));
 

      ng = insertNodeNullGeometry(newNode->getType(),store_manifolds);
      newNode->setManifoldGeometry(ng,false);
      assert((!store_manifolds)||(ng->getType()>=0));
      assert((!store_manifolds)||(ng->getType()<=ndims));


      int maxT = std::max(newNode->getType(),cur->getType());
      bool store=false;

      if ((maxT==ndims_net)&&(store_arcsGeom&(1<<0))) store=true;
      else if ((maxT==1)&&(store_arcsGeom&(1<<2))) store=true;
      else if (store_arcsGeom&(1<<1)) store=true;
      
      arc_it newArc=addArc(cur,newNode);
      
      setArcNullGeometry(newArc,store);
      nadd[newNode->getType()]++;
      count++;
      
      boundary[cur]=newNode;      
    }
  
  for (nodeMap_it it=boundary.begin();it!=boundary.end();it++)
    {
      node_it cur = it->first;
      node_it img = it->second;

      arc_it a=cur->getArc();
      
      do {
	node_it n=a->getOtherNode(cur);
	
	if ((!n->isBoundary())||
	    (n->getType()>cur->getType()))
	  {
	    a = a->getNextRef(cur);
	    continue;
	  }
	
	arc_it newArc= addArc(boundary[n],img);
	arcGeom_it ag = a->getGeometry();

	bool store=false;
	int maxT = std::max(newArc->endpoints().first->getType(),newArc->endpoints().second->getType());
	if ((maxT==ndims_net)&&(store_arcsGeom&(1<<0))) store=true;
	else if ((maxT==1)&&(store_arcsGeom&(1<<2))) store=true;
	else if (store_arcsGeom&(1<<1)) store=true;

	if (ag==arcsGeom_end())
	  {
	    setArcNullGeometry(newArc,store);
	  }
	else
	  {
	    if (!store) setArcNullGeometry(newArc,store);
	    else newArc->setManifoldGeometry(duplicateArcGeometry(ag));
	  }
	//assert(newArc->getGeometry()!=arcsGeom_end());
	
	narcs++;
	a = a->getNextRef(cur);
      } while(a!=arcs_end());            
    }
      
    
  printf("done.\n");
  printf("   Dummy nodes : ");
  for (i=0;i<=ndims_net;i++)
    {     
      printf("%ld %s%s",nadd[i],nodeType2Str(i).c_str(),(ndims==i)?".\n":", ");
    }
  printf("   Dummy arcs : ");
  printf("%ld from boundary, %ld at infinity.\n",count,narcs);
    
}



inline void NDcomplex::
removeBoundaries(bool removeAll)
{
  long nRemoved=0;
    printf ("Removing boundaries ... ");fflush(0);
    
    if (removeAll) printf ("(full) ");fflush(0);
    
    for (node_it it=nodes_begin();it != nodes_end();)
    {
      bool remove=true;

      if (it->isOut()) remove=true;
      else if (it->isBoundary())
	{
	  if (!removeAll)
	    {
	      std::vector<node_it> nei;
	      std::vector<node_it>::iterator nei_it;
	      
	      it->getNeighbors(it,arcs_end(),back_inserter(nei));
	      for (nei_it=nei.begin();nei_it!=nei.end();nei_it++)
		{
		  if ((!(*nei_it)->isOut())&&
		      (!(*nei_it)->isBoundary())&&
		      (!(*nei_it)->isInfinite())) remove=false;
		}
	    }	  
	}
      else remove=false;

      //if (!it->isOut()) remove=false;
      if (remove) 
	{
	  node_it pn= it->getPNode();
	  if ((pn!=nodes_end())&&(pn!=it))
	    pn->setPNode(pn);
	  
	  it=removeNode(it,false);
	  nRemoved++;
	}
      else it++;

    }
    
    printf("done. (%ld nodes removed)\n",nRemoved);
}

inline bool NDcomplex::
haveBoundaries()
{
  node_it it;

  for (it=nodes_begin();it != nodes_end();it++)
    {
      if (it->isBoundary()) break;
    }

  return (it!=nodes_end());
}

inline void NDcomplex::
removeOutNodes()
{
  long nRemoved=0;
  printf ("Removing out nodes ... ");fflush(0);
  
  
  for (node_it it=nodes_begin();it != nodes_end();)
    {
      bool remove=false;
      
      if (it->isOut()) remove=true;
      
      //if (!it->isOut()) remove=false;
      if (remove) 
	{
	  node_it pn= it->getPNode();
	  if ((pn!=nodes_end())&&(pn!=it))
	    pn->setPNode(pn);
	  
	  it=removeNode(it,false);
	  nRemoved++;
	}
      else it++;
      
    }
  
  printf("done. (%ld nodes removed)\n",nRemoved);
}

NDcomplex::node_it NDcomplex::
getParentNode(node_it node, bool useRatio)
{
  if (!getCyclesState()[0])
    {
      std::vector<int> wh(1,0);
      fprintf(stderr,"WARNING: asked for parent node but 0-cycles are not available.\n");
      fprintf(stderr,"         will compute them ...\n");
      computePersistencePairs(wh);
    }

  bool down=(node->getType()==0);
  if ((!down)&&(node->getType()!=ndims_net)) return nodes_end();
  
  if (node->getPNode()==node) return nodes_end();
  double p_ref=(useRatio)?(node->persistenceR()):node->persistence();
  node_it n=node;
      
  while(true)
    {
      node_it res=(down)?n->getRefDown():n->getRefUp();
      
      if (res==nodes_end()) return nodes_end();
      if (res->isDeleted())
	{
	  n=res;
	  continue;
	}
      if (res->getPNode()==res) return res;
 
      double p=(useRatio)?(res->persistenceR()):res->persistence();
      if (p<p_ref) n=res;
      else return res;
    }
 
  // return nodes_end();
}

inline std::pair<NDcomplex::node_it,NDcomplex::node_it> NDcomplex::
getReferenceNodes(arc_it arc)
{
  std::vector<bool> cyclesState = getCyclesState(); 
  if (cyclesState[1]) return getReferenceNodes_cycles(arc);

  std::vector<int> wh(1,1);
  fprintf(stderr,"WARNING: asked for reference nodes but 1-cycles are not available.\n");
  fprintf(stderr,"         will compute them ...\n");
  computePersistencePairs(wh);
  return getReferenceNodes_cycles(arc);

  //return getReferenceNodes_neighbors(arc);
}

inline std::pair<NDcomplex::node_it,NDcomplex::node_it> NDcomplex::
getReferenceNodes_cycles(arc_it arc)
{
  node_it nd=arc->endpoints().first;
  node_it nu=arc->endpoints().second;
  int reverse=false;
  std::pair<NDcomplex::node_it,NDcomplex::node_it> result;

  if (nd->getType()>nu->getType()) 
    {
      std::swap(nd,nu);
      reverse=true;
    }
  
  if (nd->getType()==0)
    result=std::make_pair(nd->getRefUp(),nu->getRefUp());
  else if (nu->getType()==ndims_net)
    result=std::make_pair(nd->getRefDown(),nu->getRefDown());
  else 
    result=std::make_pair(nodes.end(),nodes.end());

  if (reverse) std::swap(result.first,result.second);

  return result;
}

inline bool NDcomplex::
computePersistencePairs(const std::vector<int> compCycles)
{
  int i;
  //havePPairs=false;
  //haveCycles.clear();
  
  if (havePPairs)
    { 
      for (i=0;i<compCycles.size();i++)
	{
	  if (haveCycles.size()<compCycles[i]) break;
	  if (!haveCycles[compCycles[i]]) break;
	}

      if (i==compCycles.size())
	{
	  printf("Computing persistence pairs and cycles ... SKIPPED.\n");
	  return false;
	}
    }
 
  computePersistencePairs_matrix(compCycles);
  
  havePPairs=true;
  haveCycles.resize(ndims_net+1);
 
  for (i=0;i<compCycles.size();i++)
    haveCycles[compCycles[i]]=true;
    
  return true;
}


inline sparseZ2Matrix* NDcomplex::
buildBoundaryMatrix(std::vector<node_it> &nodeTab, char curType, bool reverse, bool skipAlreadyPaired)
{
  std::map<node_it,long,node::compareItLess> nodeID;
  long i,j;
  long nnodes = nodeTab.size();
  sparseZ2Matrix *Mp = new sparseZ2Matrix(nnodes,nnodes);
  sparseZ2Matrix &M = *Mp;
  
  if (reverse)
    {
      for (i=0;i<nodeTab.size();i++)
	nodeID.insert(std::make_pair(nodeTab[i],nnodes-i-1));
    }
  else
    {
      for (i=0;i<nodeTab.size();i++)
	nodeID.insert(std::make_pair(nodeTab[i],i));
    }
  
  std::vector<node_it> nei;
  Z2set<long> rem;
  Z2set<long>::iterator rem_it;
  std::vector<node_it>::iterator nei_it;
  bool ninf=false;
  int delta;
  
  if (reverse) delta=1;
  else delta=-1;
   
  for (j=0;j<nnodes;j++)
    {
      i=(reverse)?nnodes-j-1:j;
      node_it n = nodeTab[j];
      if ((n->getType()==curType)&&(!n->isInfinite()))
	{
	  if (skipAlreadyPaired)
	    {
	      if (n->getPNode() != nodes.end()) 
		{
		  if (n->getPNode()->getType()!=n->getType()+delta)
		    continue;
		}
	    }

	  nei.clear();
	  n->getNeighbors(n,arcs_end(),back_inserter(nei),n->getType()+delta);
	  for (nei_it=nei.begin();nei_it!=nei.end();nei_it++)
	    {
	      node_it n_it = *nei_it;
	      
	      if (n_it->isInfinite()) continue;
	      if (n_it->isOut()!=n->isOut()) continue;
	      
	      bool res;
	      long nid=nodeID[n_it];
	      //#pragma omp critical
	      res = M.set(nid,i);
	      if (res) rem.insert(nid);
	    }
	  
	  if (rem.size())
	    {
	      //#pragma omp critical
	      for (rem_it=rem.begin();rem_it!=rem.end();rem_it++)
		M.set(*rem_it,i,false);
	      rem.clear();
	    }
	}
    }

  return Mp;
}

inline bool NDcomplex::
setCyclesFromBoundaryMatrix(int cycle_type,sparseZ2Matrix &M, const std::vector<node_it> &nodeTab, bool reverse)
{     
  long ntot=nodeTab.size();
  typedef sparseZ2Matrix::colType colType;
  typedef sparseZ2Matrix::colItType colType_it;
  int i=0;

  if (cycle_type==0)
    {
      
      for (long ii=0;ii<ntot;ii++)
	{
	  long i2=(reverse)?ntot-ii-1:ii;     
	  long ll=M.low(i2); 
	  if (ll<0) continue;
	  node_it n=nodeTab[ii];
	  long hh=M.high(i2); 	  

	  long l2=(reverse)?ntot-ll-1:ll;	  
	  node_it low=nodeTab[l2];
	  long h2=(reverse)?ntot-hh-1:hh;	  
	  node_it high=nodeTab[h2];
	  
	  if (low==high) continue;

	  if (low->getType()>n->getType())
	    low->setRefUp(high);
	  else
	    low->setRefDown(high);	  
	}    
    }

  if (cycle_type==1)
    {
      boundaryMatrix(&M).canonize();

      for (long ii=0;ii<ntot;ii++)
	{
	  long i2=(reverse)?ntot-ii-1:ii;      
	  colType col=M.getCol(i2);
	  if (M.colIsVoid(col)) continue;
	  node_it n=nodeTab[ii];

	  for (colType_it it=col->begin();it!=col->end();it++)
	    {
	      long ll=**it;
	      long l2=(reverse)?ntot-ll-1:ll;	  
	      node_it p=nodeTab[l2];
	  
	      if (p->getType()>n->getType())
		{
		  node_it ref_down=p->getRefDown();
		  if (ref_down==nodes.end()) p->setRefDown(n);
		  else if (ref_down->getVal()>n->getVal()) p->setRefDown(n);
		}
	      else
		{
		  node_it ref_up=p->getRefUp();
		  if (ref_up==nodes.end()) p->setRefUp(n);
		  else if (ref_up->getVal()<n->getVal()) p->setRefUp(n);
		}
	    }
	}

      if (debug_dump)
	{
	  FILE *f;
	  char full_fname[255];
	  char dbg_fname[255];
      
	  if (reverse)
	    strcpy(dbg_fname,"cycles_up");
	  else
	    strcpy(dbg_fname,"cycles_down");

	  sprintf(full_fname,"%s.i",dbg_fname);
	  f=fopen(full_fname,"w");
	  fprintf(f,"func get_%s(img,x0,y0,x1,y1) \n{\n",dbg_fname);
	  fprintf(f,"all=[];ndims=%d;\n",ndims);
	  fprintf(f,"d2p=[1,3,2];\n");
	  fprintf(f,"if (is_void(x0)) x0=double(0); else x0=double(x0);\n");
	  fprintf(f,"if (is_void(x1)) x1=double(0); else x1=double(x1);\n");
	  fprintf(f,"if (is_void(y0)) y0=double(1); else y0=double(y0);\n");
	  fprintf(f,"if (is_void(y1)) y1=double(1); else y1=double(y1);\n");
	  fprintf(f,"dims=dimsof(img)(2:);dx=(x1-x0)/dims(1);dy=(y1-y0)/dims(2);\n");
      
	  for (long ii=0;ii<ntot;ii++)
	    {
	      long i2=(reverse)?ntot-ii-1:ii;      
	      colType col=M.getCol(i2);
	      if (M.colIsVoid(col)) continue;
	      node_it n=nodeTab[ii];
	  
	      fprintf(f,"c=[[%ld,%ld,%ld,%ld]];\n",
		      (long)n->getPNode()->getCell().type(),(long)n->getPNode()->getCell().id(),
		      (long)n->getCell().type(),(long)n->getCell().id());
	  
	      for (colType_it it=col->begin();it!=col->end();it++)
		{
		  long ll=**it;
		  long l2=(reverse)?ntot-ll-1:ll;	  
		  node_it p=nodeTab[l2];
	      
		  long ct;
		  std::vector<node_it> nei;
		  p->getNeighbors(p,arcs_end(),back_inserter(nei),(reverse)?(p->getType()+1):(p->getType()-1));
		  for (ct=0;ct<nei.size();ct++)
		    {
		      fprintf(f,"grow,c,[[%ld,%ld,%ld,%ld]];\n",
			      (long)nei[ct]->getCell().type(),(long)nei[ct]->getCell().id(),
			      (long)p->getCell().type(),(long)p->getCell().id());
		    }
	      
		}
	      fprintf(f,"d=[c(2,)/d2p(1+c(1,)),c(4,)/d2p(1+c(3,))];\n");
	      fprintf(f,"d=[x0+dx*(d(,1)%%dims(1)),y0+dy*(d(,1)/dims(1)),x0+dx*(d(,2)%%dims(1)),y0+dy*(d(,2)/dims(1))];\n");
	  
	      fprintf(f,"grow,all,&d;\n");//for regular grid only
	    }
	  fprintf(f,"return all;\n}\n");
	  fprintf(f,"func plot_%s(all,w=,width=,delta=) \n{\n",dbg_fname);
	  fprintf(f,"if (is_void(w)) w=indgen(numberof(all));\n");
	  fprintf(f,"col=char(random(3,numberof(w))*255);\n");
	  fprintf(f,"for (i=1;i<=numberof(w);i++) {c=*all(w(i));\n");
	  fprintf(f,"c+=cos([c(,1)*c(,2)*1.E10,c(,1)+c(,2)*1.E10,c(,3)*c(,4)*1.E10,c(,3)+c(,4)*1.E10]+100*i)*delta;\n");
	  fprintf(f,"pldj,c(,1),c(,2),c(,3),c(,4),color=col(,i),width=width;}\n}\n");  
	  fclose(f);
	}
    }

  return true;
}

inline void NDcomplex::
propagateReferenceFrom1Cycles()
{
  std::list<arc_it> upOut;
  std::list<arc_it> downOut;
  for (arc_it it=arcs_begin();it!=arcs.end();it++)
    {
      std::pair<node_it,node_it> p=it->endpoints();
      node_it up=p.first;
      node_it down=p.second;
      if (up->getType()<down->getType()) std::swap(up,down);

      if (down->getType()==0)
	{
	  if (up->getRefUp()==nodes.end()) 
	    {
	      downOut.push_back(it);
	      continue;
	    }

	  if (down->getRefUp()==nodes.end()) 
	    down->setRefUp(up->getRefUp());
	  else if (down->getRefUp()->getVal()<up->getRefUp()->getVal())
	    down->setRefUp(up->getRefUp());
	}

      if (up->getType()==ndims_net)
	{
	  if (down->getRefDown()==nodes.end()) 
	    {
	      upOut.push_back(it);
	      continue;
	    }

	  if (up->getRefDown()==nodes.end()) 
	    up->setRefDown(down->getRefDown());
	  else if (up->getRefDown()->getVal()>down->getRefDown()->getVal())
	    up->setRefDown(down->getRefDown());
	}	    
    }
  long len=upOut.size();
  while (len)
    {
      for (std::list<arc_it>::iterator it=upOut.begin();it!=upOut.end();)
	{
	  std::pair<node_it,node_it> p=(*it)->endpoints();
	  node_it up=p.first;
	  node_it down=p.second;
	  if (up->getType()<down->getType()) std::swap(up,down);
	  if (up->getRefDown()==nodes.end()) 
	    {
	      if (down->getRefDown()==nodes.end())
		it++;
	      else
		{
		  up->setRefDown(down->getRefDown());
		  it=upOut.erase(it);
		}
	      continue;
	    }
	  down->setRefDown(up->getRefDown());
	  it=upOut.erase(it);	      
	}
	  
      if (len!=upOut.size()) len=upOut.size();
      else {len=0;fprintf(stderr,"\nWARNING: Isolated components (up 1-cycles) ...\n");}
    }

  len=downOut.size();
  while (len)
    {
      for (std::list<arc_it>::iterator it=downOut.begin();it!=downOut.end();)
	{
	  std::pair<node_it,node_it> p=(*it)->endpoints();
	  node_it up=p.first;
	  node_it down=p.second;
	  if (up->getType()<down->getType()) std::swap(up,down);
	  if (down->getRefUp()==nodes.end()) 
	    {
	      if (up->getRefUp()==nodes.end()) 
		it++;
	      else
		{
		  down->setRefUp(up->getRefUp());
		  it=downOut.erase(it);
		}
	      continue;
	    }
	  up->setRefUp(down->getRefUp());
	  it=downOut.erase(it);	      
	}
      if (len!=downOut.size()) len=downOut.size();
      else {len=0;fprintf(stderr,"\nWARNING: Isolated components (down 1-cycles) ...\n");}
    }

}

inline long NDcomplex::
setPairsFromBoundaryMatrix(const sparseZ2Matrix &M, const std::vector<node_it> &nodeTab, bool reverse)
{
  long npb=0;
  long ntot=nodeTab.size();
  
  for (long ii=0;ii<ntot;ii++)
    {
      long i2=(reverse)?ntot-ii-1:ii;     
      long ll=M.low(i2);  
      if (ll<0) continue;
      node_it n=nodeTab[ii];

      long l2=(reverse)?ntot-ll-1:ll;	  
      node_it p=nodeTab[l2];
       
      // This should NEVER happen ... 
      if ((n->getPNode() != nodes.end())||(p->getPNode() != nodes.end()))
	{	      
	  npb++;
	  printf("WARNING: undetermined pairing:\n");	      
	      
	  if (n->getPNode() != nodes.end())
	    {
	      if (p->getPNode() != nodes.end()) continue;
	      if (n->persistence()<n->persistence(p)) continue;
	      node_it tmp_node = n->getPNode();
	      tmp_node->setPNode(tmp_node);
	    }
	      
	  if (p->getPNode() != nodes.end())
	    {
	      if (p->persistence()<p->persistence(p)) continue;
	      node_it tmp_node = p->getPNode();
	      tmp_node->setPNode(tmp_node);
	    }
	      
	  n->flagsAnd(~(NODE_FLAG_NEGATIVE|NODE_FLAG_POSITIVE));
	  p->flagsAnd(~(NODE_FLAG_NEGATIVE|NODE_FLAG_POSITIVE));
	}

      n->setPNode(p);
      p->setPNode(n);
	  
      if (n->getType()>p->getType())
	{
	  n->flagsOr(NODE_FLAG_POSITIVE);
	  p->flagsOr(NODE_FLAG_NEGATIVE);
	}
      else
	{
	  p->flagsOr(NODE_FLAG_POSITIVE);
	  n->flagsOr(NODE_FLAG_NEGATIVE);
	}
	  
      //printf("[%d,%d] : %s - %s\n",nodeID[n->uniqueID()],nodeID[p->uniqueID()],
      //n->getInfo(true).c_str(),p->getInfo(true).c_str());
	  
    }
  return npb;
}
/*
template <class Arg, class Result>
class index2Type : public std::unary_function<Arg, Result> {
private:
  std::vector<NDcomplex::node_it> *nodeVec;
  
public:
  
  index2Type(std::vector<NDcomplex::node_it> *nodeV) {
    nodeVec=nodeV;
  }
  
  Result operator() (Arg id) {return (*nodeVec)[id]->getType();}
};
*/

inline void NDcomplex::
computePersistencePairs_matrix(const std::vector<int> compCycles_p)
{
  long nInNodes;
  long i,j;
  long nnodes = num_nodes();
  std::vector<node_it> nodeTab;
  int npb=0; 
  char txt[1024];
  std::vector<int> compCycles;

  struct timeval wc_start,wc_stop;
  gettimeofday(&wc_start, NULL);

  for (i=0;i<compCycles_p.size();i++)
    {
      if (compCycles_p[i]>=haveCycles.size()) compCycles.push_back(compCycles_p[i]);
      else if (!haveCycles[compCycles_p[i]]) compCycles.push_back(compCycles_p[i]);
    } 
    
  strcpy(txt,"will compute");
  if (!havePPairs) sprintf(txt,"%s persistence pairs",txt); 

  if (compCycles.size())
    {
      i=0;
      if (!havePPairs) sprintf(txt,"%s and (%d",txt,compCycles[i++]);
      else sprintf(txt,"%s (%d",txt,compCycles[i++]);
      for (i=1;i<compCycles.size();i++)
	sprintf(txt,"%s,%d",txt,compCycles[i]);
      sprintf(txt,"%s)-cycles",txt);
    }
  if (havePPairs&&(compCycles.size()==0)) 
    {
      printf("Computing persistence pairs ... skipped.\n");fflush(0);
      return;
    }  

  int nth=1;
#ifdef USE_THREADS
  nth=glob_num_threads;
#endif
  nth=1;
  printf("%s (%dT):\n",txt,nth);fflush(0);
  printf("    Setup  ...");fflush(0);
  //printf("(setup) ");fflush(0); 
 
  //std::sort(nodeTab.begin(),nodeTab.end(),node::compareItLess());
  
  nodeTab.resize(nnodes);
  
  bool resetRef1=false;
  bool resetRef0=false;
  for (i=0;i<compCycles.size();i++) 
    {
      if (compCycles[i]==1) resetRef1=true;
      if (compCycles[i]==0) resetRef0=true;
    }

  
  i=0;
  for (node_it n = nodes.begin();n!=nodes.end();n++,i++)
    {
      nodeTab[i]=n;
      if (!havePPairs)
	{
	  n->setPNode(nodes.end());
	  n->flagsAnd(~(NODE_FLAG_POSITIVE|NODE_FLAG_NEGATIVE));     
	}      
      
      //assert(n->getRefUp()==nodes.end());
      //assert(n->getRefDown()=nodes.end());

      if (resetRef1)
	{
	  if (n->getType()!=ndims_net) n->setRefUp(nodes.end());
	  //else n->setRefDown(nodes.end());
	  if (n->getType()!=0) n->setRefDown(nodes.end());
	  //else n->setRefUp(nodes.end());
	}
      if (resetRef0)
	{
	  if (n->getType()==ndims_net) n->setRefUp(nodes.end());
	  if (n->getType()==0) n->setRefDown(nodes.end());
	}
    }

  //#ifdef USE_OPENMP
#ifdef USE_GNU_PSORT
  __gnu_parallel::sort(nodeTab.begin(),nodeTab.end(),node::compareItLess());
#else
  std::sort(nodeTab.begin(),nodeTab.end(),node::compareItLess());
#endif

  //printf("\rComputing persistence pairs (%2.2dT) ... (setup) (building) ",nth);fflush(0);
  
  
  //std::vector< std::pair<char,bool> > type_dir;
  std::vector<matrix_parms> mType;
  std::vector<bool> isSaddleSaddle;
  
  if (!havePPairs)
    {      
      mType.push_back(matrix_parms(ndims_net-1,true,true,false));
      for (j=ndims_net-1;j>0;j--)
	mType.push_back(matrix_parms(j,false,true,false));
    }

  for (i=0;i<compCycles.size();i++)
    {
      matrix_parms prm=matrix_parms(ndims_net-compCycles[i]-1,true,false,true);
      for (j=0;j<mType.size();j++)
	{
	  if ((mType[j].crit_index==prm.crit_index)&&(mType[j].up==prm.up))
	    {
	      mType[j].pairs=mType[j].cycles=true;break;
	    }
	  if ((mType[j].crit_index==prm.crit_index+1)&&(mType[j].up!=prm.up)&&(mType[j].cycles==false))
	    {
	      mType[j].cycles=true;
	      mType[j].crit_index--;
	      mType[j].up=true;
	      break;
	    }
	}
      if (j==mType.size()) mType.push_back(prm);

      prm=matrix_parms(compCycles[i]+1,false,false,true);
      for (j=0;j<mType.size();j++)
	{
	  if ((mType[j].crit_index==prm.crit_index)&&(mType[j].up==prm.up))
	    {
	      mType[j].pairs=mType[j].cycles=true;
	      break;
	    }
	  if ((mType[j].crit_index==prm.crit_index-1)&&(mType[j].up!=prm.up)&&(mType[j].cycles==false))
	    {
	      mType[j].cycles=true;
	      mType[j].crit_index++;
	      mType[j].up=false;
	      break;
	    }
	}
      if (j==mType.size()) mType.push_back(prm);
    }

  for (j=0;j<mType.size();j++) 
    {
      if (mType[j].crit_index==0) isSaddleSaddle.push_back(false);
      else if ((mType[j].crit_index==1)&&(!mType[j].up)) isSaddleSaddle.push_back(false);
      else if (mType[j].crit_index==ndims_net) isSaddleSaddle.push_back(false);
      else if ((mType[j].crit_index==ndims_net-1)&&(mType[j].up)) isSaddleSaddle.push_back(false);
      else isSaddleSaddle.push_back(true);
    }
  
 

  std::vector<int> nThreads(mType.size(),1);
 
  for (j=0;j<mType.size();j++) 
    {
      if ((mType[j].crit_index==ndims_net-1)&&(mType[j].up)) nThreads[j]=1;
      else if ((mType[j].crit_index==1)&&(!mType[j].up)) nThreads[j]=1;
      else nThreads[j]=nth;
    }
  


  printf(" done.\n");fflush(0);
  
  for (int nPass=1;nPass<=2;nPass++)
    {
      printf("    Pass %d: (building) ",nPass);
#pragma omp parallel for 
      for (j=0;j<mType.size();j++)
	{
	  bool skipPaired=true;
	 
	  if (isSaddleSaddle[j])
	    {
	      if (nPass==1) continue;
	      
	    }
	  else 
	    {
	      if (nPass==2) continue;
	     
	    }
	  
	  bool reverse=mType[j].up;
	  char curType=mType[j].crit_index;
	  sparseZ2Matrix *Mp=buildBoundaryMatrix(nodeTab,curType,reverse,skipPaired);
	  sparseZ2Matrix &M=*Mp;

#pragma omp critical
	  {
	    printf("(%s%s|%d-%d> ",
		   (mType[j].pairs)?"P":"",(mType[j].cycles)?"C":"",
		   curType,curType+(reverse?1:-1));
	    fflush(0);
	  }
      
	  boundaryMatrix(Mp).solve(curType+(reverse?0:-1),false,std::min(nThreads[j],6));	
	  if (mType[j].pairs) npb+=setPairsFromBoundaryMatrix(M,nodeTab,reverse);
	  if (mType[j].cycles) 
	    {
	      int ctype;
	      if (mType[j].up) ctype=ndims_net-curType-1;
	      else ctype=curType-1;

	    
	      setCyclesFromBoundaryMatrix(ctype,M,nodeTab,reverse);
	    }
	  delete(Mp);Mp=NULL;  
      
#pragma omp critical
	  {
	    printf("<%d-%d|%s%s) ",		 
		   curType,curType+(reverse?1:-1),
		   (mType[j].pairs)?"P":"",(mType[j].cycles)?"C":"");
	    fflush(0);
	  }
	}
      printf("\n");
    }

  std::vector<long> nPairs(ndims,0);
 
  for (node_it n = nodes.begin();n!=nodes.end();n++)
    {
      if (n->getPNode() == nodes.end()) {	
	n->setPNode(n);
      }
      if (n->getPNode()==n) nPairs[0]++;      
    }

  for (i=0;i<compCycles.size();i++) 
    {
      if (compCycles[i]==1) 
	propagateReferenceFrom1Cycles();
    }
  
  long outU=0,outD=0;
  long failU=0,failD=0;
  for (node_it n = nodes.begin();n!=nodes.end();n++)
    {
      if (n->getType()==0)
	{
	  if (n->getRefDown()==nodes.end()) outD++;
	  else if (n->getRefDown()->getType()!=0) failD++;
	}
      if (n->getType()==ndims_net)
	{
	  if (n->getRefUp()==nodes.end()) outU++;
	  else if (n->getRefUp()->getType()!=ndims_net) failU++;
	}
    }
  
  gettimeofday(&wc_stop, NULL);
  double dt  = (double)wc_stop.tv_sec + ((double)wc_stop.tv_usec)*1.E-6;
  dt -= (double)wc_start.tv_sec + ((double)wc_start.tv_usec)*1.E-6;

  printf("All done in %.2fs : %ld nodes paired, %ld free, %d undecided.\n",dt,nnodes-nPairs[0],nPairs[0],npb);
   
}


// return next valid arc (as "it" is deleted ...)
inline NDcomplex::arc_it NDcomplex::
cancelPair(NDcomplex::arc_it it, bool &canceled, std::list<node_it> &p1,std::list<node_it> &p2,std::list<node_it> &p3,NDcomplex::arc_it cur_it, bool force)
{
  std::pair<node_it,node_it> endpts = it->endpoints();
  std::list<arc_it> a1,a2;
  arc_it new_it = cur_it;
  bool impossible = false;
  char t1=endpts.first->getType();
  char t2=endpts.second->getType();
  std::vector<arc_it> multiArcs;

  p1.clear();
  p2.clear();
  p3.clear();
  //printf ("Cancelling %s\n",it->getInfo(true).c_str());

  // critical index difference must be 1
  if (abs(t1-t2)!=1)
    {
      fprintf (stderr,"Illegal critical cancellation !!!\n");
      exit(0);
      return ++new_it;
    }
  
  if ((endpts.second->isInfinite()&&(!endpts.first->isInfinite()))||
      (endpts.first->isInfinite()&&(!endpts.second->isInfinite())))
    {
      printf ("Impossible cancellation: Infinite link.\n");
      canceled = false;
      return new_it;
    }
  //mygettime(true);

  // check all arcs connected to node->first
  // and keep those with same critical index as endpts.second
  for (arc_it cur_arc = endpts.first->getArc();
       cur_arc!=arcs.end();
       cur_arc=cur_arc->getNext(endpts.first))
    {
      node_it cur_node=cur_arc->getOtherNode(endpts.first);
      if (cur_node == endpts.second)	
	{	
	  multiArcs.push_back(cur_arc);
	  if ((!force)&&(multiArcs.size()>1))
	    {
	      canceled=false;
	      return new_it;
	    }

	  continue;
	}
      
      if (cur_node->getType() == t2)
	{
	  p1.push_back(cur_node);
	  a1.push_back(cur_arc);
	}
      else p3.push_back(cur_node);
    }
  //mygettime();

  // check all arcs connected to node->second
  // and keep those with same critical index as endpts.first
  for (arc_it cur_arc = endpts.second->getArc();
       cur_arc!=arcs.end();
       cur_arc=cur_arc->getNext(endpts.second))
    {
      node_it cur_node=cur_arc->getOtherNode(endpts.second);
      if (cur_node == endpts.first) continue;
      
      if (cur_node->getType() == t1)
	{
	  p2.push_back(cur_node);
	  a2.push_back(cur_arc);
	}
      else p3.push_back(cur_node);
    }
  //mygettime();

  //create arcs from all nodes in p1 to all nodes in p2
  std::vector<arc_it>::iterator multi_it;
  int naddarcs=0;
  for (multi_it = multiArcs.begin();multi_it!=multiArcs.end(); multi_it++)
    {
      arcGeom_it ag = (*multi_it)->getGeometry();
      std::list<arc_it>::iterator ait1 = a1.begin();
      naddarcs+=p1.size()*p2.size();

      for (std::list<node_it>::iterator nit1 = p1.begin();nit1!=p1.end();nit1++,ait1++)
	{
	  std::list<arc_it>::iterator ait2 = a2.begin();
	  for (std::list<node_it>::iterator nit2 = p2.begin();nit2!=p2.end();nit2++,ait2++)
	    {
	      arc_it arc;
	      int nmarc;
	      // This parts limits the number of arcs between crit. point
	      // to 3 maximum (we only need to know whether the number is odd or even)
	      arc = (*nit1)->findMultiArc(*nit1,*nit2,arcs_end(),nmarc,3);
	      if (nmarc>2) 
		{
		  removeArc(arc);
		  continue;
		}
	      arc = addArc(*nit2,*nit1);
	      
	      if (ag!=arcsGeom_end())
    		{
		  arcGeom_it g1 = (*ait1)->getGeometry();
		  arcGeom_it g2 = (*ait2)->getGeometry();
		  
		  //remember: the order of the params is important
		  arc->setManifoldGeometry(arcsGeom.insert(g2,ag,g1));
    		}
	    }
	} 
    }
   
    // and now the nodes manifolds ...    
    node_it node_up,node_down;
    std::list<node_it> *pa,*pb;
    nodeGeom_it ng;

    if (t2>t1)
      {
	node_down=endpts.first;
	node_up=endpts.second;
	pa=&p1;
	pb=&p2;
      }
    else
      {
	node_up=endpts.first;
	node_down=endpts.second;
	pa=&p2;
	pb=&p1;
      }

    ng=node_up->getDownGeometry();      
    if (ng!=nodesGeom_end())
      {
	pa->sort(node::compareItLess());pa->unique();
	for (std::list<node_it>::iterator nit = pa->begin();nit!=pa->end();nit++)
	  {
	    nodeGeom_it g=(*nit)->getDownGeometry();
	    (*nit)->setManifoldGeometry(nodesGeom.extend(g,ng),false);
	  }
      }

    ng=node_down->getUpGeometry();
    if (ng!=nodesGeom_end())
      {
	pb->sort(node::compareItLess());pb->unique();
	for (std::list<node_it>::iterator nit = pb->begin();nit!=pb->end();nit++)
	  {
	    nodeGeom_it g=(*nit)->getUpGeometry();
	    (*nit)->setManifoldGeometry(nodesGeom.extend(g,ng),true);
	  }
      }      
       
    // and finally we can remove the nodes
    removeNode(endpts.first,new_it);
    removeNode(endpts.second,new_it);
    //mygettime();printf("\n");
    canceled = true;
    return new_it;
}

inline void NDcomplex::
simplify(std::vector<double> level, bool useRatio)
{
  int npass=0;
  std::vector<struct timeval> time;
  std::vector<struct timeval>::iterator time_it;
  struct timeval wc_time;
  double dt;
  long N=0;
  long Nold;
  long Ntot=0;
  long Nskipped=0;
  long Nimpossible=0;
  std::set<node_it,node::compareItLess> local_imp_list;

  gettimeofday(&wc_time, NULL);time.push_back(wc_time);    
 
  printf ("    Cancelling pairs with persistance%s < [%.2e",(useRatio)?" ratio":"",level[0]);
  for (int i=1;i<ndims;i++) printf(",%.2e",level[i]);
  printf("].\n");
  
  N=0;
  do {
    N=0;
    printf("    Pairing ... ");fflush(0);
    for (node_it n= nodes_begin(); n!=nodes_end();n++)
      n->setPNode(n);
    //do {
      
    Nold=N;
    N=0;
    for (node_it n1= nodes_begin(); n1!=nodes_end();n1++)
      {
	node_it n2;
	arc_it arc;
	char t=n1->getType();
	double ld,lu;

	if (t==0)
	  {ld=-1;lu=level[0];}
	else if (t==ndims_net)
	  {ld=level[ndims_net-1];lu=-1;}
	else {ld=level[t-1];lu=level[t];}
	  
	//if ((n1->getPNode() != n1)&&(n1->getPNode()->getPNode()==n1)) continue;
	  
	n2=n1->findWeakestSimpleNeighour(n1,arcs.end(),ld,lu,useRatio);
	n1->setPNode(n2);
      }
      
    for (node_it n1= nodes_begin(); n1!=nodes_end();n1++)
      {
	node_it n2=n1->getPNode();
	if (n2==n1) continue;
	if (n2->getPNode() != n1) 
	  {
	    /*//n2->setPNode(n2);*/
	    n1->setPNode(n1);
	    continue;
	  }
	N++;
      }    
    N=0;
    updtQueue<node_it,int,node::compareItLess> cancelQueue;
    //updtQueue<node_it,int,node::compareItMore> cancelQueue;
    npass++;
    

    for (node_it n1=nodes_begin(); n1!=nodes_end();n1++)
      {
	node_it n2=n1->getPNode();

	if (n1->getType()>=n2->getType()) continue;

	double p;
	
	if (useRatio) p=n1->persistenceR();
	else p=n1->persistence();
	
	if (p > level[n1->getType()]) {	 
	  continue;
	}
	

	cancelQueue.insert(n1,n1->pairCancellationOutcome());
      }
    
    if (cancelQueue.size())
      {
	printf ("cancelling ... ");fflush(0);
	do {
	  node_it n1=cancelQueue.pop();
	  node_it n2 = n1->getPNode();
	  bool loop;
	  arc_it arc = n1->findArc(n1,n2,arcs_end(),loop);
	  std::vector<node_it> node_updt;

	  if (arc==arcs_end())
	    {			
	      //Nskipped++;			
		continue;
	    }

	  bool canceled=true;
	  cancelPair(arc,canceled,back_inserter(node_updt));
	  if (canceled)
	    {	
	      node_it n3;
	      for (std::vector<node_it>::iterator up_it=node_updt.begin();up_it != node_updt.end();up_it++)
		{
		  n1 = *up_it;
		  if (n1->getPNode() == n1)
		    {	
		      n2=n1->findWeakestSimpleNeighour(n1,arcs.end(),true);
		      if (n2->getPNode()!=n2) continue;

		      n3=n2->findWeakestSimpleNeighour(n2,arcs.end(),true);
		      if (n3!=n1) continue;

		      n1->setPNode(n2);
		      n2->setPNode(n1);

		      double p;		      
		      
		      if (n1->getType()<n2->getType()) n3=n1;
		      else n3=n2;

		      if (useRatio) p=n3->persistenceR();
		      else p=n3->persistence();

		      if (p > level[n3->getType()]) continue;

		      //int val = n3->pairCancellationOutcome();
		      cancelQueue.insert(n3,n3->pairCancellationOutcome());
		      /*
		      if (cancelQueue.updateVal(n3,val)==cancelQueue.end())
			{
			  cancelQueue.insert(n3,val);
			  //Nskipped--;
			}	
		      */
		    }
		  else 
		    {
		      n2=n1->getPNode();
		      if (n1->getType()<n2->getType()) n3=n1;
		      else n3=n2;
		      
		      cancelQueue.updateVal(n3,n3->pairCancellationOutcome());
		    }
		}
	      
	      N++;
	    }
	  else
	    {
	      n1->setPNode(n1);
	      n2->setPNode(n2);
	      
	      //printf("PB\n");
	      local_imp_list.insert(n1);
	      Nimpossible++;
	    }
	} while (cancelQueue.size()!=0);
      }
    
    if (N) {printf("done. (%ld removed)\n",N);fflush(0);}
    //else {printf("no cancellable pair left.\n");fflush(0);}
    Ntot+=N;
    
    gettimeofday(&wc_time, NULL);time.push_back(wc_time);
    time_it=time.end();time_it--;
    dt = (double)(*time_it).tv_sec + ((double)(*time_it).tv_usec)*1.E-6;
    time_it--;
    dt -= (double)(*time_it).tv_sec + ((double)(*time_it).tv_usec)*1.E-6;

    } while (N!=0);
    
  gettimeofday(&wc_time, NULL);time.push_back(wc_time);
  time_it=time.end();time_it--;
  dt = (double)(*time_it).tv_sec + ((double)(*time_it).tv_usec)*1.E-6;
  time_it=time.begin();
  dt -= (double)(*time_it).tv_sec + ((double)(*time_it).tv_usec)*1.E-6;

  printf ("no cancellable pair left.\n    Cancellation took %.2fs (%ld canceled).\n",dt,Ntot);
	    
}

inline int NDcomplex::
cancelTaggedNodes_Raw(bool forceLoops, bool solveConflicts)
{
    int NCanceled =0;
    int NCanceled_old;
    bool canceled;
    int NSkipped;
    int NFailed=0;
    int NForced=0;
    int npass=1;
    struct timeval time1;
    struct timeval time2;

    printf ("    Cancelling pairs (raw - p%d) ... ",npass);fflush(0);
    gettimeofday(&time1, NULL);
    do {
	NCanceled_old=NCanceled;
	NSkipped=0;
	NFailed=0;

	printf ("\r    Cancelling pairs (raw - p%d) ... ",npass);fflush(0);

	for (node_it n1= nodes_begin(); n1!=nodes_end();)
	{
	    node_it new_node=n1;
	    new_node++;

	    if (n1->isTagged())
	    {
		
		node_it n2 = n1->getPNode();
		bool loop;
		arc_it arc = n1->findArc(n1,n2,arcs_end(),loop);
		
		if (new_node==n2) new_node++;

		if (arc==arcs_end())
		{			
		    NSkipped++;	
		    n1=new_node;
		    continue;
		}
		
		canceled=true;

		cancelPair(arc,canceled);
		if (!canceled)
		{
		    if (forceLoops)
		    {
			cancelPair(arc,canceled,true);
			NForced++;
			NCanceled++;
		    }
		    else NFailed++;
		}
		else NCanceled++;
	    }
	    n1=new_node;
	}
	npass++;
    } while (NCanceled!=NCanceled_old);
    
    printf ("done.\n");
    gettimeofday(&time2, NULL);

    double dt = (double)time2.tv_sec + ((double)time2.tv_usec)*1.E-6;
    dt -= (double)time1.tv_sec + ((double)time1.tv_usec)*1.E-6;

    printf ("    Cancellation took %.2fs (%d C, %d S, %d Fa, %d Fo).\n",dt,NCanceled,NSkipped,NFailed,NForced);
    return 0;
}

inline int NDcomplex::
cancelTaggedNodes_Smart(bool forceLoops, bool solveConflicts)
{
  std::vector<node_it> imp_list;
  std::vector<node_it> con_list;

  cancelTaggedNodes_Smart(back_inserter(imp_list),back_inserter(con_list),forceLoops,solveConflicts);
  return 1;
}

template <class OutputIterator>
inline int NDcomplex::
cancelTaggedNodes_Smart(OutputIterator imp_list,OutputIterator con_list,bool forceLoops, bool solveConflicts)
{
  int npass=0;
  int N,Ns,Nskipped;
  int Ntot=0, Nimp=0;
    
  bool canceled;
  std::vector<node_it> local_imp_list;
   
  std::vector<struct timeval> time;
  std::vector<struct timeval>::iterator time_it;
  struct timeval wc_time;
  double dt;
    
  gettimeofday(&wc_time, NULL);time.push_back(wc_time);    
 
  do {
	
    updtQueue<node_it,int,node::compareItLess> cancelQueue;
  	
    npass++;
	
    printf ("    Cancelling pairs (smart) ... ");fflush(0);
    local_imp_list.clear();
    N=Ns=Nskipped=0;
    int i=0;
    for (node_it n1= nodes_begin(); n1!=nodes_end();n1++)
      {
	if (n1->isTagged()) 
	  {
	    node_it n2 = n1->getPNode();

	    if (n1->findArc(n1,n2,arcs_end()) == arcs_end()) continue;
		
	    if (n1->getType()>n2->getType())
	      {
		cancelQueue.insert(n1,n1->pairCancellationOutcome());
	      }
	    else
	      {
		cancelQueue.insert(n2,n2->pairCancellationOutcome());
	      }			
	  }
      }
	
    //printf("%ld tags ",cancelQueue.size());fflush(0);

    if (cancelQueue.size()) do
			      {
				node_it n1=cancelQueue.pop();
				node_it n2 = n1->getPNode();
				//bool loop;
				std::vector<node_it> node_updt;
				bool loop;
				arc_it arc = n1->findArc(n1,n2,arcs_end(),loop);
	    
				if (arc==arcs_end())
				  {			
				    Nskipped++;		
				    continue;
				  } 

				canceled=true;
	    
				cancelPair(arc,canceled,back_inserter(node_updt),forceLoops);
	    
				if (canceled)
				  {
				    if (loop) Nimp++;

				    for (std::vector<node_it>::iterator up_it=node_updt.begin();up_it != node_updt.end();up_it++)
				      {
					node_it curn = (*up_it);
					if ((curn->isTagged())||(curn->getPNode()->isTagged()))
					  {		      
					    node_it na = curn;
					    node_it nb = na->getPNode();
					    node_it nc;

					    if (na->findArc(na,nb,arcs_end()) == arcs_end()) continue;

					    if (na->getType()>nb->getType()) nc=na;
					    else nc=nb;
		      
					    int val = nc->pairCancellationOutcome();

		      
					    if (cancelQueue.updateVal(nc,val)==cancelQueue.end())
					      {
						cancelQueue.insert(nc,val);
						//Nskipped--;
					      }	
		      
					  }

				      }
				    N++;
				  }
				else // impossible cancellation (loop)
				  {	
				    //if (!forceLoops) 
				    n1->clearTag(); // only the impossible nodes get untagged, not the conflicting ones.
				    local_imp_list.push_back(n1);
				    //if (!forceLoops) *imp_list = n1;
	
				    Ns++;
				    Nimp++;
				  }
	    
			      } while (cancelQueue.size()!=0);
    //Ntot += N;

    // Impossible cancellations
	
    // THIS IS NEVER TRUE == COMMENTED
    if ((false)&&(forceLoops)&&(local_imp_list.size()!=0))
      {
	int Nrem=0;
	int Ntorem = local_imp_list.size();
	int markVal= 0;
	printf ("(forcing %d loops ",Ntorem);fflush(0);
	for (std::vector<node_it>::iterator vit= local_imp_list.begin();
	     vit!=local_imp_list.end();vit++)
	  {
	    struct timeval t1;
	    struct timeval t2;
	    struct timeval t3;
	    double dt;

	    gettimeofday(&t1, NULL);

	    node_it n1 = *vit;
	    node_it n2 = n1->getPNode();
	    //printf ("\nPair : %d/%d\n      %s / %s\n",Nrem,Ntorem,n1->getInfo(true).c_str(),n2->getInfo(true).c_str());

	    bool loop;
	    arc_it arc = n1->findArc(n1,n2,arcs_end(),loop);

	    gettimeofday(&t2, NULL);
	    dt = (double)t2.tv_sec + ((double)t2.tv_usec)*1.E-6;
	    dt -= (double)t1.tv_sec + ((double)t1.tv_usec)*1.E-6;
	    //printf (" (+%.6fs)\n     arc %s\n",dt,arc->getInfo(true).c_str());
	    canceled=true;
		
	    cancelPair(arc,canceled,true);
	    gettimeofday(&t3, NULL);
	    dt = (double)t3.tv_sec + ((double)t3.tv_usec)*1.E-6;
	    dt -= (double)t2.tv_sec + ((double)t2.tv_usec)*1.E-6;
	    //printf(" (+%.6fs)\n     Cancelled.\n",dt);
	    if (Nrem==markVal) {
	      printf(".");fflush(0);
	      markVal += (int)(Ntorem/10.);
	    }
	    N++;Nrem++;
	  }
	local_imp_list.clear();
	printf (") done. (%d removed)\n",N);
	//computePersistencePairs();
	if (N) printf ("  ");
      }
    else if (N) printf ("(%d rem.)\n",N);fflush(0);
    Ntot += N;
	
    gettimeofday(&wc_time, NULL);time.push_back(wc_time);
    time_it=time.end();time_it--;
    dt = (double)(*time_it).tv_sec + ((double)(*time_it).tv_usec)*1.E-6;
    time_it--;
    dt -= (double)(*time_it).tv_sec + ((double)(*time_it).tv_usec)*1.E-6;
    
    if (!forceLoops) N=0;
  }
  while (N!=0);

  gettimeofday(&wc_time, NULL);time.push_back(wc_time);
  time_it=time.end();time_it--;
  dt = (double)(*time_it).tv_sec + ((double)(*time_it).tv_usec)*1.E-6;
  time_it=time.begin();
  dt -= (double)(*time_it).tv_sec + ((double)(*time_it).tv_usec)*1.E-6;
  //printf("Nskipped = %d\n",Nskipped);
  //Nskipped=0;
    
  char forcedtxt[255];
  char solvedtxt[255];

  if (forceLoops) strcpy(forcedtxt,"forced ");
  else strcpy(forcedtxt,"");
    
  if (solveConflicts) strcpy(solvedtxt,"solved ");
  else strcpy(solvedtxt,"");

  int nImpConf=0;
  int NskippedInit=Nskipped;
  // conflicts resolution
  if (Nskipped) 
    {
      if (solveConflicts) 
	{
	  fprintf(stderr,"WARNING in SolveConflicts:  option is NOT available !!!!");
	  solveConflicts=0;
	}
	
      int nsolved;
      do {
	nImpConf=0;
	nsolved=0;
	
	//if (solveConflicts) {printf ("\n       solving %d conflicts ... ",Nskipped);fflush(0);}

	updtQueue<node_it,int,node::compareItLess> cancelQueue;
	for (node_it n1= nodes_begin(); n1!=nodes_end();n1++)
	  {
	    if (n1->isTagged())
	      {}
	  }
	int nfailedConf=0;
	} while ((Nskipped)&&(nsolved!=0));
      }
    
    Nskipped=0;
    for (node_it n1= nodes_begin(); n1!=nodes_end();n1++)
      if (n1->isTagged()) {*con_list=n1;Nskipped++;}

    for (std::vector<node_it>::iterator vit= local_imp_list.begin();
	 vit!=local_imp_list.end();vit++)
      *imp_list = *vit;    

    char nsolvedtxt[255];
    if (solveConflicts)
      {
	if (forceLoops)
	  sprintf(nsolvedtxt,"%d(%d f.)",NskippedInit,nImpConf);
	else
	  sprintf(nsolvedtxt,"%d/%d",NskippedInit-nImpConf,NskippedInit);
      }
    else sprintf(nsolvedtxt,"%d",Nskipped);
    

    printf ("done.\n    Cancellation took %.2fs (%d canceled, %s %sconflicts, %d %sloops).\n",
	    dt,Ntot,nsolvedtxt,solvedtxt,Nimp,forcedtxt);
    //defragNodesGeometry();
    return Ntot;  
}


inline void NDcomplex::
toAscii(const char *fname)
{
    FILE *f=NULL;
    
    printf ("Saving complex to ASCII ... ");fflush(0);
    f=fopen(fname,"w");
    for (NDcomplex::arc_it it = arcs.begin();it !=arcs.end();it++)
    {
	std::pair<NDcomplex::node_it,NDcomplex::node_it> node = it->endpoints();
	fprintf (f,"%d %ld %d %ld\n",
		 node.first->getType(),node.first->getId(),
		 node.second->getType(),node.second->getId());
    }
    fclose(f);
    printf ("done.\n");
}

/*
inline void NDcomplex::
cyclesToAscii(const char *fname)
{
    return ;
	printf ("Saving cycles to ASCII ... ");fflush(0);
	FILE *f=fopen(fname,"w");
	
	std::vector<node_it> nei;
	std::vector<int> coord;
	std::vector<int>::iterator coord_it;
	std::vector<node_it>::iterator nei_it;

	fprintf (f,"%ld\n",fullCycles.size());
	for (cycle_map_it it= fullCycles.begin();it!=fullCycles.end();it++)
	{
		coord.clear();
		for (cycle_it cycle= it->second.begin();cycle!=it->second.end();cycle++)
		{
			nei.clear();
			(*cycle)->getNeighbors((*cycle),arcs_end(),back_inserter(nei),(*cycle)->getType()+1);
			for (nei_it = nei.begin();nei_it != nei.end();nei_it++)
			{
				coord.push_back((*cycle)->getId());
				coord.push_back((*cycle)->getType());
				coord.push_back((*nei_it)->getId());
				coord.push_back((*nei_it)->getType());
				//coord.push_back((*cycle)->getVal());
				//coord.push_back((*nei_it)->getVal());
			}
		}
		fprintf (f,"%ld %d %ld %d %ld\n",
				it->first->getId(),it->first->getType(),
				it->first->getPNode()->getId(),it->first->getPNode()->getType(),
				coord.size()/4);//,it->first->getVal(),it->first->getPNode()->getVal());
		for (int i=0;i<coord.size();i+=4)
		{
			fprintf (f,"%d %d %d %d\n",
					coord[i],coord[i+1],
					coord[i+2],coord[i+3]);//,coord[i+4],coord[i+5]);

		}
	}

	fclose(f);

	

	char fname2[255];
	sprintf(fname2,"canonic_%s",fname);
	f=fopen(fname2,"w");
	fprintf (f,"%ld\n",fullCycles.size());
	for (cycle_map_it it= fullCycles.begin();it!=fullCycles.end();it++)
	{
		coord.clear();
		cycle ccycle(it->second);
		//if (it->first==nodes_end()) continue;
		//if (it->first->isInfinite()) continue;
		//if (it->first->getPNode() == nodes_end()) continue;

		//printf ("start ...");fflush(0);

		//if ((it->first!=nodes_end())&&(it->first->getPNode() != nodes_end()))
		//printf ("id/typ : %d/%d %d\n",it->first->getId(),it->first->getType(),it->first->isPositive());
		if (it->first->getType() == 1)
		    makeCycleCanonic(it->first,ccycle);//,it->first->persistence());

		//printf ("stop\n");

		for (cycle_it cycle= ccycle.begin();cycle!=ccycle.end();cycle++)
		{
			nei.clear();
			(*cycle)->getNeighbors((*cycle),arcs_end(),back_inserter(nei),(*cycle)->getType()+1);
			for (nei_it = nei.begin();nei_it != nei.end();nei_it++)
			{
				coord.push_back((*cycle)->getId());
				coord.push_back((*cycle)->getType());
				coord.push_back((*nei_it)->getId());
				coord.push_back((*nei_it)->getType());
				//coord.push_back((*cycle)->getVal());
				//coord.push_back((*nei_it)->getVal());
			}
		}
		fprintf (f,"%ld %d %ld %d %ld\n",
				it->first->getId(),it->first->getType(),
				it->first->getPNode()->getId(),it->first->getPNode()->getType(),
				coord.size()/4);//,it->first->getVal(),it->first->getPNode()->getVal());
		for (int i=0;i<coord.size();i+=4)
		{
			fprintf (f,"%d %d %d %d\n",
					coord[i],coord[i+1],
					coord[i+2],coord[i+3]);//,coord[i+4],coord[i+5]);

		}
	}

	fclose(f);

	printf ("done.\n");
}
*/


void NDcomplex::write(std::ofstream &str)
{
  std::map<node *, long> node2id;
  std::map<arc *, long> arc2id;
  node_it curnode;
  arc_it curarc;
  long i;
  arcGeom_it ag;
  nodeGeom_it ng;
  char hpp = (char)havePPairs;
  int hc=0;

  for (i=0;i<haveCycles.size();i++)
    if (haveCycles[i]) hc |= (1<<i);

  str.write((const char *)&ndims,sizeof(int));  
  str.write((const char *)&ndims_net,sizeof(int));  
  str.write((const char *)&hpp,sizeof(char));  
  str.write((const char *)&hc,sizeof(int));  

  i=arcs.size();
  str.write((const char *) &i,sizeof(long));
  i=nodes.size();
  str.write((const char *) &i,sizeof(long));

  for (curnode = nodes.begin(),i=0;curnode!=nodes.end();curnode++,i++)
    node2id.insert(std::make_pair(&(*curnode),i));
  for (curarc = arcs.begin(),i=0;curarc!=arcs.end();curarc++,i++)
    arc2id.insert(std::make_pair(&(*curarc),i));

  for (curnode = nodes.begin(),i=0;curnode!=nodes.end();curnode++,i++)
    {
      curnode->write(str,node2id,arc2id,arcs.end(),nodes.end());      
    }

  for (curarc = arcs.begin(),i=0;curarc!=arcs.end();curarc++,i++)
    {
      curarc->write(str,node2id,arc2id,arcs.end());
    }

  arc2id.clear();
  node2id.clear();
  
  std::map<node::geom_el *,long> nGeom2id;
  nodesGeom.write(str,nGeom2id);
  for (curnode = nodes.begin();curnode!=nodes.end();curnode++)
    {
      ng=curnode->getUpGeometry();
      if (ng!=nodesGeom.end()) i=nGeom2id[&(*ng)];
      else i=-1;
      str.write((const char *) &i, sizeof(long));

      ng=curnode->getDownGeometry();
      if (ng!=nodesGeom.end()) i=nGeom2id[&(*ng)];
      else i=-1;
      str.write((const char *) &i, sizeof(long));
    }
  nGeom2id.clear();
 
  std::map<arc::geom_el *,long> aGeom2id;
  arcsGeom.write(str,aGeom2id);
  for (curarc = arcs.begin();curarc!=arcs.end();curarc++)
    {
      ag=curarc->getGeometry();
      if (ag!=arcsGeom.end()) i=aGeom2id[&(*ag)];
      else i=-1;
      str.write((const char *) &i, sizeof(long));
    }
  aGeom2id.clear();

  int md=(int)currentMode;
  str.write((const char *) &md, sizeof(int));

}

void NDcomplex::read(std::ifstream &str)
{
  std::vector<node_it> id2node;
  std::vector<arc_it> id2arc;
  node_it curnode;
  arc_it curarc;
  long i,j;
  char hpp;
  int hc;

  str.read((char *)&ndims,sizeof(int));
  str.read((char *)&ndims_net,sizeof(int));
  str.read((char *)&hpp,sizeof(char));
  havePPairs=hpp;
  str.read((char *)&hc,sizeof(int));
  for (i=0;i<=ndims;i++)
    {
      if (hc&(1<<i)) haveCycles.push_back(true);
      else haveCycles.push_back(false);
    }
  
  
  str.read((char *) &i,sizeof(long));
  str.read((char *) &j,sizeof(long));  
  
  arcs.resize(i);
  nodes.resize(j);
  id2arc.resize(i);  
  id2node.resize(j);
  
  for (curnode = nodes.begin(),i=0;curnode!=nodes.end();curnode++,i++)
    id2node[i]=curnode;
  
  for (curarc = arcs.begin(),i=0;curarc!=arcs.end();curarc++,i++)
    id2arc[i]=curarc;
 
  for (curnode = nodes.begin(),i=0;curnode!=nodes.end();curnode++,i++)
    curnode->read(str,id2node,id2arc,arcs.end(),nodes.end());
  
  for (curarc = arcs.begin(),i=0;curarc!=arcs.end();curarc++,i++)
    curarc->read(str,id2node,id2arc,arcs.end());

  id2node.clear();
  id2arc.clear();

  std::vector<node::geom_it> nId2geom;
  nodesGeom.read(str,nId2geom);
  for (curnode = nodes.begin();curnode!=nodes.end();curnode++)
    {
      str.read((char *) &i,sizeof(long));
      if (i!=-1) curnode->setManifoldGeometry(nId2geom[i],true);
      else curnode->setManifoldGeometry(nodesGeom.end(),true);

      str.read((char *) &i,sizeof(long));
      if (i!=-1) curnode->setManifoldGeometry(nId2geom[i],false);
      else curnode->setManifoldGeometry(nodesGeom.end(),false);
    }
  nId2geom.clear();

  std::vector<arc::geom_it> aId2geom;
  arcsGeom.read(str,aId2geom);
  for (curarc = arcs.begin();curarc!=arcs.end();curarc++)
    {
      str.read((char *) &i,sizeof(long));
      if (i!=-1) curarc->setManifoldGeometry(aId2geom[i]);
      else curarc->setManifoldGeometry(arcsGeom.end());
    }
  aId2geom.clear();

  int md=(int)currentMode;
  str.read((char *) &md, sizeof(int));
  //currentMode = (NDcomplex::cancellationMode)md;
}


#endif
