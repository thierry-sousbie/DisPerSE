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
#ifndef __ARCS_AND_NODES_HXX__
#define __ARCS_AND_NODES_HXX__

#include <vector>
#include <set>
#include <list>
#include <cmath>
#include <functional>
#include <stack>
#include <algorithm>

#include <float.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include "manifolds.hxx"
#include "cells.hxx"

#define NODE_FLAG_BOUNDARY (1<<0)
#define NODE_FLAG_OUT (1<<1)
#define NODE_FLAG_INFINITE (1<<2)
#define NODE_FLAG_TYPEMASK ((1<<3)-1)
#define NODE_FLAG_PAIRED (1<<3)
#define NODE_FLAG_POSITIVE (1<<4)
#define NODE_FLAG_NEGATIVE (1<<5)

#define NODE_FLAG_TAG (1<<6)
//#define NODE_FLAG_TAG2 (1<<7)
#define NODE_FLAG_DELETED (1<<7)

class NDcomplex;
class NDcomplex_node;
class NDcomplex_arc;

typedef std::list<NDcomplex_arc> NDcomplex_arcs_list;
typedef std::list<NDcomplex_arc>::iterator NDcomplex_arcs_list_it;
typedef std::list<NDcomplex_arc>::reverse_iterator NDcomplex_arcs_list_rit;
typedef std::list<NDcomplex_node> NDcomplex_nodes_list;
typedef std::list<NDcomplex_node>::iterator NDcomplex_nodes_list_it;
typedef std::list<NDcomplex_node>::reverse_iterator NDcomplex_nodes_list_rit;

typedef cellIdentity<> NDcomplex_cellType;

class NDcomplex_arc {
public:
  static const double infinity;
  
  typedef NDcomplex_arc arc;
  typedef NDcomplex_node node;
  
  typedef NDcomplex_arcs_list arcs_list;
  typedef NDcomplex_arcs_list_it arc_it;
  typedef NDcomplex_arcs_list_rit arc_rit;
  typedef NDcomplex_nodes_list nodes_list;
  typedef NDcomplex_nodes_list_it node_it;
  typedef NDcomplex_nodes_list_rit node_rit;
  
  typedef NDcomplex_cellType cellType ;
  typedef arcGeometry< cellType > geom;
  typedef geom::geom_it geom_it;
  typedef geom::geometry geom_el;
 

private:
  friend class NDcomplex_node;
  friend class NDcomplex;

  geom_it geometry;

  node_it n[2];
  arc_it next[2];
  arc_it prev[2];

  inline arc_it &getNextRef(const node_it node);
  inline arc_it &getPrevRef(const node_it node);
  inline void insertBeforeFirst(arc_it &self, arc_it &a, node_it &mynode);
  inline void unLink(arc_it &self, const arc_it &null_arc);

public:
  inline std::string getInfo(bool all=false) const; 
  
  //inline double persistence() const;
  //inline double persistenceR() const;  
  inline double delta() const;
  inline double ratio() const;  
  //inline double robustness() const;
  //inline double robustnessR() const;  

  inline std::pair<node_it,node_it> getReferenceNodes() const;
  inline std::pair<node_it,node_it> endpoints() const;
  inline arc_it getNext(const node_it node) const;
  inline arc_it getPrev(const node_it node) const;
  inline int getNextId(const node_it node) const;
  inline int getPrevId(const node_it node) const;
  inline node_it getOtherNode(const node_it node) const;
  inline bool isUnique(arc_it null_arc) const;
  //inline arc_it cancelingNeighbor(arc_it self,arc_it null_arc) const;
  inline geom_it const getGeometry() {return geometry;}
  inline geom_it deleteGeometry(const geom_it null_geom)
  {geom_it it = geometry;if (geometry!=null_geom) geometry=null_geom; return it;}
  
  void setManifoldGeometry(geom_it it)
  {
    geometry = it;
  }

  NDcomplex_arc() {};

  NDcomplex_arc(node_it n0, node_it n1, const arc_it &null_arc, const geom_it null_geom)
  {
    n[0]=n0;
    n[1]=n1;
    prev[0]=prev[1]=next[0]=next[1]=null_arc;
    geometry = null_geom;
  }

  ~NDcomplex_arc()
  {
  }

  void write(std::ofstream &str,std::map<NDcomplex_node *, long> &node2id,std::map<NDcomplex_arc *, long> &arc2id, const arc_it &null_arc);
  void read(std::ifstream &str, std::vector<node_it> &id2node, std::vector<arc_it> &id2arc,const arc_it &null_arc);
  /*
    void write(std::ofstream &str) {}
  
  */
};


class NDcomplex_node {
public:
  static const double infinity;
  
  typedef NDcomplex_arc arc;
  typedef NDcomplex_node node;
  
  typedef NDcomplex_arcs_list arcs_list;
  typedef NDcomplex_arcs_list_it arc_it;
  typedef NDcomplex_arcs_list_rit arc_rit;
  typedef NDcomplex_nodes_list nodes_list;
  typedef NDcomplex_nodes_list_it node_it;
  typedef NDcomplex_nodes_list_rit node_rit;
  
  typedef NDcomplex_cellType cellType;
  typedef nodeGeometry< cellType > geom;
  typedef geom::geom_it geom_it;
  typedef geom::geometry geom_el;

  typedef cellType::idT idT;
  typedef cellType::typeT typeT;

private:
  friend class NDcomplex_arc;
  friend class NDcomplex;
  
  geom_it up;
  geom_it down;
  
  node_it ref_up;
  node_it ref_down;

  node_it p_pair;
  char type;
  char flags;
  cellType cellID;
  double val;
  double val2; 
  int numArcs;
  int numUpArcs;
  
  arc_it arc_list;

  inline void set(cellType cell, char typep, double valp=0,double valp2=0);

public:
  inline void setDeleted() {flags|=NODE_FLAG_DELETED;}
  inline bool isDeleted() const {return flags&NODE_FLAG_DELETED;}
  inline void clearDeleted() {flags&=(~NODE_FLAG_DELETED);}

  inline void setTag() {flags|=NODE_FLAG_TAG;}
  inline bool isTagged() const {return flags&NODE_FLAG_TAG;}
  inline void clearTag() {flags&=(~NODE_FLAG_TAG);}
  inline void unlink(arc_it null_arc) {numArcs=0;numUpArcs=0;arc_list=null_arc;}
    
  inline std::string getInfo(bool all=false) const; 
  inline void setType(char typep) {type=typep;}
   
  inline void setPNode(node_it node) {p_pair = node;}  
  inline node_it getPNode() const {return p_pair;}    

  inline void setRefUp(node_it node) {ref_up = node;}  
  inline void setRefDown(node_it node) {ref_down = node;} 
  inline node_it getRefUp() {return ref_up;}  
  inline node_it getRefDown() {return ref_down;} 
  

  inline int pairCancellationOutcome() const;
  inline int nArcs() const {return numArcs;}
  inline int nUpArcs() const {return numUpArcs;}
  inline char getFlags() const {return flags;}
 
  inline bool isInfinite() const {return flags&NODE_FLAG_INFINITE;}
  inline bool isOut() const {return flags&NODE_FLAG_OUT;}
  inline bool isBoundary() const {return flags&NODE_FLAG_BOUNDARY;}
  inline bool isPositive() const {return flags&NODE_FLAG_POSITIVE;}
  inline bool isNegative() const {return flags&NODE_FLAG_NEGATIVE;}
  inline double uniqueID() {return (isOut()?-1.:1.)*((double)cellID.id() +((double)cellID.type())*0.1);}
  inline idT getId() const {return cellID.id();}
  inline cellType getCell() const {return cellID;}
  inline double getVal() const {return isInfinite()?(-infinity):val;}
  inline double getVal2() const {return isInfinite()?(-infinity):val2;}
  inline typeT getType() const {return type;}

  inline arc_it getArc() const {return arc_list;}
  inline geom_it const getUpGeometry() {return up;}
  inline geom_it const getDownGeometry() {return down;}
  inline geom_it deleteUpGeometry(const geom_it null_geom)
  {geom_it it = up;if (up!=null_geom) up=null_geom; return it;}
  inline geom_it deleteDownGeometry(const geom_it null_geom)
  {geom_it it = down;if (down!=null_geom) down=null_geom; return it;}

  inline char flagsAnd(char f) {return flags&=f;}
  inline char flagsOr(char f) {return flags|=f;}

  inline bool check(node_it self, arc_it null_arc);
  template <class OutputIterator>
  inline int getNeighbors(node_it self, arc_it null_arc, OutputIterator it, char type=-1);
	
  inline node_it findWeakestSimpleNeighour(node_it self, arc_it null_arc, double downThreshold=-1, double upThreshold=-1, bool useRatio=false);

  template <class excludeClass>
  inline arc_it lowestDeltaArc(node_it self, arc_it null_arc,excludeClass exclude, char type=-1);

  arc_it findMultiArc(node_it self, node_it other, arc_it null_arc, int &nc,int nmax=-1);
  arc_it findArc(node_it self, node_it other, arc_it null_arc, bool &loop);
  arc_it findArc(node_it self, node_it other, arc_it null_arc);

  inline int countArcs(node_it self, arc_it null_arc, int type =-1) const;

  inline void insertArc(arc_it &it, node_it &self,const arc_it &null_arc);
  inline bool operator<(const node &other) const;

  void setManifoldGeometry(geom_it it, bool ascending)
  {
    if (ascending) up = it;
    else down=it;
  }

  NDcomplex_node() {};

  NDcomplex_node(const NDcomplex_node &n,const NDcomplex_arc::arc_it &null_arc,const NDcomplex_node::node_it &null_node, const geom_it null_geom)
  {
    this->cellID = n.cellID;
    this->type = n.type;
    this->val = n.val;
    this->val2 = n.val2; 
    this->flags=n.flags;

    arc_list=null_arc;
    p_pair=null_node;    
    ref_up=ref_down=null_node;
    up=down=null_geom;
    numArcs=0;
    numUpArcs=0;
  }

  NDcomplex_node(cellType cell, char typep, double valp, const NDcomplex_arc::arc_it &null_arc,const NDcomplex_node::node_it &null_node, const geom_it null_geom,double valp2 = 0, char flags = 0)
  {
    set(cell,typep,valp,valp2);
    this->flags=flags;

    arc_list=null_arc;
    p_pair = null_node;
    ref_up=ref_down=null_node;
    up=down=null_geom;
    numArcs=0;
    numUpArcs=0;	  
  }

  ~NDcomplex_node()
  {}

  static bool exchangePPairs(node_it node1, node_it node2) {
    node_it n1,n2;

    if (node1->getType() == node2->getType())
      {
	n1=node1;
	n2=node2;
      }
    else if (node1->getType() == node2->getPNode()->getType())
      {
	n1=node1;
	n2=node2->getPNode();
      }
    else 
      {
	fprintf (stderr,"ERROR in NDComplex_node::exchangePPairs : Cannot exchange pairs with different types\n");
	return false;
      }
	
    node_it tmp=n1->getPNode();

    n1->setPNode(n2->getPNode());
    n2->setPNode(tmp);
	
    n1->getPNode()->setPNode(n1);
    n2->getPNode()->setPNode(n2);

    return true;
  }

  static double persistence(const node &n1,const node &n2) {
    double r;
    
    if ((n1.getId()==n2.getId())&&(n1.getType()==n2.getType())) return infinity;
    
    if ((n1.flags&NODE_FLAG_INFINITE)||(n2.flags&NODE_FLAG_INFINITE))
      {
	if ((n1.flags&NODE_FLAG_INFINITE)&&(n2.flags&NODE_FLAG_INFINITE))
	  return 0;
	else return infinity;
      }
    
    if (n1.getType()<n2.getType())
      r=n2.getVal()-n1.getVal();
    else
      r=n1.getVal()-n2.getVal();
    
    return (r<0)?-1:r;
  }
  
  static double persistence(const node_it n1,const node_it n2) {
    return persistence(*n1,*n2);
  }
  
  double persistence(const node_it n) const {
    return persistence(*this,*n);
  }    
  
  double persistence() const {
    return persistence(*this, *(this->p_pair));
  } 
  
  static double persistenceR(const node &n1,const node &n2) {
    if ((n1.getId()==n2.getId())&&(n1.getType()==n2.getType())) return infinity;
    //if (n1.cellID == n2.cellID) return infinity;

    if ((n1.flags&NODE_FLAG_INFINITE)||(n2.flags&NODE_FLAG_INFINITE))
      {
	if ((n1.flags&NODE_FLAG_INFINITE)&&(n2.flags&NODE_FLAG_INFINITE))
	  return 1;
	else return infinity;
      }
    if ((n1.getVal()-n2.getVal())*(n1.getType()-n2.getType())<0) return 0;

    if (n1.getVal()>n2.getVal())
      return n1.getVal()/n2.getVal();
    else
      return n2.getVal()/n1.getVal();
  }

  static double persistenceR(const node_it n1,const node_it n2) {
    return persistenceR(*n1,*n2);
  }

  double persistenceR(const node_it n) const {
    return persistenceR(*this,*n);
  }

  double persistenceR() const {
    return persistenceR(*this, *(this->p_pair));
  } 
  
  static bool CompareEqual(const node &n1,const node &n2) {
    if ((n1.cellID == n2.cellID)&&
	(n1.flags==n2.flags)&&
	(n1.type==n2.type)&&
	(n1.val==n2.val)&&
	(n1.val2==n2.val2)) return true;
    return false;
  }

  static bool compareLess(const node &n1,const node &n2) {
    
    if (CompareEqual(n1,n2)) return false;
    
    if ((n1.isInfinite())&&(!n2.isInfinite())) return false;
    if ((n2.isInfinite())&&(!n1.isInfinite())) return true;
    
    if (n1.getVal()!=n2.getVal()) return n1.getVal()<n2.getVal();
    else if (n1.getType()!=n2.getType()) return n1.getType()<n2.getType();
    else if (n1.getVal2()!=n2.getVal2()) return n1.getVal2()<n2.getVal2();
    else if (n1.getId()!=n2.getId()) return n1.getId()<n2.getId();
    else return (n1.getFlags()&NODE_FLAG_TYPEMASK)<(n2.getFlags()&NODE_FLAG_TYPEMASK);
  }

  static bool compareMore(const node &n1,const node &n2) {
    //if (NDcomplex_node::CompareEqual(n1,n2)) return false;
    return NDcomplex_node::compareLess(n2,n1);
  }

  class compareItLess : public std::binary_function<node_it,node_it,bool> {
  public :
    bool operator()(const node_it& n1, const node_it& n2) const
    {
      return NDcomplex_node::compareLess(*n1,*n2);
    }
  };

  class compareItPersistenceLess : public std::binary_function<node_it,node_it,bool> {
  public :
    bool operator()(const node_it& n1, const node_it& n2) const
    {
      return comparePersistenceLess(*n1,*n2);
    }
  };

  class compareItMore : public std::binary_function<node_it,node_it,bool> {
  public :
    bool operator()(const node_it& n1, const node_it& n2) const
    {
      return NDcomplex_node::compareLess(*n2,*n1);
    }
  };
    
  static bool comparePersistenceLess(const node &n1,const node &n2) {
    if (n1.getType()!=n2.getType())
      {
	return (n1.getType()>n2.getType());
      }
    double p1=n1.persistence();
    double p2=n2.persistence();
    if (p1==p2) return NDcomplex_node::compareLess(n1,n2);
    else return p1<p2;
  }
  
  static node_it lowestPersistencePair(const node_it &n,const node_it &n1,const node_it &n2) {
    double p1=persistence(n,n1);
    double p2=persistence(n,n2);
	
    if (p1==p2) {
      if (NDcomplex_node::compareLess(*n1,*n2)) return n1;
      else return n2;
    } else {
      if (p1<p2) return n1;
      else return n2;
    }
	
  }

  void write(std::ofstream &str,std::map<NDcomplex_node *, long> &node2id,std::map<NDcomplex_arc *, long> &arc2id, const arc_it &null_arc, const node_it &null_node);
  void read(std::ifstream &str, std::vector<node_it> &id2node, std::vector<arc_it> &id2arc,const arc_it &null_arc, const node_it &null_node);

  /*
    static bool comparePersistenceMore(const node &n1,const node &n2) {
    if (n1.getType()!=n2.getType())
    {
    if (n1.getType()==0) return true;
    return (n1.getType()>n2.getType());
    }
    double p1=n1.persistence();
    double p2=n2.persistence();
    if (p1==p2) return NDcomplex_node::compareMore(n1,n2);
    else return p1>p2;
    }
  */
};



/*   ARCS IMPLEMENTATION */

inline std::string NDcomplex_arc::getInfo(bool all) const
{
  std::ostringstream stm;
  stm << "A("<<n[0]->getInfo(all)<<","<<n[1]->getInfo(all);
  if (all) stm << ","<< std::scientific <<"Delta="<< delta() << "/Ratio="<<ratio();
 
  stm<<")";
  return stm.str() ;
}


inline NDcomplex_arc::arc_it &NDcomplex_arc::getNextRef(const NDcomplex_node::node_it node)
{
  assert((n[0]==node)||(n[1]==node));
  if (n[0]==node)
    return next[0];
  else
    return next[1];
}

inline NDcomplex_arc::arc_it &NDcomplex_arc::getPrevRef(const NDcomplex_node::node_it node)
{
  assert((n[0]==node)||(n[1]==node));
  if (n[0]==node)
    return prev[0];
  else
    return prev[1];
}


inline void NDcomplex_arc::insertBeforeFirst(NDcomplex_arc::arc_it &self, NDcomplex_arc::arc_it &a, NDcomplex_node::node_it &mynode)
{
  assert((a->n[0] == mynode)||(a->n[1] == mynode));
  if (a->n[0] == mynode)
    a->next[0] = self;
  else
    a->next[1] = self;

  if (n[0]==mynode)
    prev[0]=a;
  else
    prev[1]=a;
}

inline void NDcomplex_arc::unLink(NDcomplex_arc::arc_it &self, const NDcomplex_arc::arc_it &null_arc)
{
  arc_it cur_arc;
  node_it cur_node;

  cur_node = n[0];
  cur_arc = getPrev(cur_node);
  if (cur_arc == null_arc)
    cur_node->arc_list = getNext(cur_node);
  else
    cur_arc->getNextRef(cur_node) = getNext(cur_node);

  cur_arc = getNext(cur_node);
  if (cur_arc != null_arc)
    cur_arc->getPrevRef(cur_node) = getPrev(cur_node);

  cur_node = n[1];
  cur_arc = getPrev(cur_node);
  if (cur_arc == null_arc)
    cur_node->arc_list = getNext(cur_node);
  else
    cur_arc->getNextRef(cur_node) = getNext(cur_node);

  cur_arc = getNext(cur_node);
  if (cur_arc != null_arc)
    cur_arc->getPrevRef(cur_node) = getPrev(cur_node);

  n[0]->numArcs--;
  n[1]->numArcs--;  

  if (n[1]->getType()>n[0]->getType()) n[0]->numUpArcs--;
  else n[1]->numUpArcs--;

  next[0]=next[1]=prev[0]=prev[1]=null_arc;

}

/***** ARCS PERSISTENCE *****/
inline double NDcomplex_arc::delta() const {
  return NDcomplex_node::persistence(n[0],n[1]);
}

inline double NDcomplex_arc::ratio() const {
  return NDcomplex_node::persistenceR(n[0],n[1]);
}
/*
  inline double NDcomplex_arc::robustness() const {
  return std::min(n[0]->persistence(),n[1]->persistence());
  }

  inline double NDcomplex_arc::robustnessR() const {
  return std::min(n[0]->persistenceR(),n[1]->persistenceR());
  }
*/

inline std::pair<NDcomplex_node::node_it,NDcomplex_node::node_it> NDcomplex_arc::getReferenceNodes() const {
  node_it p0,p1;

  if ((n[0]->getType() - n[0]->getPNode()->getType()) == (n[1]->getType()-n[0]->getType()))
    p0 = n[0]->getPNode();
  else p0=n[0];

  if ((n[1]->getType() - n[1]->getPNode()->getType()) == (n[0]->getType()-n[1]->getType()))
    p1 = n[1]->getPNode();
  else p1=n[1];

  return std::make_pair(p0,p1);
}

inline std::pair<NDcomplex_node::node_it,NDcomplex_node::node_it> NDcomplex_arc::endpoints() const
{
  return std::pair<node_it,node_it>(n[0],n[1]);
}

inline NDcomplex_arc::arc_it NDcomplex_arc::getNext(const NDcomplex_node::node_it node) const
{
  assert((n[0]==node)||(n[1]==node));
  if (n[0]==node)
    return next[0];
  else
    return next[1];
}

inline NDcomplex_arc::arc_it NDcomplex_arc::getPrev(const NDcomplex_node::node_it node) const
{
  assert((n[0]==node)||(n[1]==node));
  if (n[0]==node)
    return prev[0];
  else
    return prev[1];
}

inline int NDcomplex_arc::getNextId(const NDcomplex_node::node_it node) const
{
  assert((n[0]==node)||(n[1]==node));
  if (n[0]==node)
    return 0;
  else
    return 1;
}

inline int NDcomplex_arc::getPrevId(const NDcomplex_node::node_it node) const
{
  assert((n[0]==node)||(n[1]==node));
  if (n[0]==node)
    return 0;
  else
    return 1;
}

inline NDcomplex_node::node_it NDcomplex_arc::getOtherNode(const NDcomplex_node::node_it node) const
{
  assert((n[0]==node)||(n[1]==node));
  if (n[0]==node)
    return n[1];
  else
    return n[0];
}

inline bool NDcomplex_arc::isUnique(NDcomplex_arc::arc_it null_arc) const
{
  node_it node;
  node_it other;
  arc_it a;
  int N=0;

  if (n[0]->nArcs()<=n[1]->nArcs()) node=n[0];
  else node=n[1];

  a=node->arc_list;
  other = this->getOtherNode(node);

  while (a!=null_arc)
    {

      if (a->getOtherNode(node)==other)
	{
	  N++;
	  if (N>1) return false;
	}
      a = a->getNext(node);
    }

  if (N==1) return true;
  else return false;
}


/*   NODES IMPLEMENTATION */

inline int NDcomplex_node::pairCancellationOutcome() const
{
  node_it p=getPNode();
  int nc,nd;

  if (p->getType()>getType())
    {
      nc=(nUpArcs()-1);if (nc<0) nc=0;
      nd=(p->nArcs()-p->nUpArcs()-1);if (nd<0) nd=0;
    }
  else
    {
      nc=(p->nUpArcs()-1);if (nc<0) nc=0;
      nd=(nArcs()-nUpArcs()-1);if (nd<0) nd=0;
    }

  nc*=nd;
  nd=nArcs()+p->nArcs()-1;
    
  return nc-nd;    
}

// weakest links may not be between a boundary and cross a boundary
inline NDcomplex_node::node_it NDcomplex_node::findWeakestSimpleNeighour(node_it self, arc_it null_arc, double downThreshold, double upThreshold, bool useRatio)
{
  arc_it a=arc_list;
  std::vector<node_it> up;
  std::vector<node_it> down;
  std::vector<node_it>::iterator it;
  node_it nu,nd;
  long upsize;
  long downsize;
  bool cond;

  up.reserve(numUpArcs);
  down.reserve(numArcs-numUpArcs);

  while(a!=null_arc)
    {
      node_it n=a->getOtherNode(self);
   
      if ((n->getFlags()&NODE_FLAG_BOUNDARY) != (flags&NODE_FLAG_BOUNDARY)) {
	a = a->getNext(self);
	continue;
      }

      if (n->getType()>self->getType()) {
	if (upThreshold<0) cond=true;
	else if (useRatio) cond = persistenceR(n)<=upThreshold;
	else  cond = persistence(n)<=upThreshold;

	if (cond) up.push_back(n);
      }
      else {
	if (downThreshold<0) cond=true;
	else if (useRatio) cond = persistenceR(n)<=downThreshold;
	else  cond = persistence(n)<=downThreshold;
	
	if (cond) down.push_back(n);
      }
      a = a->getNext(self);
    }

  upsize=up.size();
  if (upsize) {
    if (upsize==1) nu=up[0];
    else {
      std::sort(up.begin(),up.end(),compareItLess());
      nu=up[0];
 
      bool duplicate=false;
      for (int i=1;i<upsize;i++)
	{
	  if (nu!=up[i])
	    {
	      if (duplicate) {
		nu=up[i];
		duplicate=false;
	      }
	      else break;
	    }
	  else duplicate=true;
	}
      if (duplicate) upsize=0;
    }
  } 

  downsize=down.size();
  if (downsize) {
    if (downsize==1) nd=down[0];
    else {
      std::sort(down.begin(),down.end(),compareItMore());
      nd=down[0];
      
      bool duplicate=false;
      for (int i=1;i<downsize;i++)
	{
	  if (nd!=down[i])
	    {
	      if (duplicate) {
		nd=down[i];
		duplicate=false;
	      }
	      else break;
	    }
	  else duplicate=true;
	}
      if (duplicate) downsize=0;
    }
  }

  if (upsize==0)
    {
      if (downsize==0) return self;
      return nd;
    }
  if (downsize==0) return nu;
  
  if (self->persistence(nu)<self->persistence(nd)) return nu;
  else return nd;

}

inline std::string NDcomplex_node::getInfo(bool all) const
{
  std::ostringstream stm;
  stm<<"N";
  if (isInfinite()) stm<<("*");
  else {
    if (isPositive()) stm<<("+");
    if (isNegative()) stm<<("-");
  }
  stm << "("<<getId()<<","<<(int)getType();
  if (all) {
    stm << ","<< std::scientific << val <<","<<numArcs<<"("<<numUpArcs<<")";
    stm<<","<<persistence()<<"/"<<persistenceR();
    stm<<","<<(isBoundary()?"B":"")<<(isOut()?"O":"")<<(isInfinite()?"I":"");
  }
  stm<<")";
  return stm.str();
}

inline void NDcomplex_node::set(cellType cell, char typep, double valp, double valp2)
{
  cellID = cell;
  type=typep;
  val=valp;
  val2=valp2;
}

inline bool NDcomplex_node::check(NDcomplex_node::node_it self, NDcomplex_arc::arc_it null_arc)
{
  std::set<arc*> arcs_set;
  arc_it a=arc_list;

  if (a!=null_arc)assert (a->getPrev(self) == null_arc);

  while(a!=null_arc)
    {
      //printf ("checking arcs(%ld) : %ld<->%ld \n",&*a,&*(a->n[0]),&*(a->n[1]));
      if (arcs_set.find(&*a) != arcs_set.end()) return false;
      arcs_set.insert(&*a);

      if (a->getNext(self)!=null_arc)
	assert(a->getNext(self)->getPrev(self)==a);
      if (a->getNext(a->getOtherNode(self))!=null_arc)
	assert(a->getNext(a->getOtherNode(self))->getPrev(a->getOtherNode(self))==a);

      if (a->getPrev(self)!=null_arc)
	assert(a->getPrev(self)->getNext(self)==a);
      if (a->getPrev(a->getOtherNode(self))!=null_arc)
	assert(a->getPrev(a->getOtherNode(self))->getNext(a->getOtherNode(self))==a);


      a=a->getNext(self);
    }
  return true;
}

template <class excludeClass>
inline NDcomplex_arc::arc_it
NDcomplex_node::lowestDeltaArc(node_it self, arc_it null_arc,excludeClass exclude, char type)
{
  arc_it a=arc_list;
  arc_it result;
  double p,pmin;

  while ((a!=null_arc)&&(exclude(self,a->getOtherNode(self))))
    a=a->getNext(self);

  if (a==null_arc) return null_arc;

  result=a;
  pmin=a->delta();
  a=a->getNext(self);

  while(a!=null_arc)
    {
      if (!exclude(self,a->getOtherNode(self)))
	{
	  p=a->delta();

	  if (p<pmin)
	    {
	      result=a;
	      pmin=p;
	    }
	}
      a=a->getNext(self);
    }

  return result;
}


inline int NDcomplex_node::
countArcs(NDcomplex_node::node_it self, NDcomplex_arc::arc_it null_arc, int type) const
{
   
  int c=0;
  arc_it a=arc_list;
    
  while(a!=null_arc)
    {
      if ((type<0)||(a->getOtherNode(self)->getType() == type)) c++;
      a=a->getNext(self);
    }

  return c;
}

inline void NDcomplex_node::insertArc(NDcomplex_arc::arc_it &it, NDcomplex_node::node_it &self,const NDcomplex_arc::arc_it &null_arc)
{
  if (arc_list!=null_arc)
    arc_list->insertBeforeFirst(arc_list,it,self);

  arc_list = it;
  numArcs++;
  if (it->getOtherNode(self)->getType()>self->getType()) numUpArcs++;
}

inline bool NDcomplex_node::operator<(const node &other) const
{
  //return getId()<other.getId();
  if (getVal()!=other.getVal()) return getVal()<other.getVal();
  if (getType()!=other.getType())
    return getType()<other.getType();
  else
    return getId()<other.getId();
  
}

template <class OutputIterator> 
inline int NDcomplex_node::getNeighbors(node_it self, arc_it null_arc, OutputIterator it, char type)
{
  arc_it a=arc_list;
  int n=0;

  while(a!=null_arc)
    {
      if ((type<0)||(a->getOtherNode(self)->getType()==type))
	{
	  *it = a->getOtherNode(self);
	  n++;
	}
      a=a->getNext(self);
    }

  return n;
}

NDcomplex_node::arc_it
inline NDcomplex_node::findMultiArc(node_it self, node_it other, arc_it null_arc, int &nc, int nmax)
{
  arc_it a=arc_list;
  arc_it result=null_arc;
  
  nc=0;
  if (self->nArcs()<=other->nArcs())
    {
      while(a!=null_arc)
	{
	  if (a->getOtherNode(self) == other)
	    {
	      nc++;
	      result=a;
	      if (nc==nmax) return a;
	    }
	  a=a->getNext(self);
	}
    }
  else
    {
      a=other->getArc();
      while(a!=null_arc)
	{
	  if (a->getOtherNode(other) == self)
	    {
	      nc++;
	      result=a;
	      if (nc==nmax) return a;
	    }
	  a=a->getNext(other);
	}
    }
  
  return result;
}

NDcomplex_node::arc_it
inline NDcomplex_node::findArc(node_it self, node_it other, arc_it null_arc, bool &loop)
{
  arc_it a=arc_list;
  arc_it result=null_arc;
  
  loop=false;
  if (self->nArcs()<=other->nArcs())
    {
      while(a!=null_arc)
	{
	  if (a->getOtherNode(self) == other)
	    {
	      if (result==null_arc) result=a;
	      else {loop=true;return a;}
	      
	    }
	  a=a->getNext(self);
	}
    }
  else
    {
      a=other->getArc();
      while(a!=null_arc)
	{
	  if (a->getOtherNode(other) == self)
	    {
	      if (result==null_arc) result=a;
	      else {loop=true;return a;}
	      
	    }
	  a=a->getNext(other);
	}
    }
  
  return result;
}

NDcomplex_node::arc_it
inline NDcomplex_node::findArc(node_it self, node_it other, arc_it null_arc)
{
  arc_it a=arc_list;
  //arc_it result=null_arc;
  
  if (self->nArcs()<=other->nArcs())
    {
      while(a!=null_arc)
	{
	  if (a->getOtherNode(self) == other) return a;
	  a=a->getNext(self);
	}
    }
  else
    {
      a=other->getArc();
      while(a!=null_arc)
	{
	  if (a->getOtherNode(other) == self) return a;
	  a=a->getNext(other);
	}
    }
  
  return null_arc;
}

#endif
