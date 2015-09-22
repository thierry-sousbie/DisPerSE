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
#ifndef __ARCS_MANIFOLD_HXX__
#define __ARCS_MANIFOLD_HXX__

#include <stdlib.h>

#include <list>
#include <vector>
#include <fstream>
#include <functional>
#include "mytypes.h"

#define DIST(a,b) sqrt((b[0]-a[0])*(b[0]-a[0])+(b[1]-a[1])*(b[1]-a[1]))

template <typename idType = uint>
class arcGeometry {
public:
  class geometry;
  typedef typename std::list<geometry> geom_list;
  typedef typename std::list<geometry>::iterator geom_it;

private:
  friend class geometry;
  class baseGeometry;
  typedef typename std::list<baseGeometry> geomBase_list;
  typedef typename std::list<baseGeometry>::iterator geomBase_it;
   
  class baseGeometry
  {
  private:
    friend class arcGeometry;
    typedef typename std::vector<idType>::const_iterator elements_it;
    typedef typename std::vector<idType>::const_reverse_iterator elements_rit;
    std::vector<idType> elements;

    char type;
    void clear() {elements.clear();}
    template <class InputIterator> void set(InputIterator begin_it, InputIterator end_it)
    {
      elements.assign(begin_it,end_it);	
    }

  public:
    int getType() const {return type;}
	
    baseGeometry(int m_type) {type=m_type;}
	
    template <class OutputIterator> int get(OutputIterator oit, bool reverse=false) const
    {
      int nadded=0;
	
      if (reverse)
	{
	  std::not_equal_to<elements_rit> neq;
	  elements_rit prev = elements.rend();
	  elements_rit it_end = elements.rend();
	  for (elements_rit it=elements.rbegin();neq(it,it_end);++it)
	    {	
		
	      if ((prev == elements.rend())||(*it!=*prev))
		{
		  *oit = *it;
		  prev=it;
		  nadded++;
		}
		
	    }
	}
      else
	{
	  elements_it prev = elements.end();
	  for (elements_it it=elements.begin();it!=elements.end();it++)
	    {
		
	      if ((prev == elements.end())||(*it!=*prev))
		{
		  *oit = *it;
		  prev=it;
		  nadded++;
		}
		
	    }
	}
      return nadded;
    }

    void write(std::ofstream &str)
    {
      long i = elements.size();
      str.write((const char *) &type, sizeof(char));
      str.write((const char *) &i, sizeof(long));
      str.write((const char *) &elements.front(),i*sizeof(idType));
    }

    void read(std::ifstream &str)
    {
      long i;
      str.read((char *) &type, sizeof(char));
      str.read((char *) &i, sizeof(long));
      elements.resize(i);
      str.read((char *) &elements.front(),i*sizeof(idType));
    }
  };
    
  geomBase_list baseList;

    
public:

  class geometry
  {
      
  private:
      
    friend class arcGeometry;	
    char type;
    geomBase_it elements;
    geom_it *next;
    int nrefs;
    void init() {next=NULL;nrefs=0;type=-1;}
    int clear(geom_list &lst, geomBase_list &b_lst)
    {
      bool er=false;
      nrefs--;
      if (next==NULL)
	{	
	  if (nrefs==0)
	    {
	      b_lst.erase(elements);
	      elements = b_lst.end();
	    }
	}
      else
	{
	  if (next[0]->clear(lst,b_lst)==0)
	    {lst.erase(next[0]);er=true;}
	  if (next[1]->clear(lst,b_lst)==0)
	    {lst.erase(next[1]);er=true;}
	  if (next[2]->clear(lst,b_lst)==0)	
	    {lst.erase(next[2]);er=true;}
	}
		
      return nrefs;
    }
    inline void addRef()
    {
      nrefs++;
      if (next!=NULL)
	{
	  next[0]->addRef();
	  next[1]->addRef();
	  next[2]->addRef();
	}
    }
    void set(geomBase_it it,geom_it null_geom)
    {
      elements = it;
      next=NULL;
      addRef();
      type=it->getType();
    }

    void set(const geom_it &g1,const geom_it &g2,const geom_it &g3)
    {
      next=(geom_it *)malloc(sizeof(geom_it)*3);
      next[0]=g1;
      next[1]=g2;
      next[2]=g3;
      if ((g1->getType() != g2->getType())||(g2->getType() != g3->getType()))
	{
	  fprintf (stderr,"ERROR: Arc geometry must be defined with identical type of cells \n  %d %d %d\n",g1->getType(),g2->getType(),g3->getType());
	  assert(g1->getType() == g2->getType());
	  assert(g2->getType() == g3->getType());
	  exit(0);
	}

      addRef();
      type=g1->getType();
    }
  public:
    int getType() const {return type;}
    int getNRefs() const {return nrefs;}

    geometry(const geometry& g)
    {
      if (g.next!=NULL)
	{
	  next=(geom_it *)malloc(sizeof(geom_it)*3);
	  memcpy(next,g.next,sizeof(geom_it)*3);
	}
      else elements=g.elements;
      nrefs=g.nrefs;
      type=g.type;
    }

    geometry& operator=(const geometry& g) {
      if (this != &g) { 
	if (g.next!=NULL)
	  {
	    if (next==NULL) next=(geom_it *)malloc(sizeof(geom_it)*3);
	    memcpy(next,g.next,sizeof(geom_it)*3);
	  }
	else if (next!=NULL) {
	  free(next);
	  next=NULL;
	  elements=g.elements;
	}
	nrefs=g.nrefs;
	type=g.type;
      }
      return *this;    
    }

    geometry() {init();}
    ~geometry() {if (next!=NULL) free(next);}


    int get(std::vector<idType> &result, bool reverse = false, bool prt=false)
    {	
	
      int a=result.size();	
      int nadded=0;

      if (next==NULL) 
	{
	  nadded+=elements->get(back_inserter(result),reverse);	    
	}
      else
	{
	  if (reverse)
	    {		 
	      nadded+=next[2]->get(result,reverse,prt);
	      nadded+=next[1]->get(result,!reverse,prt);
	      nadded+=next[0]->get(result,reverse,prt);
	    }
	  else
	    {
	      nadded+=next[0]->get(result,reverse,prt);
	      nadded+=next[1]->get(result,!reverse,prt);
	      nadded+=next[2]->get(result,reverse,prt);
	    }	 
	}

      int ndeleted = (a+nadded-result.size())/2;
      a = result.size() - (nadded - ndeleted);
	
      if (prt) {
	printf("\nbefore : ");fflush(0);if (next==NULL)printf("(N)");
	for (int i=0;i<result.size();i++) {
	  if (i==a) printf(" | ");
	  printf("%ld ",result[i].id());
	}
	printf("\n");
      }
	
	
	
      if (a>0) {
	// special case ... too lazy to write comments, sorry ;)
	if ((a>1)&&(result[a+0]==result[a-1-1])) //remove a-0-1 
	  {
	    int i,j;
	    for (i=a-1,j=a;j<result.size();j++,i++) result[i]=result[j];
	    a--;result.resize(result.size()-1);
	  }
	else if ((result.size()>a+1)&&(result[a+1]==result[a-0-1])) //remove a
	  {
	    int i,j;
	    for (i=a,j=a+1;j<result.size();j++,i++) result[i]=result[j];
	    result.resize(result.size()-1);
	  }

	int  lim1=a;
	if (lim1>result.size()-a) lim1=result.size()-a;
	int n1=0;
	  
	while ((n1<lim1)&&(result[a+n1] == result[a-n1-1])) n1++;
	//n1--;
	if (prt) printf("lim1 = %d, n1=%d, a=%d \n",lim1,n1,a);
	if (n1>0)
	  {
	    int i,j;
	    for (i=a-n1+1,j=a+n1;j<result.size();j++,i++) result[i]=result[j];
	    result.resize(i);
	  }
	//printf ("a=%d n1=%d lim1=%d i=%d\n",a,n1,lim1);
      }
	
	

      if (prt) {
	printf("after : ");fflush(0);
	for (int i=0;i<result.size();i++) {
	  printf("%ld ",result[i].id());
	}
	printf("\n");
      }
      return nadded;
	
    }

    template <class OutputIterator> 
    void get(OutputIterator it, bool reverse = false)
    {
      std::vector<idType> tmp;
      get(tmp,reverse);

      for (typename std::vector<idType>::iterator rit=tmp.begin();rit!=tmp.end();rit++)
	{
	  *it = *rit;
	}

      return;
    }
      

    template <class OutputIterator> void getRaw(OutputIterator it, bool reverse = false)
    {	
	
      std::vector<idType> tmp;

      if (next==NULL) 
	elements->get(it,reverse);
      else
	{
	  if (reverse)
	    {
	      next[2]->getRaw(it,reverse);
	      next[1]->getRaw(it,!reverse);
	      next[0]->getRaw(it,reverse);
	    }
	  else
	    {
	      next[0]->getRaw(it,reverse);
	      next[1]->getRaw(it,!reverse);
	      next[2]->getRaw(it,reverse);
	    }
	}
    }
	
  };

private :

  geom_list manifold;

public:
  geom_it end() {return manifold.end();}

  void init()  {}

  arcGeometry()
  {
    init();
  }

  void erase(geom_it it)
  {
    static int ct=0;
	    
    if (it->clear(manifold,baseList) == 0)
      {
	manifold.erase(it);
      }
  }

  geom_it insert( geom_it &g1, geom_it &g2, geom_it &g3)
  {
    geometry geom;
    geom_it g =manifold.insert(manifold.end(),geom);
    g->set(g1,g2,g3);
    return g;
  }

  geom_it duplicate(geom_it &g)
  {
    g->addRef();
    return g;
  }
    
  template <class InputIterator> geom_it insert(InputIterator begin_it, InputIterator end_it, int m_type)
  {
    baseGeometry geomBase(m_type);
    geometry geom;
	    
    geomBase_it gb = baseList.insert(baseList.end(),geomBase);	 
    gb->set(begin_it,end_it);

    geom_it g =manifold.insert(manifold.end(),geom);
    g->set(gb,manifold.end());

    return g;
  }
  
  bool isPure()
  {
    for (geom_it it = manifold.begin();it!=manifold.end();it++)
      if (it->next!=NULL) return false;

    return true;
  }

  void write(std::ofstream &str,std::map<geometry *,long> &geom2id) {
    long i;
    long j[3];
    //std::map<geometry *,long> geom2id;
    int n;
    
    i=0;

    for (geom_it it = manifold.begin();it!=manifold.end();it++,i++)
      geom2id.insert(std::make_pair(&(*it),i));
      
    str.write((const char *) &i, sizeof(long));
    for (geom_it it = manifold.begin();it!=manifold.end();it++)
      {	
	if (it->next == NULL) 
	  n=0; 
	else 
	  {
	    n=3;
	    j[0]=geom2id[&(*it->next[0])];
	    j[1]=geom2id[&(*it->next[1])];
	    j[2]=geom2id[&(*it->next[2])];
	  }

	str.write((const char *) &(it->type), sizeof(char));
	str.write((const char *) &(it->nrefs), sizeof(int));
	str.write((const char *) &n, sizeof(int));
	if (n) {
	  str.write((const char *) j, 3*sizeof(long));
	}
	else it->elements->write(str);
      }
      
  }

  void read(std::ifstream &str, std::vector<geom_it> &id2geom) {
    long i;
    long j[3];

    str.read((char *) &i, sizeof(long));

    baseList.clear();
    manifold.clear();

    manifold.resize(i);
    id2geom.resize(i);
    i=0;

    for (geom_it it = manifold.begin();it!=manifold.end();it++,i++)
      id2geom[i]=it;
      
    for (geom_it it = manifold.begin();it!=manifold.end();it++)
      {
	int n;
	str.read((char *) &(it->type), sizeof(char));
	str.read((char *) &(it->nrefs), sizeof(int));
	str.read((char *) &n, sizeof(int));
	
	if (n) {
	  str.read((char *) j, 3*sizeof(long));
	  it->next=(geom_it*)malloc(3*sizeof(geom_it));
	  it->next[0]=id2geom[j[0]];
	  it->next[1]=id2geom[j[1]];
	  it->next[2]=id2geom[j[2]];
	}
	else {
	  geomBase_it bit = baseList.insert(baseList.end(),baseGeometry(it->type));
	  it->elements=bit;
	  bit->read(str);
	  it->next=NULL;
	}
      }
  }

  /* This is to allow multithreading */

private :
  std::vector<geomBase_list> bpool;
  std::vector<geom_list> pool;

public:
  int createPools(int N)
  {
    if (bpool.size())
      return -1;
      
    if (N>1) {
      bpool.resize(N-1);
      pool.resize(N-1);
    }
    return 0;
  }

  long commitPools()
  {
    int i;
    long N=bpool.size();

    for (i=0;i<N;i++)
      {
	baseList.splice(baseList.end(),bpool[i]);
	manifold.splice(manifold.end(),pool[i]);
      }

    bpool.clear();
    pool.clear();
    return N;    
  }

  // this inserts a new manifold geometry
  template <class InputIterator> geom_it insertInPool(InputIterator begin_it, InputIterator end_it, int m_type, int poolID)
  {
    baseGeometry geomBase(m_type);
    geometry geom;
    geom_it g;
    geomBase_it gb;

    if (poolID) 
      {
	gb=bpool[poolID-1].insert(bpool[poolID-1].end(),geomBase);
	gb->set(begin_it,end_it);
    
	g=pool[poolID-1].insert(pool[poolID-1].end(),geom);
	g->set(gb,manifold.end());
      }
    else 
      {
	gb=baseList.insert(baseList.end(),geomBase);
	gb->set(begin_it,end_it);

	g=manifold.insert(manifold.end(),geom);
	g->set(gb,manifold.end());
      }

    return g;
  }

};

#endif
