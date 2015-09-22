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
#ifndef __NODES_MANIFOLD_HXX__
#define __NODES_MANIFOLD_HXX__

#include <stdlib.h>

#include <list>
#include <vector>
#include <map>
#include <fstream>
#include <functional>
#include "mytypes.h"

#define DIST(a,b) sqrt((b[0]-a[0])*(b[0]-a[0])+(b[1]-a[1])*(b[1]-a[1]))

template <typename idType>
class nodeGeometry {
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
	friend class nodeGeometry;
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
	
	baseGeometry(int m_type) {type=m_type;/*nrefs=0;*/}
	
	
      template <class OutputIterator> void get(OutputIterator oit, bool reverse = false) const
      {
	if (reverse)
	  {
	    std::not_equal_to<elements_rit> neq;
	    elements_rit it_end=elements.rend();
	    for (elements_rit it=elements.rbegin();neq(it,it_end);it++)
	      {
		*oit = *it;
	      }
	  }
	else
	  {
	    for (elements_it it=elements.begin();it!=elements.end();it++)
	      {
		*oit = *it;
	      }
	  }
      }

      void write(std::ofstream &str)
      {
	long i = elements.size();
	long j;

	str.write((const char *) &type, sizeof(char));
	str.write((const char *) &i, sizeof(long));
	for (j=0;j<i;j++) elements[j].write(str);
	//str.write((const char *) &elements.front(),i*sizeof(idType));
      }

      void read(std::ifstream &str)
      {
	long i,j;
	str.read((char *) &type, sizeof(char));
	str.read((char *) &i, sizeof(long));
	elements.resize(i);
	for (j=0;j<i;j++) elements[j].read(str);
	//str.read((char *) &elements.front(),i*sizeof(idType));
      }
    };
    
    geomBase_list baseList;

    
public:

    class geometry
    {
      
    private:
	friend class nodeGeometry;	
	// have to add a nnext here
       
	char type;
	geomBase_it elements;
	geom_it *next;
	int nrefs;

	void init() {next=NULL;nrefs=0;type=-1;}
	int clear(geom_list &lst, geomBase_list &b_lst)
	    {
		bool er=false;
		//printf("Nrefs=%d, next=%ld\n",nrefs,next);
		nrefs--;
		//printf ("NULL ?");fflush(0);
		if (next==NULL)
		{	
		    //printf ("YES\n");fflush(0);
		    if (nrefs==0)
		    {
			//printf ("Deleting elements\n");
			b_lst.erase(elements);
			elements = b_lst.end();
		    }
		    //printf("EL\n");
		}
		else
		{
		    //printf ("NO\n");fflush(0);
		    //printf("NEXTS %ld %ld %ld    %ld\n",&*next[0],&*next[1],&*next[2],next);
		    if (next[0]->clear(lst,b_lst)==0)
		    {lst.erase(next[0]);er=true;}
		    if (next[1]->clear(lst,b_lst)==0)
		    {lst.erase(next[1]);er=true;}

		    //free(next);next=NULL;
		    //printf("NEXT\n");
		}
		//assert((!er)||(nrefs==0));

		return nrefs;
	    }

      inline void checkAndFuse(geom_list &lst, geomBase_list &b_lst)
      {
	if (next==NULL) return;
	if ((next[0]->next!=NULL)||(next[1]->next!=NULL)) return;
	if ((next[0]->nrefs!=nrefs)||(next[1]->nrefs!=nrefs)) return;

	//printf("Simp %d\n",nrefs);

	std::vector<idType> new_el;
	get(std::back_inserter(new_el));

	int type = next[0]->elements->getType();
	next[0]->clear(lst,b_lst);lst.erase(next[0]);
	next[1]->clear(lst,b_lst);lst.erase(next[1]);
	free(next);next=NULL;

	elements = b_lst.insert(b_lst.end(),baseGeometry(type));
	elements->set(new_el.begin(),new_el.end());
      }

      inline void addRefAndFuse(geom_list &lst, geomBase_list &b_lst)
      {
	
	nrefs++;
	
	if (next!=NULL)
	  {
	    next[0]->addRefAndFuse(lst,b_lst);
	    next[1]->addRefAndFuse(lst,b_lst);
	    checkAndFuse(lst,b_lst);		
	  }
	//return nrefs;
      }

      inline void addRef()
      {
	
	nrefs++;
	//if (nrefs%(1000) ==0) {printf("nNODErefs=%d\n",nrefs);fflush(0);}
	if (next!=NULL)
	  {
	    next[0]->addRef();
	    next[1]->addRef();
	  }
	//return nrefs;
      }
      
	void set(geomBase_it it,geom_it null_geom)
	    {
		elements = it;
		next=NULL;
		addRef();
		type=it->getType();
	    }
	
      /*
	// Input Iterators elements are of type const geom_it
	template <class InputIterator> 
	geom_it set(InputIterator begin_it, InputIterator end_it)
	    {
		// A IMPLEMENTER
	    }
      */
      void set(const geom_it &g1,const geom_it &g2, bool noRecRef=false)
	    {
	     
	      next=(geom_it *)malloc(sizeof(geom_it)*2);
		
	      next[0]=g1;
	      next[1]=g2;
	      if (g1->getType() != g2->getType())
		{
		  fprintf (stderr,"ERROR: Manifold geometry must be defined with identical type of cells (%d!=%d)\n",g1->getType(),g2->getType());
		  assert(g1->getType() == g2->getType());
		  exit(0);
		}

	      //printf("%ld ",addRef());fflush(0);
	      if (!noRecRef) addRef();
	      else nrefs++;
		
	      type=g1->getType();
	    }
    public:
        int getType() const {return type;}

	geometry(const geometry& g)
	    {
		if (g.next!=NULL)
		{
		  next=(geom_it *)malloc(sizeof(geom_it)*2);
		  memcpy(next,g.next,sizeof(geom_it)*2);
		}
		else  elements=g.elements;
		nrefs=g.nrefs;
		type=g.getType();
	    }

	geometry& operator=(const geometry& g) {
	    if (this != &g) { 
		if (g.next!=NULL)
		{
		  if (next==NULL) next=(geom_it *)malloc(sizeof(geom_it)*2);
		    memcpy(next,g.next,sizeof(geom_it)*2);
		}
		else if (next!=NULL)
		{
		    free(next);
		    next=NULL;
		    elements=g.elements;
		}
		nrefs=g.nrefs;
		type=g.getType();
	    }
	    return *this;    
	}

	geometry() {init();}
	~geometry() {if (next!=NULL) free(next);}

      void get(std::set<baseGeometry*> &res)
      {
	if (next==NULL) 
	  res.insert(&(*elements));
	else
	  {
	    next[0]->get(res);
	    next[1]->get(res);		    
	  }
      }
      
      template <class OutputIterator> void get(OutputIterator it, bool reverse = false)
      {	
	std::set<baseGeometry*> result;
	get(result);
	
	for (typename std::set<baseGeometry*>::iterator set_it=result.begin();
	     set_it!=result.end();set_it++)
	  {
	    (*set_it)->get(it,reverse);
	  }
	
	/*
	//printf ("calling child ...");fflush(0);
	if (next==NULL) 
	  elements->get(it,reverse);
	else
	  {
	    next[0]->get(it,reverse);
	    next[1]->get(it,reverse);		    
	  }
	*/
      }
      
    };

private :

    geom_list manifold;

public:
    geom_it end() {return manifold.end();}

    void init()  {}

    nodeGeometry()
	{
	    init();
	}

    void erase(geom_it it)
	{
	    static int ct=0;
	    //printf("erasing %ld  (%d)\n",&*it,ct++);fflush(0);
	    if (it->clear(manifold,baseList) == 0)
	    {
		manifold.erase(it);
		//printf("DELETED\n");
	    }
	}

  //adds a manifold to the current one
  geom_it extend(geom_it &cur, geom_it &ext)
	{
	  geometry geom;
	  geom_it g=manifold.insert(manifold.end(),geom);
	  
	  g->set(cur,ext,true);
	  ext->addRefAndFuse(manifold,baseList);
	  //ext->addRef();
	  return g;
	}

  geom_it insert( geom_it &g1, geom_it &g2)
	{
	    geometry geom;
	    geom_it g =manifold.insert(manifold.end(),geom);
	    g->set(g1,g2);
	    //printf("inserted %ld\n",&*g);
	    return g;
	}
    
    // this inserts a manifold geometry from pre existing manifolds
    // InputIterators elements of type geom_it
    template <class InputIterator> 
    geom_it insert(InputIterator begin_it, InputIterator end_it)
	{
	    // A IMPLEMENTER
	  
	}

    // this inserts a new manifold geometry
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
  
  bool defrag()
  {
    long N = manifold.size();
    std::map<long,geom_it> refBy;
    typedef typename std::map<long,geom_it>::iterator refMap_it ;
    std::pair<refMap_it,bool> res;

    long Nu,Nt;

    for (geom_it it=manifold.begin(); it!=manifold.end();it++)
      {
	if (it->next!=NULL)
	  {
	    res = refBy.insert(std::make_pair((long) &(*it->next[0]),it));
	    if (!res.second) res.first->second=end();
	  
	    res =refBy.insert(std::make_pair((long) &(*it->next[1]),it));
	    if (!res.second) res.first->second=end();
	  }
      }
    Nu=Nt=0;
    for (refMap_it it = refBy.begin();it!=refBy.end();it++)
      {
	if (it->second!=end()) Nu++;
	Nt++;
      }

    printf (" (%.2f%% -%ld/%ld- unique)\n",100.*(double)Nu /(double)Nt,Nu,Nt);
    return true;
  }
  
  void write(std::ofstream &str, std::map<geometry *,long> &geom2id) {
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
	    n=2;
	    j[0]=geom2id[&(*it->next[0])];
	    j[1]=geom2id[&(*it->next[1])];
	    //j[2]=geom2id[&(*it->next[2])];
	  }

	str.write((const char *) &(it->type), sizeof(char));
	str.write((const char *) &(it->nrefs), sizeof(int));
	str.write((const char *) &n, sizeof(int));
	if (n) {
	  str.write((const char *) j, 2*sizeof(long));
	}
	else it->elements->write(str);
      }
      
  }
  void read(std::ifstream &str, std::vector<geom_it> &id2geom) {
    long i;
    long j[3];
    //std::vector<geom_it> id2geom;

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
	  str.read((char *) j, 2*sizeof(long));
	  it->next=(geom_it*)malloc(2*sizeof(geom_it));
	  it->next[0]=id2geom[j[0]];
	  it->next[1]=id2geom[j[1]];
	  //it->next[2]=id2geom[j[2]];
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

private:
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
    
	g =pool[poolID-1].insert(pool[poolID-1].end(),geom);
	g->set(gb,manifold.end());
      }
    else 
      {
	gb=baseList.insert(baseList.end(),geomBase);
	gb->set(begin_it,end_it);

	g =manifold.insert(manifold.end(),geom);
	g->set(gb,manifold.end());
     }

    return g;
  }

};


#endif
