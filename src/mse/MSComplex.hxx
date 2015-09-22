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
#ifndef __MORSE_SMALE_COMPLEX_HEADER__
#define __MORSE_SMALE_COMPLEX_HEADER__

#include <sys/time.h>
#include <sys/param.h>
#include <sys/times.h>

#include "find_unordered_maps.hxx"

#include "network_interface.hxx"
#include "NDnet_interface.hxx"
#include "NDskel_interface.hxx"

#include "NDcomplex.hxx"
#include "NDskeleton.h"
#include "NDskel_tags.h"
#include "smooth.h"
#include "NDnetwork.h"
#include "NDnetwork_tags.h"
#include "myendian.h"
#include "updtQueue.hxx"
#include "probas.h"
#include "cells.hxx"

#include "unionfind.h"

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <algorithm>
#include <vector>
#include <set>
#include <map>

#include <numeric>
#include <iostream>
#include <fstream>
#include <functional>
#include <algorithm>

#ifdef USE_OPENMP
#include <parallel/algorithm>
#endif

#ifndef DEFAULT_NTHREADS
#define DEFAULT_NTHREADS 2
#endif

#define TO_NDSKL_UP (1<<0)
#define TO_NDSKL_DOWN (1<<1)
#define TO_NDSKL_INTER (1<<2)
#define TO_NDSKL_KEEP_STRUCTURE (1<<3)
#define TO_NDSKL_ALL ((1<<3)-1)

#define MSCOMPLEX_TAG "MORSE_SMALE_COMPLEX_MSC_v1.00"
  
//template <typename TypeT = char, typename IdT = uint>
class MSComplex {
public:  

  typedef NDcomplex::cellType cellType ;
  typedef cellType::typeT typeT;
  typedef cellType::idT idT;
  
  static const cellType MSC_CRITICAL;// = cellType(0,cellType::MAXIMUM);
  static const cellType MSC_UNPAIRED;// = cellType(1,cellType::MAXIMUM);
  static const cellType MSC_BOUNDARY_CRITICAL;// = cellType(2,cellType::MAXIMUM);
  
  typedef NetworkDataInterface<cellType> NetworkType;
  
  enum BoundaryType {BType_Sphere=1, BType_Torus=2, BType_Default=3, BType_NotSet=4, BType_Natural=0};
  
private:
  
  NetworkType *network; 

  typedef NetworkType::subCellsElementType subCellsElementType;
  typedef NetworkType::subCellsType subCellsType;
  typedef NetworkType::cellsInfoType cellsInfoType;
  typedef NetworkType::newCellsElementType newCellsElementType;
  typedef NetworkType::newCellsType newCellsType;

public:
  
  typedef class MSComplex parentClass;  

  class compareCell;
  
  typedef std::vector< std::vector< cellType  > > dbVec;
  typedef std::pair<double,uint> cellInfo;
  typedef uint cell_pair;
  typedef std::map< uint, NDcomplex::node_it > node_map;
  typedef node_map::iterator node_map_it;
  typedef std::pair< uint, NDcomplex::node_it > node_map_value;
  typedef std::pair<NDcomplex::node,NDcomplex::node> node_pair;
    
  typedef NDcomplex::node_it node_it;
  typedef NDcomplex::arc_it arc_it;
    
  typedef NDcomplex::arcGeom_it arcGeom_it;
  typedef NDcomplex::nodeGeom_it nodeGeom_it;

protected:
  NDcomplex cplx;
  MSComplex *self;

private:
  friend class compareCell;
  
  bool ownNet;
  bool manifoldWithBoundary;
  BoundaryType bType;
  int NThreads; 
  bool vertexAsMaxima;
  int store_manifolds; // up0,down0,up1,down1,....
  int store_arcsGeom;

  std::map<node_it,long,NDcomplex_node::compareItLess> nd2id_map;
  
  void initIdFromNodeIt()
  {
    long i=0;
    long ct=0;
    for (node_it n=cplx.nodes_begin();n!=cplx.nodes_end();n++,i++)
      if (nd2id_map.insert(std::make_pair(n,i)).second==false) ct++;

    if (ct) printf ("WARNING : %ld identical nodes\n",ct);
  }

  void resetIdFromNodeIt()
  {
    nd2id_map.clear();
  }

  bool IdFromNodeItIsAvailable()
  {
    if (nd2id_map.size()) return true;
    else return false;
  }
  
  long getIdFromNodeIt(node_it nd)
  {
    std::map<node_it,long,NDcomplex_node::compareItLess>::iterator it;
    it=nd2id_map.find(nd);
    if (it==nd2id_map.end()) {printf("HH: %s\n",nd->getInfo(true).c_str());return -1;}
    return it->second;
  }

public:
  void write(std::ofstream &str,std::string comment = std::string("No comment")) 
  {     
    str.write(MSCOMPLEX_TAG,strlen(MSCOMPLEX_TAG)); 
    int len = strlen(comment.c_str());
    str.write((char *)&len,sizeof(int)); 
    str.write((char *)comment.c_str(),len); 

    cplx.write(str);

    char b;
    int i;

    b = (char) bType;
    str.write((const char *)&b,sizeof(char));  
    
    b=vertexAsMaxima?true:false;
    str.write((const char *)&b,sizeof(char));
   
    str.write((const char *)&store_manifolds,sizeof(int));  
    str.write((const char *)&store_arcsGeom,sizeof(int));
  }

  void read(std::ifstream &str) 
  {
    std::string s;
    char txt[1024];
    
    str.read(txt,strlen(MSCOMPLEX_TAG));
    txt[strlen(MSCOMPLEX_TAG)]='\0';
        
    if (std::string(txt)!=std::string(MSCOMPLEX_TAG))
      {
	fprintf(stderr,"\nERROR in MSComplex: file is not a regular MSC file.\n");
	fprintf(stderr,"  tag : %s <> %s\n",txt,MSCOMPLEX_TAG);
	exit(0);
      }
    
    int len;    
    str.read((char *)&len,sizeof(int)); 
    char com[len+1];
    str.read((char *)com,len); 
    com[len]='\0';
    //printf("\n comment : %s\n",com);
    
    cplx.read(str);

    char b;  
    
    str.read((char *)&b,sizeof(char)); 
    bType = (BoundaryType)b;
    if (bType==BType_Natural) manifoldWithBoundary=true;
    else manifoldWithBoundary=false;
    
    str.read((char *)&b,sizeof(char));
    bool readv=(b)?true:false;
   
    vertexAsMaxima=(bool)b;
    
    str.read((char *)&store_manifolds,sizeof(int));  
    str.read((char *)&store_arcsGeom,sizeof(int));       
  }
 
  //regular members    
protected:
  double getValue(const node_it &node) const {return network->getValue(node->getCell());}
  ulong getNFacesTotal() const {ulong r=0;for (int i=0;i<=network->getNDims();i++) r+=network->getNFaces(i); return r;}     
  float *getPosition(const node_it &node,float *pos) const {return network->getPosition(node->getCell(),pos);}
  uint getNodeGroup(const node_it &node) const {return network->getNodeGroup(node->getCell());}
  char critIndexFromCellType(char index) {return (vertexAsMaxima)?(network->getNDims(true)-index):index;}
 
  
private:
  dbVec gradient_pairs;
  enum DGradientState {DG_available, DG_onDisk, DG_empty, DG_saveState, DG_loadState, DG_none};

  void setDiscreteGradientState(DGradientState statep, bool noComp=false)
  {
    
    DGradientState state = statep;
    static DGradientState cur_state = DG_empty;
    static DGradientState saved_state = DG_none;
    static std::string fname;

    if (state == cur_state) return;

    if (state == DG_saveState) {saved_state = cur_state;return;}
    if (state == DG_loadState) {state = saved_state; saved_state = DG_none;}

    if (state == DG_available)
      {
	if (cur_state==DG_empty) {if (!noComp) computeDiscreteGradient(gradient_pairs);}
	if (cur_state==DG_onDisk) //upload
	  {
	    FILE *f=fopen(fname.c_str(),"r");
	    gradient_pairs.clear();
	    do {
	      long N;
	      long i;
	      cellType::storageType tmp;
	      fread_sw(&N,sizeof(uint),1,f,0);
	      gradient_pairs.push_back(dbVec::value_type());
	      gradient_pairs.back().resize(N);
	      dbVec::value_type cur_ar = gradient_pairs.back();
	      
	      for (i=0;i<N;i++) {
		fread_sw(&tmp,sizeof(cellType::storageType),1,f,0);
		cur_ar[i].set(tmp);
	      }
	    } while (!feof(f));
	  }
      }

    if (state == DG_onDisk)
      {
	if (fname.size()==0)
	  {
	    char tmp_fname[255];
	    struct timeval time;
	    gettimeofday(&time, NULL);
	    sprintf (tmp_fname,"%ld%ld.dg.dump.tmp",(long)time.tv_sec,(long)time.tv_usec);
	    fname.assign(tmp_fname);
	  }
	if (cur_state==DG_empty) {
	  computeDiscreteGradient(gradient_pairs);
	  dumpDGradient(fname.c_str(),false); // to disk
	  gradient_pairs.clear();
	}
	if (cur_state==DG_available) {
	  dumpDGradient(fname.c_str(),false); // to disk
	  gradient_pairs.clear();
	}
      }

    if (state == DG_empty)
      {
	if (cur_state==DG_onDisk) {remove(fname.c_str());fname.assign("");} //delete
	if (cur_state==DG_available) {
	  gradient_pairs.clear();
	  if (fname.size()!=0) {
	    remove(fname.c_str());
	    fname.assign("");
	  }
	}
      }

    if (state != DG_none) cur_state = state;
  }

    
public:

  MSComplex(NetworkType *data, bool isOwner=false,int nthreads=DEFAULT_NTHREADS)
  {
    bType = BType_NotSet;
    manifoldWithBoundary=false;
    
    ownNet = isOwner;
    setNThreads(nthreads);
    network = data;
    
    init();
  }

  ~MSComplex()
  {
    setDiscreteGradientState(DG_empty);
    if (ownNet) delete network;
    network=NULL;
  }
  

protected:
  void init() {cplx.init(network->getNDims(),network->getNDims(true));}

public:
  bool netIsOwned() {return ownNet;}
  int setNThreads(int nthreads) {return NThreads=nthreads;}
    
public:
    
  template<class T> static bool Xor(const T &a,const  T &b) {return ((a&&(!b))||((!a)&&b));}

  class compareCell: public std::binary_function <cellType,cellType,bool> {
private:
    bool True;
    bool False;
    const parentClass *p;
public:
    compareCell(const parentClass *ptr, bool Less=true)
    {
      p=ptr;
      True=Less;
      False=!True;
    }
    bool operator()(const cellType& n1, const cellType& n2) const
      {
	if (n1==n2) return false;
	
	double dmax1,dmax2;
	bool f1 = p->network->isOut(n1);
	bool f2 = p->network->isOut(n2);
	
	if ((f2)&&(!f1)) return False;
	if ((f1)&&(!f2)) return True;	
	
	if (n1.type() == n2.type()) // both cells have same type
	  {
	    
	    double avg1,avg2;
	    
	    dmax1=p->network->getValueSum(n1,avg1);
	    dmax2=p->network->getValueSum(n2,avg2);
	    
	    if (dmax1==dmax2)
	      {
		if (avg1==avg2) return (n1.id()<n2.id())?True:False;
		else return (avg1<avg2)?True:False;
	      }
	    else return (dmax1<dmax2)?True:False;	       	    
	  }
	else
	  {	    
	    if (n1.type() < n2.type()) return True;
	    return False;	 
	  }
      }
  };

private:
 
  uint getSortedCells(std::vector< cellType > &v, char type, bool withBoundary=true, bool descendingFiltration=false) const
  {
    compareCell cmp(this,!descendingFiltration);
    idT nfaces = network->getNFaces(type);
    v.clear();
    cellType cell;

    if (!withBoundary)
      {
	for (idT j=0;j<nfaces;j++)
	  {
	    cell.set(type,j);
	    if (network->isBoundary(cell)||network->isOut(cell)) continue; 
	    v.push_back(cell);	
	  }
      }
    else
      {
	for (idT j=0;j<nfaces;j++) 
	  {
	    cell.set(type,j);
	    if (!network->isOut(cell)) 
	      v.push_back(cell);
	  }  
      }
 
#ifdef USE_GNU_PSORT
    __gnu_parallel::sort(v.begin(), v.end(),cmp);
#else
    std::sort(v.begin(),v.end(),cmp);
#endif

    return v.size();
  }


  cellType getLowestPairableCofacet(cellType cell,const dbVec &pairs, bool onBoundary=false, bool descendingFiltration=false) const
  {
    int i,j,k;
    compareCell cmp(this,!descendingFiltration); // (,true) <=> 'operator<()'
    cellType newindex=MSC_CRITICAL;
	
    bool boundary = network->isBoundary(cell);

    if (network->isOut(cell)) return newindex;

    static std::vector<cellType> flist,coflist;

    if (vertexAsMaxima==descendingFiltration) network->getCofaces(cell, coflist); 
    else network->getFaces(cell, coflist);
      

    for (i=0;i<coflist.size();i++)
      {			    
	if (pairs[coflist[i].type()][coflist[i].id()] != MSC_UNPAIRED) continue; // is it available ???
	if ((onBoundary)&&(!network->isBoundary(coflist[i]))) continue;

	if (boundary)
	  {
	    //if (!network->isBoundary(type+1,coflist[i])) continue;
	  }
	else 
	  {
	    if (network->isOut(coflist[i])) continue;
	    //if (network->isBoundary(type+1,coflist[i])) continue;
	  }
	    
	// select only cofacets with exactly one available face.
	if (vertexAsMaxima==descendingFiltration) network->getFaces(coflist[i],flist);
	else network->getCofaces(coflist[i],flist);

	for (j=0,k=0;j<flist.size();j++)
	  {
	    if ((onBoundary)&&(!network->isBoundary(flist[j]))) continue;
	    if (pairs[flist[j].type()][flist[j].id()] == MSC_UNPAIRED)
	      {
		if (++k>1) break;
	      }
	  }
	if (k!=1) continue;
	    
	if (newindex==MSC_CRITICAL) newindex = coflist[i];
	else if (cmp(coflist[i],newindex)) newindex=coflist[i];
	  
      }

    return newindex;
  }

  std::vector<cellType>
  getPairsMST(std::vector<cellType> &mst, char type, bool decreasingType, double threshold, std::vector<cellType> &cells_arr)
  {  
    static std::vector<cellType> flist,coflist;
    cellType coface;
    long npairs=0;
    long i,j;
    bool dirT2B=(vertexAsMaxima==decreasingType)?true:false;
    compareCell cmpLess(this,dirT2B);    
    std::vector<cellType> component;
    UF_type *uf;
    UF_type *ufc;
    int newT=(decreasingType)?(type-1):type+1;
    
    printf ("      Computing MST ... ");fflush(0);
    
    long nT = network->getNFaces(type);
    long nNewT = network->getNFaces(newT);

    uf = uf_create(nNewT+1,true);
    ufc = uf_create(nNewT+1,true);
    
    mst.resize(nT,MSC_UNPAIRED);
      
    for (i=0;i<cells_arr.size();i++)
      {
	if (network->isOut(cells_arr[i])) continue;

	if (decreasingType) {
	  network->getFaces(cells_arr[i],flist);
	}
	else {
	  network->getCofaces(cells_arr[i],flist);
	}
	
	bool bypass=false;
	
	if ((!decreasingType)&&(network->isBoundary(cells_arr[i])))
	  {
	    static std::vector<cellType> flist2;
	    for (j=0;j<flist.size();j++)
	      if (!network->isOut(flist[j])) flist2.push_back(flist[j]);

	    if (flist2.size()==0) continue;
	    else flist.assign(2,flist2[0]);
	    bypass=true;
	  }
	
	
	long ga = uf_find(uf,flist[0].id());
	long gb = uf_find(uf,flist[1].id());
	cellType a(newT,ga);
	cellType b(newT,gb);

	
	if (bypass) gb=nNewT;
       
	if (ga!=gb)
	  {
	  //here, we have a negative cell (pair or not) 
	  // -> store the component it destroys	  

	    if (ga==nNewT) {std::swap(a,b);std::swap(ga,gb);}
	    else if (gb==nNewT) {}
	    else if (cmpLess(b,a)) {std::swap(a,b);std::swap(ga,gb);}

	  double delta;
	  
	  if (dirT2B)
	    {
	      delta = network->getValue(a)-network->getValue(cells_arr[i]);
	    }
	  else
	    {
	      delta = network->getValue(cells_arr[i])-network->getValue(a);
	    }	  
	  
	  uf_gr_union(uf,gb,ga);

	  bool keepPair=false;

	  if (delta>threshold)
	    {
	      keepPair=true;
	    }	  

	  if (network->isBoundary(cells_arr[i]) != network->isBoundary(a)) 
	    keepPair=true;

	  // forces boundary cells to be critical
	  // if (network->isBoundary(cells_arr[i])) keepPair=true;
	   
	  if (keepPair)
	    {
	      mst[cells_arr[i].id()] = a;
	    }
	  else 
	    {
	      mst[cells_arr[i].id()] = MSC_CRITICAL;
	      ga = uf_find(ufc,flist[0].id());
	      gb = uf_find(ufc,flist[1].id());
	      a = cellType(newT,ga);
	      b = cellType(newT,gb);
	      
	      if (ga==nNewT) {std::swap(a,b);std::swap(ga,gb);}
	      else if (gb==nNewT) {}
	      else if (cmpLess(b,a)) std::swap(ga,gb);

	      uf_gr_union(ufc,gb,ga);
	    }
	}
	else {}//here, we have a positive cell
	//printf("\n");
      }
    uf_free(&uf);

    printf ("(components)");fflush(0);
    
    for (i=0;i<nNewT;i++)
      {
	if (network->isOut(cellType(newT,i))) continue;
	if (uf_isOwnGroup(ufc,i))
	  component.push_back(cellType(newT,i));	
      }
    
    uf_free(&ufc);  
    printf(" done.(%ld/%ld comp.)\n",component.size(),nNewT);

    return component;
  }

  void getSaddleSaddleDG(std::vector<cellType> &cells_arr,char type,bool decreasingType,dbVec &pairs)
  {
    bool dirT2B=(vertexAsMaxima==decreasingType)?true:false;
    compareCell cmpLess(this,dirT2B);   
    //int ndims=network->getNDims();    
    int newT=(decreasingType)?(type-1):type+1;
    long i,j,k;
    std::vector<cellType> flist;
    
    long nT = network->getNFaces(type);
    long nNewT = network->getNFaces(newT);

    printf ("      Computing saddle-saddle DG ... ");fflush(0);

    for (i=0;i<cells_arr.size();i++)
      {
	if (network->isOut(cells_arr[i])) continue;
	if (pairs[cells_arr[i].type()][cells_arr[i].id()] != MSC_UNPAIRED) continue;
	if (decreasingType) 
	  network->getFaces(cells_arr[i],flist);
	else 
	  network->getCofaces(cells_arr[i],flist);

	int nf=0;
	int minid=0;
	bool isb=network->isBoundary(cells_arr[i]);
	for (j=0;j<flist.size();j++)
	  {
	    if (pairs[flist[j].type()][flist[j].id()] != MSC_UNPAIRED ) 
	      continue;
	    if ( isb != network->isBoundary(flist[j]) ) 
	      continue;
	    if (network->isOut(flist[j])) continue;

	    if (j!=nf) flist[nf]=flist[j];
	    if (nf)
	      {
		if (cmpLess(flist[nf],flist[minid])) minid=nf;
	      }
	    nf++;
	  }

	if (nf)
	  {
	    if (network->getValue(flist[minid]) == network->getValue(cells_arr[i]))
	      {
		pairs[flist[minid].type()][flist[minid].id()]=cells_arr[i];
		pairs[cells_arr[i].type()][cells_arr[i].id()]=flist[minid];
	      }
	  }
      }

    printf("done.\n");
  }

  int buildDGComponent_rec(cellType cur,bool decreasingType,const std::vector<cellType> &mst,dbVec &pairs,int ct)
  {
    std::vector<cellType> flist;
    long i,j;
    cellType c;
    
    if (decreasingType) network->getFaces(cur,flist);
    else network->getCofaces(cur,flist);
    
    
    if (pairs[flist[0].type()][flist[0].id()] == MSC_UNPAIRED) 
      c=flist[0];
    else if (pairs[flist[1].type()][flist[1].id()] == MSC_UNPAIRED) 
      c=flist[1];
    else 
      {
	fprintf(stderr,"ERROR in buildDGComponent_rec !\n");
	return -1;
      }
    
   
    pairs[c.type()][c.id()] = cur;
    pairs[cur.type()][cur.id()] = c;
    
    if (!decreasingType) network->getFaces(c,flist);
    else network->getCofaces(c,flist);

    j=0;
    for (i=0;i<flist.size();i++)
      {
	if (mst[flist[i].id()] != MSC_CRITICAL) 
	  {
	    continue;
	  }
	if (pairs[flist[i].type()][flist[i].id()] != MSC_UNPAIRED) continue;

	buildDGComponent_rec(flist[i],decreasingType,mst,pairs,ct+1);
      }
    return 0;
  }
 
  void buildDGComponent(cellType start,bool decreasingType,const std::vector<cellType> &mst,dbVec &pairs)
  {
    pairs[start.type()][start.id()]=MSC_CRITICAL;
    std::vector<cellType> flist;
    int i;

    if (!decreasingType) network->getFaces(start,flist);
    else network->getCofaces(start,flist);
    
    for (i=0;i<flist.size();i++)
      {
	if (mst[flist[i].id()] != MSC_CRITICAL)
	  {
	    continue;
	  }
	if (pairs[flist[i].type()][flist[i].id()] != MSC_UNPAIRED) continue;

	int ct=buildDGComponent_rec(flist[i],decreasingType,mst,pairs,0);

      }
  }

  void computeDiscreteGradient(dbVec &pairs, bool descendingFiltration=false)
  {
    int ndims=network->getNDims();    
    int ndims_net=network->getNDims(true);    
  
    std::vector<cellType> *cells;
    cellType coface;
    char type=-1;
    char typeL;
    bool haveBoundary=false;
    long i,oldi,j,b;
    std::vector<cellType> components;
    std::vector<cellType> mst;
    std::vector<cellType> cells_arr;

    printf ("Computing discrete gradient for %ld cells:\n",getNFacesTotal());

    fflush(0);

    pairs.clear();
    pairs.resize(ndims_net+1);
    for (int i=0;i<=ndims_net;i++)
      pairs[i].resize(network->getNFaces(i),MSC_UNPAIRED);
    
    bool reverseFiltration=false;	
    bool decreasingType;
  
    for (i=0;i<std::min(ndims_net,2);i++)
      {
	if (i==0)
	  {
	    type=ndims_net-1;
	    decreasingType=false;
	    reverseFiltration=false;
	  }
	else
	  {
	    type=1;
	    decreasingType=true;
	    reverseFiltration=(ndims_net==2)?true:false;
	  }

	printf("   Identifying (%d,%d)-cell pairs:\n",type,type+((decreasingType)?-1:1));fflush(0);   
	mst.clear();
	
	if (reverseFiltration) 
	  std::reverse(cells_arr.begin(),cells_arr.end());
	else 
	  {
	    bool dirT2B=(vertexAsMaxima==decreasingType)?true:false;
	    
	    printf ("      Sorting %ld %d-cells (%s) ... ",network->getNFaces(type),type,(dirT2B)?"desc.":"asc.");fflush(0);
	    getSortedCells(cells_arr,type,true,dirT2B);
	    printf("done.\n");
	  }
	
	components = getPairsMST(mst,type,decreasingType,0,cells_arr);
	printf ("      Computing discrete Gradient (%ld comp.) ... ",components.size());fflush(0);
	
	#pragma omp parallel for
	for (j=0;j<components.size();j++)
	  buildDGComponent(components[j],decreasingType,mst,pairs);

	printf("done.\n");
	
      }    
    
    if (ndims_net==3)
      {
	decreasingType=!decreasingType;
	reverseFiltration=!reverseFiltration;
	printf("   Identifying (%d,%d)-cell pairs:\n",type,type+((decreasingType)?-1:1));fflush(0);
	
	mst.clear();
	
	if (reverseFiltration) 
	  std::reverse(cells_arr.begin(),cells_arr.end());
	else 
	  {
	    bool dirT2B=(vertexAsMaxima==decreasingType)?true:false;
	    
	    printf ("      Sorting %ld %d-cells (%s) ... ",network->getNFaces(type),type,(dirT2B)?"desc.":"asc.");fflush(0);
	    getSortedCells(cells_arr,type,true,dirT2B);
	    printf("done.\n");
	  }

	getSaddleSaddleDG(cells_arr,type,decreasingType,pairs);
      }

    std::vector<long> ncrits(ndims_net+1,0);
    std::vector<long> nplus(ndims_net+1,0);
    for (i=0;i<=ndims_net;i++)
      for (j=0;j<pairs[i].size();j++)
	{
	  if (network->isOut(cellType(i,j))) continue;
 
	  if (pairs[i][j]==MSC_UNPAIRED) 
	    {
	      pairs[i][j]=MSC_CRITICAL;
	      nplus[i]++;
		
	    }
	  if (pairs[i][j]==MSC_CRITICAL) ncrits[i]++;
	}

    printf("   Critical cells : ");
    for (i=0;i<=ndims_net;i++)
      {
	printf("%ld(+%ld) %ld-cells%s",ncrits[i],nplus[i],i,(ndims_net==i)?".\n":", ");
      }
  }

  //int ComputeManifold(const dbVec &pairs,node_map *critPoints)
  int ComputeManifold(const dbVec &pairs, std::vector<node_map> &critPoints)
  {  
    //int type = node->getType();
    int n=0;
    char tmpM[255];
    char tmpA[255];
    int ndims = network->getNDims();
    int ndims_net = network->getNDims(true);
    std::vector<int> narcs(ndims_net,0);
  
    for (int type=0;type<=ndims_net;type++)
      {
	bool doA=false;
	bool doD=false;
       
	bool doAA=false;
	bool doDA=false;
	bool doAAG=false;
	bool doDAG=false;
       
	if (store_manifolds) {
	  if (type!=ndims_net) doA=true; 
	  if (type!=0) doD=true;
	}
	//if (store_arcsGeom) {

	if (type == ndims_net-1) {doAA=true;if (store_arcsGeom&(1<<0)) doAAG=true;}
	if (type == 1) {doDA=true;if (store_arcsGeom&(1<<2)) doDAG=true;}
	if ((type!=ndims_net)&&(type>1)) {doDA=true;if (store_arcsGeom&(1<<1)) doDAG=true;}

	if (doAA) doA=true;
	if (doDA) doD=true;
	//}
       
	if ((!doA)&&(!doD)) {
	  //printf("    * Type %d nodes: SKIPPED.\n",type);
	  printf("    * %s: SKIPPED.\n",cplx.nodeType2Str(type).c_str());
	  continue;
	}

	if (!store_manifolds) sprintf(tmpM,"no manifolds");
	else if (doA&&doD) sprintf(tmpM,"A./D. manifolds");
	else if (doA) sprintf(tmpM,"A. manifolds");
	else if (doD) sprintf(tmpM,"D. manifolds");

	if ((!doAA)&&(!doDA)) sprintf(tmpA,"no arcs");
	else if (doAA&&doDA) sprintf(tmpA,"A.%s/D.%s arcs",(doAAG)?"(G+)":"",(doDAG)?"(G+)":"");
	else if (doAA) sprintf(tmpA,"A.%s arcs",(doAAG)?"(G+)":"");
	else if (doDA) sprintf(tmpA,"D.%s arcs",(doDAG)?"(G+)":"");
      
	//printf("    * Type %d nodes: %s, %s ... ",type,tmpM,tmpA);fflush(0);
	printf("    * %s: %s, %s ... ",cplx.nodeType2Str(type).c_str(),tmpM,tmpA);fflush(0);

#ifdef USE_OPENMP
	int nth=1;
        #pragma omp parallel
	nth=omp_get_num_threads();

	cplx.createPools(nth);

	std::vector<node_it> allnodes;
	allnodes.reserve(critPoints[type].size());
	for (node_map_it nmit = critPoints[type].begin();nmit != critPoints[type].end();nmit++)
	  allnodes.push_back(nmit->second);

	long i;
	#pragma omp parallel for 
	for (i=0;i<allnodes.size();i++)
	  {
	    node_it node = allnodes[i];
	    int curThread = omp_get_thread_num(); 
	    //printf("curThread = %d\n",curThread);
#else
	for (node_map_it nmit = critPoints[type].begin();nmit != critPoints[type].end();nmit++)
	  {
	    node_it node = nmit->second;
	    int curThread=0;
#endif

	    if (doA)
	      {
		std::vector<cellType> cells;
		std::vector< std::vector<cellType>* > skl;
		bool ascending=true;

		ComputeAscendingManifold(node->getCell(),pairs,back_inserter(cells),skl);
		
		if (store_manifolds)
		  {
		    std::set<cellType> ucells;
		    ucells.insert(cells.begin(),cells.end());
		    nodeGeom_it it = cplx.insertNodeGeometry(ucells.begin(),ucells.end(),node->getType(),curThread);
		    node->setManifoldGeometry(it, ascending);
		  }
		else node->setManifoldGeometry(cplx.nodesGeom_end(), ascending);
	      
		if (doAA) 
		  {
		    for (std::vector< std::vector<cellType>* >::iterator it = skl.begin();it!=skl.end();it++)
		      {
			std::vector<cellType> &path = *(*it);
			arc_it arc;
			arcGeom_it git;
		      
			arc = cplx.addArc(node,critPoints[type+1].find(path[0].id())->second,curThread);
			/*	
				if (arc->endpoints().first->getVal()>arc->endpoints().second->getVal())
				printf("ERROR: %g > %g\n",arc->endpoints().first->getVal(),arc->endpoints().second->getVal());
			*/
			#pragma omp atomic
			narcs[type]++;
			if (!doAAG) {
			  arc->setManifoldGeometry(cplx.arcsGeom_end());
			  continue;
			} 
		      
			if (arc->endpoints().first == node)
			  {
			    std::vector<cellType>::reverse_iterator begin = path.rbegin();
			    std::vector<cellType>::reverse_iterator end = path.rend();
		      
			    git = cplx.insertArcGeometry(begin,--end,curThread);
			    arc->setManifoldGeometry(git);
			  }
			  
			else
			  {
			    std::vector<cellType>::iterator begin = path.begin();
			    std::vector<cellType>::iterator end = path.end();
		      
			    git = cplx.insertArcGeometry(begin,--end,curThread);
			    arc->setManifoldGeometry(git);
			  }
			  
		      } 
		  }
		for (std::vector< std::vector<cellType>* >::iterator it = skl.begin();it!=skl.end();it++)
		  delete *it;
	      }

	    if (doD)
	      {
		std::vector<cellType> cells;
		std::vector< std::vector<cellType>* > skl;
		bool ascending=false;

		ComputeDescendingManifold(node->getCell(),pairs,back_inserter(cells),skl);
	
		if (store_manifolds)
		  {
		  
		    std::set<cellType> ucells;
		    ucells.insert(cells.begin(),cells.end());
		    //printf("D(%d/%ld)",type,cells.size());fflush(0);
		    nodeGeom_it it = cplx.insertNodeGeometry(ucells.begin(),ucells.end(),node->getType(),curThread);
		    //if (it ==cplx.nodesGeom_end() ) printf("END\n");
		    node->setManifoldGeometry(it, ascending);
		  }
		else node->setManifoldGeometry(cplx.nodesGeom_end(), ascending);
	    
		if (doDA) 
		  {
		    for (std::vector< std::vector<cellType>* >::iterator it = skl.begin();it!=skl.end();it++)
		      {
			std::vector<cellType> &path = *(*it);
			arc_it arc;
			arcGeom_it git;
		      
			arc = cplx.addArc(node,critPoints[type-1].find(path[0].id())->second,curThread);
			/*
			  if (arc->endpoints().first->getVal()<arc->endpoints().second->getVal())
			  printf("ERROR: %g < %g\n",arc->endpoints().first->getVal(),arc->endpoints().second->getVal());
			*/
			#pragma omp atomic
			narcs[type-1]++;
			if (!doDAG) {
			  arc->setManifoldGeometry(cplx.arcsGeom_end());
			  continue;
			} 
		      
			if (arc->endpoints().first == node)
			  {
			    std::vector<cellType>::reverse_iterator begin = path.rbegin();
			    std::vector<cellType>::reverse_iterator end = path.rend();
		      
			    git = cplx.insertArcGeometry(begin,--end,curThread);
			    arc->setManifoldGeometry(git);
			  }
			  
			else
			  {
			    std::vector<cellType>::iterator begin = path.begin();
			    std::vector<cellType>::iterator end = path.end();
		      
			    git = cplx.insertArcGeometry(begin,--end,curThread);
			    arc->setManifoldGeometry(git);
			  }
			  
		      } 
		  }
		for (std::vector< std::vector<cellType>* >::iterator it = skl.begin();it!=skl.end();it++)
		  delete *it;
	      }
		  
	  }
	printf("done.\n");
	cplx.commitPools();
      } 
    
    printf ("    Computed %d min/saddle arcs, %d max/saddle arcs (%d total).\n",narcs[0],narcs[ndims_net-1],cplx.num_arcs());
    //exit(0);
    return n;
  }

  void computeMSComplex(const dbVec &pairs)
  {
    int ndims=network->getNDims();
    int ndims_net=network->getNDims(true);
    //node_map critPoints[ndims+1];
    std::vector<node_map> critPoints(ndims+1);
    int ncrit_tot=0;
    std::vector<int> ncrits(ndims+1,0);
    int type;
    uint base;
    cellType cell; 

    printf("Computing discrete Morse-Smale complex: \n");fflush(0);    
	
    for (type=0; type<=ndims_net; type++)
      {
	char critIndex = critIndexFromCellType(type);
	idT nfaces = network->getNFaces(type);
	for (idT i=0; i<nfaces; i++)
	  {
	    cell.set(type,i);
	    if (pairs[type][i]==MSC_CRITICAL)
	      {
		char flags=0;
		
		if (network->isBoundary(cell)) flags |= NODE_FLAG_BOUNDARY;
		if (network->isOut(cell)) 
		  {
		    //printf("ERROR\n");
		    flags |= NODE_FLAG_OUT;
		    continue;
		  }
		
		double val,avg;

		val=network->getValueSum(cell,avg);
		critPoints[critIndex].insert(node_map_value(i,cplx.addNode(cell,val,avg,flags,critIndex)));
		ncrits[critIndex]++;
	      }
	  }
      }

    for (int i=0;i<=ndims_net;i++) ncrit_tot+=ncrits[i];

    printf ("    %d critical points: %d min",ncrit_tot,ncrits[0]);
    for (int i=1;i<ndims_net;i++)
      printf (", %d s%d",ncrits[i],i);
    printf (", %d max.\n",ncrits[ndims_net]);


    int narcs1 =0;
    int narcs2 =0;

    ComputeManifold(pairs,critPoints);//,it==critPoints[type].end());

  }

public:
  void simplifyComplex(double persistence, bool skipLoops=false, bool solveConflicts=false)
  {
    
    printf("Simplifying complex (p>%g):\n",persistence);fflush(0);
    printf ("    ");
    cplx.clearNodeTags();
    printf ("    Tagging arcs ... ");fflush(0);
    
    for (node_it n1= cplx.nodes_begin(); n1!=cplx.nodes_end();n1++)
      {
	if (n1->isInfinite()) {continue;}	
	node_it n2 = n1->getPNode();
	if (n2==cplx.nodes_end()) {continue;}			
	if (n2==n1) {continue;}
	if (n1->persistence()>persistence) {continue;}
	if (n2->getType()>n1->getType()) {continue;}
	
	n1->setTag();
      }
    printf("done.\n");
    cplx.cancelTaggedNodes(!skipLoops, solveConflicts);
  }

private:
  
  void simplifyFromPairs(std::vector<double> level, bool useRatio, bool skipLoops=false, bool dumpConflicts=false)
  {
    int ndims=network->getNDims();
    
    printf ("Cancelling pairs with persistance%s < [%.2e",(useRatio)?" ratio":"",level[0]);
    for (int i=1;i<ndims;i++) printf(",%.2e",level[i]);
    printf("].\n");
 
    cplx.clearNodeTags();
    printf ("    Tagging arcs ... ");fflush(0);
    for (node_it n1= cplx.nodes_begin(); n1!=cplx.nodes_end();n1++)
      {
	if (n1->isInfinite()) {continue;}	
	node_it n2 = n1->getPNode();
	if (n2==cplx.nodes_end()) {continue;}			
	if (n2==n1) {continue;}
	if (n2->getType()>n1->getType()) {continue;}
	if (useRatio)
	  {
	    if (n1->persistenceR()>level[n2->getType()]) 
	      continue;
	  }
	else
	  {
	    if (n1->persistence()>level[n2->getType()]) 
	      continue;
	  }	
	n1->setTag();
      }
    printf("done.\n");    

    std::vector<node_it> imp;
    std::vector<node_it> con;

    cplx.cancelTaggedNodes_Smart(back_inserter(imp),back_inserter(con),!skipLoops,false);
       
    int constart=imp.size();
    imp.insert(imp.end(),con.begin(),con.end());
    if (imp.size()&&dumpConflicts) 
      {
	std::vector<arc_it> arcs_tab;
	std::vector<dumpGeomNode> nodes_tab;
	char fname[255];
	int i=0;

	initIdFromNodeIt();
	
	for (std::vector<node_it>::iterator it = imp.begin(); it!=imp.end();it++,i++)
	  {
	    arcs_tab.clear();
	    nodes_tab.clear();
	    
	    node_it curn = *it;
	    node_it othern=(*it)->getPNode();

	    printf ("nodes : %s %s\n",curn->getInfo(true).c_str(),othern->getInfo(true).c_str());
	     
	    int narcs=0;
	    for (arc_it arc= curn->getArc(); arc!=cplx.arcs_end(); arc=arc->getNext(curn))
	      if (arc->getOtherNode(curn) == othern) {narcs++;arcs_tab.push_back(arc);}

	    
	    for (arc_it arc= curn->getArc(); arc!=cplx.arcs_end(); arc=arc->getNext(curn))
	      if (arc->getOtherNode(curn)->getType() != othern->getType()) {narcs++;arcs_tab.push_back(arc);}

	    for (arc_it arc= othern->getArc(); arc!=cplx.arcs_end(); arc=arc->getNext(othern))
	      if (arc->getOtherNode(othern)->getType() != curn->getType()) {narcs++;arcs_tab.push_back(arc);}
	    
	    sprintf (fname,"impossible_G%d",i);
	    if (i>=constart)  sprintf (fname,"impossible_conflictG%d",i-constart);
	    printf ("Dumping %d non-cancellable arcs.\n",narcs);
	    dumpGeometry(arcs_tab.begin(),arcs_tab.end(),
			 nodes_tab.begin(),nodes_tab.end(),
			 fname,0,false,false,true);
	   
	    arcs_tab.clear();
	    nodes_tab.clear();
	    
	    if (curn->getType() > othern->getType())
	      {
		nodes_tab.push_back(dumpGeomNode(curn,dumpGeomNode::Down));

		sprintf (fname,"isurfaceD_G%d",i);
		if (i>=constart)  sprintf (fname,"isurfaceD_conflictG%d",i-constart);
		
		dumpGeometry(arcs_tab.begin(),arcs_tab.end(),
			     nodes_tab.begin(),nodes_tab.end(),
			     fname,0,false,false,true);
		nodes_tab.clear();
		nodes_tab.push_back(dumpGeomNode(othern,dumpGeomNode::Up));
		sprintf (fname,"isurfaceU_G%d",i);
		if (i>=constart)  sprintf (fname,"isurfaceU_conflictG%d",i-constart);
		
		dumpGeometry(arcs_tab.begin(),arcs_tab.end(),
			     nodes_tab.begin(),nodes_tab.end(),
			     fname,0,false,false,true);
	      }
	    else
	      {
		nodes_tab.push_back(dumpGeomNode(curn,dumpGeomNode::Up));
		sprintf (fname,"isurfaceU_G%d",i);
		if (i>=constart)  sprintf (fname,"isurfaceU_conflictG%d",i-constart);
	
		dumpGeometry(arcs_tab.begin(),arcs_tab.end(),
			     nodes_tab.begin(),nodes_tab.end(),
			     fname,0,false,false,true);
		nodes_tab.clear();
		nodes_tab.push_back(dumpGeomNode(othern,dumpGeomNode::Down));
		sprintf (fname,"isurfaceD_G%d",i);
		if (i>=constart)  sprintf (fname,"isurfaceD_conflictG%d",i-constart);
		
		dumpGeometry(arcs_tab.begin(),arcs_tab.end(),
			     nodes_tab.begin(),nodes_tab.end(),
			     fname,0,false,false,true);
	      }
	    
	    // sprintf (fname,"isurface_G%d",i);
	    // if (i>=constart)  sprintf (fname,"isurface_conflictG%d",i-constart);
	    // printf ("Dumping impossible surfaces.\n");
	    // dumpGeometry(arcs_tab.begin(),arcs_tab.end(),
	    // 	     nodes_tab.begin(),nodes_tab.end(),
	    // 	     fname);
	    
	  }


	resetIdFromNodeIt();
      }
    
    //printf ("All done.\n");
  }


 
public: 

  NDnetwork *dGradientToNDNet()
  {
    long i,j,k,l;
    long nvert=0;
    std::vector<long> cum;
    for (i=0;i<gradient_pairs.size();i++) 
      {
	nvert+=gradient_pairs[i].size();
	if (i) cum.push_back(gradient_pairs[i].size()+cum[i-1]);
	else cum.push_back(gradient_pairs[i].size());
      }

    int ndims = network->getNDims();
    NDnetwork *net = CreateNetwork(ndims, nvert, true);
    std::vector<double> x0;
    std::vector<double> delta;
    cellType cell;
    

    printf ("Building NDNetwork from discrete gradients ... ");fflush(0);

    double *vtype;
    double *vbound;
    
    network->getBoundingBox(x0,delta);
    
    for (i=0;i<ndims;i++)
      {
	net->x0[i]=x0[i];
	net->delta[i]=delta[i];	
      }
    
    i=0;j=-1;
    net->f_vertexIndex[1]=(NDNET_UINT *)calloc(nvert,sizeof(NDNET_UINT));    
    long curp=0;
    
    vtype = (double *)calloc(net->nvertex,sizeof(double));
    vbound = (double *)calloc(net->nvertex,sizeof(double));

    for (k=0;k<nvert;k++) {
      j++;
      if (j>=gradient_pairs[i].size()) 
	{
	  i++;
	  j=0;
	}
      
      vtype[k]=i;
      cell.set(i,j);
      vbound[k]=(double)network->isBoundary(cell);
	
      
      network->getPosition(cell,&(net->v_coord[ndims*k]));
      cellType val=gradient_pairs[i][j];
      if ((val!=MSC_UNPAIRED)&&(val!=MSC_CRITICAL)&&(val.type()>i))
	{
	  net->f_vertexIndex[1][curp]=k;
	  net->f_vertexIndex[1][curp+1]=cum[i]+val.id();
	  curp+=2;
	}
    }
    addNDDataArr(net,0,TYPE_TAG,&vtype);
    addNDDataArr(net,0,BOUNDARY_TAG,&vbound);

    net->haveVertexFromFace[1]=1;
    net->nfaces[1]=curp/2;
    net->f_vertexIndex[1]=(NDNET_UINT *)realloc(net->f_vertexIndex[1],curp*sizeof(NDNET_UINT));

    net->periodicity=network->getPeriodicity();

    printf ("done. (%ld pairs)\n",(long)net->nfaces[1]);
    return net;
  }

  NDnetwork *pPairsToNDNet()
  {
    //std::vector<node_it> tab;
    //std::vector<node_it> main;
    std::vector<node_it>::iterator tab_it;
    int ndims = network->getNDims();
    int ndims_net = network->getNDims(true);
    long i;
    long nnodes=cplx.num_nodes();
    //std::map<node_it,int,NDcomplex_node::compareItLess> indexof;

    //cplx.computePersistencePairs();
    printf ("Building NDNetwork from persistence pairs ... ");fflush(0);
    int verbSave=verbose;
    verbose=0;
    i=0;
  
    initIdFromNodeIt();

    NDnetwork *net = CreateNetwork(ndims, nnodes, true);
    std::vector<double> x0;
    std::vector<double> delta;

    double *vtype;
    double *ptype;   
    double *vper;
    double *vperR;
    double *pper; 
    double *pperR; 
    double *pv1;
    double *pv2;
    double *vcell;
    double *parent;
    double *parent_log;

    double *vd;
    double *vi;

    bool isPositive=false;
    node_it n;

    for (n=cplx.nodes_begin();n!=cplx.nodes_end();n++)
      if (n->getVal()<=0) break;
    if (n==cplx.nodes_end()) isPositive=true;
  

    network->getBoundingBox(x0,delta);

    for (i=0;i<ndims;i++)
      {
	net->x0[i]=x0[i];
	net->delta[i]=delta[i];
      }

    vtype = (double *)calloc(net->nvertex,sizeof(double));
    vper = (double *)calloc(net->nvertex,sizeof(double));
    if (isPositive) vperR = (double *)calloc(net->nvertex,sizeof(double));
    vd = (double *)calloc(net->nvertex,sizeof(double));
    //vi = (double *)calloc(net->nvertex,sizeof(double));
    vcell = (double *)calloc(net->nvertex,sizeof(double));
    parent = (double *)calloc(net->nvertex,sizeof(double));
    if (isPositive) parent_log = (double *)calloc(net->nvertex,sizeof(double));
    i=0;  
   
    for (n=cplx.nodes_begin();n!=cplx.nodes_end();n++,i++)
      {
	int j;
	float p[ndims];
	//node_it n=*tab_it;
	bool self=((n->getPNode()==cplx.nodes_end())||(n->getPNode()==n));

	network->getPosition(n->getCell(),p);
	for (j=0;j<ndims;j++) net->v_coord[ndims*i+j] = p[j];
	vtype[i]=n->getType();
	vcell[i]=n->getCell().getAsDouble();
	vd[i]=n->getVal();

	node_it pnode=cplx.getParentNode(n,false);
	if (pnode==cplx.nodes_end()) parent[i]=-1;
	else if (pnode->isDeleted()) parent[i]=-1;
	else parent[i]=getIdFromNodeIt(pnode);
	if (isPositive)
	  {
	    pnode=cplx.getParentNode(n,true);
	    if (pnode==cplx.nodes_end()) parent_log[i]=-1;
	    else if (pnode->isDeleted()) parent_log[i]=-1;
	    else parent_log[i]=getIdFromNodeIt(pnode);
	  }

	vper[i]=(self)?-1:n->persistence();	
	if (isPositive) vperR[i]=(self)?-1:n->persistenceR();
	//mp.insert(std::make_pair(n->uniqueID(),i));
      }
  
    long nfaces=0;
    for (n=cplx.nodes_begin();n!=cplx.nodes_end();n++)
      {	
	if ((n->getPNode()==cplx.nodes_end())||(n->getPNode()==n))
	  continue;
	nfaces++;
      }
    if (nfaces&1) fprintf(stderr,"ERROR:there is an odd number of vertice in pair !!!\n");
    nfaces/=2;
    
    addNDDataArr(net,0,TYPE_TAG,&vtype);
    addNDDataArr(net,0,CELL_TAG,&vcell);
    addNDDataArr(net,0,VALUE_TAG,&vd);
    addNDDataArr(net,0,PARENT_TAG,&parent);
    if (isPositive) addNDDataArr(net,0,PARENT_LOG_TAG,&parent_log);
    addNDDataArr(net,0,PERSISTENCE_TAG,&vper);  
    if (isPositive) addNDDataArr(net,0,PERSISTENCE_RATIO_TAG,&vperR);

    net->haveVertexFromFace[1]=1;
    net->nfaces[1]=nfaces;
    net->f_vertexIndex[1]=(NDNET_UINT *)calloc(net->nfaces[1]*2,sizeof(NDNET_UINT));

    ptype = (double *)calloc(net->nfaces[1],sizeof(double));
    pper = (double *)calloc(net->nfaces[1],sizeof(double));
    if (isPositive) pperR = (double *)calloc(net->nfaces[1],sizeof(double));

    pv1 = (double *)calloc(net->nfaces[1],sizeof(double));
    pv2 = (double *)calloc(net->nfaces[1],sizeof(double));

    i=0;     
    //for (tab_it=tab.begin();tab_it!=tab.end();tab_it++)
    for (n=cplx.nodes_begin();n!=cplx.nodes_end();n++)
      {
	//node_it n=*tab_it;
	node_it p=n->getPNode();	
	
	if (n->getType()<=p->getType()) continue;

	long idOfN=getIdFromNodeIt(n);
	long idOfP=getIdFromNodeIt(p);
	net->f_vertexIndex[1][2*i] = idOfN;//mp.find(n->uniqueID())->second;
	net->f_vertexIndex[1][2*i+1] = idOfP;//mp.find(p->uniqueID())->second;
	ptype[i]=n->getType();
	//if (p->getType()>ptype[i]) ptype[i]=p->getType();
	pper[i]=n->persistence();
	if (isPositive) 
	  {
	    pperR[i]=n->persistenceR();
	    if (pperR[i]<1) {
	      printf("Nodes %s/%s\n",n->getInfo(true).c_str(),p->getInfo(true).c_str());
	      printf(" density failed : %12.12g %12.12g\n",n->getVal(),p->getVal());
	    }
	  }

	pv1[i]=idOfN;
	pv2[i]=idOfP;
	i++;	
      }
    
    addNDDataArr(net,1,TYPE_TAG,&ptype);
    addNDDataArr(net,1,PERSISTENCE_TAG,&pper);
    if (isPositive) addNDDataArr(net,1,PERSISTENCE_RATIO_TAG,&pperR);
    addNDDataArr(net,1,UP_TAG(INDEX_TAG),&pv1);
    addNDDataArr(net,1,DOWN_TAG(INDEX_TAG),&pv2);
    
    net->periodicity=network->getPeriodicity();
    verbose=verbSave;
    resetIdFromNodeIt();

    printf ("done. (%ld pairs)\n",(long)net->nfaces[1]);
    return net;
  }

  NDskel *toNDskelWithComplement(double robustnessT, int what = TO_NDSKL_UP)
  {
    cplx.restoreRemovedArcsLinks(robustnessT);
    NDskel *skl=toNDskel(what);
    cplx.unrestoreRemovedArcsLinks();
    return skl;
  }
  
  NDskel *toNDskel(int what = TO_NDSKL_UP)
  {
    int ndims = network->getNDims();
    int ndims_net = network->getNDims(true);
    NDskel *skl;
    int i,j,k;
    int dims[ndims];
    typedef std::pair<std::pair<typeT,idT>,int> map_pair;

    initIdFromNodeIt();

    // skeleton tags
    double *nodePersistence=NULL;
    double *nodePersistenceR=NULL;
    double *nodePersistenceR_nsig=NULL;
    double *nodeDensity=NULL,*segDensity1=NULL,*segDensity2=NULL;
    double *nodeLogDensity=NULL,*segLogDensity=NULL;
    double *nodeRobustness=NULL,*segRobustness=NULL;
    double *nodeRobustnessR=NULL,*segRobustnessR=NULL;
    double *segOrient=NULL;
    double *nodeCell=NULL,*segCell1=NULL,*segCell2=NULL;
    double *segType=NULL;
    double *ppair=NULL;
    double *nodeParent=NULL;
    double *nodeParent_log=NULL;
    bool computeRobustness=false;

    printf ("Building NDskeleton from NDcomplex ... ");fflush(0);
    bool isPositive=false;
    node_it nd;
    for (nd=cplx.nodes_begin();nd!=cplx.nodes_end();nd++)
      if (nd->getVal()<=0) break;
    if (nd==cplx.nodes_end()) isPositive=true;

    for (i=0;i<ndims;i++) dims[i]=0;

    if (cplx.getCyclesState()[1]) computeRobustness=true;

    std::vector<double> x0;
    std::vector<double> delta;

    network->getBoundingBox(x0,delta);

    skl = Create_NDskel(dims,ndims,&(x0[0]),&(delta[0]),"No comments",0,0);

    skl->nnodes=cplx.num_nodes();
    
    float *nodepos=(float*)malloc(skl->nnodes*ndims*sizeof(float));
    skl->nodepos=nodepos;
    skl->Node=(NDskl_node*)calloc(skl->nnodes,sizeof(NDskl_node));

    // First sets node position and tags
    if (computeRobustness) nodeRobustness =(double *)calloc(skl->nnodes,sizeof(double));
    nodePersistence =(double *)calloc(skl->nnodes,sizeof(double));
    if (isPositive) 
      {
	nodePersistenceR =(double *)calloc(skl->nnodes,sizeof(double));
	if (computeRobustness) nodeRobustnessR =(double *)calloc(skl->nnodes,sizeof(double));
	nodePersistenceR_nsig =(double *)calloc(skl->nnodes,sizeof(double));
      }
    nodeDensity =(double *)calloc(skl->nnodes,sizeof(double));
    nodeLogDensity =(double *)calloc(skl->nnodes,sizeof(double));
    //nodeType =(double *)calloc(skl->nnodes,sizeof(double));
    nodeCell =(double *)calloc(skl->nnodes,sizeof(double));
    ppair =(double *)calloc(skl->nnodes,sizeof(double));
    nodeParent =(double *)calloc(skl->nnodes,sizeof(double));
    if (isPositive) nodeParent_log =(double *)calloc(skl->nnodes,sizeof(double));
    i=0;
    for (node_it node= cplx.nodes_begin();node!=cplx.nodes_end();node++,i++)
      {
	//if (node->isInfinite()) continue;
	//std::vector<uint> seg_cells;
	node_it pnode = node->getPNode();
	int ptype = std::min(pnode->getType(),node->getType());
	//if (pnode==node) ptype=-1;
	  
	skl->Node[i].pos = &skl->nodepos[skl->ndims*i];
	network->getPosition(node->getCell(),skl->Node[i].pos);

 	
	skl->Node[i].nnext = 0;
	skl->Node[i].type = node->getType();

	skl->Node[i].flags=0;
	if (node->isInfinite()) skl->Node[i].flags |= FLAG_NDNODE_INFINITE;
	if (node->isOut()) skl->Node[i].flags |= FLAG_NDNODE_OUT;
	if (node->isBoundary()) skl->Node[i].flags |= FLAG_NDNODE_BOUNDARY;
		
	nodeDensity[i]=getValue(node);
	nodeLogDensity[i]=(nodeDensity[i]>0)?log10(nodeDensity[i]):0;

	nodeCell[i]=node->getCell().getAsDouble();

	
	nodePersistence[i]=node->persistence();
	if (isPositive) 
	  {
	    if (node==pnode)
	      nodePersistenceR[i]=nodePersistenceR_nsig[i]=-1;
	    else
	      {
		nodePersistenceR[i]=node->persistenceR();
		nodePersistenceR_nsig[i]=sigmaFromPersistenceRatio(node->persistenceR(),ndims,ptype);
	      }
	  }
	//skl->Node[i].nnext = node->countArcs(node,cplx.arcs_end());
      }
    i=0;
    for (node_it node= cplx.nodes_begin();node!=cplx.nodes_end();node++,i++)
      {
	node_it pnode = node->getPNode();

	ppair[i]=getIdFromNodeIt(pnode);

	pnode=cplx.getParentNode(node,false);
	if (pnode==cplx.nodes_end()) nodeParent[i]=-1;
	else if (pnode->isDeleted()) nodeParent[i]=-1;
	else nodeParent[i]=getIdFromNodeIt(pnode);
	if (isPositive)
	  {
	    pnode=cplx.getParentNode(node,true);
	    if (pnode==cplx.nodes_end()) nodeParent_log[i]=-1;
	    else if (pnode->isDeleted()) nodeParent_log[i]=-1;
	    else nodeParent_log[i]=getIdFromNodeIt(pnode);
	  }

      }
  
    if (isPositive) {
      AddDataField(skl,0,nodePersistenceR,PERSISTENCE_RATIO_TAG);free(nodePersistenceR);
      AddDataField(skl,0,nodePersistenceR_nsig,PERSISTENCE_NSIG_TAG);free(nodePersistenceR_nsig);
    }
    AddDataField(skl,0,nodePersistence,PERSISTENCE_TAG);free(nodePersistence);

    if ((isPositive)&&(computeRobustness)) {
      AddDataField(skl,0,nodeRobustnessR,ROBUSTNESS_RATIO_TAG);free(nodeRobustnessR);
    }
    if (computeRobustness) {
      AddDataField(skl,0,nodeRobustness,ROBUSTNESS_TAG);free(nodeRobustness);
    }

    AddDataField(skl,0,ppair,PERSISTENCE_PAIR_TAG);free(ppair);
    AddDataField(skl,0,nodeParent,PARENT_TAG);free(nodeParent);
    if (isPositive) {AddDataField(skl,0,nodeParent_log,PARENT_LOG_TAG);free(nodeParent_log);}
    AddDataField(skl,0,nodeLogDensity,LOG_TAG(VALUE_TAG));free(nodeLogDensity);
    AddDataField(skl,0,nodeDensity,VALUE_TAG);free(nodeDensity);   
    //AddDataField(skl,0,nodeType,TYPE_TAG);free(nodeType);
    AddDataField(skl,0,nodeCell,CELL_TAG);free(nodeCell);     

    // Compute how many links to other nodes for each node
    i=0;
    for (node_it node= cplx.nodes_begin();node!=cplx.nodes_end();node++,i++)
      {
	for (arc_it arc=node->getArc();arc!=cplx.arcs_end();arc=arc->getNext(node))
	  {
	    if (arc->endpoints().first != node) continue;
	 

	    node_it other= arc->getOtherNode(node);

	    if (!(what & TO_NDSKL_KEEP_STRUCTURE))
	      {
		if ((other->getType() == ndims_net)||(node->getType() == ndims_net))
		  {if (!(what & TO_NDSKL_UP)) continue;}
		else if ((other->getType() == 0)||(node->getType() == 0))
		  {if (!(what & TO_NDSKL_DOWN)) continue;}
		else if (!(what & TO_NDSKL_INTER)) continue;

		if (arc->getGeometry()==cplx.arcsGeom_end()) continue;
	      }
	
	    k=getIdFromNodeIt(other);
	    
	    skl->Node[i].nnext++;
	    skl->Node[k].nnext++;
	  }
      }
	
    // allocate memory for nodes
    i=0;
    for (node_it node= cplx.nodes_begin();node!=cplx.nodes_end();node++,i++)
      {
	//if (node->isInfinite()) continue;
	skl->Node[i].Next  = (NDskl_node**)calloc(skl->Node[i].nnext,sizeof(NDskl_node*));
	skl->Node[i].Seg  = (NDskl_seg**)calloc(skl->Node[i].nnext,sizeof(NDskl_seg*));
	skl->Node[i].nsegs  = (int*)calloc(skl->Node[i].nnext,sizeof(int));
	skl->Node[i].nnext = 0;
      }

    // Set the segments coordinates
    i=0;
    skl->nsegs=0;

    float *segpos=NULL;

    for (node_it node= cplx.nodes_begin();node!=cplx.nodes_end();node++,i++)
      {
	//if (node->isInfinite()) continue;
	std::vector<cellType> seg_cells;
		
	for (arc_it arc=node->getArc();arc!=cplx.arcs_end();arc=arc->getNext(node))
	  {
	    if (arc->endpoints().first != node) continue;
	   
	    node_it other= arc->getOtherNode(node);
	    arcGeom_it ag=cplx.arcsGeom_end();

	    if (!(what & TO_NDSKL_KEEP_STRUCTURE))
	      {
		if ((other->getType() == ndims_net)||(node->getType() == ndims_net))
		  {if (!(what & TO_NDSKL_UP)) continue;}
		else if ((other->getType() == 0)||(node->getType() == 0))
		  {if (!(what & TO_NDSKL_DOWN)) continue;}
		else if (!(what & TO_NDSKL_INTER)) continue;

		ag = arc->getGeometry();
		if (ag==cplx.arcsGeom_end()) 
		  {
		    //fprintf(stderr,"WARNING: cannot retrieve arc geometry.\n");
		    continue;
		  }
	      }
	    else
	      {
		if ((other->getType() == ndims_net)||(node->getType() == ndims_net))
		  {if (what & TO_NDSKL_UP) ag = arc->getGeometry();}
		else if ((other->getType() == 0)||(node->getType() == 0))
		  {if (what & TO_NDSKL_DOWN) ag = arc->getGeometry();}
		else if (what & TO_NDSKL_INTER) ag = arc->getGeometry();
	      }

	    // arcGeom_it ag = arc->getGeometry();
	    // if (ag==cplx.arcsGeom_end()) continue;

	    
	    k=getIdFromNodeIt(other);
	    assert(k!=-1);

	    skl->Node[i].Next[skl->Node[i].nnext] = &(skl->Node[k]);
	    skl->Node[k].Next[skl->Node[k].nnext] = &(skl->Node[i]);

	    long addr=2*skl->ndims*skl->nsegs;

	    seg_cells.clear();
	    if (ag!=cplx.arcsGeom_end()) 
	      ag->get(back_inserter(seg_cells));

	    if (seg_cells.size()==0)
	      {
		seg_cells.push_back(node->getCell());
	      }

	    int nnew=0;
	    // we allocate seg_cells.size() new segments, which is most probably higher than necessary
	    segpos=(float*)realloc(segpos,2*sizeof(float)*skl->ndims*(skl->nsegs+seg_cells.size()));
			
	    //also allocate memory for tags
	    segDensity1 =(double *)realloc(segDensity1,(skl->nsegs+seg_cells.size())*sizeof(double));
	    segDensity2 =(double *)realloc(segDensity2,(skl->nsegs+seg_cells.size())*sizeof(double));
	    segLogDensity =(double *)realloc(segLogDensity,(skl->nsegs+seg_cells.size())*sizeof(double));
	    segType =(double *)realloc(segType,(skl->nsegs+seg_cells.size())*sizeof(double));
	    segCell1 = (double *)realloc(segCell1,(skl->nsegs+seg_cells.size())*sizeof(double));
	    segCell2 = (double *)realloc(segCell2,(skl->nsegs+seg_cells.size())*sizeof(double));
	    segOrient  = (double *)realloc(segOrient,(skl->nsegs+seg_cells.size())*sizeof(double));
	    if (computeRobustness) 
	      segRobustness  = (double *)realloc(segRobustness,(skl->nsegs+seg_cells.size())*sizeof(double));
	    if (isPositive)
	      {
		if (computeRobustness) 
		  segRobustnessR  = (double *)realloc(segRobustnessR,(skl->nsegs+seg_cells.size())*sizeof(double));
	      }
	    
	    std::pair<node_it,node_it> refNodes(std::make_pair(cplx.nodes_end(),cplx.nodes_end()));
	    if (computeRobustness)
	      refNodes=cplx.getReferenceNodes(arc);
	   

	    for (j=1;j<seg_cells.size();j++)
	      {	    
		if (seg_cells[j-1] != seg_cells[j])
		  {
		    network->getPosition(seg_cells[j-1],&segpos[addr]);
		    network->getPosition(seg_cells[j],&segpos[addr+skl->ndims]);
					
		    int l=0;
		    for (l=0;l<skl->ndims;l++) 
		      if (fabs(segpos[addr+skl->ndims+l]-segpos[addr+l])>0.5*skl->delta[l])
			{
			  if (segpos[addr+skl->ndims+l]<segpos[addr+l]) 
			    segpos[addr+skl->ndims+l]+=skl->delta[l];
			  else segpos[addr+l]+=skl->delta[l];
			}
				
		    segCell1[skl->nsegs]=seg_cells[j-1].getAsDouble();
		    segCell2[skl->nsegs]=seg_cells[j].getAsDouble();
		    segDensity1[skl->nsegs]=network->getValue(seg_cells[j-1]);
		    segDensity2[skl->nsegs]=network->getValue(seg_cells[j]);
		    segOrient[skl->nsegs]=(node->getType()<other->getType())?1:-1;
		    double tmpval = 0.5*(segDensity2[skl->nsegs]+segDensity1[skl->nsegs]);
		    segLogDensity[skl->nsegs]=(tmpval>0)?log10(tmpval):0;
		    segType[skl->nsegs]=std::min(node->getType(),arc->getOtherNode(node)->getType());
		    
		    if (computeRobustness) 
		      {
			segRobustness[skl->nsegs]=-1;
			if (isPositive) segRobustnessR[skl->nsegs]=-1;

			if (refNodes.first!=cplx.nodes_end())
			  {
			    segRobustness[skl->nsegs]=fabs(refNodes.first->getVal()-tmpval);
			    if (isPositive) 
			      segRobustnessR[skl->nsegs]=
				(refNodes.first->getVal()>tmpval)?
				(refNodes.first->getVal()/tmpval):
				(tmpval/refNodes.first->getVal());
			  }
			else
			  {
			    segRobustness[skl->nsegs]=-1;
			    if (isPositive) segRobustnessR[skl->nsegs]=-1;
			  }
			
			if (refNodes.second!=cplx.nodes_end())
			  {
			    if ((segRobustness[skl->nsegs]<0)|| 
				(fabs(refNodes.second->getVal()-tmpval)<segRobustness[skl->nsegs]))
			      segRobustness[skl->nsegs]=fabs(refNodes.second->getVal()-tmpval);
			    if (isPositive) 
			      {
				double tmpr=(refNodes.second->getVal()>tmpval)?(refNodes.second->getVal()/tmpval):(tmpval/refNodes.second->getVal());
				if ((segRobustnessR[skl->nsegs]<0)|| 
				    (tmpr<segRobustnessR[skl->nsegs]))
				  segRobustnessR[skl->nsegs]=tmpr;
			      }
			  }	
			
			if (segRobustness[skl->nsegs]==-1) 
			  segRobustness[skl->nsegs]=DBL_MAX;
			if ((isPositive)&&(segRobustnessR[skl->nsegs]==-1))
			  segRobustnessR[skl->nsegs]=DBL_MAX;
		      }
		    
		    skl->nsegs++;
		    nnew++;
		    addr+=2*skl->ndims;
		  }
	      }

	    if (seg_cells.size() !=0)
	      {

		network->getPosition(seg_cells[j-1],&segpos[addr]);
		network->getPosition(other->getCell(),&segpos[addr+skl->ndims]);
		int l=0;
		for (l=0;l<skl->ndims;l++) 
		  if (fabs(segpos[addr+skl->ndims+l]-segpos[addr+l])>0.5*skl->delta[l])
		    {
		      if (segpos[addr+skl->ndims+l]<segpos[addr+l]) segpos[addr+skl->ndims+l]+=skl->delta[l];
		      else segpos[addr+l]+=skl->delta[l];
		    }
		segCell1[skl->nsegs]=seg_cells[j-1].getAsDouble();
		segCell2[skl->nsegs]=other->getCell().getAsDouble();
		segDensity1[skl->nsegs]=network->getValue(seg_cells[j-1]);
		segDensity2[skl->nsegs]=network->getValue(other->getCell());
		segOrient[skl->nsegs]=(node->getType()<other->getType())?1:-1;
		double tmpval = 0.5*(segDensity2[skl->nsegs]+segDensity1[skl->nsegs]);
		segLogDensity[skl->nsegs]=(tmpval>0)?log10(tmpval):0;
		segType[skl->nsegs]=std::min(node->getType(),arc->getOtherNode(node)->getType());
	
		if (computeRobustness) 
		  {
		    segRobustness[skl->nsegs]=-1;
		    if (isPositive) segRobustnessR[skl->nsegs]=-1;
		    
		    if (refNodes.first!=cplx.nodes_end())
		      {
			segRobustness[skl->nsegs]=fabs(refNodes.first->getVal()-tmpval);
			if (isPositive) 
			  segRobustnessR[skl->nsegs]=
			    (refNodes.first->getVal()>tmpval)?
			    (refNodes.first->getVal()/tmpval):
			    (tmpval/refNodes.first->getVal());
		      }
		    
		    if (refNodes.second!=cplx.nodes_end())
		      {
			if ((segRobustness[skl->nsegs]<0)|| 
			    (fabs(refNodes.second->getVal()-tmpval)<segRobustness[skl->nsegs]))
			  segRobustness[skl->nsegs]=fabs(refNodes.second->getVal()-tmpval);
			if (isPositive) 
			  {
			    double tmpr=(refNodes.second->getVal()>tmpval)?(refNodes.second->getVal()/tmpval):(tmpval/refNodes.second->getVal());
			    if ((segRobustnessR[skl->nsegs]<0)|| 
				(tmpr<segRobustnessR[skl->nsegs]))
			      segRobustnessR[skl->nsegs]=tmpr;
			  }
		      }
		    
		    if (segRobustness[skl->nsegs]==-1) 
		      segRobustness[skl->nsegs]=DBL_MAX;
		    if ((isPositive)&&(segRobustnessR[skl->nsegs]==-1))
		      segRobustnessR[skl->nsegs]=DBL_MAX;
		  }

		skl->nsegs++;
		nnew++;
	      }
	    skl->Node[i].nsegs[skl->Node[i].nnext++]=nnew;
	    skl->Node[k].nsegs[skl->Node[k].nnext++]=nnew;
	  }
      }

    // setup Seg
    skl->segpos=segpos;
    skl->Seg=(NDskl_seg*)calloc(skl->nsegs,sizeof(NDskl_seg));
    for (i=0;i<skl->nsegs;i++)
      {
	skl->Seg[i].pos = &(skl->segpos[2*i*skl->ndims]);
	skl->Seg[i].index=i;
	if (i!=skl->nsegs-1) skl->Seg[i].Next=&skl->Seg[i+1];
	if (i!=0) skl->Seg[i].Prev=&skl->Seg[i-1];
      }

    for (i=0;i<skl->nnodes;i++) skl->Node[i].nnext = 0;
    for (i=0;i<skl->nnodes;i++) skl->Node[i].index = i;
    for (i=0;i<skl->nsegs;i++) skl->Seg[i].index = i;

    AddDataField(skl,1,segDensity1,SEG_P1_TAG(VALUE_TAG));free(segDensity1);
    AddDataField(skl,1,segDensity2,SEG_P2_TAG(VALUE_TAG));free(segDensity2);
    AddDataField(skl,1,segOrient,ORIENTATION_TAG);free(segOrient);
    AddDataField(skl,1,segCell1,SEG_P1_TAG(CELL_TAG));free(segCell1);
    AddDataField(skl,1,segCell2,SEG_P2_TAG(CELL_TAG));free(segCell2);
    AddDataField(skl,1,segLogDensity,LOG_TAG(VALUE_TAG));free(segLogDensity); 
    AddDataField(skl,1,segType,TYPE_TAG);free(segType);
    if (computeRobustness) 
      {
	AddDataField(skl,1,segRobustness,ROBUSTNESS_TAG);free(segRobustness);    
	if (isPositive) {AddDataField(skl,1,segRobustnessR,ROBUSTNESS_RATIO_TAG);free(segRobustnessR);}
      }
 
    // Set the segments paths between nodes    
    i=0;
    skl->nsegs=0;
    for (node_it node= cplx.nodes_begin();node!=cplx.nodes_end();node++,i++)
      {
	//if (node->isInfinite()) continue;
	for (arc_it arc=node->getArc();arc!=cplx.arcs_end();arc=arc->getNext(node))
	  {
	    if (arc->endpoints().first != node) continue;
	  
	      
	    node_it other= arc->getOtherNode(node);

	    if (!(what & TO_NDSKL_KEEP_STRUCTURE))
	      {
		if ((other->getType() == ndims_net)||(node->getType() == ndims_net))
		  {if (!(what & TO_NDSKL_UP)) continue;}
		else if ((other->getType() == 0)||(node->getType() == 0))
		  {if (!(what & TO_NDSKL_DOWN)) continue;}
		else {if (!(what & TO_NDSKL_INTER)) continue;}

		if (arc->getGeometry()==cplx.arcsGeom_end()) continue;
	      }
	   
	    k=getIdFromNodeIt(other);

	    int i1=skl->Node[i].nnext++;
	    int i2=skl->Node[k].nnext++;
		
	    skl->Node[i].Seg[i1] = &skl->Seg[skl->nsegs];
	    skl->Node[k].Seg[i2] = &skl->Seg[skl->nsegs+skl->Node[i].nsegs[i1]-1];
	    assert(skl->Node[i].nsegs[i1] == skl->Node[k].nsegs[i2]);
	    assert(skl->Node[i].Next[i1] == &skl->Node[k]);
	    assert(skl->Node[k].Next[i2] == &skl->Node[i]);

	    for (j=skl->nsegs;j<skl->nsegs+skl->Node[i].nsegs[i1];j++)
	      {
		skl->Seg[j].nodes[0]=i;
		skl->Seg[j].nodes[1]=k;
	      }
	    skl->Seg[skl->nsegs].Prev=NULL;
	    skl->Seg[skl->nsegs+skl->Node[i].nsegs[i1]-1].Next=NULL;
	    skl->nsegs+=skl->Node[i].nsegs[i1];
	  }
      }

    for (i=0;i<skl->nnodes;i++) skl->Node[i].data = &skl->nodedata[i*skl->nnodedata];
    for (i=0;i<skl->nsegs;i++) skl->Seg[i].data = &skl->segdata[i*skl->nsegdata];

    // node robustness
    if (computeRobustness) 
      {
	int robSegId = getDataFieldID(skl,1,ROBUSTNESS_TAG);
	int robNodeId = getDataFieldID(skl,0,ROBUSTNESS_TAG);
	int robRSegId = getDataFieldID(skl,1,ROBUSTNESS_RATIO_TAG);
	int robRNodeId = getDataFieldID(skl,0,ROBUSTNESS_RATIO_TAG);

	for (i=0;i<skl->nnodes;i++)	
	  {
	    skl->Node[i].data[robNodeId]=0;
	    if (isPositive) skl->Node[i].data[robRNodeId]=0;
	  }

	for (i=0;i<skl->nnodes;i++)
	  {
	    NDskl_node *node=&skl->Node[i];
	
	    for (j=0;j<node->nnext;j++)
	      {
		NDskl_seg *seg=node->Seg[j];
	   
		if (node->data[robNodeId]<seg->data[robSegId])
		  node->data[robNodeId]=seg->data[robSegId];

		if (isPositive)
		  {
		    if (node->data[robRNodeId]<seg->data[robRSegId])
		      node->data[robRNodeId]=seg->data[robRSegId];
		  }
	      }
	    //printf(" ->%g\n",node->data[robNodeId]);
	  }
      }
    resetIdFromNodeIt();
    printf ("done. (%d nodes / %d arcs = %d segs)\n",cplx.num_nodes(),cplx.num_arcs(),skl->nsegs);
    
    return skl;
  }

  
  void dumpSkl(const char *fname_p=NULL, int what=(1<<0),const char *str=NULL, int smooth=0)
  {
    NDskel *skl;
    char fname[255];
    int flags=0;

    if (what&(1<<0)) flags|=TO_NDSKL_UP;
    if (what&(1<<1)) flags|=TO_NDSKL_INTER;
    if (what&(1<<2)) flags|=TO_NDSKL_DOWN;    
    if (what&(1<<3)) flags|=TO_NDSKL_KEEP_STRUCTURE;

    //if (!(flags&TO_NDSKL_KEEP_STRUCTURE))
      {
	if ((flags&TO_NDSKL_UP)&&(!(store_arcsGeom&(1<<0)))) 
	  printf("WARNING: up arcs are not defined ! (do not load from this MSC)\n");
	if ((flags&TO_NDSKL_DOWN)&&(!(store_arcsGeom&(1<<2))))
	  printf("WARNING: down arcs are not defined ! (do not load from this MSC)\n");
	if ((flags&TO_NDSKL_INTER)&&(!(store_arcsGeom&(1<<1))))
	  {
	    printf("WARNING: inter arcs are not defined ! (do not load from this MSC)\n");
	    printf("   Option -interArcsGeom will force computing their geometry.\n");
	  }
      }

    if (fname_p==NULL) strcpy(fname,"skl.");
    else sprintf(fname,"%s.",fname_p);    
    if (str==NULL)
      {
	if (flags&TO_NDSKL_UP) sprintf(fname,"%su",fname);
	if (flags&TO_NDSKL_INTER) sprintf(fname,"%si",fname);
	if (flags&TO_NDSKL_DOWN) sprintf(fname,"%sd",fname);
	if (flags&TO_NDSKL_KEEP_STRUCTURE) sprintf(fname,"%sc",fname);
      }
    else strcat(fname,str);
    //strcat(fname,".NDskl");

    skl = toNDskel(flags);
    //skl = toNDskelWithComplement(0,flags);
    //Save_NDskel(skl,fname);
    ndskel::IO::save(skl,std::string(fname)+ndskel::IO::getExtension());
    Free_NDskeleton(&skl);    
  }

  
  void dumpPPairs(const char *fname_p=NULL, bool toNDnet = true)
  {
    if (toNDnet)
      {
	NDnetwork *net= pPairsToNDNet();
	char fname[255];
	if (net==NULL) return;

	if (fname_p==NULL)
	  strcpy(fname,"ppairs.NDnet");
	else
	  sprintf (fname,"%s.ppairs.NDnet",fname_p);
	
	ndnet::IO::save(net,std::string(fname));
	FreeNDnetwork(&net);
      }
    else
      {
	int ndims = network->getNDims();
	char fname[255];
	long i;
	if (fname_p==NULL)
	  strcpy(fname,"ppairs.ASCIIpairs");
	else
	  sprintf(fname,"%s.ASCIIpairs",fname_p);

	 printf ("Saving persistence pairs to ASCII file %s ...",fname);fflush(0);
	 std::vector<long> count(ndims,0);
	 for (node_it n= cplx.nodes_begin(); n!=cplx.nodes_end();n++)
	   {
	     node_it p=n->getPNode();
	     if (p==n) continue;
	     if (p->getType()<n->getType()) continue;
	     count[n->getType()]++;
	   }
	 FILE *f;
	 f=fopen(fname,"w");
	 fprintf(f,"PERSISTENCE PAIRS\n");
	 fprintf(f,"%d\n",ndims);
	 for (i=0;i<ndims;i++) {
	   fprintf(f,"%ld %ld\n",count[i],i);
	   for (node_it n= cplx.nodes_begin(); n!=cplx.nodes_end();n++)
	     {
	       node_it p=n->getPNode();
	       if (p==n) continue;
	       if (p->getType()<n->getType()) continue;
	       if (n->getType()!=i) continue;
	       
	       fprintf(f,"%g %g\n",n->getVal(),p->getVal());
	     }
	 }
	 fclose(f);
	 printf("done.\n");
      }
  }

  void dumpDGradient(const char *fname_p=NULL, bool toNDnet = false)
  {
    char fname[255];

    if (fname_p==NULL)
      strcpy(fname,"DGradient");
    else
      sprintf (fname,"%s",fname_p);

    setDiscreteGradientState(DG_saveState);
    setDiscreteGradientState(DG_available);

    if (toNDnet) 
      {
	sprintf(fname,"%s.NDnet",fname);
	printf ("Dumping discrete gradient (NDnet) ... ");fflush(0);
	NDnetwork *net= dGradientToNDNet();
	ndnet::IO::save(net,std::string(fname));
	FreeNDnetwork(&net);
      }
    else
      {
	printf ("Dumping discrete gradient (raw) ... ");fflush(0);
	FILE *f=fopen(fname,"w");
	long j;
	for (int i=0;i<gradient_pairs.size();i++)
	  {
	    long N=gradient_pairs[i].size();
	    cellType::storageType tmp;
	    fwrite(&N,sizeof(int),1,f);	    
	    for (j=0;j<N;j++) {
	      tmp=gradient_pairs[i][j]();
	      fwrite(&tmp,sizeof(cellType::storageType),1,f);
	    }
	  }
	
	fclose(f);
	printf ("done.\n");
	printf ("Discrete gradient was dumped to file '%s'.\n",fname);
      }
  
    setDiscreteGradientState(DG_loadState);

  }


private:

  friend class dumpGeomNode;
  class dumpGeomNode { 
  public:
    enum dumpGeomNodeType {Up=(1<<0), Down=(1<<1), Cycle=(1<<2), All=((1<<3)-1)};
    node_it node;  
    dumpGeomNodeType dumpType;

    dumpGeomNode(node_it n, dumpGeomNodeType d)
    {
      node=n;dumpType=d;
    }    
  };
  

  template <class aInputIterator, class nInputIterator>
  void dumpGeometry(aInputIterator arc_begin,aInputIterator arc_end,
		    nInputIterator node_begin,nInputIterator node_end,
		    const char* fname_p, int smooth=0, bool preserve=false, bool extended=false, bool dual=false)
  {
    char fname[255];
    int ndims=network->getNDims();   
    std::vector< std::pair<long,long> > count(ndims+1,std::make_pair(0,0));
   
    if (fname_p==NULL)
      strcpy(fname,"geometry");
    else
      sprintf (fname,"%s",fname_p);
    //strcpy(fname,fname_p);

    printf("Dumping geometry ");fflush(0);

    //cplx.computePersistencePairs();      

    subCellsType cells;
    cellsInfoType cellsInfo;
    std::vector<double> &cellsInfoData=cellsInfo.second;
    newCellsType newCells;
    cellsInfoType newCellsInfo;
    std::vector<double> &newCellsInfoData=cellsInfo.second;
    newCellsElementType newSeg(2);

    bool tagSourceCellID=IdFromNodeItIsAvailable();

    //cellsInfo.first.push_back(std::string("ascending origin"));
    //cellsInfo.first.push_back(std::string("descending origin"));

    if (tagSourceCellID) 
      cellsInfo.first.push_back(std::string(SOURCE_TAG(INDEX_TAG)));
    else
      cellsInfo.first.push_back(std::string(SOURCE_TAG(CELL_TAG)));

    for (aInputIterator it = arc_begin;it!=arc_end;it++)
      {
	arc_it arc = *it;
	std::pair<NDcomplex_node::node_it,NDcomplex_node::node_it> endpts = arc->endpoints();
	cells.push_back(endpts.first->getCell());
	if (tagSourceCellID) 
	  cellsInfoData.push_back(getIdFromNodeIt(endpts.first));
	else 
	  cellsInfoData.push_back(endpts.first->getCell().getAsDouble());

	cells.push_back(endpts.second->getCell());
	if (tagSourceCellID) 
	  cellsInfoData.push_back(getIdFromNodeIt(endpts.second));
	else 
	  cellsInfoData.push_back(endpts.second->getCell().getAsDouble());
  
	arcGeom_it g=arc->getGeometry();
	if (g==cplx.arcsGeom_end()) continue;

	//char type = g->getType();
	std::vector< cellType > ids;

	g->get(back_inserter(ids));
	if (!ids.size()) continue;
	for(int i=0;i<(int)ids.size()-1;i++)
	  {
	    newSeg.front()=ids[i];
	    newSeg.back()=ids[i+1];
	    newCells.push_back(newSeg);
	  }
   
	newSeg.front()=ids.back();
	if (newSeg.front().type()==endpts.first->getCell().type())
	  newSeg.back()=endpts.second->getCell();
	else
	  newSeg.back()=endpts.first->getCell();

	newCells.push_back(newSeg);
      }

    for (nInputIterator it = node_begin;it!=node_end;it++)
      {
	node_it node=(*it).node;

	if ((*it).dumpType & dumpGeomNode::Up) {
	  long ii=cells.size();
	  count[node->getType()].first++;
	  getNodeManifold(node,back_inserter(cells),true,extended,false);

	  double did;
	  if (tagSourceCellID) 
	    did=getIdFromNodeIt(node);
	  else
	    did=node->getCell().getAsDouble();

	  for (long i=ii;i<cells.size();i++)
	    {
	      cellsInfoData.push_back(did);
	   
	    }
	    
	  std::set<cellType> bout;
	  if (vertexAsMaxima)
	    {
	      /*
	      if (0) {
		getBoundary(cells.begin()+ii, cells.end(),bout,!true);
		cells.insert(cells.end(),bout.begin(),bout.end());
	      }
	      */
	      //std::copy(bout.begin(),bout.end(),cells.end());
	    }
	  else
	    {
	      // represent it with the voronoi dual cells
	      // if appropriate ...
	      /*
	      if (0) {
		getBoundary(cells.begin()+ii, cells.end(),bout,true);
		cells.insert(cells.end(),bout.begin(),bout.end());
	      }
	      */
	      //std::copy(bout.begin(),bout.end(),cells.end());
	      if (dual)
		{
		  for (int i=ii;i<cells.size();i++) 
		    cells[i].setType(cells[i].type()+1+ndims);	  
		}
	    }
	}

	if ((*it).dumpType & dumpGeomNode::Down) {
	  long ii=cells.size();
	  count[node->getType()].second++;
	  getNodeManifold(node,back_inserter(cells),false,extended,false);
	  
	  //double did=node->getCell().getAsDouble();
	  double did;
	  if (tagSourceCellID) 
	    did=getIdFromNodeIt(node);
	  else
	    did=node->getCell().getAsDouble();
	  

	  for (long i=ii;i<cells.size();i++)
	    {
	      cellsInfoData.push_back(did);
	    }

	  std::set<cellType> bout;
	  if (!vertexAsMaxima)
	    {
	      /*
	      if (0) {
		getBoundary(cells.begin()+ii, cells.end(),bout,!true);
		cells.insert(cells.end(),bout.begin(),bout.end());
	      }
	      */
	      //std::copy(bout.begin(),bout.end(),cells.end());
	    }
	  else
	    {
	      // represent it with the voronoi dual cells
	      // if appropriate ...
	      /*
	      if (0) {
		getBoundary(cells.begin()+ii, cells.end(),bout,true);
		cells.insert(cells.end(),bout.begin(),bout.end());
	      }
	      */
	      
	      if (dual)
		{
		  for (int i=ii;i<cells.size();i++) 
		    cells[i].setType(cells[i].type()+1+ndims);
		}
	    }
	}
      }
    
    int printed=0;
    printf("(");
    for (long i=0;i<ndims+1;i++)
      {
	if (count[i].first)
	  printf("%s%ld %ld-a",(printed++)?",":"",count[i].first,i);
	if (count[i].second)
	  printf("%s%ld %ld-d",(printed++)?",":"",count[i].second,i);
      }
    printf(") ...\n");

    network->dumpSubComplex(fname, cells ,cellsInfo, newCells, newCellsInfo, smooth, preserve);
    
  }

public:

  void dumpManifolds(const char *fname, int smooth, long which, char *whichSTR)
  {
    std::vector<arc_it> arcs_tab;
    std::vector<dumpGeomNode> nodes_tab;
    int ndims = network->getNDims();
    char newname[255];
    int ct[5] = {0,0,0,0,0};
    long what = which;    

    if (!store_manifolds)
      {
	fprintf(stderr,"ERROR: cannot dump manifolds.\n");
	fprintf(stderr,"   You are probably loading a MSC file for which they were not computed.\n");
	fprintf(stderr,"   Restart 'mse' without loading the MSC file or recompute it with option '-manifolds'.\n");
	return;
      }
     
    initIdFromNodeIt();

    if (what==0)     
      what=((long)1<<(8*sizeof(which)-1))+((long)1<<0)+((long)1<<(2*ndims+1));         

    bool preserve=what&((long)1<<((long)8*sizeof(which)-1));
    bool extended=what&((long)1<<((long)8*sizeof(which)-2));
    bool join=what&((long)1<<((long)8*sizeof(which)-3));
    bool dual=what&((long)1<<((long)8*sizeof(which)-4));

    for (node_it n= cplx.nodes_begin(); n!=cplx.nodes_end();n++)    
      {
	bool skip=true;
	float p[ndims];

	if (!join)
	  {
	    arcs_tab.clear();
	    nodes_tab.clear();
	  }

	if (what&(1<<(2*n->getType()))) 
	  {
	    nodes_tab.push_back(dumpGeomNode(n,dumpGeomNode::Up));
	    skip=false;
	  }
	if (what&(1<<(2*n->getType()+1))) 
	  {
	    nodes_tab.push_back(dumpGeomNode(n,dumpGeomNode::Down));
	    skip=false;
	  }
	if (skip) continue;

	if (!join)
	  {
	    sprintf(newname,"%s_T%d_%d",fname,n->getType(),ct[n->getType()]++);
	    dumpGeometry(arcs_tab.begin(),arcs_tab.end(),
			 nodes_tab.begin(),nodes_tab.end(),
			 newname, smooth, preserve, extended, dual);
	  }
      }

    if (join)
      {
	sprintf(newname,"%s_manifolds_%s",fname,whichSTR);
	dumpGeometry(arcs_tab.begin(),arcs_tab.end(),
		     nodes_tab.begin(),nodes_tab.end(),
		     newname, smooth, preserve, extended, dual);
      }

    resetIdFromNodeIt();
  } 

  void drawDiagram(const char *fname)
  {
    char newname[255];
    sprintf(newname,"%s.diag.dot",fname);
    FILE *f=fopen(newname,"w");
    fprintf(f,"digraph G {");
    int i=0;
    int j=0;
    
    std::vector<std::string> boxtype;
    boxtype.push_back(std::string("box"));
    boxtype.push_back(std::string("circle"));
    boxtype.push_back(std::string("diamond"));
    boxtype.push_back(std::string("triangle"));
    std::vector<std::string> color;
    color.push_back(std::string("blue"));
    color.push_back(std::string("green"));
    color.push_back(std::string("yellow"));
    color.push_back(std::string("red"));
    //std::sort(cplx.nodes_begin(),cplx.nodes_end(),NDcomplex_node::compareLess);
    std::map<double,int> nodeID;
    for (node_it n=cplx.nodes_begin(); n!=cplx.nodes_end();n++,i++)
      {
	fprintf(f,"%d [shape=%s,color=%s,style=%s];\n",i,
		boxtype[n->getType()].c_str(),
		color[n->getType()].c_str(),
		n->isInfinite()?"diagonals":(n->isBoundary()?"diagonals":"filled"));
	
	nodeID.insert(std::pair<double,int>(n->uniqueID(),i));
      }
    
    
    std::vector<node_it> nei;
    std::vector<node_it> id;
    bool linked;
    char txt[255];
    std::set<node_it,NDcomplex_node::compareItLess> uniq;
    std::vector<node_it>::iterator nei_it;
    i=0;
    for (node_it n=cplx.nodes_begin(); n!=cplx.nodes_end();n++,i++)
      {
	nei.clear();
	id.clear();
	
	if (n->getType()==0) continue;
	node_it p =n->getPNode();
	if (n->getType()<p->getType()) p=n;

	n->getNeighbors(n,cplx.arcs_end(),back_inserter(nei),n->getType()-1);
	linked=false;
	for (nei_it=nei.begin();nei_it!=nei.end();nei_it++)
	  {
	    strcpy(txt,"");
	    if (n->persistence(*nei_it)==0) {
	      sprintf(txt,"style=dotted,");
	    
	    }
	    if ((p!=n)&&(p==(*nei_it))) {
	      sprintf(txt,"color=red,");
	      linked=true;
	    }
	   
	    fprintf(f,"%d->%d [%s];\n",i,
		    nodeID[(*nei_it)->uniqueID()],
		    txt);
	   
	  }
	if ((p!=n)&&(!linked)) {
	  fprintf(f,"%d->%d [style=dashed,color=red,arrowhead=none,arrowtail=none]\n",
		  nodeID[n->uniqueID()],
		  nodeID[p->uniqueID()]);
	}

	id.clear();

      }
    fprintf(f,"}");
    fclose(f);

    char cmd[255];
    sprintf(cmd,"dot -Tps %s.diag.dot -o %s.diag.ps; gv %s.diag.ps &",fname,fname,fname);
    int err = system(cmd);
  }


  // implemented functions 

private:
  int filterDGradient(dbVec &pairs)
  {
    int ndims=network->getNDims();
    uint i,j,k;
    std::vector<cellType> cells;
    std::vector< std::vector<cellType>* > skl;
    uint count=0;
    cellType curCell;
    printf("Filtering discrete gradient ... ");fflush(0);
    //int i0,imax;
    bool firstPass=true;
    std::vector<int> indexVal;
    //i0=1;imax=pairs.size()-1;
    
    indexVal.push_back(1);
    if (ndims>2) indexVal.push_back(ndims-1);
    for (i=2;i<ndims-1;i++) indexVal.push_back(i);
    


    //#pragma omp parallel shared(pairs,count) private(i,j,k,cells,skl,curCell)
    //for (i=i0;i<imax;i++)
    for (k=0;k<indexVal.size();k++)
      {
	i=indexVal[k];
	int critIndex = critIndexFromCellType(i);
	long chunk = pairs[i].size()/(glob_num_omp_threads);

	//#pragma omp for schedule(static,chunk) nowait
	for (j=0;j<pairs[i].size();j++)
	  {	  
	    curCell.set(i,j);
	    
	    if (pairs[i][j] == MSC_CRITICAL)
	      {
		
		if (network->isOut(curCell)) continue;

		std::vector< std::vector<cellType>* >::iterator skl_it;
		cells.clear();
		skl.clear();

		
		if ((critIndex==1)&&(firstPass))
		  {
		    ComputeDescendingManifold(curCell,pairs,back_inserter(cells),skl);
		  }
		else
		  {
		    ComputeAscendingManifold(curCell,pairs,back_inserter(cells),skl);
		  }
		
				
		double d=network->getValue(curCell);
		
		for (skl_it= skl.begin();skl_it!=skl.end();skl_it++)
		  {
		    std::vector<cellType> &path = *(*skl_it);
		    double d2=network->getValue(path[0]);		    		  
		    bool cancel=false;		 

		    if (d2==d) 
		      {
			// check boundary (only cancel strictly interior or boundary pairs)
			if (network->isOut(path[0])) continue;
			if (network->isBoundary(curCell)!=network->isBoundary(path[0])) continue;
		
			// check for loops
			std::vector< std::vector<cellType>* >::iterator skl_it2;
			int nfound=0;
		     
			for (skl_it2= skl.begin();skl_it2!=skl.end();skl_it2++) {
			  std::vector<cellType> &other_path = *(*skl_it2);
			  if (other_path[0] == path[0]) nfound++;
			}

			if (nfound==1) cancel=true;
		      
		      }
		  
		    if (cancel) {
		      // 0-persistence pair, cancel the path
#pragma omp critical
#ifdef USE_OPENMP
		      if ((pairs[path[0].type()][path[0].id()]== MSC_CRITICAL)&&
			  (pairs[i][j] == MSC_CRITICAL)) 
#endif
			{
			  count++;
			  
			  int k;
			  cellType nv,new_nv;
			  
			  nv = pairs[path[1].type()][path[1].id()];
			  pairs[path[1].type()][path[1].id()]=path[0];
			  pairs[path[0].type()][path[0].id()]=path[1];

			  for (k=2;k<(int)path.size();k++)
			    {
			      new_nv = pairs[path[k].type()][path[k].id()];
			      pairs[path[k].type()][path[k].id()] = nv;
			      pairs[nv.type()][nv.id()] = path[k];
			      nv=new_nv;
			    }
			  
			  //break;
			}
#ifdef USE_OPENMP
		      else j--;
#endif	
		      break;
		    }

		  } //for loop
				
		for (skl_it = skl.begin();skl_it!=skl.end();skl_it++)
		  delete *skl_it;
	      }
	  }

	if ((critIndex==1)&&(firstPass)) 
	  {
	    k--;
	    firstPass=false;
	  }
      }

    printf("done. (%d arcs removed)\n",count);

    return count;
  }


  template <class OutputIterator>
  std::set<uint>*
  ComputeAscendingManifold(cellType oid,cellType id, const dbVec &pairs,
			   OutputIterator manifold,std::vector< std::vector<cellType>* > &skl)
  {
    std::vector<cellType> coflist;
    std::set<uint> *skl_id=NULL;

    if (id!=oid) *manifold=id;
    network->getCofaces(id,coflist);
    
    for (int i=0;i<coflist.size();i++)
      {
	if (network->isOut(coflist[i])) continue;
	cellType tmp = pairs[coflist[i].type()][coflist[i].id()];
	
	if (tmp == MSC_CRITICAL)
	  {
	    if (id==oid) *manifold=id;
	    skl.push_back(new std::vector<cellType>(1,coflist[i]));
	    if (skl_id==NULL) skl_id=new std::set<uint>();
	    skl_id->insert(skl.size()-1);
	  }
	else
	  if ((tmp.type()==id.type())&&(tmp != id))
	    {
	      if (id==oid) *manifold=id;
	      std::set<uint> *skl_r=ComputeAscendingManifold(oid,tmp,pairs,manifold,skl);
	    
	      if (skl_id==NULL) skl_id=skl_r;
	      else if (skl_r!=NULL) {
		skl_id->insert(skl_r->begin(),skl_r->end());
		delete skl_r;
	      }
	    }
      }
  
    if (skl_id!=NULL)
      {
	for (std::set<uint>::iterator it= skl_id->begin();  it!=skl_id->end(); it++)
	  skl[*it]->push_back(id);
      }
  
    if (id==oid) {
      delete skl_id;
      skl_id=NULL;
    }
  
    return skl_id;
  }

  template <class OutputIterator>
  std::set<uint>*
  ComputeAscendingManifold(cellType id ,const dbVec &pairs,
			    OutputIterator manifold,
			    std::vector< std::vector<cellType>* > &skl)
  {
    if (vertexAsMaxima)
      return ComputeDescendingManifold(id,id,pairs,manifold,skl);
    else
      return ComputeAscendingManifold(id,id,pairs,manifold,skl);
  }

  template <class OutputIterator>
  std::set<uint>*
  ComputeDescendingManifold(cellType oid, cellType id ,const dbVec &pairs,
			    OutputIterator manifold,std::vector< std::vector<cellType>* > &skl, long level=0)
  {
    std::vector<cellType> flist;
    std::set<uint> *skl_id=NULL;

    if (id!=oid) *manifold=id;
    network->getFaces(id,flist);
    
    for (int i=0;i<flist.size();i++)
      {
	if (network->isOut(flist[i])) continue;
	cellType tmp = pairs[flist[i].type()][flist[i].id()];
	
	if (tmp == MSC_CRITICAL)
	  {
	    if (id==oid) *manifold=id;
	    skl.push_back(new std::vector<cellType>(1,flist[i]));
	    if (skl_id==NULL) skl_id=new std::set<uint>();
	    skl_id->insert(skl.size()-1);
	  }
	else
	  if ((tmp.type()==id.type())&&(tmp != id))
	    {
	      if (id==oid) *manifold=id;
	  
	      std::set<uint> *skl_r=ComputeDescendingManifold(oid,tmp,pairs,manifold,skl,level+1);

	      if (skl_id==NULL) skl_id=skl_r;
	      else if (skl_r!=NULL) 
		{
		  skl_id->insert(skl_r->begin(),skl_r->end());
		  delete skl_r;
		}
	    }
      }
  
    if (skl_id!=NULL)
      {
	for (std::set<uint>::iterator it= skl_id->begin();  it!=skl_id->end(); it++)
	  skl[*it]->push_back(id);
      }

    if (id==oid) {
      delete skl_id;
      skl_id=NULL;
    }

    return skl_id;
  }

  template <class OutputIterator>
  std::set<uint>*
  ComputeDescendingManifold(cellType id ,const dbVec &pairs,
			    OutputIterator manifold,
			    std::vector< std::vector<cellType>* > &skl)
  {
    if (vertexAsMaxima)
      return ComputeAscendingManifold(id,id,pairs,manifold,skl);
    else
      return ComputeDescendingManifold(id,id,pairs,manifold,skl);
  }
 
private:
  
  template <class InputIterator>
  void getBoundary(InputIterator begin, InputIterator end, std::set<cellType> &out,bool co=false)
  {
  
    std::vector<cellType> nei;
    
    for (InputIterator it=begin; it!=end;it++)
      {
	cellType c = *it;
	
	if (co) {
	  network->getCofaces(c, nei);
	  for (std::vector<cellType>::iterator f=nei.begin();f!=nei.end();f++)
	    out.insert(*f);
	}
	else {
	  network->getFaces(c, nei);
	  for (std::vector<cellType>::iterator f=nei.begin();f!=nei.end();f++)
	    out.insert(*f);
	}
	
	
      }   
    //return out.begin();
  }

  template <class InputIterator>
  void getBoundary(InputIterator begin, InputIterator end, Z2set<cellType> &out,bool co=false)
  {
  
    std::vector<cellType> nei;
    
    for (InputIterator it=begin; it!=end;it++)
      {
	cellType c = *it;
	
	if (co) {
	  network->getCofaces(c, nei);
	  for (std::vector<cellType>::iterator f=nei.begin();f!=nei.end();f++)
	    out.insert(*f);
	}
	else {
	  network->getFaces(c, nei);
	  for (std::vector<cellType>::iterator f=nei.begin();f!=nei.end();f++)
	    out.insert(*f);
	}
	
	
      }      
  }

  template <class InputIterator, class OutputIterator>
    OutputIterator getBoundary(InputIterator begin, InputIterator end, OutputIterator out,bool co=false)
  {
  
    std::set<cellType> output;
    getBoundary(begin,end,output,co);
    for (typename std::set<cellType>::iterator it=output.begin();it!=output.end();it++)
      *out = *it;
    return out;
  }

  template <class InputIterator, class OutputIterator>
    OutputIterator getBoundary(InputIterator begin, InputIterator end, OutputIterator out, bool strict,bool co)
  {
  
    if (strict) {
      Z2set<cellType> output;  
      getBoundary(begin,end,output,co);
      for (typename Z2set<cellType>::iterator it=output.begin();it!=output.end();it++)
      *out = *it;
     
    } 
    else return getBoundary(begin,end,out,co);
    
    return out;      
  }
  
  public:
  
  template <class OutputIterator>
    void getNodeManifold(node_it node, OutputIterator out, bool up, bool extended=false, bool boundary=false)
  {
    if (extended) { getNodeExtendedManifold(node,out,up); return ;}
    
    nodeGeom_it g;
    std::vector<cellType> tmp;
      
    if (up)
      {
	g=node->getUpGeometry();
	if (g==cplx.nodesGeom_end()) return ;
      }
    else
      {
	g=node->getDownGeometry();
	if (g==cplx.nodesGeom_end()) return ;
      }
      
    g->get(back_inserter(tmp));
   

    for(std::vector< cellType >::iterator id=tmp.begin();id!=tmp.end();id++)
      *out = *id;
  
  }

    
  
  template <class OutputIterator>
  void getNodeExtendedManifold(node_it node, OutputIterator out, bool up)
  {
    std::set<cellType> output;
    getNodeExtendedManifold(node, output,up);
    for (typename std::set<cellType>::iterator it=output.begin();it!=output.end();it++)
      *out = *it;
  
  }

  void getNodeExtendedManifold(node_it node, std::set<cellType> &out, bool up)
  {
    std::vector<cellType> tmp;
    nodeGeom_it g;
  
    if (up)
      {
	g=node->getUpGeometry();
	if (g==cplx.nodesGeom_end()) return;// out.begin();
      }
    else
      {
	g=node->getDownGeometry();
	if (g==cplx.nodesGeom_end()) printf("ERROR:%s\n",node->getInfo(true).c_str());
	if (g==cplx.nodesGeom_end()) return;// out.begin();
      }

    g->get(back_inserter(tmp));
   
    //Compute the closure
    out.insert(tmp.begin(),tmp.end());

   
    int ntype = node->getType();
   
    if (up)
      {
	if (ntype<network->getNDims()-1)
	  {
	    for (arc_it arc=node->getArc();arc!=cplx.arcs_end();arc=arc->getNext(node))
	      {
		node_it other = arc->getOtherNode(node);
		int otype = other->getType();
		if (otype>ntype)
		  getNodeExtendedManifold(other,out,up);
	      }
	  }
      }
    else
      {
	if (ntype>1)
	  {
	    for (arc_it arc=node->getArc();arc!=cplx.arcs_end();arc=arc->getNext(node))
	      {
		node_it other = arc->getOtherNode(node);
		int otype = other->getType();
		if (otype<ntype)
		  getNodeExtendedManifold(other,out,up);
	      }
	  }
      }
  } 

public:

  bool computePersistencePairs(bool computeCycles=false, bool force=false)
  {
    if (force) cplx.resetPPairsState();
    std::vector<int> cycles;
    cycles.push_back(0);
    if (computeCycles) cycles.push_back(1);
    bool res=cplx.computePersistencePairs(cycles);    
    return res;
  }

  bool removeBoundaries()
  {
    if (bType==BType_Torus) 
      {
	cplx.removeBoundaries(false);  
	return true;
      }
    return false;
  }

  bool enforceBoundaryConditions(BoundaryType type=BType_Default)
  {
    printf("Enforcing boundary conditions: ");fflush(0);
    if (!cplx.haveBoundaries())
      {
	printf("no boundaries.\n");
	manifoldWithBoundary=false;
	return false;
      }

    if (type==BType_Default) type=BType_Natural;       
    if (type==BType_Sphere)
      {
	fprintf(stderr,"\nWARNING: spherical boundary conditions not implemented yet.\n");
	fprintf(stderr,"         switching to 'natural' instead ...\n");
	type=BType_Natural;
      }

    if (bType==type)
      {
	char txt[255];
	if (type==BType_Natural) sprintf(txt,"natural");
	if (type==BType_Torus) sprintf(txt,"torus");
	if (type==BType_Sphere) sprintf(txt,"sphere");
	printf("%s -> %s. SKIPPING.\n",txt,txt);
	return false;
      }
    else
      {
	char txt[255];
	char txtb[255];
	if (type==BType_Natural) sprintf(txt,"natural");
	if (type==BType_Torus) sprintf(txt,"torus");
	if (type==BType_Sphere) sprintf(txt,"sphere");
	if (bType==BType_Natural) sprintf(txtb,"natural");
	if (bType==BType_NotSet) sprintf(txtb,"notset");
	if (bType==BType_Torus) sprintf(txtb,"torus");
	if (bType==BType_Sphere) sprintf(txtb,"sphere");
	printf("%s -> %s.\n",txtb,txt);
      }


    bType=type;
    
    if (bType==BType_Natural) manifoldWithBoundary=true;
    else manifoldWithBoundary=false;
    
    //manifoldWithBoundary=false;
    
    cplx.removeOutNodes();
    
    if (!manifoldWithBoundary) 
      {	
	if (bType==BType_Torus) cplx.buildReflectiveBoundary();
	cplx.setCancellationMode(NDcomplex::noBoundary);
      }
    else 
      {
	cplx.sanitizeBoundary(vertexAsMaxima, store_manifolds, store_arcsGeom);
	cplx.setCancellationMode(NDcomplex::preserveBoundary);
      }  

    cplx.resetPPairsState();
    cplx.resetCyclesState();
    return true;
  }

  void compute(bool gFilter = true, bool store_manifolds_p=true, int store_arcsGeom_p=((1<<0)|(1<<2)), bool vertexAsMinima=false, bool descendingFiltration=false)
  {
    //dbVec gradient_pairs;
    struct timeval wc_time1,wc_time2;
    struct tms tms_time1,tms_time2;
    long time1,time2;
    double t1,t2; 
    int ndims=network->getNDims();

   
    vertexAsMaxima=!vertexAsMinima;
    if (store_manifolds_p) store_manifolds=~0;
    else store_manifolds=0;
    store_arcsGeom = store_arcsGeom_p;

    printf ("Starting Morse-Smale complex computation.\n");
  
    
#ifdef USE_OPENMP
    t1 = omp_get_wtime( );
#else
    time1=times(&tms_time1);
#endif

    network->setGetValueByMax(!vertexAsMaxima);    
   
    computeDiscreteGradient(gradient_pairs,descendingFiltration);
    network->sendMessage(NetworkType::MSG_FACE_COFACE_FREE);
    setDiscreteGradientState(DG_available,true);
  
    if (debug_dump) {
      FILE *f=fopen("grad.dat","w");
      for (int i=0;i<gradient_pairs.size();i++)
	{
	  for (int j=0;j<gradient_pairs[i].size();j++)
	    {
	      if (gradient_pairs[i][j]==MSC_CRITICAL) continue;
	      if (gradient_pairs[i][j]==MSC_UNPAIRED) continue;
	      
	      if (gradient_pairs[i][j].type()>i)
		fprintf(f,"%d %d %ld\n",i,j,gradient_pairs[i][j].id());
	      else
		fprintf(f,"%d %d %ld\n",i,j,-gradient_pairs[i][j].id()-1);
	    }
	}
      fclose(f);
    }
    
    if (debug_dump) dumpDGradient("dg",true);
    //filterDGradient(gradient_pairs);
    computeMSComplex(gradient_pairs);  
    setDiscreteGradientState(DG_empty);
    gradient_pairs.clear(); 

    if (gFilter) {
      if (ndims>2)
	{
	  std::vector<double> level(ndims);
	  level.assign(ndims,0);
	  cplx.simplify(level);	  
	}
    }
    
#ifdef USE_OPENMP
    t2 = omp_get_wtime( );
#else
    time2=times(&tms_time2);
#endif
 
#ifdef USE_OPENMP
    printf ("Morse complex was computed in %.1f s.\n",t2-t1);
#else
    printf ("Morse complex was computed in %.1f s.\n",((double)(time2-time1)));
#endif
    
  }    

  void simplify(double NSig = 0.,double cut=0.,double cutR=0., bool simplifyGrps=true, bool doSanitize=true,bool skipLoops=false,bool dumpConflicts=false)
  {    
    struct timeval wc_time1,wc_time2;
    int ndims=network->getNDims();
    std::vector<double> level(ndims);
    printf ("Starting Morse-Smale complex simplification.\n");
    gettimeofday(&wc_time1, NULL);
    
    cplx.computePersistencePairs();

    if (NSig>0) 
      {
	double levelp[ndims];
	persistenceRatioFromSigma(NSig,ndims,levelp);
	level.assign(levelp,levelp+ndims);
	printf("Sampling noise level was set to %.1f-sigma.\n",NSig);
	simplifyFromPairs(level,true,skipLoops,dumpConflicts);
      }
    
    if (cut>0)
      {
	level.assign(ndims,cut);
	simplifyFromPairs(level,false,skipLoops,dumpConflicts);
      }

    if (cutR>0)
      {
	level.assign(ndims,cutR);
	simplifyFromPairs(level,true,skipLoops,dumpConflicts);
      }

    gettimeofday(&wc_time2, NULL);
    double dt  = (double)wc_time2.tv_sec + ((double)wc_time2.tv_usec)*1.E-6;
    dt -= (double)wc_time1.tv_sec + ((double)wc_time1.tv_usec)*1.E-6;
    printf ("Morse complex was simplified in %.1f s.\n",dt);
    
  }
  
  void printStats()
  {
    printf("%s",cplx.getInfo().c_str());fflush(0);
  }


};

const MSComplex::cellType MSComplex::MSC_CRITICAL = cellType(0,cellType::MAXIMUM);
const MSComplex::cellType MSComplex::MSC_UNPAIRED = cellType(1,cellType::MAXIMUM);
const MSComplex::cellType MSComplex::MSC_BOUNDARY_CRITICAL = cellType(2,cellType::MAXIMUM);



#endif
