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
#ifndef __ND_SUB_NETWORK_HXX__
#define __ND_SUB_NETWORK_HXX__

#include "network_interface.hxx"
#include "NDnetwork.h"
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <algorithm>
#include <set>
#include <map>
#include <vector>
#include <string>

template <class networkT>
class NDsubNetwork {
private:  

  NDnetwork *refNet;
  bool ownRefNet;
  const networkT *refNetI;
  long ndims;

  std::vector< std::set< NDNET_UINT > > subCellsUnique;
  std::vector< std::vector< NDNET_UINT > > subCellsDup;

  std::vector<std::string> complementaryDataName;
  std::vector<std::map<NDNET_UINT, long> > complementaryDataUnique_map;
  std::vector<std::vector<double> > complementaryDataUnique;
  std::vector<std::vector<double> > complementaryDataDup;
  
  //std::vector<uint> newVertexCount;
  NDNET_UINT newVertexCount;
  NDNET_UINT newVertexCSCount;
  NDNET_UINT newFakeVertexCount;
  NDNET_UINT newFakeVertexCSCount;
  
  typedef typename networkT::cellT cellT;

  typedef std::multimap<std::pair<char,NDNET_UINT>, NDNET_UINT> newVertexType;
  typedef newVertexType::iterator newVertexTypeIt;
  newVertexType newVertex; // uses center of mass
  newVertexType newVertexCS; // uses circumSphere

  typedef std::vector< std::vector< newVertexTypeIt > > newCellType;
  typedef newCellType::iterator newCellTypeIt;

  std::vector< std::set<NDNET_UINT> > insertedNewCells;

  std::vector< std::vector< newVertexTypeIt > > newCell;
  std::vector< std::vector< newVertexTypeIt > > newCellCS;

  template<class TypeT, class IdT>
  std::pair<newVertexTypeIt,bool> insertAsVertexPrivate(const std::pair<TypeT,IdT> &p, bool useCOM=true, bool unique=true)
  {  
    std::pair<newVertexTypeIt,bool> res;

    if (unique)
      {
	if (useCOM) res.second=((res.first=newVertex.find(p))==newVertex.end());
	else res.second=((res.first=newVertexCS.find(p))==newVertexCS.end());
	if (!res.second) return res;
      }

    res.second=true;
    if (p.first==0) 
      {
	insert(p,useCOM,unique);
	if (useCOM) res.first = newVertex.insert(make_pair(p,newFakeVertexCount++));
	else res.first = newVertexCS.insert(make_pair(p,newFakeVertexCSCount++));
      }
    else
      { 
	if (useCOM) res.first  = newVertex.insert(make_pair(p,newVertexCount++));
	else res.first  = newVertexCS.insert(make_pair(p,newVertexCSCount++));
      }
  
    return res;

  }

  template<class TypeT, class IdT>
  bool insertSubCell(const std::pair<TypeT,IdT> &p, bool unique=true)
  {
    bool res=true;
    if ((p.first==0)||(unique)) 
      res=subCellsUnique[p.first].insert((NDNET_UINT)p.second).second;  
    else
      subCellsDup[p.first].push_back((NDNET_UINT)p.second); 
          
    return res;
  }

private:

  template<class TypeT, class IdT>
  bool insertSubCell(const std::pair<TypeT,IdT> &p, std::vector<double>::iterator compB, std::vector<double>::iterator compE,bool unique=true)
  {
    bool res=true;
    
    assert(static_cast<long>(compE-compB) == complementaryDataName.size());
   
    if ((p.first==0)||(unique)) 
      {
	res=subCellsUnique[p.first].insert((NDNET_UINT)p.second).second;  
	if (res)
	  {
	    if (compB!=compE)
	      {
		complementaryDataUnique_map[p.first][(NDNET_UINT)p.second]=complementaryDataUnique[p.first].size();
		for (std::vector<double>::iterator it=compB; it!=compE;it++)
		  complementaryDataUnique[p.first].push_back(*it);
	      }
	  }
	else
	  {
	    
	  }
      }
    else
      {
	subCellsDup[p.first].push_back((NDNET_UINT)p.second); 
	if (compB!=compE)
	  {	    
	    for (std::vector<double>::iterator it=compB; it!=compE;it++)
	      complementaryDataDup[p.first].push_back(*it);
	  }
      }
    
    return res;
  }

public:
  
  void reset()
  {
    subCellsUnique.clear();
    subCellsDup.clear();
    newVertex.clear();
    newVertexCS.clear();
    newCell.clear();
    newCellCS.clear();
    insertedNewCells.clear();
    complementaryDataName.clear();
    complementaryDataUnique_map.clear();
    complementaryDataUnique.clear();
    complementaryDataDup.clear();
    newVertexCount=0;
    newVertexCSCount=0;
    newFakeVertexCount=0;
    newFakeVertexCSCount=0;

    setRefNetwork(NULL);
  }

  long setComplementaryDataName(std::vector<std::string> &str)
  {
    complementaryDataName.assign(str.begin(),str.end());
    return str.size();
  }

  bool setRefNetwork(NDnetwork *net)
  {   
    if (ownRefNet) {
      FreeNDnetwork(&refNet);
      ownRefNet=false;
    }

    refNet=net;
    ownRefNet=false;

    if (refNet==NULL)
      {
	refNet = CreateNetwork(ndims,0,0);
	std::vector<double> x0;
	std::vector<double> delta;
	int i;
	refNetI->getBoundingBox(x0,delta);
	for (i=0;i<ndims;i++)
	  {
	    refNet->x0[i]=x0[i];
	    refNet->delta[i]=delta[i];
	  }
	refNet->periodicity=refNetI->getPeriodicity();
	strcpy(refNet->comment,"dummy");
	refNet->isSimpComplex=1;
	ownRefNet=true;
      }
    
    subCellsUnique.resize(ndims+1);
    subCellsDup.resize(ndims+1);
    insertedNewCells.resize(ndims+1);
    
    complementaryDataDup.resize(ndims+1);
    complementaryDataUnique.resize(ndims+1);
    complementaryDataUnique_map.resize(ndims+1);
  }
  
  NDsubNetwork(const networkT *netI, NDnetwork *ndnet=NULL)
  {  
    refNetI=netI;
    ndims = refNetI->getNDims();
    ownRefNet=false;
    setRefNetwork(ndnet);   
    newVertexCount=0;
    newVertexCSCount=0;
    newFakeVertexCount=0;
    newFakeVertexCSCount=0;
  }

  ~NDsubNetwork()
  {
    if (ownRefNet)
      {
	FreeNDnetwork(&refNet);
      }
  }  

  bool insert(typename networkT::subCellsType &subList, typename networkT::cellsInfoType &subListInfo, bool useCOM=true, bool unique=true)
  {
    typename networkT::subCellsType::iterator s_it;
    typename networkT::cellsInfoType::second_type::iterator si_it;
    
    long nsup=setComplementaryDataName(subListInfo.first);

    assert(subList.size()*nsup == subListInfo.second.size());

    si_it=subListInfo.second.begin();
    for (s_it=subList.begin();s_it!=subList.end();s_it++,si_it+=nsup)
      insert(s_it->getAsPair(),si_it,si_it+nsup,useCOM,unique);
  }

  template<class TypeT, class IdT>
  bool insert(const std::pair<TypeT,IdT> &p, bool useCOM=true, bool unique=true)
  {
    std::vector<double> tmp(complementaryDataName.size(),0);   
    insert(p,tmp.begin(),tmp.end(),useCOM,unique);
  }

private:

  template<class TypeT, class IdT>
  bool insert(const std::pair<TypeT,IdT> &p, std::vector<double>::iterator compB, std::vector<double>::iterator compE, bool useCOM=true, bool unique=true)
  {
    if (p.first<=ndims) return insertSubCell(p,compB,compE,unique);    
    std::vector<cellT> v;
    std::vector<cellT> v2;

    char type = p.first-ndims-1;
    
    if (unique)
      {
	if (!insertedNewCells[type].insert(p.second).second) 
	  return false;
      }
    
    if (ndims == type) return insertAsVertex(std::make_pair(type,(NDNET_UINT)p.second),useCOM,unique);

    std::vector< std::pair<char,NDNET_UINT> > pts(1+ndims - type);

    if (ndims-1 == type) {
      refNetI->getCofaces(cellT(type,p.second),v);
      long n=v.size();
           
      if (n==2) {
	pts[0]=std::make_pair(type+1,v[0].id());
	pts[1]=std::make_pair(type+1,v[1].id());
	addNewCell(pts.begin(),pts.end(),useCOM,unique);
      }
      else return false;
      
    }
    else if (ndims-2 == type) {
      int i;
      bool haveRef=false;
      refNetI->getCofaces(cellT(type,p.second),v);
      long n=v.size();
      
      
      for (i=0;i<n;i++) {
	refNetI->getCofaces(v[i],v2);
	long n2=v2.size();
		
	if (n2!=2) continue;
	if (!haveRef) {
	  haveRef=true;
	  pts[0] = std::make_pair(ndims,v2[0].id());
	}
	else if ((v2[0].id()!=pts[0].second)&&(v2[1].id()!=pts[0].second))
	  {
	    pts[1] = std::make_pair(ndims,v2[0].id());
	    pts[2] = std::make_pair(ndims,v2[1].id());
	    addNewCell(pts.begin(),pts.end(),useCOM,unique);
	  }
      }
     
      return true;
    }
    else {
      //  TO BE IMPLEMENTED !!!!!!!
      /*
      // just a hack, should actually compute voronoi cells here
      int n = FacesInFace(refNet, p.second,type,ndims,&v);
      pts.resize(2);
      pts[0] = std::make_pair(type,p.second);
      int i,j;
      for (i=0;i<n;i++) {
	pts[1]=std::make_pair(ndims,v[i]);
	addNewCell(pts.begin(),pts.end(),false);
      }
      */

    }
    
    return true;
  }
  
public:

  template<class TypeT, class IdT>
  bool insertAsVertex(const std::pair<TypeT,IdT> &p, bool useCOM=true, bool unique=true)
  {    
    insertAsVertexPrivate(p,useCOM,unique);  
  }

  template<class InputIterator>  
  int insertAsVertex(InputIterator begin, InputIterator end, bool useCOM=true, bool unique=true)
    {
    InputIterator id;
    long n=0;

    for (id=begin;id!=end;id++)
      {
	if (insertAsVertex(*id,useCOM,unique).second) n++;
      }

    return n;
  }  

  template<class InputIterator1,class InputIterator2>
  int insertAsVertex(InputIterator1 type_begin, InputIterator1 type_end,InputIterator2 id_begin, InputIterator2 id_end, bool useCOM=true, bool unique=true)
    {
    InputIterator1 type;
    InputIterator2 id;
    long n=0;

    for (type=type_begin,id=id_begin;type!=type_end;type++)
      {
	if (insertAsVertex(std::make_pair(*type,*id),useCOM,unique).second) n++;
	if (id==id_end()) {
	  fprintf(stderr,"Error in NDsubNetwork::insert, more Ids than Types\n");
	  exit(0);
	}
      }
    if (id!=id_end()) {
      fprintf(stderr,"Error in NDsubNetwork::insert, more types than Ids\n");
      exit(0);
    }

    return n;
  }



  template<class InputIterator,class T>
  int insertAsVertex(T type,InputIterator id_begin, InputIterator id_end, bool useCOM=true, bool unique=true)
    {
    
    InputIterator id;
    long n=0;

    for (id=id_begin;id!=id_end;id++)
      {
	if (insertAsVertex(std::make_pair(type,*id),useCOM,unique).second) n++;
      }

    return n;
  }  
 
  template<class InputIterator>
  bool addNewCell(InputIterator begin, InputIterator end, bool useCOM=true, bool unique=true)
  {
    std::vector< newVertexTypeIt > c;
    
    for (InputIterator it=begin;it!=end;it++)
      {	
	std::pair<newVertexTypeIt,bool> res=insertAsVertexPrivate(*it,useCOM,unique);
	c.push_back(res.first);
      }
    
    if (useCOM) newCell.push_back(c);
    else newCellCS.push_back(c);
    return true;
  }

  NDnetwork *create()
  {
 
    long i,j,k;
    std::vector<double *> cellData(ndims+1,(double*)NULL);
    std::vector<long> nSubCells(ndims+1);
    std::vector< std::vector<double*> > compData(ndims+1);
    long nCompData=complementaryDataName.size();
    
    for (i=0;i<ndims+1;i++) 
      {       
	nSubCells[i]=subCellsUnique[i].size()+subCellsDup[i].size();
	if (nCompData) compData[i].assign(nCompData,(double*)NULL);
      }
    
    NDnetwork *net = CreateNetwork(ndims,nSubCells[0],0);
    std::set< NDNET_UINT >::iterator uit;
    std::vector< NDNET_UINT >::iterator dit;
    std::vector<int> newId;

    sprintf(net->comment,"(subNet) %s",refNet->comment);
    memcpy(net->x0,refNet->x0,sizeof(double)*net->ndims);
    memcpy(net->delta,refNet->delta,sizeof(double)*net->ndims);
    net->periodicity = refNet->periodicity;   
    
    for (i=0,j=0;i<refNet->ndata;i++)
      if (refNet->data[i].type==0) j++;
    net->ndata = j;
    if (net->ndata) net->data=(NDnetwork_Data *)calloc(net->ndata,sizeof(NDnetwork_Data));

    for (i=0,j=0;i<refNet->ndata;i++)
      if (refNet->data[i].type==0) 
	{
	  strcpy(net->data[j].name,refNet->data[i].name);
	  net->data[j].type=0;
	  if (nSubCells[0])
	    net->data[j].data=(double*)malloc(sizeof(double)*nSubCells[0]);
	  else net->data[j].data=NULL;
	  j++;
	}
    
    for (i=0;i<nCompData;i++)
      if (nSubCells[0])
	compData[0][i]=(double*)malloc(sizeof(double)*nSubCells[0]);	  
    
   if (nSubCells[0]) cellData[0] = (double*)malloc(sizeof(double)*nSubCells[0]);

      //compData[0].resize(nSubCells[0]*nCompData);
    
    i=0;
    newId.assign(refNetI->getNFaces(0),-1);
    for (uit=subCellsUnique[0].begin();uit!=subCellsUnique[0].end();uit++,i++)
      {
	refNetI->getPosition(cellT(0,*uit),&net->v_coord[i*net->ndims]); // (0,i) from (0,*uit)
	for (j=0,k=0;j<refNet->ndata;j++)
	  if (refNet->data[j].type==0)
	    net->data[k++].data[i]=refNet->data[j].data[(*uit)];
			
	newId[(long) *uit] = i;

	cellData[0][i]=cellT(0,*uit).getAsDouble();
	if (nCompData) 
	  {
	    long ref=complementaryDataUnique_map[0].find(*uit)->second;
	    for (j=0;j<nCompData;j++)
	      compData[0][j][i]=complementaryDataUnique[0][ref+j];
	  }
      }
    long i0=i;
    for (dit=subCellsDup[0].begin();dit!=subCellsDup[0].end();dit++,i++)
      {
	refNetI->getPosition(cellT(0,*dit),&net->v_coord[i*net->ndims]);  // (0,i) from (0,*dit)
	for (j=0,k=0;j<refNet->ndata;j++)
	  if (refNet->data[j].type==0)
	    net->data[k++].data[i]=refNet->data[j].data[(*dit)];

	newId[(long) *dit] = i;

	cellData[0][i]=cellT(0,*dit).getAsDouble();	
	if (nCompData) 
	  {
	    for (j=0;j<nCompData;j++)
	      compData[0][j][i]=complementaryDataDup[0][(i-i0)*nCompData+j];
	  }
      }
    
    std::vector<int> addedVertex;
    
    for (i=1;i<=ndims;i++)
      {
	net->nfaces[i]=0;
	if (nSubCells[i])
	  {
	    //compData[i].resize(nSubCells[i]*nCompData);
	    for (j=0;j<nCompData;j++)
	      compData[i][j]=(double*)realloc(compData[i][j],sizeof(double)*nSubCells[i]);
	    cellData[i]=(double*)realloc(cellData[i],sizeof(double)*nSubCells[i]);
	    net->haveVertexFromFace[i]=1;
	    net->f_vertexIndex[i] = (NDNET_UINT *) malloc((i+1)*sizeof(NDNET_UINT)*nSubCells[i]);
	    net->nfaces[i] = nSubCells[i];
	  }
	
	NDNET_UINT *Id=net->f_vertexIndex[i];
	//NDNET_UINT *refId=refNet->f_vertexIndex[i];
	int tmpi;

	std::vector<cellT> v;		
	for (uit=subCellsUnique[i].begin(),j=0;uit!=subCellsUnique[i].end();uit++,j++)
	  {
	    refNetI->getVertice(cellT(i,*uit),v);
	    for (k=0;k<i+1;k++)
	      {
		NDNET_UINT curid = v[k].id();
		//NDNET_UINT curid = refId[(*uit)*(i+1)+k];
		
		//printf ("%d_",*uit);fflush(0);
		tmpi=newId[curid];
		if (tmpi == -1)
		  {
		    newId[curid] = net->nvertex+addedVertex.size();
		    addedVertex.push_back(curid);
		    tmpi=newId[curid];
		  }
		Id[j*(i+1)+k] = (NDNET_UINT)tmpi; //(i,j) from (i,*uit)		
	      }
	    cellData[i][j]=cellT(i,*uit).getAsDouble();
	    if (nCompData) 
	      {
		long ref=complementaryDataUnique_map[i].find(*uit)->second;
		for (k=0;k<nCompData;k++)
		  compData[i][k][j]=complementaryDataUnique[i][ref+k];
	      }
	  }
	long j0=j;
	for (dit=subCellsDup[i].begin();dit!=subCellsDup[i].end();dit++,j++)
	  {
	    refNetI->getVertice(cellT(i,*dit),v);
	    for (k=0;k<i+1;k++)
	      {
		NDNET_UINT curid = v[k].id();
		//NDNET_UINT curid = refId[(*dit)*(i+1)+k];
		
		//printf ("%d_",*dit);fflush(0);
		tmpi=newId[curid];
		if (tmpi == -1)
		  {
		    newId[curid] = net->nvertex+addedVertex.size();
		    addedVertex.push_back(curid);
		    tmpi=newId[curid];
		  }
		Id[j*(i+1)+k] = (NDNET_UINT)tmpi; //(i,j) from (i,*dit)
	      }
	    cellData[i][j]=cellT(i,*dit).getAsDouble();
	    if (nCompData) 
	      {
		for (k=0;k<nCompData;k++)
		  compData[i][k][j]=complementaryDataDup[i][(j-j0)*nCompData];
	      }
	  }
      }

    if (addedVertex.size())
      {
	net->v_coord = (float*) realloc(net->v_coord,sizeof(float)*net->ndims*(net->nvertex+addedVertex.size()));
	for (i=0;i<net->ndata;i++) 
	  net->data[i].data= (double*) realloc(net->data[i].data,sizeof(double)*(net->nvertex+addedVertex.size())); 
	for (i=0;i<nCompData;i++)
	  compData[0][i]=(double*)realloc(compData[0][i],sizeof(double)*(net->nvertex+addedVertex.size()));
	cellData[0]=(double*)realloc(cellData[0],sizeof(double)*(net->nvertex+addedVertex.size()));
	for (i=0;i<addedVertex.size();i++,net->nvertex++)
	  {
	    refNetI->getPosition(cellT(0,addedVertex[i]),&net->v_coord[net->nvertex*net->ndims]); //(0,net->nvertex) from (0,addedVertex[i])
	    for (j=0;j<nCompData;j++) compData[0][j][net->nvertex]=0; // replace DBL_MIN ?
	    cellData[0][net->nvertex]=cellT(0,addedVertex[i]).getAsDouble();
	    for (j=0,k=0;j<refNet->ndata;j++)
	      if (refNet->data[j].type==0)
		net->data[k++].data[net->nvertex]=refNet->data[j].data[addedVertex[i]];	    
	  }
	  
      }
    addedVertex.clear();
    
    
    // center of mass new cells
    NDNET_UINT nvert=net->nvertex;
    
    net->nvertex+=newVertexCount;
    net->v_coord = (float*) realloc(net->v_coord,sizeof(float)*net->ndims*net->nvertex);
    
    for (i=0;i<net->ndata;i++)       
      net->data[i].data= (double*) realloc(net->data[i].data,sizeof(double)*net->nvertex);    
    for (i=0;i<nCompData;i++)
      compData[0][i]=(double*)realloc(compData[0][i],sizeof(double)*(net->nvertex));
    cellData[0]=(double*)realloc(cellData[0],sizeof(double)*(net->nvertex));
    for (newVertexTypeIt vit=newVertex.begin();vit!=newVertex.end();vit++)
      {
	if ((vit->first).first != 0)
	  {	    
	    std::vector<cellT> vl;		
	    refNetI->getVertice(cellT(vit->first.first,vit->first.second),vl);
	    
	    //computes COM
	    refNetI->getPosition(cellT(vit->first.first,vit->first.second),
				 &net->v_coord[net->ndims*(nvert+vit->second)]); //(0,(nvert+vit->second)) from (vit->first.first,vit->first.second)
	 
	    for (i=0,j=0;i<refNet->ndata;i++)
	      if (refNet->data[i].type==0)
		net->data[j++].data[nvert+vit->second]=0;

	    cellData[0][nvert+vit->second]=cellT(vit->first.first,vit->first.second).getAsDouble();
	    for (k=0;k<nCompData;k++) compData[0][k][nvert+vit->second]=0;

	    for (k=0;k<vit->first.first+1;k++)
	      {
		for (i=0,j=0;i<refNet->ndata;i++)
		  {   
		    if (refNet->data[i].type==0)
		      net->data[j++].data[nvert+vit->second]+=
			refNet->data[i].data[vl[k].id()];
		  }
	      }
	    for (i=0,j=0;i<refNet->ndata;i++)
	      if (refNet->data[i].type==0)
		net->data[j++].data[nvert+vit->second]/=(vit->first.first+1);
	  }
      }    
    
    std::vector<NDNET_UINT> newCellsCount(net->ndims+1);
    for (newCellTypeIt cit=newCell.begin();cit!=newCell.end();cit++)
      {
	newCellsCount[cit->size()-1] ++;
      }
    
    for (i=1;i<net->ndims+1;i++){
      if (newCellsCount[i]) {
	net->haveVertexFromFace[i]=1;	
	net->f_vertexIndex[i] = (NDNET_UINT*)realloc(net->f_vertexIndex[i],(net->nfaces[i]+newCellsCount[i])*(i+1)*sizeof(NDNET_UINT));
	cellData[i]=(double*)realloc(cellData[i],sizeof(double)*(net->nfaces[i]+newCellsCount[i]));
	for (j=0;j<nCompData;j++)
	  compData[i][j]=(double*)realloc(compData[i][j],sizeof(double)*(net->nfaces[i]+newCellsCount[i]));
      }
    }

    for (newCellTypeIt cit=newCell.begin();cit!=newCell.end();cit++)
      {
	NDNET_UINT type=cit->size()-1;
	NDNET_UINT curnf=net->nfaces[type]; //(type,curnf) from (cit->first.first,cit->first.second)->retrouver l ID, ca c'est la liste des Vertex.
	
	for (i=0;i<type+1;i++)
	  {
	    newVertexTypeIt &vert = (*cit)[i];
	    if (vert->first.first == 0) 
	      net->f_vertexIndex[type][curnf*(type+1)+i] = newId[vert->first.second]; 
	    else
	      net->f_vertexIndex[type][curnf*(type+1)+i] = vert->second+nvert; 
	  }
	cellData[type][curnf]=0;
	for (j=0;j<nCompData;j++)
	  compData[type][j][curnf]=0;
	net->nfaces[type]++;
      }


    
    // CircumCircle new cells
    nvert=net->nvertex;
    net->nvertex+=newVertexCSCount;
    //printf("+%dCS\n",newVertexCSCount);
    net->v_coord = (float*) realloc(net->v_coord,sizeof(float)*net->ndims*net->nvertex);
    for (i=0;i<net->ndata;i++) 
      net->data[i].data= (double*) realloc(net->data[i].data,sizeof(double)*net->nvertex);
    for (i=0;i<nCompData;i++)
      compData[0][i]=(double*)realloc(compData[0][i],sizeof(double)*(net->nvertex));
    cellData[0]=(double*)realloc(cellData[0],sizeof(double)*(net->nvertex));
    for (newVertexTypeIt vit=newVertexCS.begin();vit!=newVertexCS.end();vit++)
      {
	if ((vit->first).first != 0)
	  {
	    std::vector<cellT> vl;		
	    refNetI->getVertice(cellT(vit->first.first,vit->first.second),vl);

	    double v[ndims+1][ndims];
	    double *pv[ndims+1];
	    for (i=0;i<vl.size();i++)
	      {
		float pos[ndims];
		refNetI->getPosition(vl[i],pos);
		for (j=0;j<ndims;j++) v[i][j]=pos[j];
		pv[i]=&v[i][0];
	      }

	    Simplex_CS_fc(refNet,vit->first.first,pv,&net->v_coord[net->ndims*(nvert+vit->second)]); //(0,(nvert+vit->second)) from (vit->first.first,vit->first.second)
	  
	    for (i=0,j=0;i<refNet->ndata;i++)
	      if (refNet->data[i].type==0)
		net->data[j++].data[nvert+vit->second]=0;

	    cellData[0][nvert+vit->second]=cellT(vit->first.first,vit->first.second).getAsDouble();
	    for (k=0;k<nCompData;k++) compData[0][k][nvert+vit->second]=0;


	    for (k=0;k<vit->first.first+1;k++)
	      {
		for (i=0,j=0;i<refNet->ndata;i++)
		  {   
		    if (refNet->data[i].type==0)
		      net->data[j++].data[nvert+vit->second]+=
			refNet->data[i].data[vl[k].id()];
		  }
	      }
	    for (i=0,j=0;i<refNet->ndata;i++)
	      if (refNet->data[i].type==0)
		net->data[j++].data[nvert+vit->second]/=(vit->first.first+1);
	  }
      }    

    std::vector<NDNET_UINT> newCellsCSCount(net->ndims+1);
    for (newCellTypeIt cit=newCellCS.begin();cit!=newCellCS.end();cit++)
      {
	newCellsCSCount[cit->size()-1] ++;
      }
    
    for (i=1;i<net->ndims+1;i++){
      if (newCellsCSCount[i]) {
	net->haveVertexFromFace[i]=1;
	net->f_vertexIndex[i] = (NDNET_UINT*)realloc(net->f_vertexIndex[i],(net->nfaces[i]+newCellsCSCount[i])*(i+1)*sizeof(NDNET_UINT));
	cellData[i]=(double*)realloc(cellData[i],sizeof(double)*(net->nfaces[i]+newCellsCSCount[i]));
	for (j=0;j<nCompData;j++)
	  compData[i][j]=(double*)realloc(compData[i][j],sizeof(double)*(net->nfaces[i]+newCellsCSCount[i]));
      }
    }

    for (newCellTypeIt cit=newCellCS.begin();cit!=newCellCS.end();cit++)
      {
	NDNET_UINT type=cit->size()-1;
	NDNET_UINT curnf=net->nfaces[type];  //(type,curnf) from (cit->first.first,cit->first.second)->retrouver l ID, ca c'est la liste des Vertex.

	for (i=0;i<type+1;i++)
	  {	    
	    newVertexTypeIt &vert = (*cit)[i];
	    if (vert->first.first == 0) 
	      net->f_vertexIndex[type][curnf*(type+1)+i] = newId[vert->first.second]; 
	    else
	      net->f_vertexIndex[type][curnf*(type+1)+i] = vert->second+nvert;
	  }
	cellData[type][curnf]=0;
	for (j=0;j<nCompData;j++)
	  compData[type][j][curnf]=0;
	//printf("t%d : %d\n",type,net->nfaces[type]);
	net->nfaces[type]++;
      }
    
    for (i=0;i<nCompData;i++)
      {
	for (j=0;j<ndims+1;j++)
	  {
	    if (compData[j][i]==NULL) continue;
	    addNDDataArr(net,j,complementaryDataName[i].c_str(),&(compData[j][i]));	  
	  }
      }
    for (j=0;j<ndims+1;j++)
      if (cellData[j]!=NULL)
	{
	  addNDDataArr(net,j,CELL_TAG,&cellData[j]);
	}

    printNDnetStat(net,4);
    return net;
  }


};


#endif
