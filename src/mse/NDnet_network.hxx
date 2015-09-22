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
#ifndef __NDNET_MSC_HEADER__
#define __NDNET_MSC_HEADER__

#include "network_interface.hxx"

#include "MSComplex.hxx"

#include "NDnetwork.h"
#include "mytypes.h"

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <vector>
#include <algorithm>
#include "global.h"
#include "NDsubNetwork.hxx"

#include "NDnet_interface.hxx"

//template <typename TypeT = char, typename IdT = uint>
template <class cellType = cellIdentity<> >
class NDnet_network : public NetworkDataInterface<cellType>  {

  typedef class NetworkDataInterface<cellType> ParentClass;  
  typedef typename ParentClass::messageT messageT;
  typedef typename cellType::typeT typeT;
  typedef typename cellType::idT idT;
  
  //typedef typename ParentClass::cellID cellID;
  typedef typename ParentClass::subCellsElementType subCellsElementType;
  typedef typename ParentClass::subCellsType subCellsType;
  typedef typename ParentClass::cellsInfoType cellsInfoType;
  typedef typename ParentClass::newCellsElementType newCellsElementType;
  typedef typename ParentClass::newCellsType newCellsType;
  

private:
  NDnetwork *net;
  //typeT type;
  //idT id;
  
private:
  char valueTag[255];
  double *vertexVal;

  char groupTag[255];
  double *vertexGroup;
  int periodicity;


  std::pair<long,long> setBoundaries()
  {
    cellType curCell;
    char type=getNDims(true)-1;
    long nf=getNFaces(type);
    long i;
    long nb=0;
    long nnm=0;
    unsigned char *ff[getNDims(true)+1];
    
    NDNetFlags_enable(net,0);
    ff[0]=net->v_flag;
    for (i=1;i<=getNDims(true);i++) 
      {
	NDNetFlags_enable(net,i);
	ff[i]=net->f_flag[i];
      }	  
    
    for (i=0;i<net->nvertex;i++) 
      ff[0][i] &= (~NDNETFLAG_BOUNDARY);
    
    for (type=1;type<=getNDims(true);type++)
      for (i=0;i<net->nfaces[type];i++) 
	ff[type][i]&=(~NDNETFLAG_BOUNDARY);
       
    type=getNDims(true)-1;
    nf=getNFaces(type);
    for (i=0;i<nf;i++)
      {
	std::vector<cellType> coface;	
	curCell.set(type,i);
	 
	if (isOut(curCell)) continue;	
	getCofaces(curCell,coface);

	bool isb=(coface.size()<2);
	bool isnm=(coface.size()>2);

	if (!isb)
	  {	    
	    if (isOut(coface[0])||isOut(coface[1])) isb=true;
	  }	
	else
	  {
	    if (coface.size()==0) isb=false;
	  }
	
	if (isb) 
	  {
	    std::vector<cellType> face;
	    
	    ff[type][i] |= NDNETFLAG_BOUNDARY;
	    //mask_boundary[type][i]=true;
	    nb++;
	    
	    getFaces(curCell,face);
	    for (typename std::vector<cellType>::iterator curf = face.begin(); curf!=face.end(); curf++)
	      {
		//mask_boundary[type-1][curf->id()]=true;
		ff[type-1][curf->id()] |= NDNETFLAG_BOUNDARY;
	      }
	  }
	
	if (isnm) 
	  {
	    std::vector<cellType> face;
	    ff[type][i] |= NDNETFLAG_NON_MANIFOLD;
	    nnm++;
	    getFaces(curCell,face);
	    for (typename std::vector<cellType>::iterator curf = face.begin(); curf!=face.end(); curf++)
	      {
		//mask_boundary[type-1][curf->id()]=true;
		ff[type-1][curf->id()] |= NDNETFLAG_NON_MANIFOLD;
	      }
	  }
      }
    
    for (type=this->getNDims(true)-2;type>0;type--)
      {
	nf=this->getNFaces(type);
	for (i=0;i<nf;i++)
	  {
	    curCell.set(type,i);
	    if (isBoundary(curCell))
	      {
		static std::vector<cellType> face;
		
		getFaces(curCell,face);
		for (typename std::vector<cellType>::iterator curf = face.begin(); curf!=face.end(); curf++)
		  {
		    ff[type-1][curf->id()] |= NDNETFLAG_BOUNDARY;
		    //mask_boundary[type-1][curf->id()]=true;
		  }
	      }
	    if (ff[type][i]&NDNETFLAG_NON_MANIFOLD)
	      {
		static std::vector<cellType> face;
		
		getFaces(curCell,face);
		for (typename std::vector<cellType>::iterator curf = face.begin(); curf!=face.end(); curf++)
		  {
		    ff[type-1][curf->id()] |= NDNETFLAG_NON_MANIFOLD;
		    //mask_boundary[type-1][curf->id()]=true;
		  }
	      }
	  }
      }

    return std::make_pair(nb,nnm);
  }




public:
  void freeData() {FreeNDnetwork(&net);}

  int getPeriodicity() const {return periodicity;} 

  int setPeriodicity(int p)
  {
    periodicity=p;
    return p;
  }

  void setMask(std::vector<char> &mask, bool nullIsMasked)
  {
    long i,j;
    unsigned char *ff[this->getNDims()+1];
    long nf;
    bool noOut=false;
    cellType curCell;

    if (mask.size())
      {
	if (mask.size() != net->nvertex)
	  {
	    fprintf(stderr,"ERROR in setMask: mask and data do not have the same size.\n");
	    exit(-1);
	  }
      }

    NDNetFlags_enable(net,0);
    ff[0]=net->v_flag;
    for (i=1;i<=this->getNDims();i++) 
      {
	NDNetFlags_enable(net,i);
	ff[i]=net->f_flag[i];
      }	  

    if (mask.size())
      {
	if (nullIsMasked)
	  {
	    for (i=0;i<net->nvertex;i++)
	      if (!mask[i]) 
		net->v_flag[i]|=NDNETFLAG_OUT; 
	  }
	else
	  {
	    for (i=0;i<net->nvertex;i++)
	      if (mask[i]) 
		net->v_flag[i]|=NDNETFLAG_OUT; 
	  }
      }
    else
      {	
	j=-1;
	do {
	  nf = this->getNFaces(++j);
	  for (i=0;i<nf;i++)
	    if (ff[j][i]&NDNETFLAG_OUT) break;
	    
	} while ((i==nf)&&(j<this->getNDims()));
	if ((j==this->getNDims())&&(i==nf)) noOut=true;
      }
    
    printf ("Building mask ... ");fflush(0);
    
    if (!noOut)
      for (char type=1;type<=this->getNDims(true);type++) {
	long nf = this->getNFaces(type);
	for (i=0;i<nf;i++)
	  {
	    static std::vector<cellType> vert;
	    if (ff[type][i]&NDNETFLAG_OUT) continue;
	    getVertice(cellType(type,i),vert);
	    for (j=0;j<vert.size();j++) if (isOut(vert[j])) break;
	    if (j!=vert.size()) ff[type][i]|=NDNETFLAG_OUT; 
	  }
      }
        
    std::pair<long,long> fail=setBoundaries();
    long nb=fail.first;
    long nnm=fail.second;
    printf ("done. (%d %d-faces on boundary)\n",(int)nb,net->ndims-1);

    if (nnm)
      {
	printf ("WARNING: this is not a manifold! (%ld non-manifold faces found).\n",nnm);
	printf ("I will try to fix that for you ... ");fflush(0);
	printf("done.\n");
	printf ("\n well actually I can't and I am probably going to crash soon ;)\n");fflush(0);

      }

  }

  bool isOut(cellType cell) const 
  {
    typeT type=cell.type();
    idT id=cell.id();

    if (type==0) 
      return (bool)(net->v_flag[id]&NDNETFLAG_OUT); 
    else 
      return (bool)(net->f_flag[type][id]&NDNETFLAG_OUT);
  }
  
  bool isBoundary(cellType cell) const 
  {
    typeT type=cell.type();
    idT id=cell.id();
	  
    if (type>=getNDims(true)) return false;

    if (type==0) 
      return (bool)(net->v_flag[id]&NDNETFLAG_BOUNDARY); 
    else 
      return (bool)(net->f_flag[type][id]&NDNETFLAG_BOUNDARY);
  }

  double getValueSum( cellType cell, double &sum, int getByMax)  const
  {
    typeT type=cell.type();
    idT id=cell.id();
	  
    if (type==0) {sum=vertexVal[id];return vertexVal[id];}
    if (getByMax==-1) getByMax=ParentClass::getByMax();
	  
    NDNET_UINT *vid = VERTEX_IN_FACE(net,type,id);	  
    NDNET_UINT max=vid[0];
    sum=vertexVal[max];  

    if (getByMax)
      for (NDNET_UINT i=1;i<NUM_VERTEX_IN_FACE(net,type,id);i++)
	{
	  sum+=vertexVal[vid[i]];
	  if (vertexVal[vid[i]]>vertexVal[max]) max = vid[i];
	}
    else
      for (NDNET_UINT i=1;i<NUM_VERTEX_IN_FACE(net,type,id);i++)
	{
	  sum+=vertexVal[vid[i]];
	  if (vertexVal[vid[i]]<vertexVal[max]) max = vid[i];
	}
	  
    return vertexVal[max];
  }

  double getValueAll( cellType cell, std::vector<double> &values, int getByMax)  const
  {
    typeT type=cell.type();
    idT id=cell.id();

    if (type==0) {values.resize(1);values[0]=vertexVal[id];return vertexVal[id];}
    if (getByMax==-1) getByMax=ParentClass::getByMax();

    NDNET_UINT  *vid = VERTEX_IN_FACE(net,type,id);
    NDNET_UINT max=0;
    values.resize(NUM_VERTEX_IN_FACE(net,type,id));
    values[0]=vertexVal[vid[0]];
	    
    //sum=vertexVal[max];
	    
    if (getByMax)
      for (NDNET_UINT i=1;i<NUM_VERTEX_IN_FACE(net,type,id);i++)
	{
	  //sum+=vertexVal[vid[i]];
	  values[i]=vertexVal[vid[i]];
	  if (vertexVal[vid[i]]>vertexVal[vid[max]]) max = i;
	}
    else
      for (NDNET_UINT i=1;i<NUM_VERTEX_IN_FACE(net,type,id);i++)
	{
	  //sum+=vertexVal[vid[i]];
	  values[i]=vertexVal[vid[i]];
	  if (vertexVal[vid[i]]<vertexVal[vid[max]]) max = i;
	}

    if (max!=0) std::swap(values[0],values[max]);
	    
    return values[0];
  }

  double getValue(cellType cell, int getByMax) const
  {
    typeT type=cell.type();
    idT id=cell.id();

    if (type==0) return vertexVal[id];
    if (getByMax==-1) getByMax=ParentClass::getByMax();

    NDNET_UINT *vid = VERTEX_IN_FACE(net,type,id);
    NDNET_UINT max=vid[0];
    //printf("(%d %e %e)\n",ParentClass::getByMax(),vertexVal[vid[0]],vertexVal[vid[1]]);
    if (getByMax)
      {
	for (NDNET_UINT i=1;i<NUM_VERTEX_IN_FACE(net,type,id);i++)
	  if (vertexVal[vid[i]]>vertexVal[max]) max = vid[i];
      }
    else
      {
	for (NDNET_UINT i=1;i<NUM_VERTEX_IN_FACE(net,type,id);i++)
	  if (vertexVal[vid[i]]<vertexVal[max]) max = vid[i];
      }   
    return vertexVal[max];
  }

  unsigned long getNFaces(int type) const
  {
    if (type==0) return net->nvertex;
    else return net->nfaces[type];
  }

  void getCofaces(cellType cell, std::vector<cellType> &lst) const
  {
    typeT type=cell.type();
    idT id=cell.id();
    idT i,j,k,l;
    int N=type+1;
    int N2=N+1;
	 
    if (type==getNDims()) {lst.clear();return;}

    if (!net->haveFaceFromVertex[type+1]) 
      { 
#pragma omp critical
	{
	  int v=verbose;
	  verbose=0;
	  if (!net->haveFaceFromVertex[type+1]) 
	    {
	      printf("(updt %df",type+1);fflush(0);
	      ComputeFaceFromVertex(net,type+1);
	      printf(") ");fflush(0);
	    }
	  verbose=v;
	}
      }
	  
    if (type==0) {
      NDNET_UINT *fiv = FACE_IN_VERTEX(net,type+1,id);
      NDNET_UINT NF = NUM_FACE_IN_VERTEX(net,type+1,id);
      lst.resize(NF);
      for (i=0;i<NF;i++) lst[i]=cellType(type+1,fiv[i]);
	    
      return;
    }
	  
    // that's for inlining loops ...
    else if (type==1)
      {
	NDNET_UINT *vif = VERTEX_IN_FACE(net,type,id);
	int ndm1=getNDims()-1;
	lst.clear();
	i=vif[0];
	NDNET_UINT *fiv = FACE_IN_VERTEX(net,type+1,i);
	NDNET_UINT NF = NUM_FACE_IN_VERTEX(net,type+1,i);
	for (j=0;j<NF;j++) {
	  NDNET_UINT *cov = VERTEX_IN_FACE(net,type+1,fiv[j]);
	  if (cov[0]==vif[0])
	    {
	      if (cov[1]==vif[1]) 
		lst.push_back(cellType(type+1,fiv[j]));
	      else if (cov[2]==vif[1]) 
		lst.push_back(cellType(type+1,fiv[j]));
	    }
	  else if (cov[0]==vif[1])
	    {
	      if (cov[1]==vif[0]) 
		lst.push_back(cellType(type+1,fiv[j]));
	      else if (cov[2]==vif[0]) 
		lst.push_back(cellType(type+1,fiv[j]));
	    }
	  else
	    {
	      if (cov[1]==vif[0])
		{
		  if (cov[2]==vif[1]) 
		    lst.push_back(cellType(type+1,fiv[j]));
		}
	      else if (cov[1]==vif[1])
		{
		  if (cov[2]==vif[0]) 
		    lst.push_back(cellType(type+1,fiv[j]));
		}
	    }
	  if (ndm1==1)
	    if (lst.size()==2) break;
	}
	return;
      }
	  
    else if (type==2)
      {
	NDNET_UINT *vif = VERTEX_IN_FACE(net,2,id);
	int ndm1=getNDims()-1;
	i=vif[0];
	lst.clear();
	NDNET_UINT *fiv = FACE_IN_VERTEX(net,3,i);
	NDNET_UINT NF = NUM_FACE_IN_VERTEX(net,3,i);
	for (j=0;j<NF;j++) {
	  NDNET_UINT *cov = VERTEX_IN_FACE(net,3,fiv[j]);
	  NDNET_UINT nf=1;
	  for (k=0;k<4;k++) {
	    if (cov[k]==vif[1]) nf++;
	    else if (cov[k]==vif[2]) nf++;
	    if (nf==3) break;
	  }
	  if (k!=4) lst.push_back(cellType(3,fiv[j]));
	  if (ndm1==2)
	    if (lst.size()==2) break;
	}
	return;
      }
	  
    NDNET_UINT *vif = VERTEX_IN_FACE(net,type,id);
    int ndm1=getNDims()-1;
    i=vif[0];
    lst.clear();
    NDNET_UINT *fiv = FACE_IN_VERTEX(net,type+1,i);
    NDNET_UINT NF = NUM_FACE_IN_VERTEX(net,type+1,i);
    for (j=0;j<NF;j++) {
      NDNET_UINT *cov = VERTEX_IN_FACE(net,type+1,fiv[j]);
      NDNET_UINT nf=1;
      for (k=0;k<N2;k++) {
	if (cov[k]==i) continue;
	for (l=1;l<N;l++) 
	  {if (cov[k]==vif[l]) break;}

	if (l!=N) {if (++nf == N) break;}
      }
      if (k!=N2) lst.push_back(cellType(type+1,fiv[j]));
      if (type==ndm1)
	{
	  if (lst.size()==2) break;
	}
    }
	
  }
    
  void getFaces(cellType cell, std::vector<cellType> &lst) const
  {
    typeT type=cell.type();
    idT id=cell.id();
    idT i,j,k,l;
    int N=type+1;
    int N2=N-1;

    if (type==0) {lst.clear();return;}
	  
    if (type==1)
      {
	NDNET_UINT *vif = VERTEX_IN_FACE(net,type,id);
	lst.resize(2);
	lst[0]=cellType(0,vif[0]);
	lst[1]=cellType(0,vif[1]);
	return;
      }

    if (!net->haveFaceFromVertex[type-1]) 
      {
#pragma omp critical
	{
	  int v=verbose;
	  verbose=0;
	  if (!net->haveFaceFromVertex[type-1]) 
	    {
	      printf("(updt %df",type-1);fflush(0);
	      ComputeFaceFromVertex(net,type-1);
	      printf(") ");fflush(0);
	    }
	  verbose=v;
	}
      }

    NDNET_UINT *vif = VERTEX_IN_FACE(net,type,id);
    i=vif[0];
    lst.clear();
    NDNET_UINT *fiv = FACE_IN_VERTEX(net,type-1,i);
    NDNET_UINT NF = NUM_FACE_IN_VERTEX(net,type-1,i);

    for (j=0;j<NF;j++) {
      NDNET_UINT *cov = VERTEX_IN_FACE(net,type-1,fiv[j]);
      for (k=0;k<N2;k++) {
	if (cov[k]==i) continue;
	for (l=1;l<N;l++) 
	  {if (cov[k]==vif[l]) break;}
	if (l==N) break;
      }
      if (k==N2) lst.push_back(cellType(type-1,fiv[j]));
      if (lst.size()==N2) break;
    }
	  
    i=vif[1];
    fiv = FACE_IN_VERTEX(net,type-1,i);
    NF = NUM_FACE_IN_VERTEX(net,type-1,i);
    for (j=0;j<NF;j++) {
      NDNET_UINT *cov = VERTEX_IN_FACE(net,type-1,fiv[j]);
      for (k=0;k<N2;k++) {
	if (cov[k]==i) continue;
	if (cov[k]==vif[0]) break;
	for (l=2;l<N;l++) 
	  {if (cov[k]==vif[l]) break;}
	if (l==N) break;
      }
      if (k==N2) 
	{
	  lst.push_back(cellType(type-1,fiv[j]));
	  break;
	}
    }
	 
  }
    
  void getVertice(cellType cell, std::vector<cellType> &list) const
  {
    typeT type=cell.type();
    idT id=cell.id();
	
    NDNET_UINT *v = VERTEX_IN_FACE(net,type,id);
    NDNET_UINT N = NUM_VERTEX_IN_FACE(net,type,id);

    list.resize(N);
    int i;
    for (i=0;i<N;i++) list[i]=cellType(0,v[i]);
	  
  }

  float *getPosition(cellType cell,float *pos) const
  {
    typeT type=cell.type();
	  
    if (type==0)
      {
	idT id=cell.id();
	memcpy(pos,&net->v_coord[net->ndims*id],sizeof(float)*net->ndims);
	return pos;
      }
	    
    std::vector<cellType> vlist;

    getVertice(cell,vlist);
    memcpy(pos,&net->v_coord[net->ndims*vlist[0].id()],sizeof(float)*net->ndims);
    if (periodicity)
      {		
	std::vector<float> ref(pos,&pos[net->ndims]);
	int ni=0;
	for (int i=1;i<vlist.size();i++)
	  for (int j=0;j<net->ndims;j++) 
	    {
	      float newp = net->v_coord[net->ndims*vlist[i].id()+j];
	      pos[j]+=newp;
	      if (periodicity&(1<<j))
		{
		  if (fabs(ref[j]-newp)>net->delta[j]/2)
		    {
		      if (ref[j]>newp)
			pos[j]+=net->delta[j];
		      else
			{
			  pos[j]-=net->delta[j];				    
			}
		    }
		}
	    }
      }
    else
      {
	for (int i=1;i<vlist.size();i++)
	  for (int j=0;j<net->ndims;j++) 
	    pos[j]+=net->v_coord[net->ndims*vlist[i].id()+j];
      }
    double f=1./vlist.size();
    for (int j=0;j<net->ndims;j++) pos[j]*=f;

    return pos;
  }

  uint getNodeGroup(cellType cell) const
  {
    typeT type=cell.type();
    idT id=cell.id();

    if (vertexGroup == NULL) return 0;

    if (type==0) return (uint)vertexGroup[id];
	    
    std::vector<cellType> vlist;
    getVertice(cell,vlist);
	    
    double g=vertexGroup[vlist[0].id()];
    for (int i=1;i<vlist.size();i++)
      if (vertexGroup[vlist[i].id()]!=0)
	{
	  if (vertexGroup[vlist[i].id()]>g) g = vertexGroup[vlist[i].id()];
	}
	    
    return (uint)g;
  }

  bool groupsAreDefined() const
  {
    if (vertexGroup == NULL) return false;
    return true;
  }

  bool dumpSubComplex(const char *fname,subCellsType &subList,cellsInfoType &subListInfo, newCellsType &newList,cellsInfoType &newListInfo,int smooth, bool allowDuplicates) const
  {
    NDsubNetwork<ParentClass> sub(this);
    NDnetwork *refNet = net;
    char name[255];
    typename subCellsType::iterator sit;
    typename newCellsType::iterator nit;
    newCellsElementType *nce;

    sub.setRefNetwork(refNet);
    
    sub.insert(subList,subListInfo,false,!allowDuplicates);
    
    std::vector< std::pair<typeT,idT> > tmp;
    for (nit=newList.begin();nit!=newList.end();nit++)
      {	
	tmp.clear();
	for (typename newCellsElementType::iterator it=nit->begin(); it!=nit->end();it++)
	  tmp.push_back(it->getAsPair());
	sub.addNewCell(tmp.begin(),tmp.end(),false,!allowDuplicates);
      }
    
    tmp.clear();

    NDnetwork *subnet = sub.create();
    sprintf(name,"%s.NDnet",fname);
    if (smooth) NDnet_smooth(subnet, smooth);
      
    ndnet::IO::save(subnet,std::string(name));
    FreeNDnetwork(&subnet);
    
    return true;
  }

public:
  static bool checkFileType(const char *filename, bool simplicial=false) {   
    return ndnet::IO::canLoad(std::string(filename));
  }

  int getNDims(bool network=false) const 
  {
    if (!network) return net->ndims; 
    else return net->ndims_net;
  }

  void getBoundingBox(std::vector<double> &x0, std::vector<double> &delta) const 
  {
    x0.assign(net->x0,&(net->x0[net->ndims]));
    delta.assign(net->delta,&(net->delta[net->ndims]));
  }

  bool sendMessage(messageT m)
  {
    if (m==ParentClass::MSG_FACE_COFACE_NEEDED) {
      for (int i=0;i<net->ndims;i++)
	ComputeFaceFromVertex(net,net->ndims-i);
      return true;
    }

    if (m==ParentClass::MSG_FACE_COFACE_FREE) {
      for (int i=0;i<net->ndims;i++)
	{
	  if (net->haveFaceFromVertex[i])
	    {
	      free(net->v_faceIndex[i]);net->v_faceIndex[i]=NULL;
	      free(net->v_numFaceIndexCum[i]);net->v_numFaceIndexCum[i]=NULL;
	    }
	  net->haveFaceFromVertex[i]=0;
	  return true;
	}
    }

    return false;
  }

private:
  bool markBoundaries()
  {
    long i,j;
    char type;
    unsigned char *ff[net->ndims+1];
    long nb=0;
    bool hasOut=false;
    cellType curCell;
    std::vector<long> nBound(net->ndims+1,0);
    std::vector<long> nOut(net->ndims+1,0);
	    
    printf ("Tagging boundary ... ");fflush(0);
	  
    NDNetFlags_enable(net,0);
    ff[0]=net->v_flag;
    for (i=1;i<=net->ndims;i++) 
      {
	NDNetFlags_enable(net,i);
	ff[i]=net->f_flag[i];
      }	  
	    
    for (i=0;i<net->nvertex;i++) 
      {
	ff[0][i] &= (~NDNETFLAG_BOUNDARY);
	if (isOut(cellType(0,i))) hasOut=true;
      }
    for (type=1;type<=net->ndims;type++)
      for (i=0;i<net->nfaces[type];i++) 
	ff[type][i]&=(~NDNETFLAG_BOUNDARY);
	  
    if (hasOut)
      {
	for (type=1;type<=net->ndims;type++)
	  for (i=0;i<net->nfaces[type];i++)
	    {
	      static std::vector<cellType> vert;
	      getVertice(cellType(type,i),vert);
	      for (j=0;j<vert.size();j++) if (isOut(vert[j])) break;
	      if (j!=vert.size()) {
		net->f_flag[type][i]|=NDNETFLAG_OUT; 
		nOut[type]++;
	      }
	    }
      }
	

    type=net->ndims-1;
    for (i=0;i<net->nfaces[type];i++)
      {
	std::vector<cellType> coface;	
	curCell.set(type,i);
	      
	if (isOut(curCell)) continue;
	getCofaces(curCell,coface);
	      
	if ((coface.size()<2)||(isOut(coface[0])||isOut(coface[1]))) 
	  {
	    static std::vector<cellType> face;
		  
	    ff[type][i]|=NDNETFLAG_BOUNDARY;
	    nBound[type]++;
	    nb++;
		
	    getFaces(curCell,face);
	    for (typename std::vector<cellType>::iterator curf = face.begin(); curf!=face.end(); curf++)
	      {
		if (!(ff[type-1][curf->id()]&NDNETFLAG_BOUNDARY))
		  {
		    ff[type-1][curf->id()]|=NDNETFLAG_BOUNDARY;
		    nBound[type-1]++;
		  }
	      }
		  
	  }
      }
	 	  

    for (type=net->ndims-2;type>0;type--)
      {
	if (!nBound[type+1]) break;

	for (i=0;i<net->nfaces[type];i++)
	  {
	    curCell.set(type,i);
	    if (isBoundary(curCell))
	      {
		static std::vector<cellType> face;
		      
		getFaces(curCell,face);
		for (typename std::vector<cellType>::iterator curf = face.begin(); curf!=face.end(); curf++)
		  {
		    if (!(ff[type-1][curf->id()]&NDNETFLAG_BOUNDARY))
		      {
			ff[type-1][curf->id()]|=NDNETFLAG_BOUNDARY;
			nBound[type-1]++;
		      }
			  
		  }
	      }
	  }
      }
	 	
    printf("done.\n");
    for (type=0;type<=net->ndims;type++) printf("   [ %ld / %ld ] %d-faces on boundary / outside.\n",nBound[type],nOut[type],type);
    return true;
  }

public:

  NDnet_network(NDnetwork *netp,const char *valueField=VALUE_TAG, const char *groupField=GROUP_ID_PEAK_TAG, 
		double **val=NULL, bool isOwner=false,bool getValueByMax = true) : ParentClass(isOwner,getValueByMax)
  {	

    long i,j;
    net=netp;
    net->ndims_net=0;
    for (i=net->ndims;i>0;i--) 
      {
	if (net->haveVertexFromFace[i]) 
	  {
	    net->ndims_net=i;
	    break;
	  }
      }
	
    for (i=getNDims(true);i>0;i--) 
      {
	ComputeVertexFromFace(net,i);
	  
      }
	  
    for (i=0;i<getNDims(true);i++)
      {
	    
	ComputeFaceFromVertex(net,getNDims(true)-i);
	     
      }

	   
    strcpy(valueTag,valueField);
    i = NDDataIndex(net,0,valueTag);

    if ((i<0)&&(val!=NULL)) addNDDataArr(net,0,valueTag,val);
    if ((i>=0)&&(val!=NULL)) {
      printf("Replacing Value tag in network.\n");
      addNDDataArr(net,0,valueTag,val);
    }

    i = NDDataIndex(net,0,valueTag);

    if (i<0)
      {
	fprintf(stderr,"Error in NDnet_network: field '%s' must be present in NDnetwork.\n",valueTag);
	fprintf(stderr,"  I need a scalar function to work with!\n  Try tagging the network or using option '-field'.\n");
	exit(-1);
      }
    else vertexVal = net->data[i].data; 
	    
    strcpy(groupTag,groupField);
    i = NDDataIndex(net,0,groupTag);
    if (i>=0) vertexGroup = net->data[i].data; 	    
    else vertexGroup = NULL; 
    periodicity = net->periodicity;	    

    std::vector<char> dummy;
    setMask(dummy,true); 
  
    if (debug_dump) 
      {
	double *data[this->getNDims(true)+1];
	std::vector<long> ct(net->ndims+1);

	for (i=0;i<=this->getNDims(true);i++)
	  data[i]=addNDDataArr(net,i,BOUNDARY_TAG,NULL);
		
	for (i=0;i<net->nvertex;i++)
	  {
	    if ((net->v_flag[i])&NDNETFLAG_BOUNDARY) {data[0][i]=1;ct[0]++;}
	    else if ((net->v_flag[i])&NDNETFLAG_OUT) data[0][i]=2;
	    else if ((net->v_flag[i])&NDNETFLAG_NON_MANIFOLD) data[0][i]=3;
	    else data[0][i]=0;
	  }
		
	for (j=1;j<=this->getNDims(true);j++)
	  for (i=0;i<net->nfaces[j];i++)
	    {
	      if ((net->f_flag[j][i])&NDNETFLAG_BOUNDARY) {data[j][i]=1;ct[j]++;}
	      else if ((net->f_flag[j][i])&NDNETFLAG_OUT) data[j][i]=2;
	      else if ((net->f_flag[j][i])&NDNETFLAG_NON_MANIFOLD) data[0][i]=3;
	      else data[j][i]=0;
	    }

	//std::copy(ct.begin(), ct.end(), std::ostream_iterator<long>(std::cout, "\n"));
		
	ndnet::IO::save(net,std::string("tempNET.NDnet"));
      }
    //init();
  }
    
  ~NDnet_network()
  {
    if (ParentClass::netIsOwned()) freeData();
  }

};

#endif
