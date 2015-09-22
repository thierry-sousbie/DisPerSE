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
#ifndef __SIMPLICIALGRID_MSC_HEADER__
#define __SIMPLICIALGRID_MSC_HEADER__

#include "network_interface.hxx"

#include "MSComplex.hxx"
#include "mytypes.h"
#include "NDfield.h"
#include "healpix.h"
#include <stdio.h>
#include <string.h>

#include "NDsubNetwork.hxx"
#include "simplicialGrid.hxx"
#include "gridType.hxx"

#include "NDnet_interface.hxx"

//template <typename TypeT = char, typename IdT = uint>
template <class cellType = cellIdentity<> >
class SimplicialGrid_network : public NetworkDataInterface<cellType> {

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

  typedef simplicialGrid<cellType> gridT;

  //NDfield *net;
  gridType<double> *net;
  
  double *val;
  std::vector<long> nFaces;

  long id2index[1<<8];
  long index2id[2][4][10];
  int periodicity;

  gridT grid;
  
  std::vector< std::vector<bool> > mask_out;
  std::vector< std::vector<bool> > mask_boundary;
  
  void init()
  {
    nFaces=grid.getNFaces();
  
  }

  long setBoundaries()
  {
    cellType curCell;
    char type=this->getNDims()-1;
    long nf=this->getNFaces(type);
    long i;
    long nb=0;

    mask_boundary.resize(net->ndims+1);
    for (i=0;i<=this->getNDims();i++) 
      mask_boundary[i].assign(this->getNFaces(i),false);

    for (i=0;i<nf;i++)
      {
	std::vector<cellType> coface;	
	curCell.set(type,i);

	if (isOut(curCell)) continue;
	getCofaces(curCell,coface);
	
	if ((coface.size()!=2)||(isOut(coface[0])||isOut(coface[1]))) 
	  {
	    static std::vector<cellType> face;
	    
	    mask_boundary[type][i]=true;
	    nb++;
	    
	    getFaces(curCell,face);
	    for (typename std::vector<cellType>::iterator curf = face.begin(); curf!=face.end(); curf++)
	      {
		mask_boundary[type-1][curf->id()]=true;
	      }
	  }
      }
    
    for (type=this->getNDims()-2;type>0;type--)
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
		    mask_boundary[type-1][curf->id()]=true;
		  }
	      }
	  }
      }

    return nb;
  }


public:
  int getPeriodicity() const 
  {
    return periodicity;
  } 

  int setPeriodicity(int p) 
  {
    std::vector<char> dummy;

    periodicity=p;
    setMask(dummy,false);
    
    return p;
  }
  
  void freeData() 
  {
    // Free all allocated ressources here
  }
    
  void setMask(std::vector<char> &mask, bool nullIsMasked)
  {
    long i,j,k;
    typeT type;    
    cellType curCell;
    bool fullPeriodic=true;
   
    if (mask.size())
      {
	if (mask.size() != net->nval)
	  {
	    fprintf(stderr,"ERROR in setMask: mask and data do not have the same size.\n");
	    exit(-1);
	  }
      }
    else 
      {
	for (i=0;i<net->ndims;i++)
	  if (!(periodicity&(1<<i))) break;
	if (i==net->ndims) return;
      }

    for (i=0;i<getNDims();i++)
      if (!(periodicity&(1<<i))) 
	{
	  fullPeriodic=false;
	  break;
	}

    printf ("Building mask ... ");fflush(0);

    if (!mask_out.size())
      {
	mask_out.resize(net->ndims+1);
	for (i=0;i<net->ndims+1;i++)
	  mask_out[i].assign(this->getNFaces(i),false);	
      }
     
    if (mask.size())
      {
	if (nullIsMasked)
	  for (i=0;i<mask.size();i++) mask_out[0][i] = mask_out[0][i] | (!(bool)mask[i]);
	else
	  for (i=0;i<mask.size();i++) mask_out[0][i] = mask_out[0][i] | (bool)mask[i];
      }

    for (type=1;type<=net->ndims;type++) {
      long nf = this->getNFaces(type);
      for (i=0;i<nf;i++)
	{
	  static std::vector<cellType> vert;
	  if (mask_out[type][i]) continue;
	  getVertice(cellType(type,i),vert);
	  for (j=0;j<vert.size();j++) if (isOut(vert[j])) break;
	  if (j!=vert.size()) mask_out[type][i]=true; 
	  else if (!fullPeriodic)
	    {
	      std::vector<long> coords(net->ndims);
	      std::vector<double> mn(net->ndims,DBL_MAX);
	      std::vector<double> mx(net->ndims,-DBL_MAX);
	    
	      for (j=0;j<vert.size();j++)
		{
		  grid.index2coords(vert[j].id(),coords);
		  for (k=0;k<net->ndims;k++)
		    {		      
		      if (coords[k]<mn[k]) mn[k]=coords[k];
		      if (coords[k]>mx[k]) mx[k]=coords[k];
		    }
		}
	      for (k=0;k<net->ndims;k++)
		{
		  if (periodicity&(1<<k)) continue;
		  if ((mx[k]-mn[k])>1) break;
		}
	      if (k!=net->ndims) mask_out[type][i]=true; 
	    }
	}
    }

    long nb = setBoundaries();
    printf ("done. (%d %d-faces on boundary)\n",(int)nb,net->ndims-1);

  }

  bool isOut(cellType cell) const 
  {
    typeT type=cell.type();
    idT id=cell.id();
    
    // return true if the cell is on the boundary or too close to be trusted
    if (!mask_out.size()) return false;
    return mask_out[type][id];
  }
  
  bool isBoundary(cellType cell) const 
  {
    typeT type=cell.type();
    idT id=cell.id();
    
    // return true if the cell is on the boundary or too close to be trusted
    if (!mask_boundary.size()) return false;
    return mask_boundary[type][id];
  }
 
  double getValueSum( cellType cell, double &sum, int getByMax)  const
  {
    typeT type=cell.type();
    idT id=cell.id();

    if (type==0) {sum=val[id];return val[id];}
    if (getByMax==-1) getByMax=ParentClass::getByMax();

    std::vector<cellType> r;
    getVertice(cell,r);
    double value = val[r[0].id()];
    
    int mn=0;
    int mx=0; 
    double vn=value;
    double vx=value;
    
    long i;
    sum=value;

    if (getByMax)
      {
	for (i=1;i<r.size();i++) 
	  {
	    double v=val[r[i].id()];
	    sum+=v;
	    if (v>value) {value=vx=v;mx=i;}
	 
	  }
      }
    else
      {
	for (i=1;i<r.size();i++) 
	  {
	    double v=val[r[i].id()];
	    sum+=v;
	    if (v<value) {value=vn=v;mn=i;}
	 
	  }
      }
    return value;

    float p1[6];
    float p2[6];
    double dist=0;

    getPosition(r[mn],p1);
    getPosition(r[mx],p2);
    
    for (i=0;i<net->ndims;i++) 
      {
	p1[i]-=p2[i];
	dist+=p1[i]*p1[i];
      }
    
    sum=(vx-vn)/sqrt(dist);

    return value;
  }
  
  double getValueAll( cellType cell, std::vector<double> &values, int getByMax)  const
  {
    // return the value of the cell, sum being the sum of the values for each vertex
    // sum is used to compare two different cells with identical type and value
    if (getByMax==-1) getByMax=ParentClass::getByMax();
  }
  
  double getValue( cellType cell, int getByMax) const
  {
    typeT type=cell.type();
    idT id=cell.id();

    if (type==0) return val[id];    
    if (getByMax==-1) getByMax=ParentClass::getByMax();

    std::vector<cellType> r;
    getVertice(cell,r);
    double value = val[r[0].id()];
    long i;
    if (getByMax)
      {
	for (i=1;i<r.size();i++) 
	  {
	    double v=val[r[i].id()];
	    if (v>value) value=v;
	  }
      }
    else
      {
	for (i=1;i<r.size();i++) 
	  {
	    double v=val[r[i].id()];
	    if (v<value) value=v;
	  }
      }

    return value;      
  }

  unsigned long getNFaces(int type) const
  {
    return nFaces[type];
  }

  void getCofaces(cellType cell, std::vector<cellType> &result) const
  {
    return grid.getCofaces(cell,result);
  }
  
  void getFaces(cellType cell, std::vector<cellType> &result) const
  {
    return grid.getFaces(cell,result);
   }
    
  void getVertice(cellType cell, std::vector<cellType> &result) const
  {
    return grid.getVertice(cell,result);
  }
  
  float *getPosition(cellType cell,float *pos) const
  {
    grid.getPosition(cell,pos);
    //if (pos[1]<-1.5) {printf("(%d;%ld)\n",cell.type(),cell.id());exit(0);}
  }

  uint getNodeGroup(cellType cell) const
  {
    return 0;
  }
  
  bool groupsAreDefined() const
  {
    return false;
  }
  
  bool dumpSubComplex(const char *fname,subCellsType &subList,cellsInfoType &subListInfo, newCellsType &newList,cellsInfoType &newListInfo,int smooth, bool allowDuplicates) const
  {
    NDsubNetwork<ParentClass> sub(this);
    NDnetwork *refNet = NULL;
    char name[255];
    typename subCellsType::iterator sit;
    typename newCellsType::iterator nit;
    newCellsElementType *nce;
    
    sub.setRefNetwork(refNet);
   
    sub.insert(subList,subListInfo,true,!allowDuplicates);

    std::vector< std::pair<typeT,idT> > tmp;
    for (nit=newList.begin();nit!=newList.end();nit++)
      {	
	tmp.clear();
	for (typename newCellsElementType::iterator it=nit->begin(); it!=nit->end();it++)
	  tmp.push_back(it->getAsPair());
	sub.addNewCell(tmp.begin(),tmp.end(),true,!allowDuplicates);
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
  static bool checkFileType(const char *filename,bool simplicial=false) {
    
    if (gridType<double>::isLoadable(filename))
      if (!is_healpix_map(filename)) return true;
    return false;

    return false;
    
  }

   static NetworkDataInterface<cellType> *Load(const char *filename) {
     double *data;
     float *dataf;
     long i;
     gridType<double> *field= new gridType<double>(filename);

     long nout=0;
     std::vector<char> out;
     for (i=0;i<field->nval;i++)
       if (field->val[i]!=field->val[i])
	 {
	   if (!nout) out.assign(field->nval,0);
	   out[i]=1;nout++;
	 }
     
   
    if (nout) printf ("Found %ld NaN values, a mask will be generated.\n",nout);fflush(0);
    NetworkDataInterface<cellType> *net = new SimplicialGrid_network(field,true);
    if (nout) net->setMask(out);
       
    return net;
   }

  int getNDims(bool network=false) const 
  {
    return net->ndims;
    // returns the number of dimensions
  }

  void getBoundingBox(std::vector<double> &x0, std::vector<double> &delta) const 
  {
    // set x0 and delta to the origin and size of the cubic bounding box  
    x0=net->x0;
    delta=net->delta;
  }
  
  bool sendMessage(messageT m)
  {
    
    return false;
  }

public:
  
  SimplicialGrid_network(gridType<double> *netp,bool isOwner=false,bool getValueByMax = true) : ParentClass(isOwner,getValueByMax)
  {		    
    long i;
    long du,dv,dw;
    long p1,p2,p12;

    net=netp;
   val=net->val;
    
    periodicity=~(int)0;
  
    printf ("Generating implicit simplicial complex ... ");fflush(0);
 
    for (i=0;i<net->ndims;i++)
      if (net->dims[i]&1) break;
    
    if (net->ndims==1)
      {
	grid.setType(gridT::tesselation::any1D);
      }
    else if (i!=net->ndims)
      {
	if (net->ndims==2)
	  grid.setType(gridT::tesselation::regular2D);
	else if (net->ndims==3)
	  grid.setType(gridT::tesselation::regular3D);
      }
    else
      {
	if (net->ndims==2)
	  grid.setType(gridT::tesselation::alternate2D);
	else if (net->ndims==3)
	  grid.setType(gridT::tesselation::alternate3D);
      }
       
    grid.setGrid(net->x0,net->delta,net->dims);

    // ... and don't forget to call init() once set
    init();    
    printf("done.\n");

    std::vector<char> dummy;
    setMask(dummy,true);
   
  }
  
  ~SimplicialGrid_network()
  {
    if (ParentClass::netIsOwned()) freeData();
  }
  
};


#endif
