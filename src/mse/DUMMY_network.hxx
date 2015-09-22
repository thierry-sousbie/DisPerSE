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
#ifndef __DUMMY_MSC_HEADER__
#define __DUMMY_MSC_HEADER__

// Use this as an interface between your own datatype and MSComplex class.
// Change DUMMY to whetever your datatype is and implement necessary functions.
#include "network_interface.hxx"
#include "NDnet_interface.hxx"

#include "MSComplex.hxx"
#include "mytypes.h"
#include "NDfield.h"

#include <stdio.h>
#include <string.h>

class dummyDataT
{
 
};


//template <typename TypeT = char, typename IdT = uint>
template <class cellType = cellIdentity<> >
class Dummy_network : public NetworkDataInterface<cellType> {

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
  //Dummy *net; // your datatype

public:
  int getPeriodicity() const 
  {
  } 

  int setPeriodicity(int p) 
  {
  }

  void freeData() 
  {
    // Free all allocated ressources here
  }
    
  void setMask(std::vector<char> &mask, bool nullIsMasked)
  {
    
  }

  bool isOut(cellType cell) const 
  {
    // return true if the cell is outside the valid domain
  }

  bool isBoundary(cellType cell) const 
  {

  }

  double getValueSum( cellType cell, double &sum, int getByMax)  const
  {
    // return the value of the cell, sum being the sum of the values for each vertex
    // sum is used to compare two different cells with identical type and value
    if (getByMax==-1) getByMax=ParentClass::getByMax();
  }

  double getValueAll( cellType cell, std::vector<double> &values, int getByMax)  const
  {
    // return the value of the cell, sum being the sum of the values for each vertex
    // sum is used to compare two different cells with identical type and value
    if (getByMax==-1) getByMax=ParentClass::getByMax();
  }

  double getValue( cellType cell, int getByMax) const
  {
    // return the value of the cell
    if (getByMax==-1) getByMax=ParentClass::getByMax();
  }

  unsigned long getNFaces(int type) const
  {
    // returns the total number of faces of a given type 
  }

  void getCofaces(cellType cell, std::vector<cellType> &result) const
  {
    // set result to the IDs of the cofaces of cell(type,id)
    // A=cell(type,id) is a coface of B=cell(type-1,coID) if B is a face of A.
  }
    
  void getFaces(cellType cell, std::vector<cellType> &result) const
  {
    // set result to the IDs of the faces of cell(type,id)
  }
    
  void getVertice(cellType cell, std::vector<cellType> &result) const
  {
    // set result the the IDs of the vertice of cell(type,id)
  }

  float *getPosition(cellType cell,float *pos) const
  {
    // set pos to the position of cell(type,id)
    // no need to take care of memory allocation for pos
  }

  uint getNodeGroup(cellType cell) const
  {
    // returns the group of cell(type,id)
    // 0 means no group
    // used to simplify by groups only
  }

  bool groupsAreDefined() const
  {
    // returns true if groups are defined for cells, false otherwise.
    // if "false", getNodeGroup(cellType cell) may return anything.
  }

  
  bool dumpSubComplex(const char *fname,subCellsType &subList,cellsInfoType &subListInfo, newCellsType &newList,cellsInfoType &newListInfo,int smooth, bool allowDuplicates) const
  {
    // this dumps a sucomplex of the initial network to a file
   
    //return false;

    // This implementation will work for any simplicial complex
    
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
  static bool checkFileType(const char *filename) {
    //return IsDummy(filename);
    return false;
  }

  // This is optional ...
  static NetworkDataInterface<cellType> *Load(const char *filename) {
    //return NetworkDataInterface<cellType> *net = new Dummy_network(Load_DUMMY(filename),true);
  }

  int getNDims(bool network) const 
  {
    // returns the number of dimensions
  }
  void getBoundingBox(std::vector<double> &x0, std::vector<double> &delta) const 
  {
    // set x0 and delta to the origin and size of the cubic bounding box
  }
  
  bool sendMessage(messageT m)
  {
    
    return false;
  }

public:

  Dummy_network(dummyDataT *netp,bool isOwner=false,bool getValueByMax = true) : ParentClass(isOwner,getValueByMax)
  {		    
    // set isOwner to true if netp can be freed/modified
    // nthreads is the number of threads
	    
    // All the setup goes here ...

    //net=netp;

	    

	   
  }
    
  ~Dummy_network()
  {
    if (ParentClass::netIsOwned()) freeData();
  }

};


#endif
