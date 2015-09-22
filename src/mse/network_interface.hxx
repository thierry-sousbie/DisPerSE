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
#ifndef __NETWORK_INTERFACE__HXX__
#define __NETWORK_INTERFACE__HXX__

#include "mytypes.h"
#include "cells.hxx"
#include <vector>

//template <typename TypeT = char, typename IdT = uint>
template <class cellType = cellIdentity<> >
class NetworkDataInterface {

private:
  bool isOwned;
  bool getValByMax;

protected:
  bool getByMax() const {return getValByMax;}
  bool netIsOwned() const {return isOwned;}

public:
  
  enum messageT {MSG_FACE_COFACE_NEEDED=1, MSG_FACE_COFACE_FREE=2, MSG_SIMPLICIAL_COMPLEX=3};

  /*
  typedef char TypeT;
  typedef uint IdT;
  typedef std::pair<TypeT,IdT> cellID;
  */
  typedef cellType cellT;
  typedef cellT subCellsElementType;
  typedef std::vector< cellT > subCellsType;  
  typedef std::vector< cellT > newCellsElementType;
  typedef std::vector< newCellsElementType > newCellsType;

  typedef std::pair< std::vector<std::string>,std::vector<double> > cellsInfoType;

public:
  virtual int getNDims(bool network=false) const=0;
  virtual void getBoundingBox(std::vector<double> &x0, std::vector<double> &delta) const=0;
  virtual bool sendMessage(messageT m)=0;

public:
  virtual void freeData()=0;
  virtual int getPeriodicity() const=0;
  virtual int setPeriodicity(int p) =0;
  virtual void setMask(std::vector<char> &mask, bool nullIsMasked=false) =0;
  virtual bool isOut(cellType cell) const=0;
  virtual bool isBoundary(cellType cell) const=0;
  virtual double getValue(cellType cell, int getByMax=-1) const=0;
  virtual double getValueAll(cellType cell,  std::vector<double> &values, int getByMax=-1) const=0;
  virtual double getValueSum(cellType cell, double &sum, int getByMax=-1) const=0;
  virtual unsigned long getNFaces(int type) const=0;
  
  virtual void getCofaces(cellType cell, std::vector<cellType> &list) const =0;
  virtual void getFaces(cellType cell, std::vector<cellType> &list) const =0;
  virtual void getVertice(cellType cell, std::vector<cellType> &list) const =0;
  virtual float *getPosition(cellType cell, float *pos) const =0;  
  virtual uint getNodeGroup(cellType cell) const=0;
  virtual bool groupsAreDefined() const=0;
  
  //virtual bool dumpSubComplex(const char *fname,subCellsType &subList,subCellsInfoType &subListInfo,newCellsType &newList, int smooth=0, bool allowDuplicates=false) const=0;
  virtual bool dumpSubComplex(const char *fname,subCellsType &subList,cellsInfoType &subListInfo, newCellsType &newList,cellsInfoType &newListInfo,int smooth, bool allowDuplicates) const=0;

public:

  void setGetValueByMax(bool val)
  {
    getValByMax=val;
  }
     
  NetworkDataInterface(bool isOwner = false, bool getValueByMax = true)
  {
    if (isOwner) isOwned=true;
    else isOwned=false;
    setGetValueByMax(getValueByMax);
  }

  virtual ~NetworkDataInterface() 
  {
    
  }

};

#endif
