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
#ifndef __HEALPIX_MSC_HEADER__
#define __HEALPIX_MSC_HEADER__

#include "network_interface.hxx"


#include "MSComplex.hxx"
#include "mytypes.h"

#include <stdio.h>
#include <string.h>

#include <algorithm>
#include <vector>

#include "healpix.h"

#include "NDnet_interface.hxx"

template <class cellType = cellIdentity<> >
class Healpix_network : public NetworkDataInterface<cellType> {

  typedef class NetworkDataInterface<cellType> ParentClass;
  typedef typename ParentClass::messageT messageT;
  typedef typename cellType::typeT typeT;
  typedef typename cellType::idT idT;
  
  typedef typename ParentClass::subCellsElementType subCellsElementType;
  typedef typename ParentClass::subCellsType subCellsType;
  typedef typename ParentClass::cellsInfoType cellsInfoType;
  typedef typename ParentClass::newCellsElementType newCellsElementType;
  typedef typename ParentClass::newCellsType newCellsType;

  struct healpixMap
  {
    double *hmap;
    long nside;
    char coordsys[10];
    char ordering[10];
  };

private:
  
  template <class T>
  class quad {
  public :
    T a,b,c,d;

    quad(T aa, T bb, T cc, T dd)
    {
      a=aa;b=bb;c=cc;d=dd;
    }
  };

  template <class T>
  class triplet {
  public :
    T a,b,c;

    triplet(T aa, T bb, T cc)
    {
      a=aa;b=bb;c=cc;
    }
  };

  double *hmap;
  long nside;
  long npix;
  std::vector<long> nfaces;

  std::vector< std::vector<int> > segTab;
  std::vector<int> segId;

  //std::vector< quad<int> > faceTab;
  std::vector< triplet<int> > faceTab;


  std::vector<cellType> seg2faceA; 
  std::vector<cellType> seg2faceB;  

  bool nested;
  int diags;

  void buildFacesTables()
  {
    long i,j;
    int nnei;
    int res[8];
    
    printf("* Buiding lookup tables for %ld pixels ... ",npix);fflush(0);

    nfaces.resize(3);
    nfaces[0]=npix;
    segTab.clear();
    segTab.resize(npix);   
    for (i=0;i<npix;i++)
      {
	bool wasHigh[2];
	int skipped;
	int res_nod=0;
	wasHigh[0]=wasHigh[1]=false;

	nnei = pixNeighbours(i,nside,nested,true,res,&res_nod);
	
	//if (i<2) {printf("\npix %ld: ",i);for (j=0;j<nnei;j++) printf("%d(%d) ",res[j],(res_nod&(1<<j))?1:0);printf("\n");}

	for (j=0;j<nnei+2;j++)
	  {
	    int jj=j%nnei;
	    int jjm1=(j-1)%nnei;
	    if (res[jj]>i) 
	      {	
		if (j<nnei) 
		  {
		    if (res_nod&(1<<j)) 
		      segTab[i].push_back(res[j]);
		  }

		if (res_nod&(1<<(jj))) 
		  {
		  
		    if ((nnei<8)&&(wasHigh[0])&& (res_nod&(1<<jjm1)) )
		      {
			if (j<nnei+1)
			  {
			    faceTab.push_back(triplet<int>(i,res[jjm1],res[jj]));
			  }
		      }
		    else if (wasHigh[0]&&wasHigh[1]) 
		      {
		      
		     
			faceTab.push_back(triplet<int>(i,res[(j-2)],res[jjm1]));
			faceTab.push_back(triplet<int>(i,res[jjm1],res[jj]));
		
		      
			if (res[jjm1]>i) segTab[i].push_back(res[jjm1]);
			else segTab[res[jjm1]].push_back(i);
		      
		      }	 
		  }
		wasHigh[1]=wasHigh[0];
		wasHigh[0]=true;
	      } 
	    else 
	      {
		wasHigh[1]=wasHigh[0];
		wasHigh[0]=false;
	      }
	  }
      }

    segId.resize(npix+1);
    segId[0]=0;
    for (i=0;i<npix;i++)
      segId[i+1]=segId[i]+segTab[i].size();

    nfaces[1]=segId[npix];
    nfaces[2]=faceTab.size();

    seg2faceA.resize(nfaces[1],cellType::UNSET);
    seg2faceB.resize(nfaces[1],cellType::UNSET);
    j=0;
    for (i=0;i<nfaces[2];i++)
      {
	std::vector<cellType> lst;
	cellType cell(2,i);
	long id;

	getFaces(cell,lst);
	
	
	id=lst[0].id();
	if (seg2faceA[id] == cellType::UNSET) 
	  seg2faceA[id] = cell;
	else if (seg2faceB[id] == cellType::UNSET) seg2faceB[id] = cell;
	//else nnei/=j;

	id=lst[1].id();
	if (seg2faceA[id] == cellType::UNSET) 
	  seg2faceA[id] = cell;
	else if (seg2faceB[id] == cellType::UNSET) seg2faceB[id] = cell;
	//else nnei/=j;

	id=lst[2].id();
	if (seg2faceA[id] == cellType::UNSET) 
	  seg2faceA[id] = cell;
	else if (seg2faceB[id] == cellType::UNSET) seg2faceB[id] = cell;
	//else nnei/=j;

	if (lst.size()==4)
	  {
	    id=lst[3].id();
	    if (seg2faceA[id] == cellType::UNSET) 
	      seg2faceA[id] = cell;
	    else if (seg2faceB[id] == cellType::UNSET) seg2faceB[id] = cell;
	    //else nnei/=j;
	  }
      }
    printf("done. (%ld edges, %ld faces)\n ",nfaces[1],nfaces[2]);fflush(0);
   

  }

  
 
public:
  int getPeriodicity() const 
  {
    return ~((int)0);
  } 

  int setPeriodicity(int p) 
  {
    return 0;
  }

  void freeData() 
  {
    // Free all allocated ressources here
    free(hmap);
  }
    
  void setMask(std::vector<char> &mask, bool nullIsMasked)
  {
    printf("ERROR : Masks not implemented on healpix yet !");
  }

  bool isOut(cellType cell) const 
  {
    return false; // implement masks ...
    // return true if the cell is on the boundary or too close to be trusted
  }

  bool isBoundary(cellType cell) const 
  {
    return false; //implement masks
  }

  double getValueSum( cellType cell, double &sum, int getByMax)  const
  {
    // return the value of the cell, sum being the sum of the values for each vertex
    // sum is used to compare two different cells with identical type and value
    typeT type=cell.type();
    idT id=cell.id();	 
	  
    // return the value of the cell
    if (type==0) {
      sum=hmap[id];
      return hmap[id];
    }
    if (getByMax==-1) getByMax=ParentClass::getByMax();

    std::vector<cellType> lst;
    getVertice(cell,lst);

    if (type==1)
      {
	idT a = lst[0].id();
	idT b = lst[1].id();
	sum=hmap[a]+hmap[b];
	if (getByMax)
	  return (hmap[a]>hmap[b])?hmap[a]:hmap[b];
	else
	  return (hmap[a]<hmap[b])?hmap[a]:hmap[b];
      }

    int i;
    int max=lst[0].id();
    sum=hmap[max];
    if (getByMax)
      {
	for (i=1;i<lst.size();i++)
	  {
	    idT j=lst[i].id();
	    sum+=hmap[j];
	    if (hmap[j]>hmap[max]) max = j;
	  }
      }
    else
      {
	for (i=1;i<lst.size();i++)
	  {
	    idT j=lst[i].id();
	    sum+=hmap[j];
	    if (hmap[j]<hmap[max]) max = j;
	  }
      }
	  
    return hmap[max];
  }

  double getValueAll( cellType cell, std::vector<double> &values, int getByMax)  const
  {
    // return the value of the cell, sum being the sum of the values for each vertex
    // sum is used to compare two different cells with identical type and value
    typeT type=cell.type();
    idT id=cell.id();	 
	  
    // return the value of the cell
    if (type==0) {
      values.resize(1);
      values[0]=hmap[id];
      return hmap[id];
    }
    if (getByMax==-1) getByMax=ParentClass::getByMax();

    std::vector<cellType> lst;
    getVertice(cell,lst);

    if (type==1)
      {
	idT a = lst[0].id();
	idT b = lst[1].id();	      
	//sum=hmap[a]+hmap[b];
	values.resize(2);
	if (getByMax == (hmap[a]>hmap[b]) )
	  {
	    values[0]=hmap[a];values[1]=hmap[b];
	    return hmap[a];
	  }
	else
	  {
	    values[1]=hmap[a];values[0]=hmap[b];
	    return hmap[b];
	  }
      }

    int i;
    int max=lst[0].id();
    values.resize(lst.size());
    values[0]=hmap[max];
    if (getByMax)
      {
	for (i=1;i<lst.size();i++)
	  {
	    idT j=lst[i].id();
	    values[i]=hmap[j];
	    if (hmap[j]>hmap[lst[max].id()]) max = i;
	  }
      }
    else
      {
	for (i=1;i<lst.size();i++)
	  {
	    idT j=lst[i].id();
	    values[i]=hmap[j];
	    if (hmap[j]<hmap[lst[max].id()]) max = i;
	  }
      }

    std::swap(values[0],values[max]);
	  
    return hmap[max];
  }



  double getValue( cellType cell, int getByMax) const
  {
    typeT type=cell.type();
	  
    // return the value of the cell
    if (type==0) return hmap[cell.id()];
    if (getByMax==-1) getByMax=ParentClass::getByMax();

    std::vector<cellType> lst;
    getVertice(cell,lst);

    if (type==1)
      {
	idT a = lst[0].id();
	idT b = lst[1].id();
	if (getByMax)
	  return (hmap[a]>hmap[b])?hmap[a]:hmap[b];
	else
	  return (hmap[a]<hmap[b])?hmap[a]:hmap[b];
      }

    int i;
    int max=lst[0].id();
    if (getByMax)
      {
	for (i=1;i<lst.size();i++)
	  {
	    idT j=lst[i].id();
	    if (hmap[j]>hmap[max]) max = j;
	  }
      }
    else
      {
	for (i=1;i<lst.size();i++)
	  {
	    idT j=lst[i].id();
	    if (hmap[j]<hmap[max]) max = j;
	  }
      }
	  
    return hmap[max];
  }

  unsigned long getNFaces(int type) const
  {
    return nfaces[type];
    // returns the total number of faces of a given type 
  }

  void getCofaces(cellType cell, std::vector<cellType> &lst) const
  {
	  
    typeT type=cell.type();
  
    if (type>1) return;

    idT id=cell.id();
    int i,j;
    int res[8];

    if (type==0) {
      int skipped;
	    
      int nnei = pixNeighbours(id,nside,nested,true,res,&skipped);
      lst.clear();
	    
      for (j=0;j<nnei;j++)
	{
	  if (id<res[j])
	    {
	      for (i=0;i<segTab[id].size();i++) 
		if (segTab[id][i] == res[j]) break;
	      //if (i>=segTab[id].size()) {printf("\n%d nei\n",nnei);exit(0);}
	      if (i<segTab[id].size())
		lst.push_back(cellType(1,segId[id]+i));
	    }
	  else 
	    {
	      for (i=0;i<segTab[res[j]].size();i++) 
		if (segTab[res[j]][i] == id) break;
	      //if (i>=segTab[res[j]].size()) {printf("\n%d nei\n",nnei);exit(0);}
	      if (i<segTab[res[j]].size())
		lst.push_back(cellType(1,segId[res[j]]+i));
	    }
	}
	    
    }
    else {
      lst.resize(2);
      lst[0]=seg2faceA[cell.id()];
      lst[1]=seg2faceB[cell.id()];
    }
	    
  }
    
  void getFaces(cellType cell, std::vector<cellType> &lst) const
  {
	  
    typeT type=cell.type();
    
    if (type<1) return;
	  
    idT id=cell.id();	  
    std::vector<int>::const_iterator low_it;

    if (type==1) {
      int a,b;   
      low_it = std::lower_bound(segId.begin(),segId.end(),id);
      a=(long)(low_it - segId.begin());
      if (segId[a] != id) a--;
      int delta= id-segId[a];
      if (delta>=segTab[a].size())
	{
	  do {
	    delta-=segTab[a].size();
	    a++;
	  } while (delta>=segTab[a].size());
	}
      b=segTab[a][delta];

      lst.resize(2);
      lst[0]=cellType(0,a);
      lst[1]=cellType(0,b);
    }
    else {
      std::vector<cellType> vlst;
      int a,b;
      int v0,v1,v2;
      int i;
      //int j=0;

      getVertice(cell,vlst);
      lst.resize(vlst.size());
      v0=vlst[0].id();v1=vlst[1].id();v2=vlst[2].id();

      if (v0<v1) {a=v0;b=v1;}
      else {a=v1;b=v0;}
      for (i=0;i<segTab[a].size();i++) if (segTab[a][i]==b) break;
      //if (i==segTab[a].size()) i/=j;
      lst[0]=cellType(1,segId[a]+i);

      if (v1<v2) {a=v1;b=v2;}
      else {a=v2;b=v1;}
      for (i=0;i<segTab[a].size();i++) if (segTab[a][i]==b) break;
      //if (i==segTab[a].size()) i/=j;
      lst[1]=cellType(1,segId[a]+i);

      if (vlst.size()==3) {
	if (v0<v2) {a=v0;b=v2;}
	else {a=v2;b=v0;}
	for (i=0;i<segTab[a].size();i++) if (segTab[a][i]==b) break;
	//if (i==segTab[a].size()) i/=j;
	lst[2]=cellType(1,segId[a]+i);
      }
      else {
	int v3=vlst[3].id();
	if (v2<v3) {a=v2;b=v3;}
	else {a=v3;b=v2;}
	for (i=0;i<segTab[a].size();i++) if (segTab[a][i]==b) break;
	//if (i==segTab[a].size()) i/=j;
	lst[2]=cellType(1,segId[a]+i);

	if (v0<v3) {a=v0;b=v3;}
	else {a=v3;b=v0;}
	for (i=0;i<segTab[a].size();i++) if (segTab[a][i]==b) break;
	//if (i==segTab[a].size()) i/=j;
	lst[3]=cellType(1,segId[a]+i);
      }

    }

    return;
  }
    
  void getVertice(cellType cell, std::vector<cellType> &lst) const
  {
    // set list the the IDs of the vertice of cell(type,id)
    typeT type=cell.type();
    if (type >2) return;

    if (type==0) lst.assign(1,cell);
    else if (type==1) getFaces(cell,lst);
    else {
      idT id=cell.id();
      lst.resize(3);
      lst[0]=faceTab[id].a;
      lst[1]=faceTab[id].b;
      lst[2]=faceTab[id].c;
    }
  }

  float *getPosition(cellType cell,float *pos) const
  {
    // set pos to the position of cell(type,id)
    // no need to take care of memory allocation for pos
    typeT type=cell.type();
    idT id=cell.id();
    double a,b;
	  
    if (type==0)
      {
	if (nested) pix2ang_nest(nside,id,&a,&b);
	else pix2ang_ring(nside,id,&a,&b);
	pos[0]=a;pos[1]=b;
      }
    else if (type==1)
      {
	double c,d;
	std::vector<cellType> lst;
	getVertice(cell,lst);
	if (nested) {
	  pix2ang_nest(nside,lst[0].id(),&a,&b);
	  pix2ang_nest(nside,lst[1].id(),&c,&d);
	}
	else {
	  pix2ang_ring(nside,lst[0].id(),&a,&b);
	  pix2ang_ring(nside,lst[1].id(),&c,&d);
	}
	// if segment spans over 2pi limit 
	if (fabs(d-b)>M_PI)
	  {
	    if (d>b) {
	      if (b>2*M_PI-d) d-=2*M_PI;
	      else b+=2*M_PI;
	    } else {
	      if (d>2*M_PI-b) b-=2*M_PI;
	      else d+=2*M_PI;
	    }
	  }

	pos[0]=0.5*(a+c);
	pos[1]=0.5*(b+d);
      }
    else if (type==2)
      {
	int i;
	std::vector<cellType> lst;
	getVertice(cell,lst);	      
	double c,d;
	double delta=0,bi=0;
	c=0;d=0;
	for (i=0;i<lst.size();i++)
	  {
	    if (nested) pix2ang_nest(nside,lst[i].id(),&a,&b);
	    else pix2ang_ring(nside,lst[i].id(),&a,&b);
		  
	    if (i) {
	      if(fabs(b-bi)>M_PI) b+=delta;
	    }
	    else {
	      bi=b;
	      if (b<M_PI) delta=-2*M_PI;
	      else delta=2*M_PI;
	    }
	      
	    c+=a;d+=b;
	  }
	pos[0]=c/lst.size();
	pos[1]=d/lst.size();
	if (pos[1]<0) pos[1]+=2*M_PI;
	if (pos[1]>2*M_PI) pos[1]-=2*M_PI;
      }
  }

  uint getNodeGroup(cellType cell) const
  {
    return 0;
    // returns the group of cell(type,id)
    // 0 means no group
    // used to simplify by groups only
  }

  bool groupsAreDefined() const
  {
    return false;
    // returns true if groups are defined for cells, false otherwise.
    // if "false", getNodeGroup(cellType cell) may return anything.
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
    bool result=is_healpix_map(filename);
    return result;
  }

  static NetworkDataInterface<cellType> *Load(const char *filename) {    
    healpixMap map;

    printf ("Reading HEALPIX map from file %s ...",filename);fflush(0);
    map.hmap = read_healpix_map_d(filename, &map.nside, map.coordsys, map.ordering);

    bool nest;
    if (strstr(map.ordering,"NESTED") != NULL) nest=true;
    else if (strstr(map.ordering,"RING") != NULL) nest =false;
    else 
      {
	fprintf (stderr," done. Unknown ordering (=%s). \n",map.ordering);
	fprintf (stderr,"    Trying with 'NESTED' ...");fflush(0);
	nest=true;
      }

    if (nest)
      printf (" done. (NESTED, nside=%ld)\n",map.nside);
    else
      printf (" done. (RING, nside=%ld)\n",map.nside);

    long i;
    long nout=0;
    std::vector<char> out;
    long np = nside2npix(map.nside);
    for (i=0;i<np;i++)
      {
	if (map.hmap[i]!=map.hmap[i])
	  {
	    if (!nout) out.assign(np,0);
	    out[i]=1;nout++;
	  }
      }

    if (nout) printf ("Found %ld NaN values, a mask will be generated.\n",nout);fflush(0);
    NetworkDataInterface<cellType> *net = new Healpix_network(&map,true,true);
    if (nout) net->setMask(out);

    return net;     
  }

  int getNDims(bool network=false) const 
  {
    return 2;
    // returns the number of dimensions
  }
  void getBoundingBox(std::vector<double> &x0, std::vector<double> &delta) const 
  {
    x0.assign(2,0.);
    x0[0]=-M_PI/2;

    delta.assign(2,3*M_PI/2);
    delta[1]=2*M_PI;
    // set x0 and delta to the origin and size of the cubic bounding box
  }
  
  bool sendMessage(messageT m)
  {
    
    return false;
  }

public:

  Healpix_network(healpixMap *map, bool includeDiags=true, bool isOwner=false,bool getValueByMax = true) : ParentClass(isOwner,getValueByMax)
  {
     if (strstr(map->ordering,"NESTED") != NULL) nested=true;
    else if (strstr(map->ordering,"RING") != NULL) nested=false;
    else {
      fprintf(stderr,"ERROR in Healpix_network(): map is neither RING nor NESTED (=%s).\n",map->ordering);
      exit(0);
    }
    diags=(int)includeDiags;
    nside=map->nside;
    hmap=map->hmap;
    npix = nside2npix(nside);
      
    buildFacesTables();
  }
    
  ~Healpix_network()
  {
    if (ParentClass::netIsOwned()) freeData();
    //delete map;
  }

};


#endif
