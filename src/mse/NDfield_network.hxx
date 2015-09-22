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
#ifndef __NDFIELD_MSC_HEADER__
#define __NDFIELD_MSC_HEADER__

// Use this as an interface between your own datatype and MSComplex class.
// Change NDfield to whetever your datatype is and implement necessary functions.
#include "network_interface.hxx"

#include "MSComplex.hxx"
#include "mytypes.h"

#include "NDfield.h"
#include "hypercube.h"

#include <stdio.h>
#include <string.h>

//template <typename TypeT = char, typename IdT = uint>
template <class cellType = cellIdentity<> >
class NDfield_network : public NetworkDataInterface<cellType> {

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
  //class NDfield_MSComplex ; 

    void Index2CoordsND(uint index, uint *coord) const
	{
	    idT val;idT old_val;
	    int i;
	    
	    val=1;
	    for (i=0;i<net->ndims;i++)
	    {     
		old_val=val;
		val*=net->dims[i];
		coord[i]=(index%val)/old_val;
	    }
	}
   
private:
    hypercube_type *cube;
    nfacelist_type *cubeFace;

    NDfield *net; // your datatype
    double *val;
    double ndims_inv;
    uint dx;
    int periodic;
  int periodicity;

  std::vector< std::vector<bool> > mask_out;
  std::vector< std::vector<bool> > mask_boundary;

  //typeT type;
  //idT id;

public:
    
  int getPeriodicity() const {return periodicity;} 

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
    Free_NDfield(&net);
  }

  void setMask(std::vector<char> &mask, bool nullIsMasked)
  {
    long i,j,k;
    typeT type;
    long nf;
    long nb=0;
    cellType curCell;

    if (mask.size())
      {
	if (mask.size() != net->nval)
	  {
	    fprintf(stderr,"ERROR in setMask: mask and data do not have the same size.\n");
	    exit(-1);
	  }
      }
    else if (periodicity==~(int)0) return;

    printf ("Setting mask ... ");fflush(0);

    if (!mask_out.size())
      {
	mask_out.resize(net->ndims+1);
	mask_out[0].assign(this->getNFaces(0),false);
      }

    mask_boundary.resize(net->ndims+1);
    for (i=0;i<=net->ndims;i++) mask_boundary[i].assign(this->getNFaces(i),false);
    for (i=1;i<=net->ndims;i++) mask_out[i].assign(this->getNFaces(i),false);
    
    if (mask.size())
      {
	if (nullIsMasked)
	  for (i=0;i<mask.size();i++) mask_out[0][i] = mask_out[0][i] | (!(bool)mask[i]);
	else
	  for (i=0;i<mask.size();i++) mask_out[0][i] = mask_out[0][i] | (bool)mask[i];
      }

    for (type=1;type<=net->ndims;type++) {
      nf = this->getNFaces(type);
      for (i=0;i<nf;i++)
	{
	  static std::vector<cellType> vert;
	  if (mask_out[type][i]) continue;
	  getVertice(cellType(type,i),vert);
	  for (j=0;j<vert.size();j++) if (isOut(vert[j])) break;
	  if (j!=vert.size()) mask_out[type][i]=true; 
	  else if (periodicity!=~(int)0)
	    {
	      uint coords[net->ndims];
	      std::vector<uint> mn(net->ndims,1000000000);
	      std::vector<uint> mx(net->ndims,0);
	      for (j=0;j<vert.size();j++)
		{
		  Index2CoordsND(vert[j].id(),coords);
		  for (k=0;k<net->ndims;k++)
		    {		      
		      if (coords[k]<mn[k]) mn[k]=coords[k];
		      if (coords[k]>mx[k]) mx[k]=coords[k];
		    }
		}
	      for (k=0;k<net->ndims;k++)
		{
		  if (periodicity&(1<<k)) continue;
		  if (mx[k]-mn[k]>1) break;
		}
	      if (k!=net->ndims) mask_out[type][i]=true; 
	    }
	}
    }

    type=net->ndims-1;
    nf=this->getNFaces(type);
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
    
    for (type=net->ndims-2;type>0;type--)
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

    printf ("done. (%d %d-faces on boundary)\n",(int)nb,net->ndims-1);

    /*
    std::vector<char> tmp(mask_out[0].size());
    for (i=0;i<mask_out[0].size();i++) tmp[i]=(char)mask_out[0][i];
    NDfield *fd=Create_NDfield(net->dims,net->ndims,net->fdims_index,ND_CHAR,net->x0, net->delta, (void*) &(tmp[0]),net->comment);
    Save_NDfield(fd,"out0.ND");
    for (i=0;i<mask_out[0].size();i++) tmp[i]=(char)mask_boundary[0][i];
    Save_NDfield(fd,"boundary0.ND");
    */
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

	  // return the value of the cell, sum being the sum of the values for each vertex
	  // sum is used to compare two different cells with identical type and value
	  if (type==0) {sum=val[id];return sum;}
	  if (getByMax==-1) getByMax=ParentClass::getByMax();
	  
	  nface_type *face = HC_getFace(cube,type,id,net->dims,periodic);
	  uint cube_index = HC_getCubeIndex(cube,type,id,periodic);
	  uint cube_coords[net->ndims];
	  Index2CoordsND(cube_index,cube_coords);
	  
	  nface_type *vertex = &(cubeFace[net->ndims].NFace[face->vertex[0]]);
	  uint vid = HC_getIndex(cube,vertex,cube_coords,net->dims,periodic);
	  int i;
	  
	  sum=val[vid];
	  double valmax=sum;
	  if (getByMax)
	    for (i=1;i<face->nvertice;i++)
	      {
		vertex = &(cubeFace[net->ndims].NFace[face->vertex[i]]);
		vid = HC_getIndex(cube,vertex,cube_coords,net->dims,periodic);
		if (val[vid]>valmax) valmax=val[vid];
		sum+=val[vid];
	      }
	  else
	    for (i=1;i<face->nvertice;i++)
	      {
		vertex = &(cubeFace[net->ndims].NFace[face->vertex[i]]);
		vid = HC_getIndex(cube,vertex,cube_coords,net->dims,periodic);
		if (val[vid]<valmax) valmax=val[vid];
		sum+=val[vid];
	      }
	  
	  return valmax;	    
	}
  
  double getValueAll( cellType cell, std::vector<double> &values, int getByMax)  const
	{
	  typeT type=cell.type();
	  idT id=cell.id();

	  // return the value of the cell, sum being the sum of the values for each vertex
	  // sum is used to compare two different cells with identical type and value
	  if (type==0) {values.resize(1);values[0]=val[id];return val[id];}
	  if (getByMax==-1) getByMax=ParentClass::getByMax();

	  nface_type *face = HC_getFace(cube,type,id,net->dims,periodic);
	  uint cube_index = HC_getCubeIndex(cube,type,id,periodic);
	  uint cube_coords[net->ndims];
	  Index2CoordsND(cube_index,cube_coords);
	  
	  nface_type *vertex = &(cubeFace[net->ndims].NFace[face->vertex[0]]);
	  uint vid = HC_getIndex(cube,vertex,cube_coords,net->dims,periodic);
	  int i;

	  long max=vid;

	  values.resize(face->nvertice);
	  values[0]=val[vid];

	    if (getByMax)
	      for (i=1;i<face->nvertice;i++)
		{
		  vertex = &(cubeFace[net->ndims].NFace[face->vertex[i]]);
		  vid = HC_getIndex(cube,vertex,cube_coords,net->dims,periodic);
		  if (val[vid]>val[max]) max=vid;
		  values[i]=val[vid];
		}
	    else
	      for (i=1;i<face->nvertice;i++)
		{
		  vertex = &(cubeFace[net->ndims].NFace[face->vertex[i]]);
		  vid = HC_getIndex(cube,vertex,cube_coords,net->dims,periodic);
		  if (val[vid]<val[max]) max=vid;
		  values[i]=val[vid];
		}

	    std::swap(values[0],values[max]);

	    return val[max];	    
	}

    double getValue( cellType cell, int getByMax) const
	{
	  typeT type=cell.type();
	  idT id=cell.id();

	  if (type==0) return val[id];
	  if (getByMax==-1) getByMax=ParentClass::getByMax();

	  nface_type *face = HC_getFace(cube,type,id,net->dims,periodic);
	  uint cube_index = HC_getCubeIndex(cube,type,id,periodic);
	  uint cube_coords[net->ndims];
	  Index2CoordsND(cube_index,cube_coords);
	  
	  nface_type *vertex = &(cubeFace[net->ndims].NFace[face->vertex[0]]);
	  uint vid = HC_getIndex(cube,vertex,cube_coords,net->dims,periodic);
	  int i;
	  //printf("\n%d %e ",face->nvertice,val[vid]);
	  double valmax=val[vid];
	  if (getByMax)
	    for (i=1;i<face->nvertice;i++)
	      {
		vertex = &(cubeFace[net->ndims].NFace[face->vertex[i]]);
		vid = HC_getIndex(cube,vertex,cube_coords,net->dims,periodic);
		if (val[vid]>valmax) valmax=val[vid];
		//printf("%e ",val[vid]);
	      }
	  else
	    for (i=1;i<face->nvertice;i++)
	      {
		vertex = &(cubeFace[net->ndims].NFace[face->vertex[i]]);
		vid = HC_getIndex(cube,vertex,cube_coords,net->dims,periodic);
		if (val[vid]<valmax) valmax=val[vid];
		//printf("%e ",val[vid]);
	      }
	  //printf(" -> %e \n",valmax);
	  return valmax;
	}
  
    unsigned long getNFaces(int type) const
	{
	  
	    // returns the total number of faces of a given type 
	  /*
	    if (periodic==(~0))
		return cube->NFaceList[net->ndims-type].nFacesPerCube*net->nval;
	    else
	  */
		return cube->NFaceList[net->ndims-type].nFacesPerCube*net->nval;
	}

    void getCofaces(cellType cell, std::vector<cellType> &result) const
	{
	  typeT type=cell.type();
	  idT id=cell.id();

	    // set result to the IDs of the cofaces of cell(type,id)
	    // A=cell(type,id) is a coface of B=cell(type-1,coID) if B is a face of A.
	    nface_type *face = HC_getFace(cube,type,id,net->dims,periodic);
	    uint cube_index = HC_getCubeIndex(cube,type,id,periodic);
	    uint cube_coords[net->ndims];
	    Index2CoordsND(cube_index,cube_coords);
	    result.resize(face->nmothers*2);
	    int i;
	    for (i=0;i<face->nmothers;i++)
	    {
	      result[i]=cellType(type+1,HC_getIndex(cube,face->Mothers[i],cube_coords,net->dims,periodic));
	    }

	    nface_type *new_face;
	    uint new_cube_index; 
	    HC_getSymetricConfig(cube,face,cube_coords,&new_face,&new_cube_index,net->dims,periodic);
	    
	    Index2CoordsND(new_cube_index,cube_coords);
	    
	    for (i=0;i<new_face->nmothers;i++)
	    {
	      result[face->nmothers+i]=cellType(type+1,HC_getIndex(cube,new_face->Mothers[i],cube_coords,net->dims,periodic));
	    }

	}
    
    void getFaces(cellType cell, std::vector<cellType> &result) const
	{
	  typeT type=cell.type();
	  idT id=cell.id();

	    // set result to the IDs of the faces of cell(type,id)
	    nface_type *face = HC_getFace(cube,type,id,net->dims,periodic);
	    uint cube_index = HC_getCubeIndex(cube,type,id,periodic);
	    uint cube_coords[net->ndims];
	    Index2CoordsND(cube_index,cube_coords);
	    int i;
	    result.resize(face->ndaughters);
	    for (i=0;i<face->ndaughters;i++)
	    {
		//printf ("D:%d\n",face->Daughters[i]->tag);
	      result[i]=cellType(type-1,HC_getIndex(cube,face->Daughters[i],cube_coords,net->dims,periodic));
	    }
/*
	    printf ("cell(%d,%d)=[ ",type,id);
	    for (i=0;i<net->ndims;i++) printf ("%d ",cube_coords[i]);
	    printf("] has %d faces: ",face->ndaughters);
	    for (i=0;i<result.size();i++) printf ("%u ",result[i]);
	    printf (".\n");
*/
	}
    
    void getVertice(cellType cell, std::vector<cellType> &result) const
	{
	  typeT type=cell.type();
	  idT id=cell.id();

	    // set result the the IDs of the vertice of cell(type,id)
	    if (type==0) {result.assign(1,cell);return;}
	    
	    nface_type *face = HC_getFace(cube,type,id,net->dims,periodic);
	    uint cube_index = HC_getCubeIndex(cube,type,id,periodic);
	    uint cube_coords[net->ndims];
	    Index2CoordsND(cube_index,cube_coords);
	    nface_type *vertex;
	    uint vid;
	    int i;
	    result.resize(face->nvertice);
	    for (i=0;i<face->nvertice;i++)
	    {
		vertex = &(cubeFace[net->ndims].NFace[face->vertex[i]]);
		vid = HC_getIndex(cube,vertex,cube_coords,net->dims,periodic);
		result[i]=cellType(0,vid);
	    }
	}

    float *getPosition(cellType cell,float *pos) const
	{
	  typeT type=cell.type();
	  idT id=cell.id();

	    // set pos to the position of cell(type,id)
	    // no need to take care of memory allocation for pos
	    std::vector<cellType> v;
	    getVertice(cell,v);

	    int i,j;
	    uint coordsI[net->ndims];
	    uint c0[net->ndims];
	    double coords[net->ndims];

	    Index2CoordsND(v[0].id(),coordsI);
	    for (j=0;j<net->ndims;j++) 
	    {
		c0[j] = coordsI[j];
		coords[j]=coordsI[j];
	    }	    
	    for (i=1;i<v.size();i++)
	    {
	      Index2CoordsND(v[i].id(),coordsI);
		for (j=0;j<net->ndims;j++) 
		{
		    if (abs(c0[j]-coordsI[j])<(net->dims[j]>>1))
			coords[j]+=coordsI[j];
		    else
		    {
			if (c0[j]<coordsI[j])
			    coords[j]+=(coordsI[j]-net->dims[j]);
			else
			    coords[j]+=(coordsI[j]+net->dims[j]);
		    }
		}	    
	    }
	    double ninv = 1./v.size();
	    for (j=0;j<net->ndims;j++) 
		pos[j] = net->x0[j]+coords[j]*ninv*net->delta[j]/net->dims[j];
	    
	}

    uint getNodeGroup(cellType cell) const
	{
	  typeT type=cell.type();
	  idT id=cell.id();

	    // returns the group of cell(type,id)
	    // 0 means no group
	    // used to simplify by groups only
	    return 0;
	}

    bool groupsAreDefined() const
	{
	    // returns true if groups are defined for cells, false otherwise.
	    // if "false", getNodeGroup(cellType cell) may return anything.
	    return false;
	}

  bool dumpSubComplex(const char *fname,subCellsType &subList,cellsInfoType &subListInfo, newCellsType &newList,cellsInfoType &newListInfo,int smooth, bool allowDuplicates) const
  {
    return false;
  }


public:
  static bool checkFileType(const char *filename,bool simplicial=false) { 
    if (simplicial) return false;
    if (IsNDfield(filename)) return true;
    if (IsNDfieldFromFITS(filename)) return true;
    /*
#ifdef HAVE_CFITS_IO
    else
      {
	int status = 0;
	fitsfile *fptr;
	if (!fits_open_file(&fptr, filename, READONLY, &status))
	  {
	    fits_close_file(fptr, &status);
	    return true;
	  }
	return false;
      }
#endif
    */
    return false;
  }

  static NetworkDataInterface<cellType> *Load(const char *filename) {
   
    double *data;
    float *dataf;
    NDfield *field=NULL;
    long i;

    if (IsNDfield(filename))
      field=Load_NDfield(filename);
    else if (IsNDfieldFromFITS(filename))
      field=Load_NDfieldFromFITS(filename);
    

    if (field==NULL)
      {
	fprintf(stderr,"ERROR: cannot load '%s', wrong format.\n",filename);
	fprintf(stderr,"  should have type 'NDfield' or a FITS image file.\n");
	exit(-1);
      }
    
    if ((field->datatype!=ND_FLOAT)&&(field->datatype!=ND_DOUBLE))
      {
	fprintf(stderr,"ERROR: Array data type is neither DOUBLE nor FLOAT\n");
	exit(0);
      }

    if (field->datatype==ND_FLOAT) 
      {
	dataf=(float *)field->val;
	data=NULL;
      }
    else
      {
	data=(double *)field->val;
	dataf=NULL;
      }

    long nout=0;
    std::vector<char> out;
    if (data==NULL)
      {
	for (i=0;i<field->nval;i++)
	  if (dataf[i]!=dataf[i])
	    {
	      if (!nout) out.assign(field->nval,0);
	      out[i]=1;nout++;
	    }
      }
    else
      {
	for (i=0;i<field->nval;i++)
	  if (data[i]!=data[i])
	    {
	      if (!nout) out.assign(field->nval,0);
	      out[i]=1;nout++;
	    }
      }

    if (nout) printf ("Found %ld NaN values, a mask will be generated.\n",nout);fflush(0);
    NetworkDataInterface<cellType> *net = new NDfield_network(field,true);
    if (nout) net->setMask(out);
     
    return net;
	//return NULL;      
  }

    int getNDims(bool network=false) const 
	{
	    // returns the number of dimensions
	    return net->ndims;
	}
    void getBoundingBox(std::vector<double> &x0, std::vector<double> &delta) const 
	{
	    // set x0 and delta to the origin and size of the cubic bounding box
	    x0.assign(net->x0,&(net->x0[net->ndims]));
	    delta.assign(net->delta,&(net->delta[net->ndims]));
	}
  
  bool sendMessage(messageT m)
  {
    
    return false;
  }

public:

  NDfield_network(NDfield *netp,bool isOwner=false,bool getValueByMax = true) : ParentClass(isOwner,getValueByMax)
  {		    
	    // set isOwner to true if netp can be freed/modified
	    // nthreads is the number of threads
	    
	    // All the setup goes here ...
	    net=netp;
	    Convert_NDfield(net,ND_DOUBLE);
	    val=(double *)net->val;
	    ndims_inv=1./net->ndims;
	    periodic = ~(int)0;
	    periodicity=~(int)0;
	    printf ("Generating %dD cube topology ... ",net->ndims);fflush(0);
	    cube = GenerateHypercube(net->ndims);
	    cubeFace = cube->NFaceList;
	    printf("done.\n");
	    std::vector<char> dummy;
	    setMask(dummy,true);
	    // ... and don't forget to call init() once set ... or not :)
	    //init();

	}
    
    ~NDfield_network()
    {
	FreeHypercube(cube);
	if (ParentClass::netIsOwned()) freeData();
    }

};

#endif
