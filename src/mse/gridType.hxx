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
#ifndef __GRID_TYPE_HEADER__
#define __GRID_TYPE_HEADER__

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>

#include "NDfield.h"
#include "mytypes.h"
#include "global.h"
#include "NDfield_interface.hxx"

template<typename T> struct type_name
{
  static const char* name() { return std::string("ERROR : DATATYPE NOT IMPLEMENTED").c_str(); }
};

#define DECL_TYPE_NAME(x) template<> struct type_name<x> { static const std::string name() {return std::string(#x);} }
DECL_TYPE_NAME(char);
DECL_TYPE_NAME(unsigned char);
DECL_TYPE_NAME(short);
DECL_TYPE_NAME(unsigned short);
DECL_TYPE_NAME(int);
DECL_TYPE_NAME(unsigned int);
DECL_TYPE_NAME(long);
DECL_TYPE_NAME(unsigned long);
DECL_TYPE_NAME(bool);
DECL_TYPE_NAME(float);
DECL_TYPE_NAME(double);

template <typename T=double>
class gridType 
{
public:
  typedef T dataT;
 
private:
   typedef ndfield::IO IO;
 
public:  
  dataT *val; 
  std::string comment;
  int ndims;
  long nval;
  int dataType;
  std::vector<int> dims;
  std::vector<double> x0;
  std::vector<double> delta; 
  

  static bool isLoadable(const char *filename)
  {
    return IO::canLoad(filename); 
  }
  
  bool isLoadable(const std::string fname)
  {
    return isLoadable(fname.c_str());
  }

  bool load(const char *filename)
  {
    if (val!=NULL) free(val);
    NDfield *field=IO::load(std::string(filename));
    
    if (field!=NULL)
      {
	long i;
	if (field->fdims_index!=0)
	  {
	    fprintf(stderr,"WARNING (will be ignored) :\n   File %s seems to contain coordinates, not value.\n   Maybe fdims_index is not correctly set (current value : %d, should be 0).\n",filename,field->fdims_index);
	  }
	
	if (type_name<dataT>::name()==type_name<unsigned char>::name())
	  {Convert_NDfield(field,ND_UCHAR);dataType=ND_UCHAR;}
	else if (type_name<dataT>::name()==type_name<char>::name())
	  {Convert_NDfield(field,ND_CHAR);dataType=ND_CHAR;}
	else if (type_name<dataT>::name()==type_name<unsigned short>::name())
	  {Convert_NDfield(field,ND_USHORT);dataType=ND_USHORT;}
	else if (type_name<dataT>::name()==type_name<short>::name())
	  {Convert_NDfield(field,ND_SHORT);dataType=ND_SHORT;}
	else if (type_name<dataT>::name()==type_name<unsigned int>::name())
	  {Convert_NDfield(field,ND_UINT);dataType=ND_UINT;}
	else if (type_name<dataT>::name()==type_name<int>::name())
	  {Convert_NDfield(field,ND_INT);dataType=ND_INT;}
	else if (type_name<dataT>::name()==type_name<unsigned long>::name())
	  {Convert_NDfield(field,ND_ULONG);dataType=ND_ULONG;}
	else if (type_name<dataT>::name()==type_name<long>::name())
	  {Convert_NDfield(field,ND_LONG);dataType=ND_LONG;}
	else if (type_name<dataT>::name()==type_name<float>::name())
	  {Convert_NDfield(field,ND_FLOAT);dataType=ND_FLOAT;}
	else if (type_name<dataT>::name()==type_name<double>::name())
	  {Convert_NDfield(field,ND_DOUBLE);dataType=ND_DOUBLE;}
	else 
	  {
	    printf ("ERROR: class gridType<dataT> not implemented for type dataT.\n");
	    exit(0);
	  }
	
	ndims=field->ndims;
	x0.assign(field->x0,field->x0+field->ndims);
	delta.assign(field->delta,field->delta+field->ndims);
	dims.assign(field->dims,field->dims+field->ndims);
	comment=std::string(field->comment);
	val = (dataT *) field->val;
	nval=1;
	for (i=0;i<ndims;i++) nval*=dims[i];
	//save("test.ND");
	return true;
      }
 
    fprintf(stderr,"ERROR: cannot load '%s', wrong file type.\n",filename);
    fprintf(stderr,"       Expected type : grid.\n");
    exit(-1);
    
    return false;
  }

  bool load(const std::string fname)
  {
    return load(fname.c_str());
  }

  int save(const std::string fname,std::string type=std::string(""))
  {
    if (val==NULL) return -1;
    NDfield *f=Create_NDfield(&dims[0],ndims,0,dataType,&x0[0],&delta[0],(void*)val,"");
    return IO::save(f,fname,type);
  }

  static int save(NDfield *f,const std::string fname,std::string type=std::string(""))
  {
    return IO::save(f,fname,type);
  }
  
  gridType(const std::string fname):val(NULL)
  {
    if (fname!=std::string("")) load(fname);
  }

  gridType(const char *fname):val(NULL)
  {
    if (fname!=NULL) load(fname);
    
  }

  gridType():val(NULL)
  {

  }

  virtual ~gridType()
  {
    if (val!=NULL) free(val);
  }

  dataT *steal()
  {
    dataT* v2=val;
    val=NULL;
    return v2;
  }
  
  std::vector<dataT> toVector()
  {
    std::vector<dataT> result;
    if (val!=NULL) result.assign(val,val+nval);
    return result;
  }

  void toVector(std::vector<dataT> &out)
  {
    if (val==NULL) out.clear();
    else out.assign(val,val+nval);
  }
  
};

#endif
