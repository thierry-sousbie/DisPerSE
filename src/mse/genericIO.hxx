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
#ifndef __GENERIC_IO_HXX__
#define __GENERIC_IO_HXX__

#include <string>
#include <vector>
#include "mytypes.h"
#include "global.h"

template <class T>
class genericIO_interface {
public:
  typedef T IOType;
  genericIO_interface()
  {}
  virtual ~genericIO_interface()
  {}
  
  virtual std::string getTypeStr()=0;
  virtual std::string getExtensionStr()=0;
  virtual bool canLoad(std::string fname=std::string(""))=0;
  virtual bool canSave()=0;
  virtual IOType *load(std::string fname)=0;
  virtual int save(IOType *,std::string fname)=0;
};

template <class interfaceT, class genT> 
class genericIO 
{
public:
  typedef typename interfaceT::IOType IOType;
  
private:
  typedef std::vector<interfaceT*> allT;

  static allT generateInterface()
  {
    allT all;
    genT::generate(std::back_inserter(all));   
    return all;
  }

  static void destroyInterface(allT &all)
  {
    for (unsigned long i=0;i<all.size();i++)
      delete all[i];
  }
  
  //std::vector<interfaceT *> all; 
  
public:
  /*
  genericIO()
  {
    all.clear();
    genT::generate(std::back_inserter(all));   
  }    
    
  ~genericIO()
  {
    long i;
    for (i=0;i<all.size();i++)
      delete all[i];
  }    
  */
  static std::vector<std::string> getTypeList(bool read, bool write)
  {
    allT all=generateInterface();
    std::vector<std::string> result;

    for (long i=0;i<(long)all.size();i++)
      {
	if ((read && all[i]->canLoad())||(write && all[i]->canSave())) 
	  result.push_back(all[i]->getTypeStr()); 
      }
	
    destroyInterface(all);
    return result;
  }

  static std::string getExtension(std::string type=std::string(""))
  {
    allT all=generateInterface();
    std::string result=std::string("");

    if (type==std::string(""))
      {
	for (long i=0;i<(long)all.size();i++)
	  if (all[i]->canSave()) 
	    {
	      result=all[i]->getExtensionStr();
	      break;
	    }
      }
    else
      {
	for (long i=0;i<(long)all.size();i++)
	  if (all[i]->getTypeStr() == type)
	    {
	      //if (all[i]->canSave())
	      result = all[i]->getExtensionStr();
	      break;
	    }
      }
    
    destroyInterface(all);
    return std::string(".")+result;
  }

  static bool canLoad(std::string fname)
  {
    allT all=generateInterface(); // this was not made static for thread safety ...
    bool result=false;
    
    for (unsigned long i=0;i<all.size();i++)
      if (all[i]->canLoad(fname)) 
	{
	  result=true;
	  break;
	}
    
    destroyInterface(all);
    return result;
  }

  static bool canSave(std::string type)
  {
    allT all=generateInterface(); // this was not made static for thread safety ...
    bool result=false;
    
    for (long i=0;i<(long)all.size();i++)
      if (all[i]->getTypeStr() == type)
	{
	  if (all[i]->canSave()) 
	    {
	      result=true;
	      break;
	    }
	}
    
    destroyInterface(all);
    return result;
  }

  static IOType *load(std::string fname)
  {
    allT all=generateInterface();
    IOType *result=NULL;
    int oldverbose=verbose;
    //verbose=0;

    for (unsigned long i=0;i<all.size();i++)
      if (all[i]->canLoad(fname))
	{
	  if (verbose>=1) printf("Will read %s from file '%s'.\n",genT::type().c_str(),
				 fname.c_str());
	  fflush(0);
	  result=all[i]->load(fname);
	  //printf("done.\n");
	  break;
	}

    if (result==NULL) fprintf(stderr,"ERROR: file '%s' cannot be read.\n",fname.c_str());
      

    verbose=oldverbose;
    destroyInterface(all);
    return result;
  }

  static int save(IOType *f, std::string fname, std::string type=std::string(""))
  {
    allT all=generateInterface();
    int result=-1;
    int oldverbose=verbose;
    //verbose=0;

    //printf("size : %ld\n",all.size());
    if (type==std::string(""))
      {
	for (long i=0;i<(long)all.size();i++)
	  if (all[i]->canSave()) 
	    {
	      printf("Will write %s to file '%s'.\n",genT::type().c_str(),fname.c_str());
	      fflush(0);
	      result=all[i]->save(f,fname);
	      //printf("done.\n");
	      break;
	    }
      }
    else
      {
	for (long i=0;i<(long)all.size();i++)
	  {
	    //printf("type=%s ... == %s? -> %d\n",type.c_str(),all[i]->getTypeStr().c_str(),
	    //(int)(type==all[i]->getTypeStr()));
	    if (all[i]->getTypeStr() == type)
	      {
		if (all[i]->canSave())
		  {
		    printf("Will write %s to file '%s'.\n",genT::type().c_str(),fname.c_str());
		    // printf("Saving network to file %s ...",fname.c_str());
		    fflush(0);
		    result = all[i]->save(f,fname);
		    //printf("done.\n");
		  }
		break;
	      }
	  }
      }

    verbose=oldverbose;
    destroyInterface(all);
    return result;
  }
  
};

#endif
