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
#ifndef __NDNET_INTERFACE_HXX__
#define __NDNET_INTERFACE_HXX__

#include <string>
#include <limits>

#include "genericIO.hxx"
#include "NDnet_VTK_IO.h"
#include "NDnet_PLY_IO.h"
#include "NDnetwork.h"
#include "gadget_io.h"
#include "global.h"
#include "mystring.h"

#include "NDskel_breakdown.h"
#include "NDskel_interface.hxx"

namespace ndnet {
  
  typedef genericIO_interface<NDnetwork> interface;

  class fromNDnet : public interface {
    std::string getTypeStr() {return std::string("NDnet");}
    std::string getExtensionStr() {return std::string("NDnet");}
    bool canLoad(std::string fname) 
    {
      if (fname==std::string("")) return true;
      if (IsNDnetwork(fname.c_str())==1) return true;
      return false;
    }
    bool canSave() {return true;}
    
    NDnetwork *load(std::string fname)
    {
      return Load_NDnetwork(fname.c_str());
    }
    
    int save(NDnetwork *f, std::string fname)
    {
      return Save_NDnetwork(f,fname.c_str());
    }
  };

  class fromNDnet_ascii : public interface {
    std::string getTypeStr() {return std::string("NDnet_ascii");}
    std::string getExtensionStr() {return std::string("a.NDnet");}
    bool canLoad(std::string fname) 
    {
      if (fname==std::string("")) return true;
      if (IsNDnetwork(fname.c_str())==2) return true;
      return false;
    }
    bool canSave() {return true;}
    
    NDnetwork *load(std::string fname)
    {
      return Load_NDnetwork(fname.c_str());
    }
    
    int save(NDnetwork *f, std::string fname)
    {   
      return Save_NDnetwork_ASCII(f,fname.c_str());
    }
  };

  class fromGadget : public interface {
    std::string getTypeStr() {return std::string("gadget");}
    std::string getExtensionStr() {return std::string("gad");}
    bool canLoad(std::string fname) 
    {
      if (fname==std::string("")) return true;
      if (IsGadgetFile(fname.c_str())==1) return true;
      return false;
    }
    bool canSave() {return false;}
    
    NDnetwork *load(std::string fname)
    {
      snapshot_data *g=(snapshot_data *)calloc(1,sizeof(snapshot_data));
      ReadGadget(fname.c_str(),g,FLAG_ALL);
      unsigned long npart=0;
      for (int i=0;i<6;i++) npart+=g->header.npart[i];
      NDnetwork *net=CreateNetwork(3,npart,false);
      std::copy(g->Pos,g->Pos + npart*3,net->v_coord);    
      net->ndata=4;
      net->data=(NDnetwork_Data*)calloc(net->ndata,sizeof(NDnetwork_Data));
      strcpy(net->data[0].name,"Vx");
      strcpy(net->data[1].name,"Vy");
      strcpy(net->data[2].name,"Vz");
      strcpy(net->data[3].name,"Id");
      for (int i=0;i<4;++i)
	net->data[i].data=(double*)malloc(npart*sizeof(double));

      for (unsigned long i=0;i<npart;i++)
	net->data[0].data[i]=g->Vel[i*3];
      for (unsigned long i=0;i<npart;i++)
	net->data[1].data[i]=g->Vel[i*3+1];
      for (unsigned long i=0;i<npart;i++)
	net->data[2].data[i]=g->Vel[i*3+2];
      for (unsigned long i=0;i<npart;i++)
	net->data[3].data[i]=g->Id[i];      

      freegadgetstruct(g);
      return net;
    }
    
    int save(NDnetwork *f, std::string fname)
    {   
      return -1;
    }
  };
  
  class fromPLY : public interface {
    std::string getTypeStr() {return std::string("ply");}
    std::string getExtensionStr() {return std::string("ply");}
    bool canLoad(std::string fname) 
    {
      if (fname==std::string("")) return true;
      if (IsPLYNetwork(fname.c_str())) return true;
      return false;
    }
    bool canSave() {return true;}
    NDnetwork *load(std::string fname)
    {
      return Load_NDnetworkFromPLY(fname.c_str());
    }
    
    int save(NDnetwork *f, std::string fname)
    {
      return Save_NDnetworkToPLY(f,fname.c_str(),NDNET_PLY_BIN);
    }
  };

  class fromPLY_ascii : public interface {
    std::string getTypeStr() {return std::string("ply_ascii");}
    std::string getExtensionStr() {return std::string("a.ply");}
    bool canLoad(std::string fname) {return false;}
    bool canSave() {return true;}
    NDnetwork *load(std::string fname)
    {
      return Load_NDnetworkFromPLY(fname.c_str());
    }
    
    int save(NDnetwork *f, std::string fname)
    {
      return Save_NDnetworkToPLY(f,fname.c_str(),NDNET_PLY_ASCII);
    }
  };

  class fromVTK : public interface {
    std::string getTypeStr() {return std::string("vtk");}
    std::string getExtensionStr() {return std::string("vtk");}
    bool canLoad(std::string fname) {return false;}
    bool canSave() {return true;}
    NDnetwork *load(std::string fname)
    {
      return Load_NDnetworkFromVTK(fname.c_str());
    }
    
    int save(NDnetwork *f, std::string fname)
    {
      return Save_NDnetworkToVTK(f,fname.c_str(),NDNET_VTK_BIN|NDNET_VTK_LEGACY);
    }
  };

  class fromVTK_ascii : public interface {
    std::string getTypeStr() {return std::string("vtk_ascii");}
    std::string getExtensionStr() {return std::string("a.vtk");}
    bool canLoad(std::string fname) {return false;}
    bool canSave() {return true;}
    NDnetwork *load(std::string fname)
    {
      return NULL;
    }
    
    int save(NDnetwork *f, std::string fname)
    {
      return Save_NDnetworkToVTK(f,fname.c_str(),NDNET_VTK_ASCII|NDNET_VTK_LEGACY);
    }
  };

  class fromVTK_xml : public interface {
    std::string getTypeStr() {return std::string("vtu");}
    std::string getExtensionStr() {return std::string("vtu");}
    bool canLoad(std::string fname) {return false;}
    bool canSave() {return true;}
    NDnetwork *load(std::string fname)
    {
      return NULL;
    }
    
    int save(NDnetwork *f, std::string fname)
    {
      return Save_NDnetworkToVTK(f,fname.c_str(),NDNET_VTK_BIN|NDNET_VTK_XML);
    }
  };

  class fromVTK_xml_ascii : public interface {
    std::string getTypeStr() {return std::string("vtu_ascii");}
    std::string getExtensionStr() {return std::string("a.vtu");}
    bool canLoad(std::string fname) {return false;}
    bool canSave() {return true;}
    NDnetwork *load(std::string fname)
    {
      return NULL;
    }
    
    int save(NDnetwork *f, std::string fname)
    {
      return Save_NDnetworkToVTK(f,fname.c_str(),NDNET_VTK_ASCII|NDNET_VTK_XML);
    }
  };

  class fromNDskel : public interface {
    std::string getTypeStr() {return std::string("skeleton(any)");}
    std::string getExtensionStr() {return std::string("skl");}
    bool canLoad(std::string fname) {return ndskel::IO::canLoad(fname);}
    bool canSave() {return false;}
    NDnetwork *load(std::string fname)
    {
      NDskel *skl=ndskel::IO::load(fname);
      NDnetwork *net = NDskel2NDnet(skl);
      Free_NDskeleton(&skl);

      return net;
    }
    
    int save(NDnetwork *f, std::string fname)
    {
      return -1;
    }
  };
  
  struct allTypes
  {
    template <class OutputIterator>
    static void generate(OutputIterator out)
    {
      *out = new fromNDnet();
      *out = new fromNDnet_ascii();
      *out = new fromPLY();
      *out = new fromPLY_ascii();
      *out = new fromVTK();
      *out = new fromVTK_ascii();
      *out = new fromVTK_xml();
      *out = new fromVTK_xml_ascii();
      *out = new fromGadget();
      //*out = new fromNDskel();
    }
    static std::string type()
    {
      return std::string("network");
    }
  };

  typedef genericIO<interface, allTypes> IO;
}

#endif
