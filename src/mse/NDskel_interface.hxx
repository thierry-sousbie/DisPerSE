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
#ifndef __NDSKEL_INTERFACE_HXX__
#define __NDSKEL_INTERFACE_HXX__

#include <string>
#include <limits>

#include "NDskeleton.h"
#include "NDnetwork.h"
#include "NDskel_tags.h"
#include "NDskel_VTK_IO.h"
#include "genericIO.hxx"
#include "global.h"
#include "mystring.h"

namespace ndskel {
  
  typedef genericIO_interface<NDskel> interface;

  class fromNDskel : public interface {
    std::string getTypeStr() {return std::string("NDskl");}
    std::string getExtensionStr() {return std::string("NDskl");}
    bool canLoad(std::string fname) 
    {
      if (fname==std::string("")) return true;
      if (IsNDskeleton(fname.c_str())==1) return true;
      return false;
    }
    bool canSave() {return true;}
    
    NDskel *load(std::string fname)
    {
      return Load_NDskel(fname.c_str());
    }
    
    int save(NDskel *f, std::string fname)
    {
      return Save_NDskel(f,fname.c_str());
    }
  };

   class fromNDskel_ascii : public interface {
    std::string getTypeStr() {return std::string("NDskl_ascii");}
    std::string getExtensionStr() {return std::string("a.NDskl");}
    bool canLoad(std::string fname) 
    {
      if (fname==std::string("")) return false;
      if (IsNDskeleton(fname.c_str())==2) return true;
      return false;
    }
    bool canSave() {return true;}
    
    NDskel *load(std::string fname)
    {
      return Load_NDskel(fname.c_str());
    }
    
    int save(NDskel *f, std::string fname)
    {
      return Save_ASCIIskel(f,fname.c_str());
    }
  };

  class fromVTK : public interface {
    std::string getTypeStr() {return std::string("vtk");}
    std::string getExtensionStr() {return std::string("vtk");}
    bool canLoad(std::string fname) {return false;}
    bool canSave() {return true;}
    NDskel *load(std::string fname)
    {
      return Load_NDskeletonFromVTK(fname.c_str());
    }
    
    int save(NDskel *f, std::string fname)
    {
      return Save_NDskeletonToVTK(f,fname.c_str(),NDSKEL_VTK_BIN|NDSKEL_VTK_LEGACY);
    }
  };

  class fromVTK_ascii : public interface {
    std::string getTypeStr() {return std::string("vtk_ascii");}
    std::string getExtensionStr() {return std::string("a.vtk");}
    bool canLoad(std::string fname) {return false;}
    bool canSave() {return true;}
    NDskel *load(std::string fname)
    {
      return NULL;
    }
    
    int save(NDskel *f, std::string fname)
    {
      return Save_NDskeletonToVTK(f,fname.c_str(),NDSKEL_VTK_ASCII|NDSKEL_VTK_LEGACY);
    }
  };

  class fromVTK_xml : public interface {
    std::string getTypeStr() {return std::string("vtp");}
    std::string getExtensionStr() {return std::string("vtp");}
    bool canLoad(std::string fname) {return false;}
    bool canSave() {return true;}
    NDskel *load(std::string fname)
    {
      return NULL;
    }
    
    int save(NDskel *f, std::string fname)
    {
      return Save_NDskeletonToVTK(f,fname.c_str(),NDSKEL_VTK_BIN|NDSKEL_VTK_XML);
    }
  };

  class fromVTK_xml_ascii : public interface {
    std::string getTypeStr() {return std::string("vtp_ascii");}
    std::string getExtensionStr() {return std::string("a.vtp");}
    bool canLoad(std::string fname) {return false;}
    bool canSave() {return true;}
    NDskel *load(std::string fname)
    {
      return NULL;
    }
    
    int save(NDskel *f, std::string fname)
    {
      return Save_NDskeletonToVTK(f,fname.c_str(),NDSKEL_VTK_ASCII|NDSKEL_VTK_XML);
    }
  };


  class fromSegs : public interface {
    std::string getTypeStr() {return std::string("segs_ascii");}
    std::string getExtensionStr() {return std::string("a.segs");}
    bool canLoad(std::string fname) 
    {
      return false;
    }
    bool canSave() {return true;}
    
    NDskel *load(std::string fname)
    {
      return NULL;
    }
    
    int save(NDskel *skl, std::string fname)
    {
      FILE *f;
      long i,j;
      int seg_val1=-1;
      int seg_val2=-1;
      int seg_type=-1;
      int seg_orient=-1;
      int node_val=-1;

      seg_val1=getDataFieldID(skl,1,SEG_P1_TAG(VALUE_TAG));
      if (seg_val1<0) fprintf(stderr,"WARNING(IGNORED): tag '%s' not found.\n",SEG_P1_TAG(VALUE_TAG));    
      seg_val2=getDataFieldID(skl,1,SEG_P2_TAG(VALUE_TAG));
      if (seg_val2<0) fprintf(stderr,"WARNING(IGNORED): tag '%s' not found.\n",SEG_P2_TAG(VALUE_TAG));
      seg_type=getDataFieldID(skl,1,TYPE_TAG);
      if (seg_type<0) fprintf(stderr,"WARNING(IGNORED): tag '%s' not found.\n",TYPE_TAG);
      seg_orient=getDataFieldID(skl,1,ORIENTATION_TAG);
      //if (seg_orient<0) fprintf(stderr,"WARNING: tag '%s' not found.\n",ORIENTATION_TAG);
      node_val=getDataFieldID(skl,0,VALUE_TAG);
      if (node_val<0) 
	{
	  fprintf(stderr,"ERROR: node tag '%s' not found.\n",VALUE_TAG);
	  return -1;
	}
      
      f=fopen(fname.c_str(),"w");
      if (f==NULL) return -1;
      fprintf(f,"#arc segments\n");
      fprintf(f,"#U0");
      for (i=1;i<skl->ndims;i++) fprintf(f," U%ld",i);
      for (i=0;i<skl->ndims;i++) fprintf(f," V%ld",i);
      fprintf(f," value_U value_V type boundary\n");
      fprintf(f,"#%ld %ld\n",(long)skl->ndims,(long)skl->nsegs);
      for (i=0;i<skl->nsegs;i++)
	{
	  NDskl_seg *seg=&skl->Seg[i];
	  int type;
	  int dir;
	  
	  if (seg_type<0)
	    type=(skl->Node[seg->nodes[0]].type<skl->Node[seg->nodes[1]].type)?skl->Node[seg->nodes[0]].type:skl->Node[seg->nodes[1]].type;
	  else type=(int)seg->data[seg_type];
	  /*
	  if (seg_orient<0)
	    dir=(skl->Node[seg->nodes[0]].data[node_val]<skl->Node[seg->nodes[1]].data[node_val]);
	  else dir=(int)seg->data[seg_orient];
	  */
	  if (seg_orient<0) dir=1;
	  else dir=(int)seg->data[seg_orient];

	  if (dir>0)
	    {
	      fprintf(f,"%g",seg->pos[0]);
	      for (j=1;j<2*skl->ndims;j++) fprintf(f," %g",seg->pos[j]);
	      fprintf(f," %g %g",
		      (seg_val1<0)?0:seg->data[seg_val1],
		      (seg_val2<0)?0:seg->data[seg_val2]);
	      fprintf(f," %d %d\n",type,seg->flags);
	    }
	  else // segments are always oriented in the direction of increasing value.
	    {
	      //printf("REV\n");
	      fprintf(f,"%g",seg->pos[skl->ndims]);
	      for (j=skl->ndims+1;j<2*skl->ndims;j++) fprintf(f," %g",seg->pos[j]);
	      for (j=0;j<skl->ndims;j++) fprintf(f," %g",seg->pos[j]);
	      fprintf(f," %g %g",
		      (seg_val2<0)?0:seg->data[seg_val2],
		      (seg_val1<0)?0:seg->data[seg_val1]);
	      fprintf(f," %d %d\n",type,seg->flags);
	    }
	}
	
      fclose(f);
      return 0;
    }
  };

  class fromCrits : public interface {
    std::string getTypeStr() {return std::string("crits_ascii");}
    std::string getExtensionStr() {return std::string("a.crits");}
    bool canLoad(std::string fname) 
    {
      return false;
    }
    bool canSave() {return true;}
    
    NDskel *load(std::string fname)
    {
      return NULL;
    }
    
    int save(NDskel *skl, std::string fname)
    {
      FILE *f;
      long i,j;
      int node_val=-1;
      int node_pair=-1;
            
      node_val=getDataFieldID(skl,0,VALUE_TAG);
      if (node_val<0) fprintf(stderr,"WARNING: tag '%s' not found.\n",VALUE_TAG);
      node_pair=getDataFieldID(skl,0,PERSISTENCE_PAIR_TAG);
      if (node_val<0) fprintf(stderr,"WARNING: tag '%s' not found.\n",VALUE_TAG);
     
      f=fopen(fname.c_str(),"w");
      if (f==NULL) return -1;
      fprintf(f,"#critical points\n");
      fprintf(f,"#X0");
      for (i=1;i<skl->ndims;i++) fprintf(f," X%ld",i);
      fprintf(f," value type pair_id boundary\n");
      fprintf(f,"#%ld %ld\n",(long)skl->ndims,(long)skl->nnodes);
      for (i=0;i<skl->nnodes;i++)
	{
	  NDskl_node *node=&skl->Node[i];
	  
	  fprintf(f,"%g",node->pos[0]);
	  for (j=1;j<skl->ndims;j++) fprintf(f," %g",node->pos[j]);
	  fprintf(f," %g %d %ld %d\n",
		  (node_val<0)?0:node->data[node_val],
		  node->type,
		  (node_pair<0)?i:(long)node->data[node_pair],
		  node->flags);
	}
	
      fclose(f);
      return 0;
    }
  };

  class fromNDnet : public interface {
    std::string getTypeStr() {return std::string("NDnet");}
    std::string getExtensionStr() {return std::string("NDnet");}
    bool canLoad(std::string fname) {return false;}
    bool canSave() {return true;}
    NDskel *load(std::string fname) {return NULL;}
    
    int save(NDskel *skl, std::string fname)
    {
      NDnetwork *net = NDskel2NDnet(skl);
      Save_NDnetwork(net,fname.c_str());
      FreeNDnetwork(&net);
      return 0;
    }
  };

  struct allTypes
  {
    template <class OutputIterator>
    static void generate(OutputIterator out)
    {
      *out = new fromNDskel();
      *out = new fromNDskel_ascii();
      *out = new fromNDnet();
      *out = new fromSegs();
      *out = new fromCrits();
      *out = new fromVTK();
      *out = new fromVTK_ascii();
      *out = new fromVTK_xml();
      *out = new fromVTK_xml_ascii();
    }

    static std::string type()
    {
      return std::string("skeleton");
    }
  };

  typedef genericIO<interface, allTypes> IO;
}

#endif

