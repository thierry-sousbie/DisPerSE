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
#ifndef __PD_PARAMS_HXX__
#define __PD_PARAMS_HXX__

#include "NDnetwork.h"
#include "NDnetwork_tags.h"
#include "global.h"

#include <vector>
#include <list>

#include <mgl/mgl_data.h>
#include <string>


struct PD_params
{
  std::vector<bool> showP;
  std::list<int> shownP;

  double x1;
  double y1;
  double x2;
  double y2;

  std::string x_label;
  std::string y_label;

  bool hideP;
  bool logX,logY;
  bool logH;
  bool showHisto;
  NDnetwork *net;

  mglData H;
  mglData Hx;
  mglData Hy;

  mglData v_V;
  mglData v_T;

  mglData V[4];
  mglData P[4];
  mglData Pr[4];

  mglData *x_data;
  mglData *y_data;

  int ndims;
  /*
  mglData *x_v_data;
  mglData *y_v_data;
  */
  PD_params()
  {
    net=NULL;
    x_data=y_data=NULL;
    showP.assign(4,true);
    showHisto=false;
    logH=true;
    hideP=false;
    ndims=0;
    for (unsigned long i=0;i<showP.size();i++)
      shownP.push_back(i);

    x_label=std::string("Value");
    y_label=std::string("");
  }

  ~PD_params()
  {
    FreeNDnetwork(&net);
  }

  std::pair<bool,bool> setNetwork(NDnetwork *_net)
  {
    std::string tag;
    int id;
    unsigned long i;
    bool plotR=false;
    
    double *vertV=NULL;
    tag=std::string(VALUE_TAG);
    id = NDDataIndex(_net,0,tag.c_str());
    if (id<0) 
      {
	fprintf(stderr,"WARNING: field '%s' not found.\n",tag.c_str());
      }
    else 
      {
	vertV=_net->data[id].data;
	v_V.Set(_net->data[id].data,_net->nvertex);
      }

    tag=std::string(TYPE_TAG);
    id = NDDataIndex(_net,0,tag.c_str());
    
    if (id<0) 
      {
	fprintf(stderr,"WARNING: field '%s' not found.\n",tag.c_str());
      }
    else 
      {
	ndims=(int)*std::max_element(_net->data[id].data,&_net->data[id].data[_net->nvertex]);
	v_T.Set(_net->data[id].data,_net->nvertex);
      }

    bool havePr=false;
    double *segT=NULL;
    tag=std::string(TYPE_TAG);
    id = NDDataIndex(_net,1,tag.c_str());
    if (id<0) fprintf(stderr,"WARNING: field '%s' not found.\n",tag.c_str());
    else segT=_net->data[id].data;
    
    if (segT!=NULL)
      {
	tag=std::string(PERSISTENCE_TAG);
	id = NDDataIndex(_net,1,tag.c_str());
	if (id<0) fprintf(stderr,"WARNING: field '%s' not found.\n",tag.c_str());
	else 
	  {	    
	    std::vector<double> tmp[4];
	    double *val=_net->data[id].data;
	    //v_P.Set(val,_net->nfaces[1]);
	    for (i=0;i<_net->nfaces[1];i++)  tmp[(int)segT[i]-1].push_back(val[i]);
	    for (i=0;i<4;i++) P[i].Set(&tmp[i][0],tmp[i].size());
	  }
	
	tag=std::string(PERSISTENCE_RATIO_TAG);
	id = NDDataIndex(_net,1,tag.c_str());
	if (id<0)
	  {
	    tag=std::string("persistence Ratio");
	    id = NDDataIndex(_net,1,tag.c_str());
	  }

	if (id<0) 
	  {
	    fprintf(stderr,"INFO: field '%s' not available.\n",tag.c_str());
	  }
	else 
	  {
	    std::vector<double> tmp[4];
	    double *val=_net->data[id].data;
	    //v_Pr.Set(val,_net->nfaces[1]);
	    for (i=0;i<_net->nfaces[1];i++) tmp[(int)segT[i]-1].push_back(val[i]);
	    for (i=0;i<4;i++) Pr[i].Set(&tmp[i][0],tmp[i].size());
	    havePr=true;
	  }
      

	if (vertV!=NULL)
	  {
	    std::vector<double> val(_net->nfaces[1]);
	    for (i=0;i<_net->nfaces[1];i++)
	      {
		NDNET_UINT *v=VERTEX_IN_FACE(_net,1,i);
		val[i]=std::min(vertV[v[0]],vertV[v[1]]);
	      }
	    std::vector<double> tmp[4];
	    for (i=0;i<_net->nfaces[1];i++) tmp[(int)segT[i]-1].push_back(val[i]);
	    for (i=0;i<4;i++) V[i].Set(&tmp[i][0],tmp[i].size());	   
	  }    
      }

    x_data=V;

    if (0*havePr)
      {
	y_data=Pr;
	logX=true;
	logY=true;
	plotR=true;
	//y_label=std::string("Persistence ratio");
      }
    else 
      {
	y_data=P;
	logX=false;
	logY=false;
	plotR=false;
	//y_label=std::string("Persistence");
      }    

    /*
    logX=true;
    logY=true;
    */
    FreeNDnetwork(&net);
    net=_net; 
     
    return std::make_pair(true,plotR);
  }

};


#endif
