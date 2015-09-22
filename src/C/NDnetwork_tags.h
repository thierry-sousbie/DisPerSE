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
#ifndef __NDNETWORK_TAGS_H__
#define __NDNETWORK_TAGS_H__

//#include "patches.h"
#include "NDnetwork.h"
#include "NDfield.h"

#ifdef __cplusplus
extern "C" {
#endif

#define GROUP_ID_PEAK_TAG "Peak group ID"
#define GROUP_ID_VOID_TAG "Void group ID"
#define PEAK_PATCH_INDEX_TAG "peak patch index"
#define VOID_PATCH_INDEX_TAG "void patch index"
#define PEAK_PATCH_PROBA_TAG "peak patch proba"
#define VOID_PATCH_PROBA_TAG "void patch proba"
#define PEAK_PATCH_INDEX_GUESS_TAG "guessed peak patch index"
#define VOID_PATCH_INDEX_GUESS_TAG "guessed void patch index"
#define DENSITY_TAG "Density"
#define VOLUME_TAG "Volume"
#define NPATCHES_TAG "NPatches"
  //#define VALUE_TAG "Field value"

  //#define LOG_TAG(tag) "Log " tag

    int NDDataIndex(NDnetwork *net, int type, const char *name);
    int NDSupDataIndex(NDnetwork *net, int type, const char *name);

    
    NDnetwork_Data *addNDData(NDnetwork *net, int type, const char *name, 
			      double(*compute_val)(NDnetwork* net,int type,int index,void* info),
			      void *info);
    double *addNDDataArr(NDnetwork *net, int type, const char *name, double **data_p);
    void *addNDSupDataArr(NDnetwork *net, int type, const char *name, 
			  int datasize, const char *datatype, void **supData_p);
  int AddNDDataTagFromGrid(NDnetwork *net,NDfield *field,const char *name);

    

    //  int TagNetworkPatches(NDnetwork *net,int type, NetPatch_type **vpatches,
//			  NetPatch_type **ppatches,int storeProba, int storePatches);
    int NDNetTagCritical(NDnetwork *net);
//  int TagNetworkPatches(NDnetwork *net,int type, NetPatch_type **vpatches,NetPatch_type **ppatches,int storeProba, int storePatches);

    double NDNetPatchProbaTag(NDnetwork *net, int type, int index, void *info);
    double NDNetPatchTag(NDnetwork *net, int type, int index, void *info);

    double NDLogOf(NDnetwork *net, int type, NDNET_UINT index, void *info);
    double NDVolume(NDnetwork *net, int type, NDNET_UINT index, void *info);
    double NDDensity(NDnetwork *net, int type, NDNET_UINT index, void *info);
//  double NDNetNPatchesTag(NDnetwork *net, int type, int index, void *info);

    double NDVertexDensity(NDnetwork *net, NDNET_UINT index, double *volume);
    double NDFaceDensity(NDnetwork *net, int type, NDNET_UINT index);
    

    NDnetwork_Data *SmoothNDData(NDnetwork *net, int type,const char *name, uint ntimes);


#ifdef __cplusplus
}
#endif

#endif
