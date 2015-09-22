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
#ifndef __HYPERCUBE_H__
#define __HYPERCUBE_H__

#ifdef __cplusplus
extern "C" {
#endif


#include "mytypes.h"

struct nface_str
{
    int ndims;
    INT tag;
    int nvertice;
    int *vertex;
    
    int nmothers;
    int ndaughters;
    
    struct nface_str **Mothers;
    struct nface_str **Daughters;
    
}; 

typedef struct nface_str nface_type;

struct nfacelist_str
{
    int ndims;
    int nfaces;
    int nFacesPerCube;
    nface_type *NFace;

    struct nfacelist_str *mother;
    struct nfacelist_str *daughter;
};

typedef struct nfacelist_str nfacelist_type;

typedef struct hypercube_str
{
  int ndims;
  nface_type **faceFromTag;
  nfacelist_type *NFaceList;

} hypercube_type;

    hypercube_type *GenerateHypercube(int ndims);
    unsigned int HC_getCubeIndex(hypercube_type *h, char type, unsigned int id,int periodic);
    unsigned int HC_getIndex(hypercube_type *h,const nface_type *face, const unsigned int *cube_index, const int *dims, int periodic);
    nface_type *HC_getFace(hypercube_type *h, char type, unsigned int id,const int *dims,int periodic);
    void HC_getSymetricConfig(hypercube_type *h,const nface_type *face,const unsigned int *cube_coords,nface_type **new_face,unsigned int *new_cube_index,const int *dims,int periodic);

    void FreeHypercube(hypercube_type *h);

#ifdef __cplusplus
}
#endif

#endif
