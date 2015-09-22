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
#ifndef __FOF_STRUCT_H__
#define __FOF_STRUCT_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef struct group_type_str {
    int N;//Numberof particles in the cluster
    int cut;//tells if the cluster is on the edge of the box
    int *index;//indices of the N particles 
    float *Center;//center of the cluster
    float *Spin;//spin of the cluster
    //float xmin,xmax,ymin,ymax,zmin,zmax;
    float *VMean;//Speed of the cluster
    //float VSig[6]; // xx,xy,xz,yy,yz,zz
    
    float sig_evec[9];//Velocity dispertion eigenvectors [V1x,V1y,V1z,V2x,...]
    float sig_eval[3];//Velocity dispersion eigenvalues [l1,l2,l3]
    
} group_type;

typedef struct group_list_str {
    group_type *group;// All the groups informations
    
    float *Pos;//position of the particles
    float *Vel;//velocity of the particles
  
    int NGroups;//number of groups
    int NPart;//number of particles
    int NMin;//Minimum number of particles in a group
    float ll;//linking length
    int ndims;
} group_list;

#ifdef __cplusplus
}
#endif

#endif
