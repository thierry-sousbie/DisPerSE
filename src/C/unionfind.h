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
#ifndef __UNION_FIND_H__
#define __UNION_FIND_H__

// A union-find data structure
// groups indice start at 0.

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct {
    int NNodes;
    int *parent;
    int *size;
    int preserveID;
  } UF_type;

  typedef struct {
    long NNodes;
    long *parent;
    long *size;
    int preserveID;
  } UFL_type;

  int uf_free(UF_type **uf);
  UF_type *uf_create(int N,int preserveID);
  int uf_isOwnGroup(UF_type *uf, long node);
  int uf_find(UF_type *uf, int node);
  void uf_union(UF_type *uf,int node1, int node2);
  //void uf_union_preserve(UF_type *uf,int node1, int node2);
  void uf_gr_union(UF_type *uf,int set1, int set2);
  //void uf_gr_union_preserve(UF_type *uf,int set1, int set2);

  int ufl_free(UFL_type *uf);
  UFL_type *ufl_create(long N,int preserveID);
  int ufl_isOwnGroup(UFL_type *uf, long node);
  long ufl_find(UFL_type *uf, long node);
  void ufl_union(UFL_type *uf,long node1, long node2);
  //void ufl_union_preserve(UFL_type *uf,long node1, long node2);
  void ufl_gr_union(UFL_type *uf,long set1, long set2);
  //void ufl_gr_union_preserve(UFL_type *uf,long set1, long set2);

#ifdef __cplusplus
}
#endif

#endif
