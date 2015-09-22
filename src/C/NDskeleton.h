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
#ifndef __ND_SKELETON_H__
#define __ND_SKELETON_H__

#ifdef __cplusplus
extern "C" {
#endif

#define NDSKEL_DATA_STR_SIZE 20
#define NDSKEL_TAG "NDSKEL"
#define NDSKEL_ASCII_TAG "ANDSKEL"
#define NDSKEL_MAX_DIMS 20

#define FLAG_NDNODE_BOUNDARY (1<<0)
#define FLAG_NDNODE_OUT      (1<<1)
#define FLAG_NDNODE_INFINITE (1<<2)

  /* 
#define FLAG_NDSKL_SAD      (1<<0)
#define FLAG_NDSKL_MIN      (1<<1)
#define FLAG_NDSKL_MAX      (1<<2)
#define FLAG_NDSKL_MINMAX   (1<<3)
#define FLAG_NDSKL_BIF      (1<<4)
#define FLAG_NDSKL_FAKE_SAD (1<<5)
#define FLAG_NDSKL_LOOP     (1<<31)
  */

  struct NDskl_seg_str {
    int nodes[2];
    float *pos;
    int flags;
    int index;
    double *data;
    struct NDskl_seg_str *Next;
    struct NDskl_seg_str *Prev;
  }; 
  typedef struct NDskl_seg_str NDskl_seg;
  
  struct NDskl_node_str {
    float *pos;
    int flags;
    int nnext;
    int type;
    int index;
    int *nsegs; 
    double *data;
    struct NDskl_node_str **Next;  
    struct NDskl_seg_str **Seg;
  }; 
  typedef struct NDskl_node_str NDskl_node;
  
  typedef struct NDskel_str {
    char comment[80];
    
    int ndims;
    int *dims;
    double *x0;
    double *delta;
    
    int nsegs;
    int nnodes;
    
    int nsegdata;
    int nnodedata;
    char **segdata_info;
    char **nodedata_info;
    
    float *segpos;
    float *nodepos;
    double *segdata;
    double *nodedata;
    
    NDskl_seg *Seg;
    NDskl_node *Node;
  } NDskel;

  int NDskel_SegDataIndex(NDskel *skl, const char *name);
  int NDskel_NodeDataIndex(NDskel *skl, const char *name);
  
  int NDskel_realloc(NDskel *skl,int nsegs,int nnodes,int nsegalloc,int nnodealloc, int delta);
  int NDskel_SetNNextInNode(NDskl_node *node,int nnodes);
  NDskel *Copy_NDskel(NDskel *skl);
  NDskel *Create_NDskel(int *dims,int ndims,double *x0, double *delta,const char *comment,int nnodes,int nsegs);
  int Save_NDskel(NDskel *skl,const char *filename);
  int Save_ASCIIskel(NDskel *skl,const char *filename);
  int IsNDskeleton(const char *filename);
  NDskel *Load_NDskelHeader(const char *filename);
  NDskel *Load_NDskel(const char *filename);
  int Free_NDskeleton(NDskel **skl);
  int NDskelCheckSanity(NDskel *skl,int periodic);
  NDskel *NDskel_Subregion(NDskel *skl_p, double margin, int **seglist, long *nsegs, int inplace);
  
  double ComputeSegLen(NDskel *skl, NDskl_seg* seg);
  double ComputeDistFromPrev(NDskel *skl, NDskl_seg* seg);
  double ComputeDistToNext(NDskel *skl, NDskl_seg* seg);
  double ComputeDistToSaddle(NDskel *skl, NDskl_seg* seg_p);
  int SortNDskelSegments(NDskel *skl,char *nodefieldname,int periodic);
  
  long getNDskelFilTab(NDskel *skl, NDskl_seg ***filTab, int **filSize);
  long getNDskelFilTabInfo(NDskel *skl, NDskl_seg **filTab, int nFil, char ***dataName, double ***data);
  void freeNDskelFilTab(NDskl_seg ***filTab, int **filSize);
  void freeNDskelFilTabInfo(char ***dataName, double ***data, int nData, int nFil);
  int printNDskelStat(NDskel *skl, int dec);

#ifdef HAVE_CFITS_IO
  int Save_FITSskel(NDskel *skl,const char *filename, long *naxes_p);
#endif

#ifdef __cplusplus
}
#endif

#endif
