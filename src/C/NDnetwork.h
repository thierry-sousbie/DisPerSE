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
#ifndef __ND_NETWORK_H__
#define __ND_NETWORK_H__

#include "mytypes.h"

#ifdef __cplusplus
extern "C" {
#endif

  // IMPORTANT Support for NDNET_FLOAT = double not implemented yet !!!

#define FWRITE_DUMMY_BFR_SIZE    (1<<20)
#define NDNET_UINT                UINT
#define NDNET_INT                 INT
#define NDNET_FLOAT               FLOAT
#define NDNET_IDCUMT              long

#define NDNETFLAG_OUT            (1<<0)
#define NDNETFLAG_PERIODIC_CUT   (1<<1)
#define NDNETFLAG_SHARED         (1<<3)
#define NDNETFLAG_KEEP	         (1<<4)
#define NDNETFLAG_TAG	         (1<<5)
#define NDNETFLAG_NON_MANIFOLD	 (1<<6)
#define NDNETFLAG_BOUNDARY	 (1<<7)

#define NDNETWORK_DATA_STR_SIZE  16
#define NDNETWORK_TAG            "NDNETWORK"
#define NDNETWORK_ASCII_TAG      "ANDNET"
#define NDNETWORK_ASCII_ADDTAG   "[ADDITIONAL_DATA]"
#define NDNETWORK_V_FLAG_TAG     "internal_vertex_flags"
#define NDNETWORK_F_FLAG_TAG     "internal_face_flags"

// returns the index corresponding to faces of dimension ndims
//#define NDIMS_TYPE_INDEX(net,ndims,type) {for (type=0;type<net->f_numDims;type++) {if (net->f_ndims[type]==ndims) break; } if (type>=net->f_numDims) type=-1;}

// The number and list of vertices/faces a face/vertex belongs to.
#define NUM_FACE_IN_VERTEX(net,type,vertex) (net->v_numFaceIndexCum[type][vertex+1]-net->v_numFaceIndexCum[type][vertex])
#define FACE_IN_VERTEX(net,type,vertex) (&(net->v_faceIndex[type][net->v_numFaceIndexCum[type][vertex]]))

#define NUM_VERTEX_IN_FACE(net,type,face) ((net->isSimpComplex)?(type+1):(net->f_numVertexIndexCum[type][face+1]-net->f_numVertexIndexCum[type][face]))
#define VERTEX_IN_FACE(net,type,face) (&(net->f_vertexIndex[type][(net->isSimpComplex)?((long)face*((long)type+1)):(net->f_numVertexIndexCum[type][face])]))

#define NUM_FACE_IN_FACE(net,ref_type,ref_face,type) (net->f_numFaceIndexCum[ref_type][type][ref_face+1]-net->f_numFaceIndexCum[ref_type][type][ref_face])
#define FACE_IN_FACE(net,ref_type,ref_face,type) (&(net->f_faceIndex[ref_type][type][net->f_numFaceIndexCum[ref_type][type][ref_face]]))

// this is quite useless, for speed only ...
/*
#define NUM_VERTEX_IN_FACE_SIMPLEX(net,type,face) (net->f_ndims[type]+1)
#define VERTEX_IN_FACE_SIMPLEX(net,type,face) (&(net->f_vertexIndex[type][(face*(net->f_ndims[type]+1))]))
#define NUM_VERTEX_IN_FACE_NONSIMPLEX(net,type,face) (net->f_numVertexIndexCum[type][face+1]-net->f_numVertexIndexCum[type][face])
#define VERTEX_IN_FACE_NONSIMPLEX(net,type,face) (&(net->f_vertexIndex[type][(net->f_numVertexIndexCum[type][face]])))
*/

//#define uint unsigned int 

  typedef struct
  {
    int type; // the cell-type
    char name[255];  // name of the field
    double *data;  // value for each of the nfaces[type] type-faces
  } NDnetwork_Data;

typedef struct
{
    int type; //0 for vertex, t for faces of type (i.e. dimension) t.
    char name[255];
    int datasize;
    char datatype[255];// a string to identity how data should be casted
    void *data;
} NDnetwork_SupData;

 typedef struct 
 {
    char comment[80];
    int periodicity;
    int ndims; // the number of spatial dimensions
    int ndims_net; // number of dimension of the network itself (e.g. 2 for a sphere embedded in 3D)
    int isSimpComplex;  // 1 if network is a simplicial complex (always true in disperse)
    double *x0;  // origin of the bounding box
    double *delta;  // size of the bounding box
    int indexSize; // size of NDNET_UINT type in Bytes
    int cumIndexSize; // size of NDNET_IDCUMT type in Bytes
    int floatSize; // size of NDNET_IDCUMT type in Bytes
    char dummy[160-4*3]; // dummy data reserved for future extensions
  
    NDNET_UINT nvertex;  // total number of vertices
    float *v_coord; //vertices coodinates (X_0,Y_0,Z_0,X_1,Y_1,...,Z_nvertex-1)
    
    NDNET_UINT *nfaces; // number of faces of a given type t is given by nfaces[t]
  
    int *haveVertexFromFace; // haveVertexFromFace[n] is 1 if we have an explicit definition of the n-faces (at least one type of face must be defined).
    NDNET_IDCUMT **f_numVertexIndexCum;// cumulative number of vertice in the faces f of type t, NULL for simplicial faces
    NDNET_UINT **f_vertexIndex; // list of vertices defining the n-faces is stored in f_vertexIndex[n], all vertices being enumerated for each face (the indices of the vertices in the kth n-face start at f_vertexIndex[n][(n+1)*k] )
                                // see also macro  NUM_VERTEX_IN_FACE(net,type,face) and VERTEX_IN_FACE(net,type,face)
  
    //This may be computed internally within DisPerSE but does not need to be defined explicitely
    int *haveFaceFromVertex; // haveFaceFromVertex[n] is 1 if we have an explicit list of all the n-faces that contain each vertex (used to navigate within the network)
    NDNET_IDCUMT **v_numFaceIndexCum; // cumulative number of faces of type t a vertex v belongs to
    NDNET_UINT **v_faceIndex; // indices of the faces of type t a vertex v belongs to ( the list of n-faces of vertex k starts at v_faceIndex[n][net->v_numFaceIndexCum[n][k]] and ends at v_faceIndex[n][net->v_numFaceIndexCum[n][k+1]] )
                                // see also macro  NUM_FACE_IN_VERTEX(net,type,vertex) and  FACE_IN_VERTEX(net,type,vertex)
    
    // This can become extremely memory heavy ... NOT used in DisPerSE
    int **haveFaceFromFace; // haveFaceFromFace[k][n] is 1 if we have an explicit list of all the n-faces that have a boundary/co-boundary relation with each k-face (used to navigate within the network)
    NDNET_IDCUMT ***f_numFaceIndexCum; //  cumulative number of n-faces having a boundary / co-boundary relation with each k-face f_numFaceIndexCum[k][n]
    NDNET_UINT ***f_faceIndex; // indices of the faces (similar to v_faceIndex)
                               // see also macro NUM_FACE_IN_FACE(net,ref_type,ref_face,type) and FACE_IN_FACE(net,ref_type,ref_face,type)
    
    int haveVFlags;  // do we have flags associated to each vertex ?
    int *haveFFlags;  // do we have flags associated to each n-face ?
    unsigned char *v_flag; // nvertex flag values (1 for each vertex) or NULL 
    unsigned char **f_flag; // nfaces[n] flag values (1 of each n-face) or NULL
 
    int ndata; // number of additional data fields.
    NDnetwork_Data *data; // array of all additionnal data (data in total)
    
    int nsupData;
    NDnetwork_SupData *supData;
 
 } NDnetwork;
  
  int IsNDnetwork(const char *filename);
  void FreeNDnetwork(NDnetwork **net_p);
  int ComputeFaceFromVertex(NDnetwork *net, int type);
  int ComputeVertexFromFace(NDnetwork *net, int type);
  int FacesInFace(NDnetwork *net, NDNET_UINT ref_face, int ref_type,int type, NDNET_UINT **list);
  int NDFindNeighbours(NDnetwork *net, int type, NDNET_UINT index, NDNET_UINT **list, double **interface);
  void NDnet_smooth(NDnetwork *net, int n);
  int NDNetFlags_enable(NDnetwork *net,int type);

  NDnetwork *CreateNetwork(int ndims, int nvertex, int vflags);

  int Save_NDnetwork(NDnetwork *net,const char *filename);
  NDnetwork *Load_NDnetwork(const char *filename);
  int Save_NDnetwork_ASCII(NDnetwork *net, const char *filename);
  NDnetwork *Load_NDnetwork_ASCII(const char *filename);
   
  int Simplex_COM(NDnetwork *net,int type, long index, float *result);
  int Simplex_CS(NDnetwork *net,int type, long index, float *result);
  int Simplex_CS_fc(NDnetwork *net,int type, double **pv, float *result);

  int tagBoundaries(NDnetwork *net);
  int printNDnetStat(NDnetwork *net, int dec);
  
  int getSimplexCoords(NDnetwork *net,int type, NDNET_UINT index, float *coords, int unfold);
  int setPeriodicCutFlag(NDnetwork *net);

 
#define NDNET_FSTR_VAR(NAME) fpos_t NAME;	\
  size_t s_ ## NAME;				\
  size_t rn_ ## NAME;				\
  size_t rs_ ## NAME;				\
  size_t rss_ ## NAME;

#define NDNET_FSTR_W_VAR(NAME) fpos_t NAME;	\
  size_t rn_ ## NAME;				\
  size_t rs_ ## NAME;				\

  
  typedef struct NDnet_fstruct_info_str {
    int swap;
    NDNET_FSTR_VAR(v_coord)
    NDNET_FSTR_VAR(f_vertexIndex[16])
    NDNET_FSTR_VAR(f_numVertexIndexCum[16])
    NDNET_FSTR_VAR(v_numFaceIndexCum[16])
    NDNET_FSTR_VAR(v_faceIndex[16])
    NDNET_FSTR_VAR(f_numFaceIndexCum[16][16])
    NDNET_FSTR_VAR(f_faceIndex[16][16])
    NDNET_FSTR_VAR(v_flag)
    NDNET_FSTR_VAR(f_flag[16])
    NDNET_FSTR_VAR(data[1024])
    NDNET_FSTR_VAR(supData[1024])
  } NDnet_fstruct_info;

  typedef struct NDnet_fstruct_info_str NDnet_fstruct_w_info;
  /*
  typedef struct NDnet_fstruct_w_info_str {
    NDNET_FSTR_W_VAR(v_coord)
    NDNET_FSTR_W_VAR(f_vertexIndex[16])
    NDNET_FSTR_W_VAR(f_numVertexIndexCum[16])
    NDNET_FSTR_W_VAR(v_numFaceIndexCum[16])
    NDNET_FSTR_W_VAR(v_faceIndex[16])
    NDNET_FSTR_W_VAR(f_numFaceIndexCum[16][16])
    NDNET_FSTR_W_VAR(f_faceIndex[16][16])
    NDNET_FSTR_W_VAR(v_flag)
    NDNET_FSTR_W_VAR(f_flag[16])
    NDNET_FSTR_W_VAR(data[1024])
    NDNET_FSTR_W_VAR(supData[1024])
  } NDnet_fstruct_w_info;
  */
  
  
  NDnetwork *Load_NDnetwork_header(const char *filename, NDnet_fstruct_info *info);
  void NDnetReadData(void *res,FILE *f,fpos_t *where,size_t rs, size_t rn, size_t rss, int swap);
  void NDnetWriteData(void *data,FILE *f,fpos_t *where,size_t rs, size_t rn, long count);
  
  size_t fwrite_dummy(size_t size, size_t nmemb,FILE *stream);
  NDnet_fstruct_w_info *Save_NDnetwork_header(NDnetwork *net,const char *filename, FILE **f);


#define NDNET_READ_DATA(RESULT,FF,STRCT,WHAT)				\
  NDnetReadData(RESULT,FF,&(STRCT)->WHAT,(STRCT)->rs_ ## WHAT,(STRCT)->rn_ ## WHAT,(STRCT)->rss_ ## WHAT,(STRCT)->swap)

#define NDNET_WRITE_DATA(DATA,FF,STRCT,WHAT,COUNT)				\
  NDnetWriteData(DATA,FF,&(STRCT)->WHAT,(STRCT)->rss_ ## WHAT,(STRCT)->rn_ ## WHAT, COUNT)
  
#ifdef __cplusplus
}
#endif

#endif
