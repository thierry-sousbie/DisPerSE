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
#ifndef __ND_FIELD_H__
#define __ND_FIELD_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "mytypes.h"

#define NDFIELD_MAX_DIMS 20
#define NDFIELD_TAG "NDFIELD"
#define NDFIELD_ASCII_TAG "ANDFIELD"

#define ND_CHAR   (1<<0)
#define ND_UCHAR  (1<<1)
#define ND_SHORT  (1<<2)
#define ND_USHORT (1<<3)
#define ND_INT    (1<<4)
#define ND_UINT   (1<<5)
#define ND_LONG   (1<<6)
#define ND_ULONG  (1<<7)
#define ND_FLOAT  (1<<8)
#define ND_DOUBLE (1<<9)

#define TESTPAT 0xaaaaaaaa

#define DEFINE_NDVARS(varname)\
  unsigned char *varname##_uc=NULL;\
  char *varname##_c=NULL;\
  unsigned short *varname##_us=NULL;\
  short *varname##_s=NULL;\
  unsigned int*varname##_ui=NULL;\
  int *varname##_i=NULL;\
  unsigned long *varname##_ul=NULL;\
  long *varname##_l=NULL;\
  float *varname##_f=NULL;\
  double *varname##_d=NULL;\

#define SETNDPOINTER(var,voidptr,type)  switch (type)\
    {\
    case ND_UCHAR: var##_uc = (unsigned char *)voidptr;break;	\
    case ND_CHAR: var##_c = (char *)voidptr;break;		\
    case ND_USHORT: var##_us = (unsigned short *)voidptr;break;	\
    case ND_SHORT: var##_s = (short *)voidptr;break;		\
    case ND_UINT: var##_ui = (unsigned int *)voidptr;break;	\
    case ND_INT: var##_i  = (int *)voidptr;break;		\
    case ND_ULONG:var##_ul = (unsigned long *)voidptr;break;	\
    case ND_LONG: var##_l = (long *)voidptr;break;		\
    case ND_FLOAT: var##_f = (float *)voidptr;break;		\
    case ND_DOUBLE: var##_d = (double *)voidptr;break;		\
    }\

#define SETNDFIELD_VAL(var,val,index,type)  switch (type)\
    {\
    case ND_UCHAR: var= val##_uc [index];break;\
    case ND_CHAR: var= val##_c [index];break;\
    case ND_USHORT: var= val##_us [index];break;\
    case ND_SHORT: var= val##_s [index];break;\
    case ND_UINT: var= val##_ui [index];break;\
    case ND_INT: var= val##_i [index] ;break;\
    case ND_ULONG:var= val##_ul [index];break;\
    case ND_LONG: var= val##_l [index];break;\
    case ND_FLOAT: var= val##_f [index];break;\
    case ND_DOUBLE: var= val##_d [index];break;\
    }\



typedef struct NDfield_str
{
  char comment[80];
  int dims[NDFIELD_MAX_DIMS];
  int ndims;
  int n_dims;
  int fdims_index; //index where dims are the fields dims
  int datatype;
  double x0[NDFIELD_MAX_DIMS];
  double delta[NDFIELD_MAX_DIMS];
  char dummy[160];
  void *val;

  long nval;
  int datasize;
} NDfield;

typedef struct OLD_density_grid_str
{
    int Nx,Ny,Nz;//Number of nodes along x,y and z.
    INT NNodes;//number of nodes

    float x0,y0,z0;//start coordinates
    float dx,dy,dz;//grid spacing

    FLOAT *grid;//value at every node ...

    int HasGradiant;//is gradiant computed ?
    FLOAT *grad;//gradient at every node (gx1,gy1,gz1,gx2,gy2,gz2,...)
    int HasHessian;
    FLOAT *hessian;//Hessian (xx1,xy1,xz1,yy1,yz1,zz1,xx2 ...)
    int maxindex;//index of the maximum value of grid (computed with the gradiant)

  float redshift;
  float smoothing;
} OLD_density_grid;

int Free_NDfield(NDfield **field);
int Get_NDtype(int size,int is_integer, int not_signed);
int sizeof_NDfield(int type);
NDfield *Create_NDfield(int *dims,int ndims,int fdims_index,int type,double *x0, double *delta, void *data,const char *comment);
int Init_NDfield(NDfield *field,int *dims,int ndims,int fdims_index,int type,double *x0, double *delta, void *data,char *comment);
int Free_NDfield(NDfield **field);
int Convert_NDfield(NDfield *field,int type);

int Save_NDfield(NDfield *field,const char *filename);
int IsNDfield(const char *filename);
int IsNDfieldFromFITS(const char *filename);
NDfield *Load_NDfield(const char *filename);
NDfield *Load_NDfieldFromFITS(const char *filename);

NDfield *Load_NDfieldHeader(const char *filename);
NDfield *Load_NDfieldChunk(char *filename, double *x0, double *delta, int periodic);
NDfield *Load_NDfieldChunkHeader(char *filename, double *x0, double *delta, int periodic);
int Save_NDfieldPartial(char *filename, NDfield *header, NDfield *field, int periodic);
int NDIntersection(double *x0_a,double *delta_a,double *x0_b,double *delta_b, double *x0_bbox, double *delta_bbox,double *x0_cut,double *delta_cut, int ndims, int periodic);

void Coords2IndexfND(float *coord_p, INT *index,double *x0, double *delta,int *dims, int ndims,int periodic);
int InterpolateND(void *field_p,FLOAT *result,FLOAT *pos_p, int keepdim, int periodic);
  char *print_dataType(char *dest, int type);
  int printNDfieldStat(NDfield *field, int dec);

#ifdef __cplusplus
}
#endif

#endif
