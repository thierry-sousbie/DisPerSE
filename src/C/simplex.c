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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "simplex.h"

#ifdef HAVE_GSL

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

double SimplexVolumed(double **coord,int nvert, int ndims)
{
    int i,j,k;
    double det;
    int s;
    double m_data[(nvert-1)*(nvert-1)];
    double w_data[(nvert-1)*(nvert-1)];
/*
    double w_data[3*3] = {1.,0.,0.,
			  1.,1.,0.,
			  1.,1.,1.};*/
    //double wt_data[(nvert-1)*ndims];
    gsl_permutation * perm;
    
    gsl_matrix_view m = gsl_matrix_view_array (m_data, nvert-1, nvert-1);    
    gsl_matrix_view w = gsl_matrix_view_array (w_data, nvert-1, ndims);
    //gsl_matrix_view wt = gsl_matrix_view_array (wt_data, nvert-1,ndims);

    if (nvert==2)
    {
	for (i=0,det=0;i<ndims;i++)
	    det+=(coord[1][i]-coord[0][i])*(coord[1][i]-coord[0][i]);
	return sqrt(det);
    }

    perm = gsl_permutation_alloc (nvert-1);

    for (i=0,k=0;i<nvert-1;i++)
	for (j=0;j<nvert-1;j++,k++)
	    w_data[k]=coord[i+1][j]-coord[0][j];
    
    //gsl_matrix_transpose_memcpy (wt,w);
    
    gsl_blas_dgemm (CblasNoTrans, CblasTrans,
		    1.0, &w.matrix, &w.matrix,
		    0.0, &m.matrix);
    
    gsl_linalg_LU_decomp(&m.matrix,perm,&s);
    det=gsl_linalg_LU_det(&m.matrix,s);
    
    gsl_permutation_free(perm);
    
    for (i=2,j=1;i<=nvert-1;i++) j*=i;
    
    return sqrt(fabs(det))/j;
}


double SimplexVolumef(float **coord,int nvert, int ndims)
{
    int i,j,k;
    double det;
    int s;
    double m_data[(nvert-1)*(nvert-1)];
    double w_data[(nvert-1)*(nvert-1)];
/*
    double w_data[3*3] = {1.,0.,0.,
			  1.,1.,0.,
			  1.,1.,1.};*/
    //double wt_data[(nvert-1)*ndims];
    gsl_permutation * perm;
    
    gsl_matrix_view m = gsl_matrix_view_array (m_data, nvert-1, nvert-1);    
    gsl_matrix_view w = gsl_matrix_view_array (w_data, nvert-1, ndims);
    //gsl_matrix_view wt = gsl_matrix_view_array (wt_data, nvert-1,ndims);

    if (nvert==2)
    {
	for (i=0,det=0;i<ndims;i++)
	    det+=(double)(coord[1][i]-coord[0][i])*(double)(coord[1][i]-coord[0][i]);
	return sqrt(det);
    }

    perm = gsl_permutation_alloc (nvert-1);

    for (i=0,k=0;i<nvert-1;i++)
	for (j=0;j<nvert-1;j++,k++)
	    w_data[k]=coord[i+1][j]-coord[0][j];
    
    
    //gsl_matrix_transpose_memcpy (wt,w);
    
    gsl_blas_dgemm (CblasNoTrans, CblasTrans,
		    1.0, &w.matrix, &w.matrix,
		    0.0, &m.matrix);
    
    gsl_linalg_LU_decomp(&m.matrix,perm,&s);
    det=gsl_linalg_LU_det(&m.matrix,s);
    
    gsl_permutation_free(perm);
    
    for (i=2,j=1;i<=nvert-1;i++) j*=i;
    
    return sqrt(fabs(det))/j;
}

// center must be allocated to ndims vals at least.
// coord are the vertice coords
// center is the center of the resulting ndims-sphere
// function returns its radius squared
double SimplexSphere(double **coord, int ndims,double *center)
{
    double a_data[ndims*ndims];
    double b_data[ndims];
    
    gsl_matrix_view m=gsl_matrix_view_array (a_data, ndims, ndims);
    gsl_vector_view b=gsl_vector_view_array (b_data, ndims);
    gsl_vector_view x = gsl_vector_view_array (center,ndims);
    gsl_permutation * p = gsl_permutation_alloc (ndims);

    int s;
    int i,j,k;
    double d;

    for (i=0,k=0;i<ndims;i++)
    {
	b_data[i]=0;
	for (j=0;j<ndims;j++,k++)
	{
	    a_data[k]=((double)coord[0][j]-(double)coord[i+1][j]);
	    b_data[i]+=((double)coord[0][j]*(double)coord[0][j]-(double)coord[i+1][j]*(double)coord[i+1][j]);
	}
	b_data[i]*=0.5;
    }

    gsl_linalg_LU_decomp (&m.matrix, p, &s);
    gsl_linalg_LU_solve (&m.matrix, p, &b.vector, &x.vector);      
    gsl_permutation_free (p);

    for (i=0,d=0;i<ndims;i++) 
	d+=((double)center[i]-(double)coord[0][i])*((double)center[i]-(double)coord[0][i]);

    return d;
}

#else

double SimplexVolumef(float **coord,int nvert, int ndims)
{
    int i;
    double det;
  
    if (nvert==2)
    {
	for (i=0,det=0;i<ndims;i++)
	    det+=(double)(coord[1][i]-coord[0][i])*(double)(coord[1][i]-coord[0][i]);
	return sqrt(det);
    }
    

    fprintf(stderr,"ERROR: Function 'SimplexVolume' needs GSL.\n");
    fprintf(stderr,"Please, recompile with GSL.\n");
    return -1;
}

double SimplexSphere(double **coords, int ndims,double *center)
{
    fprintf(stderr,"ERROR: Function 'SimplexSphere' needs GSL.\n");
    fprintf(stderr,"Please, recompile with GSL.\n");

    return -1;
}


#endif
