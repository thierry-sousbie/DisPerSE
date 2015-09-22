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
#include "dtfi.h"

#ifdef HAVE_GSL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// mu is the average value
unsigned int randomPoisson(double mu)
{
  static int first_call = 1;
  static const gsl_rng_type *T;
  static gsl_rng * r;
     
  if (first_call) {
    first_call=0;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
  }
  
  return gsl_ran_poisson (r, mu);
}

void randomSimplexCoords(double **coord, int ndims, double *result, double *P)
{
  static int first_call = 1;
  static const gsl_rng_type *T;
  static gsl_rng * r;
  int i,j,k;
  double u;
  //double P[50];
  
  if (first_call) {
    first_call=0;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
  }

  result[0]=gsl_rng_uniform (r);
  
  if (ndims==3)
    {
      u=gsl_rng_uniform (r);
      if (u<result[0])
	{
	  result[1]=result[0];
	  result[0]=u;
	}
      else result[1]=u;

      u=gsl_rng_uniform (r);
      if (u<result[0])
	{
	  result[2]=result[1];
	  result[1]=result[0];
	  result[0]=u;
	}
      else if (u<result[1])
	{
	  result[2]=result[1];
	  result[1]=u;
	}
      else result[2]=u;

    }
  else
    {
      for (i=1;i<ndims;i++)
	{
	  u = gsl_rng_uniform (r);
	  
	  for (j=0;j<i;j++)
	    if (u<result[j])
	      {
		result[i]=result[j];
		result[j]=u;
		break;
	      } 
	  
	  if (j==i) result[i]=u;
	  else {
	    j++;
	    for (;j<i;j++)
	      for (k=j+1;k<=i;k++)
		{
		  if (result[j]>result[k])
		    {
		      double tmp = result[k];
		      result[k]=result[j];
		      result[j]=tmp;
		    }
		}
	  }
	}
    }
    
  P[0]=result[0];
  for (i=1;i<ndims;i++) P[i]= result[i]-result[i-1];
  P[ndims]=1.-result[ndims-1];

  for (i=0;i<ndims;i++) result[i]=P[0]*coord[0][i];
  
  for (i=1;i<ndims+1;i++)
    for (j=0;j<ndims;j++) 
      result[j]+=P[i]*coord[i][j];
    
}

void genDTFIInSimplex(double **coord, double *Proba, int ngen, int ndims, double **result_p)
{
  
  if (*result_p==NULL) {
    *result_p=(double*)malloc(sizeof(double)*ngen*ndims);
  }

  double *result = *result_p;
  double P[ndims];
  int i,n;
  double keepProba;
  double pmax_inv=0;

  static int first_call = 1;
  static const gsl_rng_type *T;
  static gsl_rng * r;

  //if (ngen>100) printf("negen = %d\n",ngen);

  if (first_call) {
    first_call=0;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
  }
  
  if (Proba!=NULL) {
    pmax_inv=0;
    for (i=0;i<ndims+1;i++) 
      if (Proba[i]>pmax_inv) pmax_inv=Proba[i];
    pmax_inv=(double)1./pmax_inv;
  }

  for (n=0;n<ngen;)
    {
      randomSimplexCoords(coord,ndims,&(result[n*ndims]),P);
      if (Proba!=NULL) {
	keepProba=Proba[0]*pmax_inv*P[0];
	for (i=1;i<ndims+1;i++) keepProba+=Proba[i]*pmax_inv*P[i];
            
	if (gsl_rng_uniform (r)<=keepProba) n++;
      } else n++;
    }
  
}
     
#else

unsigned int randomPoisson(double mu)
{
  fprintf(stderr,"ERROR: Function 'randomPoisson' needs GSL.\n");
  fprintf(stderr,"Please, recompile with GSL.\n");

  return -1;
}

void randomSimplexCoords(double **coord, int ndims, double *result, double *P)
{
  fprintf(stderr,"ERROR: Function 'randomSimplexCoords' needs GSL.\n");
  fprintf(stderr,"Please, recompile with GSL.\n");

  return -1;
}

void genDTFIInSimplex(double **coord, double *Proba, int ngen, int ndims, double **result_p)
{
  fprintf(stderr,"ERROR: Function 'genDTFIInSimplex' needs GSL.\n");
  fprintf(stderr,"Please, recompile with GSL.\n");
}

#endif
