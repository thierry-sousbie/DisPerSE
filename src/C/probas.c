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
#include "probas.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define INVERSION_PRECISION 1.E-6
#define SQRT2_INV 0.707106781186547461715008466853760182857513427734375

double inverseProba(double(*f)(const double, int), double val)
{
    if (val>=1) return 1;
    if (val<=0) return -1;

    double min=0.;
    double max=1.;
    
    do {max*=2;} while (f(max,0)>val);

    double cur;
    double fcur;
    do {
	cur = (max+min)/2;
	fcur = f(cur,0);
	if (fcur>=val) min = cur;
	else max = cur;
	//printf ("%e : %e = %e\n",max-min,f(cur,0),val);
    } while (fabs((fcur-val)/fcur)>INVERSION_PRECISION);

    return (max+min)/2;
}

double sigmaProba(const double s, int inverse)
{
  if (inverse) return inverseProba(sigmaProba,s);

  return erfc(s*SQRT2_INV);
}

double ratioProba2D_1(const double r, int inverse)
{
  if (inverse) return inverseProba(ratioProba2D_1,r);
  
  double a[3] = {2,0.01,3.5};
  double u=r-1;

  return exp(-a[0]*u-a[1]*pow(u,a[2]));
}

double ratioProba2D_2(const double r, int inverse)
{
  if (inverse) return inverseProba(ratioProba2D_2,r);

  double a[2] = {0.75,0.2};
  double u=r;
  
  return pow(u,-a[0]*(a[1]*log(u)+1));
}

double ratioProba3D_1(const double r, int inverse)
{
  if (inverse) return inverseProba(ratioProba3D_1,r);
  
  double a[3] = {3.69381,0.44052,2.53805};
  double u=r-1;

  return exp(-a[0]*u-a[1]*pow(u,a[2]));
}


double ratioProba3D_2(const double r, int inverse)
{
  if (inverse) return inverseProba(ratioProba3D_2,r);

  double a[4] = {2.55393,4.,9.,1.785};
  double u=r-1;
  
  double f1= exp(-a[0]*u);
  double f2= a[1]*pow(r,-a[2]);
  double t;

  if (u==0) t=0;
  else t=1./(1.+pow(a[3]/u,14));
  
  return f2*t + (1.-t)*f1;
}

double ratioProba3D_3(const double r, int inverse)
{
  //if (inverse) return inverseProba(ratioProba3D_3,r);
  
  double a[2] = {0.448993,2.5633};
  double u=r-1;

  if (inverse) return (pow(r,-1./a[1])-1)/a[0] +1;
  else return pow(1.+a[0]*u,-a[1]);
}

double probaFromSigma(double N)
{
  return erfc(N*SQRT2_INV);
}

double sigmaFromPersistenceRatio(double r, int ndims, int type)
{
  if (ndims==3) {
    if (type==0) return sigmaProba(ratioProba3D_1(r,0),1);
    if (type==1) return sigmaProba(ratioProba3D_2(r,0),1);
    if (type==2) return sigmaProba(ratioProba3D_3(r,0),1);

    fprintf (stderr,"sigmaFromPersistenceRatio: no type %d pairs in %dD.\n",type,ndims);
    exit(0);
  }
  else if (ndims==2) {
    //printf("res(%g) = %g\n",r,ratioProba2D_2(r,0));
    if (type==0) return sigmaProba(ratioProba2D_1(r,0),1);
    if (type==1) return sigmaProba(ratioProba2D_2(r,0),1);
    
    fprintf (stderr,"sigmaFromPersistenceRatio: no type %d pairs in %dD.\n",type,ndims);
    exit(0);
  }
  else
    {
      return r;
      //fprintf (stderr,": not implemented in %d dims.\n",ndims);
      //exit(0);
    }
}

int persistenceRatioFromSigma(double N, int ndims, double *r)
{
    int i;

    if (ndims==3) {
	double p = probaFromSigma(N);
	//printf ("\n p(%e sig)=%e\n",N,p);
	r[0] = ratioProba3D_1(p,1);
	r[1] = ratioProba3D_2(p,1);
	r[2] = ratioProba3D_3(p,1);
    }
    else if (ndims==2) {
	double p = probaFromSigma(N);
	//printf ("\n p(%e sig)=%e\n",N,p);
	r[0] = ratioProba2D_1(p,1);
	r[1] = ratioProba2D_2(p,1);
    }
    else
      {
	fprintf (stderr,"persistenceRatioFromSigma: not implemented in %d dims.\n",ndims);
	exit(0);
      }

    
    for (i=0;i<ndims;i++)
	if (r[i]<1) return -1;

    return 0;
}
