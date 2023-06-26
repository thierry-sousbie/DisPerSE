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
// distances are in Mpc
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "distances.h"
#include "mystring.h"

double cosmoD_Z_STEP = 1.E-5; 
double cosmoD_local_h=0.72;
double cosmoD_local_Om=0.27;
double cosmoD_local_Ol=0.73;
double cosmoD_local_Ok=0.0;
double *cosmoD_local_w_val=NULL;
double *cosmoD_local_z_val=NULL;
char cosmoD_local_w_fname[255]="";
double cosmoD_local_Nw_val;
double cosmoD_local_distances_initialized=0;

#ifndef COSMOD_LIGHTSPEED
#define COSMOD_LIGHTSPEED  299792.458
#endif

static const double cosmoD_richardson[] = {  
  3.333333333333333333e-01, 6.666666666666666667e-02, 1.587301587301587302e-02,
  3.921568627450980392e-03, 9.775171065493646139e-04, 2.442002442002442002e-04,
  6.103888176768601599e-05, 1.525902189669642176e-05, 3.814711817595739730e-06,
  9.536752259018191355e-07, 2.384186359449949133e-07, 5.960464832810451556e-08,
  1.490116141589226448e-08, 3.725290312339701922e-09, 9.313225754828402544e-10,
  2.328306437080797376e-10, 5.820766091685553902e-11, 1.455191522857861004e-11,
  3.637978807104947841e-12, 9.094947017737554185e-13, 2.273736754432837583e-13,
  5.684341886081124604e-14, 1.421085471520220567e-14, 3.552713678800513551e-15,
  8.881784197001260212e-16, 2.220446049250313574e-16
};

#define COSMOD_MAX_COLUMNS 1+sizeof(cosmoD_richardson)/sizeof(cosmoD_richardson[0])

#define COSMOD_MAX(x,y) ( (x) < (y) ? (y) : (x) )
#define COSMOD_MIN(x,y) ( (x) < (y) ? (x) : (y) )

////////////////////////////////////////////////////////////////////////////////
//  double Rombergs_Integration_Method( double a, double h, double tolerance, //
//                             int max_cols, double (*f)(double), int *err ); //
//                                                                            //
//  Description:                                                              //
//    If T(f,h,a,b) is the result of applying the trapezoidal rule to approx- //
//    imating the integral of f(x) on [a,b] using subintervals of length h,   //
//    then if I(f,a,b) is the integral of f(x) on [a,b], then                 //
//                           I(f,a,b) = lim T(f,h,a,b)                        //
//    where the limit is taken as h approaches 0.                             //
//    The classical Romberg method applies Richardson Extrapolation to the    //
//    limit of the sequence T(f,h,a,b), T(f,h/2,a,b), T(f,h/4,a,b), ... ,     //
//    in which the limit is approached by successively deleting error terms   //
//    in the Euler-MacLaurin summation formula.                               //
//                                                                            //
//  Arguments:                                                                //
//     double a          The lower limit of the integration interval.         //
//     double h          The length of the interval of integration, h > 0.    //
//                       The upper limit of integration is a + h.             //
//     double tolerance  The acceptable error estimate of the integral.       //
//                       Iteration stops when the magnitude of the change of  //
//                       the extrapolated estimate falls below the tolerance. //
//     int    max_cols   The maximum number of columns to be used in the      //
//                       Romberg method.  This corresponds to a minimum       //
//                       integration subinterval of length 1/2^max_cols * h.  //
//     double *f         Pointer to the integrand, a function of a single     //
//                       variable of type double.                             //
//     int    *err       0 if the extrapolated error estimate falls below the //
//                       tolerance; -1 if the extrapolated error estimate is  //
//                       greater than the tolerance and the number of columns //
//                       is max_cols.                                         //
//                                                                            //
//  Return Values:                                                            //
//     The integral of f(x) from a to a +  h.                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
double cosmoD_Rombergs_Integration_Method( double a, double h, double tolerance,
                               int max_cols, double (*f)(double), int *err ) { 
   
   double upper_limit = a + h;     // upper limit of integration
   double dt[COSMOD_MAX_COLUMNS];         // dt[i] is the last element in column i.
   double integral = 0.5 * ( (*f)(a) + (*f)(a+h) );
   double x, old_h, delta;
   int j,k;

// Initialize err and the first column, dt[0], to the numerical estimate //
// of the integral using the trapezoidal rule with a step size of h.     //
 
   *err = 0;
   dt[0] = 0.5 * h *  ( (*f)(a) + (*f)(a+h) );

// For each possible succeeding column, halve the step size, calculate  //
// the composite trapezoidal rule using the new step size, and up date  //
// preceeding columns using Richardson extrapolation.                   //

   max_cols = COSMOD_MIN(COSMOD_MAX(max_cols,0),COSMOD_MAX_COLUMNS);
   for (k = 1; k < max_cols; k++) {
      old_h = h;

                 // Calculate T(f,h/2,a,b) using T(f,h,a,b) //
 
      h *= 0.5;
      integral = 0.0;
      for (x = a + h; x < upper_limit; x += old_h) integral +=  (*f)(x);
      integral = h * integral + 0.5 * dt[0];

         //  Calculate the Richardson Extrapolation to the limit //

      for (j = 0; j < k; j++) {
         delta =  integral - dt[j];
         dt[j] = integral;
         integral += cosmoD_richardson[j] * delta;
      } 

      //  If the magnitude of the change in the extrapolated estimate //
      //  for the integral is less than the preassigned tolerance,    //
      //  return the estimate with err = 0.                           //

      if ( fabs( delta ) < tolerance ) {
         return integral;
      }
      
             //  Store the current esimate in the kth column. //

      dt[k] = integral;
   }

     // The process didn't converge within the preassigned tolerance //
     // using the maximum number of columns designated.              //
     // Return the current estimate of integral and set err = -1.    //
   
   *err = -1;
   return integral;
}

//find the index of the value in tab closest to val
int cosmoD_findIndex(double *tab,int n,double val)
{
  int imin,imax;

  if (tab[0]>tab[n-1])
    {
      imin=n-1;
      imax=0;
    }
  else
    {
      imin=0;
      imax=n-1;
    }

  while (abs(imax-imin)>1)
    {
      if (tab[(imin+imax)/2]>val)
	imax=(imin+imax)/2;
      else
	imin=(imin+imax)/2;
    }

  if (fabs(tab[imin]-val)<fabs(tab[imax]-val)) 
      return imin;
  else
      return imax;
}



/*inline*/ double cosmoD_E_z(double z,double Om,double Ol,double Ok,double w)
{
  return 
    sqrt(Om*pow(1+z,3)+Ok*pow(1+z,2.)+(Ol)*pow(1+z,3*(1+w)));

}

double cosmoD_f_z(double z)
{
  double w;

  if (cosmoD_local_Nw_val==1)
    w=*cosmoD_local_w_val;
  else
    w=cosmoD_local_w_val[cosmoD_findIndex(cosmoD_local_w_val,cosmoD_local_Nw_val,z)]; 


  return 1./cosmoD_E_z(z,cosmoD_local_Om,cosmoD_local_Ol,cosmoD_local_Ok,w);    
}

//d here is the comoving distance , not luminosity distance
//dl=d*(1+z)
double cosmoD_z2d_local(double z,double h,double Om_p,double Ol_p,double Ok_p,double *w_p,double *z_p,int Nw_val)
{
    static double Om=-1;
    static double Ol=-1;
    static double Ok=-1;
    static double *wptr=NULL;
    static double z_max=0;
    static double *d=NULL;

    double a;
    int i;
    double dist;
    double w;
    int err;
    
    //recompute everything if cosmo has changed
    if ((Om!=Om_p)||(Ol!=Ol_p)||(Ok!=Ok_p)||(wptr!=w_p)) {z_max=0;Om=Om_p;Ok=Ok_p;Ol=Ol_p;wptr=w_p;}
    //compute up to max redshift, and at least by steps of dz=500*Z_STEP
    if (z>z_max)
    {
	dist = 0.5+((int)(z/(0.5)))*(0.5);

	d=(double *) realloc(d,sizeof(double)*(int)(3+dist/cosmoD_Z_STEP));
	if (z_max==0) d[0]=0;

	for (a=z_max+cosmoD_Z_STEP,i=(z_max/cosmoD_Z_STEP)+1;a<=dist+2*cosmoD_Z_STEP;a+=cosmoD_Z_STEP,i++)
	{
	  /*
	  if (Nw_val==1)
	    w=*w_p;
	  else
	    w=w_p[FindIndex(z_p,Nw_val,a)]; 
	  */

	  d[i]=cosmoD_Rombergs_Integration_Method( 0., a, 1.E-7,30,cosmoD_f_z,&err);

	  //d[i]=d[i-1]+Z_STEP/E_z(a+Z_STEP/2,Om,Ol,Ok,w);
	  //d[i]=d[i-1]+1./sqrt(Om*pow(1+a+Z_STEP/2,3)+Ok*(1+a+Z_STEP/2,2.)+(Ol)*pow(1+a+Z_STEP/2,3*(1+w)))*LIGHTSPEED*Z_STEP;
	}
	z_max=cosmoD_Z_STEP*(i-1);
    }
    
    if (z>0.) 
    {
	dist = 1-(z-cosmoD_Z_STEP*(int)(z/cosmoD_Z_STEP))/cosmoD_Z_STEP;
	dist = COSMOD_LIGHTSPEED/(100.*h) * (dist*d[(int)(z/cosmoD_Z_STEP)] + (1-dist)*d[(int)(z/cosmoD_Z_STEP)+1]);
	return dist;
    }
    

    return 0.;
}

//Inverts z2d_local function
//d here is the comoving distance , not luminosity distance
double cosmoD_d2z_local(double d,double h,double Om,double Ol, double Ok,double *w_val,double *z_val,int Nw_val)
{
    double zmin=0,zmax=1000*cosmoD_Z_STEP;
    double dmax,dmin;
    double z;
    
    while (cosmoD_z2d_local(zmax,h,Om,Ol,Ok,w_val,z_val,Nw_val)<d) {zmin=zmax;zmax+=1000*cosmoD_Z_STEP;}

    while(zmax-zmin>cosmoD_Z_STEP)
    {
	z=(zmax+zmin)/2;
	if (cosmoD_z2d_local(z,h,Om,Ol,Ok,w_val,z_val,Nw_val)<=d) zmin=z;
	else zmax=z;
    }

    dmin = cosmoD_z2d_local(zmin,h,Om,Ol,Ok,w_val,z_val,Nw_val);
    dmax = cosmoD_z2d_local(zmax,h,Om,Ol,Ok,w_val,z_val,Nw_val);

    return zmin+(zmax-zmin)*(d-dmin)/(dmax-dmin);
}

/*inline*/ double cosmoD_d2dm(double d)
{
  double Ok;
  double DH=COSMOD_LIGHTSPEED/(100.*cosmoD_local_h);

  Ok=cosmoD_local_Ok;
 
  if (Ok>1.E-5)
    return DH/sqrt(Ok)*sinh(sqrt(Ok)*d/DH);
  
  if (Ok<-1.E-5)
    return DH/sqrt(-Ok)*sin(sqrt(-Ok)*d/DH);

  return d;	
}

/*inline*/ double cosmoD_dm2d(double d)
{
  double Ok;
  double DH=COSMOD_LIGHTSPEED/(100.*cosmoD_local_h);

  Ok=cosmoD_local_Ok;

  if (Ok>1.E-5)
    return DH/sqrt(Ok)*asinh(sqrt(Ok)*d/DH);

  if (Ok<-1.E-5)
    return DH/sqrt(-Ok)*asin(sqrt(-Ok)*d/DH);

  return d;
}



double cosmoD_z2d(double z)
{
    if (!cosmoD_local_distances_initialized)
    {
	fprintf (stderr,"ERROR: Distances not initialized !!! (forgot to call InitCosmoDistances())");
	return 0;
    }
    return cosmoD_z2d_local(z,cosmoD_local_h,cosmoD_local_Om,cosmoD_local_Ol,cosmoD_local_Ok,cosmoD_local_w_val,cosmoD_local_z_val,cosmoD_local_Nw_val);
}

double cosmoD_Hz(double z)
{
  return COSMOD_LIGHTSPEED*z/cosmoD_z2d(z);
}

/*inline*/ double cosmoD_z2dm(double z)
{
  return cosmoD_d2dm(cosmoD_z2d(z));
}

double cosmoD_d2z(double z)
{
    if (!cosmoD_local_distances_initialized)
    {
	fprintf (stderr,"ERROR: Distances not initialized !!! (forgot to call InitCosmoDistances())");
	return 0;
    }
    return cosmoD_d2z_local(z,cosmoD_local_h,cosmoD_local_Om,cosmoD_local_Ol,cosmoD_local_Ok,cosmoD_local_w_val,cosmoD_local_z_val,cosmoD_local_Nw_val);
}

/*
double v2z(double v)
{
    return d2z(v/(100.*cosmoD_local_h));
}

double z2v(double z)
{  
    return cosmoD_local_h*100.*z2d(z);
}

double v2d(double v)
{
    return v/(100.*cosmoD_local_h);
}

double d2v(double d)
{
    return cosmoD_local_h*100.*d;
}
*/

/*inline*/ double cosmoD_v2z(double v)
{
  return sqrt((1.+v/COSMOD_LIGHTSPEED)/(1.-v/COSMOD_LIGHTSPEED))-1.;
}

/*inline*/ double cosmoD_z2v(double z)
{
  return (COSMOD_LIGHTSPEED*z*(2+z))/(2+2*z+z*z);
}

/*inline*/ double cosmoD_v2z_approx(double v)
{
  return v/COSMOD_LIGHTSPEED;
}

/*inline*/ double cosmoD_z2v_approx(double z)
{
  return COSMOD_LIGHTSPEED*z;
}

/*inline*/ double cosmoD_v2d(double v)
{
  return cosmoD_z2d(cosmoD_v2z(v));
}

/*inline*/ double cosmoD_d2v(double d)
{
  return cosmoD_z2v(cosmoD_d2z(d));
}

//frees everything ...
int cosmoD_reset()
{
    cosmoD_local_distances_initialized=0;
    free(cosmoD_local_w_val);
    free(cosmoD_local_z_val);
    cosmoD_local_w_val=NULL;
    cosmoD_local_z_val=NULL;

    return 0;
}

int cosmoD_initialized()
{
  return (cosmoD_local_distances_initialized!=0);
}

int cosmoD_getParamsStr(char* str)
{
  
  if (strcmp(cosmoD_local_w_fname,"")==0) {
    sprintf(str,"Om=%.2f Ol=%.2f Ok=%.2f h=%.2f w=%.2f",cosmoD_local_Om,cosmoD_local_Ol,cosmoD_local_Ok,cosmoD_local_h,cosmoD_local_w_val[0]);
  }
  else
    {
      if (strlen(cosmoD_local_w_fname)<40)
	sprintf(str,"Om=%.2f Ol=%.2f Ok=%.2f h=%.2f w=<'%s'",cosmoD_local_Om,cosmoD_local_Ol,cosmoD_local_Ok,cosmoD_local_h,cosmoD_local_w_fname);
      else
	{
	  sprintf(str,"Om=%.2f Ol=%.2f Ok=%.2f h=%.2f w=<'%s'",cosmoD_local_Om,cosmoD_local_Ol,cosmoD_local_Ok,cosmoD_local_h,&cosmoD_local_w_fname[strlen(cosmoD_local_w_fname)-39]);
	}
    }

  return 0;
}

//w is ignored if filename is not NULL
//h is reduced hubble (~0.7)
int cosmoD_init(double Om,double Ol,double Ok,double h,double w,char *Filename)
{
  cosmoD_local_Om=Om;
  cosmoD_local_Ol=Ol;
  cosmoD_local_Ok=Ok;
  cosmoD_local_h=h;
  FILE *f;
  int i;
  
  cosmoD_local_distances_initialized=1;

  if (Filename!=NULL)
    strcpy(cosmoD_local_w_fname,Filename);
  else
    strcpy(cosmoD_local_w_fname,"");
  
  //Now set value for w
  if ((Filename==NULL)||(strcmp(Filename,"")==0))
    {
      cosmoD_local_w_val=(double *) malloc(2*sizeof(double));
      cosmoD_local_z_val=(double *) malloc(2*sizeof(double));
      
      cosmoD_local_w_val[0]=cosmoD_local_w_val[1]=w;
      cosmoD_local_z_val[0]=1000;
      cosmoD_local_z_val[1]=0;
      cosmoD_local_Nw_val=2;
    }
  //For quintessence
  else
    {
      char buf[200];
      char buf1[200]; 
      char buf2[200];
      
      cosmoD_local_z_val = (double *) malloc (100000*sizeof(double));
      cosmoD_local_w_val = (double *) malloc (100000*sizeof(double));
      
      if((f = fopen(Filename, "r"))!=NULL)
	{
	  i=0;
	  char *line=NULL;
	  int linesize=0;
	  while(!feof(f))
	    {
	      Mygetline(&line,&linesize,f);

	      if (line[0]=='#') continue;
	      if (strlen(line)<=1) continue;

	      sscanf(line, "%s%s", buf1, buf2);
	      cosmoD_local_z_val[i]=atof(buf1);
	      cosmoD_local_w_val[i++]=atof(buf2);
	      //printf("%f %f\n",cosmoD_local_z_val[i-1],cosmoD_local_w_val[i-1]);
	    }
        }
      else 
	{
	  fprintf (stderr,"File %s does not exist ...\n",Filename);
	  return -1;
	}
      fclose(f);

      if (cosmoD_local_z_val[i-1]!=0) 
	{
	  cosmoD_local_z_val[i]=0;
	  cosmoD_local_w_val[i]=cosmoD_local_w_val[i-1];
	  cosmoD_local_Nw_val = i+1;
	}
      else cosmoD_local_Nw_val = i;

      cosmoD_local_w_val = (double *) realloc(cosmoD_local_w_val,cosmoD_local_Nw_val*sizeof(double));
      cosmoD_local_z_val = (double *) realloc(cosmoD_local_z_val,cosmoD_local_Nw_val*sizeof(double));
    }
  char str[1024];
  cosmoD_getParamsStr(str);
  printf("\nInitialized cosmo distances: %s.\n",str);
  return 0;
}
