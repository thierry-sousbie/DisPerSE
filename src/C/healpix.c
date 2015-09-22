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
/* This file uses functions copy/pasted from chealpix but is not part of Healpix itelf ... */

/* -----------------------------------------------------------------------------
 *
 *  Copyright (C) 1997-2005 Krzysztof M. Gorski, Eric Hivon, 
 *                          Benjamin D. Wandelt, Anthony J. Banday, 
 *                          Matthias Bartelmann, 
 *                          Reza Ansari & Kenneth M. Ganga 
 *
 *
 *  This file is part of HEALPix.
 *
 *  HEALPix is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  HEALPix is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with HEALPix; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix see http://healpix.jpl.nasa.gov
 *
 *----------------------------------------------------------------------------- */

/* Standard Includes */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <math.h>

/* Local Includes */
#include "healpix.h"

static int order_max=13;
static int jrll[12] = {2,2,2,2,3,3,3,3,4,4,4,4};
static int jpll[12] = {1,3,5,7,0,2,4,6,1,3,5,7};
static short *utab=NULL;
static short *ctab=NULL;


#ifdef HAVE_CFITS_IO

#include "fitsio.h"

void printerror( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/

    if (status)
    {
       fits_report_error(stderr, status); /* print error report */

       exit( status );    /* terminate the program, returning error status */
    }
    return;
}

int is_healpix_map(const char *infile)
{
  fitsfile *fptr=NULL;
  char ordering[10];
  int status=0;
  int hdutype, nfound, anynul;
  char     comment[FLEN_COMMENT];

  if ( fits_open_file(&fptr, infile, READONLY, &status) ) {
    return 0;
  }

  if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) ) {
    fits_close_file(fptr, &status);
    return 0;
  }
 
  if (fits_read_key(fptr, TSTRING, "ORDERING", ordering, comment, &status)) {
    fits_close_file(fptr, &status);
    return 0;
  }
  else {
    if ((strstr(ordering,"NESTED") == NULL)&&(strstr(ordering,"RING") == NULL))
      return 0;
  }

  if ( fits_close_file(fptr, &status) ) {
    return 0;
  }

  return 1;
}

float *read_healpix_map_f(const char *infile, long *nside, char *coordsys, char *ordering) {
  
  /* Local Declarations */
  long     naxes, *naxis, npix, npercol, irow;
  int      status, hdutype, nfound, anynul;
  float    nulval, *map;
  char     comment[FLEN_COMMENT];
  fitsfile *fptr;

  /* Initializations */
  status = 0;

  /* Open the file */
  if ( fits_open_file(&fptr, infile, READONLY, &status) ) {
    printerror( status );
    return (float *)NULL;
  }

  /* Move to the HDU */
  if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) ) {
    printerror( status );
    return (float *)NULL;
  }
  if (hdutype != BINARY_TBL) {
    fprintf(stderr, "%s (%d): Extension is not binary!\n", __FILE__, __LINE__);
    return (float *)NULL;
  }

  /* Read the sizes of the array */
  if ( fits_read_key_lng(fptr, "NAXIS", &naxes, comment, &status) ) {
    printerror( status );
    return (float *)NULL;
  }
  naxis = (long *)malloc(((size_t)naxes)*sizeof(long));
  if ( fits_read_keys_lng(fptr, "NAXIS", 1, naxes, naxis, &nfound, &status) 
       || nfound != naxes ) {
    printerror( status );
    return (float *)NULL;
  }

  if ( fits_read_key_lng(fptr, "NSIDE", nside, comment, &status) ) {
    printerror(status);
    return (float *)NULL;
  }
  
  npix = 12*(*nside)*(*nside);
  if ( (npix%naxis[1]) != 0 ) {
    fprintf(stderr, "%s (%d): Problem with npix.\n", __FILE__, __LINE__);
    return (float *)NULL;
  }
  npercol = npix/naxis[1];
  /*
  if (fits_read_key(fptr, TSTRING, "COORDSYS",coordsys, comment, &status)) {
    //fprintf(stderr, "WARNING: Could not find %s keyword in in file %s\n", 
    //	    "COORDSYS",infile);
    status = 0;
  }
  */
  if (fits_read_key(fptr, TSTRING, "ORDERING", ordering, comment, &status)) {
    fprintf(stderr, "WARNING: Could not find %s keyword in in file %s\n", 
	    "ORDERING",infile);
    status = 0;
  }

  /* Read the array */
  map = (float *)malloc(((size_t)npix)*sizeof(float));
  nulval = HEALPIX_NULLVAL;
  for (irow = 0; irow < naxis[1]; irow++) {
    if ( fits_read_col(fptr, TFLOAT, 1, irow+1, 1, npercol, &nulval, 
		       &(map[irow*npercol]), &anynul, &status) ) {
      printerror(status);
      return (float *)NULL;
    }
  }

  /* Close the file */
  if ( fits_close_file(fptr, &status) ) {
    printerror( status );
    return (float *)NULL;
  }

  /* Later */
  return map;
}


double *read_healpix_map_d(const char *infile, long *nside, char *coordsys, char *ordering) {
  
  /* Local Declarations */
  long     naxes, *naxis, npix, npercol, irow;
  int      status, hdutype, nfound, anynul;
  double    nulval, *map;
  char     comment[FLEN_COMMENT];
  fitsfile *fptr;

  /* Initializations */
  status = 0;

  /* Open the file */
  if ( fits_open_file(&fptr, infile, READONLY, &status) ) {
    printerror( status );
    return (double *)NULL;
  }

  /* Move to the HDU */
  if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) ) {
    printerror( status );
    return (double *)NULL;
  }
  if (hdutype != BINARY_TBL) {
    fprintf(stderr, "%s (%d): Extension is not binary!\n", __FILE__, __LINE__);
    return (double *)NULL;
  }

  /* Read the sizes of the array */
  if ( fits_read_key_lng(fptr, "NAXIS", &naxes, comment, &status) ) {
    printerror( status );
    return (double *)NULL;
  }
  naxis = (long *)malloc(((size_t)naxes)*sizeof(long));
  if ( fits_read_keys_lng(fptr, "NAXIS", 1, naxes, naxis, &nfound, &status) 
       || nfound != naxes ) {
    printerror( status );
    return (double *)NULL;
  }

  if ( fits_read_key_lng(fptr, "NSIDE", nside, comment, &status) ) {
    printerror(status);
    return (double *)NULL;
  }
  
  npix = 12*(*nside)*(*nside);
  if ( (npix%naxis[1]) != 0 ) {
    fprintf(stderr, "%s (%d): Problem with npix.\n", __FILE__, __LINE__);
    return (double *)NULL;
  }
  npercol = npix/naxis[1];
  /*
  if (fits_read_key(fptr, TSTRING, "COORDSYS",coordsys, comment, &status)) {
    //fprintf(stderr, "WARNING: Could not find %s keyword in in file %s\n", 
    //	    "COORDSYS",infile);
    status = 0;
  }
  */
  if (fits_read_key(fptr, TSTRING, "ORDERING", ordering, comment, &status)) {
    fprintf(stderr, "WARNING: Could not find %s keyword in in file %s\n", 
	    "ORDERING",infile);
    status = 0;
  }

  /* Read the array */
  map = (double *)malloc(((size_t)npix)*sizeof(double));
  nulval = HEALPIX_NULLVAL;
  for (irow = 0; irow < naxis[1]; irow++) {
    if ( fits_read_col(fptr, TDOUBLE, 1, irow+1, 1, npercol, &nulval, 
		       &(map[irow*npercol]), &anynul, &status) ) {
      printerror(status);
      return (double *)NULL;
    }
  }

  /* Close the file */
  if ( fits_close_file(fptr, &status) ) {
    printerror( status );
    return (double *)NULL;
  }

  /* Later */
  return map;
}

#else

void printerror( int status) {printf ("Healpix error: libCFITSIO no found.\n");}
int is_healpix_map(const char *infile) {return 0;}
float *read_healpix_map_f(const char *infile, long *nside, char *coordsys, char *ordering) {printf ("Healpix error: libCFITSIO not found.\n");return NULL;}
double *read_healpix_map_d(const char *infile, long *nside, char *coordsys, char *ordering) {printf ("Healpix error: libCFITSIO not found.\n");return NULL;}

#endif


long nside2npix(const long nside) { return 12*nside*nside; }
long npix2nside(const long npix) {return (long)floor(sqrt(npix/12.)+0.5);}
static int nside2order(const long nside)
{
  int m;
  for (m=0; m<=order_max; ++m)
    {
      int nstest = 1<<m;
      if (nside == nstest) return m;
      if (nside < nstest) return -1;
    }
  return -1;
}

void mk_pix2xy(int *pix2x, int *pix2y) {

  /* =======================================================================
   * subroutine mk_pix2xy
   * =======================================================================
   * constructs the array giving x and y in the face from pixel number
   * for the nested (quad-cube like) ordering of pixels
   *
   * the bits corresponding to x and y are interleaved in the pixel number
   * one breaks up the pixel number by even and odd bits
   * =======================================================================
   */

  int i, kpix, jpix, IX, IY, IP, ID;
  for (i = 0; i < 1023; i++) pix2x[i]=0;
  
  for( kpix=0;kpix<1024;kpix++ ) {
    jpix = kpix;
    IX = 0;
    IY = 0;
    IP = 1 ;//              ! bit position (in x and y)
    while( jpix!=0 ){// ! go through all the bits
      ID = (int)fmod(jpix,2);//  ! bit value (in kpix), goes in ix
      jpix = jpix/2;
      IX = ID*IP+IX;
      
      ID = (int)fmod(jpix,2);//  ! bit value (in kpix), goes in iy
      jpix = jpix/2;
      IY = ID*IP+IY;
      
      IP = 2*IP;//         ! next bit (in x and y)
    }
    
    pix2x[kpix] = IX;//     ! in 0,31
    pix2y[kpix] = IY;//     ! in 0,31
  }
  
  /* Later */
  return;
}

void mk_xy2pix(int *x2pix, int *y2pix) {
  /* =======================================================================
   * subroutine mk_xy2pix
   * =======================================================================
   * sets the array giving the number of the pixel lying in (x,y)
   * x and y are in {1,128}
   * the pixel number is in {0,128**2-1}
   *
   * if  i-1 = sum_p=0  b_p * 2^p
   * then ix = sum_p=0  b_p * 4^p
   * iy = 2*ix
   * ix + iy in {0, 128**2 -1}
   * =======================================================================
   */
  int i, K,IP,I,J,ID;
  
  for(i = 0; i < 127; i++) x2pix[i] = 0;
  for( I=1;I<=128;I++ ) {
    J  = I-1;//            !pixel numbers
    K  = 0;//
    IP = 1;//
    truc : if( J==0 ) {
      x2pix[I-1] = K;
      y2pix[I-1] = 2*K;
    }
    else {
      ID = (int)fmod(J,2);
      J  = J/2;
      K  = IP*ID+K;
      IP = IP*4;
      goto truc;
    }
  }     
  
}

void pix2ang_ring( long nside, long ipix, double *theta, double *phi) {
  /*
    c=======================================================================
    c     gives theta and phi corresponding to pixel ipix (RING) 
    c     for a parameter nside
    c=======================================================================
  */
  
  int nl2, nl4, npix, ncap, iring, iphi, ip, ipix1;
  double  fact1, fact2, fodd, hip, fihip;
  double PI=M_PI;
  //      PARAMETER (pi     = 3.1415926535897932384626434d0)
  //      parameter (ns_max = 8192) ! 2^13 : largest nside available
  
  int ns_max=8192;
  
  if( nside<1 || nside>ns_max ) {
    fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
    exit(0);
  }
  npix = 12*nside*nside;      // ! total number of points
  if( ipix<0 || ipix>npix-1 ) {
    fprintf(stderr, "%s (%d): ipix out of range: %ld\n", __FILE__, __LINE__, ipix);
    exit(0);
  }
  
  ipix1 = ipix + 1; // in {1, npix}
  nl2 = 2*nside;
  nl4 = 4*nside;
  ncap = 2*nside*(nside-1);// ! points in each polar cap, =0 for nside =1
  fact1 = 1.5*nside;
  fact2 = 3.0*nside*nside;
  
  if( ipix1 <= ncap ) {  //! North Polar cap -------------
    
    hip   = ipix1/2.;
    fihip = floor(hip);
    iring = (int)floor( sqrt( hip - sqrt(fihip) ) ) + 1;// ! counted from North pole
    iphi  = ipix1 - 2*iring*(iring - 1);
    
    *theta = acos( 1. - iring*iring / fact2 );
    *phi   = (1.*iphi - 0.5) * PI/(2.*iring);
  }
  else if( ipix1 <= nl2*(5*nside+1) ) {//then ! Equatorial region ------
    
    ip    = ipix1 - ncap - 1;
    iring = (int)floor( ip / nl4 ) + nside;// ! counted from North pole
    iphi  = (int)fmod(ip,nl4) + 1;
    
    fodd  = 0.5 * (1 + fmod((double)(iring+nside),2));//  ! 1 if iring+nside is odd, 1/2 otherwise
    *theta = acos( (nl2 - iring) / fact1 );
    *phi   = (1.*iphi - fodd) * PI /(2.*nside);
  }
  else {//! South Polar cap -----------------------------------
    
    ip    = npix - ipix1 + 1;
    hip   = ip/2.;
/* bug corrige floor instead of 1.* */
    fihip = floor(hip);
    iring = (int)floor( sqrt( hip - sqrt(fihip) ) ) + 1;//     ! counted from South pole
    iphi  = (int)(4.*iring + 1 - (ip - 2.*iring*(iring-1)));
    
    *theta = acos( -1. + iring*iring / fact2 );
    *phi   = (1.*iphi - 0.5) * PI/(2.*iring);
  }
}


void pix2ang_nest( long nside, long ipix, double *theta, double *phi) {

  /*
    c=======================================================================
    subroutine pix2ang_nest(nside, ipix, theta, phi)
    c=======================================================================
    c     gives theta and phi corresponding to pixel ipix (NESTED) 
    c     for a parameter nside
    c=======================================================================
  */
    
      int npix, npface, face_num;
      int  ipf, ip_low, ip_trunc, ip_med, ip_hi;
      int     ix, iy, jrt, jr, nr, jpt, jp, kshift, nl4;
      double z, fn, fact1, fact2;
      double piover2=0.5*M_PI;
      int ns_max=8192;
      
      static int pix2x[1024], pix2y[1024];
      //      common /pix2xy/ pix2x, pix2y      
      
      
      if( nside<1 || nside>ns_max ) {
	fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
	exit(0);
      }
      npix = 12 * nside*nside;
      if( ipix<0 || ipix>npix-1 ) {
	fprintf(stderr, "%s (%d): ipix out of range: %ld\n", __FILE__, __LINE__, ipix);
	exit(0);
      }

      /* initiates the array for the pixel number -> (x,y) mapping */
      if( pix2x[1023]<=0 ) mk_pix2xy(pix2x,pix2y);

      fn = 1.*nside;
      fact1 = 1./(3.*fn*fn);
      fact2 = 2./(3.*fn);
      nl4   = 4*nside;

      //c     finds the face, and the number in the face
      npface = nside*nside;

      face_num = ipix/npface;//  ! face number in {0,11}
      ipf = (int)fmod(ipix,npface);//  ! pixel number in the face {0,npface-1}

      //c     finds the x,y on the face (starting from the lowest corner)
      //c     from the pixel number
      ip_low = (int)fmod(ipf,1024);//       ! content of the last 10 bits
      ip_trunc =   ipf/1024 ;//       ! truncation of the last 10 bits
      ip_med = (int)fmod(ip_trunc,1024);//  ! content of the next 10 bits
      ip_hi  =     ip_trunc/1024   ;//! content of the high weight 10 bits

      ix = 1024*pix2x[ip_hi] + 32*pix2x[ip_med] + pix2x[ip_low];
      iy = 1024*pix2y[ip_hi] + 32*pix2y[ip_med] + pix2y[ip_low];

      //c     transforms this in (horizontal, vertical) coordinates
      jrt = ix + iy;//  ! 'vertical' in {0,2*(nside-1)}
      jpt = ix - iy;//  ! 'horizontal' in {-nside+1,nside-1}

      //c     computes the z coordinate on the sphere
      //      jr =  jrll[face_num+1]*nside - jrt - 1;//   ! ring number in {1,4*nside-1}
      jr =  jrll[face_num]*nside - jrt - 1;
      //      cout << "face_num=" << face_num << endl;
      //      cout << "jr = " << jr << endl;
      //      cout << "jrll(face_num)=" << jrll[face_num] << endl;
      //      cout << "----------------------------------------------------" << endl;
      nr = nside;//                  ! equatorial region (the most frequent)
      z  = (2*nside-jr)*fact2;
      kshift = (int)fmod(jr - nside, 2);
      if( jr<nside ) { //then     ! north pole region
         nr = jr;
         z = 1. - nr*nr*fact1;
         kshift = 0;
      }
      else {
	if( jr>3*nside ) {// then ! south pole region
         nr = nl4 - jr;
         z = - 1. + nr*nr*fact1;
         kshift = 0;
	}
      }
      *theta = acos(z);
      
      //c     computes the phi coordinate on the sphere, in [0,2Pi]
      //      jp = (jpll[face_num+1]*nr + jpt + 1 + kshift)/2;//  ! 'phi' number in the ring in {1,4*nr}
      jp = (jpll[face_num]*nr + jpt + 1 + kshift)/2;
      if( jp>nl4 ) jp = jp - nl4;
      if( jp<1 )   jp = jp + nl4;

      *phi = (jp - (kshift+1)*0.5) * (piover2 / nr);

}

void ang2pix_nest( const long nside, double theta, double phi, long *ipix) {

  /* =======================================================================
   * subroutine ang2pix_nest(nside, theta, phi, ipix)
   * =======================================================================
   * gives the pixel number ipix (NESTED) corresponding to angles theta and phi
   *
   * the computation is made to the highest resolution available (nside=8192)
   * and then degraded to that required (by integer division)
   * this doesn't cost more, and it makes sure that the treatement of round-off 
   * will be consistent for every resolution
   * =======================================================================
   */
  
  double z, za, z0, tt, tp, tmp;
  int    face_num,jp,jm;
  long   ifp, ifm;
  int    ix, iy, ix_low, ix_hi, iy_low, iy_hi, ipf, ntt;
  double piover2 = 0.5*M_PI, pi = M_PI, twopi = 2.0*M_PI;
  int    ns_max = 8192;
  static int x2pix[128], y2pix[128];
  static char setup_done = 0;
  
  if( nside<1 || nside>ns_max ) {
    fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
    exit(0);
  }
  if( theta<0. || theta>pi ) {
    fprintf(stderr, "%s (%d): theta out of range: %f\n", __FILE__, __LINE__, theta);
    exit(0);
  }
  if( !setup_done ) {
    mk_xy2pix(x2pix,y2pix);
    setup_done = 1;
  }
  
  z  = cos(theta);
  za = fabs(z);
  z0 = 2./3.;
  if( phi>=twopi ) phi = phi - twopi;
  if( phi<0. )    phi = phi + twopi;
  tt = phi / piover2; /* in [0,4[ */
  
  if( za<=z0 ) { /* equatorial region */
    
    /* (the index of edge lines increase when the longitude=phi goes up) */
    jp = (int)floor(ns_max*(0.5 + tt - z*0.75)); /* ascending edge line index */
    jm = (int)floor(ns_max*(0.5 + tt + z*0.75)); /* descending edge line index */
    
    /* finds the face */
    ifp = jp / ns_max; /* in {0,4} */
    ifm = jm / ns_max;
    
    if( ifp==ifm ) face_num = (int)fmod(ifp,4) + 4; /* faces 4 to 7 */
    else if( ifp<ifm ) face_num = (int)fmod(ifp,4); /* (half-)faces 0 to 3 */
    else face_num = (int)fmod(ifm,4) + 8;           /* (half-)faces 8 to 11 */
    
    ix = (int)fmod(jm, ns_max);
    iy = ns_max - (int)fmod(jp, ns_max) - 1;
  }
  else { /* polar region, za > 2/3 */
    
    ntt = (int)floor(tt);
    if( ntt>=4 ) ntt = 3;
    tp = tt - ntt;
    tmp = sqrt( 3.*(1. - za) ); /* in ]0,1] */
    
    /* (the index of edge lines increase when distance from the closest pole
     * goes up)
     */
    /* line going toward the pole as phi increases */
    jp = (int)floor( ns_max * tp          * tmp ); 

    /* that one goes away of the closest pole */
    jm = (int)floor( ns_max * (1. - tp) * tmp );
    jp = (int)(jp < ns_max-1 ? jp : ns_max-1);
    jm = (int)(jm < ns_max-1 ? jm : ns_max-1);
    
    /* finds the face and pixel's (x,y) */
    if( z>=0 ) {
      face_num = ntt; /* in {0,3} */
      ix = ns_max - jm - 1;
      iy = ns_max - jp - 1;
    }
    else {
      face_num = ntt + 8; /* in {8,11} */
      ix =  jp;
      iy =  jm;
    }
  }
  
  ix_low = (int)fmod(ix,128);
  ix_hi  =     ix/128;
  iy_low = (int)fmod(iy,128);
  iy_hi  =     iy/128;

  ipf = (x2pix[ix_hi]+y2pix[iy_hi]) * (128 * 128)+ (x2pix[ix_low]+y2pix[iy_low]);
  ipf = (long)(ipf / pow(ns_max/nside,2));     /* in {0, nside**2 - 1} */
  *ipix =(long)( ipf + face_num*pow(nside,2)); /* in {0, 12*nside**2 - 1} */
}

void ang2pix_ring( const long nside, double theta, double phi, long *ipix) {
  /*
    c=======================================================================
    c     gives the pixel number ipix (RING) 
    c     corresponding to angles theta and phi
    c=======================================================================
  */
  
  int nl2, nl4, ncap, npix, jp, jm, ipix1;
  double  z, za, tt, tp, tmp;
  int ir, ip, kshift;
  
  double piover2 = 0.5*M_PI;
  double PI=M_PI;
  double twopi=2.0*M_PI;
  double z0=2.0/3.0;
  long ns_max=8192;
  
  if( nside<1 || nside>ns_max ) {
    fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
    exit(0);
  }
  
  if( theta<0. || theta>PI) {
    fprintf(stderr, "%s (%d): theta out of range: %f\n", __FILE__, __LINE__, theta);
    exit(0);
  }
  
  z = cos(theta);
  za = fabs(z);
  if( phi >= twopi)  phi = phi - twopi;
  if (phi < 0.)     phi = phi + twopi;
  tt = phi / piover2;//  ! in [0,4)
  
  nl2 = 2*nside;
  nl4 = 4*nside;
  ncap  = nl2*(nside-1);// ! number of pixels in the north polar cap
  npix  = 12*nside*nside;
  
  if( za <= z0 ) {
    
    jp = (int)floor(nside*(0.5 + tt - z*0.75)); /*index of ascending edge line*/
    jm = (int)floor(nside*(0.5 + tt + z*0.75)); /*index of descending edge line*/
    
    ir = nside + 1 + jp - jm;// ! in {1,2n+1} (ring number counted from z=2/3)
    kshift = 0;
    if (fmod(ir,2)==0.) kshift = 1;// ! kshift=1 if ir even, 0 otherwise
    
    ip = (int)floor( ( jp+jm - nside + kshift + 1 ) / 2 ) + 1;// ! in {1,4n}
    if( ip>nl4 ) ip = ip - nl4;
    
    ipix1 = ncap + nl4*(ir-1) + ip ;
  }
  else {
    
    tp = tt - floor(tt);//      !MOD(tt,1.d0)
    tmp = sqrt( 3.*(1. - za) );
    
    jp = (int)floor( nside * tp * tmp );// ! increasing edge line index
    jm = (int)floor( nside * (1. - tp) * tmp );// ! decreasing edge line index
    
    ir = jp + jm + 1;//        ! ring number counted from the closest pole
    ip = (int)floor( tt * ir ) + 1;// ! in {1,4*ir}
    if( ip>4*ir ) ip = ip - 4*ir;
    
    ipix1 = 2*ir*(ir-1) + ip;
    if( z<=0. ) {
      ipix1 = npix - 2*ir*(ir+1) + ip;
    }
  }
  *ipix = ipix1 - 1;// ! in {0, npix-1}
  
}



int xyf2ring (int ix, int iy, int nside_,int face_num)
{
  int nl4 = 4*nside_;
  int jr = (jrll[face_num]*nside_) - ix - iy  - 1;
  int npix_ = nside2npix(nside_);
  int nr, kshift, n_before;
  int ncap_ = 2*nside_*(nside_-1);

  if (jr<nside_)
    {
      nr = jr;
      n_before = 2*nr*(nr-1);
      kshift = 0;
    }
  else if (jr > 3*nside_)
    {
      nr = nl4-jr;
      n_before = npix_ - 2*(nr+1)*nr;
      kshift = 0;
    }
  else
    {
      nr = nside_;
      n_before = ncap_ + (jr-nside_)*nl4;
      kshift = (jr-nside_)&1;
    }
  
  int jp = (jpll[face_num]*nr + ix - iy + 1 + kshift) / 2;
  if (jp>nl4)
    jp-=nl4;
  else
    if (jp<1) jp+=nl4;
  
  return n_before + jp - 1;
}


void ring2xyf (int pix, long nside_,int *ix, int *iy, int *face_num)
{
  int iring, iphi, kshift, nr;

  int nl2 = 2*nside_;
  int ncap_ = nl2*(nside_-1);
  int order_;// = nside2order(nside_);
  int npix_;// = nside2npix(nside_);

  order_ = nside2order(nside_);
  npix_ =nside2npix(nside_);

  if (pix<ncap_) // North Polar cap
    {
      iring = (int)(0.5*(1+sqrt(1+2*pix))); //counted from North pole
    iphi  = (pix+1) - 2*iring*(iring-1);
    kshift = 0;
    nr = iring;
    *face_num=0;
    int tmp = iphi-1;
    if (tmp>=(2*iring))
      {
      *face_num=2;
      tmp-=2*iring;
      }
    if (tmp>=iring) ++(*face_num);
    }
  else if (pix<(npix_-ncap_)) // Equatorial region
    {
    int ip = pix - ncap_;
    if (order_>=0)
      {
      iring = (ip>>(order_+2)) + nside_; // counted from North pole
      iphi  = (ip&(4*nside_-1)) + 1;
      }
    else
      {
      iring = (ip/(4*nside_)) + nside_; // counted from North pole
      iphi  = (ip%(4*nside_)) + 1;
      }
    kshift = (iring+nside_)&1;
    nr = nside_;
    int ire = iring-nside_+1;
    int irm = nl2+2-ire;
    int ifm, ifp;
    if (order_>=0)
      {
      ifm = (iphi - ire/2 + nside_ -1) >> order_;
      ifp = (iphi - irm/2 + nside_ -1) >> order_;
      }
    else
      {
      ifm = (iphi - ire/2 + nside_ -1) / nside_;
      ifp = (iphi - irm/2 + nside_ -1) / nside_;
      }
    if (ifp == ifm) // faces 4 to 7
      *face_num = (ifp==4) ? 4 : ifp+4;
    else if (ifp<ifm) // (half-)faces 0 to 3
      *face_num = ifp;
    else // (half-)faces 8 to 11
      *face_num = ifm + 8;
    }
  else // South Polar cap
    {
    int ip = npix_ - pix;
    iring = (int)(0.5*(1+sqrt(2*ip-1))); //counted from South pole
    iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));
    kshift = 0;
    nr = iring;
    iring = 2*nl2-iring;
    *face_num=8;
    int tmp = iphi-1;
    if (tmp>=(2*nr))
      {
      *face_num=10;
      tmp-=2*nr;
      }
    if (tmp>=nr) ++(*face_num);
    }

  int irt = iring - (jrll[*face_num]*nside_) + 1;
  int ipt = 2*iphi- jpll[*face_num]*nr - kshift -1;
  if (ipt>=nl2) ipt-=8*nside_;

  *ix =  (ipt-irt) >>1;
  *iy =(-(ipt+irt))>>1;
}

void nest2xyf (int pix, int nside_, int *ix, int *iy, int *face_num)
{
  int order_ =  nside2order(nside_);
  int npface_ = nside_*nside_;

  *face_num = pix>>(2*order_);
  pix2xy(pix&(npface_-1),ix,iy);
}

int xyf2nest (int ix, int iy, int nside_, int face_num)
{
  int order_ =  nside2order(nside_);

  return (face_num<<(2*order_))+xy2pix(ix,iy);
}

void fill_tab(void)
{
  int m;

  ctab = malloc (0x100*sizeof(short));
  utab = malloc (0x100*sizeof(short));

  for (m=0; m<0x100; ++m)
    {
      ctab[m] =
	(m&0x1 )       | ((m&0x2 ) << 7) | ((m&0x4 ) >> 1) | ((m&0x8 ) << 6)
	| ((m&0x10) >> 2) | ((m&0x20) << 5) | ((m&0x40) >> 3) | ((m&0x80) << 4);
      utab[m] =
	(m&0x1 )       | ((m&0x2 ) << 1) | ((m&0x4 ) << 2) | ((m&0x8 ) << 3)
	| ((m&0x10) << 4) | ((m&0x20) << 5) | ((m&0x40) << 6) | ((m&0x80) << 7);
    }
}

void pix2xy (int pix, int *x, int *y)
{
  int raw = (pix&0x5555) | ((pix&0x55550000)>>15);
  if (ctab==NULL) fill_tab();

  *x = ctab[raw&0xff] | (ctab[raw>>8]<<4);
  raw = ((pix&0xaaaa)>>1) | ((pix&0xaaaa0000)>>16);
  *y = ctab[raw&0xff] | (ctab[raw>>8]<<4);
}

int xy2pix (int x, int y)
{
  if (utab==NULL) fill_tab();
  return utab[x&0xff] | (utab[x>>8]<<16) | (utab[y&0xff]<<1) | (utab[y>>8]<<17);
}


 int pixNeighbours(int pix, long nside, int nested, int withDiag, int *result, int *nonDiagId)
{
  
  static const int xoffset[] = { -1,-1, 0, 1, 1, 1, 0,-1 };
  static const int yoffset[] = {  0, 1, 1, 1, 0,-1,-1,-1 };
  static const int facearray[][12] =
    { {  8, 9,10,11,-1,-1,-1,-1,10,11, 8, 9 },   // S
      {  5, 6, 7, 4, 8, 9,10,11, 9,10,11, 8 },   // SE
      { -1,-1,-1,-1, 5, 6, 7, 4,-1,-1,-1,-1 },   // E
      {  4, 5, 6, 7,11, 8, 9,10,11, 8, 9,10 },   // SW
      {  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11 },   // center
      {  1, 2, 3, 0, 0, 1, 2, 3, 5, 6, 7, 4 },   // NE
      { -1,-1,-1,-1, 7, 4, 5, 6,-1,-1,-1,-1 },   // W
      {  3, 0, 1, 2, 3, 0, 1, 2, 4, 5, 6, 7 },   // NW
      {  2, 3, 0, 1,-1,-1,-1,-1, 0, 1, 2, 3 } }; // N
  static const int swaparray[][12] =
    { {  0,0,0,0,0,0,0,0,3,3,3,3 },   // S
      {  0,0,0,0,0,0,0,0,6,6,6,6 },   // SE
      {  0,0,0,0,0,0,0,0,0,0,0,0 },   // E
      {  0,0,0,0,0,0,0,0,5,5,5,5 },   // SW
      {  0,0,0,0,0,0,0,0,0,0,0,0 },   // center
      {  5,5,5,5,0,0,0,0,0,0,0,0 },   // NE
      {  0,0,0,0,0,0,0,0,0,0,0,0 },   // W
      {  6,6,6,6,0,0,0,0,0,0,0,0 },   // NW
      {  3,3,3,3,0,0,0,0,0,0,0,0 } }; // N

  int ix, iy, face_num;
  long nside_ = nside;
  int m,i;
  int skipped=10;
  int nnei=8;
  
  if (!nested) 
    ring2xyf(pix,nside_,&ix,&iy,&face_num); 
  else
    nest2xyf(pix,nside_,&ix,&iy,&face_num);

   const int nsm1 = nside_-1;
   if ((ix>0)&&(ix<nsm1)&&(iy>0)&&(iy<nsm1))
     {
     if (!nested)
       for (m=0; m<8; ++m)
         result[m] = xyf2ring(ix+xoffset[m],iy+yoffset[m],nside_,face_num);
     else
       for (m=0; m<8; ++m)
         result[m] = xyf2nest(ix+xoffset[m],iy+yoffset[m],nside_,face_num);
     }
   else
     {
     for (i=0,nnei=0; i<8; ++i)
       {
       int x=ix+xoffset[i];
       int y=iy+yoffset[i];
       int z;
       int nbnum=4;
       if (x<0)
         { x+=nside_; nbnum-=1; }
       else if (x>=nside_)
         { x-=nside_; nbnum+=1; }
       if (y<0)
         { y+=nside_; nbnum-=3; }
       else if (y>=nside_)
         { y-=nside_; nbnum+=3; }

       int f = facearray[nbnum][face_num];
       if (f>=0)
         {
         if (swaparray[nbnum][face_num]&1) x=nside_-x-1;
         if (swaparray[nbnum][face_num]&2) y=nside_-y-1;
         if (swaparray[nbnum][face_num]&4) {z=x;x=y;y=z;}
         result[nnei] = (!nested) ? xyf2ring(x,y,nside_,f) : xyf2nest(x,y,nside_,f);
	 nnei++;
         }
       else
	 {
	   skipped=i;
	   //result[i] = -1;
	   //printf ("Error:i=%d\n",i);
	 }
       }
     }

   (*nonDiagId)=(1<<0)|(1<<2)|(1<<4)|(1<<6);
   if (nnei<8) {
     if (skipped==1) (*nonDiagId)=(1<<0)|(1<<1)|(1<<3)|(1<<5);
     else if (skipped==3) (*nonDiagId)=(1<<0)|(1<<2)|(1<<3)|(1<<5);
     else if (skipped==5) (*nonDiagId)=(1<<0)|(1<<2)|(1<<4)|(1<<5);
     
   }
 
   if (!withDiag)
     {
       
       if (nnei==8)
	 {
	   result[1]=result[2];
	   result[2]=result[4];
	   result[3]=result[6];
	 }
       else
	 {
	   int m=0;
	   int s=(*nonDiagId);
	   for (i=0; i<8; i++)
	     if (s&(1<<i)) result[m++]=result[i];
	 }
       nnei=4;
     }  

   return nnei;
}
