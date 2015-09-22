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

#ifndef __CHEALPIX_H__
#define __CHEALPIX_H__

#ifdef __cplusplus
extern "C" {
#endif

/* -------------------- */
/* Constant Definitions */
/* -------------------- */

#ifndef HEALPIX_NULLVAL
#define HEALPIX_NULLVAL (-1.6375e30)
#endif /* HEALPIX_NULLVAL */

/* --------------------- */
/* Function Declarations */
/* --------------------- */
int is_healpix_map(const char *infile);

/* CXX pixel operations */
/* -------------------- */
void ring2xyf (int pix, long nside_,int *ix, int *iy, int *face_num);
int xyf2ring (int ix, int iy, int nside_,int face_num);
void nest2xyf (int pix, int nside_, int *ix, int *iy, int *face_num);
int xyf2nest (int ix, int iy, int nside_, int face_num);
void pix2xy (int pix, int *x, int *y);
int xy2pix (int x, int y);

/* pixel operations */
/* ---------------- */
void ang2pix_nest(const long nside, double theta, double phi, long *ipix);
void ang2pix_ring(const long nside, double theta, double phi, long *ipix);

void pix2ang_nest(long nside, long ipix, double *theta, double *phi);
void pix2ang_ring(long nside, long ipix, double *theta, double *phi);

void nest2ring(long nside, long ipnest, long *ipring);
void ring2nest(long nside, long ipring, long *ipnest);

void mk_pix2xy(int *pix2x, int *pix2y);
void mk_xy2pix(int *x2pix, int *y2pix);

long nside2npix(const long nside);
long npix2nside(const long pix  );

void ang2vec(double theta, double phi,   double *vec);
void vec2ang(double *vec, double *theta, double *phi);

void vec2pix_nest(const long nside, double *vec, long *ipix);
void vec2pix_ring(const long nside, double *vec, long *ipix);

void pix2vec_nest(long nside, long ipix, double *vec);
void pix2vec_ring(long nside, long ipix, double *vec);

  int pixNeighbours(int pix, long nside, int nested, int withDiag, int *result, int *nonDiagId);
  
 

/* FITS operations */
/* --------------- */

void printerror (int) ;

float *read_healpix_map_f (const char *, long *, char *, char *) ;
double *read_healpix_map_d (const char *, long *, char *, char *) ;

//int write_healpix_map( float *, long , char *, char ,char *) ;

long get_fits_size(char *, long *, char * ) ;


/* ------------------ */
/* end of header file */
/* ------------------ */

#ifdef __cplusplus
}
#endif

#endif /* __CHEALPIX_H__ */
