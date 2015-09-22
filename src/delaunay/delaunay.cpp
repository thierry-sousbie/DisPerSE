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
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include <iostream>
#include <fstream>

#include "gridType.hxx"
#include "box.hxx"
#include "sampledDataInput.hxx"

#include "gadget_io.h"
#include "mystring.h"
#include "NDfield.h"
#include "NDnetwork.h"
#include "mystring.h"
#include "distances.h"

#include "NDnet_interface.hxx"
#include "NDnet_gather.hxx"

#define GLOBAL_DEFINITION
#include "global.h"

using namespace std;

#define PERIODIC_MARGIN_SIZE 0.1
#define NONPERIODIC_MARGIN_SIZE 0.06
#define DEFAULT_MAX_ANGULAR_SIZE 5.0

#define MARGIN_SAFETY_FACTOR 5

void Usage(char *fname)
{
  int i;
  fprintf(stderr,"delaunay version %s\n",VERSION_STR);
  fprintf (stderr,"\nUsage:\n  %s <fname> [-outName <fname>] [-outDir <dir>]\n",CutName(fname));
#if PERIODICITY!=1
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
#if NDIMS==3
  fprintf (stderr,"   [-subbox <x0> <y0> <z0> <dx0> <dy0> <dz0>]\n");
#elif NDIMS==2
  fprintf (stderr,"   [-subbox <x0> <y0> <dx0> <dy0>]\n");
#endif
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   [-periodic] [-minimal]\n");
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   [-blocks <NChunks> <NThreads>]\n");
  //for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  //fprintf (stderr,"   [-periodic] [-margin <M=%.2f(p)/%.2f(np)>] [-btype <t=mirror>]\n",PERIODIC_MARGIN_SIZE,NONPERIODIC_MARGIN_SIZE);
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   [-margin <M>] [-btype <t=void>]\n");
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   [-mask <fname.ND>]\n");// [-mangle <mangle mask>]\n");
  //for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  //fprintf (stderr,"   [-smooth <N=0>]\n");
#ifdef HAVE_CFITS_IO
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   [-angmask <fname.fits> [<maximum angular size (degrees) = %.2f>]]\n",DEFAULT_MAX_ANGULAR_SIZE);
#endif
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   [-radialDensity <A> <Dr> <B>]\n");
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   [-subSample < 0<s<1 >]\n");
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   [-cosmo <Om=%.2f Ol=%.2f Ok=%.2f h=%.2f w=%.2f>]\n",OMEGAM_DEFAULT,OMEGAL_DEFAULT,OMEGAK_DEFAULT,HUBBLE_DEFAULT,W_DEFAULT);
  fprintf (stderr,"\n");
  fprintf (stderr,"   * '-minimal' prevents explicit computation of intermediate cells\n");
  fprintf (stderr,"      such as segments and triangles in 3D.\n");
  fprintf (stderr,"   * '-blocks' computes the tesselation in chunks, using NThreads in parallel.\n");
  fprintf (stderr,"     This is also usefull for significantly lowering memory consumption\n");
  fprintf (stderr,"     when NChunks >> NThreads.\n");
  fprintf (stderr,"     INCOMPATIBLE with '-mask', '-btype smooth', '-angmask', '-radialDensity'.\n");
  //for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   * Available boundary types: mirror, periodic, smooth, void.\n");
  //for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   * Margin is expressed as a fraction of the box (or subbox) size.\n");
  fprintf (stderr,"   * '-angmask' : apply a Healpix type angular mask.\n");
  fprintf (stderr,"      A hole smaller than max. angular size is filled, while a larger one \n");
  fprintf (stderr,"      is considered outside the boundary. Default value (%.2f deg)\n",DEFAULT_MAX_ANGULAR_SIZE);
  fprintf (stderr,"      is adequate for SDSS. \n");
  fprintf (stderr,"   * '-radialDensity' : correct the mean number of galaxy N as a function\n");
  fprintf (stderr,"     of distance d: \n");
  fprintf (stderr,"             N(d) = A * d^2 * exp(-(d/Dr)^B).\n");
#endif
  
  printf ("\n");
  exit(0);
}

int main(int argc, char **argv)
{  
  verbose=1;

  int N;
  long i,j,k,l;

  double *x0=NULL;
  double *delta=NULL;
  std::vector<double> x0_vec;
  std::vector<double> delta_vec;
  std::vector<double> x0_full;
  std::vector<double> delta_full;


  std::vector<double> sub_x0;
  std::vector<double> sub_delta;
  int npart;
  float *pos;
  int ndims; 
  char fname[255];
  char fullfname[255];
  char outname[255];
  char outdir[255];
  char maskfname[255];
  char manglefname[255];
  char angMaskfname[255];
  double dt;
  NDnetwork *net;
  struct timeval wc_time1,wc_time2;

  int Opt_fname=0;
  int Opt_mangle=0;
  int Opt_angMask=0;
  int Opt_mask=0;
  int Opt_periodic=0;
  //int Opt_margin_set=0;
  int Opt_nSmooth=0;
  double Opt_margin=NONPERIODIC_MARGIN_SIZE;  
  int Opt_guard=0;
  int Opt_smoothBoundaries=0;
  int Opt_radialDensity=0;
  double Opt_subSample=1.0;
  int Opt_outFName=0;
  int Opt_outDir=0;
  int Opt_computeAllFaces=1;  
  double radialDensity_A=1.;
  double radialDensity_Dr=21300;
  double radialDensity_B=1.55;

  bool margin_set=false;
  double borderSize=0;
  double big_borderSize=0;
  double Opt_maxAngSize=DEFAULT_MAX_ANGULAR_SIZE;
  int return_val=0;
  int Opt_blocks=0;
  int nblocks=1;
  int nthreads=1;

  bConditions::BoundaryType BType = bConditions::BType_None;  

  if (argc==1) Usage(argv[0]);
  printf("\n");  
  for (i=1;i<argc;)
    {
	if (argv[i][0]!='-')
	{
	  strcpy(fullfname,argv[i]);
	  Opt_fname=1;
	  i++;
	}
	else if (i==1)
	  {
	    fprintf (stderr,"First argument must be a filename.\n");
	    Usage(argv[0]);
	  }
	else if (!strcmp(argv[i],"-cosmo"))
	  {
	    double val[6];
	    
	    if (argc<i+6)
	    {
	      fprintf(stderr,"Invalid arguments for option '-cosmo'.\n");
	      Usage(argv[0]);
	    }
	    
	    for (j=1;j<5;j++) val[j-1]=atof(argv[i+j]);
	    if (!isFile(argv[i+5]))
	      cosmoD_init(val[0],val[1],val[2],val[3],atof(argv[i+5]),NULL);
	    else
	      cosmoD_init(val[0],val[1],val[2],val[3],-1,argv[i+5]);

	    i+=6;
	  }
	else if (!strcmp(argv[i],"-blocks"))
	  {
	    Opt_blocks=1;
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect two integers as argument of '-blocks'\n");
		Usage(argv[0]);
	    }
	    nblocks=atoi(argv[i++]);
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect two integers as argument of '-blocks'\n");
		Usage(argv[0]);
	    }
	    nthreads=atoi(argv[i++]);	    
	  }
	else if (!strcmp(argv[i],"-outName"))
	  {
	    Opt_outFName=1;
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a filename as argument of '-outName'\n");
		Usage(argv[0]);
	    }
	    strcpy(outname,argv[i]);
	    i++;
	  }
	else if (!strcmp(argv[i],"-outDir"))
	  {
	    Opt_outDir=1;
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a filename as argument of '-outDir'\n");
		Usage(argv[0]);
	    }
	    strcpy(outdir,argv[i]);
	    i++;
	  }
#if PERIODICITY!=1
	else if (!strcmp(argv[i],"-subbox"))
	{
	  i++;	  
	  for (j=0;j<NDIMS;j++)
	    {
	      if (argc==i+j)
		{
		  fprintf (stderr,"I Expect %d numbers as argument of '-subbox'\n",NDIMS*2);
		  Usage(argv[0]);
		}
	      sub_x0.push_back(atof(argv[i+j]));
	    }
	  i+=NDIMS;

	  for (j=0;j<NDIMS;j++)
	    {
	      if (argc==i+j)
		{
		  fprintf (stderr,"I Expect %d numbers as argument of '-subbox'\n",NDIMS*2);
		  Usage(argv[0]);
		}
	      sub_delta.push_back(atof(argv[i+j]));
	    }
	  i+=NDIMS;
	  
	  //i++;
	}
	else if (!strcmp(argv[i],"-mask"))
	{
	    Opt_mask=1;	    
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a filename as argument of '-mask'\n");
		Usage(argv[0]);
	    }
	    strcpy(maskfname,argv[i]);
	    i++;
	}
	/*
	else if (!strcmp(argv[i],"-smooth"))
	{
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect an integer as argument of '-smooth'\n");
		Usage(argv[0]);
	    }
	    Opt_nSmooth = atoi(argv[i]);
	    i++;
	}
	*/
	else if (!strcmp(argv[i],"-subSample"))
	{
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a float within [0,1.] as argument of '-subSample'\n");
		Usage(argv[0]);
	    }
	    Opt_subSample = atof(argv[i]);
	    if ((Opt_subSample<0)||(Opt_subSample>1))
	    {
		fprintf (stderr,"I Expect a float within [0,1.] as argument of '-subSample'\n");
		Usage(argv[0]);
	    }
	    i++;
	}
	else if (!strcmp(argv[i],"-periodic"))
	{
	    Opt_periodic=1;
	    BType = bConditions::BType_Periodic;
	    i++;
	}
	else if (!strcmp(argv[i],"-minimal"))
	{
	  Opt_computeAllFaces=0;
	  i++;
	}
	else if (!strcmp(argv[i],"-mangle"))
	{
	    Opt_mangle=1;
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a filename as argument of '-mangle'\n");
		Usage(argv[0]);
	    }
	    strcpy(manglefname,argv[i]);
	    i++;
	}
	else if (!strcmp(argv[i],"-angmask"))
	{
	    Opt_angMask=1;
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a 'fits' filename as argument of '-angmask'\n");
		Usage(argv[0]);
	    }
	    strcpy(angMaskfname,argv[i]);
	    i++;

	    if ((argc!=i)&&(argv[i][0]!='-'))
	    {
	      Opt_maxAngSize = atof(argv[i]);
	      i++;
	    }
	    
	}
	else if (!strcmp(argv[i],"-radialDensity"))
	  {
	    Opt_radialDensity = 1;
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect 3 arguments for option '-radialDensity'\n");
		Usage(argv[0]);
	    }
	    radialDensity_A = atof(argv[i]);i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect other arguments for option '-radialDensity'\n");
		Usage(argv[0]);
	    }
	    radialDensity_Dr = atof(argv[i]);i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect other arguments for option '-radialDensity'\n");
		Usage(argv[0]);
	    }
	    radialDensity_B = atof(argv[i]);i++;
	  }

	else if (!strcmp(argv[i],"-margin"))
	{
	  //Opt_margin_set=1;
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a number as argument of '-margin'\n");
		Usage(argv[0]);
	    }
	    Opt_margin = atof(argv[i]);
	    /*
	    if (Opt_margin>=0.5)
	      {
		fprintf (stderr,"Margin must be < 0,5. Value set to ~0,5.\n");
		Opt_margin=0.49999999999;
	      }
	    */
	    i++;
	    margin_set=true;
	}
	else if (!strcmp(argv[i],"-btype"))
	{
	    
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a string argument for '-btype'.\n");
		fprintf (stderr,"  available types: periodic, mirror, smooth, void.\n");
		Usage(argv[0]);
	    }
	    if (strcmp(argv[i],"mirror")==0)
	      BType = bConditions::BType_Mirror;
	    else if (strcmp(argv[i],"periodic")==0)
	      BType = bConditions::BType_Periodic;
	    else if (strcmp(argv[i],"smooth")==0) {
	      BType = bConditions::BType_None;
	      Opt_smoothBoundaries=1;
	    }
	    else if (strcmp(argv[i],"void")==0)
	      BType = bConditions::BType_None;
	    else if (strcmp(argv[i],"none")==0)
	      BType = bConditions::BType_None;
	    else
	      {
		fprintf (stderr,"%s is not a valid boundary type.\n",argv[i]);
		fprintf (stderr,"  available types: periodic, mirror, smooth, void.\n");
		Usage(argv[0]);
	      }
	    i++;
	}
#endif
	else 
	{
	    printf ("\nWhat is %s ???\n",argv[i]);
	    Usage(argv[0]);
	}
    }
  
  if (!Opt_fname) Usage(argv[0]);

  if (!Opt_outFName) 
    {
      if (Opt_subSample<1.)
	sprintf(outname,"%s_sub%.2f",CutName(fullfname),Opt_subSample*100);
      else
	strcpy(outname,CutName(fullfname));    
    }

  if (Opt_outDir) {
    char tmp[255];
    if (outdir[strlen(outdir)-1]!='/') sprintf(outdir,"%s/",outdir);
    sprintf(tmp,"%s%s",outdir,outname);
    strcpy(outname,tmp);
  }

  if ((Opt_periodic)&&(BType != bConditions::BType_Periodic))
    {
      fprintf(stderr,"Periodic boundary type was enforced by option'-periodic'.\n");
      BType = bConditions::BType_Periodic;
    }

  //strcpy(fname,CutName(fullfname));

  sampledDataInput data(fullfname);

  if (Opt_subSample<1.) 
    {
      printf ("Subsampling (f=%.2g) ... ",Opt_subSample);fflush(0);
      data.subSample(Opt_subSample);  
      printf ("done.\n");
    }
  
  if (Opt_mangle)
    {
      if (data.getCoordSystem() != sampledDataInput::cs_spherical)
	{
	  fprintf(stderr,"ERROR: option 'mangle' can only be used in a spherical coordinates system.\n");
	  Usage(argv[0]);
	}
      if (system("drangle > /dev/null")) 
	{
	  fprintf(stderr,"ERROR: 'mangle' is not correctly installed.\n");
	  Usage(argv[0]);
	}
    }
    

  if (!data.isLoaded()) {
    fprintf(stderr,"ERROR: Unknown or invalid data type.\n");
    return -1;
  }
  
  /*
  if ((margin_set)&&(BType == bConditions::BType_None))
    {
      printf("WARNING: boundary conditons are set to 'void'. Option '-margin' will be ignored.");
    }
  */
  pos=data.getPosFPtr();
  npart=data.getNPart();
  ndims=data.getNDims();
  data.getBBox(&x0,&delta);
  //printf("size = %d / %d\n",NDIMS,ndims);
  x0_full.assign(x0,x0+NDIMS);
  delta_full.assign(delta,delta+NDIMS);

  //if ((!margin_set)&&(Opt_periodic)) Opt_margin=PERIODIC_MARGIN_SIZE;
  
  bConditions boundary(BType,x0,delta,NDIMS);

  
  if (sub_x0.size())
    {
      //if (Opt_periodic) printf ("Option '-periodic' was set with a subbox.\n");
      Opt_periodic=0;

      for (i=0;i<ndims;i++)
	{
	  x0[i]=sub_x0[i];
	  delta[i]=sub_delta[i];
	}
    }
  
  x0_vec.assign(x0,x0+NDIMS);
  delta_vec.assign(delta,delta+NDIMS);

  if (!margin_set) 
    {
      double vol = delta_full[0];
      for (i=1;i<ndims;i++) vol*=delta_full[i];
      vol/=npart;      
      borderSize=big_borderSize=MARGIN_SAFETY_FACTOR*pow(vol,1./ndims);
      Opt_margin=borderSize/delta[0];
      for (i=1;i<ndims;i++)
	if (Opt_margin > borderSize/delta[i]) 
	  Opt_margin=borderSize/delta[i];
      printf("Margin set to %g (actual size is ~%g).\n",Opt_margin,borderSize);
    }
  else
    {
      borderSize=0;
      for (i=0;i<ndims;i++)
	if (delta[i]*Opt_margin > borderSize) 
	  borderSize=delta[i]*Opt_margin;
      
      for (i=0;i<ndims;i++)
	if (delta_full[i]*Opt_margin > big_borderSize) 
	  big_borderSize=delta_full[i]*Opt_margin;
    }

  if (NDIMS!=ndims)
  {
      fprintf (stderr,"ERROR: Trying to load %dD data, progam was compiled for %dD.\n",ndims,NDIMS);
      exit(0);
  }  

  
  
  gettimeofday(&wc_time1, NULL);

//#if PERIODICITY!=1
    

  if (Opt_blocks)
    {
      if (Opt_angMask) 
	{fprintf(stderr,"ERROR: option '-blocks' is incompatible with '-angmask'.\n");exit(-1);}
      if (Opt_smoothBoundaries) 
	{fprintf(stderr,"ERROR: option '-blocks' is incompatible with '-btype smooth'.\n");exit(-1);}
      if (Opt_mask) 
	{fprintf(stderr,"ERROR: option '-blocks' is incompatible with '-mask'.\n");exit(-1);}
      if (Opt_radialDensity)
	{fprintf(stderr,"ERROR: option '-blocks' is incompatible with '-radialDensity'.\n");exit(-1);}

      printf ("Tesselating %d particles (%dD) using block decomposition:\n",npart,ndims);fflush(0);

      ThreadedDelaunay::params threadedP;
      threadedP.boundary=&boundary;
      threadedP.pos=pos;
      threadedP.mass=data.getMassFPtr();     
      threadedP.npart=npart;
      threadedP.ndims=ndims;
      threadedP.x0=x0;
      threadedP.delta=delta;
      threadedP.borderSize=borderSize;
      threadedP.periodic=Opt_periodic;
      threadedP.outname=outname;
      threadedP.margin=Opt_margin;
      threadedP.checkBoundaries=(BType != bConditions::BType_None);
      threadedP.computeAllFaces = (Opt_computeAllFaces!=0)?true:false;
      
      ThreadedDelaunay TD;
      TD.compute(&threadedP,nblocks,nthreads);
      data.clear();
      char GFName[255];
      sprintf(fname,"%s%s",outname,ndnet::IO::getExtension().c_str());
      sprintf(GFName,"%s_G%s",outname,ndnet::IO::getExtension().c_str());      
      //NDnetGather gather(fname);
      //gather.write(GFName);
      NDnetGather(fname).write(GFName);
      
      exit(0);
    }
  
  printf ("Tesselating %d particles (%dD) ... ",npart,ndims);fflush(0);

  Delaunay *T = new Delaunay(x0,delta,Opt_periodic);

  //Delaunay T(x0,delta,Opt_periodic);

  std::vector<Delaunay::Point> pts;
  std::vector<long> index;
  std::vector<INT> true_index;

  std::vector<double> x0v;
  std::vector<double> deltav;
  std::vector<double> x0v_old;
  std::vector<double> deltav_old;
 
  long nnew=0;
  long old_nnew=0;
  bool needRecomp;
  int nPasses=0;
  
  do {
    needRecomp=false;
    nPasses++;

    x0v_old=x0v;
    deltav_old=deltav;
    x0v.assign(x0,x0+NDIMS);
    deltav.assign(delta,delta+NDIMS);
    
    for (i=0;i<ndims;i++)
      {
	x0v[i]-=borderSize;
	deltav[i]+=2*borderSize;
      }

    for (i=0,k=0;i<npart;i++) {
      std::vector<double> P(NDIMS);
      std::vector<double> result;      
      for (j=0;j<NDIMS;j++) P[j]=pos[i*NDIMS+j];
      bool direct = boundary.toSubBox(P, result, x0v, deltav);
    
      for (j=0;j<result.size()/NDIMS;j++)
	{

	  if (nPasses>1)
	    {
	      for (l=0;l<ndims;l++)
		{
		  if ((result[j*NDIMS+l]<x0v_old[l])||
		      (result[j*NDIMS+l]>x0v_old[l]+deltav_old[l]))
		    break;
		}
	      if (l==ndims) continue;
	    }

#if NDIMS==3
	  pts.push_back(Delaunay::Point(result[j*NDIMS],result[j*NDIMS+1],result[j*NDIMS+2]));
#else
	  pts.push_back(Delaunay::Point(result[j*NDIMS+0],result[j*NDIMS+1]));
#endif
	  if ((j==0)&&direct)
	    {
	      true_index.push_back(i);
	      index.push_back(i);
	    }
	  else
	    {
	      true_index.push_back(i);
	      index.push_back(npart+nnew);
	      nnew++;
	    }
	}
      k++;
    }
    
    if (nPasses==1) {printf("(+%ld in boundary)",nnew);fflush(0);}
    else {printf ("Adding %ld guard particles ... ",nnew-old_nnew);fflush(0);}
    old_nnew=nnew;
    T->insert(pts.begin(),pts.end(),data.getMassFPtr(),NULL,&(index[0]),&(true_index[0]));

    gettimeofday(&wc_time2, NULL);
    dt  = (double)wc_time2.tv_sec + ((double)wc_time2.tv_usec)*1.E-6;
    dt -= (double)wc_time1.tv_sec + ((double)wc_time1.tv_usec)*1.E-6;
    printf(" done. (%.2fs elapsed)\n",dt);fflush(0);


    if (BType != bConditions::BType_None)
      {
	printf ("Identifying boundaries ... ");fflush(0);
	std::pair<int,double> bcheck=T->CheckBoundaries(borderSize);
	printf ("done.\n");

	if (bcheck.first)
	  {
	    //printf("WARNING: circumsphere test failed for %d simplexes.\n",bcheck.first);
	    //printf("  This means that the margin size was not correctly estimated.\n");
	    //printf("  Current margin is %e  = %.2f%%/%.2f%% of box/subbox size.\n",
	    //   borderSize, Opt_margin*100*(borderSize/big_borderSize), Opt_margin*100);
	    double suggestedBorder=(bcheck.second*1.001);
	    if (suggestedBorder>2*borderSize) suggestedBorder=2*borderSize;
	    double suggestedMargin=suggestedBorder/borderSize*Opt_margin;	
	    //printf("  Will try ~%e.  (option '-margin %.5f')\n",
	    //   suggestedBorder,suggestedMargin);
	    needRecomp=true;

	    borderSize = suggestedBorder;
	    Opt_margin = suggestedMargin;
	  }
	//else printf("Circumsphere test was successful.\n");
      }
    else printf ("Circumsphere test was skipped.\n");
    
    pts.clear();
    index.clear();
    true_index.clear();

  } while (needRecomp);

  /*
#else //PERIODICITY==1
#if NDIMS==3
  if ((delta[0]!=delta[1])||(delta[1]!=delta[2])) {
    printf("\nSORRY: Periodic boundary conditions work only on cubic boxes. Use 'Delaunay_3D -periodic' instead.\n");
    exit(0);
  }
  Iso_cuboid domain(x0[0],x0[1],x0[2],x0[0]+delta[0],x0[1]+delta[1],x0[2]+delta[2]); 
#else //NDIMS!=3
  if (delta[0]!=delta[1]) {
    printf("\nSORRY: Periodic boundary conditions work only on cubic boxes. Use 'Delaunay_xD -periodic' instead.\n");
    exit(0);
  }
  double d=delta[0];
  if (delta[1]>d) d=delta[1];
  Iso_cuboid domain(x0[0],x0[1],-d,x0[0]+delta[0],x0[1]+delta[1],d);
#endif //NDIMS==3
  Delaunay T(x0,delta,domain);
  T.insert(pos,npart,data.getMassFPtr());
#endif
  */

#ifdef HAVE_CFITS_IO
  if (Opt_angMask) {
    printf ("Computing density ... ");fflush(0);
    T->SetValueToDensity();
    printf ("done.\n");

    int ninside = T->applyAngMask(angMaskfname,std::back_inserter(pts),Opt_maxAngSize);

    for (j=0;j<pts.size();j++)
      {
	if (j<ninside) true_index.push_back(NEWPARTICLE_INSIDE_INDEX);
	else true_index.push_back(NEWPARTICLE_OUTSIDE_INDEX);
	index.push_back(npart+nnew);
	nnew++;
      }

    printf("Adding %ld particles ... ",pts.size());fflush(0);
    if (pts.size())
      T->insert(pts.begin(),pts.end(),NULL,NULL,&(index[0]),&(true_index[0]));
    printf("done.\n");
    
    pts.clear();
    index.clear();
    true_index.clear();
  }
#endif
  
  if (Opt_smoothBoundaries) 
    {
      printf ("Computing density ... ");fflush(0);
      T->SetValueToDensity();
      printf ("done.\n");
	  
      T->GenerateSmoothBoundaries(borderSize,std::back_inserter(pts));
	  
      for (j=0;j<pts.size();j++)
	{
	  true_index.push_back(NEWPARTICLE_OUTSIDE_INDEX);
	  index.push_back(npart+nnew);
	  nnew++;
	}
	  
      printf("Adding %ld particles ... ",pts.size());fflush(0);
      if (pts.size())
	T->insert(pts.begin(),pts.end(),NULL,NULL,&(index[0]),&(true_index[0]));
      printf("done.\n");
	  
      pts.clear();
      index.clear();
      true_index.clear();
    }
            
   
  printf ("Computing density ... ");fflush(0);
  T->SetValueToDensity();
  if (Opt_radialDensity) {
    printf ("correcting ...");
    T->correctNofD(radialDensity_A,radialDensity_Dr,radialDensity_B);
  }
  printf ("done.\n");

  //if (Opt_mangle) T->Mangle(data.getPosFPtr(true),manglefname);

  if (Opt_mask) 
    {     
      T->applyMask(gridType<char>(maskfname).steal());
      /*
      NDfield *f = Load_NDfield(maskfname);
      Convert_NDfield(f,ND_CHAR);
      T->applyMask((char *)(f->val));
      Free_NDfield(&f);
      */
    }

  gettimeofday(&wc_time2, NULL);
  dt  = (double)wc_time2.tv_sec + ((double)wc_time2.tv_usec)*1.E-6;
  dt -= (double)wc_time1.tv_sec + ((double)wc_time1.tv_usec)*1.E-6;
  printf("All done in %.2f s.\n",dt);fflush(0);

  net = T->ToNDnetwork(Opt_periodic,Opt_computeAllFaces);
  delete T;
  /*
  if (Opt_nSmooth) {
    printf("Smoothing %d times ... ",Opt_nSmooth);fflush(0);
    SmoothNDData(net,0,VALUE_TAG, Opt_nSmooth);
    printf("done.\n");
  }
  */
  sprintf(fname,"%s%s",outname,ndnet::IO::getExtension().c_str());
  if (data.getComment().size()<80) strcpy(net->comment,data.getComment().c_str());
  ndnet::IO::save(net,std::string(fname));   
  
  printf("All done.\n");
  
  printf ("\nNetwork was saved as : %s\n",fname);
 
  //printf ("\nResulting complex: %s\n",fname);
  printNDnetStat(net,3);

  printf("\n");
  FreeNDnetwork(&net);
  return return_val;
}
