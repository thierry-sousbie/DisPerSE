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
//#include <sys/stat.h>

#include <iterator>
#include <fstream>
#include <algorithm>
#include <vector>

#include "fof_io.h"
#include "NDnetwork.h"
#include "mystring.h"
#include "smooth.h"
#include "NDnetwork_tags.h"

#include "NDnet_network.hxx"
#include "NDfield_network.hxx"
#include "Healpix_network.hxx"
#include "SimplicialGrid_network.hxx"
#include "DUMMY_network.hxx"

#include "gridType.hxx"
#include "NDnet_interface.hxx"
#include "PDV_remoteInterface.hxx"

#define GLOBAL_DEFINITION
#include "global.h"

using namespace std;

void Usage(char *fname)
{
  int i;
  fprintf(stderr,"mse version %s\n",VERSION_STR);
  fprintf(stderr,"\nUsage:\n  %s <network filename> [-field <fname>]\n",CutName(fname));
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf(stderr,"   [-outName <fname>] [-noTags] [-outDir <dir>] \n");
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf(stderr,"   [-periodicity <val>]\n");
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf(stderr,"   [-mask <fname.ND>[~]]\n");// [-nonSimplicial]\n"); 
#if defined (USE_OPENMP) || defined (USE_THREADS)
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf(stderr,"   [-nthreads <N=$OMP_NUM_THREADS>]\n");
#endif
  printf("\n");

  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf(stderr,"   [-nsig <n1, n2, ...>] [-cut <l1, l2, ...>]\n");// [-cutR <l1, l2, ...>]\n");
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf(stderr,"   [-interactive [<full/path/to/pdview>]]\n");
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf(stderr,"   [-forceLoops] [-keepLoops]\n");
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf(stderr,"   [-robustness] [-no_robustness]\n");
  printf("\n");

  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf(stderr,"   [-manifolds] [-interArcsGeom] [-no_arcsGeom]\n");
  printf("\n");

  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf(stderr,"   [-ppairs] [-ppairs_ASCII] \n"); 
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf(stderr,"   [-upSkl] [-downSkl] [-interSkl]\n");   
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf(stderr,"   [-dumpArcs <CUID>]\n");   
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf(stderr,"   [-dumpManifolds [<JEP0123456789ad>]] \n");
  //for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  //fprintf(stderr,"   [-dumpConflicts]\n");
  printf("\n"); 

  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf(stderr,"   [-compactify <type=natural>]\n"); 
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf(stderr,"   [-vertexAsMinima] [-descendingFiltration] \n");  
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," "); 
  fprintf(stderr,"   [-no_saveMSC] [-loadMSC <fname>]\n");
  printf("\n");
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");  
  fprintf(stderr,"   [-no_gFilter]\n");
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf(stderr,"   [-diagram] [-smooth <Ntimes=0>] \n");  
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf(stderr,"   [-debug] [-dispSLevel]\n");
 

  fprintf (stderr,"NOTE: the function whose complex is computed may be specified with '-field'\nor tagged within the input network as 'field_value'.\n");
   printf ("\n\n");
  //fprintf (stderr,"   * '-simplicial' forces decomposition into a simplicial complex.\n");
  fprintf (stderr,"   * '-periodicity' sets periodic boundary conditions (PBC) when appropriate.\n");
  fprintf (stderr,"        parameter is a serie of 0/1 enabling periodic boundary conditions\n");
  fprintf (stderr,"        along the corresponding direction.\n");
  fprintf (stderr,"        Example: '-periodicity 0' sets non periodic boundaries\n");
  fprintf (stderr,"                 '-periodicity 101' sets PBC along dims 1 and 3 (x and z)\n");
  fprintf (stderr,"   * '-interactive' uses 'pdview' to interactively select persistence threshold.\n");
  fprintf (stderr,"   * '-nsig' sets persistence ratio threshold for DTFE type densities.\n");
  fprintf (stderr,"        This sets a limit on persistence ratios in terms of 'n-sigmas'\n"); 
  fprintf (stderr,"        Use this for DTFE densities in delaunay tesselations.\n"); 
  fprintf (stderr,"   * '-cut' sets persistence threshold.\n");
  fprintf (stderr,"        There will be distinct outputs for each given threshold\n"); 
  fprintf (stderr,"   * '-field' may be used to set/replace the function value\n");
  fprintf (stderr,"   * '-compactify': can be 'natural' (default), 'sphere' or 'torus'\n");
  fprintf (stderr,"        * 'natural' is the default value, and almost always the best choice.\n");
  fprintf (stderr,"        * 'torus' creates reflective (i.e periodic) boundaries.\n");
  fprintf (stderr,"        * 'sphere' links boundaries to a cell at infinity.\n");  
  fprintf (stderr,"   * '-mask': Mask must be and array.\n");
  fprintf (stderr,"        By default, non-null values correspond to masked pixels.\n");
  fprintf (stderr,"        Adding a trailing '~' reverses the convention.\n");
  //fprintf (stderr,"   * '-upSkl'/'-downSkl'/'-interSkl' \n");
  fprintf (stderr,"   * '-dumpArcs' saves arcs geometry (may be called several times).\n");
  fprintf (stderr,"        U(p): arcs leading to maxima.\n");
  fprintf (stderr,"        D(own): arcs leading to minima.\n");
  fprintf (stderr,"        I(nter): other arcs.\n");
  fprintf (stderr,"        C(onnect): keeps at least the connectivity information for all arcs.\n");
  fprintf (stderr,"        Ex: CU dumps geometry of arcs from maxima, only connectivity for others.\n");
  fprintf (stderr,"   * '-dumpManifolds' saves manifolds geometry (may be called several times).\n");
  fprintf (stderr,"        J(oin): join all the the manifolds in one single file\n");
  fprintf (stderr,"        E(xtended): compute extended manifolds\n");
//fprintf (stderr,"        B: add the boundary to the  manifolds\n");
  fprintf (stderr,"        P(reserve): do not merge infinitely close submanifolds.\n");
  fprintf (stderr,"        D(ual): compute dual cells geometry when appropriate.\n");
  fprintf (stderr,"        0123456789: specifies the critical index (for opt. a/d)\n");
  fprintf (stderr,"        a/d: compute the Ascending and/or Descending manifold.\n");
  fprintf (stderr,"        Ex: JD0a2ad3d dumps to a single file the ascending manifolds\n");
  fprintf (stderr,"         of critical index 0 and 2, and the descending ones for critical\n");
  fprintf (stderr,"         index 2 and 3. Cells are replaced by their dual where appropriate.\n");
  
  

//for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
//fprintf(stderr,"    [-groups <grps filename>] [-simple] [-skipSanitize] [-explicit_overwrite]\n");  

  printf ("\n");
  exit(0);
}


int main(int argc, char **argv)
{
  group_list *groups=NULL;
  NDskel *skl=NULL;
  int ndims=0;
    
  int i,j;

  char filename[255];
  char grps_fname[255];
  char outname[255];
  char outdir[255];
  char tmpname[255];
  char MSC_fname[255];
  char scalarField_fname[255];
  char maskfname[255];
  char PDVfname[255];
    
  std::vector<double> Opt_NSig;
  std::vector<double> Opt_Cut;
  std::vector<double> Opt_CutR;

  int Opt_FNameTags=1;
  int Opt_skl=1;    
  int Opt_up=0;
  int Opt_down=0;
  int Opt_inter=0;
  int Opt_ppairs=0;
  int Opt_ppairsASCII=0;
  int Opt_manifolds=0;
  int Opt_noCycles=0;
  int Opt_robustness=-1;
  int Opt_arcsGeom=1;
  int Opt_hasgrps=0;
  int Opt_skipLoops=0;
  int Opt_hasnetwork=0;
  int Opt_smooth=0;
  int Opt_explicit_overwrite=0;
  int Opt_simpleObjects=0;
  int Opt_sanitize=1;
  int Opt_gFilter=1;
  int Opt_querySigma=0;
  int Opt_simplicial=1;
  int Opt_periodic=0xffffffff;
  int Opt_interactive=0;

  int Opt_dumpManifolds=0;
  long Opt_whichManifolds[50];
  char Opt_whichManifoldsStr[50][255];

  int Opt_dumpArcs=0;
  long Opt_whichArcs[50];
  char Opt_whichArcsStr[50][255];
  int Opt_arcsGeomFlags=(1<<0)|(1<<2); // default to up and down
    
  int Opt_dumpConflicts=0;
  int Opt_withBoundary=0;
  int Opt_scalarField=0;
    
  int Opt_loadMSC=0;
  int Opt_saveMSC=1;
  int Opt_outFName=0;
  int Opt_outDir=0;
  int Opt_mask=0;
  int Opt_diag=0;
  int Opt_vertexAsMinima=0;
  int Opt_descendingFiltration=0;
  int Opt_ompNThreads=0;
  MSComplex::BoundaryType Opt_btype=MSComplex::BType_Natural;
  double *groupid;
  bool netValueFromDTFE=false;

  PDV_remote PDV(false);

    debug_dump=0;
    verbose=2;
    if (argc==1) Usage(argv[0]);
    for (i=1;i<argc;)
    {
      
      if (argv[i][0]!='-')
	{
	  Opt_hasnetwork=1;
	  strcpy(filename,argv[i]);
	  i++;
	}
      
      else if (i==1)
	{
	  if (!strcmp(argv[i],"-dispSLevel"))
	    {
	      Opt_querySigma=1;
	      i++;	    
	    }
	  else
	    {
	      fprintf (stderr,"First argument must be a filename.\n");
	      Usage(argv[0]);
	    }
	}
	else if (!strcmp(argv[i],"-groups"))
	{
	    Opt_hasgrps=1;
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a filename as argument of '-groups'\n");
		Usage(argv[0]);
	    }
	    strcpy(grps_fname,argv[i]);
	    Opt_simpleObjects=1;
	    i++;
	}
	else if (!strcmp(argv[i],"-field"))
	  {
	    i++;
	    Opt_scalarField=1;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a filename as argument of '-field'\n");
		Usage(argv[0]);
	    }
	    strcpy(scalarField_fname,argv[i]);
	    
	    i++;
	  }
	else if (!strcmp(argv[i],"-periodicity"))
	{
	  i++;
	  if ((i==argc)||(argv[i][0]=='-'))
	    {
	      printf ("\noption '-periodicity' needs an argument.\n");
	      Usage(argv[0]);
	    }
	  Opt_periodic=0;
	  for (j=0;j<strlen(argv[i]);j++)
	    {
	      if (argv[i][j]-'0')
		Opt_periodic |= (1<<j);
	    }
	  
	  i++;
	}
	else if (!strcmp(argv[i],"-smooth"))
	{
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a number as argument of '-smooth'\n");
		Usage(argv[0]);
	    }
	    Opt_smooth = atoi(argv[i]);
	    i++;
	}
	else if (!strcmp(argv[i],"-nsig"))
	  {
	    bool loopPrms;	  
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	      {
		fprintf (stderr,"I Expect a number as argument of '-nsig'\n");
		Usage(argv[0]);
	      }
	    do {
	      loopPrms=false;
	      Opt_NSig.push_back(atof(argv[i]));
	      i++;
	      if (i<argc)
		{
		  if (argv[i][0]!='-') loopPrms=true;
		}
	    } while(loopPrms);
	}
	else if (!strcmp(argv[i],"-cut"))
	{
	  bool loopPrms;	  
	  i++;
	  if ((argc==i)||(argv[i][0]=='-'))
	    {
	      fprintf (stderr,"I Expect a value as argument of '-cut'\n");
	      Usage(argv[0]);
	    }
	  do {
	    loopPrms=false;
	    Opt_Cut.push_back(atof(argv[i]));
	    i++;
	    if (i<argc)
	      {
		if (argv[i][0]!='-') loopPrms=true;
	      }
	  } while(loopPrms);
	}
	else if (!strcmp(argv[i],"-cutR"))
	{
	  bool loopPrms;	  
	  i++;
	  if ((argc==i)||(argv[i][0]=='-'))
	    {
	      fprintf (stderr,"I Expect a value as argument of '-cutR'\n");
	      Usage(argv[0]);
	    }
	  do {
	    loopPrms=false;
	    Opt_CutR.push_back(atof(argv[i]));
	    i++;
	    if (i<argc)
	      {
		if (argv[i][0]!='-') loopPrms=true;
	      }
	  } while(loopPrms);
	}
	else if (!strcmp(argv[i],"-nocut"))
	  {
	    i++;
	    //Opt_Cut.push_back(-1);
	  }
	else if (!strcmp(argv[i],"-explicit_overwrite"))
	{
	    Opt_explicit_overwrite=1;
	    i++;
	}
	else if (!strcmp(argv[i],"-simple"))
	{
	    Opt_simpleObjects=1;
	    i++;
	}
	else if (!strcmp(argv[i],"-noTags"))
	  {
	    Opt_FNameTags=0;
	    i++;
	  }
	else if (!strcmp(argv[i],"-withBoundary"))
	{
	  Opt_withBoundary=1;
	  i++;
	}
	else if (!strcmp(argv[i],"-dumpManifolds"))
	{
	  int cur=Opt_dumpManifolds;
	  Opt_whichManifolds[cur]=0;
	  Opt_dumpManifolds+=1;
	    i++;
	    if ((argc!=i)&&(argv[i][0]!='-')) 
	      {
		int curD=-1;
		for (j=0;j<strlen(argv[i]);j++)
		  {
		    if (isdigit(argv[i][j])) curD=argv[i][j]-'0';
		    if ((argv[i][j]=='a')||(argv[i][j]=='d')) 
		      {
			if (curD<0) {
			  fprintf (stderr,"Error: Invalid argument for '-dumpManifolds'.\n");
			  Usage(argv[0]);
			}
			if (argv[i][j]=='a') Opt_whichManifolds[cur]|=((long)1<<(2*curD));
			if (argv[i][j]=='d') Opt_whichManifolds[cur]|=((long)1<<(2*curD+1));
		      }
		    if (argv[i][j]=='P') Opt_whichManifolds[cur]|=((long)1<<(8*sizeof(Opt_whichManifolds[cur])-1));
		    if (argv[i][j]=='E') Opt_whichManifolds[cur]|=((long)1<<(8*sizeof(Opt_whichManifolds[cur])-2));
		    if (argv[i][j]=='J') Opt_whichManifolds[cur]|=((long)1<<(8*sizeof(Opt_whichManifolds[cur])-3));
		    if (argv[i][j]=='D') Opt_whichManifolds[cur]|=((long)1<<(8*sizeof(Opt_whichManifolds[cur])-4));
		  }
		strcpy(Opt_whichManifoldsStr[cur],argv[i]);
		i++;
	      }
	}
	else if (!strcmp(argv[i],"-dumpArcs"))
	{
	  int cur=Opt_dumpArcs;
	  Opt_whichArcs[cur]=0;
	  Opt_dumpArcs+=1;
	  i++;
	  if ((argc!=i)&&(argv[i][0]!='-')) 
	    {
	      for (j=0;j<strlen(argv[i]);j++)
		{		  
		  if ((argv[i][j]=='u')||(argv[i][j]=='U')) Opt_whichArcs[cur]|=1<<0;
		  if ((argv[i][j]=='i')||(argv[i][j]=='I')) Opt_whichArcs[cur]|=1<<1;
		  if ((argv[i][j]=='d')||(argv[i][j]=='D')) Opt_whichArcs[cur]|=1<<2;
		  if ((argv[i][j]=='c')||(argv[i][j]=='C')) Opt_whichArcs[cur]|=1<<3;
		}
		strcpy(Opt_whichArcsStr[cur],argv[i]);
		i++;
	    }
	}
	else if (!strcmp(argv[i],"-nonSimplicial"))
	{
	    Opt_simplicial=0;
	    i++;
	}
	else if (!strcmp(argv[i],"-dumpConflicts"))
	{
	    Opt_dumpConflicts=1;
	    i++;
	}
	else if (!strcmp(argv[i],"-skipSanitize"))
	{
	    Opt_sanitize=0;
	    i++;
	}
	else if (!strcmp(argv[i],"-manifolds"))
	{
	    Opt_manifolds=1;
	    i++;
	}
        else if (!strcmp(argv[i],"-no_cycles"))
	{
	    Opt_noCycles=1;
	    i++;
	}
	else if (!strcmp(argv[i],"-forceLoops"))
	{
	    Opt_skipLoops=0;
	    i++;
	}
        else if (!strcmp(argv[i],"-keepLoops"))
	{
	    Opt_skipLoops=1;
	    i++;
	}
        else if (!strcmp(argv[i],"-robustness"))
	{
	    Opt_robustness=1;
	    i++;
	}
        else if (!strcmp(argv[i],"-no_robustness"))
	{
	    Opt_robustness=0;
	    i++;
	}      
	else if (!strcmp(argv[i],"-no_gFilter"))
	{
	    Opt_gFilter=0;
	    i++;
	}
	else if (!strcmp(argv[i],"-no_arcsGeom"))
	{
	    Opt_arcsGeom=0;
	    i++;
	}
	else if (!strcmp(argv[i],"-upSkl"))
	{
	    Opt_up=1;
	    i++;
	}
	else if (!strcmp(argv[i],"-downSkl"))
	{
	    Opt_down=1;
	    i++;
	}
	else if (!strcmp(argv[i],"-interSkl"))
	{
	    Opt_inter=1;
	    i++;
	}
	else if (!strcmp(argv[i],"-interArcsGeom"))
	{
	  Opt_arcsGeomFlags|=(1<<1);
	  i++;
	}
	else if (!strcmp(argv[i],"-ppairs"))
	{
	    Opt_ppairs=1;
	    i++;
	}
	else if (!strcmp(argv[i],"-ppairs_ASCII"))
	{
	    Opt_ppairsASCII=1;
	    i++;
	}
	
	else if (!strcmp(argv[i],"-no_sklDump"))
	{
	    Opt_skl=0;
	    i++;
	}	
	else if (!strcmp(argv[i],"-vertexAsMinima"))
	{
	    Opt_vertexAsMinima=1;
	    i++;
	}
	else if (!strcmp(argv[i],"-descendingFiltration"))
	{
	    Opt_descendingFiltration=1;
	    i++;
	}
	else if (!strcmp(argv[i],"-dispSLevel"))
	  {
	    Opt_querySigma=1;
	    i++;	    
	  }
	else if (!strcmp(argv[i],"-loadMSC"))
	  {
	    //Opt_saveMSC=0;
	    Opt_loadMSC=1;
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a filename as argument of '-loadMSC'\n");
		Usage(argv[0]);
	    }
	    strcpy(MSC_fname,argv[i]);
	    i++;
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
	else if (!strcmp(argv[i],"-interactive"))
	  {
	    Opt_interactive=1;
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	      {
		//fprintf (stderr,"I Expect a filename as argument of '-outDir'\n");
		//Usage(argv[0]);
	      }
	    else
	      {
		Opt_interactive=2;
		strcpy(PDVfname,argv[i]);
		i++;
		if ((argc==i)||(argv[i][0]=='-'))
		  {
		    fprintf (stderr,"Only one optional argument is accepted for '-interactive'\n");
		    Usage(argv[0]);
		  }
	      }
	  }
	else if (!strcmp(argv[i],"-no_saveMSC"))
	  {
	    Opt_saveMSC=0;
	    i++;	    
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
	    if (maskfname[strlen(maskfname)-1]=='~')
	      {
		maskfname[strlen(maskfname)-1]='\0';
		Opt_mask=-1;
	      }
	    i++;
	}
	else if (!strcmp(argv[i],"-compactify"))
	{
	  //Opt_compactify=1;	    
	  i++;
	  if ((argc==i)||(argv[i][0]=='-'))
	    {
	      fprintf (stderr,"I Expect an argument to option '-compactify'\n");
	      Usage(argv[0]);
	    }
	  if (strcmp(argv[i],"sphere")==0)
	    Opt_btype = MSComplex::BType_Sphere;
	  else if (strcmp(argv[i],"torus")==0)
	    Opt_btype = MSComplex::BType_Torus;
	  else if (strcmp(argv[i],"natural")==0)
	    Opt_btype = MSComplex::BType_Natural;
	  else 
	    {
	      fprintf (stderr,"\nERROR: %s is not a valid argument for '-compactify'.\n",argv[i]);
	      Usage(argv[0]);
	    }
	  
	  i++;
	}
	else if (!strcmp(argv[i],"-diagram"))
	  {
	    Opt_diag=1;
	    i++;
	  }
	else if (!strcmp(argv[i],"-debug"))
	  {
	    debug_dump=1;
	    i++;
	  }
	else if (!strcmp(argv[i],"-nthreads"))
	{
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a number as argument of '-ompNThreads'\n");
		Usage(argv[0]);
	    }
	    Opt_ompNThreads = atoi(argv[i]);
	    i++;
	}
	else 
	{
	    printf ("\nWhat is %s ???\n",argv[i]);
	    Usage(argv[0]);
	}
    }

    glob_num_threads=1;
    glob_num_omp_threads=1;
   
    if (Opt_ompNThreads) {
#ifdef USE_OPENMP
      omp_set_num_threads(Opt_ompNThreads);
#pragma omp parallel
      glob_num_omp_threads=omp_get_num_threads();

#else
      fprintf(stderr,"WARNING: option '-nThreads' ignored.\n");
      fprintf(stderr,"   Add command 'set (USE_THREADS 1)' in CMakeList.txt and recompile\n");
      fprintf(stderr,"   in order to enable openMP/pthreads support.\n");
#endif
    }
    else 
      {
#ifdef USE_OPENMP
	#pragma omp parallel
	glob_num_omp_threads=omp_get_num_threads();
#endif
      }
#ifdef USE_THREADS
#ifdef USE_OPENMP
    glob_num_threads=glob_num_omp_threads;
#else
    glob_num_threads=Opt_ompNThreads;
#endif
#endif

    printf ("\n****** MSE v%s (%d bits) ******\n",VERSION_STR,(sizeof(long)==8)?64:32);

    if (Opt_interactive)
      {
	char path_[255];
	char *path=&path_[0];

	printf ("Looking for 'pdview' ... ");fflush(0);

	getPath(argv[0],&path);	
	if (Opt_interactive==1) 
	  PDV.findBinary(std::string(path));
	else 
	  {
	    if (!PDV.findBinary(std::string(PDVfname))) 
	      PDV.findBinary(std::string(path));
	  }
	
	if (!PDV.binaryFound()) 
	  {
	    fprintf(stderr,"\n  ERROR: could not open binary for 'pdview' (-interactive).\n");
	    fprintf(stderr,"  Try specifying the location: '-interactive <full/path/to/pdview>'\n");
	    exit(-1);
	  }
	else printf ("found.\n");

      }
    
    if (Opt_Cut.size()!=Opt_NSig.size())
      {
	if ((Opt_Cut.size()>1)&&(Opt_NSig.size()>1))
	  {
	    fprintf(stderr,"\n FATAL ERROR: you gave an inconsistent numberof of '-cut' / '-nsig'\n");
	    fprintf(stderr,"  There should be as many of each, or one of them  must appear 0 or 1 time.\n");
	    exit(0);
	  }
	else
	  {
	    //if (Opt_Cut.size()==0) Opt_Cut.push_back(0);
	    if (Opt_Cut.size()==1)
	      for (i=1;i<Opt_NSig.size();i++) Opt_Cut.push_back(Opt_Cut[0]);

	    //if (Opt_NSig.size()==0) Opt_NSig.push_back(0);
	    if (Opt_NSig.size()==1)
	      for (i=1;i<Opt_Cut.size();i++) Opt_NSig.push_back(Opt_NSig[0]);
	  }
      }
    else if ((Opt_NSig.size()==0)&&(Opt_Cut.size()==0))
      {
	//Opt_NSig.push_back(0);
	//Opt_Cut.push_back(0);
      }

    if ((Opt_CutR.size()!=Opt_NSig.size())||(Opt_Cut.size()!=Opt_NSig.size()))
      {
	if ((Opt_CutR.size()>1)+(Opt_NSig.size()>1)+(Opt_Cut.size()>1) >1)
	  {
	    fprintf(stderr,"\n FATAL ERROR: you gave an inconsistent numberof of '-cut' / '-nsig' / 'cutR'\n");
	    fprintf(stderr,"  There should be as many of each, or one of them  must appear 0 or 1 time.\n");
	    exit(0);
	  }
	else
	  {
	    long nmax=std::max(Opt_CutR.size(),Opt_NSig.size());
	    nmax=std::max(nmax,(long)Opt_Cut.size());

	    if (Opt_CutR.size()==0) Opt_CutR.push_back(0);
	    if (Opt_CutR.size()==1)
	      for (i=1;i<nmax;i++) Opt_CutR.push_back(Opt_CutR[0]);
	    
	    if (Opt_NSig.size()==0) Opt_NSig.push_back(0);
	    if (Opt_NSig.size()==1)
	      for (i=1;i<nmax;i++) Opt_NSig.push_back(Opt_NSig[0]);

	    if (Opt_Cut.size()==0) Opt_Cut.push_back(0);
	    if (Opt_Cut.size()==1)
	      for (i=1;i<nmax;i++) Opt_Cut.push_back(Opt_Cut[0]);
	  }
      }
   
    if (!Opt_outFName) strcpy(outname,CutName(filename));    

    if (Opt_outDir) {
      char tmp[255];
      if (outdir[strlen(outdir)-1]!='/') sprintf(outdir,"%s/",outdir);
      sprintf(tmp,"%s%s",outdir,outname);
      strcpy(outname,tmp);
    }

    if (Opt_querySigma)
      {
	double perRatio[50];
	for (i=0;i<Opt_NSig.size();i++)
	  {
	    persistenceRatioFromSigma(Opt_NSig[i],3,perRatio);
	    printf("\n%10.10e\n%10.10e\n%10.10e\n",perRatio[0],perRatio[1],perRatio[2]);
	  }
	exit(0);
      }

    if ((!Opt_arcsGeom)&&(Opt_skl))
    {
	printf ("warning: Skeleton cannot be computed without arc geometry.\n");
	Opt_skl=0;
    }

    if (Opt_smooth<0) Opt_smooth=-Opt_smooth;
    
    if ((!Opt_hasnetwork))//||(!Opt_hasgrps))
    {
	fprintf (stderr,"I need a network to work with ...\n");
	Usage(argv[0]);
    }
  
    MSComplex::NetworkType *data;

    if (NDnet_network<MSComplex::cellType>::checkFileType(filename))
    {
      NDnetwork *net = ndnet::IO::load(std::string(filename));
      if (!strcmp(net->comment,DELAUNAY_TESSELATION_TAG)) 
	netValueFromDTFE=true;
	
	// tag particles with their group index 
        //printf ("*********************************\n");
	//printf ("******* Preparing network *******\n");
	if (Opt_hasgrps) 
	{
	    groups = Load_FOF(grps_fname);
	    groupid = Groups2Id_d(groups);
	    if (groups->NPart != net->nvertex) {
	      groupid = (double *)realloc(groupid,net->nvertex*sizeof(double));
	      memset(&groupid[groups->NPart],0,sizeof(double) * (net->nvertex-groups->NPart));
	    }
	    addNDDataArr(net,0,GROUP_ID_PEAK_TAG,&groupid);
	    freeFOF(&groups);
	}
	
	if (Opt_scalarField)
	  {
	    double *field_val = gridType<double>(scalarField_fname).steal();
	    netValueFromDTFE=false;
	    data = new NDnet_network<MSComplex::cellType>(net,VALUE_TAG,GROUP_ID_PEAK_TAG,&field_val,true);
	  }
	else 
	  data = new NDnet_network<MSComplex::cellType>(net,VALUE_TAG,GROUP_ID_PEAK_TAG,NULL,true);
	
    }
    else if (Healpix_network<MSComplex::cellType>::checkFileType(filename,Opt_simplicial))
      data = Healpix_network<MSComplex::cellType>::Load(filename);
    else if (NDfield_network<MSComplex::cellType>::checkFileType(filename,Opt_simplicial))
      data = NDfield_network<MSComplex::cellType>::Load(filename);
    else if (SimplicialGrid_network<MSComplex::cellType>::checkFileType(filename,Opt_simplicial))
      data = SimplicialGrid_network<MSComplex::cellType>::Load(filename); 
    else if (Dummy_network<MSComplex::cellType>::checkFileType(filename))
      data = Dummy_network<MSComplex::cellType>::Load(filename); 
    else
      {
	fprintf (stderr,"ERROR: the format of file '%s' is unknown.\n",filename);
	return (-1);
      }
    printf ("*********************************\n");
    data->setPeriodicity(Opt_periodic);
    ndims=data->getNDims();
    
    if (Opt_mask) 
      { 
	std::vector<char> mymask=gridType<char>(maskfname).toVector();
	
	if (Opt_mask>0) 
	  data->setMask(mymask,false);
	else 
	  data->setMask(mymask,true);
      }

    MSComplex msc(data,true);
    
    if (Opt_up) 
      {
	strcpy(Opt_whichArcsStr[Opt_dumpArcs],"up");
	Opt_whichArcs[Opt_dumpArcs++]=(1<<0);
      }
    if (Opt_inter)
      {
	strcpy(Opt_whichArcsStr[Opt_dumpArcs],"inter");
	Opt_whichArcs[Opt_dumpArcs++]=(1<<1);
      }
    if (Opt_down)
      {
	strcpy(Opt_whichArcsStr[Opt_dumpArcs],"down");
	Opt_whichArcs[Opt_dumpArcs++]=(1<<2);
      }
    
    if (Opt_arcsGeom)
      {
	for (i=0;i<Opt_dumpArcs;i++)
	  {
	    if (Opt_whichArcs[i]&(1<<0)) Opt_arcsGeomFlags|=(1<<0);
	    if (Opt_whichArcs[i]&(1<<1)) Opt_arcsGeomFlags|=(1<<1);
	    if (Opt_whichArcs[i]&(1<<2)) Opt_arcsGeomFlags|=(1<<2);
	  }
      }
    else Opt_arcsGeomFlags=0;

    if (Opt_dumpManifolds)
      {
	if (!Opt_manifolds)
	  {
	    printf("INFORMATION: enforcing manifolds computation.\n");
	    Opt_manifolds=1;
	  }
      }
       
    if (Opt_loadMSC) {
      std::ifstream ifile(MSC_fname, ios::binary);
      printf ("Reading MSC from file '%s' ... ",MSC_fname);fflush(0);
      msc.read(ifile);
      ifile.close();
      printf("done.\n");
    }
    else msc.compute(Opt_gFilter,Opt_manifolds,Opt_arcsGeomFlags,Opt_vertexAsMinima,Opt_descendingFiltration);

    data->sendMessage(MSComplex::NetworkType::MSG_FACE_COFACE_FREE);
     
    if (Opt_saveMSC&&(!Opt_loadMSC)) {
      sprintf(tmpname,"%s.MSC",outname);
      std::ofstream ofile(tmpname, ios::binary);
      std::string cmd("cmd line: ");
      for (i=0;i<argc;i++) cmd+=std::string(argv[i])+std::string(" ");
      printf ("Writing MSC to file '%s'... ",tmpname);fflush(0);
      msc.write(ofile,cmd);
      ofile.close();
      printf("done.\n");
    }
    
    //if (!Opt_withBoundary) 
    bool ebs = msc.enforceBoundaryConditions(Opt_btype);

    bool computeCycles;
    if (Opt_robustness<0)
      {
	if (ndims<3) computeCycles=false;
	else computeCycles=false;
      }
    else computeCycles=Opt_robustness;

    bool cpps = msc.computePersistencePairs(computeCycles);
       
    //if (!Opt_gFilter) msc.simplifyComplex(0);
   
    if (Opt_saveMSC&&(ebs||cpps)) {
      sprintf(tmpname,"%s.MSC",outname);
      char mvcmd[255];
      int ret;
      sprintf(mvcmd,"mv -f %s %s.backup",tmpname,tmpname);
      ret=system(mvcmd);
      std::ofstream ofile(tmpname, ios::binary);
      std::string cmd("cmd line: ");
      for (i=0;i<argc;i++) cmd+=std::string(argv[i])+std::string(" ");
      printf ("Writing MSC to file '%s'... ",tmpname);fflush(0);
      msc.write(ofile,cmd);
      ofile.close();
      sprintf(mvcmd,"rm -f %s.backup",tmpname);
      ret=system(mvcmd);
      printf("done.\n");
    }

    printf ("*********** Information *********\n");
    msc.printStats();


    if (Opt_interactive)
      {
	printf ("******* Dumping p. pairs ********\n");
	msc.dumpPPairs(outname);
	printf ("*********************************\n");
      }
  
    if ((Opt_interactive)&&(!Opt_NSig.size()))
      {	
	char cmd[1024];
	printf("Starting pdview ...");fflush(0);
 	Opt_interactive=0;
	sprintf(cmd,"%s.ppairs.NDnet",outname);
	std::pair<double,double> res=PDV.getThreshold(cmd,0,0,netValueFromDTFE);
	
	if ((res.first>0)||(res.second>0))
	  {
	    Opt_NSig.push_back(res.first);	    
	    Opt_Cut.push_back(res.second);
	    Opt_CutR.push_back(0);
	    if (res.first>0) printf("done. (-nsig %g)\n",res.first);
	    else if (res.second>0) printf("done. (-cut %g)\n",res.second);
	  }
	else printf(" done. (no threshold)\n");
      }

    if ((Opt_dumpArcs)&&(Opt_NSig.size()==0))
      {
	Opt_NSig.push_back(0);
	Opt_Cut.push_back(0);
	Opt_CutR.push_back(0);
      }
    else if ((Opt_dumpManifolds)&&(Opt_NSig.size()==0))
      {
	Opt_NSig.push_back(0);
	Opt_Cut.push_back(0);
	Opt_CutR.push_back(0);
      }
    
    for (i=0;i<Opt_NSig.size();i++)
      {
	if (Opt_interactive)
	  {
	    char cmd[1024];
	    printf("Starting pdview ...");fflush(0);
	    Opt_interactive=0;
	    sprintf(cmd,"%s.ppairs.NDnet",outname);
	    std::pair<double,double> res=PDV.getThreshold(cmd,Opt_NSig[i],Opt_Cut[i],netValueFromDTFE);
	    
	    if ((res.first>0)||(res.second>0))
	      {
		Opt_NSig[i]=res.first;	    
		Opt_Cut[i]=res.second;
		Opt_CutR[i]=0;
		if (res.first>0) printf("done. (-nsig %g)\n",res.first);
		else if (res.second>0) printf("done. (-cut %g)\n",res.second);
	      }
	    else printf(" done. (no threshold)\n");
	    	  
	  }
	
	printf ("****** Simplifying complex ******\n");
	msc.simplify(Opt_NSig[i], Opt_Cut[i],  Opt_CutR[i],
		     Opt_simpleObjects, Opt_sanitize, 
		     Opt_skipLoops, Opt_dumpConflicts);

	if (msc.removeBoundaries())
	  {
	    if (Opt_NSig.size()>1)
	      {
		printf("WARNING : I cannot compute several simplifications with this type of boundary conditions.\n");
		printf("          Restart the program for each level of simplification ...\n"); 
		Opt_NSig.resize(1);
		Opt_Cut.resize(1);
	      }
	  }
	printf ("*********** Information *********\n");
	msc.printStats();
	printf ("*********************************\n");
	char outname_cut[1024];

	strcpy(outname_cut,outname);

	if (Opt_FNameTags) {	  
	  strcpy(tmpname,outname_cut);
	  if (Opt_Cut[i]!=0)
	    sprintf(outname_cut,"%s_c%.3g",tmpname,Opt_Cut[i]);
	  if (Opt_CutR[i]!=0)
	    sprintf(outname_cut,"%s_cr%.3g",tmpname,Opt_CutR[i]);
	  if (Opt_NSig[i]!=0)
	    sprintf(outname_cut,"%s_s%.3g",tmpname,Opt_NSig[i]);
	}

	if (Opt_skl)
	  {
	    printf ("******* Dumping skeleton ********\n");
	    for (i=0;i<Opt_dumpArcs;i++)
	      msc.dumpSkl(outname_cut, Opt_whichArcs[i],Opt_whichArcsStr[i], Opt_smooth);
	    printf ("*********************************\n");
	  }
	if (Opt_ppairs)
	  {
	    printf ("******* Dumping p. pairs ********\n");
	    msc.dumpPPairs(outname_cut);
	    printf ("*********************************\n");
	  }
	
	if (Opt_ppairsASCII)
	  {
	    printf ("*** Dumping p. pairs (ASCII) ****\n");
	    msc.dumpPPairs(outname_cut,false);
	    printf ("*********************************\n");
	  }
	
	if (Opt_diag) msc.drawDiagram(outname);
	
	if (Opt_dumpManifolds) 
	  {
	    for (j=0;j<Opt_dumpManifolds;j++)
	      {
		printf ("******* Dumping manifolds *******\n");
		msc.dumpManifolds(outname_cut, Opt_smooth, 
				  Opt_whichManifolds[j],
				  Opt_whichManifoldsStr[j]);
		printf ("*********************************\n");
	      }
	  }
      }
    printf ("*********** ALL DONE ************\n");
    
}
