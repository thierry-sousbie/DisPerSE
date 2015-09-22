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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "mystring.h"

#include <QApplication>
#include "PDV_interface.hxx"

#define GLOBAL_DEFINITION
#include "global.h"


void Usage(char *fname)
{
  fprintf(stderr,"pdview version %s\n",VERSION_STR);
  fprintf (stderr,"\nUsage: %s <fname.ppairs.NDnet> [-dtfe] [-nsig <val>] [-cut <val>]\n",CutName(fname));
  fprintf (stderr,"\nnb: Only binary NDnet type files are accepted (option -ppairs in 'mse').\n");
  fprintf (stderr,"Options:\n");
  fprintf (stderr,"   * '-dtfe' : set default parameters for DTFE type networks.\n");
  fprintf (stderr,"               Use that for delaunay tesselations produced with 'delaunay_nD'.\n");
  fprintf (stderr,"   * '-nsig' : set default persistence ratio threshold and enforce '-dtfe'.\n");
  fprintf (stderr,"               Use that for delaunay tesselations produced with 'delaunay_nD'.\n");
  fprintf (stderr,"   * '-cut' : set default persistence threshold.\n");
  
  exit(0);
}

int main(int argc, char **argv)
{ 
  long i;
  char fname[255]; 
  bool Opt_fname=false;
  bool Opt_useRatio=false;
  bool Opt_nsig=false;
  bool Opt_cut=false;
  double nsig=0;
  double cut=0;

  if (argc==1) Usage(argv[0]);
  for (i=1;i<argc;)
    {
      if (argv[i][0]!='-')
	{
	  strcpy(fname,argv[i]);
	  Opt_fname=true;
	  i++;
	}
      else if (!strcmp(argv[i],"-dtfe"))
	{
	  Opt_useRatio=true;
	  i++;
	}   
      else if (!strcmp(argv[i],"-nsig"))
	{
	  i++;
	  if ((i==argc)||(argv[i][0]=='-'))
	    {
	      printf ("\noption '-nsig' needs an argument.\n");
	      Usage(argv[0]);
	    }
	  Opt_nsig=true;
	  nsig= atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-cut"))
	{
	  i++;
	  if ((i==argc)||(argv[i][0]=='-'))
	    {
	      printf ("\noption '-nsig' needs an argument.\n");
	      Usage(argv[0]);
	    }
	  Opt_cut=true;
	  cut= atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-sayHello"))
	{
	  printf("Hi there.\n");
	  exit(0);
	}  
      else Usage(argv[0]);
    }

  if (!Opt_fname)
    {
      fprintf(stderr,"I need a file name !!!!\n");
      Usage(argv[0]);
    }
  
  QApplication app(argc, argv);
  app.setApplicationName("pdview");
  
  PDV_interface pdv(0);  
  pdv.resize(800,600);
  pdv.show();

  if (Opt_nsig) Opt_useRatio=true;
  if (Opt_fname) pdv.loadNetwork(fname,Opt_useRatio);
  
  if (Opt_nsig) pdv.setNSig(nsig);
  else if (Opt_cut) pdv.setCut(cut);

  return app.exec();
}
