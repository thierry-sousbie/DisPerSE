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
#include <string>
#include "mystring.h"
#include "distances.h"


#include "NDfield_interface.hxx"

#define GLOBAL_DEFINITION
#include "global.h"

void Usage(char *fname)
{
  int i;
  fprintf(stderr,"fieldconv version %s\n",VERSION_STR);
  fprintf (stderr,"\nUsage: %s <NDfield filename> \n",CutName(fname));
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   [-outName <output filename>] [-outDir <dir>]\n");
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   [-cosmo <Om=%.2f Ol=%.2f Ok=%.2f h=%.2f w=%.2f>]\n",OMEGAM_DEFAULT,OMEGAL_DEFAULT,OMEGAK_DEFAULT,HUBBLE_DEFAULT,W_DEFAULT);
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   [-info] [-to <format>] \n");
  printf ("\n");

  fprintf(stderr,"Accepted file formats:\n");
  std::vector<std::string> lst=ndfield::IO::getTypeList(true,false);
  fprintf(stderr,"  Input: ");
  if (lst.size()) fprintf(stderr,"%s",lst[0].c_str());
  for (i=1;i<lst.size();i++) fprintf(stderr,", %s",lst[i].c_str());
  fprintf(stderr,".\n");
  lst=ndfield::IO::getTypeList(false,true);
  fprintf(stderr,"  Ouput: ");
  if (lst.size()) fprintf(stderr,"%s",lst[0].c_str());
  for (i=1;i<lst.size();i++) fprintf(stderr,"%s%s",(i%6)?", ":",\n         ",lst[i].c_str());
  fprintf(stderr,".\n\n");
  
  exit(0);
}

int main(int argc,char **argv)
{
  long i,j;
  char fName[255];
  char outname[255];
  char outdir[255];
  char fname[255];
  char tmpname[255];
  
  int Opt_hasFName=0;
  int Opt_outName=0;
  int Opt_outDir=0;
  int Opt_info=0;

  std::vector<std::string> Opt_to;

  if (argc==1) Usage(argv[0]);
  for (i=1;i<argc;)
    {
      if (argv[i][0]!='-')
	{
	  strcpy(fName,argv[i]);
	  Opt_hasFName=1;
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
      else if (!strcmp(argv[i],"-info"))
	{
	  i++;
	  Opt_info=1;
	}
      else if (!strcmp(argv[i],"-outName"))
	{
	  i++;
	  if ((i==argc)||(argv[i][0]=='-'))
	    {
	      printf ("\noption '-outName' needs an argument.\n");
	      Usage(argv[0]);
	    }
	  Opt_outName=1;
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
      else if (!strcmp(argv[i],"-to"))
	{
	  i++;
	  if ((argc==i)||(argv[i][0]=='-'))
	    {
	      fprintf (stderr,"I Expect a file type as argument to '-to'\n");
	      Usage(argv[0]);
	    }
	  Opt_to.push_back(std::string(argv[i]));
	  if (!ndfield::IO::canSave(Opt_to.back()))
	    {
	      fprintf (stderr,"Error: format '%s' is unknown or read-only.\n",argv[i]);
	      Usage(argv[0]);
	    }
	  i++;
	}
      else 
	{
	  printf ("\nWhat is %s ???\n",argv[i]);
	  Usage(argv[0]);
	}
    }

  verbose=2;
  bool needSave=false;

  if (!Opt_outName) strcpy(outname,CutName(fName));
  else needSave=true;
  if (Opt_outDir) {
    needSave=true;
    strcpy(tmpname,outdir);
    if (outdir[strlen(outdir)-1]!='/') sprintf(outdir,"%s/",tmpname);
    sprintf(tmpname,"%s%s",outdir,outname);
    strcpy(outname,tmpname);
  }

  NDfield *field=ndfield::IO::load(std::string(fName));
  
  if (Opt_info) {
    printf("field statistics:\n");
    printNDfieldStat(field,3);
  }
  
  if (Opt_to.size())
    {
      for (i=0;i<Opt_to.size();i++)
	ndfield::IO::save(field,std::string(outname)+ndfield::IO::getExtension(Opt_to[i]),Opt_to[i]);
    }
  else if (needSave) 
    ndfield::IO::save(field,std::string(outname)+ndfield::IO::getExtension());

  Free_NDfield(&field);
  return 0;
}
