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

#include "smooth.h"
#include "NDskeleton.h"
#include "mystring.h"
#include "NDnetwork.h"
#include "NDnetwork_tags.h"
#include "distances.h"

#include "NDnet_interface.hxx"
#include "sampledDataInput.hxx"
#include "NDnet_unperiodize.hxx"

#define GLOBAL_DEFINITION
#include "global.h"

void Usage(char *fname)
{
  int i;
  fprintf(stderr,"netconv version %s\n",VERSION_STR);
  fprintf (stderr,"\nUsage: %s <NDnet filename> \n",CutName(fname));
  for (i=0;i<(int)strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   [-outName <output filename>] [-outDir <dir>]\n");
  for (i=0;i<(int)strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   [-addField <filename> <field name>] [-smooth <Ntimes>] \n");
  for (i=0;i<(int)strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   [-smoothData <vertex field name> <Ntimes>] \n");
  for (i=0;i<(int)strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   [-noTags] [-toRaDecZ] [-toRaDecDist] \n");
  for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   [-cosmo <Om=%.2f Ol=%.2f Ok=%.2f h=%.2f w=%.2f>]\n",OMEGAM_DEFAULT,OMEGAL_DEFAULT,OMEGAK_DEFAULT,HUBBLE_DEFAULT,W_DEFAULT);
  for (i=0;i<(int)strlen(CutName(fname));i++) fprintf (stderr," ");
  fprintf (stderr,"   [-unperiodize] [-info] [-to <format>] \n");
  printf ("\n");
  fprintf(stderr,"Accepted file formats:\n");
  std::vector<std::string> lst=ndnet::IO::getTypeList(true,false);
  fprintf(stderr,"  Input: ");
  if (lst.size()) fprintf(stderr,"%s",lst[0].c_str());
  for (i=1;i<(int)lst.size();i++) fprintf(stderr,", %s",lst[i].c_str());
  fprintf(stderr,".\n");
  lst=ndnet::IO::getTypeList(false,true);
  fprintf(stderr,"  Ouput: ");
  if (lst.size()) fprintf(stderr,"%s",lst[0].c_str());
  for (i=1;i<(int)lst.size();i++) fprintf(stderr,"%s%s",(i%6)?", ":",\n         ",lst[i].c_str());
  fprintf(stderr,".\n\n");
  
  exit(0);
}

int main(int argc,char **argv)
{
  long i,j;
  char fName[255];
  char outname[255];
  char outdir[255];
  char sklfname[255];
  char tmpname[255];
  
  std::vector<std::string> addFieldFile;
  std::vector<std::string> addFieldName;
  
  std::vector<int> smoothDataN;
  std::vector<std::string> smoothDataFieldName;

  int Opt_addField=0;
  int Opt_smoothData=0;
  int Opt_hasFName=0;
  int Opt_outName=0;
  int Opt_outDir=0;
  int Opt_skelTag=0;
  int Opt_smooth=0;  
  int Opt_info=0;
  int Opt_noTags=0;
  int Opt_toRaDecZ=0;
  int Opt_toRaDecDist=0;
  int Opt_unperiodize=0;
  int Opt_save=0;

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
      else if (!strcmp(argv[i],"-skelTag"))
	{
	  i++;
	  if ((i==argc)||(argv[i][0]=='-'))
	    {
	      printf ("\noption '-skelTag' needs an argument.\n");
	      Usage(argv[0]);
	    }
	  Opt_skelTag=1;
	  strcpy(sklfname,argv[i]);
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
      else if (!strcmp(argv[i],"-unperiodize"))
	{
	  i++;
	  Opt_unperiodize=1;
	}
      else if (!strcmp(argv[i],"-info"))
	{
	  i++;
	  Opt_info=1;
	}
      else if (!strcmp(argv[i],"-noTags"))
	{
	  Opt_noTags=1;
	  i++;
	}
      else if (!strcmp(argv[i],"-addField"))
	{
	  Opt_addField++;
	  i++;
	  if ((argc==i)||(argv[i][0]=='-'))
	    {
	      fprintf (stderr,"I Expect a 2 argument for '-addField'\n");
	      Usage(argv[0]);
	    }
	  addFieldFile.push_back(std::string(argv[i++]));
	  if ((argc==i)||(argv[i][0]=='-'))
	    {
	      fprintf (stderr,"I Expect a 2 argument for '-addField'\n");
	      Usage(argv[0]);
	    }
	  addFieldName.push_back(std::string(argv[i++]));
	}
      else if (!strcmp(argv[i],"-smoothData"))
	{
	  Opt_smoothData++;
	  i++;
	  if ((argc==i)||(argv[i][0]=='-'))
	    {
	      fprintf (stderr,"I Expect a 2 argument for '-smoothData'\n");
	      Usage(argv[0]);
	    }
	  smoothDataFieldName.push_back(std::string(argv[i++]));
	  if ((argc==i)||(argv[i][0]=='-'))
	    {
	      fprintf (stderr,"I Expect a 2 argument for '-smoothData'\n");
	      Usage(argv[0]);
	    }
	  smoothDataN.push_back(atoi(argv[i++]));
	}
      else if (!strcmp(argv[i],"-toRaDecZ"))
	{
	  if (Opt_toRaDecDist)
	    {
	      fprintf(stderr,"ERROR: options '-toRaDecZ' and '-toRaDecDist' are NOT compatible.\n");
	      Usage(argv[0]);
	    }
	  Opt_toRaDecZ=1;
	  i++;
	}
      else if (!strcmp(argv[i],"-toRaDecDist"))
	{
	  if (Opt_toRaDecZ)
	    {
	      fprintf(stderr,"ERROR: options '-toRaDecZ' and '-toRaDecDist' are NOT compatible.\n");
	      Usage(argv[0]);
	    }
	  Opt_toRaDecDist=1;
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
	  if (!ndnet::IO::canSave(Opt_to.back()))
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

  NDnetwork *net=ndnet::IO::load(std::string(fName));

  if (Opt_info) {
    printf("Network statistics:\n");
    printNDnetStat(net,3);
  }

  if (!Opt_outName) { 
    //Opt_save=1;
    strcpy(outname,CutName(fName));
  } else Opt_save=1;

  if (Opt_outDir) {
    Opt_save=1;
    char tmp[255];
    strcpy(tmpname,outdir);
    if (outdir[strlen(outdir)-1]!='/') sprintf(outdir,"%s/",tmpname);
    sprintf(tmp,"%s%s",outdir,outname);
    strcpy(outname,tmp);
  }
  /*
  if (Opt_skelTag)
    {
      Opt_info=0;
      NDskel *skl=Load_NDskel(sklfname);
      sprintf(outname,"%s.tgd",outname);
      
      Free_NDskeleton(&skl);
    }
  */
  if (Opt_smooth) 
    {
      Opt_save=1;
      strcpy(tmpname,outname);
      if (!Opt_noTags)  sprintf(outname,"%s.S%3.3d",tmpname,Opt_smooth);
      printf("Smoothing network geometry %d times ... ",Opt_smooth);fflush(0);
      NDnet_smooth(net, Opt_smooth);
      printf("done.\n");
    }

  if (Opt_addField)
    {
      Opt_save=1;
      for (i=0;i<Opt_addField;i++)
	{
	  NDfield *f=ndfield::IO::load(addFieldFile[i].c_str());
	  if (f->ndims==1)
	    {
	      if (f->nval==net->nvertex)
		{
		  Convert_NDfield(f,ND_DOUBLE);
		  addNDDataArr(net,0,addFieldName[i].c_str(),(double**)&f->val);
		}
	      
	      for (j=0;j<=net->ndims;j++)
		if ((net->nfaces[j])&&(f->nval == net->nfaces[j])) 
		  {
		    Convert_NDfield(f,ND_DOUBLE);
		    addNDDataArr(net,j,addFieldName[i].c_str(),(double**)&f->val);
		  }
	    }
	  else if (f->ndims == net->ndims)
	    {
	      if (f->fdims_index)
		{
		  fprintf(stderr,"ERROR: file '%s' is of type coords (fdims_index=1).\n",addFieldFile[i].c_str());
		  fprintf(stderr,"       I don t know how to tag a network with that ...\n");
		}
	      AddNDDataTagFromGrid(net,f,addFieldName[i].c_str());
	    }
	  Free_NDfield(&f);
	}
      //fprintf(stderr,"\n SORRY: NOT IMPLEMENTED YET :).\n");
    }

    if (Opt_smoothData)
    {
      Opt_save=1;
      int old_verbose=verbose;
      verbose=1;
      strcpy(tmpname,outname);
      if (!Opt_noTags)  sprintf(outname,"%s.SD",tmpname); 
      
      for (i=0;i<Opt_smoothData;i++)
	{
	  printf("Smoothing '%s' over network %d times ... ",smoothDataFieldName[i].c_str(),smoothDataN[i]);
	  fflush(0);
	  SmoothNDData(net,0,smoothDataFieldName[i].c_str(),smoothDataN[i]);
	  printf("done.\n");
	}
      verbose=old_verbose;
    }

    if (Opt_unperiodize)
      {
	Opt_save=1;
	strcpy(tmpname,outname);
	if (!Opt_noTags)  sprintf(outname,"%s.unPer",tmpname); 
	printf("Unperiodizing the network ... ");fflush(0);
	net=NetworkUnperiodizer(net).get();
	printf("done.\n");
      }

    if (Opt_toRaDecZ)
    {
      Opt_save=1;
      strcpy(tmpname,outname);
      if (!Opt_noTags) sprintf(outname,"%s.RaDecZ",tmpname);
      if (net->ndims!=3)
	{
	  fprintf(stderr,"ERROR converting to (Ra,Dec,z) : ndims = %d\n",net->ndims);
	  exit(0);
	}
      for (i=0;i<net->nvertex;i++)
	{
	  sampledDataInput::cartesian2spherical<float>(&net->v_coord[i*net->ndims],&net->v_coord[i*net->ndims],net->ndims);
	  net->v_coord[i*net->ndims+2]=survey_dist2red(net->v_coord[i*net->ndims+2]);
	}
    }
      
  if (Opt_toRaDecDist)
    {
      Opt_save=1;
      strcpy(tmpname,outname);
      if (!Opt_noTags) sprintf(outname,"%s.RaDecDist",tmpname);
      for (i=0;i<net->nvertex;i++)
	sampledDataInput::cartesian2spherical<float>(&net->v_coord[i*net->ndims],&net->v_coord[i*net->ndims],net->ndims);
    }
  
  if (Opt_to.size())
    {
      for (i=0;i<(long)Opt_to.size();i++)
	ndnet::IO::save(net,std::string(outname)+ndnet::IO::getExtension(Opt_to[i]),Opt_to[i]);
    }
  else if (Opt_save)
    ndnet::IO::save(net,std::string(outname)+ndnet::IO::getExtension());
    
}
