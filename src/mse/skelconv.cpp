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
#include <limits>

#include "smooth.h"
#include "NDskeleton.h"
#include "NDskel_tags.h"
#include "NDskel_breakdown.h"
//#include "NDskel_assemble.h"
#include "mystring.h"
#include "distances.h"

#include "sampledDataInput.hxx"
#include "NDskel_assemble.hxx"
#include "NDskel_trim.hxx"
#include "NDskel_interface.hxx"
#include "NDnet_interface.hxx"

#define GLOBAL_DEFINITION
#include "global.h"

void Usage(char *fname)
{
  int i;
  fprintf(stderr,"skelconv version %s\n",VERSION_STR);
    fprintf (stderr,"\nUsage: %s <skeleton filename> \n",CutName(fname));
    for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
    fprintf (stderr,"   [-noTags] [-outName <output filename>] [-outDir <dir>]\n");//[-toASCII]
    for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
    fprintf (stderr,"   [-smooth <Ntimes=0>] [-breakdown]\n");
    for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
    fprintf (stderr,"   [-assemble [[<field_name>] <threshold>] <angle>]\n");
    for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
    fprintf (stderr,"   [-trimAbove [<field_name>] <threshold>]\n");
    for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
    fprintf (stderr,"   [-trimBelow [<field_name>] <threshold>] \n");
    for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
    fprintf (stderr,"   [-toRaDecZ] [-toRaDecDist]\n");

    for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
    fprintf (stderr,"   [-rmBoundary] [-rmOutside]\n");
#ifdef HAVE_CFITS_IO
    for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
    fprintf (stderr,"   [-toFITS [<Xres>] [<Yres>] [<Zres>]]\n");
#endif
    for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
    fprintf (stderr,"   [-cosmo <Om=%.2f Ol=%.2f Ok=%.2f h=%.2f w=%.2f>]\n",OMEGAM_DEFAULT,OMEGAL_DEFAULT,OMEGAK_DEFAULT,HUBBLE_DEFAULT,W_DEFAULT);
    for (i=0;i<strlen(CutName(fname));i++) fprintf (stderr," ");
    fprintf (stderr,"   [-info] [-addField <filename> <field_name>] [-to <format>] \n");
    printf ("\n");

#ifdef HAVE_CFITS_IO
    fprintf (stderr,"   * '-toFITS' samples filaments to a FITS image, pixel values correspond to\n");
    fprintf (stderr,"       filament ID. If dimensions are ommited, one pixel is one unit distance.\n");
    printf ("\n");
#endif
    
    fprintf(stderr,"Accepted file formats:\n");
    std::vector<std::string> lst=ndskel::IO::getTypeList(true,false);
    fprintf(stderr,"  Input: ");
    if (lst.size()) fprintf(stderr,"%s",lst[0].c_str());
    for (i=1;i<lst.size();i++) fprintf(stderr,", %s",lst[i].c_str());
    fprintf(stderr,".\n");
    lst=ndskel::IO::getTypeList(false,true);
    fprintf(stderr,"  Ouput: ");
    if (lst.size()) fprintf(stderr,"%s",lst[0].c_str());
    for (i=1;i<lst.size();i++) fprintf(stderr,"%s%s",(i%6)?", ":",\n         ",lst[i].c_str());
    fprintf(stderr,".\n\n");
    //printf("Note: skeleton files can also be loaded/converted as networks with 'netconv'.\n");
    exit(0);
}


int main(int argc,char **argv)
{
  verbose=2;

  long i,j;
  char fName[255];
  char outname[255];
  char outdir[255];
  char tmpname[255];

  std::vector<std::string> addFieldFile;
  std::vector<std::string> addFieldName;
  std::vector<long> FITSRes;

  int Opt_toFITS=0;  
  int Opt_info=0;
  int Opt_addField=0;
  int Opt_hasFName=0;
  int Opt_outName=0;
  int Opt_outDir=0;
  int Opt_toNDskel=1;
  int Opt_smooth=0;
  int Opt_smoothN=0;
  int Opt_toRaDecZ=0;
  int Opt_toRaDecDist=0;
  int Opt_assemble=0;
  int Opt_nTrim=0;
  int Opt_rmBoundary=0;
  int Opt_rmOutside=0;
  int Opt_trim[255];
  double ass_threshold=0.;
  double trim_threshold[255];
  std::string trim_field[255];
  bool trim_above[255];
  std::string ass_field=std::string(VALUE_TAG);
  double ass_angle=45;  

  int Opt_breakdown=0;
  int Opt_noTags=0;

  int count=1;

  std::vector<std::string> Opt_to;
  
  for (i=0;i<255;i++) 
    {
      trim_field[i]=std::string(VALUE_TAG);
      trim_threshold[i]=0;
      trim_above[i]=false;
      Opt_trim[i]=0;
    }

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
      else if (!strcmp(argv[i],"-noTags"))
	{
	  Opt_noTags=1;
	  i++;
	}
#ifdef HAVE_CFITS_IO
      else if (!strcmp(argv[i],"-toFITS"))
	{
	    Opt_toFITS=1;
	    i++;
	    while ((argc!=i)&&(argv[i][0]!='-'))
	    {
	      FITSRes.push_back(atof(argv[i]));
	      i++;
	    }
	}
#endif     
      else if (!strcmp(argv[i],"-assemble"))
	{
	  Opt_assemble=count++;
	  i++;
	  if (argc==i)
	    {
	      fprintf (stderr,"I Expect at least 1 argument to '-assemble'\n");
	      Usage(argv[0]);
	    }
	  
	  if (isalpha(argv[i][0]))
	    {
	      ass_field = std::string(argv[i]);
	      i++;
	    }
	  
	  if ((argc==i)||((argv[i][0]=='-')&&(!isdigit(argv[i][1]))))
	    {
	      fprintf (stderr,"I Expect at least 1 numeric argument to '-assemble'\n");
	      Usage(argv[0]);
	    }

	  ass_threshold = atof(argv[i]);
	  i++;
	  if ((argc==i)||((argv[i][0]=='-')&&(!isdigit(argv[i][1]))))
	    {
	      ass_angle = ass_threshold;
	      ass_threshold = -std::numeric_limits<double>::max();
	    }
	  else if (argv[i][0]=='-')
	    {	      
	      if (isdigit(argv[i][1]))
		{
		  ass_angle = atof(argv[i]);
		  i++;
		}
	    }
	  else 
	    {
	      ass_angle = atof(argv[i]);
	      i++;
	    }
	}
      else if (!strcmp(argv[i],"-info"))
	{
	  i++;
	  Opt_info=1;
	}
      else if ((!strcmp(argv[i],"-trimBelow"))||(!strcmp(argv[i],"-trimAbove"))||(!strcmp(argv[i],"-trim")))
	{
	  if (!strcmp(argv[i],"-trimAbove")) trim_above[Opt_nTrim]=true;
	  Opt_trim[Opt_nTrim]=count++;
	  i++;
	  if (argc==i)
	    {
	      fprintf (stderr,"I Expect at least one argument to '-trim'\n");
	      Usage(argv[0]);
	    }

	  if (argv[i][0]=='-')
	    {	      
	      if (isdigit(argv[i][1]))
		{
		  trim_threshold[Opt_nTrim] = atof(argv[i]);
		  i++;
		}
	      else
		{
		  fprintf (stderr,"Invalid arguments for '-trim'\n");
		  Usage(argv[0]);
		}
	    }
	  else
	    {
	      if (!isdigit(argv[i][0]))
		{
		  trim_field[Opt_nTrim] = std::string(argv[i]);
		  i++;
		}
		
	      if (argc==i)
		{
		  fprintf (stderr,"Invalid arguments for '-trim'\n");
		  Usage(argv[0]);
		}

	      trim_threshold[Opt_nTrim] = atof(argv[i]);
	      i++;
	    }
	  Opt_nTrim++;
	}
      else if (!strcmp(argv[i],"-breakdown"))
	{
	    Opt_breakdown=count++;
	    i++;
	}
      
      else if (!strcmp(argv[i],"-rmOutside"))
	{
	    Opt_rmOutside=count++;
	    i++;
	}
      else if (!strcmp(argv[i],"-rmBoundary"))
	{
	    Opt_rmBoundary=count++;
	    i++;
	}
      
      else if (!strcmp(argv[i],"-toRaDecZ"))
	{
	  if (Opt_toRaDecDist)
	    {
	      fprintf(stderr,"ERROR: options '-toRaDecZ' and '-toRaDecDist' are NOT compatible.\n");
	      Usage(argv[0]);
	    }
	  Opt_toRaDecZ=count++;
	  i++;
	}
      else if (!strcmp(argv[i],"-toRaDecDist"))
	{
	  if (Opt_toRaDecZ)
	    {
	      fprintf(stderr,"ERROR: options '-toRaDecZ' and '-toRaDecDist' are NOT compatible.\n");
	      Usage(argv[0]);
	    }
	  Opt_toRaDecDist=count++;
	  i++;
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
      else if (!strcmp(argv[i],"-smooth"))
	{
	    i++;
	    if ((argc==i)||(argv[i][0]=='-'))
	    {
		fprintf (stderr,"I Expect a number as argument of '-smooth'\n");
		Usage(argv[0]);
	    }
	    Opt_smoothN = atoi(argv[i]);
	    Opt_smooth=count++;
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
	  if (!ndskel::IO::canSave(Opt_to.back()))
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

  if (!Opt_outName) strcpy(outname,CutName(fName));
  
  if (Opt_outDir) {
    strcpy(tmpname,outdir);
    if (outdir[strlen(outdir)-1]!='/') sprintf(outdir,"%s/",tmpname);
    sprintf(tmpname,"%s%s",outdir,outname);
    strcpy(outname,tmpname);
  }

  /*
  sprintf(outname,"%s.NDskl",outname);
  else sprintf(outname,"%s.ASCIIskl",outname);
  */
  bool needToSave=false;
  NDskel *skl=ndskel::IO::load(fName);
  int ct;
  for (ct=1;ct<count;ct++)
    {
      if (Opt_smooth==ct) 
	{
	  needToSave=true;
	  strcpy(tmpname,outname);
	  if (!Opt_noTags) sprintf(outname,"%s.S%3.3d",tmpname,Opt_smoothN);
	  printf("Smoothing skeleton %d times ... ",Opt_smoothN);fflush(0);
	  for (i=0;i<Opt_smoothN;i++)
	    skl = MatrixSmoothSkeleton(skl,1);
	  printf("done.\n");
	}
      
      if (Opt_breakdown==ct)
	{
	  needToSave=true;
	  strcpy(tmpname,outname);
	  if (!Opt_noTags) sprintf(outname,"%s.BRK",tmpname);
	  printf("Breaking down skeleton ... ");fflush(0);
	  NDskel *skl2=NDskelBreakdown(skl);
	  printf("done.\nBroken down skeleton has %d(+%d) nodes and %d(%d) segments.\n",
		 skl2->nnodes,skl2->nnodes-skl->nnodes,
		 skl2->nsegs,skl2->nsegs-skl->nsegs);
	  Free_NDskeleton(&skl);
	  skl=skl2;
	}

      if (Opt_assemble==ct)
	{
	  needToSave=true;
	  strcpy(tmpname,outname);
	  if (!Opt_noTags) sprintf(outname,"%s.ASMB",tmpname);
	  //printf("Assembling skeleton ... ");fflush(0);
	  NDskel *sklT=NDskelAssembleT(skl,ass_field).build(ass_threshold, ass_angle);//NDskelAssemble(skl,ass_threshold,ass_angle);   
	  NDskel *skl2=NDskelTrimT(sklT).trimBelow(ass_threshold,ass_field);
	  printf("Assembled skeleton: %d(%d) nodes and %d(%d) segments left.\n",
		 skl2->nnodes,skl2->nnodes-skl->nnodes,
		 skl2->nsegs,skl2->nsegs-skl->nsegs);
	  Free_NDskeleton(&sklT);
	  Free_NDskeleton(&skl);
	  skl=skl2;
	}
      for (i=0;i<Opt_nTrim;i++)
	if (Opt_trim[i]==ct)
	  {
	    needToSave=true;
	    strcpy(tmpname,outname);
	    if (!Opt_noTags) sprintf(outname,"%s.TRIM",tmpname);
	    //printf("Assembling skeleton ... ");fflush(0);
	    //NDskel *skl2=NDskelAssembleT(skl,trim_field[i]).trim(trim_threshold[i]);
	    NDskel *skl2;

	    if (trim_above[i])
	      skl2=NDskelTrimT(skl).trimAbove(trim_threshold[i],trim_field[i]);
	    else
	      skl2=NDskelTrimT(skl).trimBelow(trim_threshold[i],trim_field[i]);

	    printf("Trimmed skeleton: %d(%d) nodes and %d(%d) segments left.\n",
		   skl2->nnodes,skl2->nnodes-skl->nnodes,
		   skl2->nsegs,skl2->nsegs-skl->nsegs);
	    Free_NDskeleton(&skl);
	    skl=skl2;
	  }
      
      if (Opt_rmBoundary==ct)
	  {
	    needToSave=true;
	    strcpy(tmpname,outname);
	    if (!Opt_noTags) sprintf(outname,"%s.rmB",tmpname);
	    NDskel *skl2=NDskelTrimT(skl).trimBoundary(true);
	    printf(" Trimmed boundary: %d(%d) nodes and %d(%d) segments left.\n",
		   skl2->nnodes,skl2->nnodes-skl->nnodes,
		   skl2->nsegs,skl2->nsegs-skl->nsegs);
	    Free_NDskeleton(&skl);
	    skl=skl2;
	  }
      
      if (Opt_rmOutside==ct)
	  {
	    needToSave=true;
	    strcpy(tmpname,outname);
	    if (!Opt_noTags) sprintf(outname,"%s.rmO",tmpname);
	    NDskel *skl2=NDskelTrimT(skl).trimBoundary(false);
	    printf(" Trimmed boundary (out): %d(%d) nodes and %d(%d) segments.\n",
		   skl2->nnodes,skl2->nnodes-skl->nnodes,
		   skl2->nsegs,skl2->nsegs-skl->nsegs);
	    Free_NDskeleton(&skl);
	    skl=skl2;
	  }
      
      if (Opt_toRaDecZ==ct)
	{
	  needToSave=true;
	  strcpy(tmpname,outname);
	  if (!Opt_noTags) sprintf(outname,"%s.RaDecZ",tmpname);
	  if (skl->ndims!=3)
	    {
	      fprintf(stderr,"ERROR converting to (Ra,Dec,z) : ndims = %d\n",skl->ndims);
	      exit(0);
	    }
	  for (i=0;i<skl->nnodes;i++)
	    {
	      sampledDataInput::cartesian2spherical<float>(skl->Node[i].pos,skl->Node[i].pos,skl->ndims);
	      skl->Node[i].pos[2]=survey_dist2red(skl->Node[i].pos[2]);
	    }
	  for (i=0;i<skl->nsegs;i++)
	    {
	      sampledDataInput::cartesian2spherical<float>(skl->Seg[i].pos,skl->Seg[i].pos,skl->ndims);
	      sampledDataInput::cartesian2spherical<float>(skl->Seg[i].pos+skl->ndims,skl->Seg[i].pos+skl->ndims,skl->ndims);
	      skl->Seg[i].pos[2]=survey_dist2red(skl->Seg[i].pos[2]);
	      skl->Seg[i].pos[2+3]=survey_dist2red(skl->Seg[i].pos[2+3]);
	    }
	}
      
      if (Opt_toRaDecDist==ct)
	{
	  needToSave=true;
	  strcpy(tmpname,outname);
	  if (!Opt_noTags) sprintf(outname,"%s.RaDecDist",tmpname);
	  for (i=0;i<skl->nnodes;i++)
	    sampledDataInput::cartesian2spherical<float>(skl->Node[i].pos,skl->Node[i].pos,skl->ndims);
	  for (i=0;i<skl->nsegs;i++)
	    {
	      sampledDataInput::cartesian2spherical<float>(skl->Seg[i].pos,skl->Seg[i].pos,skl->ndims);
	      sampledDataInput::cartesian2spherical<float>(skl->Seg[i].pos+skl->ndims,skl->Seg[i].pos+skl->ndims,skl->ndims);
	    }
	}
    }

  if (Opt_addField)
    {
      needToSave=true;
      for (i=0;i<Opt_addField;i++)
	{
	  NDfield *f=ndfield::IO::load(addFieldFile[i].c_str());
	  addDataFieldFromNDfield(skl,f,addFieldName[i].c_str(),0);
	  Free_NDfield(&f);
	}
    }
  
  if (Opt_info) {
    printf("Skeleton statistics:\n");
    printNDskelStat(skl,3);
  }

#ifdef HAVE_CFITS_IO
  if (Opt_toFITS)
    {
      if ((FITSRes.size()!=0)&&(FITSRes.size()!=skl->ndims))
	{
	  fprintf(stderr,"ERROR in 'toFITS': FITS image is %dD, but %ld dimensions were specified.\n",skl->ndims,FITSRes.size());
	}
      else 
	{
	  if (FITSRes.size()==0)
	    Save_FITSskel(skl,(std::string(outname)+std::string(".fits")).c_str(),NULL);
	  else
	    Save_FITSskel(skl,(std::string(outname)+std::string(".fits")).c_str(),&FITSRes[0]);
	}
      
    }
#endif

  if (Opt_to.size())
    {
      for (i=0;i<Opt_to.size();i++)
	ndskel::IO::save(skl,std::string(outname)+ndskel::IO::getExtension(Opt_to[i]),Opt_to[i]);
    }
  else if (needToSave) 
    ndskel::IO::save(skl,std::string(outname)+ndskel::IO::getExtension());

  Free_NDskeleton(&skl);
}
