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
#ifndef __NDFIELD_INTERFACE_HXX__
#define __NDFIELD_INTERFACE_HXX__

#include <string>
#include <limits>

#include "NDfield.h"
#include "NDfield_VTK_IO.h"
#include "NDnetwork.h"
#include "genericIO.hxx"
#include "global.h"
#include "mystring.h"
#include "distances.h"

//#include "sampledDataInput.hxx"
#include "asciiSurvey.h"

#ifdef HAVE_SDL
#include "SDL/SDL.h"
#include "SDL/SDL_image.h"
#endif 

namespace ndfield {

  typedef genericIO_interface<NDfield> interface;
  
  class fromNDfield : public interface {
    std::string getTypeStr() {return std::string("NDfield");}
    std::string getExtensionStr() {return std::string("ND");}
    bool canLoad(std::string fname) 
    {
      if (fname==std::string("")) return true;
      if (IsNDfield(fname.c_str())) return true;
      return false;
    }
    bool canSave() {return true;}
    
    NDfield *load(std::string fname)
    {
      return Load_NDfield(fname.c_str());
    }
    
    int save(NDfield *f, std::string fname)
    {
      return Save_NDfield(f,fname.c_str());
    }
  };
  
  class fromASCII : public interface {
    std::string getTypeStr() {return std::string("NDfield_ascii");}
    std::string getExtensionStr() {return std::string("a.ND");}
    bool canLoad(std::string fname) 
    {
      
      if (fname==std::string("")) return true;
      FILE *f = fopen(fname.c_str(),"r");
      if (f==NULL) return false;
      char *line=NULL;
      int n;
      Mygetline(&line,&n,f);      
      fclose(f);
      
      if (strstr(line,NDFIELD_ASCII_TAG)!=NULL) 
	{
	  free(line);
	  return true;
	}
      free(line);
      return false;
    }
    bool canSave() {return true;}
    
    NDfield *load(std::string fname)
    {
      long i;      
      FILE *f = fopen(fname.c_str(),"r");
      if (f==NULL) return NULL;
      char *line=NULL;
      int n;
      long dummy;
      char comment[80];

      strcpy(comment,"");
      
      Mygetline(&line,&n,f);  
      if (strstr(line,NDFIELD_ASCII_TAG)==NULL) 
	{
	  free(line);
	  return NULL;
	}   
      
      bool coord=false;
      char *str[1000];
      
      if (strstr(line,"COORDS")!=NULL)
	coord=true;
      
      Mygetline(&line,&n,f);
    
      if ((strstr(line,"[")==NULL)||(strstr(line,"]")==NULL))
	{
	  fprintf(stderr,"ERROR trying to load file '%s'.\n",fname.c_str());
	  fprintf(stderr,"   Not a valid ASCII NDfield format.\n");
	  fprintf(stderr,"   --> could not read dimensions.\n");
	  free(line);
	  return NULL;
	}

      if (strstr(line,"BBOX"))
	{
	  fprintf(stderr,"ERROR trying to load file '%s'.\n",fname.c_str());
	  fprintf(stderr,"   bounding box must be defined AFTER dimensions.\n");	  
	  free(line);
	  exit(-1);
	}
      if (line[0]=='#')
	{
	  fprintf(stderr,"ERROR trying to load file '%s'.\n",fname.c_str());
	  fprintf(stderr,"   comments must be defined AFTER dimensions.\n");	  
	  free(line);
	  exit(-1);
	}

      char *c=strstr(line,"[");
      std::vector<int> dims;
      while (*c!=']')
	{
	  if (!isdigit(*c)) c++;
	  else
	    {
	      int v;
	      sscanf(c,"%d",&v);
	      dims.push_back(v);
	      while(isdigit(*c)) c++;
	    }
	};

      int ndims=coord?dims[0]:dims.size();
      long nval=1;
      for (i=0;i<(long)dims.size();i++) nval*=dims[i];
      std::vector<double> x0;
      std::vector<double> delta;
      dummy=Mygetline(&line,&n,f);
      if (line[0]=='#')
	{
	  if (dummy>80) line[79]='\0';
	  strcpy(comment,line);
	  dummy=Mygetline (&line,&n,f);
	}

      if ((c=strstr(line,"BBOX"))!=NULL)
	{
	  int ntok=str2tok(c,"[ ],;|",6,0,0,str);
	  if (ntok != ndims*2+1)
	    {
	      fprintf(stderr,"ERROR trying to load file '%s'.\n",fname.c_str());
	      fprintf(stderr,"  Bad definition of the bounding box, should have ndims=%d\n",ndims);
	      exit(-1);
	    }
	  for (i=1;i<ntok;i++)
	    {
	      if ((long)x0.size()<ndims) x0.push_back(strtod(str[i],NULL));
	      else delta.push_back(strtod(str[i],NULL));
	    }
	  Mygetline(&line,&n,f);
	}
      
      char dimstr[256];
      sprintf(dimstr,"[%d",dims[0]);
      for (i=1;i<(long)dims.size();i++) sprintf(dimstr,"%s,%d",dimstr,dims[i]);
      strcat(dimstr,"]");
      if (verbose>=1) printf("Reading %s %s from NDfied ASCII file '%s' ... ",dimstr,coord?"coords":"grid",fname.c_str());
      fflush(0);
      
      double *tab=(double*)malloc(sizeof(double)*nval);
      long nread=0;
      do {
	
	//Mygetline(&line,&n,f);
	//printf ("line : %s\n",line);
	if (feof(f)) break;
	//printf ("line : %s\n",line);
	int ntok=str2tok(line," ,;|[]",6,'#',0,str);
	if ((nread+ntok) > nval) 
	  {
	    ntok=nval-nread;
	    fprintf(stderr,"WARNING: Too many values to read,stopping.\n");
	  }
	for (i=0;i<ntok;i++) tab[nread++]=strtod(str[i],NULL);
	Mygetline(&line,&n,f);
      } while (nread<nval);
      free(line);

      if (nread!=nval)
	{
	  fprintf(stderr,"ERROR: too few values to read, aborting. (%ld/%ld)\n",nread,nval);
	  free(tab);
	  return NULL;
	}

      if (!x0.size())
	{
	  if (!coord)
	    {
	      x0.assign(ndims,0);
	      delta.assign(dims.begin(),dims.end());
	    }
	  else
	    {
	      long j;
	      x0.assign(ndims,std::numeric_limits<double>::max());
	      delta.assign(ndims,-std::numeric_limits<double>::max());
	      for (i=0;i<nval;i+=dims[0])
		for (j=0;j<dims[0];j++)
		  {
		    if (tab[i+j]<x0[j]) x0[j]=tab[i+j];
		    if (tab[i+j]>delta[j]) delta[j]=tab[i+j];
		  }
	      
	      for (j=0;j<dims[0];j++) delta[j]-=x0[j];
	      
	    }
	}
      NDfield *field=Create_NDfield(&dims[0],ndims,coord?1:0,ND_DOUBLE,&x0[0],&delta[0],(void*)(tab),comment);

      if (verbose>=1) printf("done. \n");
      return field;
    }
    
    int save(NDfield *field, std::string fname)
    {
      FILE *f;
      long i,j;
      f=fopen(fname.c_str(),"w");
      if (f==NULL) return -1;
      if (field->fdims_index) fprintf(f,"%s COORDS\n",NDFIELD_ASCII_TAG);
      else fprintf(f,"%s\n",NDFIELD_ASCII_TAG);
      fprintf(f,"[%ld",(long)field->dims[0]);
      for (i=1;i<field->n_dims;i++) fprintf(f," %ld",(long)field->dims[i]);
      fprintf(f,"]\n");
      fprintf(f,"#%s\n",field->comment);
      fprintf(f,"BBOX [%g",field->x0[0]);
      for (i=1;i<field->ndims;i++) fprintf(f,",%g",field->x0[i]);
      fprintf(f,"] [%g",field->delta[0]);
      for (i=1;i<field->ndims;i++) fprintf(f,",%g",field->delta[i]);
      fprintf(f,"]\n");
      
      DEFINE_NDVARS(val);
      SETNDPOINTER(val,field->val,field->datatype);

      double tmp=-1;
      if (field->fdims_index)
	{
	  //printf("%d %d\n",field->dims[0],field->dims[1]);
	  for (i=0;i<field->dims[1];i++)
	    {
	      SETNDFIELD_VAL(tmp,val,i++,field->datatype);
	      fprintf(f,"%g",tmp);
	      for (j=1;j<field->dims[0];j++) 
		{
		  SETNDFIELD_VAL(tmp,val,i++,field->datatype);
		  fprintf(f," %g",tmp);
		}
	      fprintf(f,"\n");
	    }
	}
      else
	{
	  for (i=0;i<field->nval;i++)
	    {
	      SETNDFIELD_VAL(tmp,val,i,field->datatype);
	      fprintf(f,"%g\n",tmp);
	    }
	}
      
      fclose(f);
      return 0;
      //return Save_NDfield(f,fname.c_str());
    }
  };

  class fromfits : public interface {
    std::string getTypeStr() {return std::string("fits");}
    std::string getExtensionStr() {return std::string("fits");}
    bool canLoad(std::string fname) 
    {
      if (fname==std::string("")) return true;
      if (IsNDfieldFromFITS(fname.c_str())) return true; 
      return false;
    }
    bool canSave() {return false;}

    NDfield *load(std::string fname)
    {
      return Load_NDfieldFromFITS(fname.c_str());
    }

    int save(NDfield *f, std::string fname)
    {
      return -1;
    }
  };

#ifdef HAVE_SDL
  class fromSDL : public interface {
  private:
    std::string SDL_imageType(SDL_RWops *rwop)
    {
      std::string result;
   
      if (IMG_isBMP(rwop)) result=std::string("BMP");
      else if (IMG_isCUR(rwop)) result=std::string("CUR");
      else if (IMG_isGIF(rwop)) result=std::string("GIF");
      else if (IMG_isICO(rwop)) result=std::string("ICO");
      else if (IMG_isJPG(rwop)) result=std::string("JPG");
      else if (IMG_isLBM(rwop)) result=std::string("LBM");
      else if (IMG_isPCX(rwop)) result=std::string("PCX");
      else if (IMG_isPNG(rwop)) result=std::string("PNG");
      else if (IMG_isPNM(rwop)) result=std::string("PNM");
      //else if (IMG_isTGA(rwop)) result=std::string("TGA"); //does not work
      else if (IMG_isTIF(rwop)) result=std::string("TIF");
      else if (IMG_isXCF(rwop)) result=std::string("XCF");
      else if (IMG_isXPM(rwop)) result=std::string("XPM");
      else if (IMG_isXV(rwop)) result=std::string("XV");
      //exit(0);
      return result;
    }

    Uint32 getpixel(SDL_Surface *surface, int x, int y)
    {
      int bpp = surface->format->BytesPerPixel;
      Uint8 *p = (Uint8 *)surface->pixels + y * surface->pitch + x * bpp;
    
      switch (bpp) {
      case 1:
	return *p;
      
      case 2:
	return *(Uint16 *)p;
      
      case 3:
	if (SDL_BYTEORDER == SDL_BIG_ENDIAN)
	  return p[0] << 16 | p[1] << 8 | p[2];
	else
	  return p[0] | p[1] << 8 | p[2] << 16;
      
      case 4:
	return *(Uint32 *)p;
      
      default:
	return 0;     
      } 
    }

  public:

    std::string getTypeStr() {return std::string("SDL-image");}
    std::string getExtensionStr() {return std::string("SDL");}
    bool canLoad(std::string fname) 
    {
      if (fname==std::string("")) return true;

      SDL_RWops *rwop=SDL_RWFromFile(fname.c_str(), "rb");
      if (SDL_imageType(rwop) != std::string(""))
	{
	  SDL_FreeRW(rwop);
	  return true;
	}
      SDL_FreeRW(rwop); 
      return false;
    }
    bool canSave() {return false;}

    NDfield *load(std::string fname)
    {       
      SDL_RWops *rwop=SDL_RWFromFile(fname.c_str(), "rb"); 
      std::string IMGtype=SDL_imageType(rwop);

      if (IMGtype != std::string(""))
	{
	  NDfield *f=NULL;
	  long i,j,k;
	  char ftype[10];
	  strcpy(ftype,IMGtype.c_str());
	  SDL_Surface *image=IMG_LoadTyped_RW(rwop, 1, ftype);

	  if ( !image )
	    {
	      printf ( "IMG_Load: %s\n", IMG_GetError () );
	      
	    }

	  printf ("Loading [%d,%d] array from file %s ... ",image->w,image->h,fname.c_str());fflush(0);
	  	
	  double *val = (double *)malloc(image->w*image->h*sizeof(double));
	  //double norm = 1./(255.*(0.2989+0.5870+0.1140));
	  for (j=image->h-1,k=0;j>=0;j--)
	    for (i=0;i<image->w;i++,k++)	  
	      {
		Uint8 r,g,b;
		Uint32 p=getpixel(image,i,j);
		SDL_GetRGB(p,image->format,&r,&g,&b);
		val[k] = (double)((0.2989*(double)r + 0.5870*(double)g + 0.1140*(double)b));
	      }
	  int dims[2];
	  double x0[2];
	  double delta[2];
	  x0[0]=x0[1]=0;
	  delta[0]=image->w;
	  delta[1]=image->h;
	  dims[0]=image->w;
	  dims[1]=image->h;
	  f=Create_NDfield(dims,2,0,ND_DOUBLE,x0,delta,(void*)val,"SDL_image");
	 
	  printf("done.\n");
	  return f;
	}
      
      SDL_FreeRW(rwop);
      return NULL;
    }

    int save(NDfield *f, std::string fname)
    {
      return -1;
    }
  };
#endif
  
  class fromASCIISurvey : public interface {
    std::string getTypeStr() {return std::string("survey_ascii");}
    std::string getExtensionStr() {return std::string(".survey");}
    bool canLoad(std::string fname) 
    {
      if (fname==std::string("")) return true;
      if (isAsciiSurvey(fname.c_str())) return true;
      return false;
    }
    bool canSave() {return false;}

    NDfield *load(std::string fname)
    {
      asciiSurvey *s=readAsciiSurvey(fname.c_str());
      int dims[2];
      double *p=NULL;
      double x0[s->ndims];
      double delta[s->ndims];
      
      asciiSurveyGetCoords_d(s,&p);
      dims[0]=s->ndims;
      dims[1]=s->npart;

      long i,j;
      for (i=0;i<s->ndims;i++) x0[i]=delta[i]=p[i];
      for (i=0;i<s->npart;i++)
	{
	  for (j=0;j<s->ndims;j++)
	    {
	      if (x0[j]>p[i*s->ndims+j]) x0[j]=p[i*s->ndims+j];
	      if (delta[j]<p[i*s->ndims+j]) delta[j]=p[i*s->ndims+j];
	    }
	}
      double vol=1.;
      
      for (j=0;j<s->ndims;j++) {
	delta[j] -= x0[j];
	vol*=delta[j];		
      }
	    
      double dec = 0.5*pow(vol/s->npart,1./s->ndims);
      //double dec = 0.1 * pow(vol,1./ndims);
      
      for (j=0;j<s->ndims;j++) {
	delta[j]+=dec;
	x0[j]-=dec/2;
      }
      
      char tmpc[255];
      char cp[255];
      cosmoD_getParamsStr(cp);
      sprintf(tmpc,"%s",cp);
      //sprintf(tmpc,"from ASCII survey (om=%2.2f ol=%2.2f H0=%2.2f)",OMEGAM_DEFAULT,OMEGAL_DEFAULT,HUBBLE_DEFAULT);

      NDfield *f=Create_NDfield(dims,s->ndims,1,ND_DOUBLE,x0,delta,(void*)p,tmpc);         
      freeSurvey(&s);

      return f;
    }

    int save(NDfield *f, std::string fname)
    {
      return -1;
    }
  };

  class dummy : public interface {
    std::string getTypeStr() {return std::string("dummy");}
    std::string getExtensionStr() {return std::string("dum");}
    bool canLoad(std::string fname) 
    {
      if (fname==std::string("")) return true;
      return false;
    }
    bool canSave() {return false;}

    NDfield *load(std::string fname)
    {
      return NULL;
    }

    int save(NDfield *f, std::string fname)
    {
      return -1;
    }
  };
  
  class fromVTK : public interface {
    std::string getTypeStr() {return std::string("vtk");}
    std::string getExtensionStr() {return std::string("vtk");}
    bool canLoad(std::string fname) {return false;}
    bool canSave() {return true;}
    NDfield *load(std::string fname)
    {
      return Load_NDfieldFromVTK(fname.c_str());
    }
    
    int save(NDfield *f, std::string fname)
    {
      if (f->fdims_index)
	{
	  fprintf(stderr,"ERROR writing VTK file: NDfield has type 'coords', should be a 'grid'.\n");
	  fprintf(stderr," Try converting to 'NDnet' format, and then use 'netconv'.\n");
	  return -1;
	}
      return Save_NDfieldToVTK(f,fname.c_str(),NDFIELD_VTK_BIN|NDFIELD_VTK_LEGACY);
    }
  };

  class fromVTK_ascii : public interface {
    std::string getTypeStr() {return std::string("vtk_ascii");}
    std::string getExtensionStr() {return std::string("a.vtk");}
    bool canLoad(std::string fname) {return false;}
    bool canSave() {return true;}
    NDfield *load(std::string fname)
    {
      return NULL;
    }
    
    int save(NDfield *f, std::string fname)
    {
      if (f->fdims_index)
	{
	  fprintf(stderr,"ERROR writing VTK file: NDfield has type 'coords', should be a 'grid'.\n");
	  fprintf(stderr," Try converting to 'NDnet' format, and then use 'netconv'.\n");
	  return -1;
	}
      return Save_NDfieldToVTK(f,fname.c_str(),NDFIELD_VTK_ASCII|NDFIELD_VTK_LEGACY);
    }
  };

  class fromVTK_xml : public interface {
    std::string getTypeStr() {return std::string("vti");}
    std::string getExtensionStr() {return std::string("vti");}
    bool canLoad(std::string fname) {return false;}
    bool canSave() {return true;}
    NDfield *load(std::string fname)
    {
      return NULL;
    }
    
    int save(NDfield *f, std::string fname)
    {
      if (f->fdims_index)
	{
	  fprintf(stderr,"ERROR writing VTK file: NDfield has type 'coords', should be a 'grid'.\n");
	  fprintf(stderr," Try converting to 'NDnet' format, and then use 'netconv'.\n");
	  return -1;
	}
      return Save_NDfieldToVTK(f,fname.c_str(),NDFIELD_VTK_BIN|NDFIELD_VTK_XML);
    }
  };

  class fromVTK_xml_ascii : public interface {
    std::string getTypeStr() {return std::string("vti_ascii");}
    std::string getExtensionStr() {return std::string("a.vti");}
    bool canLoad(std::string fname) {return false;}
    bool canSave() {return true;}
    NDfield *load(std::string fname)
    {
      return NULL;
    }
    
    int save(NDfield *f, std::string fname)
    {
      if (f->fdims_index)
	{
	  fprintf(stderr,"ERROR writing VTK file: NDfield has type 'coords', should be a 'grid'.\n");
	  fprintf(stderr," Try converting to 'NDnet' format, and then use 'netconv'.\n");
	  return -1;
	}
      return Save_NDfieldToVTK(f,fname.c_str(),NDFIELD_VTK_ASCII|NDFIELD_VTK_XML);
    }
  };

  class fromNDnet : public interface {
    std::string getTypeStr() {return std::string("NDnet");}
    std::string getExtensionStr() {return std::string("NDnet");}
    bool canLoad(std::string fname) {return false;}
    bool canSave() {return true;}
    NDfield *load(std::string fname) {return NULL;}
    
    int save(NDfield* field, std::string fname)
    {
      NDnetwork *net=NULL;    
      long i;

      if (field->fdims_index)
	{
	  net=CreateNetwork(field->dims[0],field->dims[1],0);
	  memcpy(net->x0,field->x0,sizeof(double)*net->ndims);
	  memcpy(net->delta,field->delta,sizeof(double)*net->ndims);      
	  DEFINE_NDVARS(d);
	  SETNDPOINTER(d,field->val,field->datatype);
	  for (i=0;i<(long)net->nvertex*net->ndims;i++)
	    SETNDFIELD_VAL(net->v_coord[i],d,i,field->datatype);
	}
      else
	{
	  fprintf(stderr,"ERROR writing NDnet file: NDfield has type 'grid', should be 'coords'.\n");
	  fprintf(stderr," Convertion from 'grid' type to 'NDnet' is not implemented yet.\n");
	  return -1;
	
	}
     
      Save_NDnetwork(net,fname.c_str());
      FreeNDnetwork(&net);
      return 0;
    }
  };


  struct allTypes
  {
    template <class OutputIterator>
    static void generate(OutputIterator out)
    {
      *out = new fromNDfield();
      *out = new fromASCII();
      *out = new fromNDnet();
      *out = new fromfits();      
      *out = new fromASCIISurvey();
#ifdef HAVE_SDL
      *out = new fromSDL();
#endif
      *out = new fromVTK();
      *out = new fromVTK_ascii();
      *out = new fromVTK_xml();
      *out = new fromVTK_xml_ascii();
    }
    static std::string type()
    {
      return std::string("field");
    }
  };

  typedef genericIO<interface, allTypes> IO;
}

#endif
