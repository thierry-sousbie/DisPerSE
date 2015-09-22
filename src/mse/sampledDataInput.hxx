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
#ifndef __SAMPLED_INPUT_HXX__
#define __SAMPLED_INPUT_HXX__

#include <string>
#include <algorithm>
#include <vector>

#include "gadget_io.h"
#include "asciiSurvey.h"
#include "mystring.h"
#include "NDfield.h"
#include "mytypes.h"

#include "NDfield_interface.hxx"

class sampledDataInput {
public:
    enum fTypes {ft_notset, ft_unknown, ft_ascii, ft_gadget, ft_ndfield};
    enum cSystems {cs_notset, cs_unknown, cs_cartesian, cs_spherical };

    typedef std::pair<fTypes,cSystems> fileProperties;

private:

  typedef std::string string;
  std::string fname;
  //bool  valid;

  fTypes fileType;
  cSystems coordSystem;

  bool loaded;

  asciiSurvey *S;
  snapshot_data *G;
  NDfield *F;
  std::string comment;

  std::vector<float> posv;
  std::vector<float> posv_s;


  float *pos; // always cartesian
  float *pos_s; // always spherical

  float *vel;
  float *mass;
  long npart;
  int ndims;

  std::vector<double> x0;
  std::vector<double> delta;

  std::vector<double> x0s;
  std::vector<double> deltas; 

public:
  // template<class OutputIterator> 
  //OutputIterator getSubBox((std::vector<double> &x0, std::vector<double> &delta, double margin, OutputIterator it)

  std::string getComment() {return comment;}

  template <class T>
  static void spherical2cartesian(T *src, T *dest, int ndims)
  {
    if ((ndims!=2)&&(ndims!=3))
      {
	fprintf(stderr,"ERROR in spherical2cartesian: spherical coordinates only in 2D and 3D.\n");
	exit(-1);
      }

    double ra = src[0];
    double dec = src[1];	   
    double dist;
    
    if (ndims==3) dist=src[2];
    else dist=1.;

    double dcr=dist*cos(dec*DEG2RAD);

    dest[0] = dcr*cos(ra*DEG2RAD);
    dest[1] = dcr*sin(ra*DEG2RAD);
    if (ndims == 3) dest[2] = dist*sin(dec*DEG2RAD);
    /*
    T d[3];
    T s[3];
    memcpy(s,dest,sizeof(T)*ndims);
    
    cartesian2spherical(s,d,ndims);
    for (int i=0;i<3;i++) printf (" %g==%g\n",src[i],d[i]);
    */
    
  }
  
  template <class T>
  static void cartesian2spherical(T *src, T *dest, int ndims,bool radians=false)
  {
    //double pi = 3.141592653589793116;
	    
    if (ndims!=3)
      {
	fprintf(stderr,"ERROR in spherical2cartesian: spherical coordinates only in 2D and 3D.\n");
	exit(-1);
      }
    
    double x = src[0];
    double y = src[1];	   
    double z = src[2];
    double d;
    d=sqrt(x*x+y*y+z*z);

    dest[2]=d;

    if (d==0) dest[1]=0;
    else dest[1]=asin(z/d);
       
    if ((y==0)||(x==0)) dest[0]=0;
    else dest[0]=atan2(y,x);

    //dest[0]+=PI;
    if (dest[0]<0) dest[0]+=2.*PI;  
    
    if (!radians) {
      dest[0]*=RAD2DEG;
      dest[1]*=RAD2DEG;
    }
    /*
    T dd[3];
    T ss[3];
    
    memcpy(ss,dest,sizeof(T)*ndims);
    
    spherical2cartesian(ss,dd,ndims);
    if (x!=dd[0]) printf (" %g==%g\n",x,dd[0]);
    if (y!=dd[1]) printf (" %g==%g\n",y,dd[1]);
    if (z!=dd[2]) printf (" %g==%g\n",z,dd[2]);
    */
    //for (int i=0;i<3;i++) if (src[i]!=dd[i]) printf (" %g==%g\n",src[i],dd[i]);
    /*
    if (!radians) {
      dest[0]*=RAD2DEG;
      dest[1]*=RAD2DEG;
    }
    */

    /*
    dest[0]=acos(z/d)*180./pi -90.;
    dest[1]=atan2(y,x)*180./pi;
    if (dest[1]<0) dest[1]+=360.;
    */
    
    
  }
  
    

private:

  void setCartesianFromSpherical() 
  {
    if ((ndims!=2)&&(ndims!=3))
      {
	fprintf(stderr,"ERROR in setCartesianFromSpherical(): spherical coordinates only in 2D and 3D.\n");
	exit(-1);
      }
    
    int i=0;
    posv.resize(ndims*npart);
    pos = &(posv[0]);
    
    for (i=0;i<npart;i++)
      spherical2cartesian(&pos_s[i*ndims],&pos[i*ndims],ndims);
    
  }

  void setSphericalFromCartesian() 
  {
    if ((ndims!=2)&&(ndims!=3))
      {
	fprintf(stderr,"ERROR in setCartesianFromSpherical(): spherical coordinates only in 2D and 3D.\n");
	exit(-1);
      }
    
    int i=0;
    posv_s.resize(ndims*npart);
    pos_s = &(posv_s[0]);
    
    for (i=0;i<npart;i++)
      spherical2cartesian(&pos[i*ndims],&pos_s[i*ndims],ndims);
    
  }

 
  
    int guessBBox(std::vector<double> &x0, std::vector<double> &delta, bool spherical=false)
	{
	    if (!loaded) return -1;
	    float *p;
	    
	    if (spherical)
	    {	 
	      if (pos_s==NULL) setSphericalFromCartesian(); 
	      
	      p=pos_s;
	      x0s.resize(ndims);
	      deltas.resize(ndims);
	      int i,j;
	      for (i=0;i<ndims;i++) x0s[i]=deltas[i]=p[i];
	      for (i=0;i<npart;i++)
		{
		  
		  for (j=0;j<ndims;j++)
		    {
		      
		      if (x0s[j]>p[i*ndims+j]) x0s[j]=p[i*ndims+j];
		      if (deltas[j]<p[i*ndims+j]) deltas[j]=p[i*ndims+j];
		    }
		}
	      return 1;
	    }
	    
	    if (pos==NULL) setCartesianFromSpherical(); 
	    
	    p=pos;
	    x0.resize(ndims);
	    delta.resize(ndims);
	    int i,j;
	    for (i=0;i<ndims;i++) x0[i]=delta[i]=p[i];
	    for (i=0;i<npart;i++)
	    {
	      for (j=0;j<ndims;j++)
		{
		  if (x0[j]>p[i*ndims+j]) x0[j]=p[i*ndims+j];
		  if (delta[j]<p[i*ndims+j]) delta[j]=p[i*ndims+j];
		}
	    }

	    
	    double vol=1.;
		
	    for (j=0;j<ndims;j++) {
	      delta[j] -= x0[j];
	      vol*=delta[j];		
	    }
	    
	    double dec = 0.5*pow(vol/npart,1./ndims);
	    //double dec = 0.1 * pow(vol,1./ndims);
	      
	      for (j=0;j<ndims;j++) {
		delta[j]+=dec;
		x0[j]-=dec/2;
	      }
	       
	    return 1;
	}

public:  

  static fTypes checkFileType(const std::string filename) {
	fTypes type;
	
	if (filename == string()) type=ft_notset;
	else if (IsGadgetFile(filename.c_str())) type=ft_gadget;
	else if (isAsciiSurvey(filename.c_str())) type=ft_ascii;
	else if (ndfield::IO::canLoad(filename)) type=ft_ndfield;
	else type=ft_unknown;

	return type;
    }

  static fileProperties checkFileProperties(const std::string filename) {
    cSystems system;
    fTypes type = checkFileType(filename);
    int oldverbose=verbose;
    verbose=0;

    if (filename == std::string()) system=cs_notset;
    else if (type==ft_ndfield)
      {
	system=cs_cartesian;
	NDfield *f;
	bool head=false;
	if ((head=IsNDfield(filename.c_str())))
	  f=Load_NDfieldHeader(filename.c_str());
	else
	  f=ndfield::IO::load(filename);
	     
	std::string comment(f->comment);
	std::transform(comment.begin(), comment.end(), comment.begin(), ::toupper);
	if (comment.find("SPHERICAL") != string::npos)
	  system=cs_spherical;
	    
	if (head) free(f);
	else Free_NDfield(&f);

      }
    else if (type==ft_gadget) system=cs_cartesian;
    else if (type==ft_ascii)
      {
	system=cs_cartesian;
	asciiSurvey *s=readAsciiSurvey_header(filename.c_str());
	if ((s->ra!=-1)&&(s->dec!=-1)) system=cs_spherical;	    
	freeSurvey(&s);
      }
    else system=cs_unknown;
    //printf("hello\n");
    verbose=oldverbose;

    return fileProperties(type,system);
  }


  bool isLoaded() const {return loaded;}
  bool hasMass() const {return mass!=NULL;}
  
  bool isValid() const 
  {
    if ((fileType==ft_notset)||(fileType==ft_unknown)) 
      return false;
    int oldverbose=verbose;
    verbose=0;
    if (fileType==ft_ndfield)
      {
	NDfield *tmp;
	bool head=false;
	if ((head=IsNDfield(fname.c_str())))
	  tmp=Load_NDfieldHeader(fname.c_str());
	else
	  tmp=ndfield::IO::load(fname);

	if (tmp->fdims_index==0)
	  {
	    fprintf(stderr,"WARNING (will be ignored) :\n   File %s seems to contain array values, not coordinates.\n   Maybe fdims_index is not correctly set (current value : %d, should be 1).\n",fname.c_str(),tmp->fdims_index);
	  }
	//if (tmp->ndims!=tmp->dims[0]) {free(tmp);return false;}

	if (head) free(tmp);
	else Free_NDfield(&tmp);
	
      }
    verbose=oldverbose;
    return true;
  }

  fTypes getType() const {return fileType;}
    
  cSystems getCoordSystem() const { return coordSystem;}
    
  void clear()
  {
    posv.clear();
    posv_s.clear();
    if (loaded==true) 
      {
	loaded = false;
	
	if (fileType==ft_gadget) {
	  freegadgetstruct(G);G=NULL;
	}
	else if (fileType==ft_ndfield) {
	  Free_NDfield(&F);
	}
	else if (fileType==ft_ascii) {	  
	  freeSurvey(&S);
	}
      }
  }
    
  fTypes setFname(const string filename, bool reset=true)
  {
    if ((loaded)&&(reset)) clear();
    
    fname = filename;
    fileProperties tc=checkFileProperties(fname);
    fileType = tc.first;
    coordSystem = tc.second;
    return fileType;
  }


  bool read() {
    
    if (loaded==true) return true;
    if (!isValid()) return false;
    loaded=true;
    if (fileType==ft_gadget) {
      G=(snapshot_data *)calloc(1,sizeof(snapshot_data));
      ReadGadget(fname.c_str(),G,FLAG_ALL);
      comment=std::string("from gadget file");
      pos=G->Pos;
      pos_s=NULL;
      npart=G->header.npartTotal[1];
      vel=G->Vel;
      mass=G->Mass;
      ndims=3;
      if (G->header.BoxSize>0)
	{
	  x0.assign(3,0);
	  delta.assign(3,G->header.BoxSize);
	}
      else
	{
	  
	}
     
      return true;
    }
    else if (fileType==ft_ndfield) {
      
      F=ndfield::IO::load(fname);
      comment=std::string(F->comment);
      Convert_NDfield(F,ND_FLOAT);
      vel=NULL;
      npart=F->dims[1];
      ndims=F->dims[0];
      if (F->delta[0] != 0) {
	if (coordSystem == cs_cartesian)
	  {
	    pos=(float *)F->val;
	    pos_s=NULL;
	    x0.assign(F->x0,&(F->x0[ndims]));
	    delta.assign(F->delta,&(F->delta[ndims]));
	  }
	else
	  {
	    pos_s=(float *)F->val;
	    //setCartesianFromSpherical();
	    x0s.assign(F->x0,&(F->x0[ndims]));
	    deltas.assign(F->delta,&(F->delta[ndims]));
	    //guessBBox(x0,delta);
	  }
      }
      else
	{
	  if (coordSystem == cs_cartesian)
	  {
	    pos=(float *)F->val;
	    pos_s=NULL;
	  }
	else
	  {
	    pos_s=(float *)F->val;
	  }
	}
            
      return true;
    }
    else if (fileType==ft_ascii) {	  
      S=readAsciiSurvey(fname.c_str());
      char tag[256];
      if (!cosmoD_initialized())
	sprintf(tag,"No cosmological conversion.");
      else
	{cosmoD_getParamsStr(tag);tag[255]='\0';}
      comment=std::string(tag);
      npart=S->npart;
      ndims=S->ndims;
      mass=S->m;
      if (coordSystem == cs_cartesian) {
	pos_s=NULL;
	pos=S->pos;

      }
      else{
	pos_s=S->pos;

      }
      vel=S->vel;
      
      return true;
    }

    loaded=false;
    return false;
  }

  sampledDataInput(const string filename=string(), bool load=true)
  {

    G=NULL;
    mass=NULL;
    loaded=false;
    setFname(filename);
    pos_s=pos=NULL;
    
    if (load) 
      {
	if (!read())
	  {
	    fprintf(stderr,"ERROR in sampledDataInput: could not read file.");
	    exit(-1);
	  }
      }
  }

  ~sampledDataInput() {
      clear();      
  }
  
  long getNPart() {if (loaded) return npart; else return 0;}
  int getNDims() {if (loaded) return ndims; else return 0;}

  float *getPosFPtr(bool spherical=false) 
  {
    if (!loaded) 
      {
	fprintf(stderr,"WARNING in getPosFPtr: data is not loaded");
	return NULL;
      } 
    /*
    if ((spherical)&&(coordSystem!=cs_spherical)) 
      {
	fprintf(stderr,"WARNING in getPosFPtr: trying to get spherical cordinates from a cartesian data file.\n");
	return NULL;
      }
    */
    if (spherical)
      {
	if (pos_s==NULL) setSphericalFromCartesian();
	return pos_s;
      }
    else
      {
	if (pos==NULL) setCartesianFromSpherical();
	return pos;
      }
    
  }

  void subSample(double factor)
  {
    static const long long RDMX = RAND_MAX;
    static const long long RDMX_LL = ((long long)RAND_MAX)*((long long)RAND_MAX);
    long newNPart=npart*factor;
    
    for (long i=newNPart; i<npart; ++i)
      {	
	double r=double((long long)rand()+RDMX*rand())/(RDMX_LL);
	if (r < double(newNPart)/double(i+1))
	  {
	    long index=((long long)rand()+RDMX*rand())%newNPart;	    
	    std::copy(&pos[ndims*i],&pos[ndims*(i+1)],&pos[index*ndims]);
	    if (vel!=NULL) std::copy(&vel[ndims*i],&vel[ndims*(i+1)],&vel[index*ndims]);
	    if (mass!=NULL) mass[index]=mass[i];
	  }
      }    
    
    npart=newNPart;
  }

  float *getVelFPtr() {if (loaded) return vel; else return NULL;}

  float *getMassFPtr() {if (loaded) return mass; else return NULL;}

  void getBBox(double **x0p, double **deltap, bool spherical=false) 
  {
    if (loaded==false) return;

    if ((*x0p)==NULL) *x0p = (double *)malloc(sizeof(double)*ndims);
    if ((*deltap)==NULL) *deltap = (double *)malloc(sizeof(double)*ndims);
	    
    if (spherical)
      {
	/*
	  if (coordSystem!=cs_spherical)
	  {
	  fprintf(stderr,"ERROR in getBBox: coord system is not spherical.");
	  return;
	  }
	*/
	if (pos_s==NULL) setSphericalFromCartesian();
	if (x0s.size()==0) guessBBox(x0s,deltas,true);

	int i;
	for (i=0;i<ndims;i++)
	  {
	    (*x0p)[i]=x0s[i];
	    (*deltap)[i]=deltas[i];
	  }
	return;
      }
    else
      {
	if (pos==NULL) setCartesianFromSpherical();
	if (x0.size()==0) guessBBox(x0,delta,false);

	int i;
	for (i=0;i<ndims;i++)
	  {
	    (*x0p)[i]=x0[i];
	    (*deltap)[i]=delta[i];
	  }
      }
  }

};

#endif
