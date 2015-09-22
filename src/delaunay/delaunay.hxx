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
#ifndef __MY_OWN_DELAUNAY_HXX__
#define __MY_OWN_DELAUNAY_HXX__

#include <stdio.h>
#include <math.h>
#include <vector>
#include <CGAL/Random.h>
//#include "NDsubNetwork.hxx"
#include "network_interface.hxx"

#include "delaunay_setup.h"
#include "vertex_info.hxx"
#include "sampledDataInput.hxx"
#include "boundary.hxx"
#include "box.hxx"

#include "simplex.h"
#include "dtfi.h"
#include "NDnetwork.h"
#include "NDnetwork_tags.h"
#include "healpix.h"
#include "global.h"
#include "mytypes.h"

#include "NDnet_interface.hxx"
#include "find_unordered_maps.hxx"

using namespace std;

template < class Tr = Dt>
class MyDelaunay
  : public Tr
{
private:
    double MinDistToBBox_square(double *pt);

  template <class T>
  int IntersectionFlags(double *pt, double r, 
			std::vector<T> &x0, std::vector<T> &delta, 
			double borderSize=0);

    int periodicity;
    bool valueIsSet;

    vector<double> x0;
    vector<double> delta;

public: 
    typedef Tr                                   Tr_Base;
    typedef typename Tr_Base::Geom_traits        Geom_traits;
    typedef typename Tr_Base::Point              Point;
    typedef typename Tr_Base::Vertex_handle      Vertex_handle; 

  typedef typename Tr_Base::Locate_type Locate_type;
    
#if PERIODICITY==0

    typedef typename Tr_Base::Finite_vertices_iterator Finite_vertices_iterator;   
    typedef typename Tr_Base::Finite_edges_iterator    Finite_edges_iterator;
    
#if NDIMS==2 
    typedef typename Tr_Base::Finite_faces_iterator FCell_it;
  typedef typename Tr_Base::Face_handle Face_handle;
#else
    typedef typename Tr_Base::Cell_handle        Cell_handle; 
    typedef typename Tr_Base::Finite_cells_iterator    Finite_cells_iterator;    
    typedef typename Tr_Base::Finite_facets_iterator   Finite_facets_iterator;
    typedef typename Tr_Base::Finite_cells_iterator  FCell_it; 
    typedef typename Tr_Base::Cell_circulator Cell_circulator;
#endif

#else // PERIODICITIY==1

    typedef typename Tr_Base::Vertex_iterator Finite_vertices_iterator;   
    typedef typename Tr_Base::Edge_iterator    Finite_edges_iterator;
    typedef typename Tr_Base::Cell_iterator    Finite_cells_iterator;    
    typedef typename Tr_Base::Facet_iterator   Finite_facets_iterator;
    typedef typename Tr_Base::Cell_iterator  FCell_it; 
    typedef typename Tr_Base::Iso_cuboid Iso_cuboid;

#endif // PERIODICITIY

  using Tr_Base::number_of_vertices;
  using Tr_Base::geom_traits;
    
  bool SetBBox(const std::vector<double> &x0_p, const std::vector<double> &delta_p);
  bool SetBBox(double *x0_p, double *delta_p, int n);

  std::pair<int,double> CheckBoundaries(double borderSize); 
  void TagBoundaries(); 

  template < class OutputIterator> 
  int GenerateSmoothBoundaries(double borderSize, OutputIterator out);

  /*   
       int CheckBoundaries() {return CheckBoundaries(std::vector<double>(),std::vector<double>(),0);}
  */
  double CheckGuards();
  template <class T> int Mangle(T* coord, const char *fname);
  template < class OutputIterator> 
  int applyAngMask(const char *fname, OutputIterator out, double maxHoleSize=-1, int nSamples=5);
  NDnetwork *ToNDnetwork(int periodic, bool allFaces=true);
  void correctNofD(double A, double Dr, double B);
  bool SetValueToDensity();
  void applyMask(char *mask);

  template <class T> 
  INT insert(const T *pos, INT N, float *mass=NULL, double *val=NULL, long *index=NULL, INT *true_index=NULL);

  template <class inputIterator> 
  INT insert(inputIterator p_begin, inputIterator p_end, float *mass=NULL, double *val=NULL, long *index=NULL, INT *true_index=NULL); 
   
#if PERIODICITY!=0
  
    MyDelaunay(vector<double> x0=vector<double>(), vector<double> delta=vector<double>(),Iso_cuboid cub=Iso_cuboid(), int period = ~0,const Geom_traits& traits = Geom_traits())     
	       
	:Tr_Base(cub,traits)
	{
	    periodicity=period;
	    SetBBox(x0,delta);
	    valueIsSet=false;
	}

    MyDelaunay(double *x0=NULL,double *delta=NULL,Iso_cuboid cub=Iso_cuboid(), int ndims=NDIMS,
	       int period = ~0,const Geom_traits& traits = Geom_traits()) 
	:Tr_Base(cub,traits)
	{
	    periodicity=period;
	    SetBBox(x0,delta,ndims);
	    valueIsSet=false;
	}    

#else    

    MyDelaunay(vector<double> x0=vector<double>(), vector<double> delta=vector<double>(),
	       int period = 0,const Geom_traits& traits = Geom_traits()) 
	:Tr_Base(traits)
	{
	    periodicity=period;
	    SetBBox(x0,delta);
	    valueIsSet=false;
	}

    MyDelaunay(double *x0=NULL,double *delta=NULL,
	       int period = 0,const Geom_traits& traits = Geom_traits()) 
	:Tr_Base(traits)
	{
	    periodicity=period;
	    SetBBox(x0,delta,NDIMS);
	    valueIsSet=false;
	}    
#endif
};

template < class Tr >
bool MyDelaunay<Tr>::
SetBBox(const vector<double> &x0_p, const vector<double> &delta_p)
{
    int i;

    if ((x0_p.size()>NDIMS)||(delta_p.size()>NDIMS)||
	(x0_p.size()!=delta_p.size())) 
	return false;
    
    x0=vector<double>(x0_p.begin(),x0_p.end());
    delta=vector<double>(delta_p.begin(),delta_p.end());
    for (i=x0.size();i<NDIMS;i++) {x0.push_back(0);delta.push_back(0);}
    
    return true;
}

template < class Tr >
bool MyDelaunay<Tr>::
SetBBox(double *x0_p,double *delta_p, int n)
{
    int i;

    if (n>NDIMS) 
    {
	fprintf (stderr,"ERROR in SetBBox : ndims>%d\n",NDIMS);
	exit(0);
    }
    for (i=0;i<n;i++) {x0.push_back(x0_p[i]);delta.push_back(delta_p[i]);}
    for (i=n;i<NDIMS;i++) {x0.push_back(0);delta.push_back(0);}
    
    return true;
}

template < class Tr >
double MyDelaunay<Tr>::
MinDistToBBox_square(double *pt)
{
    double d=0;
    double t;
    
    for (int i=0;i<NDIMS;i++)
    {
	t=(pt[i]-x0[i]);
	if ((t<x0[i])||(t>x0[i]+delta[i])) return 0;
	if (t>0.5*delta[i]) t=delta[i]-t;
	d+=t*t;
    }

    return d;
}

template < class Tr >
template < class T >
int MyDelaunay<Tr>::
IntersectionFlags(double *pt, double r,std::vector<T> &x0p, std::vector<T> &deltap, double borderSize)
{
    int flg=0;
    if (borderSize != 0)
    {
	for (int i=0;i<NDIMS;i++)
	{
	    if (pt[i]+r >= x0p[i]+deltap[i]+borderSize) flg |= (1<<(2*i+1));
	    if (pt[i]-r < x0p[i]-borderSize) flg |= (1<<(2*i));
	}
    }
    else
    {
	for (int i=0;i<NDIMS;i++)
	{
	    if (pt[i]+r >= x0p[i]+deltap[i]) flg |= (1<<(2*i+1));
	    if (pt[i]-r < x0p[i]) flg |= (1<<(2*i));
	}
    }

    return flg;
}

template < class Tr> 
template < class inputIterator> 
INT MyDelaunay<Tr>::
insert(inputIterator p_begin, inputIterator p_end, float *mass, double *val, long *index, INT *true_index)
{
  long i,j,k;
  vertex_info *nfo;
  Vertex_handle v= Vertex_handle();
  long N;

  j=this->number_of_vertices();   

  std::map<Point,INT> pm;
  
  //printf ("(Pre) ");fflush(0);
  i=0;
  for (inputIterator p = p_begin; p!=p_end; p++) 
    {
      pm.insert(std::make_pair(*p,i));
      i++;
    }
  N=i;

  //printf ("(Tess) ");fflush(0);
#if PERIODICITY == 1
    long ni = Tr_Base::insert(p_begin,p_end,true);
#else
    long ni = Tr_Base::insert(p_begin,p_end);
#endif

    //printf ("(Post) ");fflush(0);
    
    long Nadded=0;
    for (Finite_vertices_iterator vertex=this->finite_vertices_begin();vertex!=this->finite_vertices_end();vertex++)
    {
	typename std::map<Point,INT>::iterator result;
    	if ((result=pm.find(vertex->point())) != pm.end())
	{
	    nfo = &vertex->info();

	    if (index==NULL)
	      nfo->index=result->second + j;
	    else 
	      nfo->index=index[result->second];

#if PERIODICITY != 1
	    nfo->boundary_flags=0;
	    if (true_index==NULL)
	      nfo->true_index=nfo->index;
	    else
	      nfo->true_index=true_index[result->second];
#endif
	    if (true_index==NULL) 
	      i=nfo->index;
	    else
	      i=(long)true_index[result->second];
	    //i=(long)result->second;printf("i=%ld\n",i);
	    //i=nfo->index;printf("i=%ld\n",i);
	    if (val==NULL) nfo->value=-1;
	    else nfo->value=val[i];
	    if (mass==NULL) nfo->mass=1.;
	    else nfo->mass=mass[i];
	    Nadded++;
	}
    }
       
    pm.clear();

    if (val!=NULL) valueIsSet=true;

    if (N+j!=this->number_of_vertices())
      {
	printf("\nWARNING: Only %ld/%ld point are unique.\n",this->number_of_vertices()-j,(long)N);
	printf("         Some points have identical coordinates !\n");
      }

    return N;
}

template < class Tr> 
template < class T> 
INT MyDelaunay<Tr>::
insert(const T *pos_p, INT N, float *mass, double *val, long *index, INT *true_index)
{
    long i,j,k;
    vertex_info *nfo;
    Vertex_handle v= Vertex_handle();

    j=this->number_of_vertices();    

    std::vector<Point> pv;
    std::map<Point,INT> pm;
    pv.resize(N);
    printf ("(Pr) ");fflush(0);
    for (i=0;i<N;i++) 
    {
#if NDIMS==3
	pv[i]=Point(pos_p[NDIMS*i],pos_p[NDIMS*i+1],pos_p[NDIMS*i+2]);
#elif NDIMS==2 
#if PERIODICITY == 1
	pv[i]=Point(pos_p[NDIMS*i],pos_p[NDIMS*i+1],0);
#else
	pv[i]=Point(pos_p[NDIMS*i],pos_p[NDIMS*i+1]);
#endif
#else
	fprintf (stderr,"NOT IMPLEMENTED in %d-D",NDIMS);
#endif
	pm.insert(std::make_pair(pv[i],i));
    }

    printf ("(T) ");fflush(0);
#if PERIODICITY == 1
    long ni = Tr_Base::insert(pv.begin(),pv.end(),true);
#else
    long ni = Tr_Base::insert(pv.begin(),pv.end());
#endif
    printf ("(Po) ");fflush(0);
    pv.clear();
    long Nadded=0;
    for (Finite_vertices_iterator vertex=this->finite_vertices_begin();vertex!=this->finite_vertices_end();vertex++)
    {
	typename std::map<Point,INT>::iterator result;
    	if ((result=pm.find(vertex->point())) != pm.end())
	{
	    nfo = &vertex->info();

	    if (index==NULL)
	      nfo->index=result->second + j;
	    else 
	      nfo->index=index[result->second];

#if PERIODICITY != 1
	    nfo->boundary_flags=0;
	    if (true_index==NULL)
	      nfo->true_index=nfo->index;
	    else
	      nfo->true_index=true_index[result->second];
#endif
	    i=(long)result->second;
	    if (val==NULL) nfo->value=-1;
	    else nfo->value=val[i];
	    if (mass==NULL) nfo->mass=1.;
	    else nfo->mass=mass[i];
	    Nadded++;
	}
    }
       
    pm.clear();

    if (val!=NULL) valueIsSet=true;

    if (N+j!=this->number_of_vertices())
      {
	printf("WARNING: Only %ld/%ld point are unique.\n",this->number_of_vertices()-j,(long)N);
	printf("         This means that some points have identical coordinates.\n");
      }

    return N;
}

#if PERIODICITY==0
template <class Tr>
template < class OutputIterator> 
int MyDelaunay<Tr>::
applyAngMask(const char *fname, OutputIterator out, double maxHoleSize, int nSamples)
{
  double *hmap;
  long nside;
  char coordsys[10];
  char ordering[10];
  bool nested=true;
  int i,j;
  double delta = 1./(nSamples+1);
  double x[NDIMS+1];

  double v[NDIMS+1][NDIMS];
  double *pv[NDIMS+1];
  double r2;
  double center[NDIMS];
  double maxAngSize;
  int Ninside=0;
  long ipix;
  std::vector<Point> outPoints;
  /*
  NDnetwork *net = this->ToNDnetwork();
  NDsubNetwork sub(net);
  std::vector< std::pair<char,uint> > net_cell(4);
  */  
    
  if (maxHoleSize>0) maxAngSize=tan(maxHoleSize*DEG2RAD);
  else maxAngSize=0;

  for (i=0;i<NDIMS+1;i++) pv[i]=&v[i][0];
  printf ("Applying angular mask :\n");
  /*
  printf("  Preparing ... ");fflush(0);

  for (Finite_vertices_iterator vertex=this->finite_vertices_begin();
       vertex!=this->finite_vertices_end();vertex++)
    vertex->info().boundary_flags &= (~BOUNDARY_FLAG_DUMMY);

  printf("done.\n");
  */
  
  #ifdef HAVE_CFITS_IO
  printf ("  Reading HEALPIX map from file %s ...",fname);fflush(0);
  hmap = read_healpix_map_d(fname, &nside, coordsys, ordering);
  
  if (strstr(ordering,"NESTED") != NULL) nested=true;
  else if (strstr(ordering,"RING") != NULL) nested =false;
  else fprintf (stderr,"Unknown ordering %s. \n",ordering);
  if (nested)
    printf (" done. (NESTED)\n");
  else
    printf (" done. (RING)\n");
  long nfound=0;

#if NDIMS==2 && PERIODICITY==0 
  int ntot = this->number_of_faces();
#elif NDIMS==3 || PERIODICITY==1
#if PERIODICITY==0
  int ntot = this->number_of_finite_cells();
#else
  int ntot = this->number_of_cells();
#endif
#endif
 
  printf ("  Checking cells ... ");fflush(0);
  
  FILE *f;FILE *g;
  f=fopen("test.dat","w");g=fopen("testO.dat","w");
  fprintf(f,"# ra dec dist\n");fprintf(g,"# ra dec dist\n");
  

  int nsofar=0;
  #if NDIMS==2  
    for (FCell_it fc=this->finite_faces_begin();fc!=this->finite_faces_end();fc++)
#elif NDIMS==3 
    for (FCell_it fc=this->finite_cells_begin();fc!=this->finite_cells_end();fc++)
#endif
    {
      Point p1,p2;
      int i1,i2;
      double dx;

      if ( (nsofar & ((1<<15)-1)) == 0)
	{printf ("\r  Checking cells ... (%.2f%%) ",100.*(double)nsofar/(double)ntot);fflush(0);}
      
      for (i1=0;i1<NDIMS+1;i1++)
	for (i2=i1+1;i2<NDIMS+1;i2++)
	  {

	    p1 = fc->vertex(i1)->point();
	    p2 = fc->vertex(i2)->point();

	    dx=delta;
	    for (i=0,dx=delta;i<nSamples;i++,dx+=delta)
	      {
		x[0]=p1.x() + dx * (p2.x()-p1.x());
		x[1]=p1.y() + dx * (p2.y()-p1.y());
#if NDIMS==3
		x[2]=p1.z() + dx * (p2.z()-p1.z());
#endif
		sampledDataInput::cartesian2spherical(x,x,NDIMS,true);
		
		long ipix;
		
		if (nested)
		  ang2pix_nest(nside,PI - (x[0]+PI/2),x[1],&ipix);
		else
		  ang2pix_ring(nside,PI - (x[0]+PI/2),x[1],&ipix);
		
		if (hmap[ipix]<=0.75) break;
	      }
	    
	    if (i!=nSamples) {
	      i1=NDIMS+1;
	      break;
	    }
	  }
       
      //one edge at least is out
      if (i2!=NDIMS+1)
      	{
	  double rho_avg;
	  double rho[NDIMS+1];
	  
	  nfound++;
	  
	  rho_avg=0;
	  for (i=0;i<NDIMS+1;i++)
	    {
	      v[i][0] = CGAL::to_double(fc->vertex(i)->point().x());
	      v[i][1] = CGAL::to_double(fc->vertex(i)->point().y());
#if NDIMS==3
	      v[i][2] = CGAL::to_double(fc->vertex(i)->point().z());
#endif
	      rho[i] =   fc->vertex(i)->info().value;
	      rho_avg += rho[i];
	    }
	  rho_avg /= NDIMS+1;
	  
	  bool isInside=false;
	  double volume=SimplexVolumed(pv,NDIMS+1,NDIMS);
	  int ngen= randomPoisson(volume*rho_avg);
	  double *result=NULL;	  
	  
	  if (maxHoleSize>0) {
	    double dist2 = 0;
	    r2 = SimplexSphere(pv,NDIMS,center);
	    for (i=0;i<NDIMS;i++) dist2+=center[i]*center[i];
	    isInside=(sqrt(r2/dist2)<maxAngSize);
	  }
	  
	  genDTFIInSimplex(pv,rho,ngen,NDIMS,&result);

	  for (i=0;i<ngen*NDIMS;i+=NDIMS) {
	    sampledDataInput::cartesian2spherical(&result[i],x,NDIMS,true);
	    if (nested)
	      ang2pix_nest(nside,PI - (x[0]+PI/2),x[1],&ipix);
	    else
	      ang2pix_ring(nside,PI - (x[0]+PI/2),x[1],&ipix);
	    
	    if (hmap[ipix] < ((double)random())/RAND_MAX) 
	      {
		Point p;
#if NDIMS==3
		p=Point(result[i],result[i+1],result[i+2]);
#else
		p=Point(result[i],result[i+1]);
#endif
		if (isInside) {
		  *out=p;
		  Ninside++;		  
		}
		else outPoints.push_back(p);
		if (isInside) fprintf(f,"%g %g %g\n",x[0]*RAD2DEG,x[1]*RAD2DEG,x[2]);
		else fprintf(g,"%g %g %g\n",x[0]*RAD2DEG,x[1]*RAD2DEG,x[2]);
		
	      }

	  }

	  free(result);
	} // finished generating dist in this cell
      nsofar++;
    } // finished checking all cells
    
    fclose(f);fclose(g);
  for (i=0;i<outPoints.size();i++) *out = outPoints[i];
  printf ("\r  Checking cells ... done. (~%.2f%% outside mask)\n",100.*(double)nfound/(double)ntot);
 
  #else
  fprintf(stderr,"ERROR: library CFitsio was not found.\n");
  fprintf(stderr,"   I need it to read angular mask. GO DOWNLOAD IT !!!!! \n");
  fprintf(stderr,"   And don't forget to recompile or I'll complain again.\n");
  exit(-1);
  #endif


  return Ninside;
}

template < class Tr> 
template < class T> 
int MyDelaunay<Tr>::
Mangle(T* coord, const char *mask_fname)
{
  FILE *fi;
  FILE *fo;
  int i,j;

  char tmpin[255];
  char tmpout[255];
  char cmd[255];

  double r2;
  double v[NDIMS+1][NDIMS];
  double *pv[NDIMS+1];
  double center[NDIMS];
  double center_sph[NDIMS];
  int ntrue=0;
  vertex_info *nfo;
  double pi = 3.141592653589793116;

  //printf ("Mangling: ");fflush(0);
  printf ("Preparing input file for drangle ...");fflush(0);
  for (int i=0;i<NDIMS+1;i++) pv[i]=&v[i][0];

  strcpy(tmpin,"drangle.tmp.input");
  strcpy(tmpout,"drangle.tmp.output");
  sprintf(cmd,"drangle -ud,r %s %s %s",mask_fname,tmpin,tmpout);
  fi=fopen(tmpin,"w");

#if NDIMS==2  
  for (FCell_it fc=this->finite_faces_begin();fc!=this->finite_faces_end();fc++)
#elif NDIMS==3 
    for (FCell_it fc=this->finite_cells_begin();fc!=this->finite_cells_end();fc++)
#endif
      {
	ntrue=0;
       	for (int i=0;i<NDIMS+1;i++)
	  {
	    nfo=&fc->vertex(i)->info();
	    if (nfo->index==nfo->true_index) ntrue++;

	    v[i][0] = CGAL::to_double(fc->vertex(i)->point().x());
	    v[i][1] = CGAL::to_double(fc->vertex(i)->point().y());
#if NDIMS==3
	    v[i][2] = CGAL::to_double(fc->vertex(i)->point().z());
#endif
	  }
	if (!ntrue) continue;
	
	//double pavg[NDIMS];
	//for (j=0;j<NDIMS;j++) pavg[j]=0;
	int nsuccess=0;
	for (int i=0;i<NDIMS+1;i++)
	  {
	    nfo=&fc->vertex(i)->info();
	    if (nfo->index!=nfo->true_index) continue;
	    
	    INT id = nfo->index;
	    //for (j=0;j<NDIMS;j++)
	    //pavg[j]+=coord[id*NDIMS+j];
	  }
	//for (j=0;j<NDIMS;j++) pavg[j]/=ntrue;
	      
	double angle;
	r2 = SimplexSphere(pv,NDIMS,center);
	
#if NDIMS==2  
	angle=sqrt(r2);
#elif NDIMS==3 
	sampledDataInput::cartesian2spherical(center,center,NDIMS);
	angle=atan2(sqrt(r2),center[2])*180./pi;
#endif
	fprintf(fi,"%e %e %e\n",center[0],center[1],fabs(angle));
      } 
  fclose(fi);
  printf("done.\n");
  printf("Executing cmd : %s\n",cmd);
  if (!system(cmd))
    {
      fprintf(stderr,"ERROR using mangle. \n");
      exit(0);
    }
  //remove(fi);
  //fclose(fo);
  
}

template < class Tr >
double MyDelaunay<Tr>::
CheckGuards()
{
  vertex_info *nfo;
  long nfail=0;
  long ntot=0;
  double v[NDIMS];
  int i,j;

#if NDIMS==2 && PERIODICITY==0 
    for (FCell_it fc=this->finite_faces_begin();fc!=this->finite_faces_end();fc++)
#elif NDIMS==3 || PERIODICITY==1
    for (FCell_it fc=this->finite_cells_begin();fc!=this->finite_cells_end();fc++)
#endif
    {
      int nguard=0;
      int nout=0;
      for (i=0;i<NDIMS+1;i++)
	{
	  nfo = &fc->vertex(i)->info();
	  if (nfo->true_index == NEWPARTICLE_OUTSIDE_INDEX) nguard++;
	}
      if (nguard)
	{
	  for (i=0;i<NDIMS+1;i++)
	    {
	      v[0] = CGAL::to_double(fc->vertex(i)->point().x());
	      v[1] = CGAL::to_double(fc->vertex(i)->point().y());
	  
#if NDIMS==3
	      v[2] = CGAL::to_double(fc->vertex(i)->point().z());
#endif
	      for (j=0;j<NDIMS;j++)
		{		  
		  if ((v[j]<x0[j])||(v[j]>x0[j]+delta[j])) break;
		}
	      if (j!=NDIMS) nout++;
	    }
	  
	  if (nout!=NDIMS+1) nfail++;
	}
      
      ntot++;
    }

    return (double)nfail/(double)ntot;
}

template < class Tr >
template < class OutputIterator> 
int MyDelaunay<Tr>::
GenerateSmoothBoundaries(double borderSize, OutputIterator out)
{  
  double v[NDIMS+1][NDIMS];
  double *pv[NDIMS+1];

  std::vector<float> bx;
  std::vector<Point> pts;
  std::vector<long> index;
  std::vector<INT> true_index;
  std::vector<double> x0v;
  std::vector<double> deltav;
  double fmax=1.E13;
  int i,j;
  long ngentot=0;
  

  //pts.clear();
  //index.clear();
  //true_index.clear();
  x0v.assign(x0.begin(),x0.end());
  deltav.assign(delta.begin(),delta.end());
  for (i=0;i<NDIMS;i++)
    {
      //printf("bdbd : %lg\n",borderSize/deltav[i]);
      if (fmax>borderSize/deltav[i]) fmax = borderSize/deltav[i];
      x0v[i]-=borderSize;
      deltav[i]+=2*borderSize;
    }
  for (i=0;i<NDIMS+1;i++) pv[i]=&v[i][0];

  //printf ("fmax=%f\n",fmax);
  if (fmax>0.01)
    generateGuardBox(borderSize,x0v,deltav,std::back_inserter(bx));
  else generateBox(100,x0v,deltav,std::back_inserter(bx),true);

  printf ("Building smooth boundaries (%ld guards): ",bx.size()/NDIMS);fflush(0);

      for (j=0;j<bx.size()/NDIMS;j++)
      {
#if NDIMS==3
	pts.push_back(Point(bx[j*NDIMS],bx[j*NDIMS+1],bx[j*NDIMS+2]));
#else
	pts.push_back(Point(bx[j*NDIMS+0],bx[j*NDIMS+1]));
#endif
	true_index.push_back(NEWPARTICLE_OUTSIDE_INDEX);
	index.push_back(NEWPARTICLE_OUTSIDE_INDEX);
      }
      printf (" tesselation ... ");fflush(0);
      insert(pts.begin(),pts.end(),NULL,NULL,&(index[0]),&(true_index[0]));
      printf (" generating ... ");fflush(0);
      
      
	FILE *f;
	f=fopen("test_smooth.dat","w");
	fprintf(f,"# px py pz\n");
      

#if NDIMS==2  
      for (FCell_it fc=this->finite_faces_begin();fc!=this->finite_faces_end();fc++)
#elif NDIMS==3 
	for (FCell_it fc=this->finite_cells_begin();fc!=this->finite_cells_end();fc++)
#endif
	  {
	    double rho_avg;
	    double rho[NDIMS+1];

	    for (i=0;i<NDIMS+1;i++)
	      if (fc->vertex(i)->info().index == NEWPARTICLE_OUTSIDE_INDEX) break;
	    if (i==NDIMS+1) continue;

	    rho_avg=0;
	    for (i=0;i<NDIMS+1;i++)
	      {
		v[i][0] = CGAL::to_double(fc->vertex(i)->point().x());
		v[i][1] = CGAL::to_double(fc->vertex(i)->point().y());
#if NDIMS==3
		v[i][2] = CGAL::to_double(fc->vertex(i)->point().z());
#endif
		if (fc->vertex(i)->info().index == NEWPARTICLE_OUTSIDE_INDEX)
		  rho[i]=0;
		else 
		  rho[i] =   fc->vertex(i)->info().value;
		
		rho_avg += rho[i];
	      }
	    rho_avg /= NDIMS+1;

	    if (rho_avg==0) continue;

	    double *result=NULL;
	    double volume=SimplexVolumed(pv,NDIMS+1,NDIMS);
	    int ngen =  randomPoisson(volume*rho_avg);
	    ngentot+=ngen;
	    
	    //if (ngen>5000) printf("\n WARNING:Will gen %d particles\n",ngen);

	    genDTFIInSimplex(pv,rho,ngen,NDIMS,&result);
	    for (i=0;i<ngen*NDIMS;i+=NDIMS) {
#if NDIMS==3
	      *out=Point(result[i],result[i+1],result[i+2]);
#else
	      *out=Point(result[i],result[i+1]);
#endif
	      fprintf(f,"%g %g %g\n",result[0],result[1],result[2]); 
	    }
	    free(result);
	  }
      fclose(f);

      
      //printf("(%d new)",nin);fflush(0);
      printf ("cleaning up ... ");fflush(0);
      
      for (Finite_vertices_iterator vertex=this->finite_vertices_begin();vertex!=this->finite_vertices_end();vertex++)
	{
	  if (vertex->info().index==NEWPARTICLE_OUTSIDE_INDEX)
	    Tr_Base::remove(vertex);
	}
      

      printf("done. (+%ld particles)\n",ngentot);
      return ngentot;
}


template < class Tr >
void MyDelaunay<Tr>::
TagBoundaries()
{
  double v[NDIMS+1][NDIMS];
  double *pv[NDIMS+1];
  int bflags[NDIMS+1];
  long i;
   
  for (i=0;i<NDIMS+1;i++) pv[i]=&v[i][0];

  for (Finite_vertices_iterator vertex=this->finite_vertices_begin();vertex!=this->finite_vertices_end();vertex++)
    {
      vertex->info().boundary_flags=0;
    }

#if NDIMS==2  
  for (FCell_it fc=this->finite_faces_begin();fc!=this->finite_faces_end();fc++)
#elif NDIMS==3 
    for (FCell_it fc=this->finite_cells_begin();fc!=this->finite_cells_end();fc++)
#endif
      {
       	for (i=0;i<NDIMS+1;i++)
	  {
	    v[i][0] = CGAL::to_double(fc->vertex(i)->point().x());
	    v[i][1] = CGAL::to_double(fc->vertex(i)->point().y());

#if NDIMS==3
	    v[i][2] = CGAL::to_double(fc->vertex(i)->point().z());
#endif
	    bflags[i]=IntersectionFlags(v[i],0,x0,delta);
	    fc->vertex(i)->info().boundary_flags |= bflags[i];
	    
	    if (fc->vertex(i)->info().true_index == NEWPARTICLE_OUTSIDE_INDEX)
	      fc->vertex(i)->info().boundary_flags |= ~0; 
	  }
      }
}

template < class Tr >
 std::pair<int,double> MyDelaunay<Tr>::
CheckBoundaries(double borderSize)
{
    double v[NDIMS+1][NDIMS];
    double *pv[NDIMS+1];
    double center[NDIMS];
    int bflags[NDIMS+1];
    int boundary_flags;
    double r2;
    double dist;
    //int ncells=0;
    long i;
    long nfailed=0;
    double suggested_size=0;
    
    for (i=0;i<NDIMS+1;i++) pv[i]=&v[i][0];

    for (Finite_vertices_iterator vertex=this->finite_vertices_begin();vertex!=this->finite_vertices_end();vertex++)
      {
	vertex->info().boundary_flags=0;
      }

#if NDIMS==2  
    for (FCell_it fc=this->finite_faces_begin();fc!=this->finite_faces_end();fc++)
#elif NDIMS==3 
    for (FCell_it fc=this->finite_cells_begin();fc!=this->finite_cells_end();fc++)
#endif
      {
       	for (i=0;i<NDIMS+1;i++)
	  {
	    v[i][0] = CGAL::to_double(fc->vertex(i)->point().x());
	    v[i][1] = CGAL::to_double(fc->vertex(i)->point().y());

#if NDIMS==3
	    v[i][2] = CGAL::to_double(fc->vertex(i)->point().z());
#endif
	    bflags[i]=IntersectionFlags(v[i],0,x0,delta);
	    fc->vertex(i)->info().boundary_flags |= bflags[i];
	    
	    if (fc->vertex(i)->info().true_index == NEWPARTICLE_OUTSIDE_INDEX)
	      fc->vertex(i)->info().boundary_flags |= ~0; 
	  }
	
	for (i=0;i<NDIMS+1;i++)
	  if (fc->vertex(i)->info().unsafe!=-1) break;
	if (i==NDIMS+1) continue;
	
	r2 = SimplexSphere(pv,NDIMS,center);
	
	for (i=0;i<NDIMS+1;i++) if (!bflags[i]) break;
	if ((i!=NDIMS+1)&&(borderSize>0)) 
	  {
	    boundary_flags  = IntersectionFlags(center,sqrt(r2),x0,delta,borderSize);
	    if (boundary_flags) 
	      {
		for (i=0;i<NDIMS;i++) 
		  fc->vertex(i)->info().unsafe=1;
		nfailed++;	  
	    	    
		for (i=0;i<NDIMS;i++)
		  {
		    if (x0[i] - suggested_size > center[i]-sqrt(r2)) 
		      suggested_size = x0[i] - (center[i]-sqrt(r2));
		    if (x0[i] + delta[i] + suggested_size < center[i]+sqrt(r2)) 
		      suggested_size = (center[i]+sqrt(r2)) - (x0[i] + delta[i]);
		  }
	      }
	}
    }
    
    for (Finite_vertices_iterator vertex=this->finite_vertices_begin();vertex!=this->finite_vertices_end();vertex++)
      {
	if (vertex->info().unsafe==0) vertex->info().unsafe=-1;
      }

    return std::make_pair((int)nfailed,suggested_size);
}

template < class Tr >
void MyDelaunay<Tr>::
applyMask(char *mask)
{
  vertex_info *nfo;

  for (Finite_vertices_iterator vertex=this->finite_vertices_begin();vertex!=this->finite_vertices_end();vertex++)
    {
      nfo=&vertex->info();
      if (nfo->index != nfo->true_index) continue; 
      if (!mask[nfo->index]) nfo->boundary_flags = ~0;
    }
}


#else // PERIODICITY==1
template < class Tr>
 template <class T>
int MyDelaunay<Tr>::
applyAngMask(const char *fname, OutputIterator out, double maxHoleSize, int nSamples)
{
  return 0;
}

template < class Tr> 
template < class T> 
int MyDelaunay<Tr>::
Mangle(T* coord, const char *mask_fname)
{
  return 0;
}

template < class Tr >
double MyDelaunay<Tr>::
CheckGuards()
{
  return 0;
}

template < class Tr >
template < class OutputIterator> 
int MyDelaunay<Tr>::
GenerateSmoothBoundaries(double borderSize, OutputIterator out)
{
  return 0;
}
/*
template < class Tr >
template < class T> 
std::pair<int,double>  MyDelaunay<Tr>::
CheckBoundaries(std::vector<T> x0, std::vector<T> delta,double borderSize)
{
  return std::pair<int,double>();
}
*/
template < class Tr >
void MyDelaunay<Tr>::
applyMask(char *mask)
{
  //return ;
}

#endif // PERIODICITY

inline double NofD_eval(double d, double A, double Dr, double B)
{
  return A*d*d*exp(-pow(d/Dr,B));
}

template < class Tr >
void MyDelaunay<Tr>::
correctNofD(double A, double Dr, double B)
{
  double d;
  double dmax;
  double maxval_inv;

  dmax = Dr * pow(2./B,1./B);
  maxval_inv = 1./NofD_eval(dmax,A,Dr,B);
  
  for (Finite_vertices_iterator vertex=this->finite_vertices_begin();vertex!=this->finite_vertices_end();vertex++)
    {
      d=vertex->point().x()*vertex->point().x()+vertex->point().y()*vertex->point().y();
#if NDIMS==3
      d+=vertex->point().z()*vertex->point().z();
#endif
      d=sqrt(d);
      vertex->info().value /= (NofD_eval(d,A,Dr,B) * maxval_inv);
      //printf ("N(%e) = %e\n",d,NofD_eval(d,A,Dr,B) * maxval_inv);
    }
}

template < class Tr >
bool MyDelaunay<Tr>::
SetValueToDensity()
{
    double v[NDIMS+1][NDIMS];
    double *pv[NDIMS+1];

    //printf ("Computing density ... ");fflush(0);

    for (int i=0;i<NDIMS+1;i++) pv[i]=&v[i][0];
#if PERIODICITY==0    
    for (Finite_vertices_iterator vertex=this->finite_vertices_begin();vertex!=this->finite_vertices_end();vertex++)
    {
	vertex->info().value = 0;	
    }
#else
    for (Finite_vertices_iterator vertex=this->vertices_begin();vertex!=this->vertices_end();vertex++)
    {
	vertex->info().value = 0;	
    }
#endif
    //FILE *f=fopen("test.dat","w");
#if NDIMS==2 && PERIODICITY==0 
    for (FCell_it fc=this->finite_faces_begin();fc!=this->finite_faces_end();fc++)
#elif NDIMS==3 || PERIODICITY==1
    for (FCell_it fc=this->finite_cells_begin();fc!=this->finite_cells_end();fc++)
#endif
    {
	double vol;
	for (int i=0;i<NDIMS+1;i++)
	    {
		Point p =fc->vertex(i)->point();
		v[i][0] = CGAL::to_double(p.x());
		v[i][1] = CGAL::to_double(p.y());
#if NDIMS==3
		v[i][2] = CGAL::to_double(p.z());	
#endif	
	    }

#if PERIODICITY==1

	int region=0;
	for (int i=0;i<NDIMS+1;i++)
	{
	    for (int j=0;j<NDIMS;j++)
	    {
		if (fabs(v[i][j]-x0[j])>2*delta[j]/3) region|=(1<<(2*j+1));
		if (fabs(v[i][j]-x0[j])<delta[j]/3) region|=(1<<(2*j));
	    }
	}
	for (int j=0;j<NDIMS;j++)
	{
	    if (region & ((1<<(2*j+1)) | (1<<(2*j))))
	    {
		for (int i=0;i<NDIMS+1;i++)
		{
		    if (fabs(v[i][j]-x0[j])>delta[j]/2) v[i][j]-=delta[j];
		}
	    }
	}
	
#endif   

	vol=SimplexVolumed(pv,NDIMS+1,NDIMS);
	for (int i=0;i<NDIMS+1;i++)
	{
	    fc->vertex(i)->info().value+=vol;
	}
	//fprintf (f,"%d %e\n",vertex->info().index,vertex->info().value);
    }

   

#if PERIODICITY==0    
    for (Finite_vertices_iterator vertex=this->finite_vertices_begin();vertex!=this->finite_vertices_end();vertex++)
#else
    for (Finite_vertices_iterator vertex=this->vertices_begin();vertex!=this->vertices_end();vertex++)
#endif
      //for (Finite_vertices_iterator vertex=this->finite_vertices_begin();vertex!=this->finite_vertices_end();vertex++)
      {
    	vertex->info().value = vertex->info().mass * (NDIMS+1)/vertex->info().value;
	//fprintf (f,"%d %e\n",vertex->info().index,vertex->info().value);
    	//vertex->info().value = log(vertex->info().value);
      }
    //fclose(f);
    //printf ("done.\n");
    

    return true;
}


template < class Tr >
NDnetwork *MyDelaunay<Tr>::
ToNDnetwork(int periodic, bool allFaces)
{
    long ndims = NDIMS;
    NDnetwork *net;
    long index;
    long i,j,k;
    vertex_info nfo;
    //double *val;
    bool useIndexId;
    long useIndexMapFrom=-1;
    
    typedef my_unordered_map<long,INT>::type myMapT;
    myMapT index2id;
    
#if NDIMS==2 && PERIODICITY==0 
    long ncells = this->number_of_faces();
#elif NDIMS==3 || PERIODICITY==1
#if PERIODICITY==0
    long ncells = this->number_of_finite_cells();
#else
    long ncells = this->number_of_cells();
#endif
#endif

    if (verbose>0) printf ("Converting tesselation to NDnetwork ... ");fflush(0);
    this->TagBoundaries();
    net=CreateNetwork(ndims,0,0);
    
    strcpy(net->comment,DELAUNAY_TESSELATION_TAG);
    for (i=0;i<ndims;i++) net->x0[i]=this->x0[i];
    for (i=0;i<ndims;i++) net->delta[i]=this->delta[i];
    net->nvertex=this->number_of_vertices();
    net->isSimpComplex=1;

    net->nfaces[ndims] = ncells;
    net->haveVertexFromFace[NDIMS]=1;
    net->f_vertexIndex[ndims]=(NDNET_UINT *)malloc(sizeof(NDNET_UINT)*(size_t)net->nfaces[ndims]*(size_t)(ndims+1));
   
    NDNetFlags_enable(net,0);
    NDNetFlags_enable(net,ndims);
    
    long n_trueVert=0;
    net->v_coord = (float *)malloc(sizeof(float)*(size_t)net->nvertex*(size_t)NDIMS);
    
    long maxindex=0;
    
    for (Finite_vertices_iterator vertex=this->finite_vertices_begin();vertex!=this->finite_vertices_end();vertex++)
      {
	if (vertex->info().index>maxindex) maxindex=vertex->info().index;
      }
    //printf("maxindex : %ld; nvert : %ld\n",(long)maxindex,(long)net->nvertex);
    if (maxindex<net->nvertex) useIndexId=true;
    else 
      {
	useIndexId=false;
#ifdef HAVE_TR1
	index2id.rehash(std::ceil((double)net->nvertex / index2id.max_load_factor()));
#endif
      }
    i=-1;

    for (Finite_vertices_iterator vertex=this->finite_vertices_begin();vertex!=this->finite_vertices_end();vertex++)
    {
	Point p = vertex->point();
	nfo=vertex->info();

#if PERIODICITY!=1
	if ((periodic)&&(nfo.index!=nfo.true_index)) continue;
#endif
	
	n_trueVert++;
	
	if (useIndexId) i=nfo.index;
	else index2id.insert(std::make_pair(nfo.index,++i));
	
	net->v_coord[i*NDIMS] = p.x();
	net->v_coord[i*NDIMS+1] = p.y();
#if NDIMS==3
	net->v_coord[i*NDIMS+2] = p.z();
#endif

#if PERIODICITY!=1
	if ((!periodic)&&(nfo.boundary_flags)) net->v_flag[i]|=NDNETFLAG_OUT;
#endif
    }
    
    if (n_trueVert != net->nvertex)
      {
	net->nvertex=n_trueVert;
	net->v_flag=(unsigned char*) realloc(net->v_flag,sizeof(unsigned char)*(size_t)net->nvertex);
	net->v_coord=(float *)realloc(net->v_coord,sizeof(float)*(size_t)net->nvertex*(size_t)NDIMS);
      }
    index=0;
    //FreeNDnetwork(&net);
    std::vector<bool> keepVertex(net->nvertex,false);    
    std::vector<long> curVertexIndex(NDIMS+1);
    long n_trueFace=0;
#if NDIMS==2 && PERIODICITY!=1
    for (FCell_it fc=this->finite_faces_begin();fc!=this->finite_faces_end();fc++)
#elif NDIMS==3 || PERIODICITY==1 
    for (FCell_it fc=this->finite_cells_begin();fc!=this->finite_cells_end();fc++)
#endif
    {	
      long minid=net->nvertex;
      long minvertex=0;
      bool inside=false;

      net->f_flag[ndims][(index)/(NDIMS+1)]=0;
      
      for (i=0;i<NDIMS+1;i++)
	{
	  nfo = fc->vertex(i)->info();
	  
#if PERIODICITY!=1
	  if ((minid>nfo.true_index)&&(nfo.true_index>=0)) {minid=nfo.true_index;minvertex=i;}
	  if (periodic) curVertexIndex[i] = net->f_vertexIndex[ndims][index++] = nfo.true_index;
	  else {
	    //if (nfo.boundary_flags) net->f_flag[ndims][(index)/(NDIMS+1)] |= NDNETFLAG_OUT;
	    
	    if (useIndexId) curVertexIndex[i] = net->f_vertexIndex[ndims][index++] = nfo.index;
	    else curVertexIndex[i] = net->f_vertexIndex[ndims][index++] = index2id[nfo.index];
	  }
#else	    
	  curVertexIndex[i] = net->f_vertexIndex[ndims][index++] = nfo.index;
#endif
	  if (!(net->v_flag[curVertexIndex[i]]&NDNETFLAG_OUT)) inside=true;
	}

      if (inside) 
	{
	  for (i=0;i<NDIMS+1;i++) 
	    keepVertex[curVertexIndex[i]]=true;
	}
      else
	{
	  index-=(NDIMS+1);
	  continue;
	}

      nfo = fc->vertex(minvertex)->info();

#if PERIODICITY!=1
      if (!periodic) n_trueFace++;
      else if (nfo.index == nfo.true_index) n_trueFace++;
      else index-=(NDIMS+1);
#else
      n_trueFace++;
#endif
    }
    
    if (n_trueFace != net->nfaces[ndims])
      {
	net->nfaces[ndims]=n_trueFace;
	net->f_flag[ndims]=(unsigned char*) realloc(net->f_flag[ndims],sizeof(unsigned char)*(size_t)net->nfaces[ndims]);
	net->f_vertexIndex[ndims]=(NDNET_UINT *)realloc(net->f_vertexIndex[ndims],sizeof(NDNET_UINT)*(size_t)net->nfaces[ndims]*(size_t)(ndims+1));
      }

    
    for (i=0;i<net->nvertex;i++) 
      if (!keepVertex[i]) break;
    
    if (i!=net->nvertex)
      {
	std::vector<long> id2id(net->nvertex);
	j=0;
	for (i=0;i<net->nvertex;i++)
	  {
	    if (keepVertex[i]) 
	      {
		if (i!=j) 
		  {
		    memcpy(&net->v_coord[ndims*j],&net->v_coord[ndims*i],ndims*sizeof(float));
		    net->v_flag[j]=net->v_flag[i];
		  }
		id2id[i]=j++;
	      }
	    else 
	      {
		id2id[i]=-1;
	      }
	  }
	net->nvertex=j;
	net->v_flag=(unsigned char*) realloc(net->v_flag,sizeof(unsigned char)*(size_t)net->nvertex);
	net->v_coord=(float *)realloc(net->v_coord,sizeof(float)*(size_t)net->nvertex*(size_t)NDIMS);
	
	for (i=0;i<net->nfaces[ndims]*(NDIMS+1);i++)
	  net->f_vertexIndex[ndims][i]= id2id[net->f_vertexIndex[ndims][i]];

	if (useIndexId)
	  {
	    /*
	    i=0;
	    for (Finite_vertices_iterator vertex=this->finite_vertices_begin();vertex!=this->finite_vertices_end();vertex++)
	      {
		nfo=vertex->info();
		index2id.insert(std::make_pair(nfo.index,nfo.index));
	      }
	    */
	    useIndexId=false;
	    //long nmaps=0;
	    //printf("HELLO\n");
	    useIndexMapFrom=net->nvertex+1;
	    for (Finite_vertices_iterator vertex=this->finite_vertices_begin();vertex!=this->finite_vertices_end();vertex++)
	      {
		nfo=vertex->info();

		if ((id2id[nfo.index]>=0)&&(id2id[nfo.index]!=nfo.index))
		  index2id.insert(std::make_pair(nfo.index,id2id[nfo.index]));
		//nmaps++;
		
		if ((id2id[nfo.index]<0)||(id2id[nfo.index]!=nfo.index))
		  {
		    if (useIndexMapFrom > nfo.index)
		      useIndexMapFrom = nfo.index;
		  }
	      }
	    /*
#ifdef HAVE_TR1
	    index2id.rehash(std::ceil((double)nmaps / index2id.max_load_factor()));
#endif

	    for (Finite_vertices_iterator vertex=this->finite_vertices_begin();vertex!=this->finite_vertices_end();vertex++)
	      {
		nfo=vertex->info();
		if ((id2id[nfo.index]>=0)&&(id2id[nfo.index]!=nfo.index))
		  index2id.insert(std::make_pair(nfo.index,id2id[nfo.index]));
	      }
	    */
	    //printf("HI\n");
	  }
	else
	  {
	    for (myMapT::iterator it=index2id.begin();it!=index2id.end();)
	      {
		if (id2id[it->second]<0) 
		  {
		    myMapT::iterator newIT=it;
		    newIT++;
		    index2id.erase(it);
		    it=newIT;
		  }
		else
		  {
		    it->second = id2id[it->second];
		    it++;
		  }
	      }
	  }	
      }
    keepVertex.clear();
    // FROM HERE WE CAN PROBABLY COMMENT EVERYTHING OUT   
    if (allFaces) {
#if NDIMS==2
    /*
    ComputeVertexFromFace(net,1);    
    free(net->v_faceIndex[2]);
    free(net->v_numFaceIndexCum[2]);    
    net->haveFaceFromVertex[2]=0;
    NDNetFlags_enable(net,1);
    */
#elif NDIMS==3
    index=0;
#if PERIODICITY==0
    net->nfaces[2] = (uint)this->number_of_finite_facets();
#else
    net->nfaces[2] = (uint)this->number_of_facets();
#endif
    net->haveVertexFromFace[2]=1;
    net->f_vertexIndex[2]=(NDNET_UINT *)malloc(sizeof(NDNET_UINT)*(size_t)net->nfaces[2]*(size_t)(2+1));
    NDNetFlags_enable(net,2);
    uint n_trueFacet=0;
    index=0;
#if PERIODICITY==0
    for (Finite_facets_iterator f=this->finite_facets_begin(); f!=this->finite_facets_end();f++)
#else
      for (Finite_facets_iterator f=this->facets_begin(); f!=this->facets_end();f++)
#endif
    {
      int minid=net->nvertex;
      int minvertex=0;
      bool inside=true;

      net->f_flag[2][(index)/(2+1)]=0;
      for (i=0;i<NDIMS+1;i++)
	{
	  if (f->second ==i) continue;
	  if (!inside) {index++;continue;}

	  nfo = f->first->vertex(i)->info();
#if PERIODICITY!=1
	  if ((minid>nfo.true_index)&&(nfo.true_index>=0)) {minid=nfo.true_index;minvertex=i;}
	  if (periodic) net->f_vertexIndex[2][index++] = nfo.true_index;
	  else {
	    //if (nfo.boundary_flags) net->f_flag[2][(index)/(2+1)] |= NDNETFLAG_OUT;
	    if (useIndexId) net->f_vertexIndex[2][index++] = nfo.index;
	    else 
	      {
		INT curindex=nfo.index;
		
		if (curindex>=useIndexMapFrom)
		  {
		    myMapT::iterator it=index2id.find(curindex);
		    if (it==index2id.end()) {inside=false;index++;continue;}
		    curindex=it->second;
		  }
		net->f_vertexIndex[2][index++] = curindex;
	      }
	  }
#else	    
	  net->f_vertexIndex[2][index++] = nfo.index;
#endif
	}

      if (!inside)
	{
	  index-=2+1;
	  continue;
	}
      nfo = f->first->vertex(minvertex)->info();
      
#if PERIODICITY!=1
      if (!periodic) n_trueFacet++;
      else if (nfo.index == nfo.true_index) n_trueFacet++;
      else index-=(2+1);
#else
      n_trueFacet++;
#endif
    }
      
    if (n_trueFacet != net->nfaces[2])
      {
	net->nfaces[2]=n_trueFacet;
	
	net->f_flag[2]=(unsigned char*) realloc(net->f_flag[2],sizeof(unsigned char)*(size_t)net->nfaces[2]);
	net->f_vertexIndex[2]=(NDNET_UINT *)realloc(net->f_vertexIndex[2],sizeof(NDNET_UINT)*(size_t)net->nfaces[2]*(size_t)(2+1));
      }


    //#endif // NDIMS==3    


    index=0;
#if PERIODICITY==0
    net->nfaces[1] = (uint)this->number_of_finite_edges();
#else
    net->nfaces[1] = (uint)this->number_of_edges();
#endif
    net->haveVertexFromFace[1]=1;
    net->f_vertexIndex[1]=(NDNET_UINT *)malloc(sizeof(NDNET_UINT)*(size_t)net->nfaces[1]*(size_t)(1+1));
    NDNetFlags_enable(net,1);
    uint n_trueSeg=0;
    for (Finite_edges_iterator f=this->finite_edges_begin(); f!=this->finite_edges_end();f++)
    {
      int minid=net->nvertex;
      int minvertex=0;

      net->f_flag[1][(index)/(1+1)]=0;
      
      nfo = f->first->vertex(f->second)->info();
#if PERIODICITY!=1
      if ((minid>nfo.true_index)&&(nfo.true_index>=0)) {minid=nfo.true_index;minvertex=0;}
      if (periodic) net->f_vertexIndex[1][index++] = nfo.true_index;
      else {
	//if (nfo.boundary_flags) net->f_flag[1][(index)/(1+1)] |= NDNETFLAG_OUT;
	if (useIndexId) net->f_vertexIndex[1][index++] = nfo.index;
	else 
	  {
	    INT curindex=nfo.index;
		
	    if (curindex>=useIndexMapFrom)
	      {
		myMapT::iterator it=index2id.find(curindex);
		if (it==index2id.end()) {continue;}
		curindex=it->second;
	      }
	    net->f_vertexIndex[1][index++] = curindex;
	  }
      }
#else	    
      net->f_vertexIndex[1][index++] = nfo.index;
#endif

      nfo = f->first->vertex(f->third)->info();
#if PERIODICITY!=1
      if ((minid>nfo.true_index)&&(nfo.true_index>=0)) {minid=nfo.true_index;minvertex=1;}
      if (periodic) net->f_vertexIndex[1][index++] = nfo.true_index;
      else {
	//if (nfo.boundary_flags) net->f_flag[1][(index)/(1+1)] |= NDNETFLAG_OUT;
	if (useIndexId) net->f_vertexIndex[1][index++] = nfo.index;
	else 
	  {	    
	    INT curindex=nfo.index;
		
	    if (curindex>=useIndexMapFrom)
	      {
		myMapT::iterator it=index2id.find(curindex);
		if (it==index2id.end()) {index--;continue;}
		curindex=it->second;
	      }
	    net->f_vertexIndex[1][index++] = curindex;
	  }
      }
#else	    
      net->f_vertexIndex[1][index++] = nfo.index;
#endif
      
      if (minvertex==0)
	nfo = f->first->vertex(f->second)->info();
      else 
	nfo = f->first->vertex(f->third)->info();

#if PERIODICITY!=1
      if (!periodic) n_trueSeg++;
      else if (nfo.index == nfo.true_index) n_trueSeg++;
      else index-=(1+1);
#else
      n_trueSeg++;
#endif
    }

    if (n_trueSeg != net->nfaces[1])
      {
	net->nfaces[1]=n_trueSeg;
	net->f_flag[1]=(unsigned char*) realloc(net->f_flag[1],sizeof(unsigned char)*(size_t)net->nfaces[1]);
	net->f_vertexIndex[1]=(NDNET_UINT *)realloc(net->f_vertexIndex[1],sizeof(NDNET_UINT)*(size_t)net->nfaces[1]*(size_t)(1+1));
      }
    
#endif // NDIMS==3
    } // if allfaces
    // DOWN TO HERE 
 

    double *val1=(double *)malloc(sizeof(double)*(size_t)net->nvertex);  
    double *val2=(double *)malloc(sizeof(double)*(size_t)net->nvertex);  
    double *val3=(double *)malloc(sizeof(double)*(size_t)net->nvertex);  
    double *val4=NULL;
    double *val5=NULL;
    if (!useIndexId) val4=(double *)malloc(sizeof(double)*(size_t)net->nvertex);
    if (!useIndexId) val5=(double *)malloc(sizeof(double)*(size_t)net->nvertex);

    for (Finite_vertices_iterator vertex=this->finite_vertices_begin();vertex!=this->finite_vertices_end();vertex++)
    {      
      nfo=vertex->info();
#if PERIODICITY!=1
      if ((periodic)&&(nfo.index != nfo.true_index)) continue;
#endif
      
      if (useIndexId) 
	{
	  //printf("HIHO\n");
	  val1[nfo.index] = nfo.value;
	  val2[nfo.index] = log(nfo.value);
	  val3[nfo.index] = nfo.mass;
	}
      else 
	{
	  //printf("HIHO\n");
	  INT curindex=nfo.index;

	  if (curindex>=useIndexMapFrom)
	    {
	      myMapT::iterator it=index2id.find(curindex);
	      if (it==index2id.end()) continue;
	      curindex=it->second;
	    }

	  val1[curindex] = nfo.value;
	  val2[curindex] = log(nfo.value);
	  val3[curindex] = nfo.mass;
	  val4[curindex] = nfo.index;
	  val5[curindex] = nfo.true_index;
	}
    }
    
    addNDDataArr(net,0,VALUE_TAG,&val1);
    addNDDataArr(net,0,LOG_TAG(VALUE_TAG),&val2);
    addNDDataArr(net,0,MASS_TAG,&val3);
    if (val4!=NULL) addNDDataArr(net,0,INDEX_TAG,&val4);
    if (val5!=NULL) addNDDataArr(net,0,TRUE_INDEX_TAG,&val5);
 
#if PERIODICITY!=0
    net->periodicity=~0;
#endif
    if (periodic) net->periodicity=~0;
    
    #if NDIMS==2
    if (allFaces)
      {
	ComputeVertexFromFace(net,1);    
	free(net->v_faceIndex[2]);
	free(net->v_numFaceIndexCum[2]);    
	net->haveFaceFromVertex[2]=0;
	
	NDNetFlags_enable(net,1);
      }
    #endif
    
    if (verbose>0) printf ("done.\n");
    return net;
}

typedef MyDelaunay<> Delaunay ;

#endif
