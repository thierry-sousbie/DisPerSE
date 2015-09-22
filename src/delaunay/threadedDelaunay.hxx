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
#ifndef __THREADED_DELAUNEY_HXX__
#define __THREADED_DELAUNEY_HXX__

#include <unistd.h>
#include <sys/wait.h>
#include <limits>

#include "delaunay.hxx"
#include "boxPartition.hxx"
#include "NDnet_interface.hxx"

#include "global.h"

#ifdef USE_THREADS
//#include "pthreadWrapper.h"
//#include <errno.h>
//#include <pthread.h>
#endif

class ThreadedDelaunay
{
public:

  struct returnType
  {
    std::vector<long> nfaces;
  };

  struct params {
    bConditions *boundary;
    float *pos;
    float *mass;
    long npart;
    int ndims;
    double *x0,*delta;
    double borderSize; 
    int periodic;  
    char *outname;
    double margin;
    bool checkBoundaries;
    bool computeAllFaces;
  };

  returnType compute(params *p, int nchunks_target, int nthreads)
  {
    long i,j;
    int nDone=0;
    int curId=0;   
    returnType ret;
    
    std::vector<double> x0(p->x0,p->x0+p->ndims);
    std::vector<double> delta(p->delta,p->delta+p->ndims);
    
    int nchunks = partition.build(nchunks_target,x0,delta);

    std::vector<threadParams> tp(std::min(nthreads,nchunks));
    std::vector<pthread_t> thread(std::min(nthreads,nchunks));
    for (curId=0;curId<tp.size();curId++)
      {
	tp[curId]=newThreadParams(p,curId,nchunks);
	printf("Chunk %d/%d started: ",1+tp[curId].chunk,nchunks);
	for (j=0;j<p->ndims;j++)
	  printf("[%g,%g]",tp[curId].x0[j],tp[curId].x0[j]+tp[curId].delta[j]);
	printf("\n");
	//printf("starting thread for chunk %d ... ",tp[curId].chunk);fflush(0);
	int thread_ret=pthread_create(&thread[curId],NULL,
				      run,
				      (void*)(&(tp[curId])));
	//printf("done.\n");
      }

    ret.nfaces.resize(p->ndims+1);
    
    do {
      
      for (i=0;i<tp.size();i++)
	{
	  if (tp[i].done==1)
	    {
	      pthread_join(thread[i],NULL);
	   	      
	      for (j=0;j<ret.nfaces.size();j++) ret.nfaces[j]+=tp[i].ret.nfaces[j];
	      //printf("done.\n");
	      if (curId<nchunks)
		{	  
		  tp[i]=newThreadParams(p,curId++,nchunks);
		  
		  printf("Chunk %d/%d started: ",1+tp[i].chunk,nchunks);
		  for (j=0;j<p->ndims;j++)
		    printf("[%g,%g]",tp[i].x0[j],tp[i].x0[j]+tp[i].delta[j]);
		  printf("\n");
		  int thread_ret=pthread_create(&thread[i],NULL,
						run,
						(void*)(&(tp[i])));		  
		  
		}
	      else  tp[i].done=-1;
	      nDone++;
	    }
	}

      usleep(0.02 * 1.E6);
      
    } while (nDone!=nchunks);
    
    return ret;
  }
  
private:

  boxPartition partition;
  int npartitions;
  int nDone;

  struct threadParams
  {
    params *p;
    std::vector<double> x0;
    std::vector<double> delta;
    int done;    
    pthread_t thread;
    int chunk;
    int nchunks;
    returnType ret;
  };

  threadParams newThreadParams(params *p, int index,int nchunks)
  {
    threadParams tp;
    unsigned long i;

    tp.p=p;
    tp.done=0;
    partition.getX0(index,std::back_inserter(tp.x0));
    partition.getDelta(index,std::back_inserter(tp.delta));
    /*
    for (i=0;i<tp.x0.size();i++)
      {
	tp.x0[i]-=p->borderSize;
	tp.delta[i]+=2*p->borderSize;
      }
    */
    tp.chunk=index;
    tp.nchunks=nchunks;
    return tp;
  }

  static void *run(void *prms)
  {
    threadParams *tp=static_cast<threadParams *>(prms);  
    params *params = tp->p;
    returnType &ret=tp->ret;

    std::vector<Delaunay::Point> pts;
    std::vector<long> index;
    std::vector<INT> true_index;
    long nnew=0;
    long i,j,k,l;
    int ndims=params->ndims;
    
    std::vector<double> x0v;
    std::vector<double> deltav;
    std::vector<double> x0v_old;
    std::vector<double> deltav_old;
    
    long old_nnew=0;
    bool needRecomp;
    int nPasses=0;

    double borderSize = params->borderSize;
    double margin = params->margin;

    Delaunay *T = new Delaunay(tp->x0,tp->delta,0*params->periodic);

    std::vector<double> x0f=params->boundary->getX0();
    std::vector<double> deltaf=params->boundary->getDelta();
    
    do {
      needRecomp=false;
      nPasses++;

      x0v_old=x0v;
      deltav_old=deltav;
      x0v=tp->x0;
      deltav=tp->delta;
    
      for (i=0;i<ndims;i++)
	{
	  x0v[i]-=borderSize;
	  deltav[i]+=2*borderSize;
	}

      for (i=0,k=0;i<params->npart;i++) {
	std::vector<double> P(NDIMS);
	std::vector<double> result;      
	for (j=0;j<NDIMS;j++) P[j]=params->pos[i*NDIMS+j];
	bool direct = params->boundary->toSubBox(P, result, x0v,deltav);
      
	for (j=0;j<result.size()/ndims;j++)
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
		long dec=0;
		if (params->periodic) dec=i;
		else
		  {
		    for (l=0;l<NDIMS;l++)
		      {
			if ((result[j*NDIMS+l]<x0f[l])||
			    (result[j*NDIMS+l]>=x0f[l]+deltaf[l]))
			  dec|=(1<<l);
		      }
		    
		    dec = dec*params->npart+i;
		  }
		index.push_back(dec);
		nnew++;
	      }
	  }
	k++;
      }

       old_nnew=nnew;
      T->insert(pts.begin(),pts.end(),params->mass,NULL,&(index[0]),&(true_index[0]));
      
      if (params->checkBoundaries)
	{
	  std::pair<int,double> bcheck=T->CheckBoundaries(borderSize);
	
	  if (bcheck.first)
	    {
	      //printf("-%d-:WARNING: circumsphere test failed for %d simplexes.\n",tp->chunk,bcheck.first);
	      //printf("-%d-:  This means that the margin size was not correctly estimated.\n",tp->chunk);
	      //printf("-%d-:  Current margin is %e  = %.2f%% of box size.\n",tp->chunk,
	      //borderSize,  margin*100);
	      double suggestedBorder=(bcheck.second*1.001);
	      if (suggestedBorder>2*borderSize) suggestedBorder=2*borderSize;
	      double suggestedMargin=suggestedBorder/borderSize*margin;	
	      //printf("-%d-:  Will try ~%e.  (option '-margin %.5f')\n",tp->chunk,
	      //suggestedBorder,suggestedMargin);
	      needRecomp=true;

	      borderSize = suggestedBorder;
	      margin = suggestedMargin;
	    } 
	  //else printf("-%d-:Circumsphere test was successful.\n",tp->chunk);
	 
	}

      pts.clear();
      index.clear();
      true_index.clear();
    } while (needRecomp);

    T->SetValueToDensity();
    NDnetwork *net = T->ToNDnetwork(0*params->periodic,params->computeAllFaces);
    //net->periodicity=params->periodic;
    delete T;
    char fname[255];
    sprintf(fname,"%s%s_%5.5d",params->outname,ndnet::IO::getExtension().c_str(),tp->chunk);
    analyse(net,prms);
    if (params->periodic) for (i=0;i<ndims;i++) net->periodicity|=(1<<i);
    sprintf(net->comment,"File %d in %d",tp->chunk,tp->nchunks);
    ndnet::IO::save(net,std::string(fname));
    
    FreeNDnetwork(&net);

    tp->done=1;  
    return NULL;
  }


  static void analyse(NDnetwork* net, void *prms)
  {
    threadParams *tp=static_cast<threadParams *>(prms);  
    params *params = tp->p;
    returnType &ret=tp->ret;
    
    long i,j,k,l;    
  
    std::vector<double> x0(net->x0,net->x0+net->ndims);
    std::vector<double> delta(net->delta,net->delta+net->ndims);
    std::vector<double> x0f=params->boundary->getX0();
    std::vector<double> deltaf=params->boundary->getDelta();
    printf("X0 = %f %f\n",x0[0],x0f[0]);
    ret.nfaces.resize(net->ndims+1);
  
    for (i=0;i<net->nvertex;i++)
      {
	float *p=&net->v_coord[i*net->ndims];

	if (!(net->v_flag[i]&NDNETFLAG_OUT))
	  {
	    net->v_flag[i]|=NDNETFLAG_KEEP;
	    ret.nfaces[0]++;
	  }
	else if (!params->periodic)
	  {
	    if (!params->boundary->inBBox(p))
	      {
		net->v_flag[i]|=NDNETFLAG_KEEP;
		ret.nfaces[0]++;
	      }
	  }
	
      }
    
    for (j=1;j<=net->ndims;j++)
      for (i=0;i<net->nfaces[j];i++)
	{
	  NDNET_UINT *vertex = VERTEX_IN_FACE(net,j,i);
	  int flags=0;
	  for (k=0;k<j+1;k++)
	    {
	      if (net->v_flag[vertex[k]]&NDNETFLAG_KEEP) flags|=(1<<k);
	    }

	  if ( flags==((1<<(j+1))-1) )
	    {
	      net->f_flag[j][i]|=NDNETFLAG_KEEP;
	      ret.nfaces[j]++;
	    }
	  
	  if (!flags) continue;

	  int ref=0;
	  float *p_ref=&net->v_coord[vertex[0]*net->ndims];

	  for (k=1;k<j+1;k++)
	    {
	      float *p=&net->v_coord[vertex[k]*net->ndims];
	  
	      for (l=0;l<net->ndims;l++)
		{
		  if (p_ref[l]>p[l])
		    {
		      p_ref=p;
		      ref=k;
		      break;
		    }
		  else if (p_ref[l]!=p[l]) break;
		}
	    }
	  if (k!=j+1) continue;

	  bool inside_b = params->boundary->inBBox(p_ref,x0,delta);

	  if (inside_b) 
	    {
	      net->f_flag[j][i]|=NDNETFLAG_KEEP;
	      ret.nfaces[j]++;
	      for (k=0;k<j+1;k++) net->v_flag[vertex[k]]|=NDNETFLAG_SHARED;
	    }
	}
  }

};


#endif
