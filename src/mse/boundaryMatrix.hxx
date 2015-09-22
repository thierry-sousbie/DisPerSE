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
#ifndef BOUNDARY_MATRIX_HXX_
#define BOUNDARY_MATRIX_HXX_

#include "sparseZ2matrix.hxx"
#include "global.h"

#ifdef USE_OPENMP
#include <parallel/algorithm>
#endif

#ifdef USE_THREADS
#include "pthreadBarrier.hxx"
//#include "pthreadWrapper.h"
//#include <errno.h>
//#include <pthread.h>
#endif

//#define NOKILLING

class boundaryMatrix {
  
public:
  typedef sparseZ2Matrix matrixT;  
  typedef matrixT::idType idType;

private:  
  matrixT *bM;
  idType M_nCol;
  std::vector<idType> idx;

  struct solveBM_args{
    boundaryMatrix *This;    
    int threadID;
    int curType;
    int *NThreads;
    long imin;
    long imax;
    
    sparseZ2Matrix *M; 
    std::vector<int> *eq_col;
    bool skipLastCol;
    std::vector<int> *curindex;
    std::vector<bool> *waiting;   
  };

#ifdef HAVE_PTHREADS
#ifdef USE_THREADS
  pthread_mutex_t* curindex_mutex;
  pthread_cond_t* curindex_threshold_cv;  
  pthread_mutex_t* barrier_mutex;
  pthread_barrier_t* barrier;
#endif
#endif

public:

  matrixT *get()
  {
    return bM;
  }

  void init(matrixT *M)
  {
    bM=M;
    M_nCol=M->nCols();
  }  
  
  boundaryMatrix(matrixT *M)
  {
    init(M);
  }

  boundaryMatrix()
  {
    bM=NULL;
  }

  ~boundaryMatrix()
  {
    
  }

  void canonize()
  {
    long i,j;
    int ncol=bM->nCols();
    std::vector<int> low;
    std::vector<int> colOf;
    sparseZ2Matrix &M=*bM; 

    typedef matrixT::colType colType;
    typedef matrixT::colItType colType_it;

    low.assign(ncol,-1);
    colOf.assign(ncol,-1);
    
    for (i=0;i<ncol;i++)
      {
	long li=M.low(i);
	low[i]=li;
	if (li>=0) colOf[li]=i;
      }
    
    for (i=0;i<ncol;i++)
      {
	if (low[i]<0) continue;
	colType col=M.getCol(i);
	colType_it it=col->begin();
	while(it!=col->end())
	  {
	    long j=**it;
	    if ((j==low[i])||(colOf[j]<0))
	      {
		it++;
		continue;
	      }
	    
	    colType_it it2=M.addCol(i,colOf[j],it);

	    if (colOf[j]<i) it=it2;
	    else it=col->begin();
	  }
      }

  }

  void solve(int curType=0,bool skipLastCol=true,int nthreads=-1)
  {
    long i,j;
    sparseZ2Matrix &M=*bM; 
    
    if (bM==NULL) 
      {
	fprintf(stderr,"ERROR in boundaryMatrix: matrix was not initialized.\n");
	exit(0);
      }
    int ncol=bM->nCols();
    
    std::vector<int> eq_col(ncol,-1);

#ifdef HAVE_PTHREADS   
#ifdef USE_THREADS 
    
    pthread_t thread[glob_num_threads];
    int thread_ret[glob_num_threads];
    int nCurTh;
    solveBM_args pv[glob_num_threads]; 
    idx.clear();

    // Setup threads parameters
    //printf("nthreads = %d\n",nthreads);
    if (nthreads>1)
      {
	int nth =(nthreads<1)?glob_num_threads:nthreads;
	idx.clear();
	for (i=0;i<ncol;i++)
	  {
	    if (M.getCol(i)==M.nullCol()) continue;
	    idx.push_back(i);
	  }
	
	int nval = idx.size();
	if (!nval) return;
	
	nCurTh=nth;
	std::vector<int> imin(nth);
	std::vector<int> imax(nth);
	std::vector<int> curindex(nth+1);    
	std::vector<bool> waiting(nth,false);
	
	curindex_mutex=(pthread_mutex_t*)malloc((nth+1)*sizeof(pthread_mutex_t));
	curindex_threshold_cv=(pthread_cond_t*)malloc((nth+1)*sizeof(pthread_cond_t)); 
	barrier_mutex=(pthread_mutex_t*)malloc(1*sizeof(pthread_mutex_t));
	barrier=(pthread_barrier_t *)malloc(1*sizeof(pthread_barrier_t));
	
	for (i=0;i<=nth;i++) {
	  pthread_mutex_init(&curindex_mutex[i], NULL);
	  pthread_cond_init(&curindex_threshold_cv[i], NULL);
	}
	pthread_mutex_init(barrier_mutex,NULL);
			
	for (i=0;i<nth;i++) {
	  imin[i]=((long)i*(long)nval)/(nth+1);
	  imax[i]=((long)(i+1)*(long)nval)/(nth+1);
	  curindex[i]=imin[i];
	} 
	curindex[nth]=ncol+10;

	for (i=0;i<nth;i++) {	    
	  solveBM_args p ={this,(int)i,curType,&nCurTh,imin[i],imax[i],
			   bM,&eq_col,skipLastCol,&curindex,&waiting};//&(eq_col[i]),
	  pv[i]=p;
	}
	
	for (i=nth-1;i>=0;i--) 
	  thread_ret[i] = pthread_create(&thread[i],NULL,
					 solveBoundaryMatrix,
					 (void*)(&(pv[i])));

	for (i=nth-1;i>=0;i--) {
	  //printf("Waiting for threads %ld ... ",i);fflush(0);
	  pthread_join(thread[i],NULL);
	  //$printf("joined!\n");
	}
	
	for (i=0;i<idx.size();i++)
	  {
	    int li = M.low(idx[i]);
	    if (li>=0) M.clearCol(li);
	  }

	for (i=0;i<=nth;i++) {
	  pthread_mutex_destroy(&curindex_mutex[i]);
	  pthread_cond_destroy(&curindex_threshold_cv[i]);
	}
	free(curindex_mutex);
	free(curindex_threshold_cv);
	
	pthread_mutex_destroy(barrier_mutex);
	free(barrier_mutex);
	free(barrier);
	return;
      }
#endif
#endif
    //#else
    nthreads=1;
    //#endif

    if (nthreads==1)
      {
	//printf("(solve) ");fflush(0);           
	long ntot = ncol;
	if (skipLastCol) ntot--;
	
	for (i=0;i<ntot;i++)
	  {
	    long li=M.low(i);
	    if (li<0) continue;
	    
	    while (eq_col[li]>=0) {
	      M.addCol(i,eq_col[li]);
	      li=M.low(i);
	      if (li<0) break;
	    };
	    
	    if (li>=0) eq_col[li] = i;
	  } 
	
	if (skipLastCol) M.clearCol(ncol-1);
	//M.addCol(ncol-1,ncol-1);
      }
  }

  
private:
  
  static void *solveBoundaryMatrix(void *prms)
  {
    long i=0,j,k;
    solveBM_args *p = static_cast<solveBM_args *>(prms);    
    sparseZ2Matrix &M=*(p->M); 
    long ncol = M.nCols();
    //std::vector<int> &eq_col=*(p->eq_col);
    std::vector<int> &curindex=*(p->curindex);
    std::vector<bool> &waiting=*(p->waiting); 
    int &threadID = p->threadID;
    boundaryMatrix *This=p->This;  
    int delay=-1;
    int nth;

    int curType=p->curType;
    std::vector<idType> &idx=This->idx;
#ifdef HAVE_PTHREADS
#ifdef USE_THREADS
    pthread_mutex_t* curindex_mutex=This->curindex_mutex;
    pthread_cond_t* curindex_threshold_cv=This->curindex_threshold_cv;
    pthread_mutex_t* barrier_mutex=This->barrier_mutex;
    pthread_barrier_t* barrier=This->barrier;
#endif
#endif
    static std::vector<int> deltaP_a;
    static std::vector<int> deltaN_a;
    
    sparseZ2Matrix::colListType Mlist;
    std::vector<int> eq_col;
    std::vector<sparseZ2Matrix::colType> eq_colM;
    
    char fname[255];
    FILE *f=NULL;
    struct timeval t_new;
    double dt,dti;

#pragma omp critical
    if (debug_dump) 
      {
	static int fcount=0;
	sprintf(fname,"timing_%2.2d_%2.2d_%2.2d.dat",curType,threadID,fcount++);

	f=fopen(fname,"w");
	gettimeofday(&t_new, NULL);
	dt=(double)t_new.tv_sec + ((double)t_new.tv_usec)*1.E-6;
	fprintf(f,"%ld %16.16g 0\n",(long)idx[p->imin],dt);
      }


    //printf("Started thread %d\n",threadID);
    
    //struct timeval wc_time1,wc_time2;
    //gettimeofday(&wc_time1, NULL);
    if (threadID==0)
      {
	deltaP_a.assign(*p->NThreads+1,1);
	deltaN_a.assign(*p->NThreads+1,-1);
	eq_col.assign(ncol,-1);
	for (k=p->imin;k<p->imax;k++)
	  {
	    i=idx[k];
	    if (p->skipLastCol)
	      if (i==ncol-1) {
		M.clearCol(i);
		continue;
	      }
	    
	    long li=M.low(i);
	    if (li<0) continue;
	    
	    while (eq_col[li]>=0) {
	      M.addCol(i,eq_col[li]);
	      li=M.low(i);
	      if (li<0) break;
	    }
	    if (li>=0) eq_col[li] = i;
	    //continue;	  
	  }
      }
    else
      {
	eq_colM.assign(ncol,M.nullCol());
	for (k=p->imin;k<p->imax;k++)
	  {
	    i=idx[k];
	    if (p->skipLastCol)
	      if (i==ncol-1) continue;
	    
	    long li=M.low(i);	  
	    if (li<0) continue;
	    
	    while (eq_colM[li]!=M.nullCol()) {
	      M.addCol(i,eq_colM[li]);
	      li=M.low(i);
	      if (li<0) break;
	    }
	    if (li>=0) {
	      //eq_colM[li] = M.getCol(i);
	      Mlist.push_front(*M.getCol(i));
	      eq_colM[li]=Mlist.begin();
	    }
	    //continue;	  
	  }
      }
#ifdef HAVE_PTHREADS    
#ifdef USE_THREADS        
    pthread_mutex_lock(barrier_mutex);
    if (!waiting[threadID])
      {
	for (j=0;j<(*p->NThreads);j++) if (j!=threadID) waiting[j]=true;
	pthread_barrier_init(barrier,NULL,(*p->NThreads));
	//printf("barrier set by thread %d\n",threadID);
      }
    pthread_mutex_unlock(barrier_mutex);
    curindex[threadID]=p->imax; 
    //printf("Barrier reached by thread %d\n",threadID);
    pthread_barrier_wait(barrier);

    int &deltaP=deltaP_a[threadID];
    int &deltaN=deltaN_a[threadID];

    if (!waiting[threadID]) {
      pthread_barrier_destroy(barrier);
      //printf("Barrier destroyed\n");
    }
    waiting[threadID]=false;

#endif 
#endif

    gettimeofday(&t_new, NULL);
    dt=(double)t_new.tv_sec + ((double)t_new.tv_usec)*1.E-6;
    dti=dt;
    if (debug_dump) fprintf(f,"%ld %16.16g 0\n",i,dt);
   
    const int len = idx.size();

    for (k=p->imax;k<len;k++) 
      {
	i=idx[k];

	if (debug_dump) {
	  if (threadID!=0)
	    if (((k-p->imax) % 100) == 0) {
	      gettimeofday(&t_new, NULL);
	      dt=(double)t_new.tv_sec + ((double)t_new.tv_usec)*1.E-6;
	      fprintf(f,"%ld %16.16g 0\n",i,dt);
	    }
	}
#ifdef HAVE_PTHREADS
#ifdef USE_THREADS      
	curindex[threadID]=i; 
	
	if (curindex[threadID+deltaP]<=curindex[threadID]+2)
	  {
	    pthread_mutex_lock(&curindex_mutex[threadID]);	  
	    if (curindex[threadID+deltaP]<ncol) 
	      {
		waiting[threadID]=true;
		//printf("--Thread %d is waiting at %ld/%d/%ld (%d t)\n",threadID,i,curindex[threadID+deltaP],ncol,(*p->NThreads));//
		pthread_cond_wait(&curindex_threshold_cv[threadID], &curindex_mutex[threadID]);	  
		waiting[threadID]=false;
#ifndef NOKILLING
		if (threadID==0) {
			
		  deltaP+=deltaP_a[threadID+deltaP];
		  deltaN_a[threadID+deltaP]=-deltaP;
		}
	
#endif
#endif
		
	      }
	    pthread_mutex_unlock(&curindex_mutex[threadID]);
	  }

	
	if ((threadID>0)&&(waiting[threadID+deltaN]))
	  {
	   
	    if (delay==-1) delay=1000;
#ifndef NOKILLING
	    if (threadID+deltaN==0)
	      {
	
		pthread_mutex_lock(&curindex_mutex[threadID+deltaN]);
		pthread_cond_signal(&curindex_threshold_cv[threadID+deltaN]);      
		i=ncol;k=len;
		curindex[threadID]=i;
		pthread_mutex_unlock(&curindex_mutex[threadID+deltaN]);
	
		continue;
	      }
	    else 
#endif
	      if (curindex[threadID+deltaN]<=curindex[threadID]-delay)
		{
		  pthread_mutex_lock(&curindex_mutex[threadID+deltaN]);
		  pthread_cond_signal(&curindex_threshold_cv[threadID+deltaN]);
		  pthread_mutex_unlock(&curindex_mutex[threadID+deltaN]);
		}
	  }
	
#endif
	if (threadID==0)
	  {
	   
	    if (p->skipLastCol)
	      if (i==ncol-1) {
		M.clearCol(i);
		continue;
	      }
	    long li=M.low(i);
	    if (li<0) continue;	    
	    
	    while (eq_col[li]>=0) {
	      M.addCol(i,eq_col[li]);
	      li=M.low(i);	    
	      if (li<0) break;
	    }
	    if (li>=0) eq_col[li] = i;

	    if (debug_dump) {
	      {
		gettimeofday(&t_new, NULL);
		dt=(double)t_new.tv_sec + ((double)t_new.tv_usec)*1.E-6;
		fprintf(f,"%ld %16.16g %ld\n",i,dt,li);
	      }
	    }
	  }
	else
	  {
	    if (p->skipLastCol)
	      if (i==ncol-1) continue;
	    long li=M.low(i);	  
	    if (li<0) continue;
	    //int count=0;
	    while (eq_colM[li]!=M.nullCol()) {

	      M.addCol(i,eq_colM[li]);
	  
	      li=M.low(i);
	   
	      if (li<0) break;
	    }
	
	    if (li>=0) {
	    
	      if (i!=ncol-1)
		{
		  Mlist.push_front(*M.getCol(i));
		  eq_colM[li]=Mlist.begin();
		}
	    }
	  }
      }

    if (debug_dump) {
      gettimeofday(&t_new, NULL);
      dt=(double)t_new.tv_sec + ((double)t_new.tv_usec)*1.E-6;
      fprintf(f,"%ld %16.16g 0\n",i,dt);
    }

#ifdef HAVE_PTHREADS
#ifdef USE_THREADS  
    curindex[threadID]=ncol+10; 
    (*p->NThreads)--;
    if (threadID>0)
      {
	pthread_mutex_lock(&curindex_mutex[threadID+deltaN]);
	pthread_cond_signal(&curindex_threshold_cv[threadID+deltaN]);
	pthread_mutex_unlock(&curindex_mutex[threadID+deltaN]);
      }
  
#endif
#endif

    
    gettimeofday(&t_new, NULL);
    dt=(double)t_new.tv_sec + ((double)t_new.tv_usec)*1.E-6;
    if (debug_dump) {
      fprintf(f,"%ld %16.16g 0\n",i,dt);
      fclose(f);
    }
#ifdef HAVE_PTHREADS
#ifdef USE_THREADS  
    pthread_exit(NULL);
#endif
#endif

    return NULL;
  }
  

};

#endif
