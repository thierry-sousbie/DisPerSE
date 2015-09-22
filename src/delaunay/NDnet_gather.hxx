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
#ifndef __NDNET_GATHER_HXX__
#define __NDNET_GATHER_HXX__

#include <stdio.h>
#include <math.h>
#include <vector>

#include "NDnetwork.h"
#include "NDnetwork_tags.h"
#include "global.h"
#include "mytypes.h"

#include "find_unordered_maps.hxx"

#define GATHER_NDNET_READ_WRITE_COPY_SELECT(TYPE,BSIZE,FLAG,FIN,INFO_IN,FOUT,INFO_OUT,WHAT) \
  {									\
    TYPE *__macr_data =new TYPE [ (INFO_IN)->rn_ ## WHAT ];		\
    NDNET_READ_DATA(__macr_data,FIN,INFO_IN,WHAT);			\
    long __macr_ii=0;							\
    long __macr_count=0;						\
    for (__macr_ii=0 ;__macr_ii < FLAG.size() ;__macr_ii++)		\
      {									\
	if (FLAG[__macr_ii])						\
	  {								\
	    if (__macr_count==__macr_ii)				\
	      {								\
		__macr_count++;						\
		continue;						\
	      }								\
	    memcpy(&__macr_data[__macr_count*BSIZE],&__macr_data[__macr_ii*BSIZE],sizeof(TYPE)*BSIZE); \
	    __macr_count++;						\
	  }								\
      }									\
    NDNET_WRITE_DATA(__macr_data,FOUT,INFO_OUT,WHAT,__macr_count*BSIZE); \
    delete __macr_data;							\
  }									\




class NDnetGather
{
  typedef my_unordered_map<long,INT>::type myMapT;

  std::string in_fname;    
  NDnetwork *net_dummy;
  std::vector<NDnetwork*> net_in;
  std::vector<NDnet_fstruct_info> net_info_in;
  std::vector<long> start_index;
  myMapT new_index;
  int nfiles;

public:

  void init(std::string fname)
  {
    long i,j;    
    int ftot=-1;

    in_fname=fname;

    nfiles=0;
    char name[255];
    printf("Gathering '%s_#####' network files ... (%d found)",fname.c_str(),0);fflush(0);
    start_index.assign(1,0);    

    while (true)
      {	
	sprintf(name,"%s_%5.5d",fname.c_str(),nfiles);
	if (access(name,F_OK) == -1) break;
	printf("\rGathering '%s_#####' network files ... (%d found)",fname.c_str(),nfiles+1);fflush(0);
	//printf("%d: %s\n",nfiles,name);
	NDnet_fstruct_info info;
	NDnetwork *net=Load_NDnetwork_header(name, &info);	

	if (!nfiles)
	  {	    
	    int fid=-1;
	    int nsc = sscanf(net->comment,"File %d in %d",&fid,&ftot);
	    char test[255];
	    sprintf(test,"File %d in %d",fid,ftot);
	    if (!strcmp(net->comment,test))
	      {
		
	      }
	      

	    if (net_dummy!=NULL) 
	      FreeNDnetwork(&net_dummy);
	    net_dummy=CreateNetwork(net->ndims,0,0); 
	    memcpy(net_dummy->x0,net->x0,net->ndims*sizeof(double));
	    memcpy(net_dummy->delta,net->delta,net->ndims*sizeof(double));
	    net_dummy->periodicity=net->periodicity;
	    
	  }
	else
	  {
	    int fid=-1;
	    int locftot=-1;
	    int nsc = sscanf(net->comment,"File %d in %d",&fid,&locftot);
	    if ((nsc!=2)||(locftot != ftot))
	      {
		FreeNDnetwork(&net);
		break;
	      }
	    for (i=0;i<net->ndims;i++)
	      {
		if (net->x0[i]<net_dummy->x0[i])
		  {
		    double mx=net_dummy->x0[i]+net_dummy->delta[i];
		    net_dummy->x0[i]=net->x0[i];
		    net_dummy->delta[i]=mx-net_dummy->x0[i];
		  }
		if (net->x0[i]+net->delta[i] > net_dummy->x0[i]+net_dummy->delta[i])
		  {
		    net_dummy->delta[i]=(net->x0[i]+net->delta[i])-net_dummy->x0[i];
		  }
	      }
	  }
	
	net_in.push_back(net);
	net_info_in.push_back(info);
	
	int tagID=NDDataIndex(net,0,INDEX_TAG);
	if (tagID<0)
	  {
	    fprintf(stderr,"ERROR in NDnetGather: tag %s was not found in file %s, cannot continue ...\n",INDEX_TAG,name);
	    exit(-1);
	  }
	
	if (!net->haveVFlags)
	  {
	    fprintf(stderr,"ERROR in NDnetGather: vertices do not have flags, invalid network, cannot continue ...\n");
	    exit(-1);
	  }

	std::vector<double> index(net->nvertex);
	std::vector<unsigned char> vflag(net->nvertex);
	
	FILE *f=fopen(name,"r");
	NDNET_READ_DATA(&index[0],f,&info,data[tagID]);
	NDNET_READ_DATA(&vflag[0],f,&info,v_flag);
	
	start_index.push_back(start_index.back());
	//printf("nvertex = %ld\n",(long)net->nvertex);
	if (!nfiles) net_dummy->haveVFlags=1;
	for (i=0;i<net->nvertex;i++) 
	  {
	    if (vflag[i]&NDNETFLAG_KEEP) 
	      {		
		if (vflag[i]&NDNETFLAG_SHARED)
		  new_index.insert(std::make_pair(index[i],start_index.back()));
		
		start_index.back()++;
		net_dummy->nvertex++;
	      }
	  }
	
	for (i=0;i<=net->ndims;i++) 
	  {
	    if (!net->haveVertexFromFace[i]) continue;
	    if (!net->haveFFlags[i]) continue;
	    if (!nfiles) net_dummy->haveVertexFromFace[i]=1;
	    if (!nfiles) net_dummy->haveFFlags[i]=1;
	    std::vector<unsigned char> fflag(net->nfaces[i]);
	    NDNET_READ_DATA(&fflag[0],f,&info,f_flag[i]);
	    for (j=0;j<net->nfaces[i];j++)
	      {
		if (fflag[j]&NDNETFLAG_KEEP) 
		  net_dummy->nfaces[i]++;
	      }
	  }

	if (!nfiles) 
	  {
	    net_dummy->ndata=net->ndata;
	    net_dummy->data=(NDnetwork_Data*)calloc(net->ndata,sizeof(NDnetwork_Data));
	  }
	  
	for (i=0;i<net->ndata;i++)
	  {
	    net_dummy->data[i].type=net->data[i].type;
	    strcpy(net_dummy->data[i].name, net->data[i].name);
	  }

	fclose(f);
	/*
	for (i=1;i<=net->ndims;i++)  printf("nfaces[%ld] : %ld\n",(long)i,(long)net->nfaces[i]);
	
	printf("Size : %ld, map : %ld\n",start_index.back(),new_index.size());
	*/
	//FreeNDnetwork(&net);
	nfiles++;
      }

    if (!nfiles)
      {
	fprintf(stderr,"\nERROR in NDnetGather: files %s_????? not found.\n",fname.c_str());
	exit(-1);
      }
    printf("\rGathering '%s_#####' network files ... done. (%d found)\n",fname.c_str(),nfiles);fflush(0);
    sprintf(net_dummy->comment,"Gathered from %d chuncks",nfiles);
	
    //printf("nvertex : %ld\n",(long)net_dummy->nvertex);
    //for (i=1;i<=net_dummy->ndims;i++)  printf("nfaces[%ld] : %ld\n",(long)i,(long)net_dummy->nfaces[i]);
  }

  NDnetGather(std::string fname=std::string())
  {
    net_dummy=NULL;
    init(fname);    
  }

  ~NDnetGather()
  {
    long i;

    if (net_dummy!=NULL) 
      FreeNDnetwork(&net_dummy);
    for (i=0;i<net_in.size();i++)
      FreeNDnetwork(&net_in[i]);
  }
  /*
  NDnetwork *toNDnet(std::string fname)
  {

  }
  */
  void write(std::string fname)
  {
    printf("Writing gathered network to file '%s':\n",fname.c_str());fflush(0);
    FILE *fout;
    NDnet_fstruct_w_info *info_out = 
      Save_NDnetwork_header(net_dummy,fname.c_str(),&fout);

    char name[255];
    int nf=0;
    long i,j,k;

    while (nf<nfiles)
      {	
	sprintf(name,"%s_%5.5d",in_fname.c_str(),nf);
	if (nf==0) {printf("  adding file '%s' (%d/%d) ... ",name,nf+1,nfiles);fflush(0);}
	else {printf("\r  adding file '%s' (%d/%d) ... ",name,nf+1,nfiles);fflush(0);}
	
	NDnetwork *net=net_in[nf];
	NDnet_fstruct_info *info_in=&net_info_in[nf];
	std::vector< std::vector<bool> > selection;
	std::vector<bool> shared;
	std::vector<NDNET_UINT> newId;
	selection.resize(net_dummy->ndims+1);
	
	FILE *f=fopen(name,"r");
	
	if (net_dummy->haveVFlags)
	  {
	    std::vector<unsigned char> vflag(net->nvertex);
	    selection[0].resize(net->nvertex);
	    shared.resize(net->nvertex);
	    NDNET_READ_DATA(&vflag[0],f,info_in,v_flag);
	    for (i=0;i<net->nvertex;i++) 
	      {
		selection[0][i]=(vflag[i]&NDNETFLAG_KEEP)?true:false;
		shared[i]=(vflag[i]&NDNETFLAG_SHARED)?true:false;
	      }
	    vflag.clear();
	    j=start_index[nf];
	    
	    int tagID=NDDataIndex(net,0,INDEX_TAG);
	    std::vector<double> index(net->nvertex);
	    NDNET_READ_DATA(&index[0],f,info_in,data[tagID]);
	    newId.resize(net->nvertex);
	    for (i=0;i<net->nvertex;i++) 
	      {
		if (selection[0][i]) newId[i]=j++;
		else if (shared[i]) 
		  {
		    //if (new_index.find(index[i])==new_index.end()) printf("ERROR ...\n");
		    newId[i]=new_index[index[i]];
		  }
		else newId[i]=-1;
	      }
	    
	    
	    GATHER_NDNET_READ_WRITE_COPY_SELECT(unsigned char,1,selection[0],f,info_in,fout,info_out,v_flag);
	  }
	for (i=1;i<=net_dummy->ndims;i++)
	  {
	    if (net_dummy->haveFFlags[i])
	      {
		std::vector<unsigned char> fflag(net->nfaces[i]);
		selection[i].resize(net->nfaces[i]);
		NDNET_READ_DATA(&fflag[0],f,info_in,f_flag[i]);
		for (j=0;j<net->nfaces[i];j++) selection[i][j]=(fflag[j]&NDNETFLAG_KEEP)?true:false;
		fflag.clear();
		
		GATHER_NDNET_READ_WRITE_COPY_SELECT(unsigned char,1,selection[i],f,info_in,fout,info_out,f_flag[i]);
	      }
	  }
	
	GATHER_NDNET_READ_WRITE_COPY_SELECT(float,net->ndims,selection[0],f,info_in,fout,info_out,v_coord);

	for (i=0;i<net_dummy->ndata;i++)
	  {
	    GATHER_NDNET_READ_WRITE_COPY_SELECT(double,1,selection[net_dummy->data[i].type],f,info_in,fout,info_out,data[i]);
	  }
		
	for (i=1;i<=net->ndims;i++)
	  {
	    if (net->nfaces[i])
	      {
		NDNET_UINT *data =new NDNET_UINT [net->nfaces[i]*(i+1)];	
		NDNET_READ_DATA(data,f,info_in,f_vertexIndex[i]);
		long count=0;						
		for (j=0 ;j < net->nfaces[i] ;j++)		
		  {									
		    if (selection[i][j])						
		      {								
			if (count!=j) memcpy(&data[count*(i+1)],&data[j*(i+1)],sizeof(NDNET_UINT)*(i+1)); 
			for (k=0;k<i+1;k++) 
			  {
			    //if (newId[data[count*(i+1) + k]] <0) printf("AYYYEEE\n");
			    data[count*(i+1) + k] = newId[data[count*(i+1) + k]];
			  }
			count++;						
		      }								
		  }									
		NDNET_WRITE_DATA(data,fout,info_out,f_vertexIndex[i],(i+1)*count); 
		delete data;							
	      }									
	  }

	nf++;
      }
    free(info_out);
    printf("all done.\n");
    printf ("\nNetwork was saved as: %s\n",fname.c_str());
    printNDnetStat(net_dummy,3);
  }
};

#endif
