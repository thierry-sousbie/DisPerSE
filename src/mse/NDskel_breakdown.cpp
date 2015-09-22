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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <list>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <limits>

#include "NDnetwork_tags.h"

#include "NDskel_breakdown.h"
#include "NDskel_breakdown_tools.hxx"
#include "global.h"

long filamentID(int a, int b)
{
  //return (long)a + (((long)b)<<31);
  
  if (b<a)
    return (long)a + (((long)b)<<31);
  else
    return (long)b + (((long)a)<<31);
  
}

NDskel *NDskelBreakdown(NDskel *skl)
{  
  long i,j;  
  std::vector<NDskl_seg *> segs;
  compareSegPosEqual cmpEq(skl->ndims);
  compareSignatureEq cmpSigEq;

  segs.reserve(skl->nsegs);
  for (NDskl_seg *seg=skl->Seg;seg!=skl->Seg+skl->nsegs;seg++)
    segs.push_back(seg);

  std::sort(segs.begin(),segs.end(),compareSegPosLess(skl->ndims));

  std::vector<int> newSegIndex(skl->nsegs,-1);
  std::vector< std::vector<long> > filSig(1); 
  std::vector< std::vector<long>* > fils;
  std::vector< int > filId;

  int robustness_id=getDataFieldID(skl,1,ROBUSTNESS_TAG);
  int robustness_ratio_id=getDataFieldID(skl,1,ROBUSTNESS_RATIO_TAG);

  int newNSegs=0;   
  newSegIndex[segs[0]->index]=newNSegs++;
  filSig[newNSegs-1].push_back(filamentID(segs[0]->nodes[0],segs[0]->nodes[1]));

  NDskl_seg *lastSeg=segs[0];

  for (i=1;i<(long)segs.size();i++)
    {     
      if (cmpEq(segs[i],segs[i-1]))
	{
	  bool replace=false;
	 
	  if (robustness_id>=0)
	    {
	      double nval=segs[i]->data[robustness_id];
	      double oval=lastSeg->data[robustness_id];
	      if (nval==oval)
		{
		  nval=segs[i]->data[robustness_ratio_id];
		  oval=lastSeg->data[robustness_ratio_id];
		}
	      replace=(oval<nval);
	    }
	  replace=false;
	  if (!replace)
	    {
	      newSegIndex[segs[i]->index]=-1;
	      filSig[newNSegs-1].push_back(filamentID(segs[i]->nodes[0],segs[i]->nodes[1]));
	    }
	  else
	    {
	      newSegIndex[lastSeg->index]=-1;
	      newSegIndex[segs[i]->index]=newNSegs-1;
	      filSig[newNSegs-1].push_back(filamentID(segs[i]->nodes[0],segs[i]->nodes[1]));
	      std::swap(filSig[newNSegs-1].front(),filSig[newNSegs-1].back());
	      lastSeg=segs[i];
	    }
	}
      else 
	{	  
	  newSegIndex[segs[i]->index]=newNSegs++;
	  filSig.push_back(std::vector<long>(1,filamentID(segs[i]->nodes[0],segs[i]->nodes[1])));
	  lastSeg=segs[i];
	}
    }
  segs.clear();
  //printf ("segments : %d --> %d !!!\n",skl->nsegs,newNSegs);
  

  fils.reserve(filSig.size());
  for (i=0;i<(long)filSig.size();i++)
    fils.push_back(&filSig[i]);

  std::sort(fils.begin(),fils.end(),compareSignatureLess());
  
  int newNFils=0;
  std::vector<long> *fil_ref=&filSig[0];

  filId.resize(filSig.size());
  filId[(long)(fils[0]-fil_ref)] = newNFils++;
  for (i=1;i<(long)fils.size();i++)
    {
      if (! cmpSigEq(fils[i],fils[i-1])) newNFils++;
      filId[(long)(fils[i]-fil_ref)]=newNFils-1;
    }
  filSig.clear();
  fils.clear();
  //printf ("--> %d filaments ?\n",newNFils);

  //std::set<newNodeT> newNodes;
  //typename std::set<newNodeT>::iterator nn;

  skeletonBuildT newSkl;
  
  for (i=0;i<skl->nnodes;i++)
    {
      NDskl_node *node=&(skl->Node[i]);
     
      skeletonBuildT::newNodeIt cn=newSkl.addNode(node,skl->ndims);
      
      for (j=0;j<node->nnext;j++) 
	{	  
	  NDskl_node *nextnode=node->Next[j];
	  if (nextnode->index>i) continue;
	  if (cn==newSkl.addNode(nextnode,skl->ndims)) continue;
	
	  //here we start a new filament

	  NDskl_seg *seg=node->Seg[j];
	  //float *pos;
	  int lastFilId;
	  int reverse;
	  //int ct=0;

	  if (seg->nodes[0]==i) reverse=0; else reverse=1;
	  //if (reverse) continue;
	 
	  if (newSegIndex[seg->index]>=0) 
	    {
	      //printf("_");
	      skeletonBuildT::newNodeIt nn=newSkl.addNode(node,skl->ndims);
	      newSkl.addFilament(); //filId[newSegIndex[seg->index]]
	      newSkl.setFilamentSrc(nn);
	      newSkl.AddSegToFilament(seg);
	      lastFilId=filId[newSegIndex[seg->index]];
	      //printf("id : %d->",lastFilId);
	    }
	  else
	    {
	      //printf("~");
	      lastFilId=-1;
	    }

	  NDskl_seg *next = (reverse)?seg->Prev:seg->Next;
	  while (next!=NULL)
	    {
	      
	      seg=next;
	      next = (reverse)?seg->Prev:seg->Next;
	      
	      int newFilId=-1;
	      if (newSegIndex[seg->index]>=0) newFilId=filId[newSegIndex[seg->index]];

	      if (newFilId != lastFilId)
		{
		  skeletonBuildT::newNodeIt nn=newSkl.addNode(seg,skl->ndims,reverse);
		  
		  if (lastFilId>=0)
		    newSkl.setFilamentDst(nn);
		  
		  if (newSegIndex[seg->index]>=0) 
		    {
		      newSkl.addFilament();
		      newSkl.setFilamentSrc(nn);
		    }
		  
		  lastFilId=newFilId;
		}

	      if (newSegIndex[seg->index]>=0) newSkl.AddSegToFilament(seg);
	     
	    }
	  
	  if (newSegIndex[seg->index]>=0) 
	    {
	      skeletonBuildT::newNodeIt nn=newSkl.addNode(nextnode,skl->ndims);
	      newSkl.setFilamentDst(nn);
	    }

	}
    }
  
  return newSkl.toNDskel(skl);
}


NDnetwork *NDskel2NDnet(NDskel *skl_)
{     
  NDskel *skl=skl_;
  long i,j,k;
  
  long cur=0;
  long curseg=0;

  long nnodeDataSup=1;
  long nsegDataSup=3;
  double *segData[skl->nsegdata+nsegDataSup];
  double *nodeData[skl->nnodedata+nnodeDataSup];  

  long nfil=0;
  NDskl_seg **filTab=NULL;
  int *filSize=NULL;
  nfil=getNDskelFilTab(skl,&filTab,&filSize);

  long nvert=skl->nnodes;
  for (i=0;i<nfil;i++) nvert+=filSize[i]-1;
  NDnetwork *net=CreateNetwork(skl->ndims,nvert,0);

  for (i=0;i<skl->nsegdata+nsegDataSup;i++)
    segData[i]=(double*)calloc(skl->nsegs,sizeof(double));
  for (i=0;i<skl->nnodedata+nnodeDataSup;i++)
    {
      nodeData[i]=(double*)calloc(nvert,sizeof(double));
      for (j=0;j<nvert;j++) nodeData[i][j]=-std::numeric_limits<double>::max();
    }
   
  for (i=0;i<skl->nnodes;i++)
    {
      memcpy(&net->v_coord[i*skl->ndims],skl->Node[i].pos,sizeof(float)*net->ndims);
      for (j=0;j<skl->nnodedata;j++)
	nodeData[j][i]=skl->Node[i].data[j];
      for (j=skl->nnodedata;j<skl->nnodedata+nnodeDataSup;j++)
	nodeData[j][i]=skl->Node[i].type;
    }

  net->nfaces[1]=skl->nsegs;
  net->haveVertexFromFace[1]=1;
  net->f_vertexIndex[1]=(NDNET_UINT*)malloc(sizeof(NDNET_UINT)*2*net->nfaces[1]);
  cur=skl->nnodes;
     
  memcpy(net->x0,skl->x0,sizeof(double)*skl->ndims);
  memcpy(net->delta,skl->delta,sizeof(double)*skl->ndims);
  memcpy(net->comment,skl->comment,sizeof(char)*80);  
  
  for (i=0;i<nfil;i++)
    {
      NDskl_seg *seg=filTab[i];
      NDskl_node *node=&skl->Node[seg->nodes[0]];
      NDskl_node *next=&skl->Node[seg->nodes[1]];
      
      for (k=0;k<skl->nsegdata;k++) 
	segData[k][curseg>>1]=seg->data[k];
      segData[k++][curseg>>1]=i;
      segData[k++][curseg>>1]=(node->type<next->type)?node->index:next->index;
      segData[k++][curseg>>1]=(node->type>next->type)?node->index:next->index;
      
      net->f_vertexIndex[1][curseg++]=node->index;
      
      if (seg->Next!=NULL)
	{
	  
	  do {
	    net->f_vertexIndex[1][curseg++]=cur;
	    seg=seg->Next;
	    memcpy(&net->v_coord[net->ndims*cur],seg->pos,sizeof(float)*net->ndims);	   
	    
	    for (k=0;k<skl->nsegdata;k++) 
	      segData[k][curseg>>1]=seg->data[k];
	    segData[k++][curseg>>1]=i;
	    segData[k++][curseg>>1]=(node->type<next->type)?node->index:next->index;
	    segData[k++][curseg>>1]=(node->type>next->type)?node->index:next->index;
	    
	    net->f_vertexIndex[1][curseg++]=cur++;
	  } while(seg->Next!=NULL);
	  	  
	}
      
      net->f_vertexIndex[1][curseg++]=next->index;    
      
    }
   freeNDskelFilTab(&filTab,&filSize);

  int oldverbose=verbose;
  verbose=0;

  for (k=0;k<skl->nnodedata;k++) 
    addNDDataArr(net,0,skl->nodedata_info[k],&nodeData[k]);
  addNDDataArr(net,0,CRITICAL_INDEX_TAG,&nodeData[k++]);

  for (k=0;k<skl->nsegdata;k++) 
    addNDDataArr(net,1,skl->segdata_info[k],&segData[k]);
  addNDDataArr(net,1,ARC_ID_TAG,&segData[k++]);
  addNDDataArr(net,1,DOWN_TAG(INDEX_TAG),&segData[k++]);
  addNDDataArr(net,1,UP_TAG(INDEX_TAG),&segData[k++]);

  net->haveVFlags=1;
  net->v_flag=(unsigned char *) calloc (net->nvertex,sizeof(unsigned char));
  for (i=0;i<skl->nnodes;i++) net->v_flag[i]=(unsigned char)skl->Node[i].flags;

  verbose=oldverbose;

  return net;
}
  
