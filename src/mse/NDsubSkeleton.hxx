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
#ifndef __NDSKEL_SUBSKELETON_HXX__
#define __NDSKEL_SUBSKELETON_HXX__

#include "NDskeleton.h"
#include "NDskel_tags.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <list>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include "mystring.h"
#include "global.h"

class skeletonBuildT;

class newNodeT {
  friend class skeletonBuildT;
private:
  
  NDskl_node *ndref;
  NDskl_seg *sgref;
  int segExtId;
  std::vector<float> pos;
  mutable int index; 
 
public:
  long getNDims() const {return pos.size();}
  float getCoord(int i) const {return pos[i];}
  int getType() const {return (ndref==NULL)?pos.size()+1:ndref->type;}
  int getFlags() const {return (ndref==NULL)?sgref->flags:ndref->flags;}
  double* getData() const {return (ndref==NULL)?sgref->data:ndref->data;}

  template <typename T>
  newNodeT(NDskl_node *nd, T ndims)
  {
    pos.assign(nd->pos,nd->pos+ndims);
    ndref=nd;  
    sgref=NULL;
    segExtId=-1;
  }

  template <typename T>
  newNodeT(NDskl_seg *sg, T ndims, int extId)
  {
    if (extId==0) 
      pos.assign(sg->pos,sg->pos+ndims);
    else if (extId==1) 
      pos.assign(sg->pos+ndims,sg->pos+2*ndims);
    else {fprintf(stderr,"ERROR in newNodeT, invalid extId\n");exit(0);}
    segExtId=extId;
    sgref=sg;
    ndref=NULL;
  }

  int fromSeg() const
  {
    return segExtId;
  }

  bool operator<(const newNodeT &b) const
  {
    if (this->pos.size()!=b.pos.size()) return (this->pos.size()<b.pos.size());

    for (unsigned long i=0;i<this->pos.size();i++)
      {
	if (this->pos[i]!=b.pos[i]) return (this->pos[i]<b.pos[i]);
      }
    return false;
  }

  bool operator==(const newNodeT &b) const
  {
    if (this->pos.size()!=b.pos.size()) return false;

    for (unsigned long i=0;i<this->pos.size();i++)
      {
	if (this->pos[i]!=b.pos[i]) return false;
      }
    return true;
  }

};

class skeletonBuildT
{
public:
  typedef std::set<newNodeT>::iterator newNodeIt;

  template <class outputIterator>
  static int getFilament(NDskl_seg *seg, outputIterator out)
  {
    bool dir=(seg->Prev==NULL);
    NDskl_seg *next=seg;
    do {
      *out = next;
      if (dir) next=next->Next;
      else next=next->Prev;      
    } while (next!=NULL);
    
    return (dir)?+1:-1;
  }
  
private:
  std::set<newNodeT> nodes;
  std::vector < std::vector<NDskl_seg *> > fil;
  std::vector < std::pair<newNodeIt,newNodeIt> > filExt; 

  void swapSeg(NDskl_seg* seg, int ndims, std::vector< std::pair<int,int> > &dataSwap, int orientation_ID=-1)
  {
    for (long i=0;i<ndims;i++) 
      std::swap(seg->pos[i],seg->pos[ndims+i]);
    
    for (unsigned long i=0;i<dataSwap.size();i++) 
      std::swap(seg->data[dataSwap[i].first],seg->data[dataSwap[i].second]);
    if (orientation_ID>=0) seg->data[orientation_ID]*=-1;
  }  
 
public:

  long nDims() {return (nodes.begin()==nodes.end())?0:nodes.begin()->getNDims();}
  long nNodes() {return nodes.size();}
  long nFils() {return fil.size();}
  long nSegs() 
  {
    unsigned long i=0;
    long ct=0;
    for (i=0;i<fil.size();i++) ct+=fil[i].size();
    return ct;
  }

  skeletonBuildT()
  {
    
  }

  void clear()
  {
    nodes.clear();
    fil.clear();
    filExt.clear();
  }

  template <typename T>
  newNodeIt addNode(NDskl_node *nd, T ndims)
  {
    //printf ("adding node %d: (%f,%f,%f)\n",nd->index,nd->pos[0],nd->pos[1],nd->pos[2]);
    std::pair<newNodeIt,bool> it=nodes.insert(newNodeT(nd,ndims));
    if (it.second) it.first->index=nodes.size()-1;
    return it.first;
  }

  template <typename T>
  newNodeIt addNode(NDskl_seg *sg, T ndims, int extId)
  {
    std::pair<newNodeIt,bool> it=nodes.insert(newNodeT(sg,ndims,extId));
    if (it.second) it.first->index=nodes.size()-1;
    return it.first;
  }

  void setFilamentSrc(newNodeIt n, int filId=-1)
  {
    //printf ("SRC : %ld",(filId<0)?filExt.size()-1:(long)filId);
    //if (n->ndref!=NULL) printf ("->%d",n->ndref->index);
    //printf("\n");
    if (filId<0) filExt.back().first=n;
    else filExt[filId].first=n;
  }

  void setFilamentDst(newNodeIt n, int filId=-1)
  {
    //printf ("DST : %ld",(filId<0)?filExt.size()-1:(long)filId);
    //if (n->ndref!=NULL) printf ("->%d",n->ndref->index);
    //printf("\n");
    if (filId<0) filExt.back().second=n;
    else filExt[filId].second=n;
  }

  long addFilament(newNodeIt src, newNodeIt dst)
  {
    if (src == dst) 
      {
	//printf ("*FIL : %ld, (%g,%g) -> (%g,%g)\n",filExt.size(),src->pos[0],src->pos[1],dst->pos[0],dst->pos[1]);
	return fil.size()-1;
      }
    fil.push_back(std::vector<NDskl_seg *>());
    filExt.push_back(std::make_pair(nodes.end(),nodes.end()));
    setFilamentSrc(src);
    setFilamentDst(dst);
    return fil.size()-1;
  }

  long addFilament()
  {    
    //printf ("*FIL : %ld\n",filExt.size());
    fil.push_back(std::vector<NDskl_seg *>());
    filExt.push_back(std::make_pair(nodes.end(),nodes.end()));
    return fil.size()-1;
  }

  
/*
  void setFilamentExt(std::pair<newNodeIt,newNodeIt> &n, int filId=-1)
  {
    if (filId<0) filExt.back()=n;
    else filExt[filId]=n;
  }
*/
  void AddSegToFilament(NDskl_seg *seg,int filId=-1)
  {
    long id=filId;
    if (id<0) id=fil.size()-1;
    fil[id].push_back(seg);    
  }

  template <class inputIT>
  void AddSegToFilament(inputIT start, inputIT stop,int filId=-1)
  {
    for (inputIT it=start;it!=stop;it++)
      {
	AddSegToFilament(*it,filId);
      }
  }

  NDskel *toNDskel(NDskel *model)
  {
    char comment[256];
    int orientation_id=getDataFieldID(model,1,ORIENTATION_TAG);
    long i,j,k,l;
    sprintf(comment,"%s (RB)",model->comment);
    NDskel *skl=Create_NDskel(model->dims,model->ndims,model->x0,model->delta,comment,nNodes(),nSegs());

    skl->nnodedata=model->nnodedata;
    skl->nodedata_info=(char**)malloc(sizeof(char *) * skl->nnodedata);
    for (i=0;i<skl->nnodedata;i++)
      {
	skl->nodedata_info[i]=(char*)malloc(1+strlen(model->nodedata_info[i]));
	strcpy(skl->nodedata_info[i],model->nodedata_info[i]);
      }
    skl->nodedata=(double*)malloc(sizeof(double)*skl->nnodedata*skl->nnodes);

    std::vector<int> nid2sid;    
    for (i=0;i<model->nnodedata;i++)
      {
	nid2sid.push_back(NDskel_SegDataIndex(model,model->nodedata_info[i]));
	if (nid2sid.back()<0)
	  {
	    char test[255];
	    sprintf(test,"%s_p1",model->nodedata_info[i]);
	    nid2sid.back()=NDskel_SegDataIndex(model,test);
	    sprintf(test,"%s_p2",model->nodedata_info[i]);
	    nid2sid.push_back(NDskel_SegDataIndex(model,test));
	  }
	else
	  {
	    nid2sid.push_back(NDskel_SegDataIndex(model,model->nodedata_info[i]));
	  }
      }

    std::vector< std::pair<int,int> > segDataSwap;
    for (i=0;i<model->nsegdata;i++)
      {
	char *w=strstr(model->segdata_info[i],"_p1");
	if (w!=NULL)
	  {
	    char tmp[255];
	    strReplace(model->segdata_info[i],tmp,"_p1","_p2");
	    int id=NDskel_SegDataIndex(model,tmp);
	    if ((id>=0)&&(i<id))
	      {
		segDataSwap.push_back(std::make_pair(i,id));
	      }
	  }
      }
    newNodeIt newNode = nodes.begin();
    for (i=0;i<skl->nnodes;i++,newNode++)
      {
	int index = newNode->index;
	NDskl_node *node=&skl->Node[index];
	for (j=0;j<skl->ndims;j++)
	  node->pos[j]=newNode->getCoord(j);
	node->type=newNode->getType();
	node->flags=newNode->getFlags();
	node->index=index;	

	node->data=&skl->nodedata[skl->nnodedata*index];
	double *data=newNode->getData();
	int segExt=newNode->fromSeg();
	if (segExt<0)
	  {
	    memcpy(node->data,data,sizeof(double)*skl->nnodedata);
	  }
	else
	  {
	    for (j=0;j<skl->nnodedata;j++)
	      if (nid2sid[2*j+segExt]<0)
		node->data[j]=-1;
	      else 
		node->data[j] = data[nid2sid[2*j+segExt]];
	  }
      }

    for (i=0;i<(long)fil.size();i++)
      {
	skl->Node[filExt[i].first->index].nnext++;
	skl->Node[filExt[i].second->index].nnext++;
      }

    for (i=0;i<skl->nnodes;i++)
      {
	NDskl_node *node= &skl->Node[i];
	node->nsegs=(int*)malloc(sizeof(int)*node->nnext);
	node->Seg=(NDskl_seg**)malloc(sizeof(NDskl_seg*)*node->nnext);
	node->Next=(NDskl_node**)malloc(sizeof(NDskl_node*)*node->nnext);
	node->nnext=0;
      }

    skl->nsegdata=model->nsegdata;
    skl->segdata_info=(char**)malloc(sizeof(char *) * skl->nsegdata);
    for (i=0;i<skl->nsegdata;i++)
      {
	skl->segdata_info[i]=(char*)malloc(1+strlen(model->segdata_info[i]));
	strcpy(skl->segdata_info[i],model->segdata_info[i]);
      }
    skl->segdata=(double*)malloc(sizeof(double)*skl->nsegdata*skl->nsegs);

    for (i=0,k=0;i<(long)fil.size();i++)
      {
	float oldPos[skl->ndims];

	NDskl_node *n1=&skl->Node[filExt[i].first->index];
	NDskl_node *n2=&skl->Node[filExt[i].second->index];
      
	n1->Next[n1->nnext]=n2;
	n2->Next[n2->nnext]=n1;	
	n1->nsegs[n1->nnext]=fil[i].size();
	n2->nsegs[n2->nnext]=fil[i].size();
	n1->Seg[n1->nnext]=&skl->Seg[k];
	n2->Seg[n2->nnext]=&skl->Seg[k+fil[i].size()-1];	
	memcpy(oldPos,n1->pos,skl->ndims*sizeof(float));
      	
	for (j=0;j<(long)fil[i].size();j++,k++)
	  {
	    bool doSwap=false;
	    NDskl_seg *seg=&skl->Seg[k];
	    memcpy(seg->pos,fil[i][j]->pos,2*skl->ndims*sizeof(float));
	  
	    for (l=0;l<skl->ndims;l++)
	      if (oldPos[l]!=seg->pos[l])
		{
		  doSwap=true;
		  break;
		}
	 
	    if (j==0) seg->Prev=NULL;
	    else seg->Prev=&skl->Seg[k-1];
	    if (j==(long)fil[i].size()-1) seg->Next=NULL;
	    else seg->Next=&skl->Seg[k+1];

	    seg->index=k;
	    seg->flags=fil[i][j]->flags;
	    seg->nodes[0]=(int)(filExt[i].first->index);
	    seg->nodes[1]=(int)(filExt[i].second->index);
	    
	    seg->data=&skl->segdata[skl->nsegdata*k];
	    memcpy(seg->data,fil[i][j]->data,sizeof(double)*skl->nsegdata);
	    
	    if (doSwap) swapSeg(seg,skl->ndims,segDataSwap,orientation_id);

	    memcpy(oldPos,&seg->pos[skl->ndims],skl->ndims*sizeof(float));
	  }

	n1->nnext++;n2->nnext++;
      }

    return skl;
  }
  
};



#endif
