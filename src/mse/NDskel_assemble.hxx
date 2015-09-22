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
#ifndef __NDSKEL_ASSEMBLE_HXX__
#define __NDSKEL_ASSEMBLE_HXX__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include "float.h"
#include <list>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <limits>

#include <assert.h>

#include "NDsubSkeleton.hxx"
//#include "NDskel_trim.hxx"

#include "mytypes.h"
#include "NDskeleton.h"
#include "global.h"


class NDskelAssembleT
{
  //public: // this has to be changed into private ...
private:
  NDskel *skl;
  int valID;
  int valID_p1;
  int valID_p2;

  bool AllowFilCrossingTh;
  bool relaxAngularCondition;
  int minNNodesInFilament;


  double maxAngle;
  
  struct curBranchT
  {   
    std::vector<int> node;
    std::vector<NDskl_seg*> fil;
    std::map<int,double> score;
    double curScore;    
  };
  
  typedef std::list<curBranchT> branchesT;
  typedef branchesT::iterator branches_it;

  branchesT branch;
  std::vector<double> bestScoreAtNode;
  std::vector<unsigned char> forbiden;
  std::set<std::pair<int,int> > forbidenPairs;
  
  bool valueThreshold;

  double segSize(NDskl_seg *s)
  {
    double dl=0;
    for (int i=0;i<skl->ndims;i++)
      {
	double dx=s->pos[i]-s->pos[i+skl->ndims];
	dl+=dx*dx;
      }
    return sqrt(dl);
  }

  double nodeDist(NDskl_node *n1,NDskl_node *n2)
  {
    double dl=0;
    for (int i=0;i<skl->ndims;i++)
      {
	double dx=fabs(n1->pos[i]-n2->pos[i]);
	if (dx > 0.5*skl->delta[i]) dx = skl->delta[i]-dx;
	dl+=dx*dx;
      }
    
    return sqrt(dl);
  }

  double nodeDistAlongFil(NDskl_seg *seg, bool cutBelowThr=true)
  {

    long i;
    double dl=0;
    NDskl_seg *next=seg;
    bool skip=false;
    bool dir=(seg->Prev==NULL);
    do {
      if (cutBelowThr)
	{
	  skip=false;
	  if (next->data[valID_p1]<valueThreshold) skip=true;
	  if (next->data[valID_p2]<valueThreshold) skip=true;
	}
      if (!skip) dl +=segSize(next);
      if (dir) next=next->Next;
      else next=next->Prev;      
    } while (next!=NULL);

    
    return dl;
  }

  double nodeAngle(float *p1, float *p2, float *p3)
  {
    double s1,s2;
    double n1=0;
    double n2=0;
    double dot=0;
    
    for (int i=0;i<skl->ndims;i++) 
      {
	double s1=(p2[i]-p1[i]);
	double s2=(p3[i]-p2[i]);
	dot+=s1*s2;n1+=s1*s1;n2+=s2*s2;
      }
    dot/=sqrt(n1*n2);

    return acos(dot)*RAD2DEG;
    
  }

  double computeScore(curBranchT &b, NDskl_seg *nextFil=NULL)
  {
    double avg=0;
    double rms=0;
    double len=0;
    double lenT=0;
    double lenv=0;
    double lenvT=0;
    int nextNode=-1;

    if (nextFil!=NULL)
      {
	if (nextFil->Prev==NULL) nextNode=nextFil->nodes[1];
	else  nextNode=nextFil->nodes[0];
      }
           
    if (!relaxAngularCondition)
      {
	std::vector<float*> p;
	if (nextNode<0)
	  {
	    if (b.node.size()>2)
	      {
		if (nodeAngle(skl->Node[b.node[b.node.size()-3]].pos,
			      skl->Node[b.node[b.node.size()-2]].pos,
			      skl->Node[b.node[b.node.size()-1]].pos)>maxAngle) return -std::numeric_limits<double>::max();
	      }
	  }
	else if (b.node.size()>1)
	  {
	    if (nodeAngle(skl->Node[b.node[b.node.size()-2]].pos,
			  skl->Node[b.node[b.node.size()-1]].pos,
			  skl->Node[nextNode].pos)>maxAngle) return -std::numeric_limits<double>::max();
	  }
	  
	if ((nextNode>=0)&&(skl->Node[nextNode].type<=skl->ndims))
	  p.push_back(skl->Node[nextNode].pos);
	for (int i=b.node.size()-1;i>=0;i--)
	  {
	    if (skl->Node[b.node[i]].type<=skl->ndims)
	      p.push_back(skl->Node[b.node[i]].pos);
	    if (p.size()==3) break;
	  }
	if (p.size()==3)
	  {
	    if (nodeAngle(p[2],p[1],p[0])>maxAngle) return -std::numeric_limits<double>::max();
	  }
      }
      
    long nn;
    nn=0;
    for (int i=0;i<b.node.size();i++)
      {
	NDskl_node * node = &skl->Node[b.node[i]];
	  	  
	//if (node->data[valID]<0) return -1;
	if (i!=b.node.size()-1) 
	  {
	    NDskl_node * node2 = &skl->Node[b.node[i+1]];
	    double dl = nodeDist(node,node2);
	    double dlT = nodeDistAlongFil(b.fil[i]);
	    len+=dl;
	    lenT+=dlT;
	    lenv+=0.5*(node->data[valID]+node2->data[valID])*dl;
	    lenvT+=0.5*(node->data[valID]+node2->data[valID])*dlT;
	  }
	avg+=node->data[valID];	  
	nn++;
      }
      
    if (nextFil==NULL)
      avg/=nn;
    else
      {	  	  
	nn++;
	NDskl_node * node = &skl->Node[b.node.back()];
	NDskl_node * node2 = &skl->Node[nextNode];
	double dl = nodeDist(node,node2);
	double dlT = nodeDistAlongFil(nextFil);
	len+=dl;
	lenT+=dlT;
	lenv+=0.5*(node->data[valID]+node2->data[valID])*dl;
	lenvT+=0.5*(node->data[valID]+node2->data[valID])*dlT;
	avg=(avg+node2->data[valID])/(nn);
      }
      
    for (int i=0;i<b.node.size();i++)
      {
	NDskl_node * node = &skl->Node[b.node[i]];
	rms+=(node->data[valID]-avg)*(node->data[valID]-avg);
      }

    if (nextNode<0)
      rms=sqrt(rms)/nn;
    else
      rms=sqrt(rms+(skl->Node[nextNode].data[valID]-avg)*(skl->Node[nextNode].data[valID]-avg))/nn;     

    return lenT;
  }
  
  double computeScore(branches_it b, NDskl_seg *nextFil=NULL)
  {
    return computeScore(*b,nextFil);
  } 

  int ComputeBranches_rec(branches_it cur, int depth, int &maxDepth)
  {
    if (depth==maxDepth) return maxDepth;  
    NDskl_node *curNode=&skl->Node[cur->node.back()];

    if ((depth)&&(forbiden[curNode->index]>=2)) return depth;
    if ((!AllowFilCrossingTh)&&(forbiden[curNode->index]>=2)) return depth;

   
    NDskl_node *prevNode;
    branches_it killMe = branch.end();
    int Nin=0;
    int maxD=0;
    branches_it old_cur=cur;
    branches_it tmpCur=branch.insert(branch.end(),*cur);

    if (cur->node.size()==1)
      prevNode=NULL;
    else
      prevNode=&skl->Node[cur->node[cur->node.size()-2]];
    
    for (int i=0;i<curNode->nnext;i++)
      {
	if ((!depth)&&(forbiden[curNode->Next[i]->index]>=2)) continue;//return depth;
	if ((!AllowFilCrossingTh)&&(forbiden[curNode->Next[i]->index]>=2)) continue;

	if (forbiden[curNode->Next[i]->index])
	  {
	    std::pair<int,int> p(curNode->index,curNode->Next[i]->index);
	    if (p.first>p.second) std::swap(p.first,p.second);
	    if (forbidenPairs.find(p) != forbidenPairs.end()) continue;
	  }

	curBranchT newCur;
	
	if (curNode->Next[i] == prevNode) {/*printf("Prev\n");*/continue;}

	// prevents from looping
	if (old_cur->score.find(curNode->Next[i]->index)!=old_cur->score.end()) 
	  {continue;}
	    
	double curScore = computeScore(old_cur,curNode->Seg[i]); //compute score here;
	if (curScore==-std::numeric_limits<double>::max()) {/*printf("Score\n");*/continue;}
	

	cur=branch.insert(branch.end(),*tmpCur);
	bestScoreAtNode[curNode->Next[i]->index]=curScore;
	cur->node.push_back(curNode->Next[i]->index);
	cur->fil.push_back(curNode->Seg[i]);
	cur->score.insert(std::make_pair(curNode->Next[i]->index,curScore));
	cur->curScore=curScore;
	
	int dp = ComputeBranches_rec(cur,depth+1,maxDepth);
	if (dp>maxD) maxD=dp;
	
	Nin++;
      }
    
    branch.erase(tmpCur);
    if (Nin)
      {
	cur=old_cur;     
	for (branches_it it=(++cur);it!=branch.end();it++)
	  {
	 
	    if (it->curScore > cur->curScore)
	      cur=it;
	  }

	if (cur!=old_cur) 
	  *old_cur=*cur;

	branch.erase(++old_cur,branch.end());

      }
    
    
    if (depth>maxD) return depth;
    else return maxD;
  }

  void init(NDskl_node *startNode)
  {   
    branch.clear();  
    bestScoreAtNode.clear();

    if (startNode!=NULL)
      {	
	//forbiden.assign(skl->nnodes,false);
	branch.assign(1,curBranchT());
	branch.back().node.push_back(startNode->index);
	branch.back().score[startNode->index]=(double)computeScore(branch.back());
	branch.back().curScore=branch.back().score[startNode->index];
	bestScoreAtNode.assign(skl->nnodes,-std::numeric_limits<double>::max());
	bestScoreAtNode[startNode->index]=branch.back().curScore;
      }
  }

  
  int compute(double threshold, double angle, int maxDepth=40)  
  {       
    branchesT result;
    long i,j;
    long nloops=0;

    forbiden.assign(skl->nnodes,0);
    valueThreshold = threshold;
    for (i=0;i<skl->nnodes;i++)
      {
	if (skl->Node[i].data[valID]<threshold) forbiden[i]=2;
      }
    maxAngle=angle;
    j=0;
    relaxAngularCondition=false;
    while (true) {
      //if (angle<0) break;
      printf ("\rAssembling skeleton: %ld filament(s)",nloops++);fflush(0);
      curBranchT best;
      for (i=0;i<skl->nnodes;i++)
	{
	  if ((forbiden[i]>=2)&&(!AllowFilCrossingTh)) continue;
	  init(&skl->Node[i]);
	  ComputeBranches_rec(branch.begin(),0,maxDepth);	
	  if (branch.front().node.size()<minNNodesInFilament) continue;

	  if (best.node.size())
	    {
	      if (branch.front().curScore>best.curScore)
		best = branch.front();
	    }
	  else
	    best = branch.front();
	}
    
      if ((best.node.size()<minNNodesInFilament)||(best.curScore==-std::numeric_limits<double>::max()))
	{
	  break;
	}
      result.insert(result.end(),best);
  
      for (i=0;i<best.node.size();i++) 
	{
	  forbiden[best.node[i]]=1;
	  if (i) 
	    {
	      if (best.node[i-1]<best.node[i])
		forbidenPairs.insert(std::make_pair(best.node[i-1],best.node[i]));
	      else
		forbidenPairs.insert(std::make_pair(best.node[i],best.node[i-1]));
	    }
	}
    }

    branch.assign(result.begin(),result.end());     
    
    printf ("\rAssembling skeleton ... done. (found %ld filaments)\n",branch.size());
    return 0;
  }

  int ComputeBranches(NDskl_node *node, int maxDepth=40)
  {   
    init(node);
    return ComputeBranches_rec(branch.begin(),0,maxDepth);	
  }

  void reset()
  {
    branch.clear();  
    bestScoreAtNode.clear();
    forbiden.clear();
    skl=NULL;   
  }

  int getValueFieldID()
  {
    return valID;
  }

  void init(NDskel *skl_, std::string valueField = std::string(VALUE_TAG))
  {
    reset();

    skl=skl_;
    valID=NDskel_NodeDataIndex(skl,valueField.c_str());
    valID_p1=NDskel_SegDataIndex(skl,(valueField+std::string("_p1")).c_str());
    valID_p2=NDskel_SegDataIndex(skl,(valueField+std::string("_p2")).c_str());

    if ((valID_p1<0)||(valID_p2<0))
      {
	valID_p1=NDskel_SegDataIndex(skl,valueField.c_str());
	valID_p2=valID_p1;
      }

    if ((valID<0)||(valID_p1<0)||(valID_p2<0))
      {
	fprintf(stderr,"ERROR in NDskelAssembleT, data field '%s' is absent for nodes and/or segments.\n",valueField.c_str());
	exit(0);
      }
  }  
  
public:  

  // TODO: make this another class, used by NDskel_assemble !
  NDskel *trim(double threshold,int maxDepth=40)
  {
    return build(threshold,-1,maxDepth);
  }

  NDskel *build(double angle,int maxDepth=40)
  {
    return build(-std::numeric_limits<double>::max(),angle,maxDepth);
  }
 
  NDskel *build(double threshold, double angle,int maxDepth=40)
  {
    long i,j;

    //AllowFilCrossingTh=true;
    minNNodesInFilament=2;

    if (angle>=0) 
      compute(threshold,angle,maxDepth);
    else
      {
	// add all filaments above threshold
	for (i=0;i<skl->nnodes;i++)
	  {
	
	    for (j=0;j<skl->Node[i].nnext;j++)
	      {
		if (skl->Node[i].Next[j]->index<skl->Node[i].index) continue;
	
		curBranchT tmp;
		tmp.node.push_back(skl->Node[i].index);
		tmp.node.push_back(skl->Node[i].Next[j]->index);
		tmp.fil.push_back(&skl->Seg[i]);
		branch.push_back(tmp);
	      }
	  }
      }

    printf ("Processing ...");fflush(0);
    skeletonBuildT newSkl;
    std::list<NDskl_seg *> fil;
    std::vector<int> fil_dir;
    for (branches_it it=branch.begin();it!=branch.end();it++)
      {
	//skeletonBuildT::newNodeIt n1=newSkl.addNode(&skl->Node[it->node.front()],skl->ndims);	
	skeletonBuildT::newNodeIt na;
	skeletonBuildT::newNodeIt nb;
	std::list<NDskl_seg *>::iterator fil_beg;
	std::list<NDskl_seg *>::iterator fil_it;

	fil.clear();
	fil_dir.clear();
	for (i=0;i<it->node.size()-1;i++)
	  {
	    long size=fil.size();
	    int dir = skeletonBuildT::getFilament(it->fil[i],std::back_inserter(fil));
	    for (j=size;j<fil.size();j++) fil_dir.push_back(dir);
	  }
	if (!fil.size()) continue;
	
	na=newSkl.addNode(&skl->Node[it->node.front()],skl->ndims);
	nb=newSkl.addNode(&skl->Node[it->node.back()],skl->ndims);
	if (na!=nb) 
	  {
	    newSkl.addFilament(na,nb);
	    newSkl.AddSegToFilament(fil.begin(),fil.end());
	  }		
      }
    
    NDskel* ret=newSkl.toNDskel(skl);
    reset();

    printf("done.\n");
    return ret;
  }  
 
  NDskelAssembleT(NDskel *skl_, std::string valueField = std::string(VALUE_TAG))
  {
    init(skl_,valueField);
  }
  
  
 
};

#endif
