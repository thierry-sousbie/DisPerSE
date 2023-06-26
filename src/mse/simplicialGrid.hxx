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
#include "cells.hxx"

template <class cellType = cellIdentity<> >
class simplicialGrid
{
 public:
  typedef cellType cellT;
  
  struct tesselation
  {
    enum type{empty,any1D,any2D,any3D,regular2D,alternate2D,regular3D,alternate3D};
  };

  typedef typename tesselation::type tesselationType;

  typedef typename cellType::typeT typeT;
  typedef typename cellType::idT idT;
  
 private:
  static const int NDIMS_MAX = 4;
  //static const int NFACEPERVOX_MAX = 4;

  static const long IDMASK = (1<<20)-1;
  static const long D[6];// initialization below
  static const long P1 = (1<<6)<<20;
  static const long P2 = (1<<7)<<20;
  static const long P12 = P1|P2;

  struct SimCell {
    long i[NDIMS_MAX];
    long id;
    
    template <class T>
    SimCell(std::vector<T> pi, long pid)
    {
      int j;
      for (j=0;j<pi.size();j++)
	i[j]=pi[j];
      id=pid;
    }

    template <class T>
    SimCell(T begin, T end, long pid)
    {
      int j=0;
      for (T it=begin;it!=end;it++)
	i[j++]=*it;
      id=pid;
    }

  };

  tesselationType curT;
  std::vector<long> id2index;
  std::vector<long> index2id[2][NDIMS_MAX+1];
  bool needParity;
  long curNFaces[2][NDIMS_MAX+1];
  long crtNFaces[NDIMS_MAX+1];
  std::vector<long> nFacesPerVox;
  std::vector<long> nFaces;
  int ndims;

  bool gridIsSet;
  std::vector<double> x0,dx,deltax,deltaxh;
  std::vector<long> nu,nuv;
  long nVox;

  std::vector<std::vector<std::vector<long> > > faceTable[NDIMS_MAX+1][2];
  std::vector<std::vector<std::vector<long> > > cofaceTable[NDIMS_MAX+1][2];

  template <class Iterator>
  long getBoundaryConfiguration(Iterator begin, Iterator end) const
  {
    long res=0;
    long u;
    int i=0;

    for (Iterator it=begin;it!=end;it++,i++)
      {
	u=*it;
	if (*it<0) u+=nu[i];
	else if (*it>nu[i]) u-=nu[i];

	if (u==0) res |= (1<<0)<<(2*i);
	else if (u==nu[i]-1) res |= (1<<1)<<(2*i);
      }

    return res;
  }
  
  int parityType(long id) const
  {
    int res=0;

    if (id2index[id]&P1) res|=(1<<0);
    if (id2index[id]&P2) res|=(1<<1);

    return res;
  }

  template <class Iterator>
  int parityType(Iterator begin, Iterator end) const
  {
    long u;
    int i=0;
    long res=0;
    
    for (Iterator it=begin;it!=end;it++,i++)
      {
	u=*it;
	if (*it<0) u+=nu[i];
	else if (*it>nu[i]) u-=nu[i];

	res+=u;
      }
    res=res&1;

    if (res) return (1<<0);
    else return (1<<1);
  }
  
  std::vector<SimCell> findSimilarCfg(cellType cell)
  {
    long type=cell.type();
    long ncomb;
    long ii[ndims+1];
    long i,j,c;
    long vindex,cid;
    std::vector<long> u,v;
    std::vector<SimCell> result;

    ncomb=1;
    for (i=0;i<ndims;i++)
      ncomb*=3;

    getReferenceId(cell,vindex,cid,u);
    v.resize(ndims);
        
    for (c=0;c<(1<<(1<<ndims));c++)
      {
	if (id2index[c]<0) continue;
	
	for (i=0,j=0;i<(1<<ndims);i++) if (c&(1<<i)) j++;
	if (j!=type+1) continue;
	
	for (i=0;i<ndims;i++) ii[i]=-1;
	ii[0]--;
	for (i=0;i<ncomb;i++)
	  {
	    j=0;ii[j]++;
	    while (ii[j]>1)
	      {	       
		ii[j]=-1;
		ii[++j]++;
	      };
	    for (j=0;j<ndims;j++) v[j]=u[j]+ii[j];

	    long index=getIndexFromId(v.begin(),v.end(),type,c);
	    if (index == cell.id())
	      {
		if (parityType(v.begin(),v.end()) & parityType(c))
		  result.push_back(SimCell(ii,ii+ndims,c));
	      }
	  }
      }

    return result;
  }

  void computeFaces(cellType cell, std::vector<cellType> &result) 
  {	   
    typeT type=cell.type();	  
    if (type==0) {result.resize(0);return;}
    result.resize(type+1);
    
    long vid[5];
    long vindex,cid;
    std::vector<long> u;
    getReferenceId(cell,vindex,cid,u);
    long i,j=0;
    for (i=0;i<5;i++) vid[i]=0;
    for (i=0;i<(1<<ndims);i++)
      if (cid&(1<<i)) vid[j++]=1<<i;
    /*
    printf("|%ld[",cell.id());
    for (i=0;i<(1<<ndims);i++) if (cid&(1<<i)) printf("%ld",i);
    printf("]");fflush(0);
    */ 
    //printf ("[(%ld,%ld):%ld/ %ld,%ld,%ld]",cell.id(),vindex,cid,u[0],u[1],u[2]);

    if (type==1) {
      result[0]=cellType(0,getIndexFromId(u.begin(),u.end(),0,vid[0]));
      result[1]=cellType(0,getIndexFromId(u.begin(),u.end(),0,vid[1]));
    }
    else if (type==2) {
      result[0]=cellType(1,getIndexFromId(u.begin(),u.end(),1,vid[0]+vid[1]));
      result[1]=cellType(1,getIndexFromId(u.begin(),u.end(),1,vid[0]+vid[2]));
      result[2]=cellType(1,getIndexFromId(u.begin(),u.end(),1,vid[1]+vid[2]));
      //printf("->[%ld,%ld,%ld]",result[0].id(),result[1].id(),result[2].id());
    }
    else if (type==3) {
      result[0]=cellType(2,getIndexFromId(u.begin(),u.end(),2,vid[0]+vid[1]+vid[2]));
      result[1]=cellType(2,getIndexFromId(u.begin(),u.end(),2,vid[0]+vid[1]+vid[3]));
      result[2]=cellType(2,getIndexFromId(u.begin(),u.end(),2,vid[0]+vid[2]+vid[3]));
      result[3]=cellType(2,getIndexFromId(u.begin(),u.end(),2,vid[1]+vid[2]+vid[3]));
    }
    else if (type==4) {
      result[0]=cellType(3,getIndexFromId(u.begin(),u.end(),3,vid[0]+vid[1]+vid[2]+vid[3]));
      result[1]=cellType(3,getIndexFromId(u.begin(),u.end(),3,vid[0]+vid[1]+vid[2]+vid[4]));
      result[2]=cellType(3,getIndexFromId(u.begin(),u.end(),3,vid[0]+vid[1]+vid[3]+vid[4]));
      result[3]=cellType(3,getIndexFromId(u.begin(),u.end(),3,vid[0]+vid[2]+vid[3]+vid[4]));
      result[4]=cellType(3,getIndexFromId(u.begin(),u.end(),3,vid[1]+vid[2]+vid[3]+vid[4]));
    }
  }  

  void computeCofaces(cellType cell, std::vector<cellType> &result) 
  {
    long i,j,k;
    typeT type=cell.type();
    long vindex,cid;
    std::vector<long> u,v;
    std::set<cellType> tmpr;

    if (type==ndims) {result.resize(0);return;}

    v.resize(ndims);
    getReferenceId(cell,vindex,cid,u);
    std::vector<SimCell> sim=findSimilarCfg(cell);
    //printf("grid :");
    for (i=0;i<sim.size();i++)
      {
	for (j=0;j<ndims;j++) v[j]=u[j]+sim[i].i[j];
	long vindex2=getIndexFromId(v.begin(),v.end(),0,1);
	//printf(" %ld-%ld (",vindex2,cell.id());
	for (j=vindex2*nFacesPerVox[type+1];j<(vindex2+1)*nFacesPerVox[type+1];j++)
	  {
	    std::vector<cellType> r;
	    computeFaces(cellType(type+1,j),r);
	   
	    for (k=0;k<r.size();k++)
	      {
	
		if (r[k]==cell) tmpr.insert(cellType(type+1,j));
	      }
	  }
	//printf(")");
      }
    //printf("\n");
    result.assign(tmpr.begin(),tmpr.end());
  }  

  
  void buildFaceCofaceTable()
  {
    long pref[5]={0,1,2,-1,0};
    long pid[NDIMS_MAX+1];
    std::vector<long> pu(ndims);
    long i,j,k,l,t,c,p;

    long ncomb=1;
    for (i=0;i<ndims;i++)
      ncomb*=4;
        
    for (t=0;t<=ndims;t++)
      {

	for (i=0;i<NDIMS_MAX+1;i++) pid[i]=0;
	pid[0]--;
	for (i=0;i<ncomb;i++)
	  {
	    j=0;pid[j]++;
	    while (pid[j]>3)
	      {
		pid[j]=0;
		pid[++j]++;
	      };

	    
	    long vid=0;
	    for (j=0;j<ndims;j++)
	      {
		pu[j]=pref[pid[j]];
		if (pu[j]<0) pu[j]+=nu[j];
		vid+=nuv[j]*pu[j];
	      }
	    
	    int parity=parityType(pu.begin(),pu.end())-1;
	    long btype=getBoundaryConfiguration(pu.begin(),pu.end());
	    
	    for (c=0;c<nFacesPerVox[t];c++)
	      {
		if (faceTable[t][parity][btype][c].size()) continue;
		if (cofaceTable[t][parity][btype][c].size()) continue;
		
		std::vector<cellType> rc,rf;
		cellType cell(t,vid*nFacesPerVox[t]+c);	    		    
		
		computeCofaces(cell,rc);
		computeFaces(cell,rf);
		
		faceTable[t][parity][btype][c].clear();
		cofaceTable[t][parity][btype][c].clear();
		
		for (l=0;l<rf.size();l++) 
		  faceTable[t][parity][btype][c].push_back(rf[l].id()-vid*nFacesPerVox[t-1]);
		for (l=0;l<rc.size();l++) 
		  cofaceTable[t][parity][btype][c].push_back(rc[l].id()-vid*nFacesPerVox[t+1]);
	      }
	  }
      }
  }
  
  public:

  std::vector<long> getNFaces(bool perVox=false)
  {
    if (perVox) return nFacesPerVox;
    else return nFaces;
  }

  template <class T>
  void index2coords(long vindex, std::vector<T> &u) const
  {
    int i;
    long cum=0;
    u.resize(ndims);

    for (i=ndims-1;i>=0;i--)
      {
	u[i]=(vindex-cum)/nuv[i];
	if (i) cum+=u[i]*nuv[i];
      }
  }

  template <class T>
  void index2position(long vindex, std::vector<T> &x) const
  {
    int i;
    std::vector<long> u;
    index2coords(vindex,u);
    x.resize(ndims);
    for (i=0;i<ndims;i++) x[i]=u[i]*dx[i]+x0[i];
  }

  void getReferenceId(cellType cell, long &vindex, long &cid,std::vector<long> &u) const
  {
    long type = cell.type();
    vindex = cell.id()/nFacesPerVox[type];
    index2coords(vindex,u);
    long ref=vindex*nFacesPerVox[type];
    cid = index2id[parityType(u.begin(),u.end())-1][type][cell.id()-ref];
  }

  template <class Iterator>
  long getIndexFromId(Iterator begin, Iterator end, char type, long id) const
  {
    long dd;
    long u;
    int i=0;
    long r1=0;
    long r2=(id2index[id]&IDMASK);
    
    for (Iterator it=begin;it!=end;it++,i++)
      {
	if (id2index[id]&D[i]) dd=nFacesPerVox[type]*nuv[i];
	else dd=0;

	u=*it;
	if (*it<0) u+=nu[i];
	else if (*it>nu[i]) u-=nu[i];

	if (dd&&(u==nu[i]-1)) dd=-u*dd;

	r1+=u*nuv[i];
	r2+=dd;	
      }
 
    return r1*nFacesPerVox[type]+r2;
  }

  long getIndexFromId(long u, long v, long w, char type, long id) const
  {
    std::vector<long> tmp;
    tmp.push_back(u);tmp.push_back(v);tmp.push_back(w);
    return getIndexFromId(tmp.begin(),tmp.end(),type,id);
  }

  void getCofaces(cellType cell, std::vector<cellType> &result) const
  {
    long i;
    typeT type=cell.type();
    long vindex,cid;
    std::vector<long> u;
    
    if (type==ndims) {result.resize(0);return;}
    getReferenceId(cell,vindex,cid,u);
 
    int parity=parityType(u.begin(),u.end())-1;
    long btype=getBoundaryConfiguration(u.begin(),u.end());
    long index = cell.id()-vindex*nFacesPerVox[type];
    
    const std::vector<long> &ft=cofaceTable[type][parity][btype][index];
    
    result.resize(ft.size());
    vindex*=nFacesPerVox[type+1];
    for (i=0;i<ft.size();i++) result[i]=cellType(type+1,vindex+ft[i]); 
  }
  
  void getFaces(cellType cell, std::vector<cellType> &result) const
  {
    long i;
    typeT type=cell.type();
    long vindex,cid;
    std::vector<long> u;

    if (type==0) {result.resize(0);return;}
    getReferenceId(cell,vindex,cid,u);
    
    int parity=parityType(u.begin(),u.end())-1;
    long btype=getBoundaryConfiguration(u.begin(),u.end());
    long index=cell.id()-vindex*nFacesPerVox[type];
    
    const std::vector<long> &ft=faceTable[type][parity][btype][index];
    
    result.resize(ft.size());
    vindex*=nFacesPerVox[type-1];
    for (i=0;i<ft.size();i++) result[i]=cellType(type-1,vindex+ft[i]);
  }

  void getVertice(cellType cell, std::vector<cellType> &result) const
  {
    long i;  
    char type=cell.type();

    if (type==0) {result.assign(1,cell);return;}
  
    long vindex,cid;
    std::vector<long> u;
    long j=0;

    getReferenceId(cell,vindex,cid,u);
    result.resize(type+1);  
    
    for (i=0;i<(1<<ndims);i++) 
      if (cid&(1<<i)) result[j++]=cellType(0,getIndexFromId(u.begin(),u.end(),0,(1<<i)));
  }

  float *getPosition(cellType cell,float *pos) const
  {
    long i,j;
    std::vector<cellType> v;
    std::vector<double> x(ndims);
    std::vector<double> xs(ndims);
    std::vector<double> xr(ndims);
     
    getVertice(cell,v);
    
    index2position(v[0].id(),xr);
    //if (cell.id() == 32513) printf(" [%g %g] ",xr[0],xr[1]);
    x=xr;

    for (i=1;i<v.size();i++)
      {
	index2position(v[i].id(),xs);
	//if (cell.id() == 32513) printf(" [%g %g] ",xs[0],xs[1]);
	for (j=0;j<xs.size();j++)
	  {
	    if (fabs(xs[j]-xr[j])>deltaxh[j])
	      {
		if (xr[j]<x0[j]+deltaxh[j]) xs[j]-=deltaxh[j]*2;
		else xs[j]+=deltaxh[j]*2;
	      }
	    x[j]+=xs[j];
	  }
      }//if (cell.id() == 32513) printf("\n");
    double inv=1./v.size();
    for (i=0;i<x.size();i++) pos[i]=x[i]*inv;

    return pos;
  }


private:
  
  void setRegular(int ndims)
  {
    if (ndims==1)
      {
	addVertice(1);
	addSimplex(0,1);
      }
    else if (ndims==2)
      {
	addVertice(2);
	addSimplex(0,1,2);
	addSimplex(1,2,3);
      }
    else if (ndims==3)
      {
	addVertice(3);
	addSimplex(0,7,1,3);
	addSimplex(0,7,3,2);
	addSimplex(0,7,2,6);
	addSimplex(0,7,6,4);
	addSimplex(0,7,4,5);
	addSimplex(0,7,5,1);
      }
  }

  /*
  void setRegular2D()
  {
    addVertice(2);
    addSimplex(0,1,2);
    addSimplex(1,2,3);
  }
  */

  void setAlternate2D()
  {
    addVertice(2);
    addSimplexP1(0,1,2);addSimplexP2(0,1,3);
    addSimplexP1(1,2,3);addSimplexP2(0,2,3);
  }

  void setAlternate3D()
  {
    addVertice(3);
    /*
    addSimplexP1(0,3);addSimplexP2(1,2);
    addSimplexP1(0,5);addSimplexP2(1,4);
    addSimplexP1(0,6);addSimplexP2(2,4);
    addSimplexP1(3,5);addSimplexP2(1,7);
    addSimplexP1(3,6);addSimplexP2(2,7);
    addSimplexP1(5,6);addSimplexP2(4,7);
    
    addSimplexP1(0,1,3);addSimplexP2(0,1,2);
    addSimplexP1(0,1,5);addSimplexP2(0,1,4);
    addSimplexP1(0,2,3);addSimplexP2(1,2,3);
    addSimplexP1(0,2,6);addSimplexP2(0,2,4);
    addSimplexP1(0,4,5);addSimplexP2(1,4,5);
    addSimplexP1(0,4,6);addSimplexP2(2,4,6);
    addSimplexP1(1,3,5);addSimplexP2(1,3,7);
    addSimplexP1(2,3,6);addSimplexP2(2,3,7);
    addSimplexP1(3,5,7);addSimplexP2(1,5,7);
    addSimplexP1(3,6,7);addSimplexP2(2,6,7);
    addSimplexP1(4,5,6);addSimplexP2(4,5,7);
    addSimplexP1(5,6,7);addSimplexP2(4,6,7);

    addSimplexP1(0,3,5);addSimplexP2(2,1,4);
    addSimplexP1(0,3,6);addSimplexP2(2,1,7);
    addSimplexP1(0,5,6);addSimplexP2(2,4,7);
    addSimplexP1(3,5,6);addSimplexP2(1,4,7);
    */
    addSimplexP1(0,1,3,5);addSimplexP2(2,0,1,4);
    addSimplexP1(0,2,3,6);addSimplexP2(2,1,3,7);
    addSimplexP1(0,3,5,6);addSimplexP2(2,1,4,7);
    addSimplexP1(0,4,5,6);addSimplexP2(2,4,6,7);
    addSimplexP1(3,5,6,7);addSimplexP2(1,4,5,7);
    /* 
    int i,j;
    printf("NFACES : \n");
    for (i=0;i<ndims+1;i++) 
      printf("   %d:(%ld<>%ld)\n",i,curNFaces[0][i],curNFaces[1][i]);
    printf("\n");
    int cnt[4]={0,0,0,0};
    
    for (i=0;i<256;i++)
      {
	if (id2index[i]>=0)
	  {
	    //printf("%2.2d: ",ct++);
	    int ct=0;
	    for (j=0;j<8;j++) if (i&(1<<j)) ct++;
	    if (ct!=2) continue;
	    //if (!((id2index[i]&P12)==P12)) continue;
	    for (j=0;j<8;j++) if (i&(1<<j)) {printf(" %d ",j);}
	    if (ct)
	      {
		cnt[ct-1]++;
		printf(": %d -> %ld(%d,%d)\n",cnt[ct-1],id2index[i]&IDMASK,
		       (id2index[i]&P1)!=0,(id2index[i]&P2)!=0);
	      }
	  }
      }
    */
      /*
    exit(0);
    */
  }

  long countVert()
  {
    int i,c;
    c=0;
    //if (ndims==0) return 0;
    for (i=0;i<(1<<ndims);i++) if (id2index[1<<i]!=-1) c++;
    return c;
  }

  long addVertex(bool du, bool dv, bool dw=0, bool dx=0)
  {
    long c=countVert();
 
    if (c==(1<<ndims))
      {
	ndims++;
	id2index.resize(1<<(1<<ndims),-1);
      }

    if (c>=(1<<NDIMS_MAX)) {
      fprintf(stderr,"ERROR in simpicialGrid, max dims =4.\n");
      exit(0);
    }

    if (id2index[(1<<c)]!=-1) {
      assert(0);
      fprintf(stderr,"ERROR in simpicialGrid, same vertex was added twice ! (%d,%d,%d,%d).\n",du,dv,dw,dx);
      exit(0);
    }

    if (c==0) {
      index2id[0][0].push_back(1<<0);
      index2id[1][0].push_back(1<<0);
    }

    if ((du==0)&&(dv==0)&&(dw==0)&&(dx==0))
      {
	curNFaces[0][0]++;
	curNFaces[1][0]++;
      }
    id2index[1<<c]=du*D[0]+dv*D[1]+dw*D[2]+dx*D[3] + P12;     
    return c;
  }  

  long addSimplexP(std::vector<long> p, long prty=0)
  {
    switch (p.size())
      {
      case 2:
	return addSimplexP(p[0],p[1],-1,-1,-1,prty);
	break;
      case 3:
	return addSimplexP(p[0],p[1],p[2],-1,-1,prty);
	break;
      case 4:
	return addSimplexP(p[0],p[1],p[2],p[3],-1,prty);
	break;
      case 5:
	return addSimplexP(p[0],p[1],p[2],p[3],p[4],prty);
	break;
      default:
	fprintf(stderr,"ERROR in addSimplex.\n");exit(0);
      }
    return -1;
  }
  
  long addSimplexP(long p1, long p2, long p3=-1, long p4=-1,long p5=-1, long prty=0)
  {
    long parity;
    long id=(1<<p1)+(1<<p2);
    int i,j;
    std::vector<long> sym;    
    std::vector<long> sym2;    
    std::vector<long> pts;  
 
    if (prty==0) parity=P12;
    else if (prty&1) parity=P1;
    else parity=P2;

    pts.push_back(p1);
    pts.push_back(p2);   

    if (p3!=-1) {
      id+=(1<<p3);
      pts.push_back(p3);
    }
    if (p4!=-1) {
      id+=(1<<p4);
      pts.push_back(p4);
    }
    if (p5!=-1) {
      id+=(1<<p5);
      pts.push_back(p5);
    }

    int type=pts.size()-1;
    long newIndex;
    
    if (id2index[id]>=0)
      {	
	id2index[id]|=parity;
	return id2index[id]&IDMASK;
      }
    
    if (pts.size()>2)
      {
	for (i=0;i<pts.size();i++)
	  {
	    std::vector<long> pts2;

	    for (j=0;j<pts.size();j++)
	      if (i!=j) pts2.push_back(pts[j]);

	     addSimplexP(pts2,prty);	
	  }
      }   

    newIndex = crtNFaces[type]++;
    newIndex |= parity; 
    id2index[id]=newIndex;

    return id2index[id]&IDMASK;
  }

  void build()
  {
    std::map<long,long> eq[6];
    std::map<long,long> ideq[6];
    typedef typename std::map<long,long>::iterator map_it;

    long i,j,k;

    for (i=0;i<1<<(1<<ndims);i++)
      {
	if (id2index[i]<0) continue;

	std::vector<int> pts;

	for (j=0;j<(1<<ndims);j++)
	  if (i&(1<<j)) pts.push_back(j);

	int type=pts.size()-1;
	std::vector<long> sym; 

	for (j=0;j<ndims;j++)
	  {
	    for (k=0;k<pts.size();k++)
	      if (!(id2index[1<<pts[k]]&D[j])) break;
	    if (k==pts.size()) sym.push_back(j);
	  }

	if (sym.size())
	  {
	    long id2=0;

	    for (j=0;j<sym.size();j++) id2index[i] |= D[sym[j]];

	    for (k=0;k<pts.size();k++)
	      {		
		long val=id2index[1<<pts[k]];
	
		for (j=0;j<sym.size();j++)
		  val &= ~D[sym[j]];

		for (j=0;j<(1<<ndims);j++) 
		  if (id2index[1<<j]==val) {id2|=(1<<j);break;}
	      }

	    if (id2index[id2]<0) {
	      printf("ERRORRRR \n");
	      continue;
	    }

	    eq[type].insert(std::make_pair(id2index[i]&IDMASK,id2index[id2]&IDMASK));
	    ideq[type].insert(std::make_pair(i,id2));
	    
	  }
      }

    for (i=0;i<1<<(1<<ndims);i++)
      {
	if (id2index[i]<0) continue;	
	//if ((id2index[i]&P12)!=P12) continue;

	int type=-1;
	for (j=0;j<(1<<ndims);j++)
	  if (i&(1<<j)) type++;

	long id=id2index[i]&IDMASK;

	map_it it=eq[type].find(id);
	if (it==eq[type].end()) continue;
	/*
	map_it it2=eq[type].find(it->second);
	while (it2!=eq[type].end())
	  {
	    it=it2;
	    it2=eq[type].find(it->second);
	  }
	*/
	id2index[i] = (id2index[i]&(~IDMASK)) + it->second;	
      }

    

    for (i=0;i<ndims+1;i++) 
      {
	eq[i].clear();
	crtNFaces[i]=0;
      }
    
    for (i=0;i<1<<(1<<ndims);i++)
      {
	if (id2index[i]<0) continue;		
	if ((id2index[i]&P12)!=P12) continue;

	int type=-1;
	for (j=0;j<(1<<ndims);j++)
	  if (i&(1<<j)) type++;	 

	long id=id2index[i]&IDMASK;
	if (eq[type].find(id) == eq[type].end())
	  eq[type].insert(std::make_pair(id,crtNFaces[type]++));
	
	id2index[i] = (id2index[i]&(~IDMASK)) + eq[type][id];
      }

   

    for (i=0;i<ndims+1;i++)
      {
	//printf("crt[%ld] = %ld\n",i,crtNFaces[i]);
	curNFaces[0][i]=curNFaces[1][i]=crtNFaces[i];
	eq[i].clear();
      }

    for (i=0;i<1<<(1<<ndims);i++)
      {
	if (id2index[i]<0) continue;		
	if (id2index[i]&P1) continue;

	int type=-1;
	for (j=0;j<(1<<ndims);j++)
	  if (i&(1<<j)) type++;
	
	if (ideq[type].find(i) != ideq[type].end()) continue;	
	
	long id=id2index[i]&IDMASK;

	if (eq[type].find(id) == eq[type].end())
	  eq[type].insert(std::make_pair(id,curNFaces[0][type]++));
	id2index[i] = (id2index[i]&(~IDMASK)) + eq[type][id];
      }

   
    for (i=0;i<ndims+1;i++) eq[i].clear();

    for (i=0;i<1<<(1<<ndims);i++)
      {
	if (id2index[i]<0) continue;	
	if (id2index[i]&P2) continue;

	int type=-1;
	for (j=0;j<(1<<ndims);j++)
	  if (i&(1<<j)) type++;
	
	if (ideq[type].find(i) != ideq[type].end()) continue;
	/*
	  {
	    id2index[i] = (id2index[i]&(~IDMASK)) + (id2index[ideq[type][i]]&(~IDMASK));
	    continue;	
	  }
	*/
	long id=id2index[i]&IDMASK;

	if (eq[type].find(id) == eq[type].end())
	  eq[type].insert(std::make_pair(id,curNFaces[1][type]++));
	id2index[i] = (id2index[i]&(~IDMASK)) + eq[type][id];
      }

 
    
    for (i=0;i<1<<(1<<ndims);i++)
      {
	if (id2index[i]<0) continue;	
	if ((id2index[i]&P12)==P12) continue;

	int type=-1;
	for (j=0;j<(1<<ndims);j++)
	  if (i&(1<<j)) type++;	
	
	if (ideq[type].find(i) != ideq[type].end())
	  {
	    /*
	    for (j=0;j<8;j++) if (i&(1<<j)) {printf(" %ld ",j);}
	    printf(" ---> ");
	    for (j=0;j<8;j++) if (ideq[type][i]&(1<<j)) {printf(" %ld ",j);}
	    printf("\n");
	    */
	    id2index[i] = (id2index[i]&(~IDMASK)) + (id2index[ideq[type][i]]&(IDMASK));
	    continue;	
	  }
      }
  
    for (i=0;i<2;i++)
      for (j=0;j<ndims+1;j++)
	index2id[i][j].resize(curNFaces[i][j]);

    for (i=0;i<1<<(1<<ndims);i++)
      {
	if (id2index[i]<0) continue;
	for (j=0;j<ndims;j++) if (id2index[i]&D[j]) break;
	if (j!=ndims) continue;
	//printf("HELLO\n");
	int type=-1;
	for (j=0;j<(1<<ndims);j++)
	  if (i&(1<<j)) type++;	
	
	int p=parityType(i);
	//printf("%ld : %d ->(%d,%ld)\n",i,p,type,id2index[i]&IDMASK);
	if (p&1) index2id[0][type][id2index[i]&IDMASK] = i;
	if (p&2) index2id[1][type][id2index[i]&IDMASK] = i;
      }

    
    /*******************   
    {
      printf("--------------------\n");
    int cnt[4]={0,0,0,0};
    for (i=0;i<(1<<(1<<ndims));i++)
      {
	if (id2index[i]>=0)
	  {
	    //printf("%2.2d: ",ct++);
	    int ct=0;
	    for (j=0;j<(1<<ndims);j++) if (i&(1<<j)) ct++;
	    if (ct!=3) continue;
	    //if (!((id2index[i]&P12)==P12)) continue;
	    for (j=0;j<(1<<ndims);j++) if (i&(1<<j)) {printf(" %ld ",j);}
	    if (ct)
	      {
		cnt[ct-1]++;
		printf(": %d -> %ld(%d,%d)\n",cnt[ct-1],id2index[i]&IDMASK,
		       (id2index[i]&P1)!=0,(id2index[i]&P2)!=0);
	      }
	  }
      }
    
    long l;
    j=2;
    //for (j=0;j<ndims+1;j++)
	{	
	  for (i=0;i<2;i++)
	    {
	      printf("dim %ld / P%ld:\n",j,i);
	      //k=1;
	      for (k=0;k<index2id[i][j].size();k++)
		{
		  printf("%ld : ",k);
		  printf (" [");
		  for (l=0;l<(1<<ndims);l++) if (index2id[i][j][k]&(1<<l)) {printf(" %ld ",l);}
		  printf("] \n");
		}
	      printf("\n");
	    }
	}
      exit(0);
    }
**********************/
  }
  
 public: 

  long addVertice(int ndims)
  {
    addVertex(0,0);
    addVertex(1,0);
    if (ndims>1)
      {
	addVertex(0,1);
	addVertex(1,1);
      }
    if (ndims>2)
      {
	addVertex(0,0,1);
	addVertex(1,0,1);
	addVertex(0,1,1);
	addVertex(1,1,1);
      }
    if (ndims>3)
      {
	addVertex(0,0,0,1);
	addVertex(1,0,0,1);
	addVertex(0,1,0,1);
	addVertex(1,1,0,1);
	addVertex(0,0,1,1);
	addVertex(1,0,1,1);
	addVertex(0,1,1,1);
	addVertex(1,1,1,1);
      }
    /*
    addSimplex(0,1);
    addSimplex(0,2);
    addSimplex(1,3);
    addSimplex(2,3);
    if (ndims>2) {
      addSimplex(0,4);
      addSimplex(1,5);
      addSimplex(2,6);
      addSimplex(3,7);
      addSimplex(4,5);
      addSimplex(4,6);
      addSimplex(5,7);
      addSimplex(6,7);
    }

    if (ndims>3)
      {
	addSimplex(8,9);
	addSimplex(8,10);
	addSimplex(9,11);
	addSimplex(10,11);
	addSimplex(8,12);
	addSimplex(9,13);
	addSimplex(10,14);
	addSimplex(11,15);
	addSimplex(12,13);
	addSimplex(12,14);
	addSimplex(13,15);
	addSimplex(14,15);

	addSimplex(0,8);
	addSimplex(1,9);
	addSimplex(2,10);
	addSimplex(3,11);
	addSimplex(4,12);
	addSimplex(5,13);
	addSimplex(6,14);
	addSimplex(7,15);
      }
    */
    return 0;
  }

  long addSimplex(long p1, long p2, long p3=-1, long p4=-1, long p5=-1)
  {
    return addSimplexP(p1,p2,p3,p4,p5,0);
  }

  long addSimplexP1(long p1, long p2, long p3=-1, long p4=-1, long p5=-1)
  {
    return addSimplexP(p1,p2,p3,p4,p5,1);
  }

  long addSimplexP2(long p1, long p2, long p3=-1, long p4=-1, long p5=-1)
  {
    return addSimplexP(p1,p2,p3,p4,p5,2);
  }
  
  void reset()
  {
    curT=tesselation::empty;
    id2index.resize(1<<(1<<0),-1);
    int i,j;
    for (i=0;i<2;i++)
      for (j=0;j<NDIMS_MAX+1;j++)
	index2id[i][j].clear();
    //index2id.clear();
    ndims=0;
    needParity=false;
    for (i=0;i<NDIMS_MAX+1;i++)
      {
	curNFaces[0][i]=curNFaces[1][i]=0;
	crtNFaces[i]=0;
      }

    resetGrid();
  }

  void setType(tesselationType t)
  {
    reset();
    if (t==tesselation::empty) return;
    
    switch (t)
      {
      case tesselation::any1D:
	setRegular(1);
	break;

      case tesselation::any2D:
      case tesselation::regular2D:
	setRegular(2);
	break;
	
      case tesselation::alternate2D:
	setAlternate2D();
	break;

      case tesselation::any3D:
      case tesselation::regular3D:
	setRegular(3);
	break;

      case tesselation::alternate3D:
	setAlternate3D();
	break;
	
      default:
	printf("ERROR, invalid type in simplicialGrid::setType()");
	exit(0);
      }
    build();

    /*****
        int i,j;
    printf("NFACES : \n");
    for (i=0;i<ndims+1;i++) 
      printf("   %d:(%ld<>%ld)\n",i,curNFaces[0][i],curNFaces[1][i]);
    printf("\n");
    int cnt[4]={0,0,0,0};
    
    for (i=0;i<256;i++)
      {
	if (id2index[i]>=0)
	  {
	    //printf("%2.2d: ",ct++);
	    int ct=0;
	    for (j=0;j<8;j++) if (i&(1<<j)) ct++;
	    if (ct!=2) continue;
	    //if (!((id2index[i]&P12)==P12)) continue;
	    for (j=0;j<8;j++) if (i&(1<<j)) {printf(" %d ",j);}
	    if (ct)
	      {
		cnt[ct-1]++;
		printf(": %d -> %ld(%d,%d)\n",cnt[ct-1],id2index[i]&IDMASK,
		       (id2index[i]&P1)!=0,(id2index[i]&P2)!=0);
	      }
	  }
      }

    exit(0);
    *****/
    curT=t;
  }

  template <typename FT, typename IT>
  void setGrid(std::vector<FT> _x0,
	       std::vector<FT> _deltax,
	       std::vector<IT> _nu)
  /*
  {
    std::vector<double> x0;
    std::vector<double> deltax;
    std::vector<long> nu;

    x0.assign(_x0.begin(),_x0.end());
    deltax.assign(_deltax.begin(),_deltax.end());
    nu.assign(_nu.begin(),_nu.end());
  }
  
  void setGrid(std::vector<double> _x0=std::vector<double>(),
	       std::vector<double> _deltax=std::vector<double>(),
	       std::vector<long> _nu=std::vector<long>())
  */
  {
    int i,j,k;

    if (_x0.size()==0) 
      {
	gridIsSet=false;
	x0.clear();
	deltax.clear();
	dx.clear();
	nu.clear();
	deltaxh.clear();
	nVox=0;
	nuv.clear();

	nFacesPerVox.clear();
	nFaces.clear();

	for (i=0;i<ndims+1;i++)
	  for (j=0;j<2;j++)
	    {
	      faceTable[i][j].clear();
	      cofaceTable[i][j].clear();
	    }
	
	gridIsSet=false;
	return;
      }

    x0.assign(_x0.begin(),_x0.end());
    deltax.assign(_deltax.begin(),_deltax.end());
    nu.assign(_nu.begin(),_nu.end());
    /*
    x0=_x0;
    deltax=_deltax;
    nu=_nu;
    */
    deltaxh=deltax;
    dx=deltax;
    for (i=0;i<deltax.size();i++)
      {
	deltaxh[i]/=2;
	dx[i]=deltax[i]/nu[i];
      }
    nVox=1;
    
    nuv=nu;
    nuv[0]=1;
    for (i=0;i<nu.size();i++)
      {
	nVox*=nu[i];
	if (i) nuv[i]=nuv[i-1]*nu[i-1];
      }

    
    
    /*
    nFacesPerVox.resize(ndims+1);

    for (i=0;i<(1<<(1<<ndims));i++)
      {
	int cv=0;
	if (id2index[i]<0) continue;
	for (j=0;j<(1<<ndims);j++) if (id2index[i]&(1<<j)) cv++;

	if ((id2index[i]&IDMASK)+1>nFacesPerVox[cv]) 
	  nFacesPerVox[cv]=1+(id2index[i]&IDMASK);
      }
    */
    //nFacesPerVox.assign(curNFaces[0],&curNFaces[0][ndims+2]);printf ("s:%ld\n",nFacesPerVox.size());
    nFacesPerVox.clear();
    nFaces.resize(ndims+1);
    for (i=0;i<=ndims;i++)
      {
	nFacesPerVox.push_back(curNFaces[0][i]);//printf ("s:%ld\n",nFacesPerVox.back());
	if (curNFaces[0][i]!=curNFaces[1][i]) 
	  {	    
	    fprintf(stderr,"ERROR in simplicalGrid: odd/even voxels do not have the same number of faces !\n");
	    assert(!(curNFaces[0][i]!=curNFaces[1][i]));
	    exit(0);
	  }
	
	nFaces[i]=nFacesPerVox[i]*nVox;
      }
    //printf("%ld %ld\n",nFaces[0],nFaces[1]);
    

    for (i=0;i<ndims+1;i++)
      for (j=0;j<2;j++)
	{
	  faceTable[i][j].resize((1<<(2*ndims))-1);
	  cofaceTable[i][j].resize((1<<(2*ndims))-1);
	  for (k=0;k<faceTable[i][j].size();k++)
	    {
	      faceTable[i][j][k].resize(nFacesPerVox[i]);
	      cofaceTable[i][j][k].resize(nFacesPerVox[i]);
	    }
	}

    
    buildFaceCofaceTable();
    gridIsSet=true;
    /*
    printf("NFACES : \n");
    for (i=0;i<ndims+1;i++) 
      printf("   %d:(%ld<>%ld) pv:%ld to:%ld\n",i,
	     curNFaces[0][i],curNFaces[1][i],
	     nFacesPerVox[i],nFaces[i]);
    printf("\n");
    */
    
  }

  void resetGrid()
  {
    setGrid(std::vector<double>(),std::vector<double>(),std::vector<long>());
  }
  
  
 public:

  simplicialGrid(tesselationType t=tesselation::empty,
		 std::vector<double> _x0=std::vector<double>(),
		 std::vector<double> _deltax=std::vector<double>(),
		 std::vector<long> _nu=std::vector<long>())
    {
      setType(t);
      setGrid(_x0,_deltax,_nu);
    }
  
  ~simplicialGrid()
    {
      
    }
  
   
};

template <class cellType>
const long simplicialGrid<cellType>::D[6] = {(1<<0)<<20,(1<<1)<<20,(1<<2)<<20,(1<<3)<<20,(1<<4)<<20,(1<<5)<<20};
