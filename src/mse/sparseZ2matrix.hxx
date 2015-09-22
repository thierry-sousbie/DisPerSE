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
#ifndef SPARSE_Z2_MATRIX_HXX_
#define SPARSE_Z2_MATRIX_HXX_

#include <list>
#include <set>
#include <string>
#include <vector>
#include <algorithm>

#ifdef USE_THREADS
#include "pthreadBarrier.hxx"
#endif

class sparseZ2Matrix {
public:
  enum initValue{zero, unity};
  typedef unsigned int idType;  
  
#ifdef USE_THREADS
private:  
  pthread_mutex_t colLock_mutex;
  typedef std::vector<bool> lockType;
  lockType lockedCol;
public:
  typedef lockType::iterator lockItType;

#endif

private:
  typedef idType rowEl_type;  
  typedef std::vector<rowEl_type> row_type;  
  typedef row_type::iterator row_it_type;

  typedef std::set<row_it_type> colListEl_type;
  typedef colListEl_type::iterator colListEl_it_type;
  typedef std::list<colListEl_type> colList_type; 
  typedef colList_type::iterator colList_it_type;

  typedef colList_it_type colEl_type; 
  typedef std::vector<colEl_type> col_type;  
  typedef col_type::iterator col_it_type;

  typedef std::pair<colList_it_type,colListEl_it_type> element_it;
  
  bool orderedRow;
  row_type row;
  col_type col;
  colList_type colList;
  colListEl_type dummyColListEl;

  bool allocAllCol;

public:
  typedef colList_type colListType;
  typedef colEl_type colType;
  typedef colListEl_it_type colItType; 

  bool colIsVoid(const colType col) const 
  {
    if (col == colList.end()) return true;
    if (col->begin()==col->end()) return true;
    return false;
  }
  
  
private:
  colListEl_it_type end() const {return dummyColListEl.end();} 

  element_it getRef(idType rowId, idType colId)
  {
    colEl_type curCol = col[colId];
    
    if (curCol==colList.end()) return element_it(colList.end(),end());
    
    row_it_type ref = row.begin()+rowId;
    colListEl_it_type cit = curCol->find(ref);
    
    if (cit==curCol->end()) return element_it(curCol,end());
    else return element_it(curCol,cit);
  }

  void insert(idType rowId,idType colId)
  {
    colList.push_front(colListEl_type());
    colList.front().insert(row.begin()+rowId);
    col[colId]=colList.begin();
  }

  void insert(idType rowId, colList_it_type it)
  {
    it->insert(row.begin()+rowId);
  }

  void insertCol(idType colId, colListEl_type &val)
  {
     colList.push_front(val);
     col[colId]=colList.begin();
  }

  void erase(element_it &it,idType colId)
  {
    it.first->erase(it.second);

    if (it.first->size()==0) 
      {
	if (allocAllCol) return;
	col[colId]=colList.end();
	colList.erase(it.first);
      }
  }

  void eraseCol(idType colId)
  {
    if (allocAllCol) return;
    colEl_type it = col[colId];
    if (it==colList.end()) return;
    col[colId]=colList.end();
    colList.erase(it);
  }

public:

  sparseZ2Matrix(idType Nrow, idType Ncol, initValue val = zero, bool AllocateCols = false,bool locksInit=false)
  {
    idType i;
    
    row.resize(Nrow);
    for (i=0;i<Nrow;i++) row[i]=i;
    col.assign(Ncol,colList.end());
    orderedRow=true;
    allocAllCol=AllocateCols;
    if (allocAllCol)
      {
	for (i=0;i<Ncol;i++)
	  {
	    colList.push_front(colListEl_type());
	    col[i]=colList.begin();
	  }
      }
    if (val==unity)
      {
	for (i=0;i<Nrow;i++) insert(i,i);
      }
#ifdef USE_THREADS
    pthread_mutex_init(&colLock_mutex, NULL);
    if (locksInit) initLocks();
#endif
  }

  ~sparseZ2Matrix()
  {
#ifdef USE_THREADS
    pthread_mutex_destroy(&colLock_mutex);
#endif
  }

  long nRows() const {return row.size();}
  long nCols() const {return col.size();}

  void initLocks()
  {
#ifdef USE_THREADS
    lockedCol.assign(col.size(),false);
#endif
  }

  // returns whether the value was already set to the same value
  bool set(idType rowId, idType colId, bool value=true)
  {
    element_it ref=getRef(rowId, colId);

    if (ref.first==colList.end())
      {
	if (!value) return true;
	insert(rowId,colId);
	return false;
      }
    
    if (ref.second==end())
      {
	if (!value) return true;
	insert(rowId,ref.first);
	return false;
      }

    if (value) return true;
    erase(ref,colId);
    return false;  
  }

  bool get(idType rowId, idType colId)
  {
    if (getRef(rowId, colId).second == end()) return false;
    return true;
  }

  sparseZ2Matrix::colType getCol(idType colId) const
  {
    return col[colId];
  }
  
  sparseZ2Matrix::colType nullCol()
  {
    return colList.end();
  }

  void clearCol(idType colId)
  {  
    eraseCol(colId);
  }

  colItType addCol(idType refColId, idType addColId, colItType cur)
  {
    colEl_type rCol = col[refColId];
    colEl_type aCol = col[addColId];
    
    if (aCol==colList.end()) return cur;
    if (rCol==colList.end()) 
      {
	insertCol(refColId, *aCol);
	return cur;
      }
    
    colListEl_it_type r=rCol->begin();
    colListEl_it_type a=aCol->begin();
    
    for (;a!=aCol->end();a++)
      {
	for (;*r<*a;r++) 
	  if (r==rCol->end()) break;
	
	if (r==rCol->end())
	  {
	    rCol->insert(a,aCol->end());
	    break;
	  }
	
	if (*r==*a) 
	  {
	    colListEl_it_type tmp=r++;
	    if (cur==tmp) cur=r;
	    rCol->erase(tmp);
	  }
	else rCol->insert(*a);
      }
    
   
    if (rCol->size()==0) 
      eraseCol(refColId);

    return cur;
  }

  void addCol(idType refColId, colType aCol)
  {
    colEl_type rCol = col[refColId];
    
    if (aCol==colList.end()) return;
    if (rCol==colList.end()) 
      {
	insertCol(refColId, *aCol);
	return;
      }
    
    colListEl_it_type r= rCol->begin();
    colListEl_it_type a= aCol->begin();
    
    for (;a!=aCol->end();a++)
      {
	for (;*r<*a;r++) 
	  if (r==rCol->end()) break;
	
	if (r==rCol->end())
	  {
	    rCol->insert(a,aCol->end());
	    break;
	  }
	
	if (*r==*a) 
	  {
	    colListEl_it_type tmp=r++;
	    rCol->erase(tmp);
	  }
	else rCol->insert(*a);
      }
    
   
    if (rCol->size()==0) eraseCol(refColId);    
  }

  void addCol(idType refColId, idType addColId)
  {
    addCol(refColId,col[addColId]);
  }


  void swapRow(idType a, idType b)
  {
    rowEl_type tmp = row[a];
    orderedRow=false;
    row[a]=row[b];
    row[b]=tmp;
  }

  void swapCol(idType a, idType b)
  {
    colEl_type tmp = col[a];
    col[a]=col[b];
    col[b]=tmp;
  }

  void reverse()
  {
    colEl_type tmp;
    idType N=(col.size()/2);
    idType i,j;

    for (i=0,j=col.size()-1; i<N;) swapCol(i++,j--);
  }

  long low(idType colId) const
  {
    colEl_type tmp = col[colId];
    if (tmp==colList.end()) return -1;
    if (tmp->size()==0) return -1;

    if (orderedRow) return (long)(*(*tmp->rbegin()));
    
    long row=-1;
    for (colListEl_it_type it=tmp->begin(); it!=tmp->end(); it++)
      {
	if (*(*it)>row) row=*(*it);
      }
    return row;
  }

  long high(idType colId) const
  {
    colEl_type tmp = col[colId];
    if (tmp==colList.end()) return -1;
    if (tmp->size()==0) return -1;

    if (orderedRow) return (long)(*(*tmp->begin()));
    
    long row=-1;
    for (colListEl_it_type it=tmp->begin(); it!=tmp->end(); it++)
      {
	if (*(*it)<row) row=*(*it);
      }
    return row;
  }

  void print(std::string txt = std::string("M")) const
  {
    int c,r,i,j;
    std::string text(txt);
    std::vector<bool> mat(row.size()*col.size(),false);
    
    for (c=0;c<col.size();c++)
      {
	if (col[c] != colList.end())
	  {
	    for (colListEl_it_type it=col[c]->begin(); it!=col[c]->end(); it++)
	      {
		mat[*(*it)+row.size()*c]=true;
	      }
	  }
      }
    
    printf("\n");    
    //for (i=0;i<text.length();i++) printf(" ");

    for (i=0;i<3;i++) {
      for (j=0;j<text.length()+4;j++) printf(" ");
      printf("    |");
      for (c=0;c<col.size();c++) printf("%d",(c/(int)(pow(10,2-i)))%10);
      printf("|\n");
    }
    
    for (i=0;i<text.length()+4;i++) printf(" ");
    printf("    -");
    for (c=0;c<col.size();c++) printf("-");
    printf("-\n");


    for (r=0;r<row.size();r++)
      {
	if (r==row.size()/2) printf(" %s = |%3.3d|",text.c_str(),r);
	else {
	  for (i=0;i<text.length();i++) printf(" ");
	  printf("    |%3.3d|",r);
	}

	for (c=0;c<col.size();c++)
	  printf("%d",(int)mat[r+c*row.size()]); 

	printf("|\n");
      }
    for (i=0;i<text.length()+4;i++) printf(" ");
    printf("    \\");
    for (c=0;c<col.size();c++) printf(" ");
    printf("/\n");
    
  }

  bool lockCol(long colId)
  {
#ifdef USE_THREADS
    pthread_mutex_lock(&colLock_mutex);
    bool res=lockedCol[colId];
    if (res) {
      pthread_mutex_unlock(&colLock_mutex);
      return false;
    }
    lockedCol[colId]=true;
    pthread_mutex_unlock(&colLock_mutex);
    return true;

#endif
    return true;
  }

  bool lockCol(long colId1, long colId2)
  {
#ifdef USE_THREADS
    pthread_mutex_lock(&colLock_mutex);
    bool res1=lockedCol[colId1];
    bool res2=lockedCol[colId2];
    if (res1||res2) {
      pthread_mutex_unlock(&colLock_mutex);
      return false;
    }
    lockedCol[colId1]=true;
    lockedCol[colId2]=true;
    pthread_mutex_unlock(&colLock_mutex);
    return true;
#endif
    return true;
  }

  void unlockCol(long colId)
  {
#ifdef USE_THREADS
    pthread_mutex_lock(&colLock_mutex);
    lockedCol[colId]=false;
    pthread_mutex_unlock(&colLock_mutex);
#endif
  }

  void unlockCol(long colId1,long colId2)
  {
#ifdef USE_THREADS
    pthread_mutex_lock(&colLock_mutex);
    lockedCol[colId1]=false;
    lockedCol[colId2]=false;
    pthread_mutex_unlock(&colLock_mutex);
#endif
  }

};

#endif
