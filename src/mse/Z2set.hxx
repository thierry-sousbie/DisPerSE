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
/*
 * Z2set.hxx
 *
 * Behaviour is identical to std::set except that inserting an element 
 * twice results in its removal
 *
 */

#ifndef Z2SET_HXX_
#define Z2SET_HXX_

#include <set>


template < class Key, class Compare = std::less<Key>,
           class Allocator = std::allocator<Key> >
class Z2set{
private:
  std::set<Key,Compare,Allocator> std_set;
    
public:
  typedef std::set<Key,Compare,Allocator> Base;
    
  typedef typename Base::iterator iterator;
  typedef typename Base::const_iterator const_iterator;
  typedef typename Base::reverse_iterator reverse_iterator;
  typedef typename Base::const_reverse_iterator const_reverse_iterator;
  typedef typename Base::size_type size_type;
    
  Z2set() {
	
  }
  ~Z2set() {
	
  }
    
  iterator begin () {return std_set.begin();};
  const_iterator begin () const {return std_set.begin();};
    
  iterator end () {return std_set.end();};
  const_iterator end () const {return std_set.end();};
    
  reverse_iterator rbegin() {return std_set.rbegin();};
  const_reverse_iterator rbegin() const {return std_set.rbegin();};
    
  reverse_iterator rend() {return std_set.rend();};
  const_reverse_iterator rend() const {return std_set.rend();};
    
  void clear () {std_set.clear();}
    
  size_type count ( const Key& x ) const {return std_set.count();}

  void erase ( iterator position ) {std_set.erase(position);}
  size_type erase ( const Key& x ) {return std_set.erase(x);}
  void erase ( iterator first, iterator last ) {std_set.erase(first,last);}

  iterator find ( const Key& x ) const {return std_set.find(x);}

  size_type size() const {return std_set.size();}

  std::pair<iterator,bool>
  insert ( const Key& x )
  {
    std::pair<iterator,bool> val;
    val = std_set.insert(x);
    if (val.second == false)
      {
	iterator next = val.first;
	next++;
	std_set.erase(val.first);
	val.first = next;
	//val.first = std_set.end();
	return val;
      }
    else return val;
  }

  iterator
  insert ( iterator position, const Key& x )
  {
    if (*position == x) {
      std_set.erase(position);
      return std_set.end();
    }
    else return std_set.insert(x);
  }

  template <class InputIterator>
  void
  insert ( InputIterator first, InputIterator last )
  {
    for (InputIterator it=first;it!=last;it++)
      insert(*it);
  }

};

#endif /* Z2SET_HXX_ */
