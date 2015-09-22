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
/************************************************************************/
/* An updtQueue is a class that stores pairs of <Key,value> and allows  */
/* retieving keys in order of their increasing corresponding values.    */
/* Pairs can be inserted, deleted and updated dynamically.              */
/* All operations are performed in logarithmic time, except retrieving  */
/* the lowest value key, which is a constant time operation.            */ 
/* Note that Keys are unique, but their values are not.                 */
/************************************************************************/

#ifndef UPDTQUEUE_HXX_
#define UPDTQUEUE_HXX_

#include <set>
#include <utility>

template < class Key, class Value, class CompareLess = std::less<Key> >
class updtQueue {

private:
    typedef std::pair<Key,Value> kv_pair;
    typedef std::pair<Value,Key> vk_pair;    

    class compareLessKV : public std::binary_function<kv_pair,kv_pair,bool> {
    public :
	bool operator()(const kv_pair& p1, const kv_pair& p2) const 
	    {	
		return CompareLess()(p1.first,p2.first); 
	    }
    };

    class compareLessVK : public std::binary_function<vk_pair,vk_pair,bool> {
    public :
	bool operator()(const vk_pair& p1, const vk_pair& p2) const
	    {
		//return CompareLess()(p1.second,p2.second);
		if (p1.first==p2.first) return CompareLess()(p1.second,p2.second);
		else return p1.first<p2.first;
	    }
    };

    typedef std::set<kv_pair,compareLessKV> kv_set;
    typedef std::set<vk_pair,compareLessVK> vk_set;
    
    typedef typename kv_set::iterator kv_iterator;
    typedef typename vk_set::iterator vk_iterator;
    
    typedef vk_set Base;

    kv_set kv;
    vk_set vk;

    typedef typename Base::iterator iterator;
    typedef typename Base::const_iterator const_iterator;
    typedef typename Base::reverse_iterator reverse_iterator;
    typedef typename Base::const_reverse_iterator const_reverse_iterator;
    typedef typename Base::size_type size_type;
    typedef typename Base::value_type value_type;

    kv_iterator kv_find(const Key& x) const {
	return kv.find(make_pair(x,Value()));

    }

    vk_iterator vk_find(const Key& x, const Value& v) const {
	return vk.find(make_pair(v,x));
    }

    vk_iterator vk_find(const Key& x) const {
	//return vk.find(make_pair(Value(),x));
	kv_iterator result = kv_find(x);

	if (result==kv.end()) return vk.end();

	return vk.find(make_pair(result->second,x));
    }

public:
    updtQueue() {}
    ~updtQueue() {}

    iterator begin () {return vk.begin();}
    const_iterator begin () const {return vk.begin();}

    iterator end () {return vk.end();}
    const_iterator end () const {return vk.end();}

    reverse_iterator rbegin() {return vk.rbegin();};
    const_reverse_iterator rbegin() const {return vk.rbegin();};
    
    reverse_iterator rend() {return vk.rend();};
    const_reverse_iterator rend() const {return vk.rend();};
    
    void clear ( ) {vk.clear();kv.clear();}

    iterator find (const Key& x) const {return vk.find(x);}

    size_type size() const {return vk.size();}

    Key pop() {
	vk_iterator vk_it = vk.begin();
	Key result = vk_it->second;

	kv.erase(make_pair(vk_it->second,vk_it->first));
	vk.erase(vk_it);

	return result;
    }

    iterator top ( ) const {
	return vk.begin();
    }

    iterator bottom ( ) const {
	iterator it = vk.end();

	if (size() == 0) return it;

	return --it;
    }

    iterator updateVal(const Key& x, Value val)
	{
	    kv_iterator kv_it = kv_find(x);

	    if (kv_it == kv.end()) return end();

	    vk_iterator vk_it = vk_find(x,kv_it->second);
	    
	    if (val == vk_it->first) return vk_it;
	    
	    kv_pair kvp = *kv_it;
	    vk_pair vkp = *vk_it;
	    
	    kvp.second = val;
	    vkp.first = val;

	    kv.erase(kv_it);
	    vk.erase(vk_it);

	    kv.insert(kvp);
	    iterator result = vk.insert(vkp).first;
	    /*
	    kv.insert(kv_it,kvp);
	    iterator result = vk.insert(vk_it,vkp);
	    
	    kv.erase(kv_it);
	    vk.erase(vk_it);
	    */
	    
	    return result;
	}

    std::pair<iterator,bool> insert(const Key& x, Value val)
	{	  
 
	    kv_iterator kv_it = kv_find(x);
	    
	    if (kv_it != kv.end())
	    {	
		return make_pair(updateVal(x,val),false);
	    }

	    kv.insert(make_pair(x,val));	    
	    return vk.insert(make_pair(val,x));	    
	}
    

};


#endif
