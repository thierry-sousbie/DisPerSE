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
#include <algorithm>

#ifdef USE_TBB
#include "tbb/parallel_invoke.h" 
#endif

namespace parallel {

template< class _Type >
inline void merge_ptr( const _Type* a_start, const _Type* a_end, const _Type* b_start, const _Type* b_end, _Type* dst )
{
    while( a_start < a_end && b_start < b_end ) {
        if ( *a_start <= *b_start ) *dst++ = *a_start++;    // if elements are equal, then a[] element is output
        else *dst++ = *b_start++;
    }
    while( a_start < a_end ) *dst++ = *a_start++;
    while( b_start < b_end ) *dst++ = *b_start++;
}

// This version is borrowed from "Introduction to Algorithms" 3rd edition, p. 799.
template< class _Type >
inline long my_binary_search( _Type value, const _Type* a, long left, long right )
{
    long low  = left;
    long high = __max( left, right + 1 );
    while( low < high )
    {
        long mid = ( low + high ) / 2;
        if ( value <= a[ mid ] ) high = mid;
        else low  = mid + 1; 
	// because we compared to a[mid] and the value was larger than a[mid].
	// Thus, the next array element to the right from mid is the next possible
	// candidate for low, and a[mid] can not possibly be that candidate.
    }
    return high;
}

#ifdef USE_TBB
// From Dr.Dobbs
template< class _Type >
inline void merge( _Type* t, long p1, long r1, long p2, long r2, _Type* a, long p3 )
{
    long length1 = r1 - p1 + 1;
    long length2 = r2 - p2 + 1;
    if ( length1 < length2 ) {
        std::swap(      p1,      p2 );
        std::swap(      r1,      r2 );
        std::swap( length1, length2 );
    }
    if ( length1 == 0 ) return;
    if (( length1 + length2 ) <= 16 )
        merge_ptr( &t[ p1 ], &t[ p1 + length1 ], &t[ p2 ], &t[ p2 + length2 ], &a[ p3 ] );
    else {
        long q1 = ( p1 + r1 ) / 2;
        long q2 = my_binary_search( t[ q1 ], t, p2, r2 );
        long q3 = p3 + ( q1 - p1 ) + ( q2 - p2 );
        a[ q3 ] = t[ q1 ];
        tbb::parallel_invoke(
        //Concurrency::parallel_invoke(
           [&] { merge( t, p1,     q1 - 1, p2, q2 - 1, a, p3     ); },
           [&] { merge( t, q1 + 1, r1,     q2, r2,     a, q3 + 1 ); }
        );
    }
}

template< class _Type, class Compare>
inline void mergeSort_rec( _Type* src, long l, long r, _Type* dst, bool srcToDst, Compare cmp=std::iterator_traits<_Type>::value_type::operator<())
{
    if ( r == l ) {    // termination/base case of sorting a single element
        if ( srcToDst )  dst[ l ] = src[ l ];    // copy the single element from src to dst
        return;
    }
    if (( r - l ) <= 64 && !srcToDst ) {     // 32 or 64 or larger seem to perform well
      std::sort( src + l, src + r + 1,cmp);    
      // want to do dstToSrc, can just do it in-place, just sort the src, no need to copy
      //stable_sort( src + l, src + r + 1 );  // STL stable_sort can be used instead, but is slightly slower than Insertion Sort
        return;
    }
    long m = (( r + l ) / 2 );
    tbb::parallel_invoke(             // Intel's     Threading Building Blocks (TBB)
    //Concurrency::parallel_invoke(       // Microsoft's Parallel Pattern Library  (PPL)
    [&] { mergeSort_rec( src, l,     m, dst, !srcToDst,cmp ); },      // reverse direction of srcToDst for the next level of recursion
    [&] { mergeSort_rec( src, m + 1, r, dst, !srcToDst,cmp ); }       // reverse direction of srcToDst for the next level of recursion
    );
    if ( srcToDst ) merge( src, l, m, m + 1, r, dst, l );
    else            merge( dst, l, m, m + 1, r, src, l );
}

template< class _Type , class Compare>
inline void mergeSort( _Type* src, long size, _Type* dst , Compare cmp=std::iterator_traits<_Type>::value_type::operator<())
{
  if (size<2) return;

  if (dst==void)
    {
      _type dst(*src);
      mergeSort_rec(src,0,size-1,dst,false,cmp);
    }
  else mergeSort_rec(src,0,size-1,dst,true,cmp);
}

#endif

  /*
template< class _Type , class Compare>
inline void mergeSort( _Type &src, Compare cmp=std::iterator_traits<_Type>::value_type::operator<())
{
  if (src.size()<2) return;
  mergeSort_rec(&(src[0]),src.size(),void,cmp);
}

template< class _Type , class Compare>
inline void mergeSort( _Type &src, _Type &dst, Compare cmp=std::iterator_traits<_Type>::value_type::operator<())
{
  if (src.size()<2) {
    dst.assign(src.begin(),src.end());
    return;
  }

  dst.resize(src.size());
  mergeSort_rec(&(src[0]),src.size(),&(dst[0]),cmp);
}
  */

template <class RandomAccessIterator, class Compare>
inline void sort(RandomAccessIterator first, RandomAccessIterator last, Compare cmp=std::iterator_traits<RandomAccessIterator>::value_type::operator<())
{

#ifdef USE_TBB
  mergeSort(&(*first), (long)(last-first), void, cmp);
#elif defined(_PARALLEL_ALGORITHM)
  __gnu_parallel::sort(v.begin(), v.end(),cmp);
#else
  std::sort(v.begin(),v.end(),cmp);
#endif

}


}
