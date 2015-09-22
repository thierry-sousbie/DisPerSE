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
#ifndef __MY__DELAUNAY_SETUP_H__
#define __MY__DELAUNAY_SETUP_H__

#include "vertex_info.hxx"

#if NDIMS==2
#if PERIODICITY==0

// NDIMS=2 and non-periodic
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_vertex_base_with_info_2<vertex_info,K>  Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>        Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds>           Dt;

//typedef Dt::Finite_vertices_iterator Finite_vertices_iterator;
//typedef Dt::Vertex_handle            Vertex_handle;
typedef Dt::Point                    Point;

#define FINITE_CELLS_BEGIN(what) ( what##->finite_faces_begin())
#define FINITE_CELLS_END(what) ( what##->finite_faces_end())

#else

// NDIMS=2 and periodic
// We actually use 3D tesselation ...
/*
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_triangulation_filtered_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Periodic_3_triangulation_ds_vertex_base_3.h>
#include <CGAL/Periodic_3_triangulation_ds_cell_base_3.h>
#include <CGAL/Periodic_3_triangulation_hierarchy_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_3_triangulation_filtered_traits_3<K> Pt;

typedef CGAL::Periodic_3_triangulation_ds_vertex_base_3<> Tvb;
typedef CGAL::Triangulation_vertex_base_3<Pt,Tvb> Vb;
typedef CGAL::Triangulation_vertex_base_with_info_3<vertex_info,Pt,Vb>   VbI;
typedef CGAL::Triangulation_hierarchy_vertex_base_3<VbI>  VbhI;

typedef CGAL::Periodic_3_triangulation_ds_cell_base_3<> Tcb;
typedef CGAL::Triangulation_cell_base_3<Pt,Tcb> Cb;

typedef CGAL::Triangulation_data_structure_3<VbI,Cb> Tds;

typedef CGAL::Periodic_3_Delaunay_triangulation_3<Pt,Tds> Dt_b;
//typedef CGAL::Periodic_3_triangulation_hierarchy_3<Dt_b>  Dt;

typedef Dt_b Dt;

typedef Dt::Point          Point;
typedef Dt::Iso_cuboid     Iso_cuboid;


#define FINITE_CELLS_BEGIN(what) (what##->finite_cells_begin())
#define FINITE_CELLS_END(what) (what##->finite_cells_end())
*/
#endif

#elif NDIMS==3

#if PERIODICITY==0

// NDIMS=3 and non-periodic
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
//#include <CGAL/Triangulation_hierarchy_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Location_policy.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_vertex_base_with_info_3<vertex_info,K>   Vb;
typedef Vb Vbh;
//typedef CGAL::Triangulation_hierarchy_vertex_base_3<Vb>  Vbh;
typedef CGAL::Triangulation_data_structure_3<Vbh>        Tds;

//typedef CGAL::Delaunay_triangulation_3<K,Tds,CGAL::Fast_location> Dt_b;
typedef CGAL::Delaunay_triangulation_3<K,Tds> Dt_b;
typedef Dt_b Dt;
//typedef CGAL::Triangulation_hierarchy_3<Dt_b>       Dt;

//typedef Dt::Finite_vertices_iterator Finite_vertices_iterator;
//typedef Dt::Vertex_handle            Vertex_handle;
typedef Dt::Point                    Point;

#define FINITE_CELLS_BEGIN(what) (what##->finite_cells_begin())
#define FINITE_CELLS_END(what) (what##->finite_cells_end())

#else

// NDIMS=3 and periodic
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_triangulation_filtered_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Periodic_3_triangulation_ds_vertex_base_3.h>
#include <CGAL/Periodic_3_triangulation_ds_cell_base_3.h>
#include <CGAL/Periodic_3_triangulation_hierarchy_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_3_triangulation_filtered_traits_3<K> Pt;

typedef CGAL::Periodic_3_triangulation_ds_vertex_base_3<> Tvb;
typedef CGAL::Triangulation_vertex_base_3<Pt,Tvb> Vb;
typedef CGAL::Triangulation_vertex_base_with_info_3<vertex_info,Pt,Vb>   VbI;
typedef CGAL::Triangulation_hierarchy_vertex_base_3<VbI>  VbhI;

typedef CGAL::Periodic_3_triangulation_ds_cell_base_3<> Tcb;
typedef CGAL::Triangulation_cell_base_3<Pt,Tcb> Cb;

typedef CGAL::Triangulation_data_structure_3<VbI,Cb> Tds;

typedef CGAL::Periodic_3_Delaunay_triangulation_3<Pt,Tds> Dt_b;
//typedef CGAL::Periodic_3_triangulation_hierarchy_3<Dt_b>  Dt;

typedef Dt_b Dt;

typedef Dt::Point          Point;
typedef Dt::Iso_cuboid     Iso_cuboid;


#define FINITE_CELLS_BEGIN(what) (what##->finite_cells_begin())
#define FINITE_CELLS_END(what) (what##->finite_cells_end())

#endif

#elif NDIMS>=4

// NDIMS>3, not fully implemented yet
#include <CGAL/Homogeneous_d.h>
#include <CGAL/Gmpzf.h>
#include <CGAL/Delaunay_d.h>

typedef CGAL::Gmpzf RT;
typedef CGAL::Homogeneous_d<RT> Kernel;
typedef CGAL::Delaunay_d<Kernel> Delaunay_d;

typedef Delaunay_d<NDIMS> Dt
typedef Dt::Point_d Point;
//typedef Delaunay_d::Simplex_handle Simplex_handle;
//typedef Delaunay_d::Vertex_handle Vertex_handle;

#endif
/*
#include <CGAL/Orthogonal_incremental_neighbor_search.h>

#if NDIMS==2

#include <CGAL/Search_traits_2.h>
typedef CGAL::Search_traits_2<K> SearchTraits;

#elif NDIMS==3

#include <CGAL/Search_traits_3.h>
typedef CGAL::Search_traits_3<K> SearchTraits;

#else

#include <CGAL/Search_traits_d.h>
typedef CGAL::Search_traits_d<K> SearchTraits;

#endif

typedef CGAL::Orthogonal_incremental_neighbor_search<SearchTraits> NN_incremental_search;
typedef NN_incremental_search::iterator NN_iterator;
typedef NN_incremental_search::Tree Tree;

*/


#endif
