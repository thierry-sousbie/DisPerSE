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
#ifdef GLOBAL_DEFINITION
#define GLOBAL
#else 
#define GLOBAL extern
#endif

#ifndef VERSION_STR
#define STRINGIFY1(x)  #x
#define STRINGIFY(x)  STRINGIFY1(x)
#define VERSION_STR STRINGIFY(VER_MAJOR) "." STRINGIFY(VER_MINOR) "." STRINGIFY(VER_BUILD) 
#endif

#ifndef VALUE_TAG
#define VALUE_TAG "field_value"
#endif

#ifndef MASS_TAG
#define MASS_TAG "mass"
#endif

#ifndef PARENT_TAG
#define PARENT_TAG "parent_index"
#endif

#ifndef PARENT_LOG_TAG
#define PARENT_LOG_TAG "parent_log_index"
#endif

#ifndef PERSISTENCE_PAIR_TAG
#define PERSISTENCE_PAIR_TAG "persistence_pair"
#endif

#ifndef PERSISTENCE_TAG
#define PERSISTENCE_TAG "persistence"
#endif

#ifndef PERSISTENCE_RATIO_TAG
#define PERSISTENCE_RATIO_TAG "persistence_ratio"
#endif

#ifndef PERSISTENCE_NSIG_TAG
#define PERSISTENCE_NSIG_TAG "persistence_nsigmas"
#endif

#ifndef ROBUSTNESS_TAG
#define ROBUSTNESS_TAG "robustness"
#endif

#ifndef ROBUSTNESS_RATIO_TAG
#define ROBUSTNESS_RATIO_TAG "robustness_ratio"
#endif

#ifndef ARC_ID_TAG
#define ARC_ID_TAG "arc_id"
#endif

#ifndef FILAMENT_ID_TAG
#define FILAMENT_ID_TAG "filament_id"
#endif
        
#ifndef ORIENTATION_TAG
#define ORIENTATION_TAG "orientation"
#endif

#ifndef TYPE_TAG
#define TYPE_TAG "type"
#endif

#ifndef INDEX_TAG
#define INDEX_TAG "index"
#endif

#ifndef TRUE_INDEX_TAG
#define TRUE_INDEX_TAG "true_index"
#endif

#ifndef CELL_TAG
#define CELL_TAG "cell"
#endif

#ifndef BOUNDARY_TAG
#define BOUNDARY_TAG "boundary"
#endif

#ifndef DELAUNAY_TESSELATION_TAG
#define DELAUNAY_TESSELATION_TAG "delaunay_tesselation"
#endif

#ifndef CRITICAL_INDEX_TAG
#define CRITICAL_INDEX_TAG "critical_index"
#endif

#ifndef LENGTH_TAG
#define LENGTH_TAG "length"
#endif

#ifndef LOG_TAG
#define LOG_TAG(tag) "log_" tag
#endif

#ifndef UP_TAG
#define UP_TAG(tag) "up_" tag
#endif

#ifndef DOWN_TAG
#define DOWN_TAG(tag) "down_" tag
#endif

#ifndef SEG_P1_TAG
#define SEG_P1_TAG(tag) tag "_p1"
#endif

#ifndef SEG_P2_TAG
#define SEG_P2_TAG(tag) tag "_p2"
#endif


#ifndef SOURCE_TAG
#define SOURCE_TAG(tag) "source_" tag
#endif

GLOBAL int verbose;
GLOBAL int debug_dump;
GLOBAL int glob_num_threads;
GLOBAL int glob_num_omp_threads;

#undef GLOBAL
