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
#ifndef __SKELEX_IO
#define __SKELEX_IO

#include <stdlib.h>
#include "gadget_struct.h"
#include "myendian.h"

#ifdef __cplusplus
extern "C" {
#endif
    
#define FLAG_POS  ((int)1<<0)
#define FLAG_VEL  ((int)1<<1)
#define FLAG_MASS ((int)1<<2)
#define FLAG_ID   ((int)1<<3)
#define FLAG_TYPE ((int)1<<4)
#define FLAG_GAS  ((int)1<<5)
#define FLAG_SUP  ((int)1<<6)
#define FLAG_ALL  (((int)1<<6)-1)
#define FLAG_DM   ((int)FLAG_POS|FLAG_VEL|FLAG_MASS|FLAG_ID)
//if the file is not of the same endian type as the computer
//this is done automatically
#define FLAG_SWAPENDIAN ((int)1<<31)
    
#define SKIP(filename) {int skpdummy;skpdummy=fread(&skpdummy,sizeof(skpdummy),1,filename);}

    int IsGadgetFile(const char *);
    int ReadGadget(const char*,snapshot_data*,int);
    void freegadgetstruct(snapshot_data *snap);
    int WriteGadget(const char *fname, snapshot_data *gadget, int *selection, int selection_size);
    
#ifdef __cplusplus
}
#endif


#endif
