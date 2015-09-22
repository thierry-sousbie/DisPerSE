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
#include <stdlib.h>
#include <stdio.h>
#include "fof_io.h"

group_list *Load_FOF(char *FileName)
{
    group_list *glist;
    group_type *gr;
    int i,j;
    float tmp[20];
    FILE *f;
    int ret;

    if ((f=fopen(FileName,"r"))==NULL)
    {
	fprintf(stderr,"ERROR: Cannot open file %s .\n",FileName);
	return NULL;
    }

    printf ("Loading FOF groups from file %s ... ",FileName);

    glist = calloc(1,sizeof(group_list));

    ret=fread(&i,sizeof(int),1,f);
    ret=fread(&glist->ndims,sizeof(int),1,f);
    ret=fread(&glist->NGroups,sizeof(int),1,f);
    ret=fread(&glist->NPart,sizeof(int),1,f);
    ret=fread(&glist->NMin,sizeof(int),1,f);
    ret=fread(&glist->ll,sizeof(float),1,f);
    ret=fread(&i,sizeof(int),1,f);
    
    glist->group = calloc(glist->NGroups,sizeof(group_type));
    for (j=0;j<glist->NGroups;j++)
    {
	gr = &glist->group[j];
	ret=fread(&i,sizeof(int),1,f);
	ret=fread(&gr->N,sizeof(int),1,f);
	ret=fread(&gr->cut,sizeof(int),1,f);
	
	gr->Center = calloc (sizeof(float),glist->ndims);
	gr->Spin = calloc (sizeof(float),glist->ndims);
	gr->VMean = calloc (sizeof(float),glist->ndims);
	
	ret=fread(gr->Center,sizeof(float),glist->ndims,f);
	ret=fread(gr->Spin,sizeof(float),glist->ndims,f);
	ret=fread(gr->VMean,sizeof(float),glist->ndims,f);
	ret=fread(gr->sig_evec,sizeof(float),9,f);
	ret=fread(gr->sig_eval,sizeof(float),3,f);
	ret=fread(&i,sizeof(int),1,f);

	gr->index = calloc(sizeof(int),gr->N);
	
	ret=fread(&i,sizeof(int),1,f);
	ret=fread(gr->index,sizeof(int),gr->N,f);
	ret=fread(&i,sizeof(int),1,f);
    }
    fclose(f);
    printf ("done.\n");
    
    return glist;
}

int *Groups2Id(group_list *list)
{
    int *arr;
    int i,j;
    group_type *gr;

    arr=(int *)malloc(sizeof(int)*list->NPart);
    for (i=0;i<list->NPart;i++)
	arr[i]=0;

    for (i=0;i<list->NGroups;i++)
    {
	gr = &list->group[i];
	for (j=0;j<gr->N;j++)
	    arr[gr->index[j]]=i+1;
    }

    return arr;
}

double *Groups2Id_d(group_list *list)
{
    double *arr;
    int i,j;
    group_type *gr;

    arr=(double *)malloc(sizeof(double)*list->NPart);
    for (i=0;i<list->NPart;i++)
	arr[i]=0;
    
    for (i=0;i<list->NGroups;i++)
    {
	gr = &list->group[i];
	for (j=0;j<gr->N;j++)
	    arr[gr->index[j]]=i+1;
    }

    return arr;
}

void freeFOF(group_list **list)
{
    int i;
    group_type *gr;

    if ((*list)->Pos!=NULL) free((*list)->Pos);
    if ((*list)->Vel!=NULL) free((*list)->Vel);
    for (i=0;i<(*list)->NGroups;i++)
    {
	gr = &(*list)->group[i];
	if (gr->index!=NULL) free(gr->index);
	if (gr->Center!=NULL) free(gr->Center);
	if (gr->index!=NULL) free(gr->Spin);
	if (gr->VMean!=NULL) free(gr->VMean);
    }

    free(*list);
    *list=NULL;
}
