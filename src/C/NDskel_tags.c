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
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "NDskel_tags.h"

#include "global.h"

//which = 0 -> node
//which = 1 -> seg
int getDataFieldID(NDskel *skl,int which, const char *name)
{
  int index=-1;
  long i;

  if (which == 0)
    {
      for (i=0;i<skl->nnodedata;i++)
	{
	  if (!strcmp(skl->nodedata_info[i],name)) index = i;
	}
      
    }
  else
    {
      for (i=0;i<skl->nsegdata;i++)
	{
	  if (!strcmp(skl->segdata_info[i],name)) index = i;
	}
     
    }
  
   return index;
}

//which = 0 -> node
//which = 1 -> seg
//val should have 1 value for every segment / node
//val[i] corresponds to skl->segpos[2*i] and skl->segpos[2*i+1] or skl->nodedata[i]
//returns the index of the new datafield
int AddDataField(NDskel *skl,int which,double *val,const char *name)
{
  int index=-1;
  long i,j;
  int old_ndata;

  if (which == 0)
    {
      old_ndata=skl->nnodedata;
      for (i=0;i<skl->nnodedata;i++)
	{
	  if (!strcmp(skl->nodedata_info[i],name)) index = i;
	}
      if (index == -1) index =skl->nnodedata++;

      if (old_ndata != skl->nnodedata)
	{
	  //new=1;
	  skl->nodedata = 
	    realloc (skl->nodedata,(long)skl->nnodedata*sizeof(double)*skl->nnodes);
	  skl->nodedata_info = 
	    realloc (skl->nodedata_info,(long)skl->nnodedata*sizeof(char *));

	  if (old_ndata!=0)
	    {
	      for (i=(skl->nnodes-1)*skl->nnodedata,j=(skl->nnodes-1)*old_ndata;
		   i>=0;i-=skl->nnodedata,j-=old_ndata)
		{
		  memmove(&skl->nodedata[i],&skl->nodedata[j],old_ndata*sizeof(double));
		}
	    }
	  if (index>=old_ndata) 
	    {
	      skl->nodedata_info[index]=malloc(NDSKEL_DATA_STR_SIZE*sizeof(char));
	      strcpy(skl->nodedata_info[index],name);
	    }
	}
      if (val!=NULL)
	for (i=0,j=0;i<skl->nnodes*skl->nnodedata;i+=skl->nnodedata,j++)
	  skl->nodedata[i+index] = val[j];
      else
	for (i=0,j=0;i<skl->nnodes*skl->nnodedata;i+=skl->nnodedata,j++)
	  skl->nodedata[i+index] = 0;

      for (i=0;i<skl->nnodes;i++)
	skl->Node[i].data = &skl->nodedata[skl->Node[i].index*skl->nnodedata];

    }
  else
    {
      old_ndata=skl->nsegdata;
      for (i=0;i<skl->nsegdata;i++)
	{
	  if (!strcmp(skl->segdata_info[i],name)) index = i;
	}
      if (index == -1) index =skl->nsegdata++;

      if (old_ndata != skl->nsegdata)
	{
	  //new=1;
	  skl->segdata = 
	    realloc (skl->segdata,(long)skl->nsegdata*sizeof(double)*skl->nsegs);
	  skl->segdata_info = 
	    realloc (skl->segdata_info,(long)skl->nsegdata*sizeof(char *));
	  if (old_ndata!=0)
	    {
	      for (i=(skl->nsegs-1)*skl->nsegdata,j=(skl->nsegs-1)*old_ndata;i>=0;
		   i-=skl->nsegdata,j-=old_ndata)
		{
		  memmove(&skl->segdata[i],&skl->segdata[j],old_ndata*sizeof(double));
		}
	    }
	  if (index>=old_ndata) 
	    {
	      skl->segdata_info[index]=malloc(NDSKEL_DATA_STR_SIZE*sizeof(char));
	      strcpy(skl->segdata_info[index],name);
	    }
	}

      if (val!=NULL)
	for (i=0,j=0;i<skl->nsegs*skl->nsegdata;i+=skl->nsegdata,j++)
	  skl->segdata[i+index] = val[j];
      else
	for (i=0,j=0;i<skl->nsegs*skl->nsegdata;i+=skl->nsegdata,j++)
	  skl->segdata[i+index] = 0;

      for (i=0;i<skl->nsegs;i++)
	skl->Seg[i].data = &skl->segdata[skl->Seg[i].index*skl->nsegdata];
    }

  return index;
}

int ComputeNodeTag(NDskel *skl,const char *name,double (*comp)(NDskel *,NDskl_node*))
{
	int index;
	int i;

	//val=calloc(skl->nnodes,sizeof(double));
	index=AddDataField(skl,0,NULL,name);
	//free(val);

	for (i=0;i<skl->nnodes;i++)
		skl->Node[i].data[index] = comp(skl,&(skl->Node[i]));

	return index;
}

int ComputeSegTag(NDskel *skl,const char *name,double (*comp)(NDskel *,NDskl_seg*))
{
	int index;
	int i;

	//val=calloc(skl->nsegs,sizeof(double));
	index=AddDataField(skl,1,NULL,name);
	//free(val);

	for (i=0;i<skl->nsegs;i++)
		skl->Seg[i].data[index] = comp(skl,&(skl->Seg[i]));

	return index;
}


int AddTag(NDskel *skl,NDfield *field,const char *name,int periodic)
{
  long i,j,k;
  FLOAT tmp;
  INT index;
  FLOAT pos[skl->ndims];
  float fpos[skl->ndims];
  double tag[(skl->nsegs>skl->nnodes)?skl->nsegs:skl->nnodes];
  char segname[255];
  DEFINE_NDVARS(d);

  if (field->datatype&(ND_CHAR|ND_UCHAR|ND_SHORT|ND_USHORT|ND_INT|ND_UINT|ND_LONG))
    {
      if (verbose>1) printf ("* Adding discrete '%s' data field to skeleton ...",name);fflush(0);
      /*
	for (i=0;i<skl->ndims;i++)
	{
	skl->x0[i] -= 0.5* skl->delta[i]/(skl->dims[i]-1);
	}
      */
      SETNDPOINTER(d,field->val,field->datatype);

      for (i=0,j=0;i<skl->nsegs*2;i+=2,j++)
	{

	  for (k=0;k<skl->ndims;k++) fpos[k] = skl->segpos[skl->ndims*i+k];
	  Coords2IndexfND(fpos, &index,field->x0,field->delta,field->dims, field->ndims,periodic);
	  SETNDFIELD_VAL(tag[j],d,index,field->datatype);
	}
      sprintf (segname,SEG_P1_TAG("%s"),name);
      AddDataField(skl,1,tag,segname);

      for (i=1,j=0;i<skl->nsegs*2;i+=2,j++)
	{
	  for (k=0;k<skl->ndims;k++) fpos[k] = skl->segpos[skl->ndims*i+k];
	  Coords2IndexfND(fpos, &index,field->x0,field->delta,field->dims, field->ndims,periodic);
	  SETNDFIELD_VAL(tag[j],d,index,field->datatype);
	}
      sprintf (segname,SEG_P2_TAG("%s"),name);
      AddDataField(skl,1,tag,segname);

      for (i=0,j=0;i<skl->nnodes;i++,j++)
	{
	  for (k=0;k<skl->ndims;k++) fpos[k] = skl->nodepos[skl->ndims*i+k];
	  Coords2IndexfND(fpos, &index,field->x0,field->delta,field->dims, field->ndims,periodic);
	  SETNDFIELD_VAL(tag[j],d,index,field->datatype);
	}
      AddDataField(skl,0,tag,name);
      /*
	for (i=0;i<skl->ndims;i++)
	{
	skl->x0[i] += 0.5* skl->delta[i]/(skl->dims[i]-1);
	}
      */
    }
  else
    {
      if (verbose>1) printf ("* Adding interpolated '%s' data field to skeleton ...",name);fflush(0);

      for (i=0,j=0;i<skl->nsegs*2;i+=2,j++)
	{
	  //for (k=0;k<skl->ndims;k++) pos[k] = (skl->segpos[skl->ndims*i+k]-skl->x0[k])*skl->dims[k]/skl->delta[k];
	  for (k=0;k<skl->ndims;k++) pos[k] = skl->segpos[skl->ndims*i+k];
	  InterpolateND(field,&tmp,pos,-1,periodic);tag[j]=tmp;

	  //if (j<10) printf ("%d: %f!=%f\n",j,tag[j],tmp);
	}
      sprintf (segname,SEG_P1_TAG("%s"),name);
      AddDataField(skl,1,tag,segname);

      for (i=1,j=0;i<skl->nsegs*2;i+=2,j++)
	{
	  //for (k=0;k<skl->ndims;k++) pos[k] = (skl->segpos[skl->ndims*i+k]-skl->x0[k])*skl->dims[k]/skl->delta[k];
	  for (k=0;k<skl->ndims;k++) pos[k] = skl->segpos[skl->ndims*i+k];
	  InterpolateND(field,&tmp,pos,-1,periodic);tag[j]=tmp;
	  //if (j<10) printf ("%d: %f!=%f\n",j,tag[j],tmp);
	}
      sprintf (segname,SEG_P2_TAG("%s"),name);
      AddDataField(skl,1,tag,segname);

      for (i=0,j=0;i<skl->nnodes;i++,j++)
	{
	  //for (k=0;k<skl->ndims;k++) pos[k] = (skl->nodepos[skl->ndims*i+k]-skl->x0[k])*skl->dims[k]/skl->delta[k];
	  for (k=0;k<skl->ndims;k++) pos[k] = skl->nodepos[skl->ndims*i+k];
	  InterpolateND(field,&tmp,pos,-1,periodic);tag[j]=tmp;
	}
      AddDataField(skl,0,tag,name);
    }

  if (verbose>1) printf (" done.\n");
  return 0;
}

int AddTagFromVal(NDskel *skl,NDfield *segs_field,NDfield *nodes_field,const char *name,int periodic)
{
  double *tag;
  char segname[255];

  if (segs_field!=NULL)
    {
      if ((segs_field->ndims==1)&&(segs_field->dims[0]==skl->nsegs))
	{
	  if (verbose>1) printf ("* Adding '%s' data field to skeleton segments ...",name);
	  Convert_NDfield(segs_field,ND_DOUBLE);

	  tag = (double *) segs_field->val;
	  sprintf (segname,"%s",name);
	  AddDataField(skl,1,tag,segname);
	  if (verbose>1) printf("done.\n");
	}
      else if ((segs_field->ndims==2)&&
	       (segs_field->dims[1]==2)&&
	       (segs_field->dims[0]==skl->nsegs))
	{
	  Convert_NDfield(segs_field,ND_DOUBLE);
		
	  tag = (double *) segs_field->val;
	  sprintf (segname,SEG_P1_TAG("%s"),name);
	  if (verbose>1) printf ("* Adding '%s' data field to skeleton segments ...",segname);
	  AddDataField(skl,1,tag,segname);
	  if (verbose>1) printf("done.\n");

	  tag = &(tag[skl->nsegs]);
	  sprintf (segname,SEG_P2_TAG("%s"),name);
	  if (verbose>1) printf ("* Adding '%s' data field to skeleton segments ...",segname);
	  AddDataField(skl,1,tag,segname);
	  if (verbose>1) printf("done.\n");
	}
      else
	{
	  fprintf (stderr,"The field used for segments tag does not have correct dimensions: [%d,%d] or [%d]\n",skl->nsegs,2,skl->nsegs);
	  return -1;
	}
    }

  if (nodes_field!=NULL)
    {
      if ((nodes_field->ndims!=1)&&(nodes_field->dims[0]==skl->nnodes))
	{
	  if (verbose>1) printf ("* Adding '%s' data field to skeleton nodes ...",name);
	  Convert_NDfield(nodes_field,ND_DOUBLE);
	  tag = (double *) nodes_field->val;	  
	  AddDataField(skl,0,tag,name);
	  if (verbose>1) printf("done.\n");
	}
      else
	{
	  fprintf (stderr,"The field used for nodes tag does not have correct dimensions: [%d]\n",skl->nnodes);
	  return -1;
	}
    }

  //printf (" done.\n");
  return 0;
}

int addDataFieldFromNDfield(NDskel *skl,NDfield *f, const char *name, int periodic)
{
  int done=1;
  if (f->ndims==1)
    {
      if (f->dims[0] == skl->nnodes)
	AddTagFromVal(skl,NULL,f,name,periodic);
      else if (f->dims[0] == skl->nsegs)
	AddTagFromVal(skl,f,NULL,name,periodic);
      else done=0;
    }
  else if ((f->ndims==2)&&(f->dims[1]==2)&&(f->dims[0]==skl->nsegs))
    {
      AddTagFromVal(skl,f,NULL,name,periodic);	  
    }
  else if (f->ndims==skl->ndims)
    {
      AddTag(skl,f,name,periodic);
    }
  else done=0;
  
  if (!done)
    {
      fprintf(stderr,"ERROR adding data field from NDfield.\n");
      fprintf(stderr,"   input field has wrong dimensions.\n");
      fprintf(stderr,"   should be a grid, with dims [%ld] [%ld] [%ld,%ld] or ndims=%ld.\n",
	      (long)skl->nnodes,(long)skl->nsegs,(long)skl->nsegs,(long)2,(long)skl->ndims);
      return -1;
    }

  return 0;
}

