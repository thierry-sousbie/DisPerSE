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
#include <string.h>
#include <math.h>
#include <assert.h>

#include "NDskel_VTK_IO.h"
#include "NDskeleton.h"
#include "NDskel_tags.h"
#include "global.h"
#include "myendian.h"

int IsVTKSkeleton(const char *fname)
{
  return 0;
}

NDskel *Load_NDskeletonFromVTK(const char *fname)
{
  NDskel *skl=NULL;

  return skl;
}

int Save_NDskeletonToVTK(NDskel *skl,const char *fname,int type)
{
  int xml=type&NDSKEL_VTK_XML;
  int binary=type&NDSKEL_VTK_BIN;
  long i,j,k;
  FILE *f;
  char format[255];
  long offset=0;
  unsigned int size[50];
  int sizeindex=0;
  int segDataIndex[skl->nnodedata][2];
  //int segDataIsShared[skl->nsegdata];
  
  NDskl_seg **filTab=NULL;  
  int *filSize=NULL;
  double **filData=NULL;
  char **filDataInfo=NULL;
  int nfil=getNDskelFilTab(skl,&filTab,&filSize);
  int nfilData=getNDskelFilTabInfo(skl,filTab,nfil,&filDataInfo,&filData);
  long npoints= skl->nnodes; 
  long ncells=(long)nfil+skl->nnodes;

  for (i=0;i<nfil;i++) npoints+=filSize[i]-1;

  //for (i=0;i<skl->nsegdata;i++) segDataIsShared[i]=0;
  for (i=0;i<skl->nnodedata;i++)
    {
      char tmp[255];
      sprintf(tmp,"%s_p1",skl->nodedata_info[i]);
      segDataIndex[i][0]=getDataFieldID(skl,1,tmp);
      sprintf(tmp,"%s_p2",skl->nodedata_info[i]);
      segDataIndex[i][1]=getDataFieldID(skl,1,tmp);
      if ((segDataIndex[i][0]<0)&&(segDataIndex[i][1]<0))
	{
	  sprintf(tmp,"%s",skl->nodedata_info[i]);
	  segDataIndex[i][0]=getDataFieldID(skl,1,tmp);
	  segDataIndex[i][1]=segDataIndex[i][0];
	}
      if ((segDataIndex[i][0]>=0)&&(segDataIndex[i][1]>=0))
	{
	  //segDataIsShared[segDataIndex[i][0]]=1;
	  //segDataIsShared[segDataIndex[i][1]]=1;
	}
    }
 
  if (xml)
    sprintf(format,"%s",(binary)?"appended":"ascii");
  else
    sprintf(format,"%s",(binary)?"BINARY":"ASCII");

  f=fopen(fname,"w");

  if (xml)
    {
      fprintf(f,"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
      fprintf(f,"  <PolyData>\n");
      fprintf(f,"    <Piece NumberOfPoints=\"%ld\" NumberOfVerts=\"%ld\" NumberOfLines=\"%ld\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n",(long)npoints,(long)skl->nnodes,(long)nfil);
      
    }
  else
    {
      fprintf(f,"# vtk DataFile Version 2.0\n");
      fprintf(f,"%s\n",skl->comment);
      fprintf(f,"%s\n",format);
      fprintf(f,"DATASET POLYDATA\n");     
    }

  if (xml) 
    {
      fprintf(f,"<Points>\n");
      if (binary)
	fprintf(f,"<DataArray NumberOfComponents=\"%ld\" Name=\"%s\" type=\"Float32\" format=\"%s\" offset=\"%ld\" />\n",
	      (long)3,"coords",format,offset);  
      else
	fprintf(f,"<DataArray NumberOfComponents=\"%ld\" Name=\"%s\" type=\"Float32\" format=\"%s\" >\n",
		(long)3,"coords",format);  
      
      if (binary)
	{
	  size[sizeindex]=sizeof(float)*3*npoints;
	  offset+=size[sizeindex++]+sizeof(int);
	}
      
    }
  else fprintf(f,"POINTS %ld float\n",(long)npoints);
  
  float v[3];
  for (i=0;i<3;i++) v[i]=0;
  for (i=0;i<skl->nnodes;i++)
    {
      memcpy(v,skl->Node[i].pos,skl->ndims*sizeof(float));
      if (binary) 
	{
	  if (!xml) 
	    fwriteBE(v,sizeof(float),3,f);
	}
      else 
	fprintf(f,"%g %g %g\n",v[0],v[1],v[2]);
    }
  for (i=0;i<nfil;i++)
    {
      NDskl_seg *seg=filTab[i];
      if (seg->Next!=NULL)
	{
	  do {
	    seg=seg->Next;
	    memcpy(v,seg->pos,skl->ndims*sizeof(float));
	    if (binary)
	      {
		if (!xml) fwriteBE(v,sizeof(float),3,f);
	      }
	    else 
	      fprintf(f,"%g %g %g\n",v[0],v[1],v[2]);
	  } while (seg->Next!=NULL);
	}
    }
  if (xml) {
    if (binary)
      fprintf(f,"</Points>\n");
    else
      fprintf(f,"</DataArray>\n</Points>\n");
  }
  else if (binary) fprintf(f,"\n");

  if (xml)
    {
      
      fprintf(f,"<Verts>\n");
      if (binary)
	fprintf(f,"<DataArray type=\"Int%ld\" Name=\"%s\" format=\"%s\" offset=\"%ld\"/>\n",(long)32,"connectivity",format,offset);
      else
	fprintf(f,"<DataArray type=\"Int%ld\" Name=\"%s\" format=\"%s\">\n",(long)32,"connectivity",format);
      
      if (binary) 
	{
	  size[sizeindex]=sizeof(int)*skl->nnodes;
	  offset+=size[sizeindex++]+sizeof(int);
	}
      else
	{	  
	  for (i=0;i<skl->nnodes;i++)
	    fprintf(f," %ld\n",i);
	}
      if (!binary) fprintf(f,"</DataArray>\n");

      if (binary)
	fprintf(f,"<DataArray type=\"Int%ld\" Name=\"%s\" format=\"%s\" offset=\"%ld\"/>\n",(long)32,"offsets",format,offset);
      else
	fprintf(f,"<DataArray type=\"Int%ld\" Name=\"%s\" format=\"%s\">\n",(long)32,"offsets",format);
      
      if (binary) 
	{
	  size[sizeindex]=sizeof(int)*skl->nnodes;
	  offset+=size[sizeindex++]+sizeof(int);
	}
      else
	{	  
	  for (i=0;i<skl->nnodes;i++)
	    fprintf(f," %ld\n",i+1);
	}
      if (!binary) fprintf(f,"</DataArray>\n");
      fprintf(f,"</Verts>\n");
      
      fprintf(f,"<Lines>\n");
      if (binary)
	fprintf(f,"<DataArray type=\"Int%ld\" Name=\"%s\" format=\"%s\" offset=\"%ld\"/>\n",(long)32,"connectivity",format,offset);
      else
	fprintf(f,"<DataArray type=\"Int%ld\" Name=\"%s\" format=\"%s\">\n",(long)32,"connectivity",format);
      
      if (binary) 
	{
	  for (i=0,k=0;i<nfil;i++) k+=filSize[i]+1;
	  size[sizeindex]=sizeof(int)*k;
	  offset+=size[sizeindex++]+sizeof(int);
	}
      else
	{	  
	  for (i=0,k=skl->nnodes;i<nfil;i++) 
	    {
	      fprintf(f," %ld",(long)filTab[i]->nodes[0]);
	      for (j=0;j<filSize[i]-1;j++,k++)
		fprintf(f," %ld",k);
	      fprintf(f," %ld\n",(long)filTab[i]->nodes[1]);
	    }
	}
      if (!binary) fprintf(f,"</DataArray>\n");
      if (binary)
	fprintf(f,"<DataArray type=\"Int%ld\" Name=\"%s\" format=\"%s\" offset=\"%ld\"/>\n",(long)32,"offsets",format,offset);
      else
	fprintf(f,"<DataArray type=\"Int%ld\" Name=\"%s\" format=\"%s\">\n",(long)32,"offsets",format);
      
      if (binary) 
	{
	  size[sizeindex]=sizeof(int)*nfil;
	  offset+=size[sizeindex++]+sizeof(int);
	}
      else
	{
	  j=0;
	  for (i=0;i<nfil;i++)
	    {
	      j+=filSize[i]+1;
	      fprintf(f,"%ld\n",j);		
	    }
	}	  

      if (!binary) fprintf(f,"</DataArray>\n");
      fprintf(f,"</Lines>\n");
    }
  else 
    {
      if (!binary)
	{	  
	  fprintf(f,"VERTICES %ld %ld\n",(long)skl->nnodes,(long)2*skl->nnodes);
	  for (i=0;i<skl->nnodes;i++)
	    fprintf(f,"1 %ld\n",i);
	  for (i=0,k=0;i<nfil;i++) k+=filSize[i]+2;
	  fprintf(f,"LINES %ld %ld\n",(long)nfil,k);
	  for (i=0,k=skl->nnodes;i<nfil;i++) 
	    {
	      fprintf(f,"%ld",(long)filSize[i]+1);
	      fprintf(f," %ld",(long)filTab[i]->nodes[0]);
	      for (j=0;j<filSize[i]-1;j++,k++)
		fprintf(f," %ld",k);
	      fprintf(f," %ld\n",(long)filTab[i]->nodes[1]);
	    }
	  /*
	  fprintf(f,"CELL_TYPES %ld\n",ncells);
	  for (i=0;i<skl->nnodes;i++)
	    fprintf(f,"%d\n",(int)1); //vertex
	  for (i=0;i<nfil;i++) 
	    fprintf(f,"%d\n",(int)4); //polyline
	  */
	}
      else
	{
	  int tmpui;
	  
	  fprintf(f,"VERTICES %ld %ld\n",(long)skl->nnodes,(long)2*skl->nnodes);
	 
	  for (i=0;i<skl->nnodes;i++)
	    {
	      tmpui=1;fwriteBE(&tmpui,sizeof(unsigned int),1,f);
	      tmpui=i;fwriteBE(&tmpui,sizeof(unsigned int),1,f);
	    }
	  for (i=0,k=0;i<nfil;i++) k+=filSize[i]+2;
	  fprintf(f,"\nLINES %ld %ld\n",(long)nfil,k);
	  for (i=0,k=skl->nnodes;i<nfil;i++) 
	    {
	      tmpui=filSize[i]+1;fwriteBE(&tmpui,sizeof(unsigned int),1,f);
	      tmpui=filTab[i]->nodes[0];fwriteBE(&tmpui,sizeof(unsigned int),1,f);
	      for (j=0;j<filSize[i]-1;j++,k++)
		{
		  tmpui=k;
		  fwriteBE(&tmpui,sizeof(unsigned int),1,f);
		}
	      tmpui=filTab[i]->nodes[1];fwriteBE(&tmpui,sizeof(unsigned int),1,f);
	    }
	  fprintf(f,"\n");
	}
    }
   
  // POINT DATA
  int found=0;
  for (i=0;i<skl->nnodedata;i++)
    {
      double *d=skl->nodedata;
      char name[255];
      strcpy(name,skl->nodedata_info[i]);
      
      if (!found)
	{
	  int nfields=skl->nnodedata;
	  if (xml)
	    {
	      fprintf(f,"<PointData>\n");
	    }
	  else
	    {
	      fprintf(f,"POINT_DATA %ld\n",(long)npoints);
	      fprintf(f,"FIELD point_field %ld\n",(long)nfields+1);
	    }
	  found=1;
	}

      if (xml)
	{
	  if (binary)
	    {		  
	      fprintf(f,"<DataArray Name=\"%s\" type=\"Float64\" format=\"%s\" offset=\"%ld\"/>\n",name,format,offset);
	      size[sizeindex]=sizeof(double)*npoints;
	      offset+=size[sizeindex++]+sizeof(int);
	    }
	  else
	    {
	      fprintf(f,"<DataArray Name=\"%s\" type=\"Float64\" format=\"%s\" >\n",name,format);
	    }
	}
      else fprintf(f,"%s 1 %ld double\n",name,npoints);

      if ((!binary)||(!xml))
	{
	  for (j=0;j<skl->nnodes;j++) 
	    {
	      if (binary) fwriteBE(&d[j*skl->nnodedata+i],sizeof(double),1,f);
	      else fprintf(f,"%g\n",d[j*skl->nnodedata+i]);	
	    }
	  for (j=0;j<nfil;j++) 
	    {
	      NDskl_seg *seg=filTab[j];
	      if (seg->Next!=NULL)
		{
		  do {
		    seg=seg->Next;
		    double val=-1;
		    if (segDataIndex[i][0]>=0) val=seg->data[segDataIndex[i][0]];
		    if (binary) fwriteBE(&val,sizeof(double),1,f);
		    else fprintf(f," %g\n",val);
		  } while (seg->Next!=NULL);
		}
	    }
	  if (xml) fprintf(f,"</DataArray>\n");
	  else if (binary) fprintf(f,"\n");
	}
    }

  if (!found)
    {
      if (xml)
	{
	  fprintf(f,"<PointData>\n");
	}
      else
	{
	  fprintf(f,"POINT_DATA %ld\n",(long)npoints);
	  fprintf(f,"FIELD point_field %ld\n",(long)1);
	}
    }
  if (xml)
    {
      if (binary)
	{		  
	  fprintf(f,"<DataArray Name=\"%s\" type=\"Float64\" format=\"%s\" offset=\"%ld\"/>\n",CRITICAL_INDEX_TAG,format,offset);
	  size[sizeindex]=sizeof(double)*npoints;
	  offset+=size[sizeindex++]+sizeof(int);
	}
      else
	{
	  fprintf(f,"<DataArray Name=\"%s\" type=\"Float64\" format=\"%s\" >\n",CRITICAL_INDEX_TAG,format);
	}
    }
  else fprintf(f,"%s 1 %ld double\n",CRITICAL_INDEX_TAG,npoints);

  if ((!xml)||(!binary))
    {
      for (i=0;i<npoints;i++)
	{
	  double val=-1;
	  if (i<skl->nnodes) val=skl->Node[i].type;
	  if (binary) fwriteBE(&val,sizeof(double),1,f);
	  else fprintf(f," %g\n",val);
	}    
    }

  if (xml) {
    if (!binary) fprintf(f,"</DataArray>\n");
    fprintf(f,"</PointData>\n");
  }
  
  filDataInfo=(char **)realloc(filDataInfo,(nfilData+1)*sizeof(char*));
  filDataInfo[nfilData]=(char*)malloc(sizeof(char)*50);
  strcpy(filDataInfo[nfilData],"flags");
  //filData=(double**)realloc(filData,(nfilData+1)*sizeof(double*));
  //filData[nfilData]=(double*)calloc(ncells,sizeof(double)); 
  //for (i=0;i<skl->nnodes;i++) filData[nfilData][i]=(double)skl->Node[i].flags;
  for (j=0;j<nfil;j++)
    {
      filData[j]=(double*)realloc(filData[j],(nfilData+1)*sizeof(double));
      int f1=skl->Node[filTab[j]->nodes[0]].flags;
      int f2=skl->Node[filTab[j]->nodes[1]].flags;
      filData[j][nfilData]=(double)(f1|f2);
      }
  nfilData++;
  
  // CELL DATA
  found=0;
  for (i=0;i<nfilData;i++)
    {
      char name[255];
      strcpy(name,filDataInfo[i]);
      
      if (!found)
	{
	  if (xml)
	    {
	      fprintf(f,"<CellData>\n");
	    }
	  else
	    {
	      fprintf(f,"CELL_DATA %ld\n",ncells);
	      fprintf(f,"FIELD cell_field %ld\n",(long)nfilData);
	    }
	  found=1;
	}
      

      if (xml)
	{
	  if (binary)
	    {		  
	      fprintf(f,"<DataArray Name=\"%s\" type=\"Float64\" format=\"%s\" offset=\"%ld\"/>\n",name,format,offset);
	      size[sizeindex]=sizeof(double)*ncells;
	      offset+=size[sizeindex++]+sizeof(int);
	    }
	  else
	    {
	      fprintf(f,"<DataArray Name=\"%s\" type=\"Float64\" format=\"%s\" >\n",name,format);
	    }
	}
      else fprintf(f,"%s 1 %ld double\n",name,ncells);

      if ((!binary)||(!xml))
	{
	  for (j=0;j<skl->nnodes;j++) 
	    {
	      double tmp=-1;
	      if (i==nfilData-1) tmp=skl->Node[j].flags;
	      if (binary) fwriteBE(&tmp,sizeof(double),1,f);
	      else fprintf(f,"%g\n",tmp);	
	    }
	  for (j=0;j<nfil;j++) 
	    {
	      if (binary) fwriteBE(&filData[j][i],sizeof(double),1,f);
	      else fprintf(f,"%g\n",filData[j][i]);	
	    }
	  
	  if (xml) fprintf(f,"</DataArray>\n");
	  else if (binary) fprintf(f,"\n");
	}
    }

   
  if (xml) {
    if (!binary) fprintf(f,"</DataArray>\n");
    fprintf(f,"</CellData>\n");
  }

    
  if (xml)
    {
      fprintf(f,"    </Piece>\n");
      fprintf(f,"  </PolyData>\n");

      if (binary)
	{
	  int sizeID=0;
	  fprintf(f,"<AppendedData encoding=\"raw\">\n_");
	  fwrite(&size[sizeID],sizeof(size[sizeID]),1,f);sizeID++;
	  float v[3];
	  for (i=0;i<3;i++) v[i]=0;
	  for (i=0;i<skl->nnodes;i++)
	    {
	      memcpy(v,skl->Node[i].pos,skl->ndims*sizeof(float));
	      fwrite(v,sizeof(float),3,f);
	    }

	  for (i=0;i<nfil;i++)
	    {
	      NDskl_seg *seg=filTab[i];
	      if (seg->Next!=NULL)
		{
		  do {
		    seg=seg->Next;
		    memcpy(v,seg->pos,skl->ndims*sizeof(float));
		    fwrite(v,sizeof(float),3,f);
		    } while (seg->Next!=NULL);
		}
	    }

	  int ct=0;
	  
	  //Verts connectivity	  
	  fwrite(&size[sizeID],sizeof(size[sizeID]),1,f);sizeID++;
	  for (ct=0;ct<skl->nnodes;ct++) 
	    fwrite(&ct,sizeof(int),1,f);
	  //Verts offset
	  fwrite(&size[sizeID],sizeof(size[sizeID]),1,f);sizeID++;
	  for (ct=1;ct<=skl->nnodes;ct++) 
	    fwrite(&ct,sizeof(int),1,f);
	  

	  //Lines connectivity	  
	  fwrite(&size[sizeID],sizeof(size[sizeID]),1,f);sizeID++;
	  for (i=0,ct=skl->nnodes;i<nfil;i++) 
	    {
	      fwrite(&filTab[i]->nodes[0],sizeof(int),1,f);
	      for (j=0;j<filSize[i]-1;j++,ct++)
		fwrite(&ct,sizeof(int),1,f);
	      fwrite(&filTab[i]->nodes[1],sizeof(int),1,f);
	    }
	  //Lines offset
	  fwrite(&size[sizeID],sizeof(size[sizeID]),1,f);sizeID++;
	  ct=0;
	  for (i=0;i<nfil;i++)
	    {
	      ct+=filSize[i]+1;
	      fwrite(&ct,sizeof(int),1,f);
	    }
	  
	  //points data
	  for (i=0;i<skl->nnodedata;i++)
	    {
	      fwrite(&size[sizeID],sizeof(size[sizeID]),1,f);sizeID++;
	      for (j=0;j<skl->nnodes;j++) 
		fwrite(&skl->nodedata[j*skl->nnodedata+i],sizeof(double),1,f);	
	      for (j=0;j<nfil;j++) 
		{
		  NDskl_seg *seg=filTab[j];
		  if (seg->Next!=NULL)
		    {
		      do {
			seg=seg->Next;
			double val=-1;
			if (segDataIndex[i][0]>=0) val=seg->data[segDataIndex[i][0]];
			fwrite(&val,sizeof(double),1,f);
		      } while (seg->Next!=NULL);
		    }
		}
	    }
	  fwrite(&size[sizeID],sizeof(size[sizeID]),1,f);sizeID++;
	  for (i=0;i<npoints;i++)
	    {
	      double val=-1;
	      if (i<skl->nnodes) val=skl->Node[i].type;
	      fwrite(&val,sizeof(double),1,f);
	    }    

	  //cells data
	  for (i=0;i<nfilData;i++)
	    {
	      double tmp=-1;	      
	      fwrite(&size[sizeID],sizeof(size[sizeID]),1,f);sizeID++;
	      for (j=0;j<skl->nnodes;j++) 
		{
		  if (i==nfilData-1) tmp=skl->Node[j].flags;
		  fwrite(&tmp,sizeof(double),1,f);	
		}
	      for (j=0;j<nfil;j++) 
		fwrite(&filData[j][i],sizeof(double),1,f);	
	    }
	  	  
	  fprintf(f,"</AppendedData>\n");
	}     
    }
  
  if (xml) fprintf(f,"</VTKFile>\n");
  fclose(f);
  
  freeNDskelFilTab(&filTab,&filSize);
  freeNDskelFilTabInfo(&filDataInfo,&filData,nfilData,nfil);
  
  return 0;
}
