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

#include "NDnet_VTK_IO.h"
#include "NDnetwork.h"
#include "NDnetwork_tags.h"
#include "global.h"
#include "myendian.h"

int IsVTKNetwork(const char *fname)
{
  return 0;
}

NDnetwork *Load_NDnetworkFromVTK(const char *fname)
{
  NDnetwork *net=NULL;

  return net;
}

int Save_NDnetworkToVTK(NDnetwork *net,const char *fname,int type)
{
  int xml=type&NDNET_VTK_XML;
  int binary=type&NDNET_VTK_BIN;
  long i,j;
  long ncells=0;
  long cell_type=0;
  FILE *f;
  char format[255];
  
  //int value_tag = NDDataIndex(net,0,VALUE_TAG);
  double *value=NULL;
  double vmin=0; 
  double dv=1;
  long offset=0;

  unsigned char vtk_cellT=0;
  NDNET_VTK_HINT size[256];
  
  double *dataPtr[256];
  long dataSize[256];
  long nData=0;
  
  int sizeindex=0;
  /*
  if (value_tag>=0) 
    {
      double vmax;
      value=net->data[value_tag].data;
      vmin=(value[0]);
      vmax=(value[0]);
      for (i=0;i<net->nvertex;i++)
	{
	  if ((value[i])<vmin) vmin=value[i];
	  if ((value[i])>vmax) vmax=value[i];
	}
      for (i=0;i<net->ndims;i++) dv+=net->delta[i];
      dv/=(net->ndims*(vmax-vmin));
    }
  */
  if (xml)
    sprintf(format,"%s",(binary)?"appended":"ascii");
  else
    sprintf(format,"%s",(binary)?"BINARY":"ASCII");

  for (i=1;i<net->ndims_net+1;i++) if (net->haveVertexFromFace[i]) cell_type=i;
  if (cell_type) ncells=net->nfaces[cell_type];

  
  if (cell_type==1) vtk_cellT=3;
  else if (cell_type==2) vtk_cellT=5;
  else if (cell_type==3) vtk_cellT=10;

  f=fopen(fname,"w");

  if (xml)
    {
      if (sizeof(NDNET_VTK_HINT)==4)
	fprintf(f,"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
      else
	fprintf(f,"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");

      fprintf(f,"  <UnstructuredGrid>\n");
      fprintf(f,"    <Piece NumberOfPoints=\"%ld\" NumberOfCells=\"%ld\">\n",(long)net->nvertex,(long)ncells);      
    }
  else
    {
      fprintf(f,"# vtk DataFile Version 2.0\n");
      fprintf(f,"%s\n",net->comment);
      fprintf(f,"%s\n",format);
      fprintf(f,"DATASET UNSTRUCTURED_GRID\n");     
    }

  if (xml) 
    {
      fprintf(f,"<Points>\n");//encoding=\"raw\"
      if (binary)
	fprintf(f,"<DataArray NumberOfComponents=\"%ld\" Name=\"%s\" type=\"Float32\" format=\"%s\" offset=\"%ld\" />\n",
	      (long)3,"coords",format,offset);  
      else
	fprintf(f,"<DataArray NumberOfComponents=\"%ld\" Name=\"%s\" type=\"Float32\" format=\"%s\" >\n",
		(long)3,"coords",format);  
      
      if (binary)
	{
	  
	  //fprintf(f,"_");
	  size[sizeindex]=sizeof(float)*3*net->nvertex;
	  offset+=size[sizeindex++]+sizeof(NDNET_VTK_HINT);
	  //fwrite(&size,sizeof(size),1,f);
	}
      
    }
  else fprintf(f,"POINTS %ld float\n",(long)net->nvertex);
  
  float v[3];
  long dimStride = net->ndims - net->ndims_net;
 
  for (i=0;i<3;i++) v[i]=0;
  for (i=0,j=0;i<net->nvertex;i++)
    {
      v[0]=net->v_coord[j++];

      if (net->ndims_net>1) 
	{
	  v[1]=net->v_coord[j++];
	  if (net->ndims_net>2) v[2]=net->v_coord[j++];
	  else if (value!=NULL) v[2]=((value[i])-vmin)*dv;
	}
      else if (value!=NULL) v[1]=((value[i])-vmin)*dv;

      j+=dimStride;

      if (binary) 
	{
	  if (!xml) fwriteBE(v,sizeof(float),3,f);
	}
      else 
	fprintf(f,"%g %g %g\n",v[0],v[1],v[2]);
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
      NDNET_IDCUMT ct;

      fprintf(f,"<Cells>\n");

      if (binary)
	fprintf(f,"<DataArray type=\"Int%ld\" Name=\"%s\" format=\"%s\" offset=\"%ld\"/>\n",sizeof(NDNET_UINT)*8,"connectivity",format,offset);
      else
	fprintf(f,"<DataArray type=\"Int%ld\" Name=\"%s\" format=\"%s\">\n",sizeof(NDNET_UINT)*8,"connectivity",format);


      if (binary) 
	{
	  //fprintf(f,"_");
	  //size=sizeof(NDNET_UINT)*(cell_type+1)*ncells;
	  size[sizeindex]=sizeof(NDNET_UINT)*(cell_type+1)*ncells;
	  offset+=size[sizeindex++]+sizeof(NDNET_VTK_HINT);	 
	  //fwrite(&size,sizeof(size),1,f);
	  //fwrite(net->f_vertexIndex[cell_type],sizeof(NDNET_UINT),(cell_type+1)*ncells,f);
	}
      else
	for (i=0;i<(cell_type+1)*ncells;i++) 
	  {
	    fprintf(f," %ld",(long)net->f_vertexIndex[cell_type][i]);
	    if ((i+1)%(cell_type+1) == 0) fprintf(f,"\n");
	  }
      //if (binary) fprintf(f,"\n");
      if (!binary) fprintf(f,"</DataArray>\n");
       if (binary)
	 fprintf(f,"<DataArray type=\"Int%ld\" Name=\"%s\" format=\"%s\" offset=\"%ld\"/>\n",sizeof(NDNET_IDCUMT)*8,"offsets",format,offset);
       else
	 fprintf(f,"<DataArray type=\"Int%ld\" Name=\"%s\" format=\"%s\">\n",sizeof(NDNET_IDCUMT)*8,"offsets",format);

      if (binary) 
	{
	  size[sizeindex]=sizeof(NDNET_IDCUMT)*ncells;
	  offset+=size[sizeindex++]+sizeof(NDNET_VTK_HINT);
	}
      else
	for (ct=(cell_type+1);ct<(cell_type+1)*(ncells+1);ct+=(cell_type+1)) 
	  fprintf(f,"%ld\n",(long)ct);
      //if (binary) fprintf(f,"\n");
      if (!binary) fprintf(f,"</DataArray>\n");
      if (binary)
	fprintf(f,"<DataArray type=\"UInt8\" Name=\"%s\" format=\"%s\" offset=\"%ld\"/>\n","types",format,offset);
      else
	fprintf(f,"<DataArray type=\"UInt8\" Name=\"%s\" format=\"%s\">\n","types",format);
      
      if (binary) 
	{
	  size[sizeindex]=sizeof(unsigned char)*ncells;
	  offset+=size[sizeindex++]+sizeof(NDNET_VTK_HINT);
	}
      else
	for (i=0;i<ncells;i++) fprintf(f,"%u\n",(unsigned int)vtk_cellT);
      //if (binary) fprintf(f,"\n");
      if (!binary) fprintf(f,"</DataArray>\n");
    }
  else 
    {
      fprintf(f,"CELLS %ld %ld\n",ncells, ncells*(cell_type+1 +1));
      unsigned int nid=cell_type+1;	  
      if (binary)
	{
	  int val[nid];
	  for (i=0;i<(cell_type+1)*ncells;i+=(cell_type+1))
	    {
	      for (j=0;j<nid;j++) val[j]=net->f_vertexIndex[cell_type][i+j];
	      fwriteBE(&nid,sizeof(unsigned int),1,f);
	      fwriteBE(val,sizeof(int),nid,f);
	    }
	  fprintf(f,"\n");
	  fprintf(f,"CELL_TYPES %ld\n",ncells);
	  int cellT=vtk_cellT;
	  for (i=0;i<ncells;i++) fwriteBE(&cellT,sizeof(int),1,f);
	}
      else
	{
	  for (i=0;i<(cell_type+1)*ncells;i+=(cell_type+1))
	    {
	      fprintf(f,"%d",nid);
	      for (j=0;j<(cell_type+1);j++) 
		fprintf(f," %ld",(long)net->f_vertexIndex[cell_type][i+j]);
	      fprintf(f,"\n");
	    }
	  fprintf(f,"CELL_TYPES %ld\n",ncells);
	  for (i=0;i<ncells;i++) fprintf(f,"%d\n",(int)vtk_cellT);
	}
    }
  //if (binary) fprintf(f,"\n");
  if (xml) fprintf(f,"</Cells>\n");
  //else fprintf(f,"\n");

  // POINT DATA
  int found=0;
  for (i=0;i<net->ndata;i++)
    {
      if (net->data[i].type!=0) continue;
      double *d=net->data[i].data;
      char name[255];
      strcpy(name,net->data[i].name);

      if (!found)
	{
	  int nfields=0;
	  for (j=0;j<net->ndata;j++) if (net->data[i].type==0) nfields++;
	  if (xml)
	    {
	      fprintf(f,"<PointData>\n");
	    }
	  else
	    {
	      fprintf(f,"POINT_DATA %ld\n",(long)net->nvertex);
	      fprintf(f,"FIELD point_field %ld\n",(long)nfields);
	    }
	  found=1;
	}

      if (xml)
	{
	  if (binary)
	    {		  
	      fprintf(f,"<DataArray Name=\"%s\" type=\"Float64\" format=\"%s\" offset=\"%ld\"/>\n",name,format,offset);
	      size[sizeindex]=sizeof(double)*net->nvertex;
	      offset+=size[sizeindex++]+sizeof(NDNET_VTK_HINT);
	      dataSize[nData]=net->nvertex;
	      dataPtr[nData++]=(void*)d;		  
	    }
	  else
	    {
	      fprintf(f,"<DataArray Name=\"%s\" type=\"Float64\" format=\"%s\" >\n",name,format);
		  
	      for (j=0;j<net->nvertex;j++) fprintf(f,"%g\n",d[j]);	
	      fprintf(f,"</DataArray>\n");
	    }
	}
      else
	{
	  for (j=0;j<strlen(name);j++) if (name[j]==' ') name[j]='_';
	  //fprintf(f,"SCALARS %s float64\n",name);
	  //fprintf(f,"LOOKUP_TABLE default\n");
	  fprintf(f,"%s 1 %ld double\n",name,(long)net->nvertex);
	  if (binary)
	    {	   
	      fwriteBE(d,sizeof(double),net->nvertex,f);
	      fprintf(f,"\n");
	    }
	  else
	    {
	      for (j=0;j<net->nvertex;j++) fprintf(f,"%g\n",d[j]);
	    }	      
	}
    }
  if ((found)&&(xml)) fprintf(f,"</PointData>\n");

  // CELL DATA
  found=0;
  for (i=0;i<net->ndata;i++)
    {
      if (net->data[i].type!=cell_type) continue;
      double *d=net->data[i].data;
      char name[255];
      strcpy(name,net->data[i].name);

      if (!found)
	{
	  int nfields=0;
	  for (j=0;j<net->ndata;j++) if (net->data[i].type==0) nfields++;
	  if (xml)
	    {
	      fprintf(f,"<CellData>\n");
	    }
	  else
	    {
	      fprintf(f,"CELL_DATA %ld\n",(long)net->nfaces[cell_type]);
	      fprintf(f,"FIELD cell_field %ld\n",(long)nfields);
	    }
	  found=1;
	}

      if (xml)
	{
	  if (binary)
	    {		  
	      fprintf(f,"<DataArray Name=\"%s\" type=\"Float64\" format=\"%s\" offset=\"%ld\"/>\n",name,format,offset);
	      size[sizeindex]=sizeof(double)*net->nfaces[cell_type];
	      offset+=size[sizeindex++]+sizeof(NDNET_VTK_HINT);
	      dataSize[nData]=net->nfaces[cell_type];
	      dataPtr[nData++]=(void*)d;		  
	    }
	  else
	    {
	      fprintf(f,"<DataArray Name=\"%s\" type=\"Float64\" format=\"%s\" >\n",name,format);
		  
	      for (j=0;j<net->nfaces[cell_type];j++) fprintf(f,"%g\n",d[j]);	
	      fprintf(f,"</DataArray>\n");
	    }
	}
      else
	{
	  for (j=0;j<strlen(name);j++) if (name[j]==' ') name[j]='_';
	  //fprintf(f,"SCALARS %s float64\n",name);
	  //fprintf(f,"LOOKUP_TABLE default\n");
	  fprintf(f,"%s 1 %ld double\n",name,(long)net->nfaces[cell_type]);
	  if (binary)
	    {	   
	      fwriteBE(d,sizeof(double),net->nfaces[cell_type],f);
	      fprintf(f,"\n");
	    }
	  else
	    {
	      for (j=0;j<net->nfaces[cell_type];j++) fprintf(f,"%g\n",d[j]);
	    }	      
	}
    }
  if ((found)&&(xml)) fprintf(f,"</CellData>\n");



  
  if (xml)
    {
      fprintf(f,"    </Piece>\n");
      fprintf(f,"  </UnstructuredGrid>\n");

      if (binary)
	{
	  fprintf(f,"<AppendedData encoding=\"raw\">\n_");
	  fwrite(&size[0],sizeof(size[0]),1,f);
	  float v[net->ndims];
	  for (i=0;i<3;i++) v[i]=0;
	  for (i=0,j=0;i<net->nvertex;i++)
	    {
	      v[0]=net->v_coord[j++];

	      if (net->ndims_net>1) 
		{
		  v[1]=net->v_coord[j++];
		  if (net->ndims_net>2) v[2]=net->v_coord[j++];
		  else if (value!=NULL) v[2]=((value[i])-vmin)*dv;
		}
	      else if (value!=NULL) v[1]=((value[i])-vmin)*dv;
	      j+=dimStride;
	      fwrite(v,sizeof(float),3,f);
	    }
	  fwrite(&size[1],sizeof(size[1]),1,f);
	  fwrite(net->f_vertexIndex[cell_type],sizeof(NDNET_UINT),(cell_type+1)*ncells,f);
	  fwrite(&size[2],sizeof(size[2]),1,f);
	  NDNET_IDCUMT ct;
	  for (ct=(cell_type+1);ct<(cell_type+1)*(ncells+1);ct+=(cell_type+1))
	    fwrite(&ct,sizeof(NDNET_IDCUMT),1,f);	
	  fwrite(&size[3],sizeof(size[3]),1,f);
	  for (i=0;i<ncells;i++) fwrite(&vtk_cellT,sizeof(unsigned char),1,f);
	  for (i=0;i<nData;i++)
	    {
	      NDNET_VTK_HINT s;
	      s=sizeof(double)*dataSize[i];
	      fwrite(&s,sizeof(s),1,f);
	      fwrite(dataPtr[i],sizeof(double),dataSize[i],f);
	    }
	  fprintf(f,"</AppendedData>\n");
	  
	}     
    }
  
  if (xml) fprintf(f,"</VTKFile>\n");
   


  fclose(f);

  return 0;
}
