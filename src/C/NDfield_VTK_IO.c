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

#include "NDfield_VTK_IO.h"
#include "NDfield.h"
#include "global.h"
#include "myendian.h"

int IsVTKNDfield(const char *fname)
{
  return 0;
}

NDfield *Load_NDfieldFromVTK(const char *fname)
{
  NDfield *field=NULL;

  return field;
}

int Save_NDfieldToVTK(NDfield *field,const char *fname,int type)
{
  int xml=type&NDFIELD_VTK_XML;
  int binary=type&NDFIELD_VTK_BIN;
  long i;
  FILE *f;
  char format[255];  
  long nx=1,ny=1,nz=1;
  double x0=0,y0=0,z0=0;
  double x1=0,y1=0,z1=0;
  //double x2=0,y2=0,z2=0;
  double dx,dy,dz;

  unsigned int size[50];
  int sizeindex=0;
  long offset=0;

  if (field->fdims_index)
    {
      fprintf(stderr,"ERROR writing VTK file: NDfield must be of type 'grid', current type is 'coords'.\n");
      return -1;
    }

  x0=x1=field->x0[0];
  //x2=field->x0[0]+field->delta[0];
  dx=field->delta[0]/field->dims[0];
  nx=field->dims[0];
  if (field->ndims>1) 
    {
      y0=y1=field->x0[1];
      //y2=field->x0[1]+field->delta[1];
      dy=field->delta[1]/field->dims[1];
      ny=field->dims[1];
    }
  else dy=dx;
  if (field->ndims>2)
    {
      z0=z1=field->x0[2];
      //z2=field->x0[2]+field->delta[2];
      dz=field->delta[2]/field->dims[2];
      nz=field->dims[2];
    }
  else dz=dy;
    
  long ncells=nx*ny*nz;

  if (xml)
    sprintf(format,"%s",(binary)?"appended":"ascii");
  else
    sprintf(format,"%s",(binary)?"BINARY":"ASCII");

  f=fopen(fname,"w");
   if (xml)
    {
      fprintf(f,"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
      /*
      fprintf(f,"  <ImageData WholeExtent=\"%g %g %g %g %g %g\" Origin=\"%g %g %g\" Spacing=\"%g %g %g\">\n",
	      x1,x2,y1,y2,z1,z2,x0,y0,z0,dx,dy,dz);
      fprintf(f,"    <Piece Extent=\"%g %g %g %g %g %g\">",x1,x2,y1,y2,z1,z2);
      */
      fprintf(f,"  <ImageData WholeExtent=\"0 %ld 0 %ld 0 %ld\" Origin=\"%g %g %g\" Spacing=\"%g %g %g\">\n",
	      nx,ny,(field->ndims<3)?0:nz,x0,y0,z0,dx,dy,dz);
      fprintf(f,"    <Piece Extent=\"0 %ld 0 %ld 0 %ld\">",nx,ny,(field->ndims<3)?0:nz);
    }
  else
    {
      fprintf(f,"# vtk DataFile Version 2.0\n");
      fprintf(f,"%s\n",field->comment);
      fprintf(f,"%s\n",format);
      fprintf(f,"DATASET STRUCTURED_POINTS\n");  
      fprintf(f,"DIMENSIONS %ld %ld %ld\n",(nx>1)?nx+1:(long)1,(ny>1)?ny+1:(long)1,(nz>1)?nz+1:(long)1);
      fprintf(f,"SPACING %g %g %g\n",dx,dy,dz);
      fprintf(f,"ORIGIN %g %g %g\n",x0,y0,z0);      
    }

   // DATA
   if (xml) 
    {
      fprintf(f,"<CellData>\n");
      if (binary)
	{	  
	  fprintf(f,"<DataArray Name=\"%s\" type=\"Float64\" format=\"%s\" offset=\"%ld\"/>\n","data",format,offset);
	  size[sizeindex]=ncells*sizeof(double);
	  offset+=size[sizeindex++]+sizeof(int);
	}
      else
	fprintf(f,"<DataArray Name=\"%s\" type=\"Float64\" format=\"%s\" >\n","data",format);
	
    }
   else
     {
       fprintf(f,"CELL_DATA %ld\n",ncells);
       fprintf(f,"FIELD FieldData %ld\n",(long)1);
       fprintf(f,"%s 1 %ld double\n","data",ncells);
     }

   if ((!binary)||(!xml))
     {
       DEFINE_NDVARS(d);
       SETNDPOINTER(d,field->val,field->datatype);

       for (i=0;i<ncells;i++)
	 {
	   double val;
	   SETNDFIELD_VAL(val,d,i,field->datatype);
	   if (binary) fwriteBE(&val,sizeof(double),1,f);
	   else fprintf(f,"%g\n",val);
	 }
     }

   if (xml)
     {
       if (!binary) fprintf(f,"</DataArray>\n");
      fprintf(f,"</CellData>\n");
     }

   if (xml)
    {      
      fprintf(f,"    </Piece>\n");
      fprintf(f,"  </ImageData>\n");

      if (binary)
	{
	  sizeindex=0;
	  fprintf(f,"<AppendedData encoding=\"raw\">\n_");
	  fwrite(&size[sizeindex],sizeof(size[sizeindex]),1,f);sizeindex++;
	  DEFINE_NDVARS(d);
	  SETNDPOINTER(d,field->val,field->datatype);

	  for (i=0;i<ncells;i++)
	    {
	      double val;
	      SETNDFIELD_VAL(val,d,i,field->datatype);
	      fwrite(&val,sizeof(double),1,f);
	    }

	  fprintf(f,"</AppendedData>\n");
	}
    }

   if (xml) fprintf(f,"</VTKFile>\n");
   fclose(f);
   return 0;
  
}
