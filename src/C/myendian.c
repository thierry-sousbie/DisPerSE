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
#include "myendian.h"

int fread_sw(void *data,size_t size,size_t nb,FILE *f,int swap)
{
  
  int ret = fread(data,size,nb,f);
  
  if (swap)
    {
      switch (size)
	{
	    case 8: Dswap8BArr(data,nb);break;
	    case 4: Dswap4BArr(data,nb);break;
	    case 2: Dswap2BArr(data,nb);break;
	}
    }
  
  return ret;
}


size_t freadBE(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t res;
  static int isLittle=-1;  
  if (isLittle<0)
    {
      int i=1;
      unsigned char *ic=(unsigned char*)&i;
      if (*ic) isLittle=1;
      else isLittle=0;
    }

  res=fread(ptr,size,nmemb,stream);
  if ((isLittle)&&(size>1))
    {
      long i,j;
      unsigned char a[16];
      unsigned char *cptr=(unsigned char*)ptr;
      for (i=0;i<nmemb*size;i+=size)
	{
	  for (j=0;j<size;j++) a[size-j-1]=cptr[i+j];
	  for (j=0;j<size;j++) cptr[i+j]=a[j];
	}
    }
  return res;
}


size_t fwriteBE(const void *ptr, size_t size, size_t nmemb,FILE *stream)
{
  size_t res;
  static int isLittle=-1;  
  if (isLittle<0)
    {
      int i=1;
      unsigned char *ic=(unsigned char*)&i;
      if (*ic) isLittle=1;
      else isLittle=0;
    }

  if ((isLittle)&&(size>1))
    {
      long i,j;
      unsigned char a[16];
      unsigned char *cptr=(unsigned char*)ptr;
      for (i=0;i<nmemb*size;i+=size)
	{
	  for (j=0;j<size;j++) a[size-j-1]=cptr[i+j];
	  for (j=0;j<size;j++) cptr[i+j]=a[j];
	}
    }

  res=fwrite(ptr,size,nmemb,stream);

  if ((isLittle)&&(size>1))
    {
      long i,j;
      unsigned char a[16];
      unsigned char *cptr=(unsigned char*)ptr;
      for (i=0;i<nmemb*size;i+=size)
	{
	  for (j=0;j<size;j++) a[size-j-1]=cptr[i+j];
	  for (j=0;j<size;j++) cptr[i+j]=a[j];
	}
    }
  return res;
}


int swapI(int val)
{
    int out;
    const char *i=(const char *)&val;
    char *o=(char *)&out;

    o[3]=i[0];
    o[2]=i[1];
    o[1]=i[2];
    o[0]=i[3];

    return out; 
}

float swapF(float val)
{
    float out;
    const char *i=(const char *)&val;
    char *o=(char *)&out;

    o[3]=i[0];
    o[2]=i[1];
    o[1]=i[2];
    o[0]=i[3];

    return out; 
}

double swapD(double val)
{
    double out;
    const char *i=(const char *)&val;
    char *o=(char *)&out;

    o[7]=i[0];
    o[6]=i[1];
    o[5]=i[2];
    o[4]=i[3];
    o[3]=i[4];
    o[2]=i[5];
    o[1]=i[6];
    o[0]=i[7];


    return out; 
}

inline void Dswap2B(void *val)
{
    char *c=(char *)val;
    char a;
    
   a=c[0];c[0]=c[1];c[1]=a; 
}

inline void Dswap4B(void *val)
{
    char *c=(char *)val;
    char a;
    
   a=c[0];c[0]=c[3];c[3]=a;
   a=c[1];c[1]=c[2];c[2]=a; 
 
}

inline void Dswap8B(void *val)
{
   char *c=(char *)val;
   char a;
    
   a=c[0];c[0]=c[7];c[7]=a;
   a=c[1];c[1]=c[6];c[6]=a;
   a=c[2];c[2]=c[5];c[5]=a;
   a=c[3];c[3]=c[4];c[4]=a;
}

void Dswap2BArr(void *val,size_t n)
{
    size_t i;
    char a;

    char *c=(char *)val;

    for (i=0;i<2*n;i+=2)
    {
	a=c[i];
	c[i]=c[i+1];
	c[i+1]=a;
    }

}


void Dswap4BArr(void *val,size_t n)
{
    size_t i;
    char a,b;

    char *c=(char *)val;

    for (i=0;i<4*n;i+=4)
    {
	a=c[i];
	b=c[i+1];
	c[i]=c[i+3];
	c[i+1]=c[i+2];
	c[i+2]=b;
	c[i+3]=a;
    }

}

void Dswap8BArr(void *val,size_t n)
{
    size_t i;
    char a,b,u,v;

    char *c=(char *)val;

    for (i=0;i<8*n;i+=8)
    {
	a=c[i];
	b=c[i+1];
	u=c[i+2];
	v=c[i+3];
	c[i]=c[i+7];
	c[i+1]=c[i+6];
	c[i+2]=c[i+5];
	c[i+3]=c[i+4];
	c[i+4]=v;
	c[i+5]=u;
	c[i+6]=b;
	c[i+7]=a;
    }

}
