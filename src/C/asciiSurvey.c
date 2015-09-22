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
#include "asciiSurvey.h"
#include "distances.h"

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <xlocale.h>
#include "mystring.h"

void freeSurvey(asciiSurvey **S)
{
  int i;

  free((*S)->m);
  free((*S)->pos);
  free((*S)->vel);
  for (i=0;i<(*S)->nfield;i++)
    {
      free((*S)->field[i]);
      free((*S)->fieldName[i]);
    }
  free((*S)->field);
  free(*S);
  *S=NULL;
  
}

int surveyFieldId(asciiSurvey *S, const char *fname)
{
  int i;
  for (i=0;i<S->nfield;i++)
    if (strcmp(fname,S->fieldName[i])) break;
  
  if (i==S->nfield) return -1;

  return i;
}

double *surveyFieldValue(asciiSurvey *S,const char *fname)
{
  return S->field[surveyFieldId(S,fname)];
}

inline double survey_red2dist(double z)
{
  if (!cosmoD_initialized())
    cosmoD_init(OMEGAM_DEFAULT,OMEGAL_DEFAULT,OMEGAK_DEFAULT,HUBBLE_DEFAULT,W_DEFAULT,NULL);
  //printf("z=%e -> %e\n",z,cosmoD_z2d(z));
  return cosmoD_z2d(z);
      
  /*
    double H0=HUBBLE_DEFAULT;
    double c=299792.458;
    double q0=0.5*Om - Ol;
    
    //return 1000.*c*(z*q0+(1-q0)*(1-sqrt(1+2*z*q0)))/(H0*q0*q0);
    return c*(z+0.5*(1-q0)*z*z)/H0;
  */
}

inline double survey_dist2red(double d)
{
  if (!cosmoD_initialized())
    cosmoD_init(OMEGAM_DEFAULT,OMEGAL_DEFAULT,OMEGAK_DEFAULT,HUBBLE_DEFAULT,W_DEFAULT,NULL);

  return cosmoD_d2z(d);

  /*
    double H0=HUBBLE_DEFAULT;
    double c=299792.458;
    double q0=0.5*Om - Ol;
    double B=c/H0;
    double A=0.5*c/H0*(1.-q0);
    double C=-d;
    double D=B*B-4.*A*C;
    
    return (-B+sqrt(D))/(2*A);
    //return 1000.*c*(z*q0+(1-q0)*(1-sqrt(1+2*z*q0)))/(H0*q0*q0);
    //return c/H0*(sqrt(1.+2.*H0/c*(1.-q0)*d)-1.);
    */
}


int isAsciiSurvey(const char *fname)
{
  char *line=NULL;
  FILE *fd;
  int i;

  fd=fopen(fname,"r");
  if (fd==NULL) return 0;
  Mygetline(&line,&i,fd);
  fclose(fd);
  if (((strstr(line,"px")!=NULL)&&(strstr(line,"py")!=NULL))||
      ((strstr(line,"ra")!=NULL)&&(strstr(line,"dec")!=NULL)))
    {
      free(line);
      return 1;
    }
  
  free(line);
  return 0;
}

asciiSurvey *readAsciiSurvey_header(const char *fname)
{
    asciiSurvey *S=calloc(1,sizeof(asciiSurvey));
    FILE *fd;
    int i;
    char *line=NULL;
    int linesize=0;
    char *token[100];
    int nsep=5;
    char sep[nsep];
  
    S->px=-1;S->py=-1;S->pz=-1;
    S->vx=-1;S->vy=-1;S->vz=-1;
    S->id=-1;
    S->ra=-1;S->dec=-1;
    S->dist=-1;S->z=-1;
    S->mass=-1;

    sep[0]='\0';
    sep[1]=' ';
    sep[2]=',';
    sep[3]=9;
    sep[4]='\n';

    if(!(fd=fopen(fname,"r")))
	return 0;
    
    Mygetline(&line,&linesize,fd);
    S->nfield=str2tok(line,sep,nsep,0,'#',token);
    S->nfield_unknown=0;
    S->field = calloc (S->nfield , sizeof(double *));

    for (i=0;i<S->nfield;i++) 
    {
	S->fieldName[i] = malloc(strlen(token[i])+1);
	strcpy(S->fieldName[i],token[i]);
    }

    for (i=0;i<S->nfield;i++)
    {
	if (!strncmp(S->fieldName[i],"px",strlen(S->fieldName[i]))) S->px=i;
	else if (!strncmp(S->fieldName[i],"py",strlen(S->fieldName[i]))) S->py=i;
	else if (!strncmp(S->fieldName[i],"pz",strlen(S->fieldName[i]))) S->pz=i;
	else if (!strncmp(S->fieldName[i],"vx",strlen(S->fieldName[i]))) S->vx=i;
	else if (!strncmp(S->fieldName[i],"vy",strlen(S->fieldName[i]))) S->vy=i;
	else if (!strncmp(S->fieldName[i],"vz",strlen(S->fieldName[i]))) S->vz=i;
	else if (!strncmp(S->fieldName[i],"id",strlen(S->fieldName[i]))) S->id=i;
	else if (!strncmp(S->fieldName[i],"ra",strlen(S->fieldName[i]))) S->ra=i;
	else if (!strncmp(S->fieldName[i],"dec",strlen(S->fieldName[i]))) S->dec=i;
	else if (!strncmp(S->fieldName[i],"dist",strlen(S->fieldName[i]))) S->dist=i;
	else if (!strncmp(S->fieldName[i],"z",strlen(S->fieldName[i]))) S->z=i;
	else if (!strncmp(S->fieldName[i],"mass",strlen(S->fieldName[i]))) S->mass=i;
	else S->nfield_unknown++;

      //printf("%s\n",fieldname[i]);
    }

    fclose(fd);
    free(line);
    return S;
}

double *asciiSurveyGetCoords_d(asciiSurvey *S,double **pos)
{
    int i;
    double d;

    if ((S->px>=0)&&(S->py>=0))
      {
	if (S->ndims==2)
	  {
	    if (*pos==NULL) *pos=(double*)malloc(sizeof(double)*S->npart*S->ndims);
	    double *p=*pos;
	    for (i=0;i<S->npart;i++) 
	      {
		p[i*S->ndims+0]=S->field[S->px][i];
		p[i*S->ndims+1]=S->field[S->py][i];
	      }
	    return p;
	  }
	else if (S->pz>=0)
	  {
	    if (*pos==NULL) *pos=(double*)malloc(sizeof(double)*S->npart*S->ndims);
	    double *p=*pos;
	    for (i=0;i<S->npart;i++) 
	      {
		p[i*S->ndims+0]=S->field[S->px][i];
		p[i*S->ndims+1]=S->field[S->py][i];
		p[i*S->ndims+2]=S->field[S->pz][i];
	      }
	    return p;
	  }
	//else {fprintf(stderr,"Error converting ascii survey to cartesian coordinates.\n");return NULL;}
      }

    if ((S->ra==-1)||(S->dec==-1)) 
      {fprintf(stderr,"Error converting ascii survey to cartesian coordinates.\n");return NULL;}
    if ((S->ndims==3)&&((S->dist==-1)&&(S->z==-1))) 
      {fprintf(stderr,"Error converting ascii survey to cartesian coordinates.\n");return NULL;}

    if (*pos==NULL) *pos=(double*)malloc(sizeof(double)*S->npart*S->ndims);
    double *p=*pos;

    for (i=0;i<S->npart;i++)
    {
	double pi = 3.141592653589793116;
	double ra=S->field[S->ra][i];
	double dec=S->field[S->dec][i];
	
	if (S->ndims==2) d=1;
	else if (S->dist!=-1) d=S->field[S->dist][i];
	else d=survey_red2dist(S->field[S->z][i]);
	//printf("ra=%g dec=%g\n",ra,dec);
	p[S->ndims*i] = d*cos(ra*pi/180.)*cos(dec*pi/180.);
	p[S->ndims*i+1] = d*sin(ra*pi/180.)*cos(dec*pi/180.);

	if (S->ndims==3) p[S->ndims*i+2] = d*sin(dec*pi/180.);
    }
    
    return p;
}

asciiSurvey *readAsciiSurvey(const char *fname)
{
  asciiSurvey *S=calloc(1,sizeof(asciiSurvey));
  
  FILE *fd;
  long i,j;
  char *line=NULL;
  int linesize=0;
  char *token[100];
  //int nfields=0;
  int nsep=5;
  char sep[nsep];
  
  S->px=-1;S->py=-1;S->pz=-1;
  S->vx=-1;S->vy=-1;S->vz=-1;
  S->id=-1;
  S->ra=-1;S->dec=-1;
  S->dist=-1;S->z=-1;
  S->mass=-1;
  S->nfield_unknown=0;

  sep[0]='\0';
  sep[1]=' ';
  sep[2]=',';
  sep[3]=9;
  sep[4]='\n';
  
  printf("reading %s (ASCII) ...",fname);fflush(0);
  if(!(fd=fopen(fname,"r")))
    return 0;
  
  Mygetline(&line,&linesize,fd);
  S->nfield=str2tok(line,sep,nsep,0,'#',token);
  S->field = calloc (S->nfield , sizeof(double *));
  //indexof=malloc(S->nfield*sizeof(int));
  for (i=0;i<S->nfield;i++) 
    {
      S->fieldName[i] = malloc(strlen(token[i])+1);
      strcpy(S->fieldName[i],token[i]);
      for (j=0;j<strlen(S->fieldName[i]);j++)
	S->fieldName[i][j]=tolower(S->fieldName[i][j]);
    }
  
  for (i=0;i<S->nfield;i++)
    {
      if (!strncmp(S->fieldName[i],"px",strlen(S->fieldName[i]))) S->px=i;
      else if (!strncmp(S->fieldName[i],"py",strlen(S->fieldName[i]))) S->py=i;
      else if (!strncmp(S->fieldName[i],"pz",strlen(S->fieldName[i]))) S->pz=i;
      else if (!strncmp(S->fieldName[i],"vx",strlen(S->fieldName[i]))) S->vx=i;
      else if (!strncmp(S->fieldName[i],"vy",strlen(S->fieldName[i]))) S->vy=i;
      else if (!strncmp(S->fieldName[i],"vz",strlen(S->fieldName[i]))) S->vz=i;
      else if (!strncmp(S->fieldName[i],"id",strlen(S->fieldName[i]))) S->id=i;
      else if (!strncmp(S->fieldName[i],"ra",strlen(S->fieldName[i]))) S->ra=i;
      else if (!strncmp(S->fieldName[i],"dec",strlen(S->fieldName[i]))) S->dec=i;
      else if (!strncmp(S->fieldName[i],"dist",strlen(S->fieldName[i]))) S->dist=i;
      else if (!strncmp(S->fieldName[i],"z",strlen(S->fieldName[i]))) S->z=i;
      else if (!strncmp(S->fieldName[i],"mass",strlen(S->fieldName[i]))) S->mass=i;
      else S->nfield_unknown++;
      
      //printf("%s\n",fieldname[i]);
    }
//printf ("(%d fields)",nfields);fflush(0);exit(0);

//printf ("%d %d %d\n",idra,iddec,idd);

  for (j=0;j<S->nfield;j++)
    S->field[j] = malloc((S->npart+10000)*sizeof(double));

  Mygetline(&line,&linesize,fd);
  do
    {
      if (S->nfield != str2tok(line,sep,nsep,0,'#',token))
	{
	  printf ("Line %d is incomplete, SKIPPED.\n",S->npart+2);
	  continue;
	}
      
      for (i=0;i<S->nfield;i++)
	{	
	  S->field[i][S->npart] = atof(token[i]);
	}
      S->npart++;

      if (S->npart%10000 == 0)
	{
	  for (j=0;j<S->nfield;j++)
	    S->field[j] = realloc(S->field[j],(S->npart+10000)*sizeof(double));
	  printf("\rreading %s (ASCII) ... (%d lines)",fname,S->npart);fflush(0);
	}
      
      Mygetline(&line,&linesize,fd);
      
    } while (!feof(fd));
 
  for (j=0;j<S->nfield;j++)
    S->field[j] = realloc(S->field[j],(S->npart)*sizeof(double));
  
  if (S->mass!=-1)
    {
      S->m = malloc(S->npart*sizeof(float));
      for (i=0;i<S->npart;i++)
	S->m[i]=S->field[S->mass][i];
    }
 
  if ((S->px!=-1)&&(S->py!=-1))
    { 
      if (S->pz!=-1) S->ndims=3;
      else S->ndims=2;

      S->pos = malloc(S->ndims*(S->npart)*sizeof(float));

      for (i=0;i<S->npart;i++)
	{
	  S->pos[S->ndims*i]=S->field[S->px][i];
	  S->pos[S->ndims*i+1]=S->field[S->py][i];
	  if (S->ndims==3) 
	    S->pos[S->ndims*i+2]=S->field[S->pz][i];
	}
    }
  else if ((S->ra!=-1)&&(S->dec!=-1))
    { 
      if ((S->dist!=-1)||(S->z!=-1)) S->ndims=3;
      else S->ndims=2;

      S->pos = malloc(S->ndims*(S->npart)*sizeof(float));

      double ra,dec,d;

      for (i=0;i<S->npart;i++)
	{
	  ra=S->field[S->ra][i];
	  dec=S->field[S->dec][i];
	  
	  if (S->ndims==2) d=1;
	  else if (S->dist!=-1) d=S->field[S->dist][i];
	  //else d=S->field[S->z][i];
	  else d=survey_red2dist(S->field[S->z][i]);
	  //printf("HELLOO\n");
	  //printf("ra2=%g dec2=%g\n",ra,dec);
	  
	  /*
	  S->pos[S->ndims*i] = d*cos(ra*pi/180.)*cos(dec*pi/180.);
	  S->pos[S->ndims*i+1] = d*sin(ra*pi/180.)*cos(dec*pi/180.);
	  if (S->ndims==3)
	    S->pos[S->ndims*i+2] = d*sin(dec*pi/180.);
	  */
	  
	  S->pos[S->ndims*i]=ra;
	  S->pos[S->ndims*i+1]=dec;
	  S->pos[S->ndims*i+2]=d;
	  
	  
	  //printf ("[%g %g %g]->[%g %g %g]\n",ra,dec,d,P->Pos[i],P->Pos[i+1],P->Pos[i+2]);
	}
    }
 
  if ((S->vx!=-1)&&(S->vy!=-1)&&((S->ndims==2)||(S->vz!=-1)))
    {
      S->vel = malloc(S->ndims*(S->npart)*sizeof(float));
      
      for (i=0;i<S->npart;i++)
	{
	  S->vel[S->ndims*i]=S->field[S->vx][i];
	  S->vel[S->ndims*i+1]=S->field[S->vy][i];
	  if (S->ndims==3) 
	    S->vel[S->ndims*i+2]=S->field[S->vz][i];
	}
    }

  free(line);
  fclose(fd);
  
  printf("\rreading %s (ASCII) ... done. (%d lines)\n",fname,S->npart);
  fflush(0);

  return S;
}
