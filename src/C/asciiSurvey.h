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
#ifndef __ASCII_SURVEY_H__
#define __ASCII_SURVEY_H__

#ifdef __cplusplus
extern "C" {
#endif

#define OMEGAM_DEFAULT (0.27)
#define OMEGAL_DEFAULT (0.73)
#define OMEGAK_DEFAULT (0.0)
#define W_DEFAULT (-1.0)
#define HUBBLE_DEFAULT (0.72)

typedef struct asciiSurvey_str { 
  int ndims;
  
  int nfield;
  int nfield_unknown;
  int npart;

  float *pos;
  float *vel;
  float *m;

  int px,py,pz;
  int vx,vy,vz;
  int id;
  int ra,dec;
  int dist,z;
  int mass;

  double** field;
  char* fieldName[255];
} asciiSurvey;

  void freeSurvey(asciiSurvey **S);
  int surveyFieldId(asciiSurvey *S, const char *fname);
  double *surveyFieldValue(asciiSurvey *S, const char *fname);
  int isAsciiSurvey(const char *fname);
  //int asciiSurveyGetSphericalCoords(asciiSurvey *S,float *pos);
  double *asciiSurveyGetCoords_d(asciiSurvey *S,double **pos);
  asciiSurvey *readAsciiSurvey(const char *fname);
  asciiSurvey *readAsciiSurvey_header(const char *fname);

  double survey_red2dist(double z);
  double survey_dist2red(double d);
    
#ifdef __cplusplus
}
#endif

#endif
