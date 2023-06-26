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
#ifndef __PDV_REMOTE_INTERFACE_HXX__
#define __PDV_REMOTE_INTERFACE_HXX__

#include <stdlib.h>
#include <stdio.h>

//#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "mystring.h"
#include "global.h"

#include <string>

class PDV_remote
{
public:

  static bool isFile(const char *fname)
  {
    int status;
    struct stat st_buf;
    
    status = stat (fname, &st_buf);
    if (status != 0) return 0;
    if (S_ISREG (st_buf.st_mode)) return 1;
    return 0;
  }

  static bool isPDV(const char *fname)
  {
    if (!isFile(fname)) return false;

    bool ret=false;
    std::string cmd(fname);
    cmd += std::string(" -sayHello");
    FILE *pdv=popen(cmd.c_str(),"r");
    char *txt=NULL;
    int nt=0;
    Mygetline (&txt,&nt,pdv);
    pclose(pdv);
    //printf("TEXT = %s\n",txt);
    if (std::string("Hi there.\n")==std::string(txt)) 
      ret=true;
    free(txt);

    return ret;
  }
  
  bool binaryFound()
  {
    return binFound;
  }

  std::string binary()
  {
    return binFName;
  }

  bool findBinary(std::string hint)
  {
    std::string tmp=find(hint);
    if (tmp==std::string(""))  
      {
	if (binFound) return false;
      }
    else 
      {
	binFName=tmp;
	binFound=true;
      }

    return binFound;
  }

  PDV_remote(bool init)
  {
    if (init)
      {
	binFName=find("./");
	if (binFName==std::string(""))  binFound=false;
      }
    else 
      {
	binFName=std::string("");
	binFound=false;
      }
  }

  PDV_remote(std::string hint=std::string("./"))
  {
    binFName=find(hint);
    if (binFName==std::string(""))  binFound=false;
  }

  std::pair<double,double> getThreshold(std::string fname, double nsig=0, double cut=0, bool dtfe=false)
  {
    if (!binFound) return std::pair<double,double>(-1,-1);

    char cmd[1024];
    std::pair<double,double> ret(0,0);
    
    if (nsig>0) 
      sprintf(cmd,"%s %s -nsig %g -dtfe",binFName.c_str(),fname.c_str(),nsig);
    else if (cut>0) 
      sprintf(cmd,"%s %s -cut %g",binFName.c_str(),fname.c_str(),cut);
    else if (dtfe) 
      sprintf(cmd,"%s %s -dtfe",binFName.c_str(),fname.c_str());
    else
      sprintf(cmd,"%s %s",binFName.c_str(),fname.c_str());
    printf("(%s)",cmd);fflush(0);
    FILE *pdv=popen(cmd,"r");
    char *txt=NULL;
    char result[255];
    int nt=0;
    strcpy(result,"");
    while(!feof(pdv))
      {
	Mygetline (&txt,&nt,pdv);
	//printf("%s\n",txt);
	if (strstr(txt,"Selected threshold")) strcpy(result,txt);
      }
    pclose(pdv);
    free(txt);

    char *tok[100];
    nt=str2tok(result," ",1,0,0,tok);
    if (nt==4)
      {
	if (!strcmp(tok[2],"-nsig")) ret.first=atof(tok[3]);
	else if (!strcmp(tok[2],"-cut")) ret.second=atof(tok[3]);
      }
    
    return ret;
  }

private:

  std::string binFName;
  bool binFound;
  
  std::string find(std::string hint=std::string("./"))
  {
    bool pdv;
    char PDVfname[255];
    char path_[255];
    char *path=&path_[0];
    
    sprintf(PDVfname,"%s",hint.c_str());
    if (isPDV(PDVfname)) return(std::string(PDVfname));

    sprintf(PDVfname,"%s/pdview",hint.c_str());
    if (isPDV(PDVfname)) return(std::string(PDVfname));

    sprintf(PDVfname,"%s/pdview/pdview",hint.c_str());
    if (isPDV(PDVfname)) return(std::string(PDVfname));

    getPath(hint.c_str(),&path);	   
 
    sprintf(PDVfname,"%spdview",path);
    if (isPDV(PDVfname)) return(std::string(PDVfname));

    sprintf(PDVfname,"%spdview/pdview",path);
    if (isPDV(PDVfname)) return(std::string(PDVfname));

    FILE *wh=popen("which pdview","r");
    int nrd=fscanf(wh,"%s",PDVfname);
    pclose(wh);
    if (isPDV(PDVfname)) return(std::string(PDVfname));

    return std::string("");   
  }

};

#endif
