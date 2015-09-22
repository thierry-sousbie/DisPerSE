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
//#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "mystring.h"

//#ifndef _GNU_SOURCE
/* Default value for line length.  */
static const int line_size = 1024;

ssize_t
Mygetdelim (char **lineptr, int *n, int delim, FILE *stream)
{
  int indx = 0;
  int c;

  /* Sanity checks.  */
  if (lineptr == NULL || n == NULL || stream == NULL)
    return -1;

  /* Allocate the line the first time.  */
  if (*lineptr == NULL)
    {
      *lineptr = malloc (line_size);
      if (*lineptr == NULL)
        return -1;
      *n = line_size;
    }

  /* Clear the line.  */
  memset (*lineptr, '\0', *n);

  while ((c = getc (stream)) != EOF)
    {
      /* Check if more memory is needed.  */
      if (indx >= *n)
        {
          *lineptr = realloc (*lineptr, *n + line_size);
          if (*lineptr == NULL)
            {
              return -1;
            }
          /* Clear the rest of the line.  */
          memset(*lineptr + *n, '\0', line_size);
          *n += line_size;
        }

      /* Push the result in the line.  */
      (*lineptr)[indx++] = c;

      /* Bail out.  */
      if (c == delim)
        {
          break;
        }
    }
  return (c == EOF) ? -1 : indx;
}
//#endif


ssize_t Mygetline (char **lineptr, int *n, FILE *stream)
{
  //#ifndef _GNU_SOURCE
  return Mygetdelim (lineptr, n, '\n', stream);
  //#else
  //return getline(lineptr,  (size_t *)n,  stream);
  //#endif

}

int isFile(const char *fname)
{
  int status;
  struct stat st_buf;
    
  status = stat (fname, &st_buf);
  if (status != 0) return 0;
  if (S_ISREG (st_buf.st_mode)) return 1;

  return 0;
}


int str2tok(char *str,const char *sep, int nsep, char comment, char ignore, char** tok)
{
    int n=0;
    int i,j;
    int len = strlen(str);
    int last=0;
    
    for (i=0;i<len;i++)
    {
	for (j=0;j<nsep;j++)
	    if (str[i]==sep[j])
	    {
		str[i]='\0';
		continue;
	    }
	    else if (comment&(str[i]==comment))
	    {
		len=i;
		continue;
	    }
    }

    for (i=0;i<len;i++)
    {
	if ((str[i]!='\0')&&(str[i]!='\n')&&(str[i]!=ignore))
	{
	    if (last == 0) 
		tok[n++] = &(str[i]);
	    last=1;
	}
	else
	    last=0;
    }
  
    return n;
}

char *strReplace(char *st,char *dest,const char *orig,const char *repl) 
{
  char *ch;
  if (!(ch = strstr(st, orig)))
    return st;
  strncpy(dest, st, ch-st);  
  dest[ch-st] = 0;
  sprintf(dest+(ch-st), "%s%s", repl, ch+strlen(orig));
  return dest;
}

char *CutName(char *Name)
{
    int i,j; 

    for (i=0,j=-1;i<strlen(Name);i++) 
	if (Name[i]=='/') j=i;

    if (j!=-1)
	return &Name[j+1];
    else
	return Name;
}

char *getPath(const char *Name,char **path_)
{
    int i,j;
    char *path;

    if (*path_ == NULL)
      {
	path = (char *)malloc(strlen(Name)*sizeof(char));
	*path_=path;
      }
    else path=*path_;

    for (i=0,j=-1;i<strlen(Name);i++) 
	if (Name[i]=='/') j=i;

    if (j!=-1)
      {
	strcpy(path,Name);	
	path[j+1]='\0';
      }
    else
      strcpy(path,"");
      
    return path;
}

int isAscii(char *fname)
{
  char cmd[1024];
  sprintf(cmd,"file -b -L %s",fname);
  FILE *f=popen(cmd,"r");
  int n=fscanf(f,"%s",cmd);
  pclose(f);
  
  if (n!=1) {
    printf ("ERROR in isAscii, could not read output of cmd '%s'\n",cmd);    
    return -1;
  }
  
  if (strstr(cmd,"ASCII")) return 1;
  else return 0;
}
