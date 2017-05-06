
/*
 *  debug.c
 *  debug messages
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 08/26/2007 06:49:02 PM CEST
 *  
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: debug.c 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/debug.c $
 *  
 */
 
 #include <stdarg.h>
 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include "debug.h"

 FILE *dbgdevice = NULL;

 int
 debugmsg( char *file, 
           int line, 
           const char *fmt, ...) {

   if (dbgdevice == NULL) {
     dbgdevice = DBGDEVICE;
   }
   
   int ret;
   va_list ap;
   va_start(ap, fmt);
#ifdef DBGNFO   
   fprintf(dbgdevice, "[%s] file: %s, line: %d: ", "SEGEMEHL", file, line);
#endif
   ret = vfprintf(dbgdevice, fmt, ap);
   va_end(ap);

   return ret; 
 }


void 
setdebugdevice(char *filename) {
  FILE *fp;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    DBG("Couldn't open file '%s'. Exit forced.\n", filename);
    exit(-1);
  }

  dbgdevice = fp;
}

int
debuglevel( char *file,
    int line,
    int level,
    const char *fmt, ... ) {

  int ret=0;
  va_list ap;

   if (dbgdevice == NULL) {
     dbgdevice = DBGDEVICE;
   }
  
   if (DBGLEVEL >= level) {

    va_start(ap, fmt);
#ifdef DBGNFO
    fprintf(dbgdevice, "[%s] file: %s, line: %d: ", "segemehl", file, line);
#endif
    ret = vfprintf(dbgdevice, fmt, ap);
    va_end(ap);
  }

  return ret;
}




