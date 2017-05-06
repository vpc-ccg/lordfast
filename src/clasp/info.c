
/*
 *  info.c
 *  nfo messages
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 08/26/2007 06:49:02 PM CEST
 *  
 *  SVN
 *  Revision of last commit: $Rev: 114 $
 *  Author: $Author: steve $
 *  Date: $Date: 2010-04-19 15:59:41 +0200 (Mon, 19 Apr 2010) $
 *
 *  Id: $Id: info.c 114 2010-04-19 13:59:41Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/info.c $
 *  
 */
 
 #include <stdarg.h>
 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <time.h>
 #include "info.h"
 #include "debug.h"

 FILE *nfodevice = NULL;


 char *timestr_r(const struct tm *timeptr) {
    static const char wday_name[7][3] = {
        "Sun", "Mon", "Tue", "Wed",
        "Thu", "Fri", "Sat"
    };

    static const char mon_name[12][3] = {
        "Jan", "Feb", "Mar", "Apr", "May", "Jun",
        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
    };

    static char result[26];

    sprintf(result, "%.3s %.3s%3d %.2d:%.2d:%.2d %d",
        wday_name[timeptr->tm_wday], mon_name[timeptr->tm_mon],
        timeptr->tm_mday, timeptr->tm_hour, timeptr->tm_min,
        timeptr->tm_sec, 1900 + timeptr->tm_year);

    return result;
 }

 int
 infomsg( char *file, 
          int line, 
          const char *fmt, ...) {

   int ret;
   va_list ap;
   time_t rawtime;
   struct tm *timeinfo;
  
   if (mute) return 0;

    time(&rawtime);
    timeinfo = localtime (&rawtime);

   if (nfodevice == NULL) {
     nfodevice = NFODEVICE;
   }
   
   va_start(ap, fmt);
#ifdef PROGNFO   
   fprintf(nfodevice, "[%s] %s: ", "SEGEMEHL", timestr_r(timeinfo));
#endif
#ifdef PROG2NFO
   fprintf(nfodevice, "[%s] %s: ", "MaFIn", timestr_r(timeinfo));
#endif
#ifdef PROG3NFO
   fprintf(nfodevice, "[%s] %s: ", "clasp", timestr_r(timeinfo));
#endif
   ret = vfprintf(nfodevice, fmt, ap);
   va_end(ap);

   return ret; 
 }


void 
setnfodevice(char *filename) {
  FILE *fp;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    DBG("Couldn't open file '%s'. Exit forced.\n", filename);
    exit(-1);
  }

  nfodevice = fp;
}

int
nfolevel( char *file,
    int line,
    int level,
    const char *fmt, ... ) {

  int ret=0;
  va_list ap;
  time_t rawtime;
  struct tm *timeinfo;
   
  if (mute) return 0;
  
  time(&rawtime);
  timeinfo = localtime (&rawtime);

   if (nfodevice == NULL) {
     nfodevice = NFODEVICE;
   }
  
   if (NFOLEVEL >= level) {

    va_start(ap, fmt);
#ifdef PROGNFO
    fprintf(nfodevice, "[%s] %s: ", "SEGEMEHL", timestr_r(timeinfo));
#endif
    ret = vfprintf(nfodevice, fmt, ap);
    va_end(ap);
  }

  return ret;
}


