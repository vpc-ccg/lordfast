 #ifndef INFO_H
 #define INFO_H

/*
 *
 *	info.h
 *  nfo messages
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 08/26/2007 07:17:44 PM CEST  
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: info.h 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/info.h $
 */

 #include <stdarg.h>
 #include <stdio.h>
 #include <string.h>

#ifndef NFOLEVEL
#define NFOLEVEL 0
#endif

#ifndef NFODEVICE
#define NFODEVICE stderr
#endif

#define NFOL(L, X, ... ) debuglevel (__FILE__, __LINE__, L, X, __VA_ARGS__)
#define NFO(X, ...) infomsg(__FILE__, __LINE__, X, __VA_ARGS__)
#define INFO(X, ...) infomsg(__FILE__, __LINE__, X, __VA_ARGS__)
#define MSG(X) infomsg(__FILE__, __LINE__, X)

extern unsigned char mute;

int infomsg(char *, int, const char *fmt, ...);
int infolevel(char *, int, int, const char *fmt, ...);

 #endif
