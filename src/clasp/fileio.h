#ifndef FILEIO_H
#define FILEIO_H

/*
 * fileio.h
 * declarations for file io
 *
 * @author Steve Hoffmann
 * @date Sat 25 Nov 2006
 *
 *  SVN
 *  Revision of last commit: $Rev: 114 $
 *  Author: $Author: steve $
 *  Date: $Date: 2010-04-19 15:59:41 +0200 (Mon, 19 Apr 2010) $
 *
 *  Id: $Id: fileio.h 114 2010-04-19 13:59:41Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/fileio.h $
 */

#ifndef ALLOCMEMORY
	#include "memory.h"
#endif

#include "stringutils.h"

char* readfile(void *, char *, unsigned long long*);
stringset_t **readcsv(void *, char *, char*, Uint *);
char **readlines(char *filename, Uint *linecount);
void writeY(char *, double  *, Uint);
void writeXYUint(char *filename, Uint *X, Uint *Y, Uint len);


#endif
