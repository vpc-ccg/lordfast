/*
 * fileio.c
 * functions to manipulate and read files
 *
 *  SVN
 *  Revision of last commit: $Rev: 114 $
 *  Author: $Author: steve $
 *  Date: $Date: 2010-04-19 15:59:41 +0200 (Mon, 19 Apr 2010) $
 *
 *  Id: $Id: fileio.c 114 2010-04-19 13:59:41Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/fileio.c $
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "basic-types.h"
#include "fileio.h"

#ifndef DIR_SEPARATOR
#define DIR_SEPARATOR '/'
#endif

#if defined (_WIN32) || defined (__MSDOS__) || defined (__DJGPP__) || \
  defined (__OS2__)
#define HAVE_DOS_BASED_FILE_SYSTEM
#ifndef DIR_SEPARATOR_2 
#define DIR_SEPARATOR_2 '\\'
#endif
#endif

/* Define IS_DIR_SEPARATOR.  */
#ifndef DIR_SEPARATOR_2
# define IS_DIR_SEPARATOR(ch) ((ch) == DIR_SEPARATOR)
#else /* DIR_SEPARATOR_2 */
# define IS_DIR_SEPARATOR(ch) \
  (((ch) == DIR_SEPARATOR) || ((ch) == DIR_SEPARATOR_2))
#endif /* DIR_SEPARATOR_2 */



char* dirname(const char *filename) {
    char *s;

    s=strrchr(filename, (int)'/');
    if(s && *s)
      *s = '\0';

    return s;
}



  char *
basename (const char *name)
{
  const char *base;

#if defined (HAVE_DOS_BASED_FILE_SYSTEM)
  /* Skip over the disk name in MSDOS pathnames. */
  if (ISALPHA (name[0]) && name[1] == ':') 
    name += 2;
#endif

  for (base = name; *name; name++)
  {
    if (IS_DIR_SEPARATOR (*name))
    {
      base = name + 1;
    }
  }
  return (char *) base;
}

char* readfile(void* space, char* filename, unsigned long long* strlen) {
  char ch, *buffer;
  FILE *fp;
  unsigned long long buffersize = MAXBUFFERSIZE;
  unsigned long long len=0;

  fp = fopen(filename, "r");
  if (fp == NULL){
    fprintf(stderr, "Opening of file %s failed. Exit forced.\n", filename);
    exit(EXIT_FAILURE);
  }

  buffer = ALLOCMEMORY(space, NULL, char, buffersize);

  while((ch=getc(fp)) != EOF) {
    if(len == buffersize-1) {
      buffersize = buffersize + MAXBUFFERSIZE;
      buffer = ALLOCMEMORY(space, buffer, char, buffersize);
    }
    len++;
    buffer[len-1]=(char)ch;	
  }
  buffer[len]='\0';
  fclose(fp);

  *strlen = len;
  return buffer;
}

char** readlines(char *filename, Uint *linecount){
  char ch, *buffer, **lines = NULL;
  FILE *fp;
  unsigned long long buffersize;
  unsigned long long len;
  *linecount = 0;
  len = 0;
  buffersize = MAXBUFFERSIZE;
  buffer = (char *) malloc(buffersize * sizeof(char));
  fp = fopen(filename, "r");
  if (fp == NULL){
    fprintf(stderr, "Opening of file %s failed. Exit forced.\n", filename);
    exit(-1);
  }
  while((ch=getc(fp)) != EOF) {
    if (len == buffersize-1) {
      buffersize = buffersize + MAXBUFFERSIZE;
      buffer = (char *) realloc(buffer, buffersize * sizeof(char));
      if (buffer == NULL){
	fprintf(stderr, "Could not allocate sufficient amount of memory for \
buffer. Exit forced.\n");
	exit(-1);
      }
    }
    if (ch == '\n' || ch == '\r'){
      /* dos based files */
      if (ch == '\n' && buffer[len - 1] == '\r'){
	len--;
      }
      if (len == 0){
	continue;
      }
      buffer[len] = '\0';
      *linecount = *linecount + 1;
      lines = (char **) realloc (lines, *linecount * sizeof(char *));
      if (lines == NULL){
	fprintf(stderr, "Could not allocate sufficient amount of memory for \
char **. Exit forced.\n");
	exit(-1);
      }
      buffer = (char *) realloc(buffer, (len + 1) * sizeof(char));
      lines[*linecount - 1] = buffer;
      buffersize = MAXBUFFERSIZE;
      buffer = (char *) malloc(buffersize * sizeof(char));
      if (buffer == NULL){
	fprintf(stderr, "Could not allocate sufficient amount of memory for \
buffer.	Exit forced.\n");
	exit(-1);
      }
      len = 0;
    }
    else {
      len++;
      buffer[len-1] = (char) ch;
    }
  }
  if (len > 0){
    buffer[len] = '\0';
    *linecount = *linecount + 1;
    lines = (char **) realloc (lines, *linecount * sizeof(char *));
    if (lines == NULL){
      fprintf(stderr, "Could not allocate sufficient amount of memory for \
char **. Exit forced.\n");
      exit(-1);
    }
    buffer = (char *) realloc(buffer, (len + 1) * sizeof(char));
    lines[*linecount-1] = buffer;
  }
  else {
    free(buffer);
  }
  fclose(fp);
  return(lines);
}

stringset_t **
readcsv(void *space, 
    char* filename, 
    char *delim, 
    Uint *linecount) {

  Uint i;
  unsigned long long contentlen;
  char *content;
  stringset_t *lines, **csv;

  content = readfile(space, filename, &contentlen);
#ifndef _CRLF_
  lines = tokensToStringset(space, "\n", content, contentlen);
#else
  lines = tokensToStringset(space, "\r\n", content, contentlen);
#endif
  FREEMEMORY(space, content);
  *linecount=lines->noofstrings;
  csv=ALLOCMEMORY(space, NULL, stringset_t *, lines->noofstrings);

  for(i=0; i < lines->noofstrings; i++) {
    csv[i] = tokensToStringset(space, delim, lines->strings[i].str, lines->strings[i].len);
  }

  destructStringset(space, lines);	
  return csv;
}

void
writeY(char *filename, double *Y, Uint len) {
  FILE *file;
  Uint i;

  file = fopen(filename, "w");
  if (file == NULL) {
    fprintf(stderr, "couldn't open %s - exit forced", filename);
    exit(-1);
  }

  for(i=0; i < len; i++) {
    fprintf(file,"%d\t%f\n", i, Y[i]);
  }

  fclose(file);
  return;
}

void
writeYUint(char *filename, Uint *Y, Uint len) {
  FILE *file;
  Uint i;

  file = fopen(filename, "w");
  if (file == NULL) {
    fprintf(stderr, "couldn't open %s - exit forced", filename);
    exit(-1);
  }

  for(i=0; i < len; i++) {
    fprintf(file,"%d\t%d\n", i, Y[i]);
  }

  fclose(file);
  return;
}

void 
writeXYUint(char *filename, Uint *X, Uint *Y, Uint len) {
  FILE *file;
  Uint i;

  file = fopen(filename, "w");
  if (file == NULL) {
    fprintf(stderr, "couldn't open %s - exit forced", filename);
    exit(-1);
  }

  for(i=0; i < len; i++) {
    fprintf(file,"%d\t%d\t%d\n", i, X[i], Y[i]);
  }

  fclose(file);
}


