/*
 * stringutils.c
 * functions to manipulate strings
 *
 *  SVN
 *  Revision of last commit: $Rev: 107 $
 *  Author: $Author: steve $
 *  Date: $Date: 2009-06-02 11:27:30 +0200 (Tue, 02 Jun 2009) $
 *
 *  Id: $Id: stringutils.c 107 2009-06-02 09:27:30Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/stringutils.c $
 */

 #include <stdlib.h>
 #include <stdio.h>
 #include <string.h>
 #include <assert.h>
 #include "stringutils.h"
 #include "basic-types.h"


void
printSubseq(char* seq, Uint start, Uint end) {
  int i;
  assert(end <= strlen(seq));
  for(i=start; i <= end; i++) {
    printf("%c", seq[i]);
  }
}

char* strrev(char *str, Uint len){
  int end = len-1;
  int start = 0;

  while (start<end) {
    str[start] ^= str[end];
    str[end] ^= str[start];
    str[start] ^= str[end];
    start++;
    end--;
  }
  return str;
}

stringset_t* tokensToStringset(void *space, char* delim, char* toTokens, 
								Uint len){
	Uint toklen;
	char* token;
	char* buffer;
	stringset_t *set;
	
	set = ALLOCMEMORY(space, NULL, stringset_t, 1);
	set->noofstrings = 0;
	set->strings = NULL;
									
	if (toTokens == NULL || len == 0)
	return set;
													
	buffer = ALLOCMEMORY(space, NULL, char, len+1);
	buffer = memcpy(buffer, toTokens, len+1);
																
	if (buffer == NULL) {
		fprintf(stderr, "copy tokenstring %s to buffer failed.\n", toTokens);
		exit(-1);	
	}
	
	token = strtok(buffer, delim);
	
	while(token != NULL) {
		
			toklen = strlen(token);
			set->noofstrings++;
			set->strings = ALLOCMEMORY(space, set->strings, string_t, set->noofstrings);
			set->strings[set->noofstrings-1].str = ALLOCMEMORY(space, NULL, char, toklen+1);
			set->strings[set->noofstrings-1].str = memcpy(set->strings[set->noofstrings-1].str, token, toklen);
			set->strings[set->noofstrings-1].str[toklen]='\0'; 
			set->strings[set->noofstrings-1].len = toklen;
			token = strtok(NULL, delim);
	}

	FREEMEMORY(space, buffer);
	return set;
}


char* strtrimquote (void *spacetab, char *toTrim, Uint *len) {
	Uint i=0;
	int start=-1;
	int end =-2;	
	char* trimmed = NULL;
	
	for(i=0; i < *len; i++) {
	  
		if(ISQUOTE(toTrim[i])) {	
			continue;
		}
		else if(start == -1) {
			start=i;
			end=i;
		}
		else 
			end=i;
	}

	if(start >= 0) {
		trimmed = ALLOCMEMORY(spacetab, NULL, char, (end-start)+2);
   		memmove(trimmed, &toTrim[start], (end-start)+1);
		trimmed[(end-start)+1]='\0';
	}
	
	*len = (end-start)+1;
	return trimmed;

}


char* strtrim(void *spacetab, char* toTrim, Uint *len) {	
	Uint i=0;
	int start=-1;
	int end =-2;	
	char* trimmed = NULL;
	
	for(i=0; i < *len; i++) {
	  
		if(ISWHITESPACE(toTrim[i])) {	
			continue;
		}
		else if(start == -1) {
			start=i;
			end=i;
		}
		else 
			end=i;
	}

	if(start >= 0) {
		trimmed = ALLOCMEMORY(spacetab, NULL, char, (end-start)+2);
   		memmove(trimmed, &toTrim[start], (end-start)+1);
		trimmed[(end-start)+1]='\0';
	}
	
	*len = (end-start)+1;
	return trimmed;
}

char* concat(void *spacetab, char* strA, char* strB, int lenA, int lenB) {

	if(strB == NULL || lenB == 0) 
	  		return strA;
	if(strA == NULL || lenA == 0)
		  		return strB;
			
	strA=ALLOCMEMORY(spacetab, strA, char, lenA+lenB+1);
	memmove(&strA[lenA], strB, lenB);
					
	strA[lenA+lenB]='\0';
					
	return strA;
}

char* concatdelim(void *spacetab, char* strA, char* strB, int lenA, int lenB, char delim) {

	if(strB == NULL || lenB == 0) 
	  		return strA;
	if(strA == NULL || lenA == 0)
		  	return strB;
			
	strA=ALLOCMEMORY(spacetab, strA, char, lenA+lenB+2);
	strA[lenA]=delim;
	memmove(&strA[lenA+1], strB, lenB);
						
	strA[lenA+lenB+1]='\0';
						
	return strA;						
}


void destructStringset(void *space, stringset_t *s) {
	Uint i;
	
	if (s->strings) {
		for(i=0; i < s->noofstrings; i++) {
			if(s->strings[i].str != NULL)	
				FREEMEMORY(space, s->strings[i].str);
		}

		FREEMEMORY(space, s->strings);
	}

	FREEMEMORY(space, s);
}


stringset_t *initStringset(void *space) {
	stringset_t *set;

	set = ALLOCMEMORY(space, NULL, stringset_t, 1);
	set->noofstrings=0;
	set->strings = NULL;

	return set;	
}

void addString(void *space, stringset_t *set, char *string, Uint len) {

	set->noofstrings++;
	set->strings=ALLOCMEMORY(space, set->strings, string_t, set->noofstrings);
	set->strings[set->noofstrings-1].str = string;
	set->strings[set->noofstrings-1].len = len;
}



/*---------------------------------- strrev ----------------------------------
 *    
 * reversing a string
 * 
 */
 
char*
strreverse(char *s, Uint len)
{	
    Uint i;
	char resc;
	
  	for(i=0; i < (len/2); i++) {
		resc = s[i];
		s[i] = s[len-1-i];
		s[len-1-i] = resc;
	}
	
	return s;
}

/* -------------------------------- my_itoa  ---------------------------------
 *    
 * just in case that there's no itoa
 * 
 */
 
char*
my_itoa (int value, char *buffer, Uint radix)
{
  	const char* base ="0123456789abcdef";
    int i=0;

	if (value == 0) {
	  buffer[0]=base[0];
	  buffer[1]='\0';
	  return buffer;
	}
	
	for (i=0; value > 0; i++, value /= radix)
	  buffer[i] = base[(value % radix)];
	
	buffer[i] ='\0';

	buffer = strreverse(buffer, i);
	return buffer;
}



/*-------------------------------- attachext ---------------------------------
 *    
 * attaches an extension to a filename
 * 
 */
 
char *
attachext (void *space, char *str, Uint l, char *ext, Uint m)
{
    char *new;

	new = ALLOCMEMORY(space, NULL, char, l+m+1);
	strncpy(new, str, l);
	new[l]='\0';
	strncat(new, ext, m);
	
	return new;
}
/*-------------------------------- attachpath --------------------------------
 *    
 * attach a path (or any other string) to a filename (or any other string) 
 * and extensions (or any other string) at the end.
 * 
 */
 
char*
attachpath (void *space, char *str, Uint l, char *path, Uint m,  
				char *ext, Uint n)
{
	char *new;
	
	new = ALLOCMEMORY(space, NULL, char, (l+m+n+1));
	strncpy(new, path, m);
	new[m] = '\0';
	strncat(new, str, l);
	strncat(new, ext, n);
	
	return new;
}

int checkmd5(unsigned char *a, unsigned char *b) { 
  return (strncmp((char*)a, (char*)b, 16));
}
