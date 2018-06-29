/**
 * vstack.c
 * implementation of a simple stack for objects of defined size
 *
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Fri Oct 10 11:37:36 CEST 2008
 */

/*
 * SVN
 * Revision of last commit: $Rev: 89 $
 * Author: $Author: steve $
 * Date: $Date: 2008-11-24 14:53:55 +0100 (Mon, 24 Nov 2008) $
 * Id: $Id: vstack.c 89 2008-11-24 13:53:55Z steve $
 * Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/vstack.c $
 */

#include <stdio.h>
#include <stdlib.h>
#include "basic-types.h"
#include "debug.h"
#include "vstack.h"

/*----------------------------- bl_vstackInit ----------------------------------
 *    
 * @brief 	init vstack
 * @author 	Christian Otto
 *   
 */
void bl_vstackInit(VStack *s, int allocelem, size_t sizeofelem){
  if (allocelem <= 0){
    DBG("vstack.c: Attempt to initialize a vstack of size %d. Exit forced.\n",
	allocelem);
    exit(-1);
  }
  if (sizeofelem <= 0){
    DBG("vstack.c: Attempt to initialize a vstack with sizeofelem %d.\
Exit forced.\n", sizeofelem);
    exit(-1);
  }
  s->stackspace = malloc(allocelem * sizeofelem);
  if (s->stackspace == NULL){
    DBG("vstack.c: Memory allocation failed. Exit forced.\n", NULL);
    exit(-1);
  }
  s->allocelem = allocelem;
  s->top = -1;
  s->sizeofelem = sizeofelem;
}

/*--------------------------- bl_vstackDestruct --------------------------------
 *    
 * @brief 	destruct vstack,
 *              remove method for elems as parameter possible
 * @author 	Christian Otto
 *   
 */
void bl_vstackDestruct(VStack *s, void (*rmv)(void*)){
  int i;
  char *p;
  if (rmv != NULL){
    p = (char *) s->stackspace;
    for(i = 0; i <= s->top; i++){
      rmv(p + (i * s->sizeofelem));
    }
  }
  free(s->stackspace);
  s->allocelem = 0;
  s->top = 0;
  s->sizeofelem = 0;  
}

/*---------------------------- bl_vstackIsEmpty --------------------------------
 *    
 * @brief 	returns if the vstack is empty
 * @author 	Christian Otto
 *   
 */
BOOL bl_vstackIsEmpty(VStack *s){
  return (s->top < 0);
}

/*----------------------------- bl_vstackPush ----------------------------------
 *    
 * @brief 	pushs elements on the top of the vstack
 * @author 	Christian Otto
 *   
 */
void bl_vstackPush(VStack *s, void *elem){
  char *p;
  if (s->top >= s->allocelem - 1){
    s->stackspace = realloc(s->stackspace,
			    s->sizeofelem * (s->allocelem + BASEINC));
    if (s->stackspace == NULL || BASEINC <= 0){
      DBG("vstack.c: Memory reallocation failed. Exit forced.\n", NULL);
      exit(-1);
    }
    s->allocelem += BASEINC;
  }
  s->top++;
  p = (char *) s->stackspace;
  memmove(p + (s->top * s->sizeofelem), elem, s->sizeofelem);
}

/*------------------------------ bl_vstackTop ----------------------------------
 *    
 * @brief 	returns top of the vstack as pointer
 * @author 	Christian Otto
 *   
 */
void* bl_vstackTop(VStack *s){
  char *p;
  if (bl_vstackIsEmpty(s)){
    return NULL;
  }
  p = (char *) s->stackspace;
  return (p + (s->top * s->sizeofelem));
}

/*------------------------------ bl_vstackTopN ---------------------------------
 *    
 * @brief 	returns Nth highest object of the vstack
 *              with N = 0,..,numofelems - 1
 * @author 	Christian Otto
 *   
 */
void* bl_vstackTopN(VStack *s, int n){
  char *p;
  if (bl_vstackIsEmpty(s) || n < 0 || n > s->top){
    return NULL;
  }
  p = (char *) s->stackspace;
  return (p + (s->top - n) * s->sizeofelem);
}

/*------------------------------ bl_vstackPop ----------------------------------
 *    
 * @brief 	pops the top of the vstack as copy 
 *              and removes it from the vstack
 * @author 	Christian Otto
 *   
 */
void* bl_vstackPop(VStack *s, void (*rmv)(void*)){
  char *p, *elem;
  if (bl_vstackIsEmpty(s)){
    return NULL;
  }
  p = (char *) s->stackspace;
  elem = (char *) malloc(s->sizeofelem);
  memmove(elem, p + (s->top * s->sizeofelem), s->sizeofelem);
  if (rmv != NULL){
    rmv(p + (s->top * s->sizeofelem));
  }
  s->top--;
  return elem;
}

/*------------------------------ bl_vstackSize ---------------------------------
 *    
 * @brief 	returns number of elements on the vstack
 * @author 	Christian Otto
 *   
 */
Uint bl_vstackSize(VStack *s){
  return (s->top + 1);
}
