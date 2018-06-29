/**
 * container.c
 * implementation of a simple container for objects of defined size
 *
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Tue Oct 14 16:31:33 CEST 2008
 */

/*
 * SVN
 * Revision of last commit: $Rev: 89 $
 * Author: $Author: steve $
 * Date: $Date: 2008-11-24 14:53:55 +0100 (Mon, 24 Nov 2008) $
 * Id: $Id: container.c 89 2008-11-24 13:53:55Z steve $
 * Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/container.c $
 */

#include <stdlib.h>
#include <stdio.h>
#include "debug.h"
#include "basic-types.h"
#include "container.h"

/*--------------------------- bl_containerInit ---------------------------------
 *    
 * @brief 	init container
 * @author 	Christian Otto
 *   
 */
void bl_containerInit(Container *c, int allocelem, size_t sizeofelem){
  if (allocelem <= 0){
    DBG("container.c: Attempt to initialize a container of size %d.\
Exit forced.\n", allocelem);
    exit(-1);
  }
  if (sizeofelem <= 0){
    DBG("container.c: Attempt to initialize a container with sizeofelem %d.\
Exit forced.\n", sizeofelem);
    exit(-1);
  }
  c->contspace = malloc(allocelem * sizeofelem);
  if (c->contspace == NULL){
    DBG("container.c: Memory allocation failed. Exit forced.\n", NULL);
    exit(-1);
  }
  c->nextfree = 0;
  c->allocelem = allocelem;
  c->sizeofelem = sizeofelem;
}

/*-------------------------- bl_containerDestruct ------------------------------
 *    
 * @brief 	destruct container,
 *              remove method for elems as parameter possible
 * @author 	Christian Otto
 *   
 */
void bl_containerDestruct(Container *c, void (*rmv) (void*)){
  int i;
  char *p;
  if (rmv != NULL){
    p = (char *) c->contspace;
    for(i = 0; i < c->nextfree; i++){
      rmv(p + (i * c->sizeofelem));
    }
  }
  free(c->contspace);
  c->nextfree = 0;
  c->allocelem = 0;
  c->sizeofelem = 0;
}

/*--------------------------- bl_containerIsEmpty ------------------------------
 *    
 * @brief 	returns if the container is empty
 * @author 	Christian Otto
 *   
 */
BOOL bl_containerIsEmpty(Container *c){
  return (c->nextfree == 0);
}

/*--------------------------- bl_containerResize -------------------------------
 *    
 * @brief 	expands the size of the container by a given value
 * @author 	Christian Otto
 *   
 */
void bl_containerResize(Container *c, int inc){
  if (inc <= 0){
    DBG("container.c: Reallocation with %d senseless. Exit forced.\n", inc);
    exit(-1);
  }
  c->contspace = realloc(c->contspace, (c->allocelem + inc) * c->sizeofelem);
  if (c->contspace == NULL){
    DBG("container.c: Memory reallocation failed. Exit forced.\n", NULL);
    exit(-1);
  }
  c->allocelem += inc;
}

/*---------------------------- bl_containerAdd ---------------------------------
 *    
 * @brief 	adds element at the end of the container
 * @author 	Christian Otto
 *   
 */
void bl_containerAdd(Container *c, void *elem){
  char *p;
  if (c->nextfree == c->allocelem){
    bl_containerResize(c, BASEINC);
  }
  p = (char *) c->contspace;
  memmove(p + (c->nextfree * c->sizeofelem), elem, c->sizeofelem);
  c->nextfree++;
}

/*---------------------------- bl_containerGet ---------------------------------
 *    
 * @brief 	returns Nth object in the container
 *              with N = 0,..,numofelems - 1
 * @author 	Christian Otto
 *   
 */
void* bl_containerGet(Container *c, int n){
  char *p;
  if (bl_containerIsEmpty(c) || n < 0 || n >= c->nextfree){
    return NULL;
  }
  p = (char *) c->contspace;
  return (p + (n * c->sizeofelem));
}

/*--------------------------- bl_containerMerge --------------------------------
 *    
 * @brief 	merges two containers together
 * @author 	Christian Otto
 *   
 */
void bl_containerMerge(Container *c, Container *s){
  int size;
  char *p;
  if (c->sizeofelem != s->sizeofelem){
    DBG("container.c: Merge of containers with different data types failed.\
Exit forced.\n", NULL);
    exit(-1);
  }
  size = s->nextfree + c->nextfree;
  if (size >= c->allocelem){
    bl_containerResize(c, s->nextfree + BASEINC);
  }
  p = (char *) c->contspace;
  memmove(p + (c->nextfree * c->sizeofelem), s->contspace,
	  s->nextfree * s->sizeofelem);
  c->nextfree = size;
}

/*----------------------------- bl_containerSize -------------------------------
 *    
 * @brief 	returns number of elements in the container
 * @author 	Christian Otto
 *   
 */
Uint bl_containerSize(Container *c){
  return c->nextfree;
}
