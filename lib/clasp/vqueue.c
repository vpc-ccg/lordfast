/**
 * vqueue.c
 * implementation of a simple queue for objects of defined size
 *
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Mon Oct 13 14:13:08 CEST 2008
 */

/*
 * SVN
 * Revision of last commit: $Rev: 89 $
 * Author: $Author: steve $
 * Date: $Date: 2008-11-24 14:53:55 +0100 (Mon, 24 Nov 2008) $
 * Id: $Id: vqueue.c 89 2008-11-24 13:53:55Z steve $
 * Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/vqueue.c $
 */

#include <stdlib.h>
#include <stdio.h>
#include "debug.h"
#include "basic-types.h"
#include "vqueue.h"

/*----------------------------- bl_vqueueInit ----------------------------------
 *    
 * @brief 	init vqueue
 * @author 	Christian Otto
 *   
 */
void bl_vqueueInit(VQueue *q, int allocelem, size_t sizeofelem){
  if (allocelem <= 0){
    DBG("vqueue.c: Attempt to initialize a vqueue of size %d. Exit forced.\n",
	allocelem);
    exit(-1);
  }
  if (sizeofelem <= 0){
    DBG("vqueue.c: Attempt to initialize a vqueue with sizeofelem %d.\
Exit forced.\n", sizeofelem);
    exit(-1);
  }
  q->queuespace = malloc(allocelem * sizeofelem);
  if (q->queuespace == NULL){
    DBG("vqueue.c: Memory allocation failed. Exit forced.\n", NULL);
    exit(-1);
  }
  q->allocelem = allocelem;
  q->numofelem = 0;
  q->enqueueindex = 0;
  q->dequeueindex = 0;
  q->sizeofelem = sizeofelem;
}

/*--------------------------- bl_vqueueDestruct --------------------------------
 *    
 * @brief 	destruct vqueue,
 *              remove method for elems as parameter possible
 * @author 	Christian Otto
 *   
 */
void bl_vqueueDestruct(VQueue *q, void (*rmv)(void*)){
  int i;
  char *p;
  if (rmv != NULL){
    p = (char *) q->queuespace;
    for(i = 0; i < q->numofelem; i++){
      rmv(p + (q->dequeueindex * q->sizeofelem));
      if (q->dequeueindex == q->allocelem - 1){
	q->dequeueindex = 0;
      }
      else {
	q->dequeueindex++;
      }
    }
  }
  free(q->queuespace);
  q->allocelem = 0;
  q->numofelem = 0;
  q->enqueueindex = 0;
  q->dequeueindex = 0;
  q->sizeofelem = 0;
}

/*---------------------------- bl_vqueueIsEmpty --------------------------------
 *    
 * @brief 	returns if the vqueue is empty
 * @author 	Christian Otto
 *   
 */
BOOL bl_vqueueIsEmpty(VQueue *q){
  return (q->numofelem == 0);
}

/*---------------------------- bl_vqueueEnqueue --------------------------------
 *    
 * @brief 	enqueues elements at the back of the vqueue
 * @author 	Christian Otto
 *   
 */
void bl_vqueueEnqueue(VQueue *q, void *elem){
  char *p;
  if (q->numofelem == q->allocelem){
    bl_vqueueResize(q);
  }
  p = (char *) q->queuespace;
  memmove(p + (q->enqueueindex * q->sizeofelem), elem, q->sizeofelem);
  q->numofelem++;
  /* implements circular data structure */
  if (q->enqueueindex == q->allocelem - 1){
    q->enqueueindex = 0;
  }
  else {
    q->enqueueindex++;
  }
}

/*---------------------------- bl_vqueueDequeue --------------------------------
 *    
 * @brief 	dequeues element from the front of the vqueue as copy 
 *              and removes it from the vqueue
 * @author 	Christian Otto
 *   
 */
void* bl_vqueueDequeue(VQueue *q, void (*rmv)(void*)){
  char *p, *elem;
  if (bl_vqueueIsEmpty(q)){
    return NULL;
  }
  p = (char *) q->queuespace;
  elem = (char *) malloc(q->sizeofelem);
  memmove(elem, p + (q->dequeueindex * q->sizeofelem), q->sizeofelem);
  if (rmv != NULL){
    rmv(p + (q->dequeueindex * q->sizeofelem));
  }
  q->numofelem--;
  /* implements circular data structure */
  if (q->dequeueindex == q->allocelem - 1){
    q->dequeueindex = 0;
  }
  else {
    q->dequeueindex++;
  }
  return elem;
}

/*---------------------------- bl_vqueueFront ----------------------------------
 *    
 * @brief 	returns the front of the queue as pointer
 *              (next element that will be dequeued)
 * @author 	Christian Otto
 *   
 */
void* bl_vqueueFront(VQueue *q){
  char *p;
  if (bl_vqueueIsEmpty(q)){
    return NULL;
  }
  p = (char *) q->queuespace;
  return (p + (q->dequeueindex * q->sizeofelem));
}

/*---------------------------- bl_vqueueFrontN ---------------------------------
 *    
 * @brief 	returns Nth nearest object to the front of the vqueue
 *              with N = 0,..,numofelems - 1
 * @author 	Christian Otto
 *   
 */
void* bl_vqueueFrontN(VQueue *q, int n){
  char *p;
  int pos;
  if (bl_vqueueIsEmpty(q) || n < 0 || n >= q->numofelem){
    return NULL;
  }
  p = (char *) q->queuespace;
  pos = (q->dequeueindex + n) % q->allocelem;
  return (p + (pos * q->sizeofelem));
}

/*---------------------------- bl_vqueueResize ---------------------------------
 *    
 * @brief 	expands the size of the vqueue to the double
 * @author 	Christian Otto
 *   
 */
void bl_vqueueResize(VQueue *q){
  char *src, *dest;
  q->queuespace = realloc(q->queuespace, q->sizeofelem * (q->allocelem * 2));
  if (q->queuespace == NULL){
    DBG("vqueue.c: Memory reallocation failed. Exit forced.\n", NULL);
      exit(-1);
  }
  /* stretch the circle to line */
  if(q->dequeueindex >= q->enqueueindex){
    src = (char *) q->queuespace + (q->sizeofelem * q->dequeueindex);
    dest = (char *) q->queuespace +
      (q->sizeofelem * (q->allocelem + q->dequeueindex));
    memmove(dest, src, (q->allocelem - q->dequeueindex) * q->sizeofelem);
    q->dequeueindex = q->dequeueindex + q->allocelem;
  }
  q->allocelem *= 2;
}

/*------------------------------ bl_vqueueSize ---------------------------------
 *    
 * @brief 	returns number of elements in the vqueue
 * @author 	Christian Otto
 *   
 */
Uint bl_vqueueSize(VQueue *q){
  return q->numofelem;
}
