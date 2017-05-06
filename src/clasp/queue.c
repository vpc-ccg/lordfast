/*
   This file is part of gdub.
   (C) 2006 Steve Hoffmann 

   gdub is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation; either version 2, or (at your
   option) any later version.

   gdub is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with gdub; see the file COPYING.  If not, write to the
   Free Software Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.	

*/
/**
 * queue.c
 * implementation of a simple queue for int
 *
 * @author Steve Hoffmann
 * @email steve@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Mon Oct 13 14:13:08 CEST 2008
 */

/*
 * SVN
 * Revision of last commit: $Rev: 72 $
 * Author: $Author: steve $
 * Date: $Date: 2008-10-28 18:14:42 +0100 (Tue, 28 Oct 2008) $
 * Id: $Id: queue.c 72 2008-10-28 17:14:42Z steve $
 * Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/queue.c $
 */

#include <stdlib.h>
#include <stdio.h>
#include "debug.h"
#include "queue.h"

/*----------------------------- bl_queueInit -----------------------------------
 *    
 * @brief 	init queue
 * @author 	Steve Hoffmann
 *   
 */
void bl_queueInit(Queue *q, int allocelem) {

  if (allocelem <= 0) {
    DBG("queue.c: Attempt to initialize a queue of size %d. Exit forced.\n",
	allocelem);
    exit(-1);
  } 

  q->queuespace = (Queueelem *) malloc(sizeof(Queueelem) * allocelem);  
  if (q->queuespace == NULL){
    DBG("queue.c: Memory allocation failed. Exit forced.\n", NULL);
    exit(-1);
  }
  q->allocelem = allocelem;
  q->numofelem = 0;
  q->enqueueindex = 0;
  q->dequeueindex = 0;
}

/*--------------------------- bl_queueDestruct ---------------------------------
 *    
 * @brief 	destruct queue
 * @author 	Steve Hoffmann
 *   
 */
void bl_queueDestruct(Queue *q) {

  free(q->queuespace);

  q->enqueueindex = 0;
  q->dequeueindex = 0;
  q->allocelem = 0;
  q->numofelem = 0;
}

/*---------------------------- bl_queueIsEmpty ---------------------------------
 *    
 * @brief 	returns if the queue is empty
 * @author 	Steve Hoffmann
 *   
 */
BOOL bl_queueIsEmpty(Queue *q) {
  return (q->numofelem == 0);
}

/*---------------------------- bl_queueEnqueue ---------------------------------
 *    
 * @brief 	enqueues elements at the back of the queue
 * @author 	Steve Hoffmann
 *   
 */
void bl_queueEnqueue(Queue *q, Queueelem elem) {

  if(q->numofelem == q->allocelem) {
    bl_queueResize(q);
  }

  q->queuespace[q->enqueueindex] = elem;
  q->numofelem++;

  /*implements circular datastructure*/
  if (q->enqueueindex == q->allocelem-1) {
    q->enqueueindex = 0;
  } else {
    q->enqueueindex++;
  }
}

/*---------------------------- bl_queueDequeue ---------------------------------
 *    
 * @brief 	dequeues element from the front of the queue
 * @author 	Steve Hoffmann
 *   
 */
Queueelem bl_queueDequeue(Queue *q) {

  Queueelem elem;		

  if(bl_queueIsEmpty(q)) {
    return 0;
  }	

  elem = q->queuespace[q->dequeueindex];
  q->numofelem--;

  /*implements circular data structure*/
  if(q->dequeueindex == q->allocelem - 1) {
    q->dequeueindex = 0;
  } else {
    q->dequeueindex++;
  }
  return elem;
}

/*---------------------------- bl_queueResize ----------------------------------
 *    
 * @brief 	expands the size of the queue to the double
 * @author 	Steve Hoffmann
 *   
 */
void bl_queueResize(Queue *q) {

  Queueelem *src;
  Queueelem *dest;
  void *ptr;

  /* resize queue to double */
  q->queuespace = (Queueelem *) realloc(q->queuespace, 
					sizeof(Queueelem) * (q->allocelem * 2));
  if (q->queuespace == NULL){
    DBG("queue.c: Memory reallocation failed. Exit forced.\n", NULL);
      exit(-1);
  }
  if (q->dequeueindex >= q->enqueueindex) {

    /* ptr arithmetics to move queue elements */
    src  = &q->queuespace[q->dequeueindex];
    dest = &q->queuespace[q->allocelem + q->dequeueindex]; 	

    ptr = memmove(dest, src,((q->allocelem)-q->dequeueindex)*sizeof(Queueelem));
    q->dequeueindex = (q->dequeueindex + q->allocelem);
  }

  q->allocelem *= 2;
}

/*------------------------------ bl_queueShow ----------------------------------
 *    
 * @brief 	prints the queue
 * @author 	Steve Hoffmann
 *   
 */
void bl_queueShow(Queue *q) {
  int i;
  Queueelem elem;

  printf("[");	
  for(i = 0; i < q->allocelem; i++) {
    elem = q->queuespace[i];
    if (i != q->enqueueindex && i != q->dequeueindex)
      printf("%d", elem);
    if (i == q->enqueueindex)
      printf("%d*", elem);
    if(i == q->dequeueindex)
      printf("%d^", elem);
    if(i+1 != q->allocelem)
      printf(",");
  }	
  printf("]\n");

}


/*------------------------------ bl_queueSize ----------------------------------
 *    
 * @brief 	returns number of elements in the queue
 * @author 	Steve Hoffmann
 *   
 */
Uint bl_queueSize(Queue *q){
  return q->numofelem;
}
