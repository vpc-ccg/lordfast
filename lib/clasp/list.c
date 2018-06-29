/**
 * list.c
 * implementation of a simple lineary linked list for object pointer
 *
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Wed Oct 15 11:39:42 CEST 2008
 */

/*
 * SVN
 * Revision of last commit: $Rev: 89 $
 * Author: $Author: steve $
 * Date: $Date: 2008-11-24 14:53:55 +0100 (Mon, 24 Nov 2008) $
 * Id: $Id: list.c 89 2008-11-24 13:53:55Z steve $
 * Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/list.c $
 */

#include <stdlib.h>
#include <stdio.h>
#include "debug.h"
#include "basic-types.h"
#include "list.h"

/*------------------------------ bl_listInit -----------------------------------
 *    
 * @brief 	init list
 * @author 	Christian Otto
 *   
 */
void bl_listInit(List *l, int allocelem, size_t sizeofelem){
  if (allocelem <= 0){
    DBG("list.c: Attempt to initialize a list of size %d.\
Exit forced.\n", allocelem);
    exit(-1);
  }
  if (sizeofelem <= 0){
    DBG("list.c: Attempt to initialize a list with sizeofelem %d.\
Exit forced.\n", sizeofelem);
    exit(-1);
  }
  l->nodes = (Listelem *) malloc(allocelem * sizeof(Listelem));  
  if (l->nodes == NULL){
    DBG("list.c: Memory allocation for nodes failed. Exit forced.\n", NULL);
    exit(-1);
  }
  l->data = malloc(allocelem * sizeofelem);  
  if (l->nodes == NULL){
    DBG("list.c: Memory allocation for data failed. Exit forced.\n", NULL);
    exit(-1);
  }
  l->first = -1;
  l->last = -1;
  l->nextfree = 0;
  l->numofelem = 0;
  l->allocelem = allocelem;
  l->sizeofelem = sizeofelem;
}

/*----------------------------- bl_listDestruct --------------------------------
 *    
 * @brief 	destruct list,
 *              remove method for elems as parameter possible
 * @author 	Christian Otto
 *   
 */
void bl_listDestruct(List *l, void (*rmv)(void*)){
  Uint cur;
  char *p;
  if (rmv != NULL){
    p = l->data;
    for (cur = l->first; cur != -1; cur = l->nodes[cur].next){
      rmv(p + (cur * l->sizeofelem));
    }
  }
  free(l->nodes);
  free(l->data);
  l->first = 0;
  l->last = 0;
  l->nextfree = 0;
  l->numofelem = 0;
  l->allocelem = 0;
  l->sizeofelem = 0;
}

/*----------------------------- bl_listIsEmpty ---------------------------------
 *    
 * @brief 	returns if the container is empty
 * @author 	Christian Otto
 *   
 */
BOOL bl_listIsEmpty(List *l){
  return (l->numofelem == 0);
}

/*----------------------------- bl_listInsert ----------------------------------
 *    
 * @brief 	adds element after an given element in the list
 *              (at beginning for cur == -1, at end for cur == l->last)
 * @author 	Christian Otto
 *   
 */
void bl_listInsert(List *l, int cur, void *elem){
  char *p;
  if (cur > l->allocelem || (cur < 0 && cur != -1)){
    return;
  }
  /* reallocation */
  if (l->nextfree >= l->allocelem){
    l->nodes = (Listelem *) realloc(l->nodes, sizeof(Listelem) *
				    (l->allocelem + BASEINC));
    if (l->nodes == NULL){
      DBG("list.c: Memory reallocation of nodes failed. Exit forced.\n", NULL);
      exit(-1);
    }
    l->data = realloc(l->data, l->sizeofelem * (l->allocelem + BASEINC));
    if (l->data == NULL){
      DBG("list.c: Memory reallocation of data failed. Exit forced.\n", NULL);
      exit(-1);
    }
    l->allocelem += BASEINC;
  }  
  p = (char *) l->data;
  /* insert data */
  memmove(p + (l->nextfree * l->sizeofelem), elem, l->sizeofelem);
  /* insert at begin (or in empty list) */
  if (cur == -1){
    l->nodes[l->nextfree].next = l->first;
    l->nodes[l->nextfree].prev = -1;
    if (l->first != -1){
      l->nodes[l->first].prev = l->nextfree;
    }
    else {      
      l->last = l->nextfree;
    }
    l->first = l->nextfree;
  }
  /* insert after elem cur */
  else {
    l->nodes[l->nextfree].prev = cur;
    l->nodes[l->nextfree].next = l->nodes[cur].next;
    /* new elem is last one */
    if (cur == l->last){
      l->last = l->nextfree;
    }
    /* otherwise */
    else {      
      l->nodes[l->nodes[l->nextfree].next].prev = l->nextfree;
    }
    l->nodes[cur].next = l->nextfree;
  }
  l->numofelem++;
  l->nextfree++;
}

/*------------------------------ bl_listUnlink ---------------------------------
 *    
 * @brief 	removes element from the list
 *              does not free
 * @author 	Christian Otto
 *   
 */
void* bl_listUnlink(List *l, Uint cur, void (*rmv)(void*)){
  char *p, *elem;
  p = (char *) l->data;
  if (cur > l->allocelem || cur < 0){
    return NULL;
  }
  elem = (char *) malloc(l->sizeofelem);
  memmove(elem, p + (cur * l->sizeofelem), l->sizeofelem);
  if (rmv != NULL){
    rmv(p + (cur * l->sizeofelem));
  }
  if (l->nodes[cur].prev != -1){
    l->nodes[l->nodes[cur].prev].next = l->nodes[cur].next;
  }
  else {
    l->first = l->nodes[cur].next;
  }
  if (l->nodes[cur].next != -1){
    l->nodes[l->nodes[cur].next].prev = l->nodes[cur].prev;
  }
  else {
    l->last = l->nodes[cur].prev;
  }
  l->nodes[cur].prev = -1;
  l->nodes[cur].next = -1;
  l->numofelem--;
  return elem;
}

/*------------------------------ bl_listSweep ----------------------------------
 *    
 * @brief 	cleans the list of all unlinked elements,
 *              implicitly sorts the nodes
 * @author 	Christian Otto
 *   
 */
void bl_listSweep(List *l){
  Uint cur, last = 0;
  Listelem *bufnodes;
  char *bufdata, *p;
  p = (char *) l->data;
  bufnodes = (Listelem *) malloc(sizeof(Listelem) * (l->numofelem + BASEINC));
  if (bufnodes == NULL){
    DBG("list.c: Memory allocation for nodes in sweep failed. Exit forced.\n",
	NULL);
    exit(-1);
  }
  bufdata = (char *) malloc(l->sizeofelem * (l->numofelem + BASEINC));
  if (bufdata == NULL){
    DBG("list.c: Memory allocation for data in sweep failed. Exit forced.\n",
	NULL);
    exit(-1);
  }
  for(cur = l->first; cur != -1; cur = l->nodes[cur].next){
    bufnodes[last].prev = last - 1;
    if (l->nodes[cur].next != -1){
      bufnodes[last].next = last + 1;
    }
    else {
      bufnodes[last].next = -1;
    }
    memmove(bufdata + (last * l->sizeofelem), p + (cur * l->sizeofelem),
	    l->sizeofelem);
  }
  free(l->nodes);
  free(l->data);
  l->nodes = bufnodes;
  l->data = bufdata;
  l->first = 0;
  l->last = l->numofelem - 1;
  l->allocelem = l->numofelem + BASEINC;
  l->nextfree = l->numofelem;
}

/*------------------------------ bl_listSize -----------------------------------
 *    
 * @brief 	returns number of elements in the list
 * @author 	Christian Otto
 *   
 */
Uint bl_listSize(List *l){
  return l->numofelem;
}
