
/**
 * bintree.c
 * a binary search tree with size $n=2^k$
 * over an search space 1...n (within the
 * bintree the priorities are from 0...n-1
 * due to using arrays starting from 0)
 * implemented as a 'static' field
 *
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Wed Jun  4 13:58:52 CEST 2008
 */ 

/* 
 * SVN
 * Revision of last commit: $Rev: 114 $
 * Author: $Author: steve $
 * Date: $Date: 2010-04-19 15:59:41 +0200 (Mon, 19 Apr 2010) $
 * Id: $Id: bintree.c 114 2010-04-19 13:59:41Z steve $
 * Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/bintree.c $
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "basic-types.h"
#include "mathematics.h"
#include "bintree.h"
#include "debug.h"

/*-------------------------- bl_bintreeInit ------------------------------------
 *    
 * @brief       initialize a binary tree with any data which is sorted from 
 *              $1 to 2^k$ with $2^k >= allocelem$ and sets father relations
 * @author      Christian Otto
 *   
 */
void bl_bintreeInit(BinTree *t, int allocelem, size_t size){
  bl_bintreeInitRange(t, 1, allocelem, size);
}

/*-------------------------- bl_bintreeInit ------------------------------------
 *    
 * @brief       initialize a binary tree with any data which is sorted from 
 *              $min to 2^k$ with $2^k >= allocelem - min + 1$ and sets father relations
 * @author      Christian Otto
 *   
 */
void bl_bintreeInitRange(BinTree *t, int min, int allocelem, size_t size){
  int base, next;
  if (allocelem <= 0){
    DBG("bintree.c: Attempt to initialize a binary tree of size %d. \
Exit forced.\n", allocelem);
    exit(-1);
  }
  /* round up to the next highest power of 2 */
  next = allocelem;
  next--;
  next |= next >> 1;
  next |= next >> 2;
  next |= next >> 4;
  next |= next >> 8;
  next |= next >> 16;
  next++;
  /* calc height of tree */
  t->allocelem = 0;
  base = next;
  while (base != 0 && t->allocelem < allocelem){
    base >>= 1;
    t->allocelem += base;
  }
  t->min = min;
  t->allocelem = next;
  t->sizeofelem = size;
  t->treespace = NULL;
  t->father = NULL;
  t->left = NULL;
  t->right = NULL;
  t->used = NULL;
  t->numofelems = -1;
  #ifndef BIN_DEFERRED
  bl_bintreeAlloc(t);
  #endif
}

void bl_bintreeAlloc(BinTree *t){
  /* Allocations */
  t->treespace = malloc(t->sizeofelem * t->allocelem);
  t->father = (int *) malloc(sizeof(int) * t->allocelem);
  t->left = (int *) malloc(sizeof(int) * t->allocelem);
  t->right = (int *) malloc(sizeof(int) * t->allocelem);
  t->used = (BOOL *) malloc(sizeof(BOOL) * (t->allocelem + 1));
  t->numofelems = 0;
  if (t->treespace == NULL ||
      t->father == NULL ||
      t->left == NULL ||
      t->right == NULL ||
      t->used == NULL){
    DBG("bintree.c: Memory allocation failed. Exit forced.\n", NULL);
    exit(-1);
  }	
	
  /* Initializations */
  t->root = MAX(0, t->allocelem/2 - 1);
  t->father[t->root]=-1;
  bl_bintreeSetFather(t, t->root, t->allocelem/2);
  if (t->allocelem > 1){
    t->father[t->allocelem-1]=t->allocelem-2;
  }
  memset(t->treespace, 0, t->allocelem * t->sizeofelem);
  memset(t->left, (int)-1, t->allocelem * sizeof(int));
  memset(t->right, (int)-1, t->allocelem * sizeof(int));
  memset(t->used, (BOOL) 0, t->allocelem);
  t->used[t->allocelem]='\0';
}

/*------------------------ bl_bintreeSetFather ---------------------------------
 *    
 * @brief       helper method for bintreeInit
 * @author      Christian Otto
 *   
 */
void bl_bintreeSetFather(BinTree *t, int current, int level){
  if (level > 1){
    level /= 2;
    t->father[current-level]=current;
    bl_bintreeSetFather(t, current-level, level);
    if (current+level < t->allocelem){
      t->father[current+level]=current;
      bl_bintreeSetFather(t, current+level, level);
    }
  }
}

/*------------------------ bl_bintreeDestruct ----------------------------------
 *    
 * @brief       destruct a binary tree and removes all the data from the tree
 *	 	using a rmv method as parameter		  
 * @author      Christian Otto
 *   
 */
void bl_bintreeDestruct(BinTree *t, void (*rmv)(void *)){
  int i;
  char *p;
  p = (char *) t->treespace;
  if (t->treespace != NULL && rmv != NULL){
    for (i = 0; i < t->allocelem; i++){
      if (t->used[i] == 1){
	rmv(p + (i * t->sizeofelem));
      }
    }
  }
  if (t->treespace != NULL){
    free(t->treespace);
    t->treespace = NULL;
  }
  if (t->father != NULL){
    free(t->father);
    t->father = NULL;
  }
  if (t->left != NULL){
    free(t->left);
    t->left = NULL;
  }
  if (t->right != NULL){
    free(t->right);
    t->right = NULL;
  }
  if (t->used != NULL){
    free(t->used);
    t->used = NULL;
  }
  t->allocelem = 0;
  t->sizeofelem = 0;
  t->numofelems = 0;
  t->min = 0;
}

/*------------------------ bl_bintreeResize ------------------------------
 *    
 * @brief       resizes the binary tree to double of the size if an element
 *		should be inserted which greater equal than allocelem
 * @author      Christian Otto
 *   
 */
void bl_bintreeResize(BinTree *t){
  int root, level;
  char *p;
  DBG("no resize possible\n", NULL);
  exit(-1);
  /* Allocations */
  t->treespace = realloc(t->treespace, t->sizeofelem * t->allocelem * 2);
  t->father = (int *) realloc(t->father, sizeof(int) * t->allocelem * 2);
  t->left = (int *) realloc(t->left, sizeof(int) * t->allocelem * 2);
  t->right = (int *) realloc(t->right, sizeof(int) * t->allocelem * 2);
  t->used = (BOOL *) realloc(t->used, sizeof(BOOL) * t->allocelem * 2);
  if (t->treespace == NULL || t->father == NULL || t->left == NULL ||
      t->right == NULL || t->used == NULL){
    DBG("bintree.c: Memory reallocation failed. Exit forced.\n", NULL);
    exit(-1);
  }

  /* Initializations */
  root = t->allocelem - 1;
  if (t->used[root] != 0){
    /* delete root without deleting the data */
    free(bl_bintreeDelete(t, root, NULL));
    t->used[root] = 1;
  }
  t->father[root] = -1;
	
  level = (root+1)/2;
  if (level > 1){
    t->father[root-level] = root;
    t->father[root+level] = root;
    bl_bintreeSetFather(t, root+level, level);
    t->father[2*t->allocelem-1] = 2 * t->allocelem - 2;

    if ((t->used[root-level] != 0) ||
	(t->left[root-level] != -1) ||
	(t->right[root-level] != -1)){
      t->left[root] = root-level;
    }
  }
  p = t->treespace;
  memset(p + (t->allocelem * t->sizeofelem), 0,
	 t->allocelem * t->sizeofelem);
  memset(t->left + t->allocelem, -1, t->allocelem * sizeof(int));
  memset(t->right + t->allocelem, -1, t->allocelem * sizeof(int));
  memset(t->used + t->allocelem, (BOOL) 0, t->allocelem);
  t->allocelem *= 2;
}

/*------------------------- bl_bintreeIsEmpty ----------------------------------
 *    
 * @brief       returns if the bintree is empty (in O(1))	  
 * @author      Christian Otto
 *   
 */
BOOL bl_bintreeIsEmpty(BinTree *t){
  if (bl_bintreeSize(t) == -1){
    return True;
  }
  return ((t->used[t->root] == 0) &&
	  (t->left[t->root] == -1) &&
	  (t->right[t->root] == -1));
}

/*-------------------------- bl_bintreeInsert ----------------------------------
 *    
 * @brief       insert an element with a sort value and data into the tree
 *		could lead to resize of the bintree (in O(log n)),
 *		possibly deletes old data using the rmv method
 * @author      Christian Otto
 *   
 */
void bl_bintreeInsert(BinTree *t, int value, void *data, void (*rmv) (void *)){
  int i, tmp;
  char *p;
  /* alloc first if deferred init method is used */
  if (bl_bintreeSize(t) == -1){
    bl_bintreeAlloc(t);
  }
  /* reducing the value by min to fit the array representation */
  value -= t->min;
  if (value < 0){
    return;
  }
  else if (value >= t->allocelem){
    while (value >= t->allocelem){
      bl_bintreeResize(t);
    }
  }
  p = (char*) t->treespace;	
  /* delete old data */
  if (t->used[value] == 1){
    if (rmv != NULL){
      rmv(p + (value * t->sizeofelem));
    }
    t->numofelems--;
  }
  /* set search path to root */
  else {
    t->used[value] = 1;
    for(i = value; t->father[i] != -1; i = t->father[i]){
      tmp = t->father[i];		
      if (tmp < i){
	t->right[tmp] = i;
	if (t->left[tmp] != -1){
	  break;
	}
      }
      else {
	t->left[tmp] = i;
	if (t->right[tmp] != -1){
	  break;
	}
      }
    }
  }
  /* copy data */	
  memmove((p + (value * t->sizeofelem)), data, t->sizeofelem);
  t->numofelems++;
}

/*---------------------------- bl_bintreeGet -----------------------------------
 *    
 * @brief       returns a pointer to the data for a given sort value (in O(1))
 *		(returns NULL pointer if it was not successful)
 * @author Christian Otto
 *   
 */
void* bl_bintreeGet(BinTree *t, int value){
  /* reducing the value by min to fit the array representation */
  value -= t->min;
  if (value < 0 || value >= t->allocelem || bl_bintreeIsEmpty(t)){
    return NULL;
  }

  if (t->used[value] == 0){
    return NULL;
  }
  else {
    char *p = (char *) t->treespace;
    return p + (value * t->sizeofelem);
  }
}

/*---------------------------- bl_bintreeDelete --------------------------------
 *    
 * @brief       delete an element from the bintree and return a copy
 *	 	of the data for the given sort value (in O(log n))
 *		(returns NULL pointer if it was not successful)
 * @author      Christian Otto
 *   
 */
void* bl_bintreeDelete(BinTree *t, int value, void (*rmv)(void *)){
  /* reducing the value by min to fit the array representation */
  value -= t->min;
  if (value < 0 || value >= t->allocelem || bl_bintreeIsEmpty(t)){
    return NULL;
  }
	
  if (t->used[value] == 0){
    return NULL;
  }
  else {
    int node, val;
    char *tmp, *p;
    /* copy data for return */
    tmp = malloc(t->sizeofelem);
    p = (char *) t->treespace;
    memmove(tmp, p + (value * t->sizeofelem), t->sizeofelem);
    /* remove the data */
    if (rmv != NULL){
      rmv(p + (value * t->sizeofelem));
    }
    /* else reset the space with 0 */
    else {
      memset(p + (value  * t->sizeofelem), 0, t->sizeofelem);
    }
    t->used[value] = 0;
    /* clean the path to the root */
    if ((t->left[value] == -1) &&
	(t->right[value] == -1) &&
	(t->father[value] != -1)){
			
      node = t->father[value];
      val = value;
      while (node != -1){
	if (node < val){
	  t->right[node] = -1;
	  if (t->left[node] != -1 || t->used[node] != 0){
	    break;
	  }
	}
	else {
	  t->left[node] = -1;
	  if (t->right[node] != -1 || t->used[node] != 0){
	    break;
	  }
	}
	val = node;
	node = t->father[node];
      }
    }
    t->numofelems--;
    return tmp;
  }
}

/*---------------------------- bl_bintreeMin -----------------------------------
 *    
 * @brief       returns the index to the minimal element in the tree according
 *              to the sort value (in O(log n)), returns -1 if bintree is empty
 *		
 * @author Christian Otto
 *   
 */
int bl_bintreeMin(BinTree *t){
  if (bl_bintreeIsEmpty(t)){
    return -1;
  }
  int node;
  node = t->root;
  while (t->left[node] != -1 || t->used[node] != 0 || t->right[node] != -1){
    if (t->left[node] != -1){
      node = t->left[node];
    }
    else if (t->used[node] != 0){
      break;
    }
    else if (t->right[node] != -1){
      node = t->right[node];
    }
  }
  /* increasing the node by min to fit the priority representation */
  return node + t->min;
}

/*------------------------------ bl_bintreeMax ---------------------------------
 *    
 * @brief       returns the index to the maximal element in the tree according
 *		to the sort value (in O(log n)), returns -1 if bintree is empty
 * @author Christian Otto
 *   
 */
int bl_bintreeMax(BinTree *t){
  if (bl_bintreeIsEmpty(t)){
    return -1;
  }
  int node;
  node = t->root;
  while (t->right[node] != -1 || t->used[node] != 0 || t->left[node] != -1){
    if (t->right[node] != -1){
      node = t->right[node];
    }
    else if (t->used[node] != 0){
      break;
    }
    else if (t->left[node] != -1){
      node = t->left[node];
    }
  }
  /* increasing the node by min to fit the priority representation */
  return node + t->min;
}

/*------------------------------ bl_bintreePred --------------------------------
 *    
 * @brief       returns the index to the maximal element in the tree which is
 *		smaller than the given value (in O(log n)),
 *		returns -1 if there is no predecessor
 * @author      Christian Otto
 *   
 */
int bl_bintreePred(BinTree *t, int value){
  /* reducing the value by min to fit the array representation */
  value -= t->min;
  if (value < 0 || bl_bintreeIsEmpty(t)){
    return -1;
  }
  else if (value >= t->allocelem){
    return bl_bintreeMax(t);
  }
  int node, father;
  node = value;
  /* going up in the tree */
  if (t->left[node] == -1){
    while (node != -1){
      father = t->father[node];
      if ((node > father) &&
	  (father != -1) &&
	  (t->left[father] != -1 || t->used[father] != 0)){
				
	if (t->used[father] != 0){
	  /* increasing the father by min to fit the priority representation */
	  return father + t->min;
	}
	else {
	  node = t->left[father];
	}
	break;
      }
      node = father;
    }
    if (node == -1){
      return -1;
    }
  }
  else {
    node = t->left[node];
  }
  /* going down in the tree (most right path) */
  while (t->right[node] != -1 || t->used[node] != 0 || t->left[node] != -1){
    if (t->right[node] != -1){
      node = t->right[node];
    }
    else if (t->used[node] != 0){
      break;
    }
    else if (t->left[node] != -1){
      node = t->left[node];
    }
  }
  /* increasing the node by min to fit the priority representation */
  return node + t->min;
}

/*------------------------------ bl_bintreeSucc --------------------------------
 *    
 * @brief       returns the index to the minimal element in the tree which is
 *		greater than the given value (in O(log n)),
 *		returns -1 if there is no successor
 * @author      Christian Otto
 *   
 */
int bl_bintreeSucc(BinTree *t, int value){
  /* reducing the value by min to fit the array representation */
  value -= t->min;
  if (value >= t->allocelem || bl_bintreeIsEmpty(t)){
    return -1;
  }
  else if (value < 0){
    return bl_bintreeMin(t);
  }
  int node, father;
  node = value;
  /* going up in the tree */
  if (t->right[node] == -1){
    while (node != -1){
      father = t->father[node];
      if ((node < father) &&
	  (t->used[father] != 0 || t->right[father] != -1)){
				
	if (t->used[father] != 0){
	  /* increasing the father by min to fit the priority representation */
	  return father + t->min;	
	}
	else {
	  node = t->right[father];
	}
	break;
      }
      node = father;
    }
    if (node == -1){
      return -1;
    }
  }
  else {
    node = t->right[node];
  }
  /* going down in the tree (most left path) */
  while (t->left[node] != -1 || t->used[node] != 0 || t->right[node] != -1){
    if (t->left[node] != -1){
      node = t->left[node];
    }
    else if (t->used[node] != 0){
      break;
    }
    else if (t->right[node] != -1){
      node = t->right[node];
    }
  }
  /* increasing the node by min to fit the priority representation */
  return node + t->min;
}

/*----------------------------- bl_bintreeSize ---------------------------------
 *    
 * @brief       returns number of elements in the bintree
 * @author      Christian Otto
 *   
 */
Uint bl_bintreeSize(BinTree *t){
  return t->numofelems;
}
