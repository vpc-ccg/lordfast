
/**
 * vebtree.c
 * an implementation of the van emde boas tree according to Johnson82 
 * over the priority domain [1, ..., N] with insert, delete, finding successor 
 * and predecessor in O(log log N) and finding max and min in constant time
 *
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Mon Jun  9 10:58:38 CEST 2008
 */

/*
 * SVN
 * Revision of last commit: $Rev: 114 $
 * Author: $Author: steve $
 * Date: $Date: 2010-04-19 15:59:41 +0200 (Mon, 19 Apr 2010) $
 * Id: $Id: vebtree.c 114 2010-04-19 13:59:41Z steve $
 * Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/vebtree.c $
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "basic-types.h"
#include "vebtree.h"
#include "container.h"
#include "debug.h"
#include "mathematics.h"

//BOOL debug = 0;

/*----------------------------- bl_vebtreeInit ---------------------------------
 *    
 * @brief 	initialize an vebtree which leafs correspond to data of any type
 *		within the priority domain [1, ..., allocelem]
 * @author 	Christian Otto
 *   
 */
void bl_vebtreeInit(VebTree *t, int allocelem, size_t size){
  bl_vebtreeInitRange(t, 1, allocelem, size);
}

/*--------------------------- bl_vebtreeInitRange ------------------------------
 *    
 * @brief 	initialize an vebtree which leafs correspond to data of any type
 *		within the priority domain [min, ..., allocelem + min - 1]
 * @author 	Christian Otto
 *   
 */
void bl_vebtreeInitRange(VebTree *t, int min, int allocelem, size_t size){
  Uint base, height;
  if (allocelem <= 0){
    DBG("vebtree.c: Attempt to initialize a vebtree of size %d.\
Exit forced.\n", allocelem);
    exit(-1);
  }
  if (min < 0){
    DBG("vebtree.c: Attempt to initialize a vebtree with negativ minimal \
priority %d. Exit forced.\n",min);
    exit(-1);
  }
  base = allocelem + 1;
  /* round up to the next highest power of 2 */
  base--;
  base |= base >> 1;
  base |= base >> 2;
  base |= base >> 4;
  base |= base >> 8;
  base |= base >> 16;
  base++;
  t->base = base;
  t->allocelem = allocelem;
  t->min = min;
	
  /* calc height of tree */
  height = 0;
  while (base != 0){
    base >>= 1;
    height++;
  }
  t->height = height;
  t->sizeofelem = size;
  /* set next of dummy element at beginning */
  t->next = -1;
  t->treespace = NULL;
  t->dataspace = NULL;
  #ifndef VEB_DEFERRED
  bl_vebtreeAlloc(t);
  #endif
}

void bl_vebtreeAlloc(VebTree *t){
  Uint i;
  Container *cont;	
  /* Allocations */
  t->treespace = (VebNode *) malloc(sizeof(VebNode) * (t->base + t->allocelem));
  t->dataspace = malloc(t->sizeofelem * t->allocelem);
  if (t->treespace == NULL || t->dataspace == NULL){
    DBG("vebtree.c: Memory allocation failed. Exit forced.\n", NULL);
    exit(-1);
  }	
	
  /* Initialization of nodes (not in O(log log n) but quite fast)*/
  memset(t->treespace, -1, sizeof(VebNode) * (t->base + t->allocelem));

  /* Initialization of leaf base */
  cont = bl_vebtreeBinSearch(t, t->base, 1);
  for (i = 0; i < bl_containerSize(cont); i++){
    int tmp = *((int *) bl_containerGet(cont, i));
    (t->treespace + tmp - 1)->left = 1;
    (t->treespace + tmp - 1)->right = 1;
    (t->treespace + tmp - 1)->lref = 1;
    (t->treespace + tmp - 1)->rref = 1;
  }
  (t->treespace)->left = t->base;
  (t->treespace)->right = t->base;
  bl_containerDestruct(cont, NULL);
  free(cont);
}

/*--------------------------- bl_vebtreeDestruct -------------------------------
 *    
 * @brief 	destruct an vebtree and removes all the data from the tree using
 *		a rmv method as parameter
 * @author 	Christian Otto
 *   
 */
void bl_vebtreeDestruct(VebTree *t, void (*rmv) (void *)){
  int i;
  char *p;
  p = t->dataspace;
  /* remove left data */
  if (t->dataspace != NULL && rmv != NULL){
    for (i = 0; i < t->allocelem; i++){
      if (bl_vebtreeIsActive(t, i + t->base)){
	rmv(p + (i * t->sizeofelem));
      }
    }
  }
  /* free data structures */
  if (t->dataspace != NULL){
    free(t->dataspace);
    t->dataspace = NULL;
  }
  if (t->treespace != NULL){
    free(t->treespace);
    t->treespace = NULL;
  }
  t->next = 0;
  t->sizeofelem = 0;
  t->height = 0;
  t->allocelem = 0;
  t->base = 0;
  t->min = 0;
}

/*--------------------------- bl_vebtreeIsEmpty --------------------------------
 *    
 * @brief 	returns if the vebtree is empty (in O(1)) which means that
 *		that there is at least one leaf active != t->base
 * @author 	Christian Otto
 *   
 */
BOOL bl_vebtreeIsEmpty(VebTree *t){
  return !(bl_vebtreeIsActive(t, 0) && ((t->treespace)->right != t->base));
}

/*--------------------------- bl_vebtreeIsActive -------------------------------
 *    
 * @brief       returns if the node value in the vebtree is active (in O(1)) 
 * @author      Christian Otto
 *   
 */
BOOL bl_vebtreeIsActive(VebTree *t, int value){
  if (t->treespace == NULL){
    return False;
  }
  return !((t->treespace + value)->left == -1 &&
	   (t->treespace + value)->right == -1 &&
	   ((t->treespace + value)->lref == -1 ||
	    (t->treespace + value)->lref == 0) &&
	   ((t->treespace + value)->rref == -1 ||
	    (t->treespace + value)->rref == 0));
}

/*---------------------------- bl_vebtreeGetNode -------------------------------
 *    
 * @brief 	returns a pointer to the node for a given node index (in O(1)),
 *		returns NULL pointer if the node does not exist
 * @author 	Christian Otto
 *   
 */
VebNode* bl_vebtreeGetNode(VebTree *t, int value){
  if (value < 0 || value >= (t->allocelem + t->base) || bl_vebtreeIsEmpty(t)){
    return NULL;
  }
  else {
    return (t->treespace + value);
  }
}

/*---------------------------- bl_vebtreeGetData -------------------------------
 *    
 * @brief 	return a pointer to the data for a given priority (in O(1)),
 *		returns NULL pointer if it was not successful
 * @author 	Christian Otto
 *   
 */
void* bl_vebtreeGetData(VebTree *t, int value){
  char *p;
  value -= (t->min - 1);
  if (value <= 0 || value > t->allocelem || bl_vebtreeIsEmpty(t) ||
      !bl_vebtreeIsActive(t, t->base + value - 1)){
    return NULL;
  }
  p = (char *) t->dataspace;
  return (p + (value - 1) * t->sizeofelem);
}

/*----------------------------- bl_vebtreeMin ----------------------------------
 *    
 * @brief 	returns the minimum priority of an element in the tree (in O(1))
 *		returns -1 if vebtree is empty
 * @author 	Christian Otto
 *   
 */
int bl_vebtreeMin(VebTree *t){
  if (bl_vebtreeIsEmpty(t)){
    return -1;
  }
  return t->next + t->min - 1;
}

/*----------------------------- bl_vebtreeMax ----------------------------------
 *    
 * @brief 	returns the maximum priority of an element in the tree (in O(1))
 *		returns -1 if vebtree is empty or only the dummy element 
 *		with priority 0 is in the tree
 * @author 	Christian Otto
 *   
 */
int bl_vebtreeMax(VebTree *t){
  int max;
  if (bl_vebtreeIsEmpty(t)){
    return -1;
  }
  max = ((t->treespace)->right - t->base);
  if (max == 0){
    return -1;
  }
  return ((t->treespace)->right - t->base + t->min - 1);
}

/*----------------------------- bl_vebtreePred ---------------------------------
 *    
 * @brief 	returns the index to the maximal element in the tree which is
 *		smaller than the given value (in O(log log n),works in O(1) if
 *		there is an element with the given value in the tree),
 *		returns -1 if there is no predecessor
 * @author 	Christian Otto
 *   
 */
int bl_vebtreePred(VebTree *t, int value){
  int f, lb, rb, ltop, rtop, lleaf, rleaf, pred;
  /* finding pred in O(1) for some conditions */
  if (bl_vebtreeIsEmpty(t)){
    return -1;
  }
  else if (value > bl_vebtreeMax(t)){
    return bl_vebtreeMax(t);
  }
  else if (value <= bl_vebtreeMin(t)){
    return -1;
  }
  value -= (t->min - 1);
  f = value + t->base;	
  /* find next right base (if node exists, this can be done in O(1)) */
  rb = bl_vebtreeLabelRuleR(t, f);
  rtop = RTOP(rb);
  rleaf = RLEAF(rtop);	
  /* if rleaf < f -> pred(f) = rleaf */
  if (rleaf < f){
    /* go to right base leaf */
    pred = RLEAF(rtop) - t->base;
    /* 
     * return right leaf if it has priority greater than 0,
     * otherwise it is just the dummy element
     */
    if (pred > 0){
      return pred + t->min - 1;
    }
    else {
      return -1;
    }
  }
  /* otherwise find next left base (if node exists, this can be done in O(1)) */
  lb = bl_vebtreeLabelRuleL(t, f);
  ltop = LTOP(lb);
  lleaf = LLEAF(ltop);
  /* if lleaf >= f -> pred(f) is rleaf of the corresponding right base */
  if (lleaf >= f){	
    /* go to corresponding right base (if ltop is not the root) */
    if (ltop != 1){
      rtop = ltop - 1;
    }
    else {
      rtop = ltop;
    }
    /* go to right base leaf */
    pred = RLEAF(rtop) - t->base;
    /* 
     * return right leaf if it has priority greater than 0,
     * otherwise it is just the dummy element
     */
    if (pred > 0){
      return pred + t->min - 1;
    }
    else {
      return -1;
    }
  }
  else {
    DBG("vebtree.c: Error in finding the predecessor (due to wrong labeling): \
%d %d %d %d %d. Exit forced.\n", f, ltop, rtop, lleaf, rleaf);
    exit(-1);
  }
}

/*----------------------------- bl_vebtreeSucc ---------------------------------
 *    
 * @brief 	returns the index to the minimal element in the tree which is 
 *		greater than the given value (in O(log log n),
 *		returns -1 if there is no successor 
 * @TODO        works in O(1) if an element with the given value exists
 * @author 	Christian Otto
 *   
 */
int bl_vebtreeSucc(VebTree *t, int value){
  int f, lb, rb, ltop, rtop, lleaf, rleaf, pred;
  /* finding pred in O(1) for some conditions */
  if (bl_vebtreeIsEmpty(t)){
    return -1;
  }
  else if (value >= bl_vebtreeMax(t)){
    return -1;
  }
  else if (value < bl_vebtreeMin(t)){
    return bl_vebtreeMin(t);
  }
  value -= (t->min - 1);
  f = value + t->base;	
  /* find next left base */
  lb = bl_vebtreeLabelRuleL(t, f);
  ltop = LTOP(lb);
  lleaf = LLEAF(ltop);	
  /* if lleaf > f -> pred(f) = lleaf */
  if (lleaf > f){
    /* go to left base leaf */
    pred = LLEAF(ltop) - t->base;
    /* 
     * returns right leaf if it has priority greater than 0,
     * otherwise it is just the dummy element
     */
    if (pred > 0){
      return pred + t->min - 1;
    }
    else {
      return -1;
    }
  }
  /* otherwise find next right base */
  rb = bl_vebtreeLabelRuleR(t, f);
  rtop = RTOP(rb);
  rleaf = RLEAF(rtop);
  /* if rleaf <= f -> pred(f) is lleaf of the corresponding left base */
  if (rleaf <= f){	
    /* go to corresponding left base (if rtop is not the root) */
    if (rtop != 1){
      ltop = rtop + 1;
    }
    else {
      ltop = rtop;
    }
    /* go to left base leaf */
    pred = LLEAF(ltop) - t->base;
    /*
     * returns right leaf if it has priority greater than 0,
     * otherwise it is just the dummy element
     */
    if (pred > 0){
      return pred + t->min - 1;
    }
    else {
      return -1;
    }
  }
  else {
    DBG("vebtree.c: Error in finding the successor (due to wrong labeling): \
%d %d %d %d %d. Exit forced.\n", f, ltop, rtop, lleaf, rleaf);
    exit(-1);
  }
}

/*----------------------------- bl_vebtreeInsert -------------------------------
 *    
 * @brief 	insert an element with a sort value and data into the tree,
 *		possibly deletes old data using the rmv method (in O(log log n))
 * @author 	Christian Otto
 *   
 */
void bl_vebtreeInsert(VebTree *t, int value, void *data, void (*rmv) (void *)){
  int d, e, b, c, f, g, h;//, neighbor;
  char *p;
  value -= (t->min - 1);
  f = value + t->base;
  if (t->treespace == NULL){
    bl_vebtreeAlloc(t);
  }
  p = t->dataspace;
  /* if value is out of range */
  if (value <= 0 || value > t->allocelem){
    return;
  }
  /* if element was already inserted -> replace old element */
  if (bl_vebtreeIsActive(t, f - 1)){
    if (rmv != NULL){
      if (bl_vebtreeGetData(t, value + t->min - 1) == NULL){
	DBG("vebtree: Error due to active node %d but no data in %d for value \
%d. Exit forced.", f-1, value+t->min - 1, value);
	exit(-1);
      }
      rmv(bl_vebtreeGetData(t, value + t->min - 1));
    }
    memmove(p + (value - 1) * t->sizeofelem, data, t->sizeofelem);
    return;
  }
  d = bl_vebtreeLabelRuleL(t, f);
  e = bl_vebtreeLabelRuleR(t, f);
  b = LTOP(d);
  c = RTOP(e);
  g = LLEAF(b);
  h = RLEAF(c);
  if (b <= 0 || c <= 0 || d <= 0 || e <= 0 || g <= 0 || h <= 0){
    DBG("vebtree.c: Label rule in insertion of %d was not successful. \
Exit forced.\n", value);
    exit(-1);
  }
  /* copy data */
  memmove(p + (value - 1) * t->sizeofelem, data, t->sizeofelem);
  if (f < g){
    /* set next */
    if (t->next == -1 || t->next > value){
      t->next = value;
    }
    /* restore proper labeling */
    bl_vebtreeIRestoreL(t, f, g, b);
  }
  else {
    /* set next */
    if (t->next == -1 || t->next > value){
      t->next = value;
    }
    /* restore proper labeling */
    bl_vebtreeIRestoreR(t, f, h, c);
  }
}

/*-------------------------- bl_vebtreeDelete ----------------------------------
 *    
 * @brief       delete an element from the vebtree and return a copy
 *	 	of the data for the given sort value (in O(log log n)),
 *		returns NULL pointer if it was not successful
 * @author      Christian Otto
 *   
 */
void* bl_vebtreeDelete(VebTree *t, int value, void (*rmv) (void *)){
  int f;
  char *p, *save;
  value -= (t->min - 1);
  if (bl_vebtreeIsEmpty(t)){
    return NULL;
  }
  p = t->dataspace;
  f = t->base + value;
  /* if val is out of range (element with priority 0 should never be deleted) */
  if (value <= 0 || value > t->allocelem){
    return NULL;
  }
  /* if element is not in the tree -> return NULL */
  if (!bl_vebtreeIsActive(t, f - 1)){
    return NULL;
  }
  /* copy data for return */
  save = (char *) malloc(t->sizeofelem);
  memmove(save, p + (value - 1) * t->sizeofelem, t->sizeofelem);
	
  /* if rmv parameter was set */
  if (rmv != NULL){
    rmv(p + (value - 1) * t->sizeofelem);
  }
  /* else reset the space with 0 */
  else {
    memset(p + (value - 1) * t->sizeofelem, 0, t->sizeofelem);
  }
  /* restore proper labeling for t */
  if (LTOP(f) < RTOP(f)){
    bl_vebtreeDRestoreL(t, f);
  }
  else {
    bl_vebtreeDRestoreR(t, f);
  }
  /* set prev and next */
  if (t->next == value){
    t->next = bl_vebtreeSucc(t, value + t->min - 1) - (t->min - 1);    
  }
  /* free if empty -> does not occur very frequently */
  #ifdef VEB_DEFERRED
  if (bl_vebtreeIsEmpty(t)){
    /* free data structures */
    free(t->dataspace);
    t->dataspace = NULL;
    free(t->treespace);
    t->treespace = NULL;
  }
  #endif
  return save;
}


/*---------------------------- bl_vebtreeBinSearch -----------------------------
 *    
 * @brief 	helper method for binary search in O(log log n), implemented 
 *		according to Johnson82, returns a list of node indexes (caution 
 *		index is from 1 to base+allocelem) on the up path from f to g
 * @author 	Christian Otto
 *   
 */
Container* bl_vebtreeBinSearch(VebTree *t, int f, int g){
  int i, j, b, ftmp;
  Container *cont;
  ftmp = f;	
  b = 0;
  j = 1;
  i = -1;
  cont = (Container *) malloc(sizeof(Container));
  bl_containerInit(cont, t->height, sizeof(int));
  bl_containerAdd(cont, &ftmp);
  while (ftmp > g){
    b += (1 << (j - 1));
    ftmp = MAX(1, (f >> b));
    j++;
    bl_containerAdd(cont, &ftmp);
  }
  i = j - 1;
  while (ftmp != g){
    if (ftmp < g){
      b -= (1 << (2 * i - j - 1));
    }
    else {
      b += (1 << (2 * i - j - 1));
    }
    /* unexpected behaviour */
    if (b <= 0){
      DBG("vebtree.c: vebtreeBinSearch: %d %d was not successful. \
Exit forced.\n", f, g);
      exit(-1);
    }
    ftmp = MAX(1, (f >> b));
    j++;
    bl_containerAdd(cont, &ftmp);
  }
  return cont;
}

/*------------------------- bl_vebtreeBinSearchRed -----------------------------
 *    
 * @brief 	helper method, binary search on the path from f to g 
 *		and stop if w is reached (in O(log log n))
 * @author 	Christian Otto
 *   
 */
Container* bl_vebtreeBinSearchRed(VebTree *t, int f, int g, int w){
  int i, j, b, ftmp, fsave;
  Container *cont;
  ftmp = f;
  fsave = f;
  b = 0;
  j = 1;
  i = -1;
  cont = (Container *) malloc(sizeof(Container));
  bl_containerInit(cont, t->height, sizeof(int));
  bl_containerAdd(cont, &ftmp);
  while (ftmp >= w){
    b += (1 << (j - 1));
    ftmp = MAX(1, (f >> b));
    j++;
    bl_containerAdd(cont, &ftmp);
    if (fsave <= w && ftmp  < w && ftmp < fsave){     /* different from paper */
      return cont;
    }
    fsave = ftmp;
  }
	
  i = j - 1;
  while (ftmp != g){
    if (ftmp < g){
      b -= (1 << (2 * i - j - 1));
    }
    else {
      b += (1 << (2 * i - j - 1));
    }
    /* unexpected behaviour */
    if (b <= 0){
      DBG("vebtree.c: vebtreeBinSearchRed: %d %d %d was not successful. \
Exit forced.\n", f, g, w);
      exit(-1);
    }		
    ftmp = MAX(1, (f >> b));
    j++;
    if (ftmp >= w){ /* has to be in it (but why?) */
      bl_containerAdd(cont, &ftmp);
    }
    if (fsave <= w && ftmp  < w && ftmp < fsave){     /* different from paper */
      return cont;
    }
    fsave = ftmp;
  }
	
  return cont;
}

/*------------------------- bl_vebtreeLabelRuleR -------------------------------
 *    
 * @brief 	helper method for insert,
 *		binary search for the first node from new leaf to the root
 *		where a right label is set (in O(log log n)),
 *              works in O(1) if f is active leaf
 * @author 	Christian Otto
 *   
 */
int bl_vebtreeLabelRuleR(VebTree *t, int f){
  int i, j, b, ftmp, flabel, frref;
  int max = -1;
  ftmp = f;
  flabel = (t->treespace + ftmp - 1)->right;
  frref = (t->treespace + ftmp - 1)->rref;
  b = 0;
  j = 1;
  i = -1;
  if (flabel > 0 && ftmp > max){
    max = ftmp;
  }
  while (flabel <= 0){
    b += (1 << (j - 1));
    ftmp = MAX(1, (f >> b));
    flabel = (t->treespace + ftmp - 1)->right;
    frref = (t->treespace + ftmp - 1)->rref;
    j++;
    if (flabel > 0 && ftmp > max){
      max = ftmp;
    }
  }
  i = j - 1;
  while (i < 0 || (2 * i - j - 1) >= 0){
    if (flabel >= 0){
      b -= (1 << (2 * i - j - 1));
    }
    else {
      b += (1 << (2 * i - j - 1));
    }
    /* unexpected behaviour */
    if (b <= 0){
      DBG("vebtree.c: vebtreeLabelRuleR: %d was not successful. \
Exit forced.\n", f);
      exit(-1);
    }
    ftmp = MAX(1, (f >> b));
    flabel = (t->treespace + ftmp - 1)->right;
    frref = (t->treespace + ftmp - 1)->rref;
    j++;	
    if (flabel > 0 && ftmp > max){
      max = ftmp;
    }		
  }
  return max;
}

/*--------------------------- bl_vebtreeLabelRuleL -----------------------------
 *    
 * @brief 	helper method for insert,
 *		binary search for the first node from new leaf to the root
 *		where a left label is set (in O(log log n)),
 *              works in O(1) if f is active leaf
 * @author 	Christian Otto
 *   
 */
int bl_vebtreeLabelRuleL(VebTree *t, int f){
  int i, j, b, ftmp, flabel, flref;
  int max = -1;
  ftmp = f;
  flabel = (t->treespace + ftmp - 1)->left;
  flref = (t->treespace + ftmp - 1)->lref;
  b = 0;
  j = 1;
  i = -1;
  if (flabel > 0 && ftmp > max){
    max = ftmp;
  }
  while (flabel <= 0){
    b += (1 << (j - 1));
    ftmp = MAX(1, (f >> b));
    flabel = (t->treespace + ftmp - 1)->left;
    flref = (t->treespace + ftmp - 1)->lref;
    j++;
    if (flabel > 0 && ftmp > max){
      max = ftmp;
    }
  }
  i = j - 1;
  while (i < 0 || (2 * i - j - 1) >= 0){
    if (flabel >= 0){
      b -= (1 << (2 * i - j - 1));
    }
    else {
      b += (1 << (2 * i - j - 1));
    }
    /* unexpected behaviour */
    if (b <= 0){
      DBG("vebtree.c: vebtreeLabelRuleL: %d was not successful. \
Exit forced.\n", f);
      exit(-1);
    }
    ftmp = MAX(1, (f >> b));
    flabel = (t->treespace + ftmp - 1)->left;
    flref = (t->treespace + ftmp - 1)->lref;
    j++;	
    if (flabel > 0 && ftmp > max){
      max = ftmp;
    }		
  }
  return max;
}

/*--------------------------- bl_vebtreeEqualRule ------------------------------
 *    
 * @brief 	helper method,
 *		binary search for the first node from the two leaves
 *		on the path to the root which is the same (in O(log log n))
 * @author 	Christian Otto
 *   
 */
int bl_vebtreeEqualRule(VebTree *t, int f, int g){
  int i, j, b, ftmp, gtmp, max = -1;
  ftmp = f;
  gtmp = g;	
  b = 0;
  j = 1;
  i = -1;
  if (ftmp == gtmp && ftmp > max){
    max = ftmp;
  }
  while (ftmp != gtmp){
    b += (1 << (j - 1));
    ftmp = MAX(1, (f >> b));
    gtmp = MAX(1, (g >> b));
    j++;
    if (ftmp == gtmp && ftmp > max){
      max = ftmp;
    }
  }
  i = j - 1;
  while (i < 0 || (2 * i - j - 1) >= 0){
    if (ftmp == gtmp){
      b -= (1 << (2 * i - j - 1));
    }
    else {
      b += (1 << (2 * i - j - 1));
    }
    /* unexpected behaviour */
    if (b <= 0){
      DBG("vebtree.c: vebtreeEqualRule: %d %d was not successful. \
Exit forced.\n", f, g);
      exit(-1);
    }
    ftmp = MAX(1, (f >> b));
    gtmp = MAX(1, (g >> b));
    j++;
    if (ftmp == gtmp && ftmp > max){
      max = ftmp;
    }
  }
  return max;
}

/*--------------------------- bl_vebtreeEqualRuleL -----------------------------
 *    
 * @brief 	helper method,
 *		binary search for the first node from the two leaves on the path
 *		to the root which is the same and which has a nonzero label on
 *              the left (in O(log log n))
 * @author 	Christian Otto
 *   
 */
int bl_vebtreeEqualRuleL(VebTree *t, int f, int g){
  int i, j, b, ftmp, gtmp, flabel, max = -1;
  ftmp = f;
  flabel = (t->treespace + ftmp - 1)->left;
  gtmp = g;	
  b = 0;
  j = 1;
  i = -1;
  if (ftmp == gtmp && flabel > 0 && ftmp > max){
    max = ftmp;
  }
  while (ftmp != gtmp || flabel <= 0){
    b += (1 << (j - 1));
    ftmp = MAX(1, (f >> b));
    flabel = (t->treespace + ftmp - 1)->left;
    gtmp = MAX(1, (g >> b));
    j++;
    if (ftmp == gtmp && flabel > 0 && ftmp > max){
      max = ftmp;
    }
  }
  i = j - 1;
  while (i < 0 || (2 * i - j - 1) >= 0){
    if (ftmp == gtmp){
      b -= (1 << (2 * i - j - 1));
    }
    else {
      b += (1 << (2 * i - j - 1));
    }
    /* unexpected behaviour */
    if (b <= 0){
      DBG("vebtree.c: vebtreeEqualRuleL: %d %d was not successful. \
Exit forced.\n", f, g);
      exit(-1);
    }
    ftmp = MAX(1, (f >> b));
    flabel = (t->treespace + ftmp - 1)->left;
    gtmp = MAX(1, (g >> b));
    j++;
    if (ftmp == gtmp && flabel > 0 && ftmp > max){
      max = ftmp;
    }
  }
  return max;
}

/*--------------------------- bl_vebtreeEqualRuleR -----------------------------
 *    
 * @brief 	helper method,
 *		binary search for the first node from the two leaves on the path
 *		to the root which is the same and which hasa nonzero label on
 *              the right (in O(log log n))
 * @author 	Christian Otto
 *   
 */
int bl_vebtreeEqualRuleR(VebTree *t, int f, int g){
  int i, j, b, ftmp, gtmp, flabel, max = -1;
  ftmp = f;
  flabel = (t->treespace + ftmp - 1)->right;
  gtmp = g;	
  b = 0;
  j = 1;
  i = -1;
  if (ftmp == gtmp && flabel > 0 && ftmp > max){
    max = ftmp;
  }
  while (ftmp != gtmp || flabel <= 0){
    b += (1 << (j - 1));
    ftmp = MAX(1, (f >> b));
    flabel = (t->treespace + ftmp - 1)->right;
    gtmp = MAX(1, (g >> b));
    j++;
    if (ftmp == gtmp && flabel > 0 && ftmp > max){
      max = ftmp;
    }
  }
  i = j - 1;
  while (i < 0 || (2 * i - j - 1) >= 0){
    if (ftmp == gtmp){
      b -= (1 << (2 * i - j - 1));
    }
    else {
      b += (1 << (2 * i - j - 1));
    }
    /* unexpected behaviour */
    if (b <= 0){
      DBG("vebtree.c: vebtreeEqualRuleR: %d %d was not successful. \
Exit forced.\n", f, g);
      exit(-1);
    }
    ftmp = MAX(1, (f >> b));
    flabel = (t->treespace + ftmp - 1)->right;
    gtmp = MAX(1, (g >> b));
    j++;
    if (ftmp == gtmp && flabel > 0 && ftmp > max){
      max = ftmp;
    }
  }
  return max;
}

/*--------------------------- bl_vebtreeIRestoreL ------------------------------
 *    
 * @brief 	helper method for insert,
 *		restore proper labeling after an insert of a new element
 *		which created a new path starting with a left fork
 *		(in O(log log n))
 * @author 	Christian Otto
 *   
 */
void bl_vebtreeIRestoreL(VebTree *t, int f, int g, int b){
  int i, p, r, s, w;
  Container *cont;
  p = bl_vebtreeEqualRule(t, f, g);
  r = 2 * p;
  s = r + 1;
  /* w in (g =>* b) | s */
  cont = bl_vebtreeBinSearchRed(t, g, b, s);
  for (i = 0; i < bl_containerSize(cont); i++){
    w = *((int *) bl_containerGet(cont, i));
    (t->treespace + w - 1)->lref--;
    if ((t->treespace + w - 1)->lref == 0){
      (t->treespace + w - 1)->left = -1;
    }
    else if (w >= s){  				      /* different from paper */
      (t->treespace + w - 1)->left = 0;
    }
  }
  bl_containerDestruct(cont, NULL);
  free(cont);
  /* w in (f =>* b) | r */
  cont = bl_vebtreeBinSearchRed(t, f, b, r);
  for (i = 0; i < bl_containerSize(cont); i++){
    w = *((int *) bl_containerGet(cont, i));
    if ((t->treespace + w - 1)->lref == -1){
      (t->treespace + w - 1)->lref = 1;
    }
    else {
      (t->treespace + w - 1)->lref++;
    }
    if (w >= b){
      (t->treespace + w - 1)->left = b;
    }
    else if ((t->treespace + w - 1)->left == -1){
      (t->treespace + w - 1)->left = 0;
    }
  }
  bl_containerDestruct(cont, NULL);
  free(cont);
  /* w in (g =>* s) */
  cont = bl_vebtreeBinSearch(t, g, s);
  for (i = 0; i < bl_containerSize(cont); i++){
    w = *((int *) bl_containerGet(cont, i));
    if ((t->treespace + w - 1)->lref == -1){
      (t->treespace + w - 1)->lref = 1;
    }
    else {
      (t->treespace + w - 1)->lref++;
    }
    if (w >= s){
      (t->treespace + w - 1)->left = s;
    }
    else if ((t->treespace + w - 1)->left == -1){
      (t->treespace + w - 1)->left = 0;
    }
  }
  bl_containerDestruct(cont, NULL);
  free(cont);
  /* w in (f =>* r) */
  // if (debug) printf("binsearch %d %d\n", f, r);
  cont = bl_vebtreeBinSearch(t, f, r);
  for (i = 0; i < bl_containerSize(cont); i++){
    w = *((int *) bl_containerGet(cont, i));
    if ((t->treespace + w - 1)->rref == -1){
      (t->treespace + w - 1)->rref = 1;
    }
    else {
      (t->treespace + w - 1)->rref++;
    }
    if (w >= r){
      (t->treespace + w - 1)->right = r;
    }
    else if ((t->treespace + w - 1)->right == -1){
      (t->treespace + w - 1)->right = 0;
    }
  }
  bl_containerDestruct(cont, NULL);
  free(cont);
	
  /* LLEAF(b) = f -> test if left of b is LLEAF */
  if ((t->treespace + b - 1)->left >= b){
    (t->treespace + b - 1)->left = f;
  }
  /* LLEAF(s) = g -> test if left of s is LLEAF */
  if ((t->treespace + s - 1)->left >= s){
    (t->treespace + s - 1)->left = g;
  }
  /* RLEAF(r) = f -> test if right of r is RLEAF */
  if ((t->treespace + r - 1)->right >= r){
    (t->treespace + r - 1)->right = f;
  }
}

/*--------------------------- bl_vebtreeIRestoreR ------------------------------
 *    
 * @brief 	helper method for insert,
 *		restore proper labeling after an insert of a new element
 *		which created a new path starting with a right fork,
 *		(in O(log log n))
 * @author 	Christian Otto
 *   
 */
void bl_vebtreeIRestoreR(VebTree *t, int f, int h, int c){
  int i, p, r, s, w;
  Container *cont;
  p = bl_vebtreeEqualRule(t, f, h);
  r = 2 * p;
  s = r + 1;
  /* w in (h =>* c) | r */
  cont = bl_vebtreeBinSearchRed(t, h, c, r);
  for (i = 0; i < bl_containerSize(cont); i++){
    w = *((int *) bl_containerGet(cont, i));
    (t->treespace + w - 1)->rref--;
    if ((t->treespace + w - 1)->rref == 0){
      (t->treespace + w - 1)->right = -1;
    }
    else if (w >= r){				      /* different from paper */
      (t->treespace + w - 1)->right = 0;
    }
  }
  bl_containerDestruct(cont, NULL);
  free(cont);
  /* w in (f =>* c) | s */
  cont = bl_vebtreeBinSearchRed(t, f, c, s);
  for (i = 0; i < bl_containerSize(cont); i++){
    w = *((int *) bl_containerGet(cont, i));
    if ((t->treespace + w - 1)->rref == -1){
      (t->treespace + w - 1)->rref = 1;
    }
    else {
      (t->treespace + w - 1)->rref++;
    }
    if (w >= c){
      (t->treespace + w - 1)->right = c;
    }
    else if ((t->treespace + w - 1)->right == -1){
      (t->treespace + w - 1)->right = 0;
    }
  }
  bl_containerDestruct(cont, NULL);
  free(cont);
  /* w in (h =>* r) */
  cont = bl_vebtreeBinSearch(t, h, r);
  for (i = 0; i < bl_containerSize(cont); i++){
    w = *((int *) bl_containerGet(cont, i));
    if ((t->treespace + w - 1)->rref == -1){
      (t->treespace + w - 1)->rref = 1;
    }
    else {
      (t->treespace + w - 1)->rref++;
    }
    if (w >= r){
      (t->treespace + w - 1)->right = r;
    }
    else if ((t->treespace + w - 1)->right == -1){
      (t->treespace + w - 1)->right = 0;
    }
  }
  bl_containerDestruct(cont, NULL);
  free(cont);
  /* w in (f =>* s) */
  cont = bl_vebtreeBinSearch(t, f, s);
  for (i = 0; i < bl_containerSize(cont); i++){
    w = *((int *) bl_containerGet(cont, i));
    if ((t->treespace + w - 1)->lref == -1){
      (t->treespace + w - 1)->lref = 1;
    }
    else {
      (t->treespace + w - 1)->lref++;
    }
    if (w >= s){
      (t->treespace + w - 1)->left = s;
    }
    else if ((t->treespace + w - 1)->left == -1){
      (t->treespace + w - 1)->left = 0;
    }
  }
  bl_containerDestruct(cont, NULL);
  free(cont);
	
  /* RLEAF(c) = f -> test if right of c is RLEAF */
  if ((t->treespace + c - 1)->right >= c){
    (t->treespace + c - 1)->right = f;
  }
  /* LLEAF(s) = f -> test if left of s is LLEAF */
  if ((t->treespace + s - 1)->left >= s){	
    (t->treespace + s - 1)->left = f;
  }
  /* RLEAF(r) = h -> test if right of r is RLEAF */
  if ((t->treespace + r - 1)->right >= r){
    (t->treespace + r - 1)->right = h;
  }
}

/*--------------------------- bl_vebtreeDRestoreL ------------------------------
 *    
 * @brief 	helper method for delete,
 *		restore proper labeling after an delete of a new element
 *		which deletes a path starting with a left fork,
 *		(in O(log log n))
 * @author 	Christian Otto
 *   
 */
void bl_vebtreeDRestoreL(VebTree *t, int f){
  int r, s, g, e, b, c, w, i;
  Container *cont;
  r = RTOP(f);
  s = r + 1;
  g = LLEAF(s);
  e = bl_vebtreeEqualRuleR(t, f, g);
  b = LTOP(f);
  c = RTOP(e);
  /* w in (f =>* r) */
  cont = bl_vebtreeBinSearch(t, f, r);
  for (i = 0; i < bl_containerSize(cont); i++){
    w = *((int *) bl_containerGet(cont, i));
    (t->treespace + w - 1)->rref--;
    if ((t->treespace + w - 1)->rref == 0){
      (t->treespace + w - 1)->right = -1;
    }
    else if (w >= r){				      /* different from paper */
      (t->treespace + w - 1)->right = 0;
    }
  }
  bl_containerDestruct(cont, NULL);
  free(cont);
  /* w in (g =>* s) */
  cont = bl_vebtreeBinSearch(t, g, s);
  for (i = 0; i < bl_containerSize(cont); i++){
    w = *((int *) bl_containerGet(cont, i));
    (t->treespace + w - 1)->lref--;
    if ((t->treespace + w - 1)->lref == 0){
      (t->treespace + w - 1)->left = -1;
    }
    else if (w >= s){				      /* different from paper */
      (t->treespace + w - 1)->left = 0;
    }
  }
  bl_containerDestruct(cont, NULL);
  free(cont);
  /* w in (f =>* b) | r */
  cont = bl_vebtreeBinSearchRed(t, f, b, r);
  for (i = 0; i < bl_containerSize(cont); i++){
    w = *((int *) bl_containerGet(cont, i));
    (t->treespace + w - 1)->lref--;
    if ((t->treespace + w - 1)->lref == 0){
      (t->treespace + w - 1)->left = -1;
    }
    else if (w >= r){				      /* different from paper */
      (t->treespace + w - 1)->left = 0;
    }
  }
  bl_containerDestruct(cont, NULL);
  free(cont);
  /* w in (g =>* b) | s */
  cont = bl_vebtreeBinSearchRed(t, g, b, s);
  for (i = 0; i < bl_containerSize(cont); i++){
    w = *((int *) bl_containerGet(cont, i));
    if ((t->treespace + w - 1)->lref == -1){
      (t->treespace + w - 1)->lref = 1;
    }
    else {
      (t->treespace + w - 1)->lref++;
    }
    if (w >= b){
      (t->treespace + w - 1)->left = b;
    }
    else if ((t->treespace + w - 1)->left == -1){
      (t->treespace + w - 1)->left = 0;
    }
  }
  bl_containerDestruct(cont, NULL);
  free(cont);	
  /* RLEAF(r) = -1 -> test if right of r is RLEAF */
  if ((t->treespace + r - 1)->right >= r){
    (t->treespace + r - 1)->right = -1;
  }
  /* LLEAF(s) = -1 -> test if left of s is LLEAF */
  if ((t->treespace + s - 1)->left >= s){
    (t->treespace + s - 1)->left = -1;
  }
  /* LLEAF(b) = g -> test if left of b is LLEAF */
  if ((t->treespace + b - 1)->left >= b){
    (t->treespace + b - 1)->left = g;
  }
}

/*--------------------------- bl_vebtreeDRestoreR ------------------------------
 *    
 * @brief 	helper method for delete,
 *		restore proper labeling after an delete of a new element
 *		which deletes a path starting with a right fork,
 *		(in O(log log n))
 * @author 	Christian Otto
 *   
 */
void bl_vebtreeDRestoreR(VebTree *t, int f){
  int r, s, h, e, b, c, w, i;
  Container *cont;
  r = LTOP(f);
  s = r - 1;
  h = RLEAF(s);
  e = bl_vebtreeEqualRuleL(t, f, h);
  b = RTOP(f);
  c = LTOP(e);
  /* w in (f =>* r) */
  cont = bl_vebtreeBinSearch(t, f, r);
  for (i = 0; i < bl_containerSize(cont); i++){
    w = *((int *) bl_containerGet(cont, i));
    (t->treespace + w - 1)->lref--;
    if ((t->treespace + w - 1)->lref == 0){
      (t->treespace + w - 1)->left = -1;
    }
    else if (w >= r){				      /* different from paper */
      (t->treespace + w - 1)->left = 0;
    }
  }
  bl_containerDestruct(cont, NULL);
  free(cont);
  /* w in (h =>* s) */
  cont = bl_vebtreeBinSearch(t, h, s);
  for (i = 0; i < bl_containerSize(cont); i++){
    w = *((int *) bl_containerGet(cont, i));
    (t->treespace + w - 1)->rref--;
    if ((t->treespace + w - 1)->rref == 0){
      (t->treespace + w - 1)->right = -1;
    }
    else if (w >= s){				      /* different from paper */
      (t->treespace + w - 1)->right = 0;
    }
  }
  bl_containerDestruct(cont, NULL);
  free(cont);
  /* w in (f =>* b) | r */
  cont = bl_vebtreeBinSearchRed(t, f, b, r);
  for (i = 0; i < bl_containerSize(cont); i++){
    w = *((int *) bl_containerGet(cont, i));
    (t->treespace + w - 1)->rref--;
    if ((t->treespace + w - 1)->rref == 0){
      (t->treespace + w - 1)->right = -1;
    }
    else if (w >= r){				      /* different from paper */
      (t->treespace + w - 1)->right = 0;
    }
  }
  bl_containerDestruct(cont, NULL);
  free(cont);
  /* w in (h =>* b) | s */
  cont = bl_vebtreeBinSearchRed(t, h, b, s);
  for (i = 0; i < bl_containerSize(cont); i++){
    w = *((int *) bl_containerGet(cont, i));
    if ((t->treespace + w - 1)->rref == -1){
      (t->treespace + w - 1)->rref = 1;
    }
    else {
      (t->treespace + w - 1)->rref++;
    }
    if (w >= b){
      (t->treespace + w - 1)->right = b;
    }
    else if ((t->treespace + w - 1)->right == -1){
      (t->treespace + w - 1)->right = 0;
    }
  }
  bl_containerDestruct(cont, NULL);
  free(cont);
  /* LLEAF(r) = -1 -> test if left of r is LLEAF */
  if ((t->treespace + r - 1)->left >= r){
    (t->treespace + r - 1)->left = -1;
  }
  /* RLEAF(s) = -1 -> test if right of s is RLEAF */
  if ((t->treespace + s - 1)->right >= s){
    (t->treespace + s - 1)->right = -1;
  }
  /* RLEAF(b) = h -> test if right of b is RLEAF */
  if ((t->treespace + b - 1)->right >= b){
    (t->treespace + b - 1)->right = h;
  }
}
