
/**
 * rangetree.c
 * an implementation of an d-dimensional semi-dynamic range tree
 * with O(n log^(d-1) n) in space, construction in O(n log^(d-1) n) in time
 * query time O(log^(d-1) n),
 * all data points must be pairwise disjunct in each dimension
 *
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Wed Aug 13 13:28:23 CEST 2008
 */

/*
 * SVN
 * Revision of last commit: $Rev: 114 $
 * Author: $Author: steve $
 * Date: $Date: 2010-04-19 15:59:41 +0200 (Mon, 19 Apr 2010) $
 * Id: $Id: rangetree.c 114 2010-04-19 13:59:41Z steve $
 * Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/rangetree.c $
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "basic-types.h"
#include "vebtree.h"
#include "bintree.h"
#include "container.h"
#include "debug.h"
#include "rangetree.h"

/*---------------------------- bl_rangetreeInit --------------------------------
 *    
 * @brief       initialize a multidimensional range tree with an array at the
 *              lowest dimension, semi-dynamic due to static initialization but
 *              dynamic insert and delete at the one dimensional array structure
 * @parameter   pos - matrix with n columns and dim rows, pos[i][j] = object at
 *              j-th position (j = 0,...,n-1) in i-th sorted array
 *              (i = 0,...,dim-1)
 * @author      Christian Otto
 *   
 */
void bl_rangetreeInit(void *t, Uint dim, Uint n, Uint** pos, size_t sizeofelem){
  Uint i, j, *initfield, height = 0, level = n, exp = 1;
  double binom = 1;
  Container points;
  rtpoint_t point, *tmppoint;
  RangeTree *rt = (RangeTree *) t;
  if (dim <= 1 || n < 1){
    DBG("rangetree.c: Attempt to initialize a range tree of dimension %d with \
%d elements. Exit forced.\n", dim, n);
    exit(-1);
  }
  if (sizeofelem <= 0){
    DBG("rangetree.c: Attempt to initialize a range tree with sizeofelem %d.\
Exit forced.\n", sizeofelem);
    exit(-1);
  }
  bl_containerInit(&points, n, sizeof(rtpoint_t));
  for(i = 0; i < n; i++){
    point.dim = dim;
    point.coord = (Uint*) malloc(dim * sizeof(Uint));
    bl_containerAdd(&points, &point);
  }
  for(i = 0; i < dim; i++){
    for(j = 0; j < n; j++){
      tmppoint = (rtpoint_t*) bl_containerGet(&points, pos[i][j]);
      /* object at pos[i][j] is in dimension i at position j */
      tmppoint->coord[i] = j;
    }
  }
  rt->space = (rtspace_t *) malloc(sizeof(rtspace_t));
  /* preallocate left and right as one chunk of memory */
  rt->space->init = (Uint **) malloc((dim + 1) * sizeof(Uint*));
  initfield = (Uint *) malloc((dim + 1) * n * sizeof(Uint));
  if (rt->space->init == NULL || initfield == NULL){  
    DBG("rangetree.c: Memory allocation of temp arrays for initialisation \
failed. Exit forced.\n", NULL);
    exit(-1);
  }
  for (i = 0; i < dim + 1; i++){
    rt->space->init[i] = initfield + i * n;
  }
  for (i = 0; i < n; i++){
    rt->space->init[0][i] = pos[0][i];
  }
  /* calc height = ceil log_2 n and exp = 2^{height + 1} */
  while (level != 0){
    level >>= 1;
    height++;
    exp *= 2;
  }
  exp *= 2;
  /* preallocate for range trees */
  rt->space->allocrt = 0;
  for (i = 1; i < dim - 1; i++){
    if (i > 1){
      binom *= (height + i - 1)/(i - 1);
    }
    rt->space->allocrt += binom * (exp - 1) - (exp - 2 * n);
  }
  rt->space->rtspace = (RangeTree *) malloc(sizeof(RangeTree) * 
					    rt->space->allocrt);
  rt->space->nextrt = 0;
  /* preallocate for vebtrees */
  if (dim > 2){
    binom *= (height + dim - 2)/(dim - 2);
  }
  rt->space->allocveb = binom * (exp - 1) - (exp - 2 * n);
  #ifndef BINTREE
  rt->space->vebspace = (VebTree *) malloc(sizeof(VebTree) *
  					   rt->space->allocveb);
  #else
  rt->space->vebspace = (BinTree *) malloc(sizeof(BinTree) *
					   rt->space->allocveb);
  #endif
  rt->space->nextveb = 0;
  /* preallocate for nodes */
  rt->space->allocnode = rt->space->allocveb + rt->space->allocrt;
  rt->space->nodespace = (rtnode_t *) malloc(sizeof(rtnode_t) *
					     rt->space->allocnode);
  rt->space->nextnode = 0;
  if (rt->space->nodespace == NULL || rt->space->rtspace == NULL ||
      rt->space->vebspace == NULL){  
    DBG("rangetree.c: Memory allocation of chunks for initialisation \
failed. Exit forced.\n", NULL);
    exit(-1);
  }
  /* fprintf(stderr, ">nodes:%d veb:%d rt:%d\n",
     rt->space->allocnode, rt->space->allocveb, rt->space->allocrt);*/

  bl_rangetreeInitDim(rt->space, rt, 1, dim, 0, n, &points, sizeofelem);

  /* reallocate 
  rt->space->rtspace = (RangeTree *) realloc(rt->space->rtspace,
					 sizeof(RangeTree) * rt->space->nextrt);
  rt->space->vebspace = (VebTree *) realloc(rt->space->vebspace, sizeof(VebTree) * rt->space->nextveb);
  rt->space->nodespace = (rtnode_t *) realloc(rt->space->nodespace, sizeof(rtnode_t) * rt->space->nextnode);*/
  /* clean up */
  free(initfield);
  free(rt->space->init);
  for(i = 0; i < n; i++){
    tmppoint = (rtpoint_t*) bl_containerGet(&points, i);
    free(tmppoint->coord);
  }
  bl_containerDestruct(&points, NULL);
}

/*--------------------------- bl_rangetreeInitDim ------------------------------
 *    
 * @brief       initializes a range tree or the array at the current dimension
 * @author      Christian Otto
 *   
 */
void bl_rangetreeInitDim(rtspace_t *space, void *t, Uint currentdim, Uint dim,
			 Uint start, Uint n, Container* points,
			 size_t sizeofelem){
  if (currentdim < dim){
    RangeTree *rt = (RangeTree*) t;
    rt->dim = currentdim;
    if (space->nextnode == space->allocnode){
      DBG("rangetree.c: Too less memory allocated for nodes. Exit forced.\n",
	  NULL);
      exit(-1);
    }
    rt->root = (rtnode_t *) space->nodespace + space->nextnode;
    space->nextnode++;
    if (rt->root == NULL){
      DBG("rangetree.c: Memory allocation of root failed. Exit forced.\n",
	  NULL);
      exit(-1);
    }
    bl_rangetreeInitLev(space, rt->root, currentdim, dim, start, n, points,
			sizeofelem);
  }
  else {
    /* if currentdim == dim, initialize vebtree as lowest dimension structure */
    int obj, max, min;
    obj = space->init[currentdim - 1][start + n - 1];
    max = ((rtpoint_t*) bl_containerGet(points, obj))->coord[currentdim-1];
    obj = space->init[currentdim - 1][start];
    min = ((rtpoint_t*) bl_containerGet(points, obj))->coord[currentdim-1];
    #ifndef BINTREE
    VebTree *vt = (VebTree*) t;
    bl_vebtreeInitRange(vt, min, max - min + 1, sizeofelem);
    #else
    BinTree *vt = (BinTree*) t;
    bl_bintreeInitRange(vt, min, max - min + 1, sizeofelem);
    #endif
  }
}

/*--------------------------- bl_rangetreeInitLev ------------------------------
 *    
 * @brief   create one level of a range tree at the current dimension,
 *          tree is recursively built from root to leaves with initialization
 *          of associated substructures
 * @author  Christian Otto
 *   
 */
void bl_rangetreeInitLev(rtspace_t *space, rtnode_t *t, Uint currentdim,
			 Uint dim, Uint start, Uint n, Container* points,
			 size_t sizeofelem){
  int i, obj, lcount, rcount, lobj, robj;
  if (n < 1){
    DBG("rangetree.c: Error in size of int array: %d elements in a row. \
Exit forced.\n", bl_containerSize(points));
    exit(-1);
  }
  if (t == NULL){
    DBG("rangetree.c: Memory allocation of node was not successful. \
Exit forced.\n", NULL);
    exit(-1);
  }
  /* create leaf in rangetree at dimension currentdim */
  if (n == 1){
    obj = space->init[currentdim-1][start];
    t->value = ((rtpoint_t*) bl_containerGet(points, obj))->coord[currentdim-1];
    t->left = NULL;
    t->right = NULL;
    /* fill array pos of one dim deeper */
    space->init[currentdim][start] = space->init[currentdim-1][start];
  }
  /* recursivly split into left and right side -> merge both */
  else {
    /* length of list for left side */
    int mid = (int) (n/2);
    if (n % 2 != 0){
      mid++;
    }
    obj = space->init[currentdim-1][start + mid - 1];
    t->value = ((rtpoint_t*) bl_containerGet(points, obj))->coord[currentdim-1];
    /* recursivly split (and merging there) */
    if (space->nextnode == space->allocnode){
      DBG("rangetree.c: Too less memory allocated for nodes. Exit forced.\n",
	  NULL);
      exit(-1);
    }
    t->left = (rtnode_t *) space->nodespace + space->nextnode;
    space->nextnode++;
    bl_rangetreeInitLev(space, t->left, currentdim, dim, start, mid,
			points, sizeofelem);

    if (space->nextnode == space->allocnode){
      DBG("rangetree.c: Too less memory allocated for nodes. Exit forced.\n",
	  NULL);
      exit(-1);
    }
    t->right = (rtnode_t *) space->nodespace + space->nextnode;
    space->nextnode++; 
    bl_rangetreeInitLev(space, t->right, currentdim, dim, start + mid, n - mid,
			points, sizeofelem);

    /* merge (start, mid) with (start + mid, n - mid) at dim deeper*/
    lcount = start; rcount = start + mid;
    lobj = -1; robj = -1;
    for (i = 0; i < n; i++){
      /* left already finished */
      if (lcount >= start + mid){
	lobj = -1;
      }
      else {
	lobj = space->init[currentdim][lcount];
      }
      /* right already finished */
      if (rcount >= start + n){
	robj = -1;
      }
      else {
	robj = space->init[currentdim][rcount];
      }
      /* take left */
      if (lobj != -1 &&
	  (robj == -1 ||
	   ((rtpoint_t *) bl_containerGet(points, lobj))->coord[currentdim] <
	   ((rtpoint_t *) bl_containerGet(points, robj))->coord[currentdim])){
	space->init[dim][i] = lobj;
	lcount++;
      }
      /* take right */
      else {
	space->init[dim][i] = robj;
	rcount++;
      }
    }
    /* copy back (in situ merge not possible) */
    memcpy(space->init[currentdim] + start, space->init[dim], sizeof(Uint) * n);
  }
  /* create associated structure */
  /* next associated structure is the vebtree */
  if (currentdim + 1 == dim){
    if (space->nextveb == space->allocveb){
      DBG("rangetree.c: Too less memory allocated for vebtrees. Exit forced.\n",
	  NULL);
      exit(-1);
    }
    #ifndef BINTREE
    t->assoc = (VebTree *) space->vebspace + space->nextveb;
    #else
    t->assoc = (BinTree *) space->vebspace + space->nextveb;
    #endif
    space->nextveb++;
  }
  /* next associated structure is another rangetree */
  else {
    if (space->nextrt == space->allocrt){
      DBG("rangetree.c: Too less memory allocated for range trees. \
Exit forced.\n", NULL);
      exit(-1);
    }
    t->assoc = (RangeTree *) space->rtspace + space->nextrt;
    space->nextrt++;
  }
  if (t->assoc == NULL){
    DBG("rangetree.c: Memory allocation of assoc was not successful. \
Exit forced.\n", NULL);
    exit(-1);
  }
  bl_rangetreeInitDim(space, t->assoc, currentdim+1, dim, start, n,
		      points, sizeofelem);
}

/*------------------------- bl_rangetreeDestruct -------------------------------
 *    
 * @brief       destructs a multidimensional range tree with an array at the
 *              lowest dimension, deletes data in the array at the root
 * @author      Christian Otto
 *   
 */
void bl_rangetreeDestruct(void *t, Uint dim, void (*rmv) (void*)){
  int i;
  if (dim <= 1){
    DBG("rangetree.c: Attempt to destruct an range tree of dimension %d. \
Exit forced.\n", dim);
    exit(-1);
  }
  RangeTree *rt = (RangeTree*) t;
  rtnode_t* root = rt->root;
  for(i = 0; i < dim - 2; i++){
    root = ((RangeTree *) root->assoc)->root;
  }
  #ifndef BINTREE
  VebTree *vt;
  vt = (VebTree *) root->assoc;  
  bl_vebtreeDestruct(vt, rmv);
  #else
  BinTree *vt;
  vt = (BinTree *) root->assoc;  
  bl_bintreeDestruct(vt, rmv);  
  #endif
  bl_rangetreeDestructDim(rt, dim);
  free(rt->space->nodespace);
  free(rt->space->rtspace);
  free(rt->space->vebspace);
  free(rt->space);
}
/*------------------------ bl_rangetreeDestructDim -----------------------------
 *    
 * @brief   helper method for bl_rangetreeDestruct
 * @author  Christian Otto
 *   
 */
void bl_rangetreeDestructDim(void *t, Uint dim){
  if (dim > 1){
    RangeTree *rt = (RangeTree*) t;
    if (rt->root != NULL){
      bl_rangetreeDestructLev(rt->root, dim);
    }
  }
  else if (dim == 1){
    #ifndef BINTREE
    VebTree *vt = (VebTree*) t;
    bl_vebtreeDestruct(vt, NULL);
    #else
    BinTree *vt = (BinTree*) t;
    bl_bintreeDestruct(vt, NULL);
    #endif
  }
  else {
    DBG("rangetree.c: Attempt to destruct range tree of dimension %d. \
Exit forced.\n", dim);
    exit(-1);
  }
}

/*------------------------ bl_rangetreeDestructLev -----------------------------
 *    
 * @brief   helper method for bl_rangetreeDestruct
 * @author  Christian Otto
 *   
 */
void bl_rangetreeDestructLev(rtnode_t *t, Uint dim){
  if (t != NULL){
    t->value = 0;
    /* if internal node */
    if (t->left != NULL && t->right != NULL){
      bl_rangetreeDestructLev(t->left, dim);
      bl_rangetreeDestructLev(t->right, dim);
    }
    if (t->assoc != NULL){
      bl_rangetreeDestructDim(t->assoc, dim-1);
    }
  }
}

/*-------------------------- bl_rangetreeInsert --------------------------------
 *    
 * @brief   inserts an element at the lowest dimension (array) of the
 *          initialized multidimensional range tree,
 *          can be done in O(log^(d-1) n * log log n),
 *          currently without the use of fractional cascading
 * @author  Christian Otto
 *   
 */
void bl_rangetreeInsert(void *t, Uint dim, Uint *coord, void *data){
  if (dim < 1){
    DBG("rangetree.c: Attempt to insert in an range tree of dimension %d. \
Exit forced.\n", dim);
    exit(-1);
  }
  bl_rangetreeInsertDim(t, 1, dim, coord, data);
}

/*------------------------- rangetreeInsertDim ---------------------------------
 *    
 * @brief   helper method for bl_rangetreeInsert
 * @author  Christian Otto
 *   
 */
void bl_rangetreeInsertDim(void *t, Uint currentdim, Uint dim,
			   Uint *coord, void *data){
  if (currentdim < dim){
    RangeTree *rt = (RangeTree*) t;
    bl_rangetreeInsertLev(rt->root, currentdim, dim, coord, data);
  }
  else {
    #ifndef BINTREE
    VebTree *vt = (VebTree*) t;
    bl_vebtreeInsert(vt, coord[currentdim-1], data, NULL);
    #else
    BinTree *vt = (BinTree*) t;
    bl_bintreeInsert(vt, coord[currentdim-1], data, NULL);
    #endif
  }
}

/*------------------------- rangetreeInsertLev ---------------------------------
 *    
 * @brief   helper method for rangetreeInsert
 * @author  Christian Otto
 *   
 */
void bl_rangetreeInsertLev(rtnode_t *t, Uint currentdim, Uint dim,
			   Uint *coord, void *data){
  if (t->left != NULL && t->right != NULL && t->assoc == NULL){
    DBG("rangetree.c: Attempt to insert element in non initialized associated \
structure. Exit forced.\n", NULL);
    exit(-1);
  }
  /* if internal node */
  if (t->left != NULL && t->right != NULL){
    if (coord[currentdim-1] <= t->value){
      bl_rangetreeInsertLev(t->left, currentdim, dim, coord, data);
    }
    else {
      bl_rangetreeInsertLev(t->right, currentdim, dim, coord, data);
    }
  }
  if (t->assoc != NULL){
    bl_rangetreeInsertDim(t->assoc, currentdim+1, dim, coord, data);
  }
}

/*--------------------------- bl_rangetreeDelete -------------------------------
 *    
 * @brief   deletes an element at the lowest dimension (array) of the
 *          initialized multidimensional range tree,
 *          can be done in O(log^(d-1) n * log log n),
 *          currently without the use of fractional cascading
 * @author  Christian Otto
 *   
 */
void bl_rangetreeDelete(void *t, Uint dim, Uint *coord, void (*rmv) (void*)){
  if (dim < 1){
    DBG("rangetree.c: Attempt to delete in an range tree of dimension %d. \
Exit forced.\n", dim);
    exit(-1);
  }
  int i;
  if (dim > 1){
    RangeTree *rt = (RangeTree*) t;
    rtnode_t* root = rt->root;
    for(i = 0; i < dim - 2; i++){
      rt = (RangeTree *) root->assoc;
      root = rt->root;
    }
    #ifndef BINTREE
    VebTree *vt;
    vt = (VebTree *) root->assoc;  
    free(bl_vebtreeDelete(vt, coord[dim-1], rmv));
    #else
    BinTree *vt;
    vt = (BinTree *) root->assoc;  
    free(bl_bintreeDelete(vt, coord[dim-1], rmv));    
    #endif
    bl_rangetreeDeleteDim(t, 1, dim, coord);
  }
  else {
    #ifndef BINTREE
    VebTree *vt;
    vt = (VebTree *) t;
    free(bl_vebtreeDelete(vt, coord[dim-1], rmv));
    #else
    BinTree *vt;
    vt = (BinTree *) t;
    free(bl_bintreeDelete(vt, coord[dim-1], rmv));
    #endif
  }
}

/*-------------------------- bl_rangetreeDeleteDim -----------------------------
 *    
 * @brief   helper method for bl_rangetreeDelete
 * @author  Christian Otto
 *   
 */
void bl_rangetreeDeleteDim(void *t, Uint currentdim, Uint dim, Uint *coord){
  if (currentdim < dim){
    RangeTree *rt = (RangeTree*) t;
    bl_rangetreeDeleteLev(rt->root, currentdim, dim, coord);
  }
  else {
    #ifndef BINTREE
    VebTree *vt = (VebTree*) t;
    free(bl_vebtreeDelete(vt, coord[currentdim-1], NULL));
    #else
    BinTree *vt = (BinTree*) t;
    free(bl_bintreeDelete(vt, coord[currentdim-1], NULL));
    #endif
  }
}

/*-------------------------- bl_rangetreeDeleteLev -----------------------------
 *    
 * @brief   helper method for bl_rangetreeDelete
 * @author  Christian Otto
 *   
 */
void bl_rangetreeDeleteLev(rtnode_t *t, Uint currentdim, Uint dim, Uint *coord){
  if (t->left != NULL && t->right != NULL && t->assoc == NULL){
    DBG("rangetree.c: Attempt to delete element from non initialized \
associated structure. Exit forced.\n", NULL);
    exit(-1);
  }
  /* if internal node */
  if (t->left != NULL && t->right != NULL){
    if (coord[currentdim-1] <= t->value){
      bl_rangetreeDeleteLev(t->left, currentdim, dim, coord);
    }
    else {
      bl_rangetreeDeleteLev(t->right, currentdim, dim, coord);
    }
  }
  if (t->assoc != NULL){
    bl_rangetreeDeleteDim(t->assoc, currentdim+1, dim, coord);
  }
}

/*--------------------------- bl_rangetreeGetData ------------------------------
 *    
 * @brief   returns the data to a given coord (index) at a given dimension
 *          in O(log n) in the d-dimensional range tree,
 *          based on the fact that at each dimension the coord (index) is
 *          unique for all points stored in the range tree
 * @author  Christian Otto
 *   
 */
void* bl_rangetreeGetData(void *t, Uint dim, Uint coord, Uint cdim){
  if (dim < 1 || cdim < 1 || cdim > dim){
    DBG("rangetree.c: Attempt to retrieve data with coordinate '%d' at \
dimension %d from an range tree of dimension %d. Exit forced.\n",
	coord, cdim, dim);
    exit(-1);
  }
  return bl_rangetreeGetDataDim(t, 1, dim, coord, cdim);
}

/*-------------------------- bl_rangetreeGetDataDim ----------------------------
 *    
 * @brief   helper method for bl_rangetreeGetData
 * @author  Christian Otto
 *   
 */
void* bl_rangetreeGetDataDim(void *t, Uint currentdim, Uint dim,
			     Uint coord, Uint cdim){
  /* if not the dimension to search in */
  if (currentdim != cdim){
    /* go to next dimension */
    if (currentdim < dim){
      RangeTree *rt = (RangeTree*) t;
      return bl_rangetreeGetDataDim(rt->root->assoc, currentdim+1,
				    dim, coord, cdim);
    }
    /* last dimension reached */
    else {
      #ifndef BINTREE
      VebTree *vt = (VebTree*) t;
      /* only one data object should remain */
      if (bl_vebtreeMin(vt) != bl_vebtreeMax(vt)){
	DBG("rangetree.c: vebtree should contain only one element.\
Exit forced.\n", NULL);
	exit(-1);
      }
      return bl_vebtreeGetData(vt, bl_vebtreeMin(vt));
      #else
      BinTree *vt = (BinTree*) t;
      /* only one data object should remain */
      if (bl_bintreeMin(vt) != bl_bintreeMax(vt)){
	DBG("rangetree.c: bintree should contain only one element.\
Exit forced.\n", NULL);
	exit(-1);
      }
      return bl_bintreeGet(vt, bl_bintreeMin(vt));
      #endif
    }
  }
  /* else -> search dimension reached */
  else {
    /* if not lowest dimension -> search in leaves */
    if (currentdim < dim){
      RangeTree *rt = (RangeTree*) t;
      return bl_rangetreeGetDataLev(rt->root, currentdim, dim, coord, cdim);
    }
    /* else -> search in vebtree */
    else {
      #ifndef BINTREE
      VebTree *vt = (VebTree*) t;
      return bl_vebtreeGetData(vt, coord);
      #else
      BinTree *vt = (BinTree*) t;
      return bl_bintreeGet(vt, coord);
      #endif
    }
  }
}

/*-------------------------- bl_rangetreeGetDataLev ----------------------------
 *    
 * @brief   helper method for bl_rangetreeGetData
 * @author  Christian Otto
 *   
 */
void* bl_rangetreeGetDataLev(rtnode_t *t, Uint currentdim, Uint dim,
			     Uint coord, Uint cdim){
  /* if internal node */
  if (t->left != NULL && t->right != NULL){
    if (coord <= t->value){
      return bl_rangetreeGetDataLev(t->left, currentdim, dim, coord, cdim);
    }
    else {
      return bl_rangetreeGetDataLev(t->right, currentdim, dim, coord, cdim);
    }
  }
  else {
    if (coord == t->value){
      return bl_rangetreeGetDataDim(t->assoc, currentdim+1, dim, coord, cdim);
    }
    else {
      return NULL;
    }
  }
}

/*---------------------------- bl_rangetreeOutput ------------------------------
 *    
 * @brief   outputs the rangetree
 * @author  Christian Otto
 *   
 */
void bl_rangetreeOutput(void *t, Uint dim){
  if (dim < 1){
    DBG("Attempt to output an range tree of dimension %d. Exit forced.", dim);
    exit(-1);
  }
  printf("\nrangetree of dim %d:\n", dim);
  Uint *tmpcoord = (Uint*) calloc(dim, sizeof(Uint));
  bl_rangetreeOutputDim(t, 1, dim, tmpcoord);
  free(tmpcoord);
  printf("\n");
}

/*--------------------------- bl_rangetreeOutputDim ----------------------------
 *    
 * @brief   helper method for bl_rangetreeOutput
 * @author  Christian Otto
 *   
 */
void bl_rangetreeOutputDim(void *t, Uint currentdim, Uint dim, Uint *tmpcoord){
  if (currentdim < dim){
    RangeTree *rt = (RangeTree*) t;
    bl_rangetreeOutputLev(rt->root, currentdim, dim, tmpcoord);
  }
  else {
    int i;
    #ifndef BINTREE
    VebTree *vt = (VebTree*) t;
    for(i = 0; i < currentdim - 1; i++){
      if (i != 0){
	printf("->");
      }
      printf("%d", tmpcoord[i]);
    }
    printf(": assoc: ");
    if (bl_vebtreeIsEmpty(vt)){
      printf("empty [%d..%d]\n", vt->min, vt->allocelem + vt->min - 1);
      return;
    }
    int value = bl_vebtreeMin(vt);
    while(value != -1){
      printf("%d ", value);
      value = bl_vebtreeSucc(vt, value);
    }
    printf(" [%d..%d]",vt->min, vt->allocelem + vt->min - 1);
    printf("\n");
    #else
    BinTree *vt = (BinTree*) t;
    for(i = 0; i < currentdim - 1; i++){
      if (i != 0){
	printf("->");
      }
      printf("%d", tmpcoord[i]);
    }
    printf(": assoc: ");
    if (bl_bintreeIsEmpty(vt)){
      printf("empty [%d..%d]\n", vt->min, vt->allocelem + vt->min - 1);
      return;
    }
    int value = bl_bintreeMin(vt);
    while(value != -1){
      printf("%d ", value);
      value = bl_bintreeSucc(vt, value);
    }
    printf(" [%d..%d]",vt->min, vt->allocelem + vt->min - 1);
    printf("\n");
    #endif
  }
}

/*--------------------------- bl_rangetreeOutputLev ----------------------------
 *    
 * @brief   helper method for bl_rangetreeOutput
 * @author  Christian Otto
 *   
 */
void bl_rangetreeOutputLev(rtnode_t *t, Uint currentdim, Uint dim,
			   Uint *tmpcoord){
  int i;
  tmpcoord[currentdim-1] = t->value;
  if (t->left != NULL && t->right != NULL){
    for(i = 0; i < currentdim; i++){
      if (i != 0)
	printf("->");
      printf("%d",tmpcoord[i]); 
    }
    printf(": left:%d, right:%d\n",((rtnode_t*)t->left)->value,
	                           ((rtnode_t*)t->right)->value);
    bl_rangetreeOutputLev(t->left, currentdim, dim, tmpcoord);
    bl_rangetreeOutputLev(t->right, currentdim, dim, tmpcoord);
  }
  if (t->assoc != NULL){
    if (t->left == NULL && t->right == NULL){
      printf("leaf:");
    }
    tmpcoord[currentdim-1] = t->value;
    bl_rangetreeOutputDim(t->assoc, currentdim+1, dim, tmpcoord);
  }
}
