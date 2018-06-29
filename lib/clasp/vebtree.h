#ifndef VEBTREE_H
#define VEBTREE_H

/**
 * vebtree.h
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
 * Id: $Id: vebtree.h 114 2010-04-19 13:59:41Z steve $
 * Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/vebtree.h $
 */
 
#include <stdlib.h>
#include "basic-types.h"
#include "container.h"

/*
 * Macros
 */
#define RTOP(f)      (((t->treespace + f - 1)->right < f)?((t->treespace + f - 1)->right):(f))
#define RLEAF(f)     (((t->treespace + f - 1)->right > f)?((t->treespace + f - 1)->right):(f))
#define LTOP(f)	     (((t->treespace + f - 1)->left < f)?((t->treespace + f - 1)->left):(f))
#define LLEAF(f)     (((t->treespace + f - 1)->left > f)?((t->treespace + f - 1)->left):(f))

//#define VEB_DEFERRED

/*
 * Typedef
 */
typedef struct {
  /* Caution: all data refers to nodes {1, ..., N} */

  /* 
   * -1              if node is no left chain
   * lefttop(node)   if q is in left chain from leaf
   * leftleaf(node)  if q is the top of the left chain
   *  0              otherwise
   */
  int left;
  int right;
  /* number of distinct left/right chains containing this node */
  int lref;
  int rref;
} VebNode;

typedef struct {
  VebNode *treespace;
  void *dataspace;
  size_t sizeofelem;
  int base;			/* first leaf in tree {1,...,N} */
  int allocelem;		/* maximal number of data leaves */
  int height;			/* length of path from leaf to root */
  int min;
  int next;                     /* next element after dummy element */
} VebTree;

/*
 * Prototypes
 */
void bl_vebtreeInit(VebTree *t, int allocelem, size_t size);
void bl_vebtreeInitRange(VebTree *t, int min, int allocelem, size_t size);
void bl_vebtreeAlloc(VebTree *t);
void bl_vebtreeDestruct(VebTree *t, void (*rmv) (void *));
BOOL bl_vebtreeIsEmpty (VebTree *t);
BOOL bl_vebtreeIsActive(VebTree *t, int value);
VebNode* bl_vebtreeGetNode(VebTree *t, int value);
void* bl_vebtreeGetData(VebTree *t, int value);
int bl_vebtreeMin(VebTree *t);
int bl_vebtreeMax(VebTree *t);
int bl_vebtreePred(VebTree *t, int value);
int bl_vebtreeSucc(VebTree *t, int value);
void bl_vebtreeInsert(VebTree *t, int value, void *, void (*rmv)(void *));
void* bl_vebtreeDelete(VebTree *t, int value, void (*rmv) (void *));

/*
 * helper methods
 */
Container* bl_vebtreeBinSearch(VebTree *t, int f, int g);
Container* bl_vebtreeBinSearchRed(VebTree *t, int f, int g, int w);
int bl_vebtreeLabelRuleR(VebTree *t, int f);
int bl_vebtreeLabelRuleL(VebTree *t, int f);
int bl_vebtreeEqualRule(VebTree *t, int f, int g);
int bl_vebtreeEqualRuleL(VebTree *t, int f, int g);
int bl_vebtreeEqualRuleR(VebTree *t, int f, int g);
void bl_vebtreeIRestoreL(VebTree *t, int f, int g, int b);
void bl_vebtreeIRestoreR(VebTree *t, int f, int h, int c);
void bl_vebtreeDRestoreL(VebTree *t, int f);
void bl_vebtreeDRestoreR(VebTree *t, int f);

#endif /* VEBTREE_H */ 
