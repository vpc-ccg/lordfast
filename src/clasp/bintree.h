#ifndef BINTREE_H
#define BINTREE_H

/**
 * bintree.h
 * a binary search tree with size $n=2^k$
 * over an search space 1...n (within the
 * bintree the priorities are from 0...n-1
 * due to using arrays starting from 0)
 * implemented as a 'static' field
 *
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Wed Jun  4 13:59:01 CEST 2008
 */
  
/*
 * SVN
 * Revision of last commit: $Rev: 114 $
 * Author: $Author: steve $
 * Date: $Date: 2010-04-19 15:59:41 +0200 (Mon, 19 Apr 2010) $
 * Id: $Id: bintree.h 114 2010-04-19 13:59:41Z steve $
 * Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/bintree.h $
 */
 
#include <stdlib.h>
#include <stdio.h>
#include "basic-types.h"

//#define BIN_DEFERRED

/*
 * Typedef
 */
typedef struct {
  void *treespace;
  int *father,
    *left,
    *right;
  BOOL *used;
  int allocelem;
  int min;
  int root;
  int numofelems;
  size_t sizeofelem;
} BinTree;

/*
 * Prototypes
 */
void bl_bintreeInit(BinTree *t, int allocelem, size_t size);
void bl_bintreeInitRange(BinTree *t, int min, int allocelem, size_t size);
void bl_bintreeAlloc(BinTree *t);
void bl_bintreeSetFather(BinTree *t, int father, int level);
void bl_bintreeDestruct(BinTree *t, void (*rmv)(void *));
BOOL bl_bintreeIsEmpty(BinTree *t);
void bl_bintreeResize(BinTree *t);
void bl_bintreeInsert(BinTree *t, int value, void *, void (*rmv)(void *));
void* bl_bintreeGet(BinTree *t, int value);
void* bl_bintreeDelete(BinTree *t, int value, void (*rmv)(void *));
int bl_bintreeMin(BinTree *t);
int bl_bintreeMax(BinTree *t);
int bl_bintreePred(BinTree *t, int value);
int bl_bintreeSucc(BinTree *t, int value);
Uint bl_bintreeSize(BinTree *t);

#endif /* BINTREE_H */
