#ifndef RANGETREE_H
#define RANGETREE_H

/**
 * rangetree.h
 * an implementation of an d-dimensional semi-dynamic range tree
 * with O(n log^(d-1) n) in space, construction in O(n log^(d-1) n) in time
 * query time O(log^(d-1) n) with fractional cascading,
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
 * Id: $Id: rangetree.h 114 2010-04-19 13:59:41Z steve $
 * Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/rangetree.h $
 */

#include <stdlib.h>
#include <stdio.h>
#include "basic-types.h"
#include "container.h"
#include "vebtree.h"
#include "bintree.h"

//#define BINTREE

/*
 * Typedef
 */
typedef struct {
  Uint value;
  void *left;  /* rtnode_t* */
  void *right; /* same as above */
  void *assoc; /* VebTree* / BinTree* or RangeTree* or data (at leaves) */
} rtnode_t;

typedef struct {
  Uint **init;
  void *nodespace;
  int nextnode;
  int allocnode;
  void *rtspace;
  int nextrt;
  int allocrt;
  void *vebspace;
  int nextveb;
  int allocveb;
} rtspace_t;

typedef struct {
  rtspace_t *space;
  rtnode_t *root;
  Uint dim;
} RangeTree;

/*
 * helper type for init
 */
typedef struct {
  Uint *coord;
  Uint dim;
} rtpoint_t;


/*
 * Prototypes
 */
void bl_rangetreeInit(void *t, Uint dim, Uint n, Uint** pos,
		      size_t sizeofelem);
void bl_rangetreeDestruct(void *t, Uint dim, void (*rmv)(void*));
void bl_rangetreeInsert(void *t, Uint dim, Uint *coord, void *data);
void bl_rangetreeDelete(void *t, Uint dim, Uint *coord, void (*rmv) (void*));
void *bl_rangetreeGetData(void *t, Uint tdim, Uint coord, Uint cdim);
void bl_rangetreeOutput(void *t, Uint dim);

/* helper methods */
void bl_rangetreeInitDim(rtspace_t *space, void *t, Uint currentdim, Uint dim,
			 Uint start, Uint n, Container* points,
			 size_t sizeofelem);
void bl_rangetreeInitLev(rtspace_t *space, rtnode_t *t, Uint currentdim, 
			 Uint dim, Uint start, Uint n, Container* points,
			 size_t sizeofelem);
void bl_rangetreeDestructDim(void *t, Uint dim);
void bl_rangetreeDestructLev(rtnode_t *t, Uint dim);
void bl_rangetreeInsertDim(void *t, Uint currentdim, Uint dim,
			   Uint *coord, void *data);
void bl_rangetreeInsertLev(rtnode_t *t, Uint currentdim, Uint dim,
			   Uint *coord, void *data);
void bl_rangetreeDeleteDim(void *t, Uint currentdim, Uint dim, Uint *coord);
void bl_rangetreeDeleteLev(rtnode_t *t, Uint currentdim, Uint dim, Uint *coord);
void *bl_rangetreeGetDataDim(void *t, Uint currentdim, Uint tdim,
	 		     Uint coord, Uint cdim);
void *bl_rangetreeGetDataLev(rtnode_t *t, Uint currentdim, Uint tdim,
			     Uint coord, Uint cdim);
void bl_rangetreeOutputDim(void *t, Uint currentdim, Uint dim, Uint *tmpcoord);
void bl_rangetreeOutputLev(rtnode_t *t, Uint currentdim, Uint dim,
			   Uint *tmpcoord);

#endif /* RANGE_TREE_H */
