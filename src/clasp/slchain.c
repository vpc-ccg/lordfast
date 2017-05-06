
/**
 * slchain.c
 * chaining of non overlapping fragments
 * with different possible gap costs
 * Note that all given fragments must be from the same reference sequence
 *
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Wed Jun  4 16:13:58 CEST 2008
 */

/* 
 * SVN
 * Revision of last commit: $Rev: 116 $
 * Author: $Author: steve $
 * Date: $Date: 2010-06-30 13:51:27 +0200 (Wed, 30 Jun 2010) $
 * Id: $Id: slchain.c 116 2010-06-30 11:51:27Z steve $
 * Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/slchain.c $ 
 */

#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include "sltypes.h"
#include "slchain.h"
#include "container.h"
#include "debug.h"
#include "mathematics.h"
#include "sort.h"
#include "vqueue.h"
#include "vebtree.h"
#include "bintree.h"
#include "rangetree.h"



/*----------------------------- bl_slExtractPoints -----------------------------
 *    
 * @brief       extracts all start and end points from the matches,
 *	        returns a sorted list according to the position in the sequence,
 *              runs in O(n log n) (consists of one sort and two linear scans),
 * @parameter   array of slmatch_t, assumed to be presorted by 
 *              start position on sequence
 * @author      Christian Otto
 *   
 */
Container* bl_slExtractPoints(slmatch_t *fragments, Uint size){
  int i, j, *space = NULL;
  Uint *sorted;
  Container *points;
  point_t point;
  slmatch_t *amatch, *bmatch;

  /* sort the indexes of the container src by sequence position of end points */
  sorted = quickSort(space, fragments, size, cmp_slmatch_end_quick, NULL);
  /* initialize data structures */
  points = (Container *) malloc(sizeof(Container));
  bl_containerInit(points, 1000, sizeof(point_t));

  /* traverse all start points according to the order */
  for (i = 0, j = 0; i < size; i++){
    amatch = fragments + i;
    bmatch = fragments + sorted[j];
    /* if end point is before sequence pos of next start point */
    while (FSTART_S(amatch) > FEND_S(bmatch)){
      point.x = FEND_S(bmatch);
      point.y = FEND_Q(bmatch);
      point.index = sorted[j];
      point.start = 0;
      bl_containerAdd(points, &point);
      bmatch = fragments + sorted[++j];
    }
    point.x = FSTART_S(amatch);
    point.y = FSTART_Q(amatch);
    point.index = i;
    point.start = 1;
    bl_containerAdd(points, &point);
  }
  /* assumes that end point is always later than start point */
  while (j < size){
    bmatch = fragments + sorted[j];
    point.x = FEND_S(bmatch);
    point.y = FEND_Q(bmatch);
    point.index = sorted[j];
    point.start = 0;
    bl_containerAdd(points, &point);
    j++;
  }
  free(sorted);
  return points;
}

/*----------------------------- bl_slClusterSop --------------------------------
 *
 * @brief       clustering of fragments for possible chains
 * @author      Christian Otto
 *
 */
void bl_slClusterLin(slmatch_t *fragments, Uint size, double eps, double lambda, int maxgap){
  Uint i, begin, length;
  int min_yst, max_yend;
  double max_score_per_pos;
  slmatch_t *current, *a, *b;
  /* nothing to do with size 0 */
  if (size == 0){
    return;
  }
  /* sorting (if already sorted -> only linear scan is required) */
  qsort(fragments, size, sizeof(slmatch_t), cmp_slmatch_qsort);

  /* initializations */
  begin = 0;
  a = fragments;
  max_score_per_pos = a->scr / (double) FLEN_Q(a);
  min_yst = FSTART_Q(a);
  max_yend = FEND_Q(a);
  for (i = 0; i < size - 1; i++){
    b = fragments + (i+1);
    /* update start and end bounds */
    if (FSTART_Q(b) < min_yst){
      min_yst = FSTART_Q(b);
    }
    if (FEND_Q(b) > max_yend){
      max_yend = FEND_Q(b);
    }
    /* update max_score_per_pos */
    if (b->scr / (double) FLEN_Q(b) > max_score_per_pos){
      max_score_per_pos = b->scr/(double) FLEN_Q(b);
    }
    /*
     * fragments are on difference reference sequences or 
     * if gap between adjacent fragments on sequence is too large or
     * if start >= end and min(G_LIN) > max_score -> new cluster
     */
    if (a->subject != b->subject || MAXGAP2(b, a) ||
	((int) FSTART_S(b) > (int) FEND_S(a) &&
	 (double) lambda * (double) DX(b, a) > (double) (max_yend - min_yst + 1) * max_score_per_pos)){
      length = i - begin + 1;

      /* no chaining necessary */
      if (length == 1){
      	current = fragments + i;
      	current->chain = (slchain_t *) malloc(sizeof(slchain_t));
      	bl_slchainInitFromMatch(current->chain, current);
      }
      /* chain cluster */
      else {
	bl_slChainLin(fragments + begin, length, eps, lambda, maxgap);
      }
      begin = i + 1;
    }
    a = b;
  }
  length = i - begin + 1;
  /* no chaining necessary */
  if (length == 1){
    current = fragments + begin;
    current->chain = (slchain_t *) malloc(sizeof(slchain_t));
    bl_slchainInitFromMatch(current->chain, current);
  }
  /* chain cluster */
  else {
    bl_slChainLin(fragments + begin, length, eps, lambda, maxgap);
  }
}

/*------------------------------ bl_slChainLin ---------------------------------
 *    
 * @brief      creates chains of non overlapping fragments from a list of
 *	       fragments with the highest score using linear gap costs (L1)
 * @parameter  fragments - the list of fragments to be chained,
 *	       lambda - cost of aligning gap with character on sequence,
 *	       eps - cost of aligning gap with character on query 
 * @author     Christian Otto
 *   
 */
void bl_slChainLin(slmatch_t *fragments, Uint size, double eps, double lambda, int maxgap){
  Uint i, *trans, *sorted, *space = NULL;
  int pred, succ;
  PairSint t;
  point_t *point;
  slmatch_t *current, *first;
  slchain_t *chain, *cand, **prev, **predmatch, **succmatch;
  #ifndef BINTREE
  VebTree veb;
  #else
  BinTree bin;
  #endif
  Container *points;
  VQueue refs;
  /* no chaining if size is zero */
  if (size == 0){
    return;
  }
  /* sorting (if already sorted -> only one scan required) */
  qsort(fragments, size, sizeof(slmatch_t), cmp_slmatch_qsort);
  /* preinitialize */
  points = bl_slExtractPoints(fragments, size);
  /* sorting by query position (required for vebtree/bintree) */
  trans = (Uint *) malloc(sizeof(Uint) * bl_containerSize(points));
  sorted = quickSort(space, points->contspace, bl_containerSize(points),
		     cmp_slmatch_trans_first_y, NULL); 
  for (i = 0; i < bl_containerSize(points); i++){
    trans[sorted[i]] = i;
  }
  /* get t.a and t.b */  
  t.a = ((point_t*) bl_containerGet(points, bl_containerSize(points) - 1))->x;
  t.b = ((point_t*) bl_containerGet(points,
				    sorted[bl_containerSize(points) - 1]))->y;
  free(sorted);
  /* initializations for chaining itself */
  prev = (slchain_t **) malloc(sizeof(slchain_t *) * size);
  memset(prev, 0, sizeof(slchain_t*) * size);
  #ifndef BINTREE
  bl_vebtreeInit(&veb, bl_containerSize(points), sizeof(slchain_t *));
  #else
  bl_bintreeInit(&bin, bl_containerSize(points), sizeof(slchain_t *));  
  #endif
  bl_vqueueInit(&refs, 2 * size, sizeof(slchain_t *));

  /* traverse all points (sorted by sequence start position) */
  for (i = 0; i < bl_containerSize(points); i++){
    point = (point_t *) bl_containerGet(points, i);

    /* point is the start point of fragment point->index (in fragments) */
    if (point->start == 1){
      /* get current match from fragments */
      current = fragments + point->index;

      /* RMQ in vebtree/bintree by y-coordinate */
      #ifndef BINTREE
      pred = bl_vebtreePred(&veb, trans[i]);
      predmatch = (slchain_t **) bl_vebtreeGetData(&veb, pred);
      #else
      pred = bl_bintreePred(&bin, trans[i]);
      predmatch = (slchain_t **) bl_bintreeGet(&bin, pred);
      #endif

      if (predmatch != NULL && !MAXGAP(current, *predmatch)){
	/* only take prev if it increases chain score */	
	if ((*predmatch)->scr >= (double) GLIN(current, *predmatch)){
	  prev[point->index] = *predmatch;
	}
	/* 
	 * addition due to clusters of local chains
	 * (differs from Abouelhoda2004)
	 */
	else if (current->scr >= (double) GLIN(current, *predmatch)){
	  first = *(slmatch_t **) bl_containerGet((*predmatch)->matches, 0);
	  if (first->chain != NULL &&
	      ((slchain_t *)first->chain)->scr < (*predmatch)->scr +
	      current->scr - GLIN(current, *predmatch)){
	    chain = (slchain_t *) malloc(sizeof(slchain_t));
	    bl_slchainInit(chain);
	    bl_slchainCopy(chain, *predmatch);
	    chain->j = FEND_Q(current) - FSTART_Q(*predmatch) + 1;
	    chain->q = FEND_S(current) - FSTART_S(*predmatch) + 1;
	    chain->scr += current->scr - (double) GLIN(current, *predmatch);
	    bl_containerAdd(chain->matches, &current);
	    bl_slchainpDestruct(&first->chain);
	    first->chain = chain;
	  }
	}
      }    
    }
    /* point is the end point of fragment point->index */
    else {
      current = fragments + point->index;
      chain = (slchain_t *) malloc(sizeof(slchain_t));
      /* 
       * chaining candidate found
       * (can be only better chain for prev or best chain for current as well)
       */
      if (prev[point->index] != NULL){
	cand = prev[point->index];
	bl_slchainInit(chain);
	chain->i = FSTART_Q(cand);
	chain->j = FEND_Q(current) - FSTART_Q(cand) + 1;
	chain->p = FSTART_S(cand);
	chain->q = FEND_S(current) - FSTART_S(cand) + 1;
	chain->scr = current->scr + cand->scr - (double) GLIN(current, cand);
	bl_containerMerge(chain->matches, cand->matches);
	bl_containerAdd(chain->matches, &current);
	if (current->subject != cand->subject){
	  DBG("slchain.c: Attempt to chain fragments of different \
reference sequences.\nExit forced.\n\n", NULL);
	  exit(-1);
	}
	chain->subject = current->subject;

	/* update best chain of cluster */
	first = *(slmatch_t **) bl_containerGet(chain->matches, 0);
	if (first->chain != NULL){
	  if (((slchain_t *) first->chain)->scr <= chain->scr){
	    slchain_t *tmp = (slchain_t *) malloc(sizeof(slchain_t));
	    bl_slchainInit(tmp);
	    bl_slchainCopy(tmp, chain);
	    bl_slchainpDestruct(&first->chain);
	    first->chain = (void *) tmp;
	  }
	}
      }
      /* no chaining candidate found */
      else {
	bl_slchainInitFromMatch(chain, current);	
	if (current->chain != NULL){
	  DBG("slchain.c: Error in updating of best chain. \
Exit forced.\n", NULL);
	  exit(-1);
	}
	else {
	  /* update best chain of current */
	  slchain_t *tmp = (slchain_t *) malloc(sizeof(slchain_t));
	  bl_slchainInit(tmp);
	  bl_slchainCopy(tmp, chain);
	  current->chain = tmp;
	}
      }
      chain->prio = chain->scr - GCLIN(chain);
      #ifndef BINTREE
      pred = bl_vebtreePred(&veb, trans[i]);
      predmatch = (slchain_t **) bl_vebtreeGetData(&veb, pred);
      if (predmatch == NULL || chain->prio >= (*predmatch)->prio){
	bl_vebtreeInsert(&veb, trans[i], &chain, NULL);
	succ = bl_vebtreeSucc(&veb, trans[i]);
	succmatch = (slchain_t **) bl_vebtreeGetData(&veb, succ);
	while(succmatch != NULL && chain->prio > (*succmatch)->prio){
	  free(bl_vebtreeDelete(&veb, succ, NULL));
	  succ = bl_vebtreeSucc(&veb, succ);
	  succmatch = (slchain_t **) bl_vebtreeGetData(&veb, succ);
	}
      }
      #else
      pred = bl_bintreePred(&bin, trans[i]);
      predmatch = (slchain_t **) bl_bintreeGet(&bin, pred);
      if (predmatch == NULL || chain->prio >= (*predmatch)->prio){
	bl_bintreeInsert(&bin, trans[i], &chain, NULL);
	succ = bl_bintreeSucc(&bin, trans[i]);
	succmatch = (slchain_t **) bl_bintreeGet(&bin, succ);
	while(succmatch != NULL && chain->prio > (*succmatch)->prio){
	  free(bl_bintreeDelete(&bin, succ, NULL));
	  succ = bl_bintreeSucc(&bin, succ);
	  succmatch = (slchain_t **) bl_bintreeGet(&bin, succ);
	}
      }
      #endif
      bl_vqueueEnqueue(&refs, &chain);
    } /* END OF else */
  } /* END OF for (i = 0; i < bl_containerSize(points); i++) */
  bl_memusage();
  free(prev);
  bl_vqueueDestruct(&refs, bl_slchainpDestruct);
  /* destruct vebtree/bintree */
  #ifndef BINTREE
  bl_vebtreeDestruct(&veb, NULL);
  #else  
  bl_bintreeDestruct(&bin, NULL);
  #endif
  /* destruct points */
  bl_containerDestruct(points, NULL);
  free(points);
  /* destruct trans */
  free(trans);

  return;
}

/*------------------------------- bl_slGetTrans --------------------------------
 *    
 * @brief       get transformation T_1 and T_2 of each point and sort by x and y
 * @return      array of four sorted arrays with index positions of points
 * @author      Christian Otto
 *   
 */
Uint** bl_slGetTrans(Container* points){
  int i, j, size = 4;
  int* space = NULL;
  point_t *point = NULL;
  /* size = dim! transformations * dim dimensions */
  Uint** trans = (Uint**) malloc(size * sizeof(Uint*));
  Uint (*cmp_trans[4])(Uint,Uint,void*,void*);
  cmp_trans[0] = &cmp_slmatch_trans_first_x;
  cmp_trans[1] = &cmp_slmatch_trans_first_y;
  cmp_trans[2] = &cmp_slmatch_trans_second_x;
  cmp_trans[3] = &cmp_slmatch_trans_second_y;
  for(i = 0; i < size; i++){
    trans[i] = quickSort(space, points->contspace, bl_containerSize(points),
			 cmp_trans[i], NULL); 
  }
  for(i = 0; i < bl_containerSize(points); i++){
    point = (point_t *) bl_containerGet(points, i);
    point->size = size;
    point->trans = (Uint*) malloc(size * sizeof(Uint));
  }
  for(i = 0; i < size; i++){
    for(j = 0; j < bl_containerSize(points); j++){
      point = (point_t *) bl_containerGet(points, trans[i][j]);
      point->trans[i] = j;
    }
  }
  return trans;
}

/*-------------------------- cmp_slmatch_qsort ---------------------------------
 *    
 * @brief      compare of two slmatch_t by seq_end
 * @author     Christian Otto
 *   
 */
int cmp_slmatch_qsort (const void * a, const void * b){
  if (((slmatch_t *) a)->subject != ((slmatch_t *) b)->subject){
    return ((slmatch_t *) a)->subject - ((slmatch_t *) b)->subject;
  }
  else {
    return (FSTART_S((slmatch_t *) a) - FSTART_S((slmatch_t *) b));
  }
}

/*------------------------ cmp_slmatch_end_quick -------------------------------
 *    
 * @brief      compare of two slmatch_t by seq_end
 * @author     Christian Otto
 *   
 */
Uint cmp_slmatch_end_quick (Uint a, Uint b, void *data, void *info){
  slmatch_t *matches = (slmatch_t *) data;
  /* sort seq end */
  if (FEND_S(matches + a) > FEND_S(matches + b)){
    return 1;
  }
  else if (FEND_S(matches + a) < FEND_S(matches + b)){
    return 2;
  }
  else{
    return 0;
  }
}

/*---------------------- cmp_slmatch_trans_first_x -----------------------------
 *    
 * @brief      compare of two point_t by T_1.x
 *             and sort end points before start points,
 *             T_1:(x, y) -> (x - y, y) (x on seq, y on qry)
 * @author     Christian Otto
 *   
 */
Uint cmp_slmatch_trans_first_x (Uint a, Uint b, void *data, void *info){
  point_t *matches = (point_t *) data;
  int ax = 0, bx = 0; 
  ax = matches[a].x - matches[a].y;
  bx = matches[b].x - matches[b].y;
  if (ax > bx){
    return 1;
  }
  else if (ax < bx){
    return 2;
  }  
  else if (matches[a].start > matches[b].start){
   return 1;
  }
  else if (matches[a].start < matches[b].start){
    return 2;
  }
  else {
    return 0;
  }
}

/*---------------------- cmp_slmatch_trans_first_y -----------------------------
 *    
 * @brief      compare of two point_t by T_1.y
 *             and sort start points before end points,
 *             T_1:(x, y) -> (x - y, y) (x on seq, y on qry)
 * @author     Christian Otto
 *   
 */
Uint cmp_slmatch_trans_first_y (Uint a, Uint b, void *data, void *info){
  point_t *matches = (point_t *) data;
  int ay, by;
  ay = matches[a].y;
  by = matches[b].y;
  if (ay > by){
    return 1;
  }
  else if (ay < by){
    return 2;
  }  
  else if (matches[a].start < matches[b].start){
    return 1;
  }
  else if (matches[a].start > matches[b].start){
    return 2;
  }
  else {
    return 0;
  }
}

/*--------------------- cmp_slmatch_trans_second_x -----------------------------
 *    
 * @brief      compare of two point_t by T_2.x
 *             and sort start points before end points,
 *             T_2:(x, y) -> (x, y - x) (x on seq, y on qry)
 * @author     Christian Otto
 *   
 */
Uint cmp_slmatch_trans_second_x (Uint a, Uint b, void *data, void *info){
  point_t *matches = (point_t *) data;
  int ax, bx;
  ax = matches[a].x;
  bx = matches[b].x;
  if (ax > bx){
    return 1;
  }
  else if (ax < bx){
    return 2;
  }  
  else if (matches[a].start < matches[b].start){
    return 1;
  }
  else if (matches[a].start > matches[b].start){
    return 2;
  }
  else {
    return 0;
  }
}

/*--------------------- cmp_slmatch_trans_second_y -----------------------------
 *    
 * @brief      compare of two point_t by T_2.y
 *             and sort end points before start points,
 *             T_2:(x, y) -> (x, y - x) (x on seq, y on qry)
 * @author     Christian Otto
 *   
 */
Uint cmp_slmatch_trans_second_y (Uint a, Uint b, void *data, void *info){
  point_t *matches = (point_t *) data;
  int ay, by;
  ay = matches[a].y - matches[a].x;
  by = matches[b].y - matches[b].x;
  if (ay > by){
    return 1;
  }
  else if (ay < by){
    return 2;
  }  
  else if (matches[a].start > matches[b].start){
    return 1;
  }
  else if (matches[a].start < matches[b].start){
    return 2;
  }
  else {
    return 0;
  }
}


/*----------------------------- bl_slClusterSop --------------------------------
 *
 * @brief       clustering of fragments for possible chains
 * @author      Christian Otto
 *
 */
void bl_slClusterSop(slmatch_t *fragments, Uint size, double eps, double lambda, int maxgap){
  Uint i, begin, length;
  int max_yst, min_yend, min_yst, max_yend, cor = 0;
  double max_score_per_pos;
  slmatch_t *current, *a, *b;

  /* nothing to do with size 0 */
  if (size == 0){
    return;
  }
  /* sorting (if already sorted -> only linear scan is required) */
  qsort(fragments, size, sizeof(slmatch_t), cmp_slmatch_qsort);

  /* initializations */
  begin = 0;
  a = fragments;
  max_score_per_pos = a->scr / (double) FLEN_Q(a);
  min_yst = FSTART_Q(a);
  max_yst = FSTART_Q(a);
  min_yend = FEND_Q(a);
  max_yend = FEND_Q(a);
  for (i = 0; i < size - 1; i++){
    b = fragments + (i+1);
    /* update start and end bounds */
    if (FSTART_Q(b) < min_yst){
      min_yst = FSTART_Q(b);
    }
    if (FSTART_Q(b) > max_yst){
      max_yst = FSTART_Q(b);
    }
    if (FEND_Q(b) < min_yend){
      min_yend = FEND_Q(b);
    }
    if (FEND_Q(b) > max_yend){
      max_yend = FEND_Q(b);
    }
    /* update cor (if necessary) and max_score_per_pos */
    if (b->scr / (double) FLEN_Q(b) > max_score_per_pos){
      max_score_per_pos = b->scr/(double) FLEN_Q(b);
    }
    if (lambda > eps && max_yst > min_yend){
      cor = (eps - lambda) * (double) (max_yst - min_yend);
    }
    /*
     * fragments are on difference reference sequences or
     * if gap between adjacent fragments on sequence is too large or
     * if dx >= dy and min(G_SOP) > max_score -> new cluster
     */
    if (a->subject != b->subject || MAXGAP2(b,a) ||
	((int) FSTART_S(b) > (int) FEND_S(a) &&
	 (int) FSTART_S(b) - (int) FEND_S(a) >= max_yst - min_yend &&
	 (double) lambda * (double) DX(b, a) + cor > (double) (max_yend - min_yst + 1) * max_score_per_pos)){
      length = i - begin + 1;
      /* no chaining necessary */
      if (length == 1){
	current = fragments + i;
	current->chain = (slchain_t *) malloc(sizeof(slchain_t));
	bl_slchainInitFromMatch(current->chain, current);
      }
      /* chain cluster */
      else {
	bl_slChainSop(fragments + begin, length, eps, lambda, maxgap);
      }
      begin = i + 1;
    }
    a = b;
  }
  length = i - begin + 1;
  /* no chaining necessary */
  if (length == 1){
    current = fragments + begin;
    current->chain = (slchain_t *) malloc(sizeof(slchain_t));
    bl_slchainInitFromMatch(current->chain, current);
  }
  /* chain cluster */
  else {
    bl_slChainSop(fragments + begin, length, eps, lambda, maxgap);
  }
}

/*------------------------------ bl_slChainSop ---------------------------------
 *    
 * @brief      creates chains of non overlapping fragments from a list of
 *	       fragments with the highest score using sum-of-pair gap costs
 * @parameter  fragments - the list of fragments to be chained,
 *	       eps - cost of aligning two anonymous characters,
 *	       lambda - cost of aligning an anonymous character with a gap
 *	       (it is demanded that lambda > 1/2 eps)
 * @author     Christian Otto
 *   
 */
void bl_slChainSop(slmatch_t *fragments, Uint size, double eps, double lambda, int maxgap){
  Uint i, **trans;
  PairSint t;
  point_t *point;
  slmatch_t *current, *first;
  slchain_t *chain, *chain2, *cand, **prev, *aprev, *bprev;
  RangeTree atrans, btrans;
  Container *points;
  VQueue refs;

  /* no chaining if size is zero */
  if (size == 0){
    return;
  }
  /* sorting (if already sorted -> only one scan required) */
  qsort(fragments, size, sizeof(slmatch_t), cmp_slmatch_qsort);
  /* preinitialize */
  points = bl_slExtractPoints(fragments, size);
  trans = bl_slGetTrans(points);
  /* initializations for chaining itself */
  prev = (slchain_t **) malloc(sizeof(slchain_t *) * size);
  memset(prev, 0, sizeof(slchain_t*) * size);
  bl_rangetreeInit(&atrans, 2, bl_containerSize(points), 
		   trans, sizeof(slchain_t *));
  bl_rangetreeInit(&btrans, 2, bl_containerSize(points), 
		   trans+2, sizeof(slchain_t *));
  bl_vqueueInit(&refs, 2 * size, sizeof(slchain_t *));
	       
  /* get t.a and t.b */
  t.a = ((point_t*) bl_containerGet(points,
				    trans[2][bl_containerSize(points)-1]))->x;
  t.b = ((point_t*) bl_containerGet(points,
				    trans[1][bl_containerSize(points)-1]))->y;
  for (i = 0; i < bl_containerSize(points); i++){
    point = (point_t *) bl_containerGet(points, i);		
    /* point is the start point of fragment point->index (in fragments) */
    if (point->start == 1){
      /* get current match from fragments */
      current = fragments + point->index;

      /* RMQ in octant O_1 and O_2 */
      aprev = (slchain_t *) bl_slChainSopRMQ(&atrans, point->trans[0],
					     point->trans[1], current,
					     eps, lambda, maxgap);
      bprev = (slchain_t *) bl_slChainSopRMQ(&btrans, point->trans[2],
					     point->trans[3], current,
					     eps, lambda, maxgap);

      /* choose the chaining candidate */
      if (aprev == NULL){
	prev[point->index] = bprev;
      }
      else if (bprev == NULL){
	prev[point->index] = aprev;
      }
      else {
	if (aprev->scr - (double) GSOP(current, aprev) >=
	    bprev->scr - (double) GSOP(current, bprev)){
	  prev[point->index] = aprev;
	}
	else {
	  prev[point->index] = bprev;
	}
      }

      /* addition due to local chaining (differ from Abouelhoda2004) */
      if (prev[point->index] != NULL){
	/* only take prev in range tree if it increases chain score */
	if (prev[point->index]->scr < (double) GSOP(current, prev[point->index])){
	  prev[point->index] = NULL;
	}
      }      
    }
    /* point is the end point of fragment point->index */
    else {
      current = fragments + point->index;
      chain = (slchain_t *) malloc(sizeof(slchain_t));
      /* 
       * chaining candidate found
       * (can be only better chain for prev or best chain for current as well)
       */
      if (prev[point->index] != NULL){
	cand = prev[point->index];
	bl_slchainInit(chain);
	chain->i = FSTART_Q(cand);
	chain->j = FEND_Q(current) - FSTART_Q(cand) + 1;
	chain->p = FSTART_S(cand);
	chain->q = FEND_S(current) - FSTART_S(cand) + 1;
	chain->scr = current->scr + cand->scr - (double) GSOP(current, cand);
	bl_containerMerge(chain->matches, cand->matches);
	bl_containerAdd(chain->matches, &current);
	if (current->subject != cand->subject){
	  DBG("slchain.c: Attempt to chain fragments of different \
reference sequences.\nExit forced.\n\n", NULL);
	  exit(-1);
	}
	chain->subject = current->subject;

	/* update best chain of cluster */
	first = *(slmatch_t **) bl_containerGet(chain->matches, 0);
	if (first->chain != NULL){
	  if (((slchain_t *) first->chain)->scr <= chain->scr){
	    slchain_t *tmp = (slchain_t *) malloc(sizeof(slchain_t));
	    bl_slchainInit(tmp);
	    bl_slchainCopy(tmp, chain);
	    bl_slchainpDestruct(&first->chain);
	    first->chain = (void *) tmp;
	  }
	}
      }
      /* no chaining candidate found */
      else {
	bl_slchainInitFromMatch(chain, current);	
	if (current->chain != NULL){
	  DBG("slchain.c: Error in updating of best chain. \
Exit forced.\n", NULL);
	  exit(-1);
	}
	else {
	  /* update best chain of current */
	  slchain_t *tmp = (slchain_t *) malloc(sizeof(slchain_t));
	  bl_slchainInit(tmp);
	  bl_slchainCopy(tmp, chain);
	  current->chain = tmp;
	}
      }
      bl_slChainSopActivate(&atrans, chain, point->trans[0],
			    point->trans[1], (double) GCSOP1(chain));
      chain2 = (slchain_t *) malloc(sizeof(slchain_t));
      bl_slchainInit(chain2);
      bl_slchainCopy(chain2, chain);
      bl_slChainSopActivate(&btrans, chain2, point->trans[2],
			    point->trans[3], (double) GCSOP2(chain2));
      bl_vqueueEnqueue(&refs, &chain);
      bl_vqueueEnqueue(&refs, &chain2);
    } /* END OF else */
  } /* END OF for (i = 0; i < bl_containerSize(points); i++) */
  bl_memusage();
  free(prev);
  bl_vqueueDestruct(&refs, bl_slchainpDestruct);
    /* destruct range trees */
  bl_rangetreeDestruct(&atrans, 2, NULL);
  bl_rangetreeDestruct(&btrans, 2, NULL);
  /* destruct transformations */
  for(i = 0; i < 4; i++){
    free(trans[i]);
  }
  free(trans);
  /* destruct points */
  for(i = 0; i < bl_containerSize(points); i++){
    point_t *tmp = (point_t*) bl_containerGet(points, i);
    if (tmp->trans != NULL){
      free(tmp->trans);
    }
  }
  bl_containerDestruct(points, NULL);
  free(points);

  return;
}



/*---------------------------- blChainSopRMQ -----------------------------------
 *    
 * @brief 	helper method for chaining with sum-of-pair gap costs,
 *              does RMQ in the range tree, currently without use of
 *              maximum field for a node or fractional cascading,
 *              algorithm like in the paper of Abouelhoda et al.
 * @return      result of two-dimensional RMQ in the range tree as
 *              slchain_t*
 * @author 	Christian Otto
 *   
 */
void* bl_slChainSopRMQ(void *data, Uint x, Uint y, slmatch_t *current, double eps, double lambda, int maxgap){
  int pred;
  double resprio;
  slchain_t *res = NULL, **tmp, *chain, *local;
  slmatch_t **first;
  RangeTree *trans = (RangeTree *) data;
  rtnode_t *node, *left;
  resprio = (-1) * DBL_MAX;
  node = trans->root;
  while (node != NULL){
    /* 
     * go right in the rangetree and check the highest value
     * in the left subtree
     */
    if (node->value <= x){
      /* the current node is not a leaf */
      if (node->left != NULL){
	left = (rtnode_t *) node->left;
	#ifndef BINTREE
	pred = bl_vebtreePred((VebTree *) left->assoc, y);
	tmp = (slchain_t **) bl_vebtreeGetData((VebTree *) left->assoc, pred);
	#else
	pred = bl_bintreePred((BinTree *) left->assoc, y);
	tmp = (slchain_t **) bl_bintreeGet((BinTree *) left->assoc, pred);
	#endif
      }
      else {
	#ifndef BINTREE
	pred = bl_vebtreePred((VebTree *) node->assoc, y);
	tmp = (slchain_t **) bl_vebtreeGetData((VebTree *) node->assoc, pred);
	#else
	pred = bl_bintreePred((BinTree *) node->assoc, y);
	tmp = (slchain_t **) bl_bintreeGet((BinTree *) node->assoc, pred);
	#endif
      }
      /*
       * check if current fragment could extend the predecessor chain
       * (even if it is will not be the best preecessor of the current fragment)
       * NOT CONSIDERED: not direct predecessors of current fragment
       */
      if (tmp != NULL && current->scr >= (double) GSOP(current, *tmp) && !MAXGAP(current, *tmp)){
	first = (slmatch_t **) bl_containerGet((*tmp)->matches, 0);
	local = (slchain_t *) (*first)->chain;
	/* 
	 * only extend if chain score stored in first fragment can be increased
	 * by use of the current fragment
	 */
	if (local->scr < current->scr + (*tmp)->scr - (double) GSOP(current, *tmp)){
	  chain = (slchain_t *) malloc(sizeof(slchain_t));
	  bl_slchainInit(chain);
	  bl_slchainCopy(chain, *tmp);
	  chain->j = FEND_Q(current) - FSTART_Q(*tmp) + 1;
	  chain->q = FEND_S(current) - FSTART_S(*tmp) + 1;
	  chain->scr += current->scr - (double) GSOP(current, *tmp);
	  bl_containerAdd(chain->matches, &current);
	  bl_slchainpDestruct(&local);
	  (*first)->chain = chain;
	}
      }
      /* check and maybe update max score */
      if (tmp != NULL && (*tmp)->prio > resprio && !MAXGAP(current, *tmp)){
	res = *tmp;
	resprio = (*tmp)->prio;
      }
      node = (rtnode_t *) node->right;
    }
    /* just go left */
    else {
      node = (rtnode_t *) node->left;
    }
  }
  return res;
}

/*-------------------------- bl_slChainSopActivate -----------------------------
 *    
 * @brief 	helper method for chaining with sum-of-pair gap costs,
 *              activation inserts a new chain in the range tree 
 *              and cleans out succ chains with a lower priority
 *		compared to the current chain
 * @author 	Christian Otto
 *   
 */
void bl_slChainSopActivate(void *data, slchain_t *current,
			   Uint x, Uint y, double gc){
  int pred, succ;
  slchain_t **predmatch, **succmatch;
  #ifndef BINTREE
  VebTree *assoc;
  #else
  BinTree *assoc;
  #endif
  RangeTree *trans = (RangeTree *) data;
  rtnode_t *node = trans->root;
  /* get priority of current chain */
  current->prio = current->scr - gc;
  while (node != NULL){
    assoc = node->assoc;    
    /* 
     * needed to search for y + 1 to find
     * if an element with y is already in it
     * (always pred(x) < x)
     */
    #ifndef BINTREE
    pred = bl_vebtreePred(assoc, y + 1);
    predmatch = (slchain_t **) bl_vebtreeGetData(assoc, pred);
    /* could also be > instead of >= (depends on purpose) */
    if (predmatch == NULL || current->prio >= (*predmatch)->prio){
      bl_vebtreeInsert(assoc, y, &current, NULL);
      succ = bl_vebtreeSucc(assoc, y);
      succmatch = (slchain_t **) bl_vebtreeGetData(assoc, succ);      
      while(succmatch != NULL && current->prio > (*succmatch)->prio){
	free((slchain_t **) bl_vebtreeDelete(assoc, succ, NULL));
	succ = bl_vebtreeSucc(assoc, succ);
	succmatch = (slchain_t **) bl_vebtreeGetData(assoc, succ);
      }
    }
    #else
    pred = bl_bintreePred(assoc, y + 1);
    predmatch = (slchain_t **) bl_bintreeGet(assoc, pred);
    /* could also be > instead of >= (depends on purpose) */
    if (predmatch == NULL || current->prio >= (*predmatch)->prio){
      bl_bintreeInsert(assoc, y, &current, NULL);
      succ = bl_bintreeSucc(assoc, y);
      succmatch = (slchain_t **) bl_bintreeGet(assoc, succ);      
      while(succmatch != NULL && current->prio > (*succmatch)->prio){
	free((slchain_t **) bl_bintreeDelete(assoc, succ, NULL));
	succ = bl_bintreeSucc(assoc, succ);
	succmatch = (slchain_t **) bl_bintreeGet(assoc, succ);
      }
    }
    #endif
    /* go right */
    if (node->value < x){
      node = (rtnode_t *) node->right;
    }
    /* go left */
    else {
      node = (rtnode_t *) node->left;
    }
  }
}
