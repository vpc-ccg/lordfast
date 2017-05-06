/**
 * sltypes.c
 * types for several slmodules
 *
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Wed Jun  4 16:13:58 CEST 2008
 */

/* 
 * SVN
 * Revision of last commit: $Rev: 114 $
 * Author: $Author: steve $
 * Date: $Date: 2010-04-19 15:59:41 +0200 (Mon, 19 Apr 2010) $
 * Id: $Id: sltypes.c 114 2010-04-19 13:59:41Z steve $
 * Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/sltypes.c $ 
 */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "sltypes.h"
#include "container.h"
#include "debug.h"

double maxmem = 0;

/*----------------------------- bl_slmatchInit ---------------------------------
 *    
 * @brief 	helper function to initialize slmatch_t
 * @author 	Christian Otto
 *   
 */
void bl_slmatchInit(slmatch_t *sl, Uint size){
  sl->i = 0;
  sl->j =0;
  sl->p = 0;
  sl->q = 0;
  sl->scr = 0;
  sl->evalue = -1;
  bl_slcountInit(&sl->mat, size);
  bl_slcountSet(&sl->mat, 0);
  bl_slcountInit(&sl->mis, size);
  bl_slcountSet(&sl->mis, 0);
  sl->ins = 0;
  sl->del = 0;
  sl->chain = NULL;
  sl->subject = 0;
  sl->idx = 0;
}

/*---------------------------- bl_slmatchDestruct ------------------------------
 *    
 * @brief 	helper function to initialize slmatch_t
 * @author 	Christian Otto
 *   
 */
void bl_slmatchDestruct(void *obj){
  slmatch_t *sl = (slmatch_t *) obj;
  sl->i = 0;
  sl->j = 0;
  sl->p = 0;
  sl->q = 0;
  sl->scr = 0;
  sl->evalue = 0;
  bl_slcountDestruct(&sl->mat);
  bl_slcountDestruct(&sl->mis);
  sl->ins = 0;
  sl->del = 0;
  sl->chain = NULL;
  sl->subject = 0;
  sl->idx = 0;
}

/*------------------------------ bl_slcountInit --------------------------------
 *    
 * @brief 	helper function to initialize slcount_t
 * @author 	Christian Otto
 *   
 */
void bl_slcountInit(slcount_t *cnt, Uint size){
  cnt->size = size;
  cnt->val = (Uint *) malloc(sizeof(Uint) * size);
}

/*------------------------------ bl_slcountInit --------------------------------
 *    
 * @brief 	helper function to initialize slcount_t
 * @author 	Christian Otto
 *   
 */
void bl_slcountDestruct(void *obj){
  slcount_t *cnt = (slcount_t *) obj;
  free(cnt->val);
  cnt->size = 0;
}

/*------------------------------- bl_slcountSet --------------------------------
 *    
 * @brief 	helper function to set slcount_t to specific value
 * @author 	Christian Otto
 *   
 */
void bl_slcountSet(slcount_t *cnt, Uint val){
  Uint i;
  for (i = 0; i < cnt->size; i++){
    cnt->val[i] = val;
  }
}

/*------------------------------- bl_slcountCopy -------------------------------
 *    
 * @brief 	helper function to copy values of slcount_t
 * @author 	Christian Otto
 *   
 */
void bl_slcountCopy(slcount_t *dest, slcount_t *src){
  Uint i;
  if (dest->size != src->size){
    DBG("seachsl.c: Copy between different sized slcount_t not \
possible. Exit forced.\n", NULL);
    exit(-1);
  }
  for (i = 0; i < src->size; i++){
    dest->val[i] = src->val[i];
  }
}

/*------------------------------- bl_slcountCmp --------------------------------
 *    
 * @brief 	helper function to compare slcount_t
 *              0 if a < b for every level,
 *              1 if a <= b for every level (at least a == b for one level),
 *              2 otherwise
 * @author 	Christian Otto
 *   
 */
Uint bl_slcountCmp(slcount_t *a, slcount_t *b){
  Uint i, ret = 0;
  if (a->size != b->size){
    DBG("seachsl.c: Comparison between different sized '%d' and '%d' slcount_t \
not possible. Exit forced.\n", a->size, b->size);
    exit(-1);
  }
  for (i = 0; i < a->size; i++){
    if (a->val[i] > b->val[i]){
      return 2;
    }
    else if (a->val[i] == b->val[i]){
      ret = 1;
    }
  }
  return ret;
}

/*------------------------------ bl_slchainInit --------------------------------
 *    
 * @brief 	helper function to initialize slchain_t
 * @author 	Christian Otto
 *   
 */
void bl_slchainInit(slchain_t *sc){
  sc->i = 0;
  sc->j = 0;
  sc->p = 0;
  sc->q = 0;
  sc->scr = 0;
  sc->prio = 0;
  sc->evalue = -1;
  sc->matches = (Container *) malloc(sizeof(Container));
  bl_containerInit(sc->matches, 10, sizeof(slmatch_t*));
  sc->subject = 0;
}

/*-------------------------- bl_slchainInitFromMatch ---------------------------
 *    
 * @brief 	helper function to initialize slchain_t from slmatch_t
 * @author 	Christian Otto
 *   
 */
void bl_slchainInitFromMatch(slchain_t *chain, slmatch_t *match){
  chain->i = match->i;
  chain->j = match->j;
  chain->p = match->p;
  chain->q = match->q;
  chain->scr = match->scr;
  chain->evalue = match->evalue;
  chain->matches = (Container *) malloc(sizeof(Container));
  bl_containerInit(chain->matches, 1, sizeof(slmatch_t*));
  bl_containerAdd(chain->matches, &match);
  chain->subject = match->subject;
}

/*------------------------------ bl_slchainCopy --------------------------------
 *    
 * @brief 	helper function to copy slchain_t
 * @author 	Christian Otto
 *   
 */
void bl_slchainCopy(slchain_t *dest, slchain_t *src){
  dest->i = src->i;
  dest->j = src->j;
  dest->p = src->p;
  dest->q = src->q;
  dest->scr = src->scr;
  bl_containerMerge(dest->matches, src->matches);
  dest->subject = src->subject;
}

/*---------------------------- bl_slchainDestruct ------------------------------
 *    
 * @brief 	helper function to destruct slchain_t
 * @author 	Christian Otto
 *   
 */
void bl_slchainDestruct(void *obj){
  slchain_t *sc = (slchain_t *) obj;
  sc->i = 0;
  sc->j = 0;
  sc->p = 0;
  sc->q = 0;
  sc->scr = 0;
  sc->prio = 0;
  sc->evalue = 0;
  if (sc->matches != NULL){
    bl_containerDestruct(sc->matches, NULL);
    free(sc->matches);
  }
  sc->matches = NULL;
  sc->subject = 0;
}

/*--------------------------- bl_slchainpDestruct ------------------------------
 *    
 * @brief 	helper function to destruct slchain_t pointer
 * @author 	Christian Otto
 *   
 */
void bl_slchainpDestruct(void *obj){
  slchain_t **sc = (slchain_t **) obj;
  if (sc != NULL){
    bl_slchainDestruct(*sc);
    free(*sc);
  }
}

/*----------------------------- bl_memusage ------------------------------------
 *
 * @brief       returns current memory usage for profiling
 *              if macro is defined
 * @author      Christian Otto
 *
 */
void bl_memusage(){
  #ifdef MAXMEM
  char buf[30];
  snprintf(buf, 30, "/proc/%u/statm", (unsigned)getpid());
  FILE* pf = fopen(buf, "r");
  if (pf) {
    unsigned size;       //  total program size (virtual)
    unsigned resident;   //  resident set size
    fscanf(pf, "%u %u", &size, &resident);
    if (maxmem < size){
      maxmem = size;
    }
  }
  fclose(pf);
  #endif
}
