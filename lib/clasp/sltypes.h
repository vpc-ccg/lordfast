#ifndef SLTYPES_H
#define SLTYPES_H

/**
 * sltypes.h
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
 * Id: $Id: sltypes.h 114 2010-04-19 13:59:41Z steve $
 * Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/sltypes.h $ 
 */

#include <stdlib.h>
#include <stdio.h>
#include "basic-types.h"
#include "container.h"

/*
 * Typedef
 */
typedef struct {
  Uint *val;
  Uint size;
} slcount_t;

typedef struct {
  /* qry_begin: i, qry_len: j */
  int i;
  int j;
  /* seq_begin: p, seq_len: q */
  // int p;
  // int q;
  long p;
  long q;
  
  double scr;
  double evalue;
  slcount_t mat;
  slcount_t mis;
  Uint ins;
  Uint del;
  void *chain;
  Uint subject;
  Uint idx;
} slmatch_t;

/* 
 * slchain_t represents the collection of multiple
 * slmatch_t*, its borders, its score and evalue
 */
typedef struct {
  /* 
   * qry_begin: i, qry_len: j
   * (i is qry_begin of first match and
   *  j is length of whole global match on qry)
   */
  int i;
  int j;
  /*
   * seq_begin: p, seq_len: q
   * (p is seq_begin of first match and
   *  q is length of whole global match on seq)
   */
  // int p;
  // int q;
  long p;
  long q;
  double scr;
  double prio;
  double evalue;
  Container *matches;
  Uint subject;
} slchain_t;

void bl_slmatchInit(slmatch_t *sl, Uint size);
void bl_slmatchDestruct(void * obj);
void bl_slcountInit(slcount_t *cnt, Uint size);
void bl_slcountDestruct(void *obj);
void bl_slcountSet(slcount_t *cnt, Uint val);
void bl_slcountCopy(slcount_t *dest, slcount_t *src);
Uint bl_slcountCmp(slcount_t *a, slcount_t *b);
void bl_slchainInit(slchain_t *);
void bl_slchainInitFromMatch(slchain_t *chain, slmatch_t *match);
void bl_slchainCopy(slchain_t *dest, slchain_t *src);
void bl_slchainDestruct(void *);
void bl_slchainpDestruct(void *);
void bl_memusage();

#endif /* SLTYPES_H */
