#ifndef CLASP_H
#define CLASP_H

/**
 * clasp.h
 * fast fragment chaining
 * using sop gap costs
 * 
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Mon Nov  2 10:06:11 CET 2009
 */

/*
 * SVN
 * Revision of last commit: $Rev: 116 $
 * Author: $Author: steve $
 * Date: $Date: 2010-06-30 13:51:27 +0200 (Wed, 30 Jun 2010) $
 * Id: $Id: clasp.h 116 2010-06-30 11:51:27Z steve $
 * Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/src/clasp.h $
 */

#include <stdio.h>
#include <stdlib.h>
#include "container.h"
#include "basic-types.h"

#include "LordFAST.h"

#define SOP     ((unsigned char) (0 << 0))
#define LIN     ((unsigned char) (1 << 0))
#define VERSION "1.1"

/* Typedef */
typedef struct {
  char *infilename;
  char *outfilename;
  FILE *dev;
  Container *fragments;
  Container *lines;
  Container *subject;
  unsigned char chainmode;
  double lambda;
  double epsilon;
  double minscore;
  int maxgap;
  Uint minfrag;
  Uint colnum;
  Uint* colorder;
  Uint* idcol;
  int idcolnum;
  BOOL outputc;
  BOOL outputf;
  BOOL outputm;
  BOOL outputorig;
} claspinfo_t;

// #ifdef __cplusplus
// extern "C" {
// #endif

int chain_seeds_clasp(Seed_t *fragment_list, uint32_t nFragment, Chain_t &bestChain);
void chain_seeds_n2(Seed_t *fragment_list, uint32_t nFragment, Chain_t &bestChain);

// #ifdef __cplusplus
// }
// #endif

#endif /* CLASP_H */
