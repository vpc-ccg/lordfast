/*
 * clasp is written by Christian Otto, Bioinformatics, University of Leipzig
 * https://www.bioinf.uni-leipzig.de/Software/clasp/
 *
 * Dynamic programming chaining is written by Ehsan Haghshenas (ehaghshe AT sfu DOT ca)
 */

#ifndef CLASP_H
#define CLASP_H

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
