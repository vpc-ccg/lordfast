#ifndef SLCHAIN_H
#define SLCHAIN_H

/**
 * slchain.h
 * chaining of non overlapping fragments
 * with different possible gap costs
 *
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Wed Jun  4 16:13:58 CEST 2008
 */

/*
 * SVN
 * Revision of last commit: $Rev: 115 $
 * Author: $Author: steve $
 * Date: $Date: 2010-04-21 11:58:46 +0200 (Wed, 21 Apr 2010) $
 * Id: $Id: slchain.h 115 2010-04-21 09:58:46Z steve $
 * Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/slchain.h $
 */

#include <stdio.h>
#include <stdlib.h>
#include "basic-types.h"
#include "container.h"

/* Macros */
#define FSTART_S(a) ((a)->p)
#define FEND_S(a)   ((a)->p+(a)->q-1)
#define FLEN_S(a)   ((a)->q)

#define FSTART_Q(a) ((a)->i)
#define FEND_Q(a)   ((a)->i+(a)->j-1)
#define FLEN_Q(a)   ((a)->j)

/* delta(a,b) = |a - b| (now ||a - b| - 1|) */
//#define D(a, b)         ((a > b)?(a - b - 1):(b - a - 1))
#define D(a, b)         ((a > b)?(a - b - 1):(b > a)?(b - a - 1):1)
/* DX on sequence, DY on qry */
#define DX(fprim, f)    (D(FSTART_S(fprim), FEND_S(f)))
#define DY(fprim, f)    (D(FSTART_Q(fprim), FEND_Q(f)))
/* SOP gap costs */
#define GSOP(fprim, f)   ((DX(fprim, f) >= DY(fprim, f))? \
			 (lambda * DX(fprim, f) + (eps - lambda) * DY(fprim, f)): \
			 (lambda * DY(fprim, f) + (eps - lambda) * DX(fprim, f)))
/* simple gap costs (L1 metric) */
#define GLIN(fprim, f)    (lambda * DX(fprim, f) + eps * DY(fprim, f))
/* geometric costs */
#define GCSOP1(f)	  (lambda * D(t.a, FEND_S(f)) + (eps - lambda) * D(t.b, FEND_Q(f)))
#define GCSOP2(f)	  (lambda * D(t.b, FEND_Q(f)) + (eps - lambda) * D(t.a, FEND_S(f)))
#define GCLIN(f)	  (lambda * D(t.a, FEND_S(f)) + eps * D(t.b, FEND_Q(f)))
/* check maxgap */
#define MAXGAP(fprim, f)  ((maxgap != -1) && ((DX(fprim, f) > maxgap) || DY(fprim, f) > maxgap))
#define MAXGAP2(fprim, f) ((maxgap != -1) && (DX(fprim, f) > maxgap))


/* Typedef */
/*
 * point_t is a temporary data structure used
 * for the chaining with all information included
 * which is needed for it
 */
typedef struct {
  int x;
  int y;
  Uint* trans;
  int size;
  Uint index;
  BOOL start;
} point_t;

/* Prototypes */
Container* bl_slExtractPoints(slmatch_t *fragments, Uint size);
Uint** bl_slGetTrans(Container *points);
void bl_slClusterLin(slmatch_t *fragments, Uint size, double eps, double lambda, int maxgap);
void bl_slChainLin(slmatch_t *fragments, Uint size, double eps, double lambda, int maxgap);
void bl_slClusterSop(slmatch_t *fragments, Uint size, double eps, double lambda, int maxgap);
void bl_slChainSop(slmatch_t *fragments, Uint size, double eps, double lambda, int maxgap);
void* bl_slChainSopRMQ(void *data, Uint x, Uint y, slmatch_t *current,
		       double eps, double lambda, int maxgap);
void  bl_slChainSopActivate(void *queue, slchain_t *sc,
			    Uint x, Uint y, double gc);
Uint cmp_slmatch_trans_first_x (Uint a, Uint b, void *data, void *info);
Uint cmp_slmatch_trans_first_y (Uint a, Uint b, void *data, void *info);
Uint cmp_slmatch_trans_second_x (Uint a, Uint b, void *data, void *info);
Uint cmp_slmatch_trans_second_y (Uint a, Uint b, void *data, void *info);
Uint cmp_slmatch_end_quick (Uint a, Uint b, void *adata, void *bdata);
int cmp_slmatch_qsort(const void *a, const void *b);

#endif /* SLCHAIN_H */
