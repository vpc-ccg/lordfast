#ifndef SORT_H
#define SORT_H

/*
 * sort.h
 * declarations for various sorting algorithms
 *
 * @author Steve Hoffmann
 * @date Mon 27 Nov 2006
 *
 *  SVN
 *  Revision of last commit: $Rev: 110 $
 *  Author: $Author: steve $
 *  Date: $Date: 2009-06-03 11:45:13 +0200 (Wed, 03 Jun 2009) $
 *
 *  Id: $Id: sort.h 110 2009-06-03 09:45:13Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/sort.h $
 */
 #include "basic-types.h"


 Uint cmp_dbl(Uint, Uint, void *, void *);
 Uint cmp_flt(Uint, Uint, void *, void *);
 Uint cmp_int(Uint, Uint, void *, void *);
 Uint cmp_Uint(Uint, Uint, void *, void *);
 Uint cmp_int_bin (Uint a, void *, void *, void *);
 Uint cmp_Uint_bin (Uint, void *, void *, void *);
 int cmp_Uint_qsort(const void *a, const void *b);
 int cmp_int_qsort(const void *a, const void *b);
 Uint binarySearch_m(void *, Uint, void*, Uint(*cmp)(Uint,void *, void *, void *), void *);
 Uint binarySearch(void *, Uint, void*, Uint(*cmp)(Uint,void *, void *, void *), void *);
 Uint *quickSort(void *, void *, Uint, Uint (*cmp)(Uint, Uint, void *, void *), void *);
 Uint compareMkstr(Uint, Uint, Uint depth, void *, void *);
 Uint compareMkstrptr(Uint, Uint, Uint depth, void *, void *);
 Uint *quickSortMultikey(void *, void*, Uint, Uint (*cmp)(Uint, Uint, Uint, void*, void*), Uint, void*);
 int cmp_PairUint_qsort(const void*, const void*);
 int cmp_PairUint_bsearch(const void*, const void*);
 int cmp_PairSint_qsort(const void*, const void*);
 int cmp_PairSint_bsearch(const void*, const void*);

#endif

