
/*
 * sort.c
 * implementation of various sorting algorithms
 *
 * @author Steve Hoffmann
 * @date Mon 27 Nov 2006
 *
 *  SVN
 *  Revision of last commit: $Rev: 110 $
 *  Author: $Author: steve $
 *  Date: $Date: 2009-06-03 11:45:13 +0200 (Wed, 03 Jun 2009) $
 *
 *  Id: $Id: sort.c 110 2009-06-03 09:45:13Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/sort.c $
 */

 #include <math.h>
 #include <string.h>
 #include "basic-types.h"
 #include "memory.h"
 #include "vstack.h"
 #include "mathematics.h"
 #include "sort.h"
 #include "debug.h"



 Uint cmp_dbl(Uint a, Uint x, void *data, void *info) {
	double *d = (double*) data;

	/*if(floor(d[a]) > floor(d[x])) return 1;
	if(floor(d[a]) < floor(d[x])) return 2;
*/
	if (fabs((double) d[a] - d[x]) <= FLT_EPSILON) return 0;
	if ((double) d[a] - d[x] >  FLT_EPSILON) return 1;
	if ((double) d[x] - d[a] >  FLT_EPSILON) return 2;

	return 0;
 }


Uint cmp_flt(Uint a, Uint x, void *data, void *info) {
	float *d = (float*) data;

	if (d[a]>d[x]) return 1;
	if (d[a]<d[x]) return 2;

	return 0;
 }
 
 Uint cmp_int(Uint a, Uint x, void *data, void *info) {
	int *d = (int*) data;

	if (d[a]>d[x]) return 1;
	if (d[a]<d[x]) return 2;

	return 0;
 }

 Uint cmp_Uint(Uint a, Uint x, void *data, void *info) {
	Uint *d = (Uint*) data;

	if (d[a]>d[x]) return 1;
	if (d[a]<d[x]) return 2;

	return 0;
 }

 Uint cmp_int_bin(Uint a, void *data, void *key, void *info) {
	int *d = (int*) data;
	int *k = (int*) key;

	if (d[a]>*k) return 1;
	if (d[a]<*k) return 2;

	return 0;
 }


 Uint cmp_Uint_bin(Uint a, void *data, void *key, void *info) {
	Uint *d = (Uint*) data;
	Uint *k = (Uint*) key;

	if (d[a]>*k) return 1;
	if (d[a]<*k) return 2;

	return 0;
 }

int cmp_Uint_qsort(const void *a, const void *b) {
  Uint  *first = (Uint*) a;
  Uint  *secnd = (Uint*) b;

  if (*first > *secnd) return 1;
  if (*first < *secnd) return -1;

  return 0;

}

 
int cmp_int_qsort(const void *a, const void *b) {
  int  *first = (int*) a;
  int  *secnd = (int*) b;

  if (*first > *secnd) return 1;
  if (*first < *secnd) return -1;

  return 0;

}

Uint binarySearch_m(void *toSearch, Uint size, void *key, 
 					Uint (*cmp)(Uint, void *, void *, void *), 
					void *info) {
	int left=0, right=size, mid=left, res;
	
	while (left<=right) {
		mid = (left+right)/2;
		res = cmp(mid, toSearch, key, info);
		
		if (res == 1) 
			right = mid-1;
		else 
		if (res == 2) 
			left = mid+1;
		else 
			return mid;
        
        if (left == right) {
          if (res == 2)
            return mid;
          else
            return mid-1;
        }
	}

    return mid;
 }

 Uint binarySearch(void *toSearch, Uint size, void *key, 
 					Uint (*cmp)(Uint, void *, void *, void *), 
					void *info) {
	int left=0, right=size, mid, res;
	
	while (left<=right) {
		mid = (left+right)/2;
		res = cmp(mid, toSearch, key, info);
		
		if (res == 1) 
			right = mid-1;
		else 
		if (res == 2) 
			left = mid+1;
		else 
			return mid;
	}

    return size+1;
 }


Uint *quickSort(void *space, void* toSort, Uint size, 
		Uint (*cmp)(Uint, Uint, void *, void*),
		void *info) {
  int left, left2, right, right2;
  PairSint ins, *lr;
  Uint i, resc, *sorted, x;
  VStack vstack;
	
  sorted = ALLOCMEMORY(space, NULL, Uint, size);
  for (i=0; i < size; i++) sorted[i]=i;
  ins.a = 0;
  ins.b = size-1;
  bl_vstackInit(&vstack, 10000, sizeof(PairSint));
  bl_vstackPush(&vstack, &ins);
   
  while (!bl_vstackIsEmpty(&vstack)){
    lr = (PairSint *) bl_vstackPop(&vstack, NULL);
    left = lr->a;
    right = lr->b;
    free(lr);
    while (left < right) {
      x=sorted[(left+right)/2];
      left2  = left;
      right2 = right;
	
      do {
	while(cmp(sorted[left2],  x, toSort, info)==2){	
	  left2++;
	}
	while(cmp(sorted[right2], x, toSort, info)==1){ 
	  right2--;
	}
			
	if(left2 <= right2) {
	  resc = sorted[right2];
	  sorted[right2]=sorted[left2];
	  sorted[left2]=resc;
	  left2++;
	  right2--;
	} 	
      } while (right2 >= left2);
			

      if ((left2-left) > (right-left2))  {		
	ins.a = left;
	ins.b = right2;
	bl_vstackPush(&vstack, &ins);
	left  = left2;
      } else {
	ins.a = left2;
	ins.b = right;
	bl_vstackPush(&vstack, &ins);
	right = right2;
      }
    }
  }
  bl_vstackDestruct(&vstack, NULL);
  return sorted;
}

Uint 
 compareMkstrptr(Uint a, Uint b, Uint depth, 
	 					void *data, void* info)
{
	char **ptr = (char**) data;
	

	if (ptr[a][depth] > ptr[b][depth]) return 1;
	if (ptr[a][depth] < ptr[b][depth]) return 2;

	return 0;
}


Uint 
 compareMkstr(Uint a, Uint b, Uint depth, 
	 					void *data, void* info)
{
	char *ptr = (char*) data;
	

	if (ptr[a+depth] > ptr[b+depth]) return 1;
	if (ptr[a+depth] < ptr[b+depth]) return 2;

	return 0;
}

/*--------------------------------- vecswap ----------------------------------
 *    
 * swaps a vector (needed by quickSortMultiKey)
 * 
 */
 
void
vecswap(int i, int j, int n, Uint *x)
{
  	while(n-- > 0) {
		SWAPUINT(x, i, j);
		i++;
		j++;
	}
}


/*---------------------------- quickSortMultikey -----------------------------
 *    
 * performs a mulitkey quicksort Bentley Sedgewick style
 * 
 * Implementation of programm ssort1, as described in: 
 * Fast Algorithms for Sorting and Searching Strings
 * Proc. of the ACM-SIAM Symposium on Discrete Algorithms, 
 * pages 360-369. http://www.cs.princeton.edu/~rs/strings/
 * 
 */

Uint *
quickSortMultikey (void *space, void* toSort, Uint size, 
		   Uint (*cmp)(Uint, Uint, Uint, void *, void*),
		   Uint sentinel, void *info) 
{
  Sint a, b, c, d, v, n, r;
  TripleSint ins;
  Uint *sorted = NULL, offset;
  Uint depth = 0;
  VStack vstack;
	
	
  if (size == 0) return NULL;

  sorted = ALLOCMEMORY(space, NULL, Uint, size);
  if (size<=1) {
    sorted[0]=0;
  }
    
  for (r=0; r < size; r++) sorted[r]=r;	
  bl_vstackInit(&vstack, 100, sizeof(TripleSint));
  n = size;
  offset=0;
	
  while (1) {
    a = rand() % n;
    SWAPUINT(sorted, offset, a+offset);
    v = sorted[offset];
    a = b = 1;
    c = d = n-1;
	
    for(;;) {
	    
      while(b<=c&&((r=cmp(sorted[b+offset],v,depth,toSort,info))==2||r==0))
	{
		  	
	  if (r==0) {
	    SWAPUINT(sorted, a+offset, b+offset);
	    a++;
	  }
	  b++;
	}
      while(b<=c&&((r=cmp(sorted[c+offset],v,depth,toSort,info))==1||r==0)) 
	{
		  		
		  
	  if (r==0) {
	    SWAPUINT(sorted, c+offset, d+offset);
	    d--;
	  }
	  c--;
	}
      if (b > c) break;
      SWAPUINT(sorted, b+offset, c+offset);
      b++;
      c--;
    }	
    r = MIN(a, (b-a));
    vecswap(offset, (b-r)+offset, r, sorted);
    r = MIN((d-c), (n-d-1));
    vecswap(b+offset, (n-r)+offset, r, sorted);
    /*sort lesser*/
    r = b-a;
    if (r > 1) {
      ins.a = offset;
      ins.b = r;
      ins.c = depth;
      bl_vstackPush(&vstack, &ins);
    }
    /*sort equal*/
    if ((a+n-d-1) > 1 && cmp(sorted[r+offset], sentinel, depth, toSort, info) != 0)
      /*if (r > 1 && sorted[r+offset]!=sentinel)*/
      { 
	ins.a = r+offset;
	ins.b = (a+n-d-1);
	ins.c = depth+1;
	bl_vstackPush(&vstack, &ins);
      }
    /*sort greater*/
    r=d-c;
    if (r > 1) {
      ins.a = (n-r)+offset;
      ins.b = r;
      ins.c = depth;
      bl_vstackPush(&vstack, &ins);
    }
		
    if (!bl_vstackIsEmpty(&vstack)){
      ins = *((TripleSint *) bl_vstackPop(&vstack, NULL));
      offset = ins.a;
      r = ins.b;
      depth = ins.c;
    } else {
      break;
    }	
  }   
  bl_vstackDestruct(&vstack, NULL);
  return sorted;
}

int
cmp_PairUint_qsort(const void* a, const void* b) {
    PairUint *i;
    PairUint *j;

    i = (PairUint*) a;
    j = (PairUint*) b;

    if (i->a > j->a) return 1;
    if (i->a < j->a) return -1;

    return 0;
}

int
cmp_PairUint_bsearch(const void* a, const void* b) {
    Uint     *key;
    PairUint *obj;

    key = (Uint*) a;
    obj = (PairUint*) b;

    if (*key > obj->a) return 1;
    if (*key < obj->a) return -1;

    return 0;
}

int
cmp_PairSint_qsort(const void* a, const void* b) {
    PairSint *i;
    PairSint *j;

    i = (PairSint*) a;
    j = (PairSint*) b;

    if (i->a > j->a) return 1;
    if (i->a < j->a) return -1;

    return 0;
}

int
cmp_PairSint_bsearch(const void* a, const void* b) {
    int     *key;
    PairSint *obj;

    key = (int*) a;
    obj = (PairSint*) b;

    if (*key > obj->a) return 1;
    if (*key < obj->a) return -1;

    return 0;
}

