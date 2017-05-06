#ifndef _MATHEMATICS_H_
#define _MATHEMATICS_H_

#include <stdlib.h>
#include <stdio.h>
#include "basic-types.h"


#define MATRIX2D(X,NCOL,I,J) 			X[(I)*(NCOL)+(J)]
#define MATRIX3D(X,DIM_M,DIM_N,I,J,K) 	X[((I)*(DIM_M)+(J))*DIM_N+(K)]
#define VECTOR(X,I) 					((X)->elements[I])


#define INITMATRIX(X,SIZE,TYPESIZE) 	initArray(X,SIZE,TYPESIZE)
#define INITMATRIX2D(X,M,N,TYPESIZE) 	initArray(X, (M)*(N), TYPESIZE)
#define INITMATRIX3D(X,M,N,L,TYPESIZE) 	initArray(X, (M)*(N)*(L), TYPESIZE)
#define INITVECTOR(V) 					(V)->elements=NULL;\
										(V)->length=0

#define RESIZEVEC(V,N)					initArray(V,N,sizeof(vectorelem))
#define LENGTHVEC(V) 					((V)->length)
#define SWEEPVEC(V)						{int m;\
										for (m=0;m<LENGTHVEC(V);m++)\
										VECTOR(V,mi)=0;}
#define APPENDVEC(S,V,E) 				appendvector(S,V,E)
#define SWAPVEC(I,J,V) 	  				{ vectorelem msv; \
						  				msv = (V)->elements[I];\
						 				(V)->elements[I] = (V)->elements[J];\
						  				(V)->elements[J] = msv; }
#define DESTRUCTVEC(S,V)                destructVector((S),(V));
#define SWAPUINT(X,A,B)					{Uint msu ;\
  										 msu= X[A]; \
										 X[A] = X[B]; \
										 X[B]=msu; }

#define REVERSEVEC(A,B,V) 				{Uint mi;\
										for(mi=0; mi<(B-A); mi++) {\
										SWAPVEC(A+mi,B-mi,V);\
										}}
#define EMPTYVEC(V)						(V->elements==NULL)
#define MAX(A,B)						((A) >= (B) ? (A) : (B))
#define MAX3(A, B, C)					MAX(MAX((A),(B)),(C))
#define MAX4(A, B, C, D)				MAX(MAX(MAX((A),(B)),(C)),(D))
#define MIN(A,B)						((A) <= (B) ? (A) : (B))
#define MAX_DOUBLE						1e200
#define ABS(X)							((X)>=0?(X):-(X))
#define OVERLAP(A,B,C,D)                (((A) >= (C) && (B) <= (D)) || \
                                        ((A) <= (D) && (B) >= (D))  || \
                                        ((B) >= (C) && (A) <= (C)))            
#define CLOSEDINTERVAL(X,A,B)           ((A) >= (X) && (X) <= (B))
#define OPENINTERVAL(X,A,B)             ((A) > (X) && (X) < (B))
#define ISELEM(X,A,B)                   ((A) <= (X) && (X) <= (B))

#ifndef FLT_EPSILON
	#define FLT_EPSILON     1.192092896e-06 
	#define DBL_EPSILON     2.2204460492503131e-016
#endif

#ifndef M_PI
    #define M_PI 3.14
#endif

#ifndef ALLOCMEMORY
	#include "memory.h"
#endif

#ifndef VECTORELEM
typedef Uint vectorelem;
#endif

typedef struct {
	vectorelem *elements;
	int length;
} vector_t;

void *initArray(void *, int, size_t);
void dumpMatrix_int(int *, int, int);
void dumpMatrix_Uint(Uint *, Uint, Uint);
void dumpMatrix_flt(float *, int, int);
void dumpMatrix_dbl(double *, Uint, Uint);
void dumpMatrix3D_int(int *, int, int, int);
double *transpose(void* space, double *, Uint, Uint);
void appendvector(void *, vector_t *, vectorelem);
void dumpVector(vector_t *);
void destructVector(void *, vector_t *);
Uint uarraymax(Uint *, Uint);
int arraymax(int *, int);
int gcd(int, int);
double power(double, int);
double uniroot(double, double, double (*f)(double, void*), double, void*);
double BLAST_Expm1(double x);
Uint fak(Uint);
double* coldel (void *, double *, Uint, Uint, Uint);
double* rowdel (void *, double *, Uint, Uint, Uint);
Uint minvecdist(void *space, vector_t *vec, Uint i, Uint j);
int* intrev(int *n, Uint len);

#endif
