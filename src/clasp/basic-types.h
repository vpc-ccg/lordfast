#ifndef BASIC_TYPES_H
#define BASIC_TYPES_H


#ifdef __CYGWIN__ 

#define CRLF '\r'

#else
#define CRLF ' '
#endif

#define MAXBUFFERSIZE 10000
#define BASEINC 10000
typedef unsigned char Uchar;
typedef unsigned int Uint;
typedef signed int Sint;
typedef unsigned char BOOL;

#define True 1
#define False 0

typedef struct {
  int  a, 
       b;
} PairSint; 


typedef struct {
  Uint  a, 
       b;
} PairUint; 


typedef struct {
  Uint a,
      b,
      c;
} TripleUint;


typedef struct {
  int a,
      b,
      c;
} TripleSint;

typedef struct {
  int a,
      b,
      c,
      d;
} QuadSint;



#endif

