#ifndef MEMORY_H
#define MEMORY_H
/*
#include "memman.h"
#include "memmac.h"
*/
#define ALLOCMEMORY(X,PTR,TYPE,SIZE) realloc(PTR,sizeof(TYPE)*(SIZE))
#define FREEMEMORY(X,PTR) free(PTR)

#endif

