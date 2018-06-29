/**
 * vstack.h
 * implementation of a simple stack for objects of defined size
 *
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Fri Oct 10 11:37:36 CEST 2008
 */

/*
 * SVN
 * Revision of last commit: $Rev: 89 $
 * Author: $Author: steve $
 * Date: $Date: 2008-11-24 14:53:55 +0100 (Mon, 24 Nov 2008) $
 * Id: $Id: vstack.h 89 2008-11-24 13:53:55Z steve $
 * Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/vstack.h $
 */

#ifndef VSTACK_H
#define VSTACK_H

#include <stdio.h>
#include <stdlib.h>
#include "basic-types.h"

#define VSTACKINC 10000

#ifndef BASEINC
#define BASEINC VSTACKINC
#endif

typedef struct{
  void* stackspace;
  int allocelem;
  int top;
  size_t sizeofelem;
} VStack;

void bl_vstackInit(VStack *s, int allocelem, size_t sizeofelem);
void bl_vstackDestruct(VStack *s, void (*rmv)(void*));
BOOL bl_vstackIsEmpty(VStack *s);
void bl_vstackPush(VStack *s, void *elem);
void* bl_vstackTop(VStack *s);
void* bl_vstackTopN(VStack *s, int n);
void* bl_vstackPop(VStack *s, void (*rmv)(void*));
Uint bl_vstackSize(VStack *s);

#endif /* VSTACK_H */
