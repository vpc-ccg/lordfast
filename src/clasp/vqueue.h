/**
 * vqueue.h
 * implementation of a simple queue for objects of defined size
 *
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Mon Oct 13 14:13:08 CEST 2008
 */

/*
 * SVN
 * Revision of last commit: $Rev: 89 $
 * Author: $Author: steve $
 * Date: $Date: 2008-11-24 14:53:55 +0100 (Mon, 24 Nov 2008) $
 * Id: $Id: vqueue.h 89 2008-11-24 13:53:55Z steve $
 * Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/vqueue.h $
 */

#ifndef VQUEUE_H
#define VQUEUE_H

#include <stdlib.h>
#include "basic-types.h"

typedef struct {
  void *queuespace;
  int enqueueindex;
  int dequeueindex;
  Uint numofelem;
  int allocelem;
  size_t sizeofelem;
} VQueue;

void bl_vqueueInit(VQueue *q, int allocelem, size_t sizeofelem);
void bl_vqueueDestruct(VQueue *q, void (*rmv)(void*));
BOOL bl_vqueueIsEmpty(VQueue *q);
void bl_vqueueResize(VQueue *q);
void bl_vqueueEnqueue(VQueue *q, void *elem);
void* bl_vqueueDequeue(VQueue *q, void (*rmv)(void*));
void* bl_vqueueFront(VQueue *q);
void* bl_vqueueFrontN(VQueue *q, int N);
Uint bl_vqueueSize(VQueue *q);

#endif /* VQUEUE_H */
