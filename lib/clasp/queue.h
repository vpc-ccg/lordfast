/*
	  This file is part of gdub.
	  (C) 2006 Steve Hoffmann 
 
	  gdub is free software; you can redistribute it and/or modify
	  it under the terms of the GNU General Public License as published
	  by the Free Software Foundation; either version 2, or (at your
	  option) any later version.
 
	  gdub is distributed in the hope that it will be useful, but
	  WITHOUT ANY WARRANTY; without even the implied warranty of
	  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	  General Public License for more details.
 
	  You should have received a copy of the GNU General Public License
	  along with gdub; see the file COPYING.  If not, write to the
	  Free Software Foundation, Inc., 59 Temple Place - Suite 330,
	  Boston, MA 02111-1307, USA.	
 
 */

/**
 * queue.c
 * implementation of a simple queue for int
 *
 * @author Steve Hoffmann
 * @email steve@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Mon Oct 13 14:13:08 CEST 2008
 */

/*
 * SVN
 * Revision of last commit: $Rev: 72 $
 * Author: $Author: steve $
 * Date: $Date: 2008-10-28 18:14:42 +0100 (Tue, 28 Oct 2008) $
 * Id: $Id: queue.h 72 2008-10-28 17:14:42Z steve $
 * Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/queue.h $
 */

#ifndef QUEUE_H
#define QUEUE_H

#include "basic-types.h"

typedef int Queueelem;

typedef struct
{
  Queueelem *queuespace;
  int 	enqueueindex,
    dequeueindex,
    allocelem,
    numofelem;
} Queue;

void bl_queueInit(Queue *q, int allocelem);
void bl_queueDestruct(Queue *q);
BOOL bl_queueIsEmpty(Queue *q);
void bl_queueResize(Queue *q);
void bl_queueEnqueue(Queue *q, Queueelem elem);
Queueelem bl_queueDequeue(Queue *q);
void bl_queueShow(Queue *q);
Uint bl_queueSize(Queue *q);

#endif 
