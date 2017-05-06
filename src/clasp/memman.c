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
 * @file memman.c
 * @author Steve Hoffmann
 * @brief functions for memory management
*/

/* 
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: memman.c 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/memman.c $
 *  
 */

#include "memman.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

void initmemoryblocks(Spacetable *st, int numberofblocks) {
	
	int i;
	
	/*alloc spacetable and blocks*/
	/*st = (Spacetable*) malloc(sizeof(Spacetable))*/
	st->numberofblocks = numberofblocks;
	st->lastalloced = 0;
	st->lastfreed=0;
	st->blocks = (Spaceblock*) malloc(st->numberofblocks*sizeof(Spaceblock));

	/*init blocks*/
	for (i=0; i < st->numberofblocks; i++) {
		st->blocks[i].spaceptr = NULL;
		st->blocks[i].fileallocated = NULL;
		st->blocks[i].lineallocated = 0;
		st->blocks[i].sizeofcell = 0;
		st->blocks[i].numberofcells = 0;
	}

}

void *allocmemory(char *file, int line, Spacetable *st, void *ptr, int size, int number) {
	
	int i;
	Spaceblock* rescueptr;
	
	/*alloc new block*/
	if (ptr==NULL) {
	
		if (st->lastfreed != 0) 
			st->lastalloced = st->lastfreed;
		
		/*find free block*/
		for(i=st->lastalloced; i < st->numberofblocks; i++) {
			if(st->blocks[i].numberofcells == 0) break;
		}

		if(i == st->numberofblocks) {

			st->numberofblocks++;
		
			/*alloc in spacetable*/
			rescueptr = (Spaceblock*) realloc(st->blocks, st->numberofblocks*sizeof(Spaceblock));
			assert(rescueptr != NULL);
			st->blocks = rescueptr;
		
			/*alloc in spaceblock*/
			st->blocks[st->numberofblocks-1].sizeofcell = size;
			st->blocks[st->numberofblocks-1].numberofcells = number;
			st->blocks[st->numberofblocks-1].fileallocated = file;
			st->blocks[st->numberofblocks-1].lineallocated = line;
			st->blocks[st->numberofblocks-1].spaceptr = (void*) malloc(number*size);	
			st->lastalloced = st->numberofblocks-1;
			st->lastfreed = 0;
			
			return st->blocks[st->numberofblocks-1].spaceptr;
		
		} else {
			
			st->blocks[i].sizeofcell = size;
			st->blocks[i].numberofcells = number;
			st->blocks[i].fileallocated = file;
			st->blocks[i].lineallocated = line;
			st->blocks[i].spaceptr = (void*) malloc(number*size);	
			st->lastalloced = i;
			st->lastfreed = 0;
			
			return st->blocks[i].spaceptr;
		}
	}


	/*resize block*/
	if(ptr != NULL) {
		
		/*get blockno*/
		for (i=0; i < st->numberofblocks; i++) {
			if (st->blocks[i].spaceptr == ptr) break;
		}
		
		assert(i < st->numberofblocks);
		st->blocks[i].sizeofcell = size;
		st->blocks[i].numberofcells = number;
		st->blocks[i].lineallocated = line;
		st->blocks[i].fileallocated = file;
		
		rescueptr = (void*) realloc(st->blocks[i].spaceptr, size*number);
		assert(rescueptr != NULL);
		
		st->blocks[i].spaceptr = rescueptr;
	
		return st->blocks[i].spaceptr;
		
	}
	
	/*a stub for the compiler*/
	return NULL;
}

void freememory(char* file, int line, Spacetable *st, void *ptr) {
		int i;
	  	
		for (i=0; i < st->numberofblocks; i++) {
			if (st->blocks[i].spaceptr == ptr && st->blocks[i].numberofcells > 0) break;
		}
	
		if (i >= st->numberofblocks) {
			printf("Attempt to free unallocated spaceblock in line %d, %s \n",line,file);
			exit(-1);
		}
		
		st->blocks[i].numberofcells = 0;
		st->blocks[i].sizeofcell = 0;
		free(st->blocks[i].spaceptr);
		st->blocks[i].spaceptr = NULL;
		
		
		return;
}

void activeblocks(Spacetable *st) {
	int i;

	for(i=0; i < st->numberofblocks; i++) {
		if (st->blocks[i].numberofcells > 0) {
			printf("# active block %d: allocated with with %d cells in file \"%s\", line %d \n", i, st->blocks[i].numberofcells, st->blocks[i].fileallocated, st->blocks[i].lineallocated);
		}
	}

}

void checkspaceleak(Spacetable *st) {
	int i;

	for(i=0; i < st->numberofblocks; i++) {
		if (st->blocks[i].numberofcells > 0){
			printf("space leak: memory for block.%d not freed \n %d cells of size %d \n", i, st->blocks[i].numberofcells, st->blocks[i].sizeofcell);
			printf("allocated in file \"%s\", line %d \n", st->blocks[i].fileallocated, st->blocks[i].lineallocated);
			break;
		}
	}

	return;
}

