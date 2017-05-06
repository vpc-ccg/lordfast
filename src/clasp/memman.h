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
* @file memman.h
 * @author Steve Hoffmann
 * @brief declarations of functions for memory management
 */

/*
 * $Log$
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: memman.h 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/memman.h $
 */

#ifndef MEMMAN_H 
#define MEMMAN_H

typedef struct
{
	void *spaceptr;
	int  sizeofcell, numberofcells;
	char* fileallocated;
	int lineallocated;

} Spaceblock;

typedef struct 
{
	int numberofblocks,
		lastalloced,
		lastfreed;
	Spaceblock *blocks;

} Spacetable;

void initmemoryblocks(Spacetable *st, int numberofblocks);
void *allocmemory(char* file, int line, Spacetable *st, void* ptr, int size, int number);
void freememory(char* file, int line, Spacetable *st, void *ptr);
void activeblocks(Spacetable *st);
void checkspaceleak(Spacetable *st);

#endif
