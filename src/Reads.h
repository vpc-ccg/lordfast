/*
 * lordFAST: sensitive and Fast Alignment Search Tool for LOng noisy Read sequencing Data
 * Copyright (C) 2018 Simon Fraser University

 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * Author: Ehsan Haghshenas (ehaghshe AT sfu DOT ca)
 */

#ifndef __READ__
#define __READ__

#include "Common.h"

typedef struct
{
    uint32_t *length;
    char *seq;
    char *qual;
    char *name;
    uint8_t *isFq;
} Read;

int readChunk(Read **seqList, unsigned int *seqListSize);
void finalizeReads();
int initRead(char *seqFile1, int maxMem);
void releaseChunk();

#endif // __READ__
