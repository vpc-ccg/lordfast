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
    int32_t hv;
    // CheckSumType checksum;
    int32_t seqInfo;

} Pair;

typedef struct
{
    // uint16_t *hits;
    uint32_t *length;
    char *seq;
    char *qual;
    // char *rseq;
    char *name;
    uint8_t *isFq;
    // uint8_t *alphCnt;
    // uint8_t *ralphCnt;
    // int32_t *hv;
    // int32_t *rhv;
    // CompressedSeq *cseq;
    // CompressedSeq *crseq;
} Read;

int readChunk(Read **seqList, unsigned int *seqListSize);
void finalizeReads();
// void getSamplingLocsInfo(int **samplingLocs, int **samplingLocsSeg, int **samplingLocsOffset, int **samplingLocsLen, int **samplingLocsLenFull, int *samplingLocsSize);// {}
//int initRead(char *seqFile1, char *seqFile2);
int initRead(char *seqFile1, int maxMem);
// void getReadIndex(ReadIndexTable ***rIndex, int **rIndexSize);// {}
void releaseChunk();

#endif // __READ__
