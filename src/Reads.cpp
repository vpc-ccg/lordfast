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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include <pthread.h>
#include "Common.h"
#include "Reads.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

pthread_t       *_r_threads;
gzFile          _r_gzfp;
kseq_t          *_r_ks;

Read            *_r_seq;
uint32_t        _r_seqCnt;
uint32_t        _r_maxSeqCnt;
uint8_t         _r_firstIteration = 1;
uint64_t        _r_maxReadMemUsage;
uint64_t        _r_readMemUsage = 0;

/**********************************************/
int initRead(char *fileName, int maxMem)
{
    _r_maxReadMemUsage = maxMem;
    _r_maxSeqCnt = maxMem / 500; // assume average read length of 500 to be safe
    _r_seq = (Read*) getMem(_r_maxSeqCnt * sizeof(Read));

    _r_gzfp = fileOpenGZ (fileName, "r");
    if (_r_gzfp == NULL)
    {
        fprintf(stderr, "[ERROR] (initRead) could not open file: %s\n", fileName);
        return 0;
    }
    _r_ks = kseq_init(_r_gzfp);

    return 1;
}

/**********************************************/
int readChunk(Read **seqList, unsigned int *seqListSize)
{
    fprintf(stderr, "Reading input... ");
    double ct = getCpuTime();
    double rt = getRealTime();
    int  size;

    _r_seqCnt = 0;
    _r_readMemUsage = 0;
    
    int i, readLen, qualLen, nameLen;

    while (kseq_read(_r_ks) >= 0)
    {
        nameLen = _r_ks->name.l;
        readLen = _r_ks->seq.l;
        qualLen = _r_ks->qual.l;
        if(qualLen == 0)
            qualLen = 1;

        size = sizeof(uint32_t) + (readLen + 1) + (qualLen + 1) + (nameLen + 1) + sizeof(uint8_t);
        _r_readMemUsage += size;
        _r_seq[_r_seqCnt].length = (uint32_t*)getMem(size);
        _r_seq[_r_seqCnt].seq  = (char*)(_r_seq[_r_seqCnt].length + 1);
        _r_seq[_r_seqCnt].qual = (char*)(_r_seq[_r_seqCnt].seq + readLen + 1);
        _r_seq[_r_seqCnt].name = (char*)(_r_seq[_r_seqCnt].qual + qualLen + 1);
        _r_seq[_r_seqCnt].isFq = (uint8_t*)(_r_seq[_r_seqCnt].name + nameLen + 1);

        strcpy(_r_seq[_r_seqCnt].name, _r_ks->name.s);
        strcpy(_r_seq[_r_seqCnt].seq, _r_ks->seq.s);
        *_r_seq[_r_seqCnt].length = readLen;

        if ( _r_ks->qual.l > 0)
        {
            strcpy(_r_seq[_r_seqCnt].qual, _r_ks->qual.s);
            *_r_seq[_r_seqCnt].isFq = 1;
        }
        else
        {
            strcpy(_r_seq[_r_seqCnt].qual, "*");
            *_r_seq[_r_seqCnt].isFq = 0;
        }
        _r_seqCnt++;

        if (_r_readMemUsage >= _r_maxReadMemUsage)
            break;
    }

    *seqList = _r_seq;
    *seqListSize = _r_seqCnt;

    if (_r_seqCnt > 0)
    {
        fprintf(stderr, "loaded %u reads in %.2f seconds (%.2f CPU seconds)\n", _r_seqCnt, getRealTime()-rt, getCpuTime()-ct);
        _r_firstIteration = 0;
    }
    else
    {
        fprintf(stderr, "no more reads\n");
    }

    if (_r_seqCnt > 0 && _r_firstIteration)
    {
        fprintf(stderr, "[WARNING] (readChunk) there is no read to map\n");
    }
    
    return _r_seqCnt;
}
/**********************************************/
void releaseChunk()
{
    int i, j;
    for (i = 0; i < _r_seqCnt; i++)
    {
        freeMem(_r_seq[i].length, 0);
    }
    memUsage -= _r_readMemUsage;
    _r_readMemUsage = 0;
}
/**********************************************/
void finalizeReads()
{
    kseq_destroy(_r_ks);
    gzclose(_r_gzfp);
    freeMem(_r_seq, _r_maxSeqCnt * sizeof(Read));
}
