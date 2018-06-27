/*
 * Copyright (c) <2008 - 2020>, University of Washington, Simon Fraser University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, 
 * are permitted provided that the following conditions are met:
 *   
 * Redistributions of source code must retain the above copyright notice, this list
 * of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or other
 *   materials provided with the distribution.
 * - Neither the name of the <ORGANIZATION> nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Author: 
 *        Faraz Hach (fhach AT sfu DOT ca)
 *        Iman Sarrafi (isarrafi AT sfu DOT ca)
 *        Ehsan Haghshenas (ehaghshe AT sfu DOT ca)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include <pthread.h>
#include "Common.h"
#include "Reads.h"

pthread_t       *_r_threads;
pthread_mutex_t _r_readIdLock;
FILE            *_r_fp;
FILE            *_r_umfp;
gzFile          _r_gzfp;

Read            *_r_seq;
uint32_t        _r_seqCnt;
uint32_t        _r_seqPos;
uint32_t        _r_maxSeqCnt = 1000000;
char            *_r_alphIndex = NULL;
char            *_r_buf;
uint32_t        *_r_buf_size;
uint32_t        *_r_buf_pos;
uint8_t         _r_fastq;
uint8_t         _r_firstIteration = 1;
uint64_t        _r_maxReadMemUsage;
uint64_t        _r_readMemUsage = 0;
uint32_t        _r_pa=0, _r_pc=0, _r_pg=0, _r_pt=0;
// uint16_t        _r_samplingInterval = 500;
// uint16_t        _r_samplingCount = 250;

/**********************************************/
void (*readBuffer)();
/**********************************************/
void readBufferGZ()
{
    (*_r_buf_size) = gzread(_r_gzfp, _r_buf, 10000000);
    (*_r_buf_pos) = 0;
}
/**********************************************/
void readBufferTXT()
{
    (*_r_buf_size) = fread(_r_buf, 1, 10000000, _r_fp);
    (*_r_buf_pos) = 0;
}
/**********************************************/
int readSeq(char *seq)
{
    int i=0, l=0;
    char cur;

    while (1)
    {
        if (*_r_buf_pos ==*_r_buf_size)
        {
            readBuffer();
            if (*_r_buf_size == 0)
                return 0;
        }

        cur = _r_buf[*_r_buf_pos];      
        (*_r_buf_pos)++;

        if ( cur == '\n')
        {
            if (l>0) i=l;
            seq[i]='\0';
            return i;
        }
        
        if (l==0 && cur == ' ')
        {
            l = i;
        }

        seq[i++]=cur;
    }
}
/**********************************************/
int getNextRead() {
    int x;
    pthread_mutex_lock(&_r_readIdLock);
        x = _r_seqPos;
        _r_seqPos++;
    pthread_mutex_unlock(&_r_readIdLock);
    return x;
}
/**********************************************/
void hash(char *str, int count, int *hvs)
{
    uint64_t windowMask = 0xffffffffffffffff >> (sizeof(uint64_t)*8 - MIN_ANCHOR_LEN*2);
    int stack=0, i, val, hv;
//  int pos=0;

    for (i=0; i < count+MIN_ANCHOR_LEN-1; i++)
    {
        val = _r_alphIndex[str[i]];
        stack++;

        if (val == 4) 
        {
            hv = 0;
            while (stack>0) {
                *(hvs++)=-1;
                stack--;
//              fprintf(stdout, "%3d => %c %3d %8llx | %d %6llx\n", i, str[i], stack, hv, pos++, *(hvs-1));
            }
        }
        else
        {
            hv = ((hv << 2) | val)&windowMask;
            if (stack == MIN_ANCHOR_LEN) {
                if (hv == _r_pa || hv == _r_pc || hv == _r_pg || hv == _r_pt)
                    *(hvs++) = -1;
                else
                    *(hvs++) = hv;
                stack --;
//              fprintf(stdout, "%3d => %c %3d %8llx | %d %6llx\n", i, str[i], stack+1, hv, pos++, *(hvs-1));
            }
//          else 
//              fprintf(stdout, "%3d => %c %3d %8llx\n", i, str[i], stack, hv);
        }
    }
}
/**********************************************/
void qGramCount(char *str, int len, uint8_t* count)
{
    int i;
    char *head=str;
    char *tail;
    uint8_t* cur = count;
    uint32_t* copy = (uint32_t*)(count+4);
    
    cur[0]=cur[1]=cur[2]=cur[3]=0;

    for (i=0; i<QGRAM_WIN_SIZE; i++)
    {
        cur[_r_alphIndex[str[i]]]++;
    }


    i=QGRAM_WIN_SIZE;

    head = str;
    tail = str+QGRAM_WIN_SIZE;
    while (i<len) 
    {
        *copy = *(copy-1);
        cur = (uint8_t*)copy++;
        cur[_r_alphIndex[*head++]]--;
        cur[_r_alphIndex[*tail++]]++;
        i++;
    }

}
/**********************************************/
void* preProcessReads(void *idp)
{
    int id = *(int*)idp;
    Read *r;
    int t, i, startPos, pos = 0, stack = 1, /*intervals,*/ k;
    char val;
    int32_t hv,sv;  

    while ((t=getNextRead())<_r_seqCnt){

        // Prcessing read t
        r = _r_seq + t;
        // fprintf(stdout,">\t%d\t%d\t%d\n", id, t,*r->length);

        // intervals = (*r->length / SAMPLING_INTERVAL ) + ((*r->length % SAMPLING_INTERVAL >= SAMPLING_COUNT+MIN_ANCHOR_LEN-1)?1:0);

        // // forward hash
        // startPos = 0;
        // for (k=0; k<intervals;k++)
        // {
        //  hash(r->seq+startPos, SAMPLING_COUNT, r->hv+k*SAMPLING_COUNT);
        //  startPos += SAMPLING_INTERVAL;
        // }
        
        // reverseComplete(r->seq, r->rseq, *r->length);

        // //reverse hash
        // startPos = 0;
        // for (k=0; k<intervals;k++)
        // {
        //  hash(r->rseq+startPos, SAMPLING_COUNT, r->rhv+k*SAMPLING_COUNT);
        //  startPos += SAMPLING_INTERVAL;
        // }

        //1-gram count forward
        // qGramCount (r->seq, *r->length, r->alphCnt);

        //1-gram count reverse
        // qGramCount (r->rseq, *r->length, r->ralphCnt);
    }
}
/**********************************************/
void preProcessReadsMT()
{
    _r_threads = (pthread_t *) getMem(sizeof(pthread_t)*THREAD_COUNT);
    int i=0; 
    for (i=0; i<THREAD_COUNT; i++)
        pthread_create(_r_threads+i, NULL, preProcessReads, THREAD_ID+i);
    
    for (i=0; i<THREAD_COUNT; i++)
        pthread_join(_r_threads[i], NULL);
    freeMem(_r_threads, sizeof(pthread_t)*THREAD_COUNT);
}

/**********************************************/
int initRead(char *fileName, int maxMem)
{
    char dummy[SEQ_MAX_LENGTH];
    char ch;
    int i, maxCnt=0;

    for (i=0; i<MIN_ANCHOR_LEN; i++)
    {
        _r_pc = (_r_pc << 2) | 1;
        _r_pg = (_r_pg << 2) | 2;
        _r_pt = (_r_pt << 2) | 3;
    }


    _r_maxReadMemUsage = maxMem;

    _r_buf = (char*) getMem(10000000);      // 10MB
    _r_buf_pos = (uint32_t*) getMem(sizeof(int));
    _r_buf_size = (uint32_t*) getMem(sizeof(int));
    *_r_buf_size = *_r_buf_pos = 0; 


    if (!seqCompressed)
    {
        _r_fp = fileOpen( fileName, "r");

        if (_r_fp == NULL)
            return 0;

        readBuffer = &readBufferTXT;
    }
    else
    {

        _r_gzfp = fileOpenGZ (fileName, "r");

        if (_r_gzfp == NULL)
        {
            return 0;
        }

        readBuffer = &readBufferGZ;
    }

    readBuffer();

    if (_r_buf[0] == '>')
        _r_fastq = 0;
    else
        _r_fastq = 1;

    _r_seq = (Read*) getMem(sizeof(Read)*_r_maxSeqCnt);

    if (!nohitDisabled)
    {
        _r_umfp = fopen(outputUnmap, "w");
    }

    _r_alphIndex = (char*) getMem(128);     // used in readChunk()
    _r_alphIndex['A'] = 0;
    _r_alphIndex['C'] = 1;
    _r_alphIndex['G'] = 2;
    _r_alphIndex['T'] = 3;
    _r_alphIndex['N'] = 4;

    return 1;
}

/**********************************************/
int readChunk(Read **seqList, unsigned int *seqListSize)
{
    double startTime=getTime();

    char seq[SEQ_MAX_LENGTH];
    char name[SEQ_MAX_LENGTH];
    char qual[SEQ_MAX_LENGTH];

    char dummy[SEQ_MAX_LENGTH];
    int  size;

    int maxCnt = 0;
    _r_seqCnt = 0;
    _r_readMemUsage = 0;
    
    int i, readLen, qualLen, nameLen, qgramLen, hvLen;

/*  typedef struct
    {
        uint16_t *hits;
        uint16_t *length;
        char *seq;
        char *rseq;
        char *qual;
        char *name;
        uint8_t *alphCnt;
        uint8_t *ralphCnt;
    } Read;
*/
    qualLen = 1;    // in case it is fasta

    while( (nameLen = readSeq(name)) )
    {
        readLen = readSeq(seq);

        if (_r_fastq)
            qualLen = readLen;
        
        qgramLen = readLen; // TODO: readLen - QGRAM_WIN_SIZE
        // hvLen = ((readLen/SAMPLING_INTERVAL)+( (readLen % SAMPLING_INTERVAL >= SAMPLING_COUNT+MIN_ANCHOR_LEN-1) ?(1) :(0)) )*SAMPLING_COUNT;

        // size = 2*sizeof(uint16_t) + 2*(readLen + 1) + (qualLen + 1) + (nameLen+1) + 2*4*qgramLen + 2*sizeof(int32_t)*hvLen;  
        size = sizeof(uint32_t) + (readLen + 1) + (qualLen + 1) + (nameLen + 1) + sizeof(uint8_t);
        // _r_seq[_r_seqCnt].hits   = (uint16_t*) getMem(size);
        _r_readMemUsage += size;
        _r_seq[_r_seqCnt].length = (uint32_t*)getMem(size);
        // _r_seq[_r_seqCnt].length = (uint16_t *)(_r_seq[_r_seqCnt].hits + 1);
        _r_seq[_r_seqCnt].seq  = (char*)(_r_seq[_r_seqCnt].length + 1);
        // _r_seq[_r_seqCnt].rseq = (char*)(_r_seq[_r_seqCnt].seq + readLen + 1);
        _r_seq[_r_seqCnt].qual = (char*)(_r_seq[_r_seqCnt].seq + readLen + 1);
        _r_seq[_r_seqCnt].name = (char*)(_r_seq[_r_seqCnt].qual + qualLen + 1);
        _r_seq[_r_seqCnt].isFq = (uint8_t*)(_r_seq[_r_seqCnt].name + nameLen + 1);
        // _r_seq[_r_seqCnt].alphCnt = (uint8_t *)(_r_seq[_r_seqCnt].name + nameLen + 1);
        // _r_seq[_r_seqCnt].ralphCnt = (uint8_t *)(_r_seq[_r_seqCnt].alphCnt + qgramLen*4);
        // _r_seq[_r_seqCnt].hv = (int32_t *)(_r_seq[_r_seqCnt].alphCnt + qgramLen);
        // _r_seq[_r_seqCnt].hv = (int32_t *)(_r_seq[_r_seqCnt].ralphCnt + qgramLen*4);
        // _r_seq[_r_seqCnt].rhv = (int32_t *)(_r_seq[_r_seqCnt].hv + hvLen);
        // _r_seq[_r_seqCnt].hits[0] = 0;


        for (i=1; i<nameLen+1; i++)
            _r_seq[_r_seqCnt].name[i-1] = name[i];

        strcpy(_r_seq[_r_seqCnt].seq, seq);
        *_r_seq[_r_seqCnt].length = readLen;


        if ( _r_fastq )
        {
            readSeq(dummy);     // comment line
            qualLen = readSeq(qual);
            strcpy(_r_seq[_r_seqCnt].qual, qual);
            *_r_seq[_r_seqCnt].isFq = 1;
        }
        else
        {
            // _r_seq[_r_seqCnt].qual = "*";
            strcpy(_r_seq[_r_seqCnt].qual, "*");
            *_r_seq[_r_seqCnt].isFq = 0;
        }
        _r_seqCnt++;


        if (_r_readMemUsage >= _r_maxReadMemUsage)
            break;
    }

    *seqList = _r_seq;
    *seqListSize = _r_seqCnt;

    _r_seqPos = 0;
    if (_r_seqCnt > 0)
    {
        // preProcessReadsMT();
        fprintf(stderr, "| *Reading Input* | %15.2f | XXXXXXXXXXXXXXX | %15.2f | XXXXXXXXXXXXXXX %15d |\n", (getTime()-startTime), getMemUsage(), _r_seqCnt );
        _r_firstIteration = 0;
    }
    else if (_r_firstIteration)
    {
        fprintf(stderr, "[readChunk] ERROR: No reads for mapping\n");
    }

    if (_r_seqCnt < _r_maxSeqCnt)       // reached end of file
        return 0;
    else
        return 1;
}
/**********************************************/
// void outputUnmapped()
// {
//  if (nohitDisabled)
//      return;

//  int i = 0;
//  for (i = 0; i < _r_seqCnt; i++)
//  {
//      if (_r_seq[i].hits[0] == 0 && _r_fastq)
//      {
//          fprintf(_r_umfp,"@%s\n%s\n+\n%s\n", _r_seq[i].name, _r_seq[i].seq, _r_seq[i].qual);
//      }
//      else if (_r_seq[i].hits[0] == 0)
//      {
//          fprintf(_r_umfp,">%s\n%s\n", _r_seq[i].name, _r_seq[i].seq);
//      }
//  }
// }
/**********************************************/
void releaseChunk()
{
    // outputUnmapped();

    int i, j;
    for (i = 0; i < _r_seqCnt; i++)
    {
        // freeMem(_r_seq[i].hits, 0);
        freeMem(_r_seq[i].length, 0);
    }
    memUsage -= _r_readMemUsage;
    _r_readMemUsage = 0;
}
/**********************************************/
void finalizeReads()
{
    if (!seqCompressed)
        fclose(_r_fp);
    else
        gzclose(_r_gzfp);
    freeMem(_r_seq, sizeof(Read)*_r_maxSeqCnt);
    freeMem(_r_alphIndex, 128);

    freeMem(_r_buf, 10000000);
    freeMem(_r_buf_pos, sizeof(int));
    freeMem(_r_buf_size, sizeof(int));

    if (!nohitDisabled)
    {
        fclose(_r_umfp);
    }
}

/**********************************************/
void getSamplingLocsInfo(int **samplingLocs, int **samplingLocsSeg, int **samplingLocsOffset, int **samplingLocsLen, int **samplingLocsLenFull, int *samplingLocsSize) {}
void getReadIndex(ReadIndexTable ***rIndex, int **rIndexSize) {}
