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
#include <sys/time.h>
#include <zlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include "Common.h"


unsigned short          SEQ_LENGTH = 0;
uint16_t                QGRAM_WIN_SIZE = 100;
unsigned short          QUAL_LENGTH = 0;
unsigned short          CMP_SEQ_LENGTH = 0;
long long               memUsage = 0;
// char                    *alphabet = "ACGTN";
char                    nVal[128];
char                    nRev[128];
char                    nHVal[128];
char                    _common_nRevVal[5];
pthread_mutex_t         _common_lock;


/**********************************************/
FILE *fileOpen(char *fileName, const char *mode)
{
    FILE *fp;
    fp = fopen (fileName, mode);
    if (fp == NULL)
    {
        fprintf(stdout, "Error: Cannot Open the file %s\n", fileName);
        fflush(stdout);
        exit(0);
    }
    return fp;
}
/**********************************************/
gzFile fileOpenGZ(char *fileName, const char *mode)
{
    gzFile gzfp;
    gzfp = gzopen (fileName, mode);
    if (gzfp == Z_NULL)
    {
        fprintf(stdout, "Error: Cannot Open the file %s\n", fileName);
        fflush(stdout);
        exit(0);
    }
    return gzfp;
}
/**********************************************/
double getTime(void)
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return t.tv_sec+t.tv_usec/1000000.0;
}

/**********************************************/
// inline char reverseCompleteChar(char c)
char reverseCompleteChar(char c)
{
    return nRev[c];
}
/**********************************************/
// inline void reverseComplete (char *seq, char *rcSeq , int length)        // TODO: efficiency check
void reverseComplete (char *seq, char *rcSeq , int length)      // TODO: efficiency check
{
    int i;
    seq+=length-1;
    for (i=0; i<length; i++)
    {
        rcSeq[i]=nRev[*(seq--)];
    }
    rcSeq[length]='\0';
}
/**********************************************/
void * getMem(size_t size)          // TODO: if malloc is unsuccessfull, return an error message
{
    
    pthread_mutex_lock(&_common_lock);
    memUsage+=size;
    pthread_mutex_unlock(&_common_lock);
    return malloc(size);
}
/**********************************************/
void freeMem(void *ptr, size_t size)
{
    pthread_mutex_lock(&_common_lock);
    memUsage-=size;
    pthread_mutex_unlock(&_common_lock);
    free(ptr);
}
/**********************************************/
double getMemUsage()
{
    return memUsage/1048576.0;
}
/**********************************************/
// inline void reverse (char *seq, char *rcSeq , int length)
void reverse (char *seq, char *rcSeq , int length)
{
    int i;
    seq += length-1;
    for (i=0; i<length; i++)
    {
        (*rcSeq++)=*(seq--);
    }
    *rcSeq='\0';
}
/**********************************************/
void stripPath(char *full, char *path, char *fileName)
{
    int i;
    int pos = -1;

    for (i=strlen(full)-1; i>=0; i--)
    {
        if (full[i]=='/')
        {
            pos = i;
            break;
        }

    }

    if (pos != -1)
    {
        sprintf(fileName, "%s%c", (full+pos+1), '\0');
        full[pos+1]='\0';
        sprintf(path,"%s%c", full, '\0');
    }
    else
    {
        sprintf(fileName, "%s%c", full, '\0');
        sprintf(path,"%c", '\0');
    }
}
/**********************************************/
inline int calculateCompressedLen(int normalLen)
{
    return (normalLen / 21) + ((normalLen%21)?1:0);
}
/**********************************************/
void compressSequence(char *seq, int seqLen, CompressedSeq *cseq)
{
    CompressedSeq val = 0;
    int i = 0, pos = 0;
    
    *cseq = 0;
    while (pos < seqLen)
    {
        *cseq = ((*cseq) << 3) | nVal[seq[pos++]];

        if (++i == 21)
        {
            i = 0;
            cseq++;
            if (pos < seqLen)   // not to write the adjacent memory in case seqLen % 21 == 0
                *cseq = 0;
        }
    }
    if (i > 0)
    {
        *cseq <<= (3*(21-i));
    }
}
/**********************************************/
void decompressSequence(CompressedSeq *cseq, int seqLen, char *seq)
{
    int i;
    int shifts = 20*3;

    for (i = 0; i < seqLen; i++)
    {
        *(seq++) = _common_nRevVal[((*cseq) >> shifts) & 7];
        shifts -= 3;
        if (shifts < 0)
        {
                cseq++;
                shifts = 20*3;
        }
    }
    *seq = '\0';
}
/**********************************************/
int hashVal(char *seq)
{
    int i=0;
    int val=0, numericVal=0;

    while(i<WINDOW_SIZE)
    {
        if (nHVal[seq[i]] == -1)
            return -1; 
        val = (val << 2) | nHVal[seq[i++]]; 
    }
    return val;
}
/**********************************************/
void initCommon()
{
    memset(nVal, 4, 128);
    nVal['A']=0;
    nVal['C']=1;
    nVal['G']=2;
    nVal['T']=3;

    memset(nRev, 'N', 128);
    nRev['A']='T';
    nRev['C']='G';
    nRev['T']='A';
    nRev['G']='C';

    memset(nHVal, -1, 128);
    nHVal['A']=0;
    nHVal['C']=1;
    nHVal['G']=2;
    nHVal['T']=3;

    _common_nRevVal[0] = 'A';
    _common_nRevVal[1] = 'C';
    _common_nRevVal[2] = 'G';
    _common_nRevVal[3] = 'T';
    _common_nRevVal[4] = 'N';
}
/**********************************************/
void reverseInPlace(char *dest, char *src, int len)
{
    int i;
    src+=len-1;
    for (i=0; i<len; i++)
    {
        *(dest++)=*(src--);
    }
    *dest='\0';
    // dest[len]='\0';
}