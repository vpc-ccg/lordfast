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
#include <sys/time.h>
#include <zlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include "Common.h"

char tableRev[128] = {
    'N','N','N','N',   'N','N','N','N',   'N','N','N','N',   'N','N','N','N', 
    'N','N','N','N',   'N','N','N','N',   'N','N','N','N',   'N','N','N','N', 
    'N','N','N','N',   'N','N','N','N',   'N','N','N','N',   'N','N','N','N', 
    'N','N','N','N',   'N','N','N','N',   'N','N','N','N',   'N','N','N','N', 
    'N','T','N','G',   'N','N','N','C',   'N','N','N','N',   'N','N','N','N', 
    'N','N','N','N',   'A','N','N','N',   'N','N','N','N',   'N','N','N','N', 
    'N','t','N','g',   'N','N','N','c',   'N','N','N','N',   'N','N','N','N', 
    'N','N','N','N',   'a','N','N','N',   'N','N','N','N',   'N','N','N','N'
};

// char                    nVal[128];
// char                    nRev[128];
// char                    nHVal[128];
// char                    _common_nRevVal[5];
pthread_mutex_t         _common_lock;


/**********************************************/
// FILE *fileOpen(char *fileName, const char *mode)
// {
//     FILE *fp;
//     fp = fopen (fileName, mode);
//     if (fp == NULL)
//     {
//         fprintf(stdout, "Error: Cannot Open the file %s\n", fileName);
//         fflush(stdout);
//         exit(0);
//     }
//     return fp;
// }
/**********************************************/
gzFile fileOpenGZ(char *fileName, const char *mode)
{
    gzFile gzfp = gzopen (fileName, mode);
    if (gzfp == Z_NULL)
    {
        fprintf(stderr, "[Error] Cannot Open the file %s\n", fileName);
        exit(EXIT_FAILURE);
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
// char reverseCompleteChar(char c)
// {
//     return tableRev[c];
// }
/**********************************************/
// inline void reverseComplete (char *seq, char *rcSeq , int length)        // TODO: efficiency check
void reverseComplement(char *seq, char *rcSeq , int length)      // TODO: efficiency check
{
    int i;
    seq+=length-1;
    for (i=0; i<length; i++)
    {
        rcSeq[i]=tableRev[*(seq--)];
    }
    rcSeq[length]='\0';
}
/**********************************************/
void* getMem(size_t size)          // TODO: if malloc is unsuccessfull, return an error message
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
// void stripPath(char *full, char *path, char *fileName)
// {
//     int i;
//     int pos = -1;

//     for (i=strlen(full)-1; i>=0; i--)
//     {
//         if (full[i]=='/')
//         {
//             pos = i;
//             break;
//         }

//     }

//     if (pos != -1)
//     {
//         sprintf(fileName, "%s%c", (full+pos+1), '\0');
//         full[pos+1]='\0';
//         sprintf(path,"%s%c", full, '\0');
//     }
//     else
//     {
//         sprintf(fileName, "%s%c", full, '\0');
//         sprintf(path,"%c", '\0');
//     }
// }
/**********************************************/
// inline int calculateCompressedLen(int normalLen)
// {
//     return (normalLen / 21) + ((normalLen%21)?1:0);
// }
/**********************************************/
// void compressSequence(char *seq, int seqLen, CompressedSeq *cseq)
// {
//     CompressedSeq val = 0;
//     int i = 0, pos = 0;
    
//     *cseq = 0;
//     while (pos < seqLen)
//     {
//         *cseq = ((*cseq) << 3) | nVal[seq[pos++]];

//         if (++i == 21)
//         {
//             i = 0;
//             cseq++;
//             if (pos < seqLen)   // not to write the adjacent memory in case seqLen % 21 == 0
//                 *cseq = 0;
//         }
//     }
//     if (i > 0)
//     {
//         *cseq <<= (3*(21-i));
//     }
// }
/**********************************************/
// void decompressSequence(CompressedSeq *cseq, int seqLen, char *seq)
// {
//     int i;
//     int shifts = 20*3;

//     for (i = 0; i < seqLen; i++)
//     {
//         *(seq++) = _common_nRevVal[((*cseq) >> shifts) & 7];
//         shifts -= 3;
//         if (shifts < 0)
//         {
//                 cseq++;
//                 shifts = 20*3;
//         }
//     }
//     *seq = '\0';
// }
/**********************************************/
// int hashVal(char *seq)
// {
//     int i=0;
//     int val=0, numericVal=0;

//     while(i < MIN_ANCHOR_LEN)
//     {
//         if (nHVal[seq[i]] == -1)
//             return -1; 
//         val = (val << 2) | nHVal[seq[i++]]; 
//     }
//     return val;
// }
/**********************************************/
// void initCommon()
// {
//     memset(nVal, 4, 128);
//     nVal['A']=0;
//     nVal['C']=1;
//     nVal['G']=2;
//     nVal['T']=3;

//     memset(nRev, 'N', 128);
//     nRev['A']='T';
//     nRev['C']='G';
//     nRev['T']='A';
//     nRev['G']='C';

//     memset(nHVal, -1, 128);
//     nHVal['A']=0;
//     nHVal['C']=1;
//     nHVal['G']=2;
//     nHVal['T']=3;

//     _common_nRevVal[0] = 'A';
//     _common_nRevVal[1] = 'C';
//     _common_nRevVal[2] = 'G';
//     _common_nRevVal[3] = 'T';
//     _common_nRevVal[4] = 'N';
// }
/**********************************************/
// void reverseInPlace(char *dest, char *src, int len)
// {
//     int i;
//     src+=len-1;
//     for (i=0; i<len; i++)
//     {
//         *(dest++)=*(src--);
//     }
//     *dest='\0';
//     // dest[len]='\0';
// }