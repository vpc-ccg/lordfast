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

pthread_mutex_t         _common_lock;
char                    tableRev[128] = {
    'N','N','N','N',   'N','N','N','N',   'N','N','N','N',   'N','N','N','N', 
    'N','N','N','N',   'N','N','N','N',   'N','N','N','N',   'N','N','N','N', 
    'N','N','N','N',   'N','N','N','N',   'N','N','N','N',   'N','N','N','N', 
    'N','N','N','N',   'N','N','N','N',   'N','N','N','N',   'N','N','N','N', 
    'N','T','N','G',   'N','N','N','C',   'N','N','N','N',   'N','N','N','N', 
    'N','N','N','N',   'A','N','N','N',   'N','N','N','N',   'N','N','N','N', 
    'N','t','N','g',   'N','N','N','c',   'N','N','N','N',   'N','N','N','N', 
    'N','N','N','N',   'a','N','N','N',   'N','N','N','N',   'N','N','N','N'
};

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
void reverseComplement(char *seq, char *rcSeq , int length)
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
void* getMem(size_t size)
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
