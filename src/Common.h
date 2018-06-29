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

#ifndef __COMMON__
#define __COMMON__

#include <stdio.h>
#include <zlib.h>
#include <stdint.h>

#include <string>
#include <sstream>

#if VERBOSITY < 1
#define LOG1(cmd)
#define LOG2(cmd)
#define LOG3(cmd)
#elif VERBOSITY < 2
#define LOG1(cmd) cmd
#define LOG2(cmd)
#define LOG3(cmd)
#elif VERBOSITY < 3
#define LOG1(cmd) cmd
#define LOG2(cmd) cmd
#define LOG3(cmd)
#elif VERBOSITY >= 3
#define LOG1(cmd) cmd
#define LOG2(cmd) cmd
#define LOG3(cmd) cmd
#endif

#define SEQ_MAX_LENGTH      50000       // maximum sequence length

typedef uint64_t CompressedSeq;
typedef uint16_t CheckSumType;

enum chainAlg_t {CHAIN_ALG_CLASP, CHAIN_ALG_DPN2, CHAIN_ALG_LISN2, CHAIN_ALG_LISNLOGN};

extern int              indexingMode;
extern int              searchingMode;
extern int              affineMode;
extern int              noSamHeader;
extern char             *seqFile;
extern char             *refFile;
extern char             outputMap[1000];
extern char             opt_commandAll[2000];
extern int              opt_outputBufferSize;
extern chainAlg_t       chainAlg;
extern double           chainReward;
extern double           chainPenalty;
extern double           gapPenalty;
extern long long        memUsage;

extern int              THREAD_COUNT;
extern int              THREAD_ID[255];
extern int              MIN_ANCHOR_LEN;     // minimum anchor length for searching
extern int              SAMPLING_COUNT;     // number of anchoring positions
extern int              MAX_MAP;
extern int              MIN_READ_LEN;
extern int              MAX_REF_HITS;

gzFile  fileOpenGZ(char *fileName, const char *mode);
void    reverseComplement(char *seq, char *rcSeq , int length);
void*   getMem(size_t size);
void    freeMem(void * ptr, size_t size);
double  getMemUsage();
void    reverse (char *seq, char *rcSeq , int length);
double  getCpuTime();
double  getRealTime();

template <typename T>
T str2type(std::string str)
{
    T n;
    std::istringstream sin(str);
    sin >> n;
    return n;
}

template <typename T>
std::string type2str(T v)
{
    std::ostringstream sout;
    sout << v;
    return sout.str();
}

#endif
