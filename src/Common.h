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

#ifndef __COMMON__
#define __COMMON__

#if SSE4==1 
#define MRSFAST_SSE4
#endif

#include <stdio.h>
#include <zlib.h>
#include <stdint.h>

#include <string>
#include <sstream>

#ifdef DODEBUG 
#define DEBUG(cmd) cmd
#else 
#define DEBUG(cmd)
#endif

#define SEQ_MAX_LENGTH		50000		// Seq Max Length
#define ALIGN_SEQ_MAX_LEN 	75000
#define CONTIG_OVERLAP		50400 		// No. of characters overlapped between contings  --  equals 2100 blocks of length 21
#define CMP_SEQ_MAX_LENGTH	10			// Compressed Seq Max Length
#define CONTIG_NAME_SIZE	200			// Contig name max size
// #define FILE_NAME_LENGTH	500			// Filename Max Length
#define MAX_SNP_PER_CHR		6000000
// #define SEED_MERGING_RANGE  (SAMPLING_COUNT+50)
#define SEED_MERGING_RANGE  500
#define MAX_QGRAM_DIS		30
#define MIN_SUP_THRESHOLD	0.4

typedef uint64_t CompressedSeq;
typedef uint16_t CheckSumType;

enum chainAlg_t { CHAIN_ALG_CLASP, CHAIN_ALG_DPN2, CHAIN_ALG_LISN2, CHAIN_ALG_LISNLOGN};

extern unsigned int		CONTIG_SIZE;
extern unsigned int		CONTIG_MAX_SIZE;
extern unsigned int		THREAD_COUNT;
extern double			MAX_MEMORY;
extern int				THREAD_ID[255];

extern unsigned char	WINDOW_SIZE;					// WINDOW SIZE for indexing/searching
extern unsigned int 	SAMPLING_INTERVAL;				// segment length
extern unsigned int 	SAMPLING_COUNT;					// seed count
extern unsigned int 	CHUNK_SIZE;
extern unsigned int 	CHUNK_OVERLAP;
extern unsigned short	SEQ_LENGTH;						// Sequence(read) length
extern uint16_t 		QGRAM_WIN_SIZE;
extern unsigned short	QUAL_LENGTH;
extern unsigned short	CMP_SEQ_LENGTH;
extern unsigned short	DISCORDANT_CUT_OFF;
extern int				SNP_QUAL_THRESHOLD;
extern unsigned int 	MAX_NUM_HITS;

extern int				indexingMode;
extern int				searchingMode;
extern int				bestMappingMode;
extern int 				affineMode;
extern int				seqCompressed;
extern int				outCompressed;
extern int				progressRep;
extern int				nohitDisabled;
extern int 				tabOutput;
extern int				noSamHeader;
extern char 			*seqFile;
extern char				*seqUnmapped;
extern char				outputMap[1000];
// extern char 			mappingOutputPath[1000];
extern char				outputUnmap[1000];
extern char 			opt_commandAll[2000];
extern int 				opt_outputBufferSize;
extern unsigned char	seqFastq;
extern int				errThreshold;
extern short			maxHits;
extern char				fileName[3][1000];
extern int				fileCnt;
extern long long		memUsage;
extern char				*alphabet;
extern chainAlg_t       chainAlg;

#pragma pack(push, 1)
typedef struct
{
	CheckSumType  checksum;
	uint32_t info;				// ReadIndex => seqInfo | GenomeIndex ==> Loc
} GeneralIndex;
#pragma pack(pop)
typedef struct
{
	int hv;
	GeneralIndex *list;
} ReadIndexTable;

typedef struct
{
	int loc;
	char alt;
} SNPLoc;

FILE	* fileOpen(char *fileName, const char *mode);
gzFile	fileOpenGZ(char *fileName, const char *mode);
double	getTime(void);
void	reverseComplete (char *seq, char *rcSeq , int length);
char reverseCompleteChar(char);
void	* getMem(size_t size);
void	freeMem(void * ptr, size_t size);
double	getMemUsage();
void 	reverse (char *seq, char *rcSeq , int length);
void 	stripPath(char *full, char *path, char *fileName);
void compressSequence(char *seq, int seqLen, CompressedSeq *cseq);
void decompressSequence(CompressedSeq *cseq, int seqLen, char *seq);
int 	calculateCompressedLen(int normalLen);
int	hashVal(char *seq);
int	checkSumVal(char *seq);
void initCommon();
void reverseInPlace(char *dest, char *src, int len);

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
