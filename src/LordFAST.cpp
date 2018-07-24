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
#include <math.h>
#include <pthread.h>

#include <algorithm>
#include <string>
#include <vector>
#include <deque>
#include <sstream>
#include <time.h>

// #include "Output.h"
// #include "PacFAST-BWT.h"
#include "LordFAST.h"
#include "BWT.h"
#include "Chain.h"

#include "ksw.h"
#include "edlib.h"

pthread_t           *_pf_threads = NULL;
pthread_mutex_t     _pf_readIdLock;
pthread_mutex_t     _pf_outputLock;

uint32_t            _pf_seqPos;
Read                *_pf_seqList;
int                 _pf_seqListSize;

// SeedList            *_pf_seedsFR;
SeedList            *_pf_seedsForward;
SeedList            *_pf_seedsReverse;
SeedList            *_pf_seedsSelected;
// uint32_t            **_pf_bwtLocs_tmp;

WinList             *_pf_topWins;
ChainList           *_pf_topChains;
// Sam_t               **_pf_topMappings;
MapInfo             *_pf_topMappings;

WinCount_t          **_pf_refWin_cnt;
uint32_t            _pf_refWin_num;
// char                *_pf_refGenSeq;
uint32_t            _pf_refGenLen;
uint64_t            _pf_refGenBWTLen;

int                 _pf_maxWin = 10; // at least 2?!
int                 _pf_maxChain = 10;
// int                 _pf_maxCandidate = 10;

uint8_t             _pf_char2int[128];
int8_t              _pf_kswMatrix[25];
int8_t              _pf_kswMatrix_clip[25];
// TODO: tune score parameters
int8_t              _pf_kswMatch = 2; // 5
int8_t              _pf_kswMismatch = 5; // 6
int8_t              _pf_kswGapOpen = 2; // 10
int8_t              _pf_kswGapExtend = 1; // 0
int8_t              _pf_kswMatch_clip = 2;
int8_t              _pf_kswMismatch_clip = 16;
int8_t              _pf_kswGapOpen_clip = 0;
int8_t              _pf_kswGapExtend_clip = 1;
char                _pf_kswCigarTable[] = "MIDNSHP=X";
// parameters for split alignments
int                 _pf_clipLen = 500;
double              _pf_clipSim = 0.75;
int                 _pf_splitLen = 80;
double              _pf_splitSim = 0.40;
double              _pf_reverseSim = 0.60;

FILE                *_pf_outFile = NULL;

void*               mapSeq(void *idp);
int                 pf_getNextRead();
void                findTopWins_coarse(uint32_t read_len, SeedList *seeds, int isRev, int readIdx, int id);
// void                findTopWins2(uint32_t read_len, SeedList *seeds, int isRev, int readIdx, int id);
void                findTopWins_fine(uint32_t read_len, SeedList *seeds, int isRev, int readIdx, float minScore, int id);
// void                findTopWins4(uint32_t read_len, SeedList *seeds, int isRev, int readIdx, int id);
void                alignWin(Win_t &win, char *query, char *query_rev, uint32_t rLen, char *qual, char* qual_rev, SamList_t &map, int id);
void                fixCigarM(std::string &cigar, std::string semiCigar);
void                fixCigar(std::string &cigar, std::string semiCigar);
void                alignChain_ksw(Chain_t &chain, char *query, int32_t readLen, int isRev, SamList_t &map);
void                alignChain_edlib(Chain_t &chain, char *query, int32_t readLen, int isRev, SamList_t &map);
void                (*alignChain)(Chain_t &chain, char *query, int32_t readLen, int isRev, SamList_t &map);

/**********************************************/
void initializeFAST()
{
    int i, j, k;

    _pf_maxWin        = MAX_MAP;
    _pf_threads       = (pthread_t*) getMem(THREAD_COUNT * sizeof(pthread_t));
    _pf_refWin_cnt    = (WinCount_t**) getMem(THREAD_COUNT * sizeof(WinCount_t*));
    // _pf_seedsFR       = (SeedList*)  getMem(THREAD_COUNT * sizeof(SeedList));
    _pf_seedsForward  = (SeedList*)  getMem(THREAD_COUNT * sizeof(SeedList));
    _pf_seedsReverse  = (SeedList*)  getMem(THREAD_COUNT * sizeof(SeedList));
    _pf_seedsSelected = (SeedList*)  getMem(THREAD_COUNT * sizeof(SeedList));
//  _pf_bwtLocs_tmp   = (uint32_t**) getMem(THREAD_COUNT * sizeof(uint32_t*));
    _pf_topWins       = (WinList*)   getMem(THREAD_COUNT * sizeof(WinList));
    _pf_topChains     = (ChainList*) getMem(THREAD_COUNT * sizeof(ChainList));
    // _pf_topMappings   = (Sam_t**)    getMem(THREAD_COUNT * sizeof(Sam_t*));
    _pf_topMappings   = new MapInfo[THREAD_COUNT];

    // _pf_refGenBWTLen = bwt_get_refGenBWTLen();
    // _pf_refWin_num = _pf_refGenBWTLen / MIN_READ_LEN; // TODO: decide about number of windows based on the minimum read length
    _pf_refGenLen = bwt_get_refGenLen();
    _pf_refWin_num = _pf_refGenLen / MIN_READ_LEN; // TODO: decide about number of windows based on the minimum read length

    for(i=0; i<THREAD_COUNT; i++)
    {
        _pf_refWin_cnt[i]         = (WinCount_t*) getMem(_pf_refWin_num * sizeof(WinCount_t));
        // _pf_seedsFR[i].list       = (Seed_t*)   getMem(SAMPLING_COUNT * MAX_REF_HITS * sizeof(Seed_t));
        _pf_seedsForward[i].list  = (Seed_t*)   getMem(SAMPLING_COUNT * MAX_REF_HITS * sizeof(Seed_t));
        _pf_seedsReverse[i].list  = (Seed_t*)   getMem(SAMPLING_COUNT * MAX_REF_HITS * sizeof(Seed_t));
        _pf_seedsSelected[i].list = (Seed_t*)   getMem(SEQ_MAX_LENGTH * sizeof(Seed_t));
//      _pf_bwtLocs_tmp[i]        = (uint32_t*) getMem(MAX_REF_HITS * sizeof(uint32_t));
        _pf_topWins[i].list       = (Win_t*)    getMem(_pf_maxWin * sizeof(Win_t));
        _pf_topChains[i].list     = (Chain_t*)  getMem(_pf_maxChain * sizeof(Chain_t));
        for(j=0; j<_pf_maxChain; j++)
        {
            _pf_topChains[i].list[j].seeds = (Seed_t*) getMem(SEQ_MAX_LENGTH * sizeof(Seed_t));
        }
        // _pf_topMappings[i]        = (Sam_t*)    getMem(_pf_maxWin * sizeof(Sam_t));
        _pf_topMappings[i].mappings  = new SamList_t[_pf_maxWin];
        // for(j=0; j<_pf_maxWin; j++)
        // {
        //     _pf_topMappings[i][j].cigar = (char*) getMem(2 * SEQ_MAX_LENGTH * sizeof(char));
        // }
    }

//  // get referenece 
//  _pf_refGenSeq = bwt_get_refGen_seq();
//  _pf_refGenLen = bwt_get_refGen_len();

    // treat all non-ACGT nucleotids as N's
    // A=0, C=1, G=2, T=3, o.w.=4
    memset(_pf_char2int, 4, 128 * sizeof(uint8_t));
    _pf_char2int['A'] = _pf_char2int['a'] = 0;
    _pf_char2int['C'] = _pf_char2int['c'] = 1;
    _pf_char2int['G'] = _pf_char2int['g'] = 2;
    _pf_char2int['T'] = _pf_char2int['t'] = 3;

    // initialize scoring matrix (from ksw.c)
    for (i=0, k=0; i<4; ++i)
    {
        for (j=0; j<4; ++j)
        {
            _pf_kswMatrix[k++] = (i==j ? _pf_kswMatch : -_pf_kswMismatch);
        }
        _pf_kswMatrix[k++] = 0; // ambiguous base
    }
    for (j=0; j<5; ++j)
        _pf_kswMatrix[k++] = 0;
    // initialize scoring matrix (from ksw.c)
    for (i=0, k=0; i<4; ++i)
    {
        for (j=0; j<4; ++j)
        {
            _pf_kswMatrix_clip[k++] = (i==j ? _pf_kswMatch_clip : -_pf_kswMismatch_clip);
        }
        _pf_kswMatrix_clip[k++] = 0; // ambiguous base
    }
    for (j=0; j<5; ++j)
        _pf_kswMatrix_clip[k++] = 0;

    // initialize output
    if(strlen(outputMap) == 0)
    {
        _pf_outFile = stdout;
    }
    else
    {
        _pf_outFile = fopen(outputMap, "w");
        if(_pf_outFile == NULL)
        {
            fprintf(stderr, "[ERROR] (initializeFAST) could not open file %s\n", outputMap);
            exit(EXIT_FAILURE);
        }
    }

    if(!noSamHeader)
        printSamHeader(_pf_outFile);

    // if(affineMode)
    // {
    //     alignChain = &alignChain_ksw;
    // }
    // else
    // {
        alignChain = &alignChain_edlib;
    // }
}
/**********************************************/
void finalizeFAST()
{
    int i, j;

    for(i=0; i<THREAD_COUNT; i++)
    {
        for(j=0; j<_pf_maxChain; j++)
        {
            freeMem(_pf_topChains[i].list[j].seeds, SEQ_MAX_LENGTH * sizeof(Seed_t));
        }
        // for(j=0; j<_pf_maxWin; j++)
        // {
        //     freeMem(_pf_topMappings[i][j].cigar, 2 * SEQ_MAX_LENGTH * sizeof(char));
        // }
        freeMem(_pf_refWin_cnt[i], _pf_refWin_num * sizeof(WinCount_t));
        // freeMem(_pf_seedsFR[i].list, SAMPLING_COUNT * MAX_REF_HITS * sizeof(Seed_t));
        freeMem(_pf_seedsForward[i].list, SAMPLING_COUNT * MAX_REF_HITS * sizeof(Seed_t));
        freeMem(_pf_seedsReverse[i].list, SAMPLING_COUNT * MAX_REF_HITS * sizeof(Seed_t));
        freeMem(_pf_seedsSelected[i].list, SEQ_MAX_LENGTH * sizeof(Seed_t));
//      freeMem(_pf_bwtLocs_tmp[i], MAX_REF_HITS * sizeof(uint32_t));
        freeMem(_pf_topWins[i].list, _pf_maxWin * sizeof(Win_t));
        freeMem(_pf_topChains[i].list, _pf_maxChain * sizeof(Chain_t));
        // freeMem(_pf_topMappings[i], _pf_maxWin * sizeof(Sam_t));
        delete[] _pf_topMappings[i].mappings;
    }

    freeMem(_pf_threads, THREAD_COUNT * sizeof(pthread_t));
    freeMem(_pf_refWin_cnt, THREAD_COUNT * sizeof(WinCount_t*));
    // freeMem(_pf_seedsFR, THREAD_COUNT * sizeof(SeedList));
    freeMem(_pf_seedsForward, THREAD_COUNT * sizeof(SeedList));
    freeMem(_pf_seedsReverse, THREAD_COUNT * sizeof(SeedList));
    freeMem(_pf_seedsSelected, THREAD_COUNT * sizeof(SeedList));
//  freeMem(_pf_bwtLocs_tmp, THREAD_COUNT * sizeof(uint32_t*));
    freeMem(_pf_topWins, THREAD_COUNT * sizeof(WinList));
    freeMem(_pf_topChains, THREAD_COUNT * sizeof(ChainList));
    // freeMem(_pf_topMappings, THREAD_COUNT * sizeof(Sam_t*));
    delete[] _pf_topMappings;

    // finalize output
    if(strlen(outputMap) > 0)
        fclose(_pf_outFile);
}
/**********************************************/
void initFASTChunk(Read *seqList, int seqListSize)
{
    int i, j;

    _pf_seqPos = 0;
    _pf_seqList = seqList;
    _pf_seqListSize = seqListSize;

    // memset(_pf_refWin_cnt[id], 0, _pf_refWin_num * sizeof(uint32_t));
    for(i = 0; i < THREAD_COUNT; i++)
        for(j = 0; j < _pf_refWin_num; j++)
            _pf_refWin_cnt[i][j].readIdx = 2 * _pf_seqListSize + 10;
            // _pf_refWin_cnt[i][j].readIdx = _pf_seqListSize + 100;

    // _pf_chainList = (Chain_t**) getMem(_pf_seqListSize * sizeof(Chain_t*));
    // _pf_chainListSize = (uint8_t*) getMem(_pf_seqListSize * sizeof(uint8_t));
    // for(i=0; i<_pf_seqListSize; i++)
    // {
    //  _pf_chainList[i] = (Chain_t*) getMem(_pf_maxWin * sizeof(Chain_t));
    //  _pf_chainListSize[i] = 0;
    // }
}
/**********************************************/
// void finalizeFASTChunk()
// {
//     int i;

//     // for(i=0; i<_pf_seqListSize; i++)
//     // {
//     //  freeMem(_pf_chainList[i], _pf_maxWin * sizeof(Chain_t));
//     // }
//     // freeMem(_pf_chainList, _pf_seqListSize * sizeof(Chain_t*));
//     // freeMem(_pf_chainListSize, _pf_seqListSize * sizeof(uint8_t));
// }
/**********************************************/
int pf_getNextRead()
{
    int x;
    pthread_mutex_lock(&_pf_readIdLock);
        x = _pf_seqPos;
        _pf_seqPos++;
    pthread_mutex_unlock(&_pf_readIdLock);
    return x;
}
/**********************************************/
void mapSeqMT()
{
    int i;

    for (i = 0; i < THREAD_COUNT; i++)
        pthread_create(_pf_threads + i, NULL, mapSeq, THREAD_ID + i);
    
    for (i = 0; i < THREAD_COUNT; i++)
        pthread_join(_pf_threads[i], NULL);

    return;
}
/*********************************************/
void printSamEntry(MapInfo &map, int readLen, int num, std::ostringstream& sout)
{
    int i, j;
    uint32_t chrBeg, chrEnd;
    char *chrName;
    int32_t chrLen;

    double bestEdit = (num > 0 ? (double)(-1*map.mappings[0].totalScore)/readLen : 1); // edit distance percentage of the best alingment
    double mapqPortion = 50.0 / (_pf_maxWin - 1);
    int x1 = 0; // number of mappings 
    int x2 = 0; // number of mappings with edit distance percentage close to that of the best mapping
    for(i = 0; i < num; i++)
    {
        if(map.mappings[i].samList.size() > 0) // mapped
        {
            x1++;
            if((double)(-1*map.mappings[i].totalScore)/readLen * 0.95 < bestEdit)
                x2++;
        }
    }
    double mapq = (x2 > 1 ? 2.1 : (_pf_maxWin - x1) * mapqPortion); // mapq 0 if there are multiple good alignments
    int32_t mapq_int;
    double lowEditBonus, betterQualBonus;
    std::string tmpStr;

    for(i = 0; i < num; i++)
    {
        if(i == 0) // first mapping should be printed anyway (mapped or un-mapped)
        {
            if(map.mappings[i].samList.size() > 0) // mapped
            {
                if(num == 1 || (num > 1 && (double)(-1*map.mappings[i].totalScore)/readLen < 0.15 && (double)(-1*map.mappings[i].totalScore)/readLen < 0.95 * (double)(-1*map.mappings[i+1].totalScore)/readLen))
                {
                    mapq_int = 60;
                }
                else 
                {
                    mapq_int = mapq + 5 * (0.2 - (double)(-1*map.mappings[i].totalScore)/readLen) / 0.2;
                }

                for(j = 0; j < map.mappings[i].samList.size(); j++)
                {
                    bwt_get_intv_info(map.mappings[i].samList[j].pos, map.mappings[i].samList[j].posEnd, &chrName, &chrLen, &chrBeg, &chrEnd);
                    //
                    map.mappings[i].samList[j].rName = chrName;
                    map.mappings[i].samList[j].rStart = chrBeg;
                    map.mappings[i].samList[j].mapQ = mapq_int;
                    tmpStr = "";
                    tmpStr += chrName; tmpStr += ",";
                    tmpStr += type2str<uint32_t>(chrBeg + 1) + ",";
                    tmpStr += (map.mappings[i].samList[j].flag & 16 ? "-," : "+,");
                    tmpStr += map.mappings[i].samList[j].cigar + ",";
                    tmpStr += type2str<int32_t>(mapq_int) + ",";
                    tmpStr += type2str<int32_t>(abs(map.mappings[i].samList[j].nmCount)) + ";";
                    map.mappings[i].samList[j].sa = tmpStr;
                }

                for(j = 0; j < map.mappings[i].samList.size(); j++)
                {
                    sout<< map.qName << "\t" // QNAME
                        << (j > 0 ? (map.mappings[i].samList[j].flag | 2048) : map.mappings[i].samList[j].flag) << "\t" // FLAG
                        << map.mappings[i].samList[j].rName << "\t" // RNAME
                        << map.mappings[i].samList[j].rStart + 1 << "\t" // POS
                        << (mapq_int >= 0 ? mapq_int : 0) << "\t" // MAPQ
                        << map.mappings[i].samList[j].cigar << "\t*\t0\t0\t" // CIGAR, RNEXT, PNEXT, TLEN
                        << (map.mappings[i].samList[j].flag & 16 ? map.seq_rev : map.seq) << "\t" // SEQ
                        << (map.mappings[i].samList[j].flag & 16 ? map.qual_rev : map.qual) << "\t" // QUAL
                        << "AS:i:" << map.mappings[i].samList[j].alnScore << "\t" 
                        << "XS:i:" << 0 << "\t"
                        << "NM:i:" << abs(map.mappings[i].samList[j].nmCount) << "\t"
                        << "MD:Z:" << map.mappings[i].samList[j].md;
                        if(strlen(readGroupId) > 0)
                            sout<< "\tRG:Z:" << readGroupId;
                        if(map.mappings[i].samList.size() > 1) // has supplementary
                        {
                            sout<< "\tSA:Z:";
                            for(int z = 0; z < map.mappings[i].samList.size(); z++)
                            {
                                if(z != j)
                                {
                                    sout<< map.mappings[i].samList[z].sa;
                                }
                            }
                        }
                        sout<< "\n";
                        // << "TS:i:" << map.mappings[i].totalScore << "\t"
                        // << "QS:i:" << map.mappings[i].samList[j].qStart << "\t"
                        // << "QE:i:" << map.mappings[i].samList[j].qEnd << "\n";
                }
            }
            else // un-mapped
            {
                sout<< map.qName  << "\t4\t*\t0\t0\t*\t*\t0\t0\t" 
                    << map.seq << "\t" << map.qual;
                if(strlen(readGroupId) > 0)
                    sout<< "\tRG:Z:" << readGroupId;
                sout << "\n";
            }
        }
        else
        {
            if(map.mappings[i].samList.size() > 0) // mapped
            {
                //
                mapq_int = mapq + 5 * (0.2 - (double)(-1*map.mappings[i].totalScore)/readLen) / 0.2;
                //
                for(j = 0; j < map.mappings[i].samList.size(); j++)
                {
                    bwt_get_intv_info(map.mappings[i].samList[j].pos, map.mappings[i].samList[j].posEnd, &chrName, &chrLen, &chrBeg, &chrEnd);
                    //
                    sout<< map.qName << "\t" // QNAME
                        << (map.mappings[i].samList[j].flag | 256)  << "\t" // FLAG
                        << chrName << "\t" // RNAME
                        << chrBeg + 1 << "\t" // POS
                        << (mapq_int >= 0 ? mapq_int : 0) << "\t" // MAPQ
                        << map.mappings[i].samList[j].cigar << "\t*\t0\t0\t" // CIGAR, RNEXT, PNEXT, TLEN
                        << (map.mappings[i].samList[j].flag & 16 ? map.seq_rev : map.seq) << "\t" // SEQ
                        << (map.mappings[i].samList[j].flag & 16 ? map.qual_rev : map.qual) << "\t" // QUAL
                        << "AS:i:" << map.mappings[i].samList[j].alnScore << "\t" 
                        << "XS:i:" << 0 << "\t"
                        << "NM:i:" << abs(map.mappings[i].samList[j].nmCount) << "\t"
                        << "MD:Z:" << map.mappings[i].samList[j].md;
                        if(strlen(readGroupId) > 0)
                            sout<< "\tRG:Z:" << readGroupId;
                        sout<< "\n";
                        // << "TS:i:" << map.mappings[i].totalScore << "\t"
                        // << "QS:i:" << map.mappings[i].samList[j].qStart << "\t"
                        // << "QE:i:" << map.mappings[i].samList[j].qEnd << "\n";
                }
            }
        }
    }
    //
    if(sout.tellp() > opt_outputBufferSize)
    {
        pthread_mutex_lock(&_pf_outputLock);
        fprintf(_pf_outFile, "%s", sout.str().c_str());
        pthread_mutex_unlock(&_pf_outputLock);
        sout.str("");
        sout.clear();
    }
}
/*********************************************/
void* mapSeq(void *idp)
{
    int id = *(int*)idp;
    int i, j;
    int t;
    Read *read;
    uint32_t readLen;
    uint32_t qualLen;
    char seq_rev[SEQ_MAX_LENGTH];
    char qual_rev[SEQ_MAX_LENGTH];
    std::ostringstream outBuffer;

    // clock_t clk;

    while ((t=pf_getNextRead())<_pf_seqListSize)
    {
        // Prcessing read t
        read = _pf_seqList + t;
        readLen = (*read->length);
        qualLen = (*read->isFq ? *read->length : 1);

        _pf_topMappings[id].qName  = read->name;
        _pf_topMappings[id].seq    = read->seq;
        _pf_topMappings[id].qual   = read->qual;

        LOG1({
            fprintf(stderr, ">%s\tlen: %d\n", read->name, readLen);
        });

        if(readLen < MIN_READ_LEN)
        {
            LOG1({
                fprintf(stderr, "\tunmapped\tshortRead\n", read->name, readLen);
            });

            _pf_topMappings[id].mappings[0].samList.clear();
            printSamEntry(_pf_topMappings[id], readLen, 1, outBuffer);
            continue;
        }

        reverseComplement(read->seq, seq_rev, *read->length);
        reverse(read->qual, qual_rev, qualLen);
        _pf_topMappings[id].seq_rev  = seq_rev;
        _pf_topMappings[id].qual_rev = qual_rev;

        // extract forward-reverse seeds
        getLocs_extend_whole_step(read->seq, readLen, SAMPLING_COUNT, _pf_seedsForward + id, _pf_seedsReverse + id);
        // getLocs_extend_whole_step2(read->seq, readLen, SAMPLING_COUNT, _pf_seedsForward + id, _pf_seedsReverse + id);
        // getLocs_extend_whole_step3(read->seq, readLen, SAMPLING_COUNT, _pf_seedsForward + id, _pf_seedsReverse + id);
        // continue; // just do seeding
        // reset number of top windows
        _pf_topWins[id].num = 0;
        findTopWins_coarse(readLen, _pf_seedsForward + id, 0, t+1, id); //find candidate paths for forward
        findTopWins_coarse(readLen, _pf_seedsReverse + id, 1, -(t+1), id); //find candidate paths for reverse

        // sort the windows based on the score!
        std::sort_heap(_pf_topWins[id].list, _pf_topWins[id].list + _pf_topWins[id].num, compareWin);
        
        LOG1({
            int ilog;
            // fprintf(stderr, "\tcandidate\tnum: %d\n", _pf_topWins[id].num);
            for(ilog = 0; ilog < _pf_topWins[id].num; ilog++)
            {
                fprintf(stderr, "\tcandidate %d:\t%u\t%u\t%c\t%f\n", ilog+1, _pf_topWins[id].list[ilog].tStart, _pf_topWins[id].list[ilog].tEnd,
                (_pf_topWins[id].list[ilog].isReverse ? '-' : '+'), _pf_topWins[id].list[ilog].score);
            }
        });

        // continue; // seeding and candidate selection

        float scoreRatio = 4;
        if(_pf_topWins[id].list[0].score >= scoreRatio * _pf_topWins[id].list[1].score) // significantly better window => coarse mode
        {
            _pf_topMappings[id].mappings[0].samList.clear();
            // _pf_topMappings[id].mappings[0].totalScore = 0;
            alignWin(_pf_topWins[id].list[0], read->seq, seq_rev, readLen, read->qual, qual_rev, _pf_topMappings[id].mappings[0], id);
            printSamEntry(_pf_topMappings[id], readLen, 1, outBuffer);
        }
        else // fine mode
        {
            _pf_topWins[id].num = 0;
            float tmp_minScore = (float)_pf_topWins[id].list[0].score / scoreRatio;
            findTopWins_fine(readLen, _pf_seedsForward + id, 0, t+_pf_seqListSize+1, tmp_minScore, id); //find candidate paths for forward
            findTopWins_fine(readLen, _pf_seedsReverse + id, 1, -(t+_pf_seqListSize+1), tmp_minScore, id); //find candidate paths for reverse

            for(i=0; i<_pf_topWins[id].num; i++)
            {
                _pf_topMappings[id].mappings[i].samList.clear();
                // _pf_topMappings[id].mappings[i].totalScore = 0;
                alignWin(_pf_topWins[id].list[i], read->seq, seq_rev, readLen, read->qual, qual_rev, _pf_topMappings[id].mappings[i], id);
            }

            // sort the windows based on the score!
            std::sort(_pf_topMappings[id].mappings, _pf_topMappings[id].mappings + _pf_topWins[id].num, compareSam);

            // print 
            printSamEntry(_pf_topMappings[id], readLen, _pf_topWins[id].num, outBuffer);
        }
    }

    if(outBuffer.tellp() > 0)
    {
        pthread_mutex_lock(&_pf_outputLock);
        fprintf(_pf_outFile, "%s", outBuffer.str().c_str());
        pthread_mutex_unlock(&_pf_outputLock);
    }

    return NULL;
}
/**********************************************/
void findTopWins_coarse(uint32_t read_len, SeedList *seeds, int isRev, int readIdx, int id)
{
    int i, winNumLimit;

    // memset(_pf_refWin_cnt[id], 0, _pf_refWin_num * sizeof(uint32_t));

    for(i=0; i<seeds->num; i++)
    {
        int32_t winId = seeds->list[i].tPos / read_len;

        // int32_t weight = (seeds->list[i].len - MIN_ANCHOR_LEN); //  979:1 4:2 1:3 16:-1
        int32_t weight = 1 + (seeds->list[i].len - MIN_ANCHOR_LEN); // 978:1 6:2 1:3 15:-1
        // int32_t weight = 1 + 2 * (seeds->list[i].len - MIN_ANCHOR_LEN); //  978:1 5:2 1:3 16:-1
        // int32_t weight = 1 + 5 * (seeds->list[i].len - MIN_ANCHOR_LEN); //  979:1 5:2 1:3 15:-1

        if(_pf_refWin_cnt[id][winId].readIdx == readIdx)
        {
            // fprintf(stderr, "@\tYES\t%d\t%d\n", _pf_refWin_cnt[id][winId].readIdx, readIdx);
            _pf_refWin_cnt[id][winId].cnt += weight;
        }
        else
        {
            // fprintf(stderr, "@\tNO\t%d\t%d\n", _pf_refWin_cnt[id][winId].readIdx, readIdx);
            _pf_refWin_cnt[id][winId].readIdx = readIdx;
            _pf_refWin_cnt[id][winId].cnt = weight;
        }
        if(winId-1 >= 0)
        {
            if(_pf_refWin_cnt[id][winId-1].readIdx == readIdx)
            {
                _pf_refWin_cnt[id][winId-1].cnt += weight;
            }
            else
            {
                _pf_refWin_cnt[id][winId-1].readIdx = readIdx;
                _pf_refWin_cnt[id][winId-1].cnt = weight;
            }
        }
    }

    winNumLimit = _pf_refGenLen / read_len + 2;
    if(winNumLimit > _pf_refWin_num)
        winNumLimit = _pf_refWin_num;

    // find top _pf_maxWin wins using a heap
    // for(i = 0; i < _pf_refWin_num; i++)
    for(i = 0; i < winNumLimit; i++)
    {
        if(_pf_refWin_cnt[id][i].readIdx == readIdx && 
            (i==0 || _pf_refWin_cnt[id][i-1].readIdx != readIdx || _pf_refWin_cnt[id][i].cnt >= _pf_refWin_cnt[id][i-1].cnt) &&
            (i==_pf_refWin_num-1 || _pf_refWin_cnt[id][i+1].readIdx != readIdx || _pf_refWin_cnt[id][i].cnt > _pf_refWin_cnt[id][i+1].cnt) )
        {
            if(_pf_topWins[id].num < _pf_maxWin) // the list has some space, push the chains 
            {
                _pf_topWins[id].list[_pf_topWins[id].num].tStart = i*read_len;
                _pf_topWins[id].list[_pf_topWins[id].num].tEnd = (i+2)*read_len-1;
                _pf_topWins[id].list[_pf_topWins[id].num].score = _pf_refWin_cnt[id][i].cnt;
                _pf_topWins[id].list[_pf_topWins[id].num++].isReverse = isRev;
                std::push_heap(_pf_topWins[id].list, _pf_topWins[id].list + _pf_topWins[id].num, compareWin);
            }
            else
            {
                if( _pf_refWin_cnt[id][i].cnt > _pf_topWins[id].list[0].score) // if the score is higher than the current smallest score, replace it with that
                {
                    std::pop_heap(_pf_topWins[id].list, _pf_topWins[id].list + _pf_topWins[id].num, compareWin); // the smallest score chain is at the last position
                    // update the last chain
                    _pf_topWins[id].list[_pf_topWins[id].num-1].tStart = i*read_len;
                    _pf_topWins[id].list[_pf_topWins[id].num-1].tEnd = (i+2)*read_len-1;
                    _pf_topWins[id].list[_pf_topWins[id].num-1].score = _pf_refWin_cnt[id][i].cnt;
                    _pf_topWins[id].list[_pf_topWins[id].num-1].isReverse = isRev;
                    std::push_heap(_pf_topWins[id].list, _pf_topWins[id].list + _pf_topWins[id].num, compareWin);
                }
            }
        }
    }
}
/**********************************************/
float calcChainScore(uint32_t rLen, uint32_t tStart, uint32_t tEnd, int isReverse, int id)
{
    uint32_t i;
    uint32_t margin = rLen >> 1;
    int64_t seedPos_Low = (int64_t)tStart - margin;
    int64_t seedPos_High = (int64_t)tEnd + margin;
    float retScore;

    if(isReverse)
    {
        _pf_seedsSelected[id].num = 0;
        for(i = 0; i < _pf_seedsReverse[id].num; i++)
        {
            if((int64_t)_pf_seedsReverse[id].list[i].tPos >= seedPos_Low && (int64_t)_pf_seedsReverse[id].list[i].tPos <= seedPos_High)
            {
                _pf_seedsSelected[id].list[_pf_seedsSelected[id].num++] = _pf_seedsReverse[id].list[i];
            }
        }

        if(chainAlg == CHAIN_ALG_CLASP)
        {
            if(seedPos_Low > 2000000000)
                for(i = 0; i < _pf_seedsSelected[id].num; i++)
                    _pf_seedsSelected[id].list[i].tPos -= 2000000000;

            chain_seeds_clasp(_pf_seedsSelected[id].list, _pf_seedsSelected[id].num, _pf_topChains[id].list[0]);
            
            if(seedPos_Low > 2000000000)
                for(i = 0; i < _pf_topChains[id].list[0].chainLen; i++)
                    _pf_topChains[id].list[0].seeds[i].tPos += 2000000000;
        }
        else if(chainAlg == CHAIN_ALG_DPN2)
        {
            chain_seeds_n2(_pf_seedsSelected[id].list, _pf_seedsSelected[id].num, _pf_topChains[id].list[0]);
        }
        
        retScore = _pf_topChains[id].list[0].score;
    }
    else
    {
        _pf_seedsSelected[id].num = 0;
        for(i = 0; i < _pf_seedsForward[id].num; i++)
        {
            if((int64_t)_pf_seedsForward[id].list[i].tPos >= seedPos_Low && (int64_t)_pf_seedsForward[id].list[i].tPos <= seedPos_High)
            {
                _pf_seedsSelected[id].list[_pf_seedsSelected[id].num++] = _pf_seedsForward[id].list[i];
            }
        }

        if(chainAlg == CHAIN_ALG_CLASP)
        {
            if(seedPos_Low > 2000000000)
                for(i = 0; i < _pf_seedsSelected[id].num; i++)
                    _pf_seedsSelected[id].list[i].tPos -= 2000000000;

            chain_seeds_clasp(_pf_seedsSelected[id].list, _pf_seedsSelected[id].num, _pf_topChains[id].list[0]);

            if(seedPos_Low > 2000000000)
                for(i = 0; i < _pf_topChains[id].list[0].chainLen; i++)
                    _pf_topChains[id].list[0].seeds[i].tPos += 2000000000;
        }
        else if(chainAlg == CHAIN_ALG_DPN2)
        {
            chain_seeds_n2(_pf_seedsSelected[id].list, _pf_seedsSelected[id].num, _pf_topChains[id].list[0]);
        }

        retScore = _pf_topChains[id].list[0].score;
    }
    return retScore;
}
/**********************************************/
// void findTopWins2(uint32_t read_len, SeedList *seeds, int isRev, int readIdx, int id)
// {
//  int i, winNumLimit;
//  float tmpScore;

//  // memset(_pf_refWin_cnt[id], 0, _pf_refWin_num * sizeof(uint32_t));

//  for(i=0; i<seeds->num; i++)
//  {
//      int32_t winId = seeds->list[i].tPos / read_len;

//      // int32_t weight = (seeds->list[i].len - MIN_ANCHOR_LEN); //  979:1 4:2 1:3 16:-1
//      int32_t weight = 1 + (seeds->list[i].len - MIN_ANCHOR_LEN); // 978:1 6:2 1:3 15:-1
//      // int32_t weight = 1 + 2 * (seeds->list[i].len - MIN_ANCHOR_LEN); //  978:1 5:2 1:3 16:-1
//      // int32_t weight = 1 + 5 * (seeds->list[i].len - MIN_ANCHOR_LEN); //  979:1 5:2 1:3 15:-1

//      if(_pf_refWin_cnt[id][winId].readIdx == readIdx)
//      {
//          // fprintf(stderr, "@\tYES\t%d\t%d\n", _pf_refWin_cnt[id][winId].readIdx, readIdx);
//          _pf_refWin_cnt[id][winId].cnt += weight;
//      }
//      else
//      {
//          // fprintf(stderr, "@\tNO\t%d\t%d\n", _pf_refWin_cnt[id][winId].readIdx, readIdx);
//          _pf_refWin_cnt[id][winId].readIdx = readIdx;
//          _pf_refWin_cnt[id][winId].cnt = weight;
//      }
//      if(winId-1 >= 0)
//      {
//          if(_pf_refWin_cnt[id][winId-1].readIdx == readIdx)
//          {
//              _pf_refWin_cnt[id][winId-1].cnt += weight;
//          }
//          else
//          {
//              _pf_refWin_cnt[id][winId-1].readIdx = readIdx;
//              _pf_refWin_cnt[id][winId-1].cnt = weight;
//          }
//      }
//  }

//  winNumLimit = _pf_refGenLen / read_len + 2;
//  if(winNumLimit > _pf_refWin_num)
//      winNumLimit = _pf_refWin_num;

//  // // put a dummy entry in _pf_topWins[id]
//  // _pf_topWins[id].list[0].tStart = 0;
//  // _pf_topWins[id].list[0].tEnd = 0;
//  // _pf_topWins[id].list[0].score = -1;
//  // _pf_topWins[id].list[0].isReverse = isRev;
//  // _pf_topWins[id].num = 1;
//  // std::push_heap(_pf_topWins[id].list, _pf_topWins[id].list + _pf_topWins[id].num, compareWin);

//  // find top _pf_maxWin wins using a heap
//  for(i = 0; i < winNumLimit; i++)
//  {
//      if(_pf_refWin_cnt[id][i].readIdx == readIdx)
//      {
//          // calculate a chain for this window
//          tmpScore = calcChainScore(read_len, i*read_len, (i+2)*read_len-1, isRev, id);
//          if(_pf_topWins[id].num < _pf_maxWin) // the list has some space, push the chains 
//          {
//              _pf_topWins[id].list[_pf_topWins[id].num].tStart = i*read_len;
//              _pf_topWins[id].list[_pf_topWins[id].num].tEnd = (i+2)*read_len-1;
//              _pf_topWins[id].list[_pf_topWins[id].num].score = tmpScore;
//              _pf_topWins[id].list[_pf_topWins[id].num++].isReverse = isRev;
//              std::push_heap(_pf_topWins[id].list, _pf_topWins[id].list + _pf_topWins[id].num, compareWin);
//          }
//          else
//          {
//              if(tmpScore > _pf_topWins[id].list[0].score) // if the score is higher than the current smallest score, replace it with that
//              {
//                  std::pop_heap(_pf_topWins[id].list, _pf_topWins[id].list + _pf_topWins[id].num, compareWin); // the smallest score chain is at the last position
//                  // update the last chain
//                  _pf_topWins[id].list[_pf_topWins[id].num-1].tStart = i*read_len;
//                  _pf_topWins[id].list[_pf_topWins[id].num-1].tEnd = (i+2)*read_len-1;
//                  _pf_topWins[id].list[_pf_topWins[id].num-1].score = tmpScore;
//                  _pf_topWins[id].list[_pf_topWins[id].num-1].isReverse = isRev;
//                  std::push_heap(_pf_topWins[id].list, _pf_topWins[id].list + _pf_topWins[id].num, compareWin);
//              }
//          }
//      }
//  }
// }
/**********************************************/
void findTopWins_fine(uint32_t read_len, SeedList *seeds, int isRev, int readIdx, float minScore, int id)
{
    int i, winNumLimit;
    float tmpScore;

    // memset(_pf_refWin_cnt[id], 0, _pf_refWin_num * sizeof(uint32_t));

    for(i=0; i<seeds->num; i++)
    {
        int32_t winId = seeds->list[i].tPos / read_len;

        // int32_t weight = (seeds->list[i].len - MIN_ANCHOR_LEN); //  979:1 4:2 1:3 16:-1
        int32_t weight = 1 + (seeds->list[i].len - MIN_ANCHOR_LEN); // 978:1 6:2 1:3 15:-1
        // int32_t weight = 1 + 2 * (seeds->list[i].len - MIN_ANCHOR_LEN); //  978:1 5:2 1:3 16:-1
        // int32_t weight = 1 + 5 * (seeds->list[i].len - MIN_ANCHOR_LEN); //  979:1 5:2 1:3 15:-1

        if(_pf_refWin_cnt[id][winId].readIdx == readIdx)
        {
            // fprintf(stderr, "@\tYES\t%d\t%d\n", _pf_refWin_cnt[id][winId].readIdx, readIdx);
            _pf_refWin_cnt[id][winId].cnt += weight;
        }
        else
        {
            // fprintf(stderr, "@\tNO\t%d\t%d\n", _pf_refWin_cnt[id][winId].readIdx, readIdx);
            _pf_refWin_cnt[id][winId].readIdx = readIdx;
            _pf_refWin_cnt[id][winId].cnt = weight;
        }
        if(winId-1 >= 0)
        {
            if(_pf_refWin_cnt[id][winId-1].readIdx == readIdx)
            {
                _pf_refWin_cnt[id][winId-1].cnt += weight;
            }
            else
            {
                _pf_refWin_cnt[id][winId-1].readIdx = readIdx;
                _pf_refWin_cnt[id][winId-1].cnt = weight;
            }
        }
    }

    winNumLimit = _pf_refGenLen / read_len + 2;
    if(winNumLimit > _pf_refWin_num)
        winNumLimit = _pf_refWin_num;

    // // put a dummy entry in _pf_topWins[id]
    // _pf_topWins[id].list[0].tStart = 0;
    // _pf_topWins[id].list[0].tEnd = 0;
    // _pf_topWins[id].list[0].score = -1;
    // _pf_topWins[id].list[0].isReverse = isRev;
    // _pf_topWins[id].num = 1;
    // std::push_heap(_pf_topWins[id].list, _pf_topWins[id].list + _pf_topWins[id].num, compareWin);

    // find top _pf_maxWin wins using a heap
    for(i = 0; i < winNumLimit; i++)
    {
        if(_pf_refWin_cnt[id][i].readIdx == readIdx && _pf_refWin_cnt[id][i].cnt > minScore &&
        (i==0 || _pf_refWin_cnt[id][i-1].readIdx != readIdx || _pf_refWin_cnt[id][i].cnt >= _pf_refWin_cnt[id][i-1].cnt) &&
        (i==_pf_refWin_num-1 || _pf_refWin_cnt[id][i+1].readIdx != readIdx || _pf_refWin_cnt[id][i].cnt > _pf_refWin_cnt[id][i+1].cnt) )
        {
            // calculate a chain for this window
            tmpScore = calcChainScore(read_len, i*read_len, (i+2)*read_len-1, isRev, id);
            if(_pf_topWins[id].num < _pf_maxWin) // the list has some space, push the chains 
            {
                _pf_topWins[id].list[_pf_topWins[id].num].tStart = i*read_len;
                _pf_topWins[id].list[_pf_topWins[id].num].tEnd = (i+2)*read_len-1;
                _pf_topWins[id].list[_pf_topWins[id].num].score = tmpScore;
                _pf_topWins[id].list[_pf_topWins[id].num++].isReverse = isRev;
                std::push_heap(_pf_topWins[id].list, _pf_topWins[id].list + _pf_topWins[id].num, compareWin);
            }
            else
            {
                if(tmpScore > _pf_topWins[id].list[0].score) // if the score is higher than the current smallest score, replace it with that
                {
                    std::pop_heap(_pf_topWins[id].list, _pf_topWins[id].list + _pf_topWins[id].num, compareWin); // the smallest score chain is at the last position
                    // update the last chain
                    _pf_topWins[id].list[_pf_topWins[id].num-1].tStart = i*read_len;
                    _pf_topWins[id].list[_pf_topWins[id].num-1].tEnd = (i+2)*read_len-1;
                    _pf_topWins[id].list[_pf_topWins[id].num-1].score = tmpScore;
                    _pf_topWins[id].list[_pf_topWins[id].num-1].isReverse = isRev;
                    std::push_heap(_pf_topWins[id].list, _pf_topWins[id].list + _pf_topWins[id].num, compareWin);
                }
            }
        }
    }
}
/**********************************************/
// void findTopWins4(uint32_t read_len, SeedList *seeds, int isRev, int readIdx, int id)
// {
//  int i, j;
//  int winNumLimit;
//  float tmpScore;
//  uint32_t lb = (uint32_t) read_len * 0.2;
//  uint32_t ub = (uint32_t) read_len * 1.5;

//  std::sort(seeds->list, seeds->list + seeds->num, compareSeed);

//  // for(i=0; i<seeds->num; i++)
//  // {
//  //  fprintf(stderr, "> %d %u %u %u\n", i, seeds->list[i].tPos, seeds->list[i].qPos, seeds->list[i].len);
//  // }

//  i = 0;
//  j = 0;
//  while(i < seeds->num || j < seeds->num)
//  {
//      // while(j+1 < seeds->num && seeds->list[j+1].tPos - seeds->list[i].tPos < ub)
//      while(j+1 < seeds->num && seeds->list[j+1].tPos - seeds->list[j].tPos < ub)
//      {
//          j++;
//      }
//      if(i == j)
//      {
//          j++;
//      }
//      // else if(seeds->list[j].tPos - seeds->list[i].tPos < ub)
//      else
//      {
//          // do chaining
//          tmpScore = calcChainScore(read_len, seeds->list[i].tPos, seeds->list[j].tPos, isRev, id);
//          if(_pf_topWins[id].num < _pf_maxWin) // the list has some space, push the chains 
//          {
//              _pf_topWins[id].list[_pf_topWins[id].num].tStart = seeds->list[i].tPos;
//              _pf_topWins[id].list[_pf_topWins[id].num].tEnd = seeds->list[j].tPos;
//              _pf_topWins[id].list[_pf_topWins[id].num].score = tmpScore;
//              _pf_topWins[id].list[_pf_topWins[id].num++].isReverse = isRev;
//              std::push_heap(_pf_topWins[id].list, _pf_topWins[id].list + _pf_topWins[id].num, compareWin);
//          }
//          else
//          {
//              if(tmpScore > _pf_topWins[id].list[0].score) // if the score is higher than the current smallest score, replace it with that
//              {
//                  std::pop_heap(_pf_topWins[id].list, _pf_topWins[id].list + _pf_topWins[id].num, compareWin); // the smallest score chain is at the last position
//                  // update the last chain
//                  _pf_topWins[id].list[_pf_topWins[id].num-1].tStart = seeds->list[i].tPos;
//                  _pf_topWins[id].list[_pf_topWins[id].num-1].tEnd = seeds->list[j].tPos;
//                  _pf_topWins[id].list[_pf_topWins[id].num-1].score = tmpScore;
//                  _pf_topWins[id].list[_pf_topWins[id].num-1].isReverse = isRev;
//                  std::push_heap(_pf_topWins[id].list, _pf_topWins[id].list + _pf_topWins[id].num, compareWin);
//              }
//          }
//          // fprintf(stderr, "chaining (%d, %d) (%u, %u) (%d)\n", i, j, seeds->list[i].tPos, seeds->list[j].tPos, seeds->list[j].tPos - seeds->list[i].tPos);
//          // i++;
//          i = j;
//      }
//      // else
//      // {
//      //  i++;
//      // }
//  }
// }

// bool compareChain(const Chain_t& c1, const Chain_t& c2)
// {
//  return c1.score > c2.score;
// }

bool compareSeed(const Seed_t& s1, const Seed_t& s2)
{
    return s1.tPos < s2.tPos;
}

bool compareWin(const Win_t& w1, const Win_t& w2)
{
    return w1.score > w2.score;
}

bool compareSam(const SamList_t& s1, const SamList_t& s2)
{
    // return (s1.samList.size() > 0 && s2.samList.size() == 0) || (s1.totalScore > s2.totalScore);
    return (s1.totalScore > s2.totalScore);
    // return (s1.totalScore > s2.totalScore || 
    //         (s1.totalScore == s2.totalScore && s1.samList.size() > 0 && s2.samList.size() > 0 && s1.samList.front().pos < s2.samList.front().pos));
}

/**********************************************/
void alignWin(Win_t &win, char *query, char *query_rev, uint32_t rLen, char *qual, char* qual_rev, SamList_t &map, int id)
{
    uint32_t i;
    uint32_t margin = rLen >> 1;
    int64_t seedPos_Low = (int64_t)win.tStart - margin;
    int64_t seedPos_High = (int64_t)win.tEnd + margin;

    if(win.isReverse)
    {
        // map.flag = map.flag | 16;
        // map.seq = query_rev;

        _pf_seedsSelected[id].num = 0;
        for(i = 0; i < _pf_seedsReverse[id].num; i++)
        {
            if((int64_t)_pf_seedsReverse[id].list[i].tPos >= seedPos_Low && (int64_t)_pf_seedsReverse[id].list[i].tPos <= seedPos_High)
            {
                _pf_seedsSelected[id].list[_pf_seedsSelected[id].num++] = _pf_seedsReverse[id].list[i];
            }
        }

        LOG2({
            int ilog;
            for(ilog = 0; ilog < _pf_seedsSelected[id].num; ilog++)
            {
                fprintf(stderr, "\tseed %d:\tqPos: %u\ttPos: %u\tlen: %u\t-\n", ilog+1, _pf_seedsSelected[id].list[ilog].qPos, _pf_seedsSelected[id].list[ilog].tPos, _pf_seedsSelected[id].list[ilog].len);
            }
        });

        if(chainAlg == CHAIN_ALG_CLASP)
        {
            if(seedPos_Low > 2000000000)
                for(i = 0; i < _pf_seedsSelected[id].num; i++)
                    _pf_seedsSelected[id].list[i].tPos -= 2000000000;

            chain_seeds_clasp(_pf_seedsSelected[id].list, _pf_seedsSelected[id].num, _pf_topChains[id].list[0]);

            LOG1({
                fprintf(stderr, "\tchain\tread: (%u, %u)\tref: (%u, %u)\t%c\tchainLen: %u\tchainScore: %f\n", 
                    _pf_topChains[id].list[0].seeds[0].qPos, _pf_topChains[id].list[0].seeds[_pf_topChains[id].list[0].chainLen - 1].qPos, 
                    _pf_topChains[id].list[0].seeds[0].tPos, _pf_topChains[id].list[0].seeds[_pf_topChains[id].list[0].chainLen - 1].tPos, 
                    (win.isReverse ? '-' : '+'), _pf_topChains[id].list[0].chainLen, _pf_topChains[id].list[0].score);
            });

            if(seedPos_Low > 2000000000)
                for(i = 0; i < _pf_topChains[id].list[0].chainLen; i++)
                    _pf_topChains[id].list[0].seeds[i].tPos += 2000000000;
        }
        else if(chainAlg == CHAIN_ALG_DPN2)
        {
            chain_seeds_n2(_pf_seedsSelected[id].list, _pf_seedsSelected[id].num, _pf_topChains[id].list[0]);

            LOG1({
                fprintf(stderr, "\tchain\tread: (%u, %u)\tref: (%u, %u)\t%c\tchainLen: %u\tchainScore: %f\n", 
                    _pf_topChains[id].list[0].seeds[0].qPos, _pf_topChains[id].list[0].seeds[_pf_topChains[id].list[0].chainLen - 1].qPos, 
                    _pf_topChains[id].list[0].seeds[0].tPos, _pf_topChains[id].list[0].seeds[_pf_topChains[id].list[0].chainLen - 1].tPos, 
                    (win.isReverse ? '-' : '+'), _pf_topChains[id].list[0].chainLen, _pf_topChains[id].list[0].score);
            });
        }
        
        // _pf_topChains[id].num = 1;

        if(_pf_topChains[id].list[0].chainLen > 1)
        {
            // fprintf(stderr, "\twinCount: %f chainLen: %u \n", win.score, _pf_topChains[id].list[0].chainLen);
            alignChain(_pf_topChains[id].list[0], query_rev, rLen, 1, map);
            map.totalScore = 0;
            for(i = 0; i < map.samList.size(); i++)
            {
                map.samList[i].alnScore = map.samList[i].nmCount + (map.samList[i].qEnd - map.samList[i].qStart);
                map.totalScore += map.samList[i].nmCount;
            }
            for(i = 0; i < map.samList.size() - 1; i++)
            {
                uint32_t diff = abs((int64_t)map.samList[i+1].pos - (int64_t)map.samList[i].posEnd) + abs((int64_t)map.samList[i+1].qStart - (int64_t)map.samList[i].qEnd);
                // fprintf(stderr, "diff: %u\n", diff);
                map.totalScore -= gapPenalty * diff;
                // if(diff > 1) map.totalScore -= (25 * log(diff));
            }
            // penalize clipping
            map.totalScore -= map.samList.front().qStart;
            map.totalScore -= (rLen - map.samList.back().qEnd);
        }
        else
        {
            // strcpy(map.cigar, "*");
            // map.flag = 4;
            // map.seq = query;
            map.totalScore = -2 * rLen;
        }
    }
    else
    {
        // map.seq = query;

        _pf_seedsSelected[id].num = 0;
        for(i = 0; i < _pf_seedsForward[id].num; i++)
        {
            if((int64_t)_pf_seedsForward[id].list[i].tPos >= seedPos_Low && (int64_t)_pf_seedsForward[id].list[i].tPos <= seedPos_High)
            {
                _pf_seedsSelected[id].list[_pf_seedsSelected[id].num++] = _pf_seedsForward[id].list[i];
            }
        }

        LOG2({
            int ilog;
            for(ilog = 0; ilog < _pf_seedsSelected[id].num; ilog++)
            {
                fprintf(stderr, "\tseed %d:\tqPos: %u\ttPos: %u\tlen: %u\t+\n", ilog+1, _pf_seedsSelected[id].list[ilog].qPos, _pf_seedsSelected[id].list[ilog].tPos, _pf_seedsSelected[id].list[ilog].len);
            }
        });

        if(chainAlg == CHAIN_ALG_CLASP)
        {
            if(seedPos_Low > 2000000000)
                for(i = 0; i < _pf_seedsSelected[id].num; i++)
                    _pf_seedsSelected[id].list[i].tPos -= 2000000000;

            chain_seeds_clasp(_pf_seedsSelected[id].list, _pf_seedsSelected[id].num, _pf_topChains[id].list[0]);

            LOG1({
                fprintf(stderr, "\tchain\tread: (%u, %u)\tref: (%u, %u)\t%c\tchainLen: %u\tchainScore: %f\n", 
                    _pf_topChains[id].list[0].seeds[0].qPos, _pf_topChains[id].list[0].seeds[_pf_topChains[id].list[0].chainLen - 1].qPos, 
                    _pf_topChains[id].list[0].seeds[0].tPos, _pf_topChains[id].list[0].seeds[_pf_topChains[id].list[0].chainLen - 1].tPos, 
                    (win.isReverse ? '-' : '+'), _pf_topChains[id].list[0].chainLen, _pf_topChains[id].list[0].score);
            });

            if(seedPos_Low > 2000000000)
                for(i = 0; i < _pf_topChains[id].list[0].chainLen; i++)
                    _pf_topChains[id].list[0].seeds[i].tPos += 2000000000;
        }
        else if(chainAlg == CHAIN_ALG_DPN2)
        {
            chain_seeds_n2(_pf_seedsSelected[id].list, _pf_seedsSelected[id].num, _pf_topChains[id].list[0]);

            LOG1({
                fprintf(stderr, "\tchain\tread: (%u, %u)\tref: (%u, %u)\t%c\tchainLen: %u\tchainScore: %f\n", 
                    _pf_topChains[id].list[0].seeds[0].qPos, _pf_topChains[id].list[0].seeds[_pf_topChains[id].list[0].chainLen - 1].qPos, 
                    _pf_topChains[id].list[0].seeds[0].tPos, _pf_topChains[id].list[0].seeds[_pf_topChains[id].list[0].chainLen - 1].tPos, 
                    (win.isReverse ? '-' : '+'), _pf_topChains[id].list[0].chainLen, _pf_topChains[id].list[0].score);
            });
        }

        // _pf_topChains[id].num = 1;
        
        if(_pf_topChains[id].list[0].chainLen > 1)
        {
            // fprintf(stderr, "\twinCount: %f chainLen: %u \n", win.score, _pf_topChains[id].list[0].chainLen);
            alignChain(_pf_topChains[id].list[0], query, rLen, 0, map);
            map.totalScore = 0;
            for(i = 0; i < map.samList.size(); i++)
            {
                map.samList[i].alnScore = map.samList[i].nmCount + (map.samList[i].qEnd - map.samList[i].qStart);
                map.totalScore += map.samList[i].nmCount;
            }
            for(i = 0; i < map.samList.size() - 1; i++)
            {
                uint32_t diff = abs((int64_t)map.samList[i+1].pos - (int64_t)map.samList[i].posEnd) + abs((int64_t)map.samList[i+1].qStart - (int64_t)map.samList[i].qEnd);
                // fprintf(stderr, "diff: %u\n", diff);
                map.totalScore -= 0.15 * diff;
                // if(diff > 1) map.totalScore -= (25 * log(diff));
            }
            // penalize clipping
            map.totalScore -= map.samList.front().qStart;
            map.totalScore -= (rLen - map.samList.back().qEnd);
        }
        else
        {
            // strcpy(map.cigar, "*");
            // map.flag = 4;
            // map.seq = query;
            map.totalScore = -2 * rLen;
        }
    }

    // fprintf(stderr, "win %u %u\n", win.tStart, win.tEnd);
    // fprintf(stderr, "chain %u %u\n", _pf_topChains[id].list[0].seeds[0].tPos, _pf_topChains[id].list[0].seeds[_pf_topChains[id].list[0].chainLen-1].tPos);
    // bwt_get_chr_boundaries(win.tStart, win.tEnd);

    // fprintf(stderr, "\nSelected seeds\n");
    // for(i = 0; i < _pf_seedsSelected[id].num; i++)
    //  fprintf(stderr, "   %u %u %u\n", _pf_seedsSelected[id].list[i].qPos, _pf_seedsSelected[id].list[i].tPos, _pf_seedsSelected[id].list[i].len);

    // fprintf(stderr, "\nChain seeds\n");
    // for(i = 0; i < _pf_topChains[id].list[0].chainLen; i++)
    //  fprintf(stderr, "   %u %u %u\n", _pf_topChains[id].list[0].seeds[i].qPos, _pf_topChains[id].list[0].seeds[i].tPos, _pf_topChains[id].list[0].seeds[i].len);
}
/**********************************************/
void convertChar2int(uint8_t *seqInt, char *seqChar, int len)
{
    for(int i=0; i<len; i++)
        seqInt[i] = _pf_char2int[seqChar[i]];
}
/**********************************************/
void reverseComplementIntStr(uint8_t *seqDst, uint8_t *seqSrc, int len)
{
    for(int i=0; i<len; i++)
        seqDst[i] = 3 - seqSrc[len - i - 1];
}
/**********************************************/
void fixCigarM(std::string &cigar, std::string semiCigar)
{
    std::istringstream sin(semiCigar);
    std::ostringstream sout;
    int n;
    char c;
    int cntM = 0;
    while(sin >> n >> c)
    {
        if(c == 'M')
        {
            cntM += n;
        }
        else
        {
            if(cntM)
            {
                sout<< cntM << "M";
                cntM = 0;
            }
            sout<< n << c;
        }
    }
    if(cntM)
    {
        sout<< cntM << "M";
    }
    cigar = sout.str();
}
/**********************************************/
void fixCigar(std::string &cigar, std::string semiCigar)
{
    std::istringstream sin(semiCigar);
    std::ostringstream sout;
    int n, cntCh = 0, operationNum = 0;
    char c, ch = 'Z';
    while(sin >> n >> c)
    {
        if(c == ch)
        {
            cntCh += n;
        }
        else
        {
            if(cntCh)
            {
                sout<< cntCh << ( operationNum == 1 && ch == 'I' ? 'S' : ch);
                cntCh = 0;
            }
            cntCh = n;
            ch = c;
            operationNum++;
        }
    }
    if(cntCh) // last CIGAR operation
    {
        sout<< cntCh << ( ch == 'I' ? 'S' : ch);
    }
    cigar = sout.str();
}
/**********************************************/
void alignChain_ksw(Chain_t &chain, char *query, int32_t readLen, int isRev, SamList_t &map)
{
    int i;
    int j;
    // uint32_t readLen = (*read->length);
    // query = (isReverse ? read->rseq : read->seq);

    uint8_t readAlnSeq[SEQ_MAX_LENGTH];
    uint8_t refAlnSeq[SEQ_MAX_LENGTH];
    uint8_t tmpAlnSeq[SEQ_MAX_LENGTH];
    uint32_t readAlnStart, refAlnStart;
    uint32_t readAlnEnd, refAlnEnd;
    int32_t readAlnLen, refAlnLen;
    int tLen, qLen;
    int cigarNum;
    uint32_t *cigar;
    int bandWidth;
    std::ostringstream soutCigar;
    std::string alnCigar;
    int alnScore = 0;
    Sam_t tmpSam;
    // int alnScore_tmp;

    // // print query
    // fprintf(stderr, "query:\t%s\n", query);
    // // print reference
    // refAlnLen = chain.seeds[chain.chainLen-1].tPos - chain.seeds[0].tPos + 200;
    // bwt_str_pac2int(chain.seeds[0].tPos - 100, refAlnLen, refAlnSeq);
    // fprintf(stderr, "  ref:\t");
    // for(j=0; j<refAlnLen; j++)
    //  fprintf(stderr, "%c", "ACGTN"[refAlnSeq[j]]);
    // fprintf(stderr, "\n");

    // set the flag
    tmpSam.flag = (isRev ? 16 : 0);
    // set sam pos
    tmpSam.pos = chain.seeds[0].tPos;

    // extend before first seed
    readAlnLen = chain.seeds[0].qPos;

    // fprintf(stderr, "readLen: %d \t before: %d\n", readLen, readAlnLen);

    if(readAlnLen > 0)
    {
        convertChar2int(tmpAlnSeq, query, readAlnLen);
        reverseComplementIntStr(readAlnSeq, tmpAlnSeq, readAlnLen);

        refAlnStart = chain.seeds[0].tPos - readAlnLen;
        refAlnLen = readAlnLen;
        bwt_str_pac2int(refAlnStart, refAlnLen, tmpAlnSeq);
        reverseComplementIntStr(refAlnSeq, tmpAlnSeq, refAlnLen);

        ksw_extend(readAlnLen, readAlnSeq, refAlnLen, refAlnSeq, 5, _pf_kswMatrix, _pf_kswGapOpen, _pf_kswGapExtend, 40, 0, 40, readAlnLen, &qLen, &tLen, 0, 0, 0);
        bandWidth = (qLen > tLen ? qLen : tLen);
        alnScore += ksw_global(qLen, readAlnSeq, tLen, refAlnSeq, 5, _pf_kswMatrix, _pf_kswGapOpen, _pf_kswGapExtend, bandWidth, &cigarNum, &cigar);

        if(qLen < readAlnLen)
        {
            soutCigar << (readAlnLen - qLen) << "S";
        }
        for (int z = cigarNum - 1; z >= 0; --z)
        {
            soutCigar << (cigar[z]>>4) << _pf_kswCigarTable[cigar[z]&0xf];
        }
        free(cigar);

        // fix sam pos
        tmpSam.pos = chain.seeds[0].tPos - tLen;
    }

    for(i=0; i<chain.chainLen-1; i++)
    {
        soutCigar << chain.seeds[i].len << "M";
        alnScore += chain.seeds[i].len * _pf_kswMatch;

        readAlnStart = chain.seeds[i].qPos + chain.seeds[i].len;
        refAlnStart = chain.seeds[i].tPos + chain.seeds[i].len;
        readAlnEnd = chain.seeds[i+1].qPos; // acutally chain.seeds[i+1].qPos - 1
        refAlnEnd = chain.seeds[i+1].tPos; // actually chain.seeds[i+1].tPos - 1
        readAlnLen = readAlnEnd - readAlnStart;
        refAlnLen = refAlnEnd - refAlnStart;

        if(readAlnLen > 0 && refAlnLen > 0)
        {
            // fprintf(stderr, "\tseed\t%u\t%u\t%u\t%u\t%u\n", chain.seeds[i].qPos, chain.seeds[i].qPos + chain.seeds[i].len - 1, 
            //  chain.seeds[i].tPos, chain.seeds[i].tPos + chain.seeds[i].len - 1, chain.seeds[i].len);
            // fprintf(stderr, "\tseed\tquery\t%.*s\n", chain.seeds[i].len, query + chain.seeds[i].qPos);
            // // print reference
            // bwt_str_pac2int(chain.seeds[i].tPos, chain.seeds[i].len, refAlnSeq);
            // fprintf(stderr, "\tseed\t  ref\t");
            // for(j=0; j<chain.seeds[i].len; j++)
            //  fprintf(stderr, "%c", "ACGTN"[refAlnSeq[j]]);
            // fprintf(stderr, "\n");

            convertChar2int(readAlnSeq, query+readAlnStart, readAlnLen);
            bwt_str_pac2int(refAlnStart, refAlnLen, refAlnSeq);

            // fprintf(stderr, "\t\t+++ aln\t(%u, %u, %u)\t(%u, %u, %u)\n", readAlnStart, readAlnEnd-1, readAlnLen, refAlnStart, refAlnEnd-1, refAlnLen);
            // // print reference
            // fprintf(stderr, "\t\t+++ ref\t");
            // for(j=0; j<refAlnLen; j++)
            //  fprintf(stderr, "%c", "ACGTN"[refAlnSeq[j]]);
            // fprintf(stderr, "\n");
            // // print reference
            // fprintf(stderr, "\t\t+++ query\t");
            // for(j=0; j<readAlnLen; j++)
            //  fprintf(stderr, "%c", "ACGTN"[readAlnSeq[j]]);
            // fprintf(stderr, "\n");

            bandWidth = (readAlnLen > refAlnLen ? readAlnLen : refAlnLen);
            alnScore += (/*alnScore_tmp =*/ ksw_global(readAlnLen, readAlnSeq, refAlnLen, refAlnSeq, 5, _pf_kswMatrix, _pf_kswGapOpen, _pf_kswGapExtend, bandWidth, &cigarNum, &cigar) );

            // fprintf(stderr, "\t\t+++ score\t%d\n", alnScore_tmp);
            // fprintf(stderr, "\t\t+++ cigar: ");
            for (int z=0; z<cigarNum; ++z)
            {
            //  fprintf(stderr, "%u%c", cigar[z]>>4, _pf_kswCigarTable[cigar[z]&0xf]);
                soutCigar << (cigar[z]>>4) << _pf_kswCigarTable[cigar[z]&0xf];
            }
            // fprintf(stderr, "\n");
            // fprintf(stderr, "\t\t+++   all: %s\n", soutCigar.str().c_str());

            free(cigar);
        }
        else // no need to do smith-waterman, either insertion or deletion!
        {
            if(readAlnLen > 0)
            {
                soutCigar << readAlnLen << "I";
                // alnScore -= _pf_kswGapOpen + (readAlnLen - 1) * _pf_kswGapExtend;
                alnScore -= _pf_kswGapOpen + readAlnLen * _pf_kswGapExtend;
            }
            else
            {
                soutCigar << refAlnLen << "D";
                // alnScore -= _pf_kswGapOpen + (refAlnLen - 1) * _pf_kswGapExtend;
                alnScore -= _pf_kswGapOpen + refAlnLen * _pf_kswGapExtend;
            }
        }
    }
    // // last seed
    // fprintf(stderr, "\tseed\t%u\t%u\t%u\t%u\t%u\n", chain.seeds[i].qPos, chain.seeds[i].qPos + chain.seeds[i].len - 1, 
    //  chain.seeds[i].tPos, chain.seeds[i].tPos + chain.seeds[i].len - 1, chain.seeds[i].len);
    // fprintf(stderr, "\tseed\tquery\t%.*s\n", chain.seeds[i].len, query + chain.seeds[i].qPos);
    // // print reference
    // bwt_str_pac2int(chain.seeds[i].tPos, chain.seeds[i].len, refAlnSeq);
    // fprintf(stderr, "\tseed\t  ref\t");
    // for(j=0; j<chain.seeds[i].len; j++)
    //  fprintf(stderr, "%c", "ACGTN"[refAlnSeq[j]]);
    // fprintf(stderr, "\n");

    soutCigar << chain.seeds[i].len << "M";
    alnScore += chain.seeds[i].len * _pf_kswMatch;
    tmpSam.posEnd = chain.seeds[i].tPos + chain.seeds[i].len - 1;

    // extend after last seed
    readAlnStart = chain.seeds[i].qPos + chain.seeds[i].len;
    readAlnLen = readLen - readAlnStart;

    // fprintf(stderr, "readLen: %d \t after: %d\n", readLen, readAlnLen);

    if(readAlnLen > 0)
    {
        refAlnStart = chain.seeds[i].tPos + chain.seeds[i].len;
        refAlnLen = readAlnLen;

        convertChar2int(readAlnSeq, query+readAlnStart, readAlnLen);
        bwt_str_pac2int(refAlnStart, refAlnLen, refAlnSeq);

        ksw_extend(readAlnLen, readAlnSeq, refAlnLen, refAlnSeq, 5, _pf_kswMatrix, _pf_kswGapOpen, _pf_kswGapExtend, 40, 0, 40, readAlnLen, &qLen, &tLen, 0, 0, 0);
        bandWidth = (qLen > tLen ? qLen : tLen);
        alnScore += ksw_global(qLen, readAlnSeq, tLen, refAlnSeq, 5, _pf_kswMatrix, _pf_kswGapOpen, _pf_kswGapExtend, bandWidth, &cigarNum, &cigar);

        for (int z=0; z<cigarNum; ++z)
        {
            soutCigar << (cigar[z]>>4) << _pf_kswCigarTable[cigar[z]&0xf];
        }
        if(qLen < readAlnLen)
        {
            soutCigar << (readAlnLen - qLen) << "S";
        }
        free(cigar);

        tmpSam.posEnd = refAlnStart + tLen - 1;
    }

    // fprintf(stderr, "score: %d\n", alnScore);
    // fprintf(stderr, "cigar: %s\n", soutCigar.str().c_str());
    fixCigar(alnCigar, soutCigar.str());
    // fprintf(stderr, "cigar: %s\n", alnCigar.c_str());
    // map.cigar = soutCigar.str();
    tmpSam.cigar = alnCigar;
    // map.cigar = (char*) malloc((alnCigar.size() + 1) * sizeof(char));
    // strcpy(map.cigar, alnCigar.c_str());
    tmpSam.alnScore = alnScore;
    // return alnScore;

    map.totalScore = alnScore;
    map.samList.push_back(tmpSam);
}
/**********************************************/
void edlibGetCigar(const unsigned char* const alignment, const int alignmentLength, const EdlibCigarFormat cigarFormat, std::deque<char> *cigar)
{
    if (cigarFormat != EDLIB_CIGAR_EXTENDED && cigarFormat != EDLIB_CIGAR_STANDARD)
    {
        return;
    }

    // Maps move code from alignment to char in cigar.
    //                        0    1    2    3
    char moveCodeToChar[] = {'=', 'I', 'D', 'X'};
    if (cigarFormat == EDLIB_CIGAR_STANDARD)
    {
        moveCodeToChar[0] = moveCodeToChar[3] = 'M';
    }

    char lastMove = 0;  // Char of last move. 0 if there was no previous move.
    int numOfSameMoves = 0;
    for (int i = 0; i <= alignmentLength; i++)
    {
        // if new sequence of same moves started
        if (i == alignmentLength || (moveCodeToChar[alignment[i]] != lastMove && lastMove != 0))
        {
            // Write number of moves to cigar string.
            int numDigits = 0;
            for (; numOfSameMoves; numOfSameMoves /= 10)
            {
                cigar->push_back('0' + numOfSameMoves % 10);
                numDigits++;
            }
            reverse(cigar->end() - numDigits, cigar->end());
            // Write code of move to cigar string.
            cigar->push_back(lastMove);
            // If not at the end, start new sequence of moves.
            if (i < alignmentLength)
            {
                // TODO: is this check necessary?!
                // Check if alignment has valid values.
                if (alignment[i] > 3)
                {
                    return;
                }
                numOfSameMoves = 0;
            }
        }
        if (i < alignmentLength)
        {
            lastMove = moveCodeToChar[alignment[i]];
            numOfSameMoves++;
        }
    }
}
/**********************************************/
void edlibGetCigarReverse(const unsigned char* const alignment, const int alignmentLength, const EdlibCigarFormat cigarFormat, std::deque<char> *cigar)
{
    if (cigarFormat != EDLIB_CIGAR_EXTENDED && cigarFormat != EDLIB_CIGAR_STANDARD)
    {
        return;
    }

    // Maps move code from alignment to char in cigar.
    //                        0    1    2    3
    char moveCodeToChar[] = {'=', 'I', 'D', 'X'};
    if (cigarFormat == EDLIB_CIGAR_STANDARD)
    {
        moveCodeToChar[0] = moveCodeToChar[3] = 'M';
    }

    char lastMove = 0;  // Char of last move. 0 if there was no previous move.
    int numOfSameMoves = 0;
    for (int i = 0; i <= alignmentLength; i++)
    {
        // if new sequence of same moves started
        if (i == alignmentLength || (moveCodeToChar[alignment[i]] != lastMove && lastMove != 0))
        {
            // Write code of move to cigar string.
            cigar->push_front(lastMove);
            // Write number of moves to cigar string.
            // int numDigits = 0;
            for (; numOfSameMoves; numOfSameMoves /= 10)
            {
                cigar->push_front('0' + numOfSameMoves % 10);
                // numDigits++;
            }
            // reverse(cigar->end() - numDigits, cigar->end());
            // If not at the end, start new sequence of moves.
            if (i < alignmentLength)
            {
                // TODO: is this check necessary?!
                // Check if alignment has valid values.
                if (alignment[i] > 3)
                {
                    return;
                }
                numOfSameMoves = 0;
            }
        }
        if (i < alignmentLength)
        {
            lastMove = moveCodeToChar[alignment[i]];
            numOfSameMoves++;
        }
    }
}
/**********************************************/
inline void edlibCigar_pushback(EdlibAlignResult &edResult, std::deque<char> *cigar)
{
    // Map integer encoded alignments to char in cigar.
    //                  0    1    2    3
    char aln2char[] = {'M', 'I', 'D', 'M'};
    // char aln2char[] = {'=', 'I', 'D', 'X'};

    for(int i = 0; i < edResult.alignmentLength; i++)
    {
        cigar->push_back(aln2char[edResult.alignment[i]]);
    }
}
/**********************************************/
inline void edlibCigar_pushfront(EdlibAlignResult &edResult, std::deque<char> *cigar)
{
    // Map integer encoded alignments to char in cigar.
    //                  0    1    2    3
    char aln2char[] = {'M', 'I', 'D', 'M'};
    // char aln2char[] = {'=', 'I', 'D', 'X'};

    for(int i = 0; i < edResult.alignmentLength; i++)
    {
        cigar->push_front(aln2char[edResult.alignment[i]]);
    }
}
/**********************************************/
inline std::string edlibCigar_toString(std::deque<char> *cigar)
{
    std::ostringstream sout;
    char ch = 0;
    int num = 0;
    int operationNum = 0;
    for(int i = 0; i < cigar->size(); i++)
    {
        if(cigar->at(i) != ch)
        {
            if(ch != 0)
            {
                sout<< num << (operationNum == 0 && ch == 'I' ? 'S' : ch);
                operationNum++;
            }
            num = 1;
            ch = cigar->at(i);
        }
        else
        {
            num++;
        }
    }
    // 
    if(num)
    {
        sout<< num << (ch == 'I' ? 'S' : ch);
    }
    // 
    return sout.str();
}
/**********************************************/
inline void edlibMD_pushback(const char *query, const char *target, EdlibAlignResult &edResult, std::deque<char> *md)
{
    if(edResult.alignmentLength <= 0)
    {
        fprintf(stderr, "[edlibMD_pushback] invalid alignment length!\n");
        exit(1);
    }

    int qIndex = 0;
    int tIndex = 0;
    for(int i = 0; i < edResult.alignmentLength; i++)
    {
        switch(edResult.alignment[i])
        {
        case 0:
            md->push_back('=');
            qIndex++;
            tIndex++;
            break;
        case 1:
            md->push_back('-');
            qIndex++;
            break;
        case 2:
            md->push_back(target[tIndex]);
            tIndex++;
            break;
        case 3:
            md->push_back(target[tIndex]);
            tIndex++;
            qIndex++;
            break;
        default:
            fprintf(stderr, "[edlibMD_pushback] invalid alignment value!\n");
            exit(1);
        }
    }
}
/**********************************************/
inline void edlibMD_pushfront(const char *query, const char *target, EdlibAlignResult &edResult, std::deque<char> *md)
{
    if(edResult.alignmentLength <= 0)
    {
        fprintf(stderr, "[edlibMD_pushback] invalid alignment length!\n");
        exit(1);
    }

    int qIndex = 0;
    int tIndex = 0;
    char tableComplement[128] = {
    'N','N','N','N',   'N','N','N','N',   'N','N','N','N',   'N','N','N','N', 
    'N','N','N','N',   'N','N','N','N',   'N','N','N','N',   'N','N','N','N', 
    'N','N','N','N',   'N','N','N','N',   'N','N','N','N',   'N','N','N','N', 
    'N','N','N','N',   'N','N','N','N',   'N','N','N','N',   'N','N','N','N', 
    'N','T','N','G',   'N','N','N','C',   'N','N','N','N',   'N','N','N','N', 
    'N','N','N','N',   'A','N','N','N',   'N','N','N','N',   'N','N','N','N', 
    'N','T','N','G',   'N','N','N','C',   'N','N','N','N',   'N','N','N','N', 
    'N','N','N','N',   'A','N','N','N',   'N','N','N','N',   'N','N','N','N'
    };

    for(int i = 0; i < edResult.alignmentLength; i++)
    {
        switch(edResult.alignment[i])
        {
        case 0:
            md->push_front('=');
            qIndex++;
            tIndex++;
            break;
        case 1:
            md->push_front('-');
            qIndex++;
            break;
        case 2:
            md->push_front(tableComplement[target[tIndex]]);
            tIndex++;
            break;
        case 3:
            md->push_front(tableComplement[target[tIndex]]);
            tIndex++;
            qIndex++;
            break;
        default:
            fprintf(stderr, "[edlibMD_pushback] invalid alignment value!\n");
            exit(1);
        }
    }
}
/**********************************************/
inline std::string edlibMD_toString(std::deque<char> *md, std::deque<char> *cigar)
{
    std::ostringstream sout;
    int num = 0;
    char ch_md;
    char ch_cigar;
    char last_move = '='; // match(=), mis-match(X), insertion(I), deletion(D)
    for(int i = 0; i < md->size(); i++)
    {
        ch_md = md->at(i);
        ch_cigar = cigar->at(i);
        if(ch_md == '=')
        {
            num++;
            last_move = '=';
        }
        else if(ch_md == '-')
        {
            // nothing
            last_move = 'I';
        }
        else if(ch_cigar == 'M')
        {
            // if(last_move != 'X')
            {
                sout<< num;
                num = 0;
            }
            sout<< ch_md;
            last_move = 'X';
        }
        else if(ch_cigar == 'D')
        {
            if(last_move != 'D')
            {
                sout<< num;
                num = 0;
                sout<< '^';
            }
            sout<< ch_md;
            last_move = 'D';
        }
    }
    // 
    sout<< num;
    return sout.str();
}
/**********************************************/
void alignChain_edlib(Chain_t &chain, char *query, int32_t readLen, int isRev, SamList_t &map)
{
    int i;
    int j;
    // for edlib
    char readAlnSeq[SEQ_MAX_LENGTH];
    char readAlnSeq_rev[SEQ_MAX_LENGTH];
    char refAlnSeq[SEQ_MAX_LENGTH];
    char refAlnSeq_rev[SEQ_MAX_LENGTH];
    uint32_t readAlnStart, refAlnStart;
    uint32_t readAlnEnd, refAlnEnd;
    int32_t readAlnLen, refAlnLen;
    // for ksw
    uint8_t readAlnSeq_ksw[SEQ_MAX_LENGTH];
    uint8_t refAlnSeq_ksw[SEQ_MAX_LENGTH];
    uint8_t readAlnSeq_ksw_rev[SEQ_MAX_LENGTH];
    uint8_t refAlnSeq_ksw_rev[SEQ_MAX_LENGTH];
    int tLen_ksw, qLen_ksw;
    // int tmp_cigarNum;
    // uint32_t *tmp_cigar;
    // int tmp_bandWidth;
    // int tmp_alnScore = 0;
    // uint32_t readAlnStart_new, refAlnStart_new;
    // uint32_t readAlnEnd_new, refAlnEnd_new;
    // int32_t readAlnLen_new, refAlnLen_new;

    int32_t editScore = 0;
    // int32_t editScore_clip = 0;
    EdlibAlignResult edResult, edResult_rev;
    std::deque<char> *edCigar = new std::deque<char>();
    std::deque<char> *edMD = new std::deque<char>();
    Sam_t tmpSam;

    uint32_t chrBeg, chrEnd;
    bwt_get_chr_boundaries(chain.seeds[0].tPos, chain.seeds[chain.chainLen-1].tPos, &chrBeg, &chrEnd);

    // fprintf(stderr, "chain %u %u\n", chain.seeds[0].tPos, chain.seeds[chain.chainLen-1].tPos);
    // fprintf(stderr, "chr %u %u\n", chrBeg, chrEnd);

    // // print query
    // fprintf(stderr, "query:\t%s\n", query);
    // // print reference
    // refAlnLen = chain.seeds[chain.chainLen-1].tPos - chain.seeds[0].tPos + 200;
    // bwt_str_pac2int(chain.seeds[0].tPos - 100, refAlnLen, refAlnSeq);
    // fprintf(stderr, "  ref:\t");
    // for(j=0; j<refAlnLen; j++)
    //  fprintf(stderr, "%c", "ACGTN"[refAlnSeq[j]]);
    // fprintf(stderr, "\n");

    // set the flag
    tmpSam.flag = (isRev ? 16 : 0);
    // set sam pos
    tmpSam.pos = chain.seeds[0].tPos;
    tmpSam.qStart = chain.seeds[0].qPos;

    // extend before first seed
    readAlnLen = chain.seeds[0].qPos;
    refAlnLen = readAlnLen + 20; // TODO: how much extra sequence?
    if(readAlnLen > 0)
    {
        if((int64_t)chain.seeds[0].tPos - refAlnLen >= (int64_t)chrBeg)
        {
            reverseComplement(query, readAlnSeq, readAlnLen);

            refAlnStart = chain.seeds[0].tPos - refAlnLen;
            bwt_str_pac2char(refAlnStart, refAlnLen, refAlnSeq);
            reverseComplement(refAlnSeq, refAlnSeq_rev, refAlnLen);

            edResult = edlibAlign(readAlnSeq, readAlnLen, refAlnSeq_rev, refAlnLen, edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH));

            LOG2({
                fprintf(stderr, "\tbeg\tqLen: %d\talnEdit: %d\talnSim: %f\n", 
                    readAlnLen, edResult.editDistance, (1 - ((float)edResult.editDistance / readAlnLen)));
            });

            if(readAlnLen > _pf_clipLen && (1 - ((float)edResult.editDistance / readAlnLen)) < _pf_clipSim)
            {
                convertChar2int(readAlnSeq_ksw_rev, query, readAlnLen);
                reverseComplementIntStr(readAlnSeq_ksw, readAlnSeq_ksw_rev, readAlnLen);

                bwt_str_pac2int(refAlnStart, refAlnLen, refAlnSeq_ksw_rev);
                reverseComplementIntStr(refAlnSeq_ksw, refAlnSeq_ksw_rev, refAlnLen);

                ksw_extend(readAlnLen, readAlnSeq_ksw, refAlnLen, refAlnSeq_ksw, 5, _pf_kswMatrix_clip, _pf_kswGapOpen_clip, _pf_kswGapExtend_clip, 40, 0, 40, readAlnLen, &qLen_ksw, &tLen_ksw, 0, 0, 0);
                
                if(qLen_ksw > 0 && qLen_ksw < readAlnLen) // perform a new edlib alignment with new coordinates
                {
                    edlibFreeAlignResult(edResult); // free old results
                    edResult = edlibAlign(readAlnSeq, qLen_ksw, refAlnSeq_rev, tLen_ksw, edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH));
                    // edlibGetCigarReverse(edResult.alignment, edResult.alignmentLength, EDLIB_CIGAR_STANDARD, edCigar);
                    edlibCigar_pushfront(edResult, edCigar);
                    edlibMD_pushfront(readAlnSeq, refAlnSeq_rev, edResult, edMD);
                    // update score
                    editScore -= edResult.editDistance;
                    // fix sam pos
                    tmpSam.pos = chain.seeds[0].tPos - edResult.endLocations[0] - 1;
                    tmpSam.qStart = chain.seeds[0].qPos - qLen_ksw;

                    edCigar->insert(edCigar->begin(), readAlnLen - qLen_ksw, 'I');
                    edMD->insert(edMD->begin(), readAlnLen - qLen_ksw, '-');
                    // update score
                    // editScore_clip -= (readAlnLen - qLen_ksw);
                }
                else // use the already calculated results
                {
                    editScore -= edResult.editDistance;
                    // edlibGetCigarReverse(edResult.alignment, edResult.alignmentLength, EDLIB_CIGAR_STANDARD, edCigar);
                    edlibCigar_pushfront(edResult, edCigar);
                    edlibMD_pushfront(readAlnSeq, refAlnSeq_rev, edResult, edMD);
                    // fix sam pos
                    tmpSam.pos = chain.seeds[0].tPos - edResult.endLocations[0] - 1;
                    tmpSam.qStart = 0;
                }
            }
            else
            {
                editScore -= edResult.editDistance;
                // edlibGetCigarReverse(edResult.alignment, edResult.alignmentLength, EDLIB_CIGAR_STANDARD, edCigar);
                edlibCigar_pushfront(edResult, edCigar);
                edlibMD_pushfront(readAlnSeq, refAlnSeq_rev, edResult, edMD);
                // fix sam pos
                tmpSam.pos = chain.seeds[0].tPos - edResult.endLocations[0] - 1;
                tmpSam.qStart = 0;
            }

            edlibFreeAlignResult(edResult);
        }
        else // not enough sequence left on the chromosome to align => soft-clip
        {
            edCigar->insert(edCigar->begin(), readAlnLen, 'I');
            edMD->insert(edMD->begin(), readAlnLen, '-');
            // update score
            // editScore_clip -= readAlnLen;
        }
    }

    int numAnchorsSoFar = 1;
    for(i=0; i<chain.chainLen-1; i++)
    {
        edCigar->insert(edCigar->end(), chain.seeds[i].len, 'M');
        // edCigar->insert(edCigar->begin(), chain.seeds[i].len, '=');
        edMD->insert(edMD->end(), chain.seeds[i].len, '=');

        readAlnStart = chain.seeds[i].qPos + chain.seeds[i].len;
        refAlnStart = chain.seeds[i].tPos + chain.seeds[i].len;
        readAlnEnd = chain.seeds[i+1].qPos; // acutally chain.seeds[i+1].qPos - 1
        refAlnEnd = chain.seeds[i+1].tPos; // actually chain.seeds[i+1].tPos - 1
        readAlnLen = readAlnEnd - readAlnStart;
        refAlnLen = refAlnEnd - refAlnStart;

        if(readAlnLen > 0 && refAlnLen > 0)
        {
            // fprintf(stderr, "\tseed\t%u\t%u\t%u\t%u\t%u\n", chain.seeds[i].qPos, chain.seeds[i].qPos + chain.seeds[i].len - 1, 
            //  chain.seeds[i].tPos, chain.seeds[i].tPos + chain.seeds[i].len - 1, chain.seeds[i].len);
            // fprintf(stderr, "\tseed\tquery\t%.*s\n", chain.seeds[i].len, query + chain.seeds[i].qPos);
            // // print reference
            // bwt_str_pac2int(chain.seeds[i].tPos, chain.seeds[i].len, refAlnSeq);
            // fprintf(stderr, "\tseed\t  ref\t");
            // for(j=0; j<chain.seeds[i].len; j++)
            //  fprintf(stderr, "%c", "ACGTN"[refAlnSeq[j]]);
            // fprintf(stderr, "\n");

            bwt_str_pac2char(refAlnStart, refAlnLen, refAlnSeq);

            // fprintf(stderr, "\t\t+++ aln\t(%u, %u, %u)\t(%u, %u, %u)\n", readAlnStart, readAlnEnd-1, readAlnLen, refAlnStart, refAlnEnd-1, refAlnLen);
            // // print reference
            // fprintf(stderr, "\t\t+++ ref\t");
            // for(j=0; j<refAlnLen; j++)
            //  fprintf(stderr, "%c", "ACGTN"[refAlnSeq[j]]);
            // fprintf(stderr, "\n");
            // // print reference
            // fprintf(stderr, "\t\t+++ query\t");
            // for(j=0; j<readAlnLen; j++)
            //  fprintf(stderr, "%c", "ACGTN"[readAlnSeq[j]]);
            // fprintf(stderr, "\n");

            edResult = edlibAlign(query+readAlnStart, readAlnLen, refAlnSeq, refAlnLen, edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH));

            LOG2({

                fprintf(stderr, "\tDP: (%d, %d)\tread: (%u, %u)\treadLen: %u\tref: (%u, %u)\trefLen: %u\tedit: %d\tsim: %f\n",
                    i, i+1, readAlnStart, readAlnEnd, readAlnLen, 
                    refAlnStart, refAlnEnd, refAlnLen, 
                    edResult.editDistance, 
                    1 - ((float)edResult.editDistance / readAlnLen));
            });

            if(abs(readAlnLen - refAlnLen) >= _pf_splitLen && (1 - ((float)edResult.editDistance / readAlnLen)) < _pf_splitSim)
            {
            //  uint8_t readAlnSeq_ksw[SEQ_MAX_LENGTH];
            //  uint8_t refAlnSeq_ksw[SEQ_MAX_LENGTH];
            //  uint8_t readAlnSeq_ksw_rev[SEQ_MAX_LENGTH];
            //  uint8_t refAlnSeq_ksw_rev[SEQ_MAX_LENGTH];
            //  int tLen_ksw, qLen_ksw;
            //  int tmp_cigarNum;
            //  uint32_t *tmp_cigar;
            //  int tmp_bandWidth;
            //  int tmp_alnScore = 0;
                uint32_t readAlnStart_new, refAlnStart_new;
                uint32_t readAlnEnd_new, refAlnEnd_new;
                int32_t readAlnLen_new, refAlnLen_new;

                // find the starting coordinate of the potential split
                convertChar2int(readAlnSeq_ksw, query+readAlnStart, readAlnLen);
                bwt_str_pac2int(refAlnStart, refAlnLen, refAlnSeq_ksw);
                // ksw_extend(readAlnLen, readAlnSeq_ksw, refAlnLen, refAlnSeq_ksw, 5, _pf_kswMatrix_clip, _pf_kswGapOpen_clip, _pf_kswGapExtend_clip, 40, 0, 1, readAlnLen, &qLen_ksw, &tLen_ksw, 0, 0, 0);
                ksw_extend2(readAlnLen, readAlnSeq_ksw, refAlnLen, refAlnSeq_ksw, 5, _pf_kswMatrix_clip, 8, 1, 4, 1, 100, 0, 200, readAlnLen, &qLen_ksw, &tLen_ksw, 0, 0, 0);
                readAlnStart_new = readAlnStart + qLen_ksw;
                refAlnStart_new = refAlnStart + tLen_ksw;

                // find the ending coordinate of the potential split
                convertChar2int(readAlnSeq_ksw_rev, query+readAlnStart, readAlnLen);
                reverseComplementIntStr(readAlnSeq_ksw, readAlnSeq_ksw_rev, readAlnLen);
                bwt_str_pac2int(refAlnStart, refAlnLen, refAlnSeq_ksw_rev);
                reverseComplementIntStr(refAlnSeq_ksw, refAlnSeq_ksw_rev, refAlnLen);
                // ksw_extend(readAlnLen, readAlnSeq_ksw, refAlnLen, refAlnSeq_ksw, 5, _pf_kswMatrix_clip, _pf_kswGapOpen_clip, _pf_kswGapExtend_clip, 40, 0, 1, readAlnLen, &qLen_ksw, &tLen_ksw, 0, 0, 0);
                ksw_extend2(readAlnLen, readAlnSeq_ksw, refAlnLen, refAlnSeq_ksw, 5, _pf_kswMatrix_clip, 8, 1, 4, 1, 100, 0, 200, readAlnLen, &qLen_ksw, &tLen_ksw, 0, 0, 0);
                readAlnEnd_new = readAlnEnd - qLen_ksw;
                refAlnEnd_new = refAlnEnd - tLen_ksw;

                refAlnLen_new = refAlnEnd_new - refAlnStart_new;
                readAlnLen_new = readAlnEnd_new - readAlnStart_new;

                LOG2({
                    fprintf(stderr, "\t\textension; read: (%u, %u)\tref: (%u, %u)\n",
                        readAlnStart_new, readAlnEnd_new,
                        refAlnStart_new, refAlnEnd_new);
                });

                // split alignment if extensions do not cross each other
                if(readAlnStart_new < readAlnEnd_new || refAlnStart_new < refAlnEnd_new)
                {
                    edlibFreeAlignResult(edResult);
                    ////////////////////// make the first part of the split alignment
                    if(readAlnStart_new > readAlnStart || refAlnStart_new > refAlnStart)
                    {
                        edResult = edlibAlign(query+readAlnStart, (readAlnStart_new - readAlnStart), refAlnSeq, (refAlnStart_new - refAlnStart), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH));
                        // edlibGetCigar(edResult.alignment, edResult.alignmentLength, EDLIB_CIGAR_STANDARD, edCigar);
                        edlibCigar_pushback(edResult, edCigar);
                        edlibMD_pushback(query+readAlnStart, refAlnSeq, edResult, edMD);
                        editScore -= edResult.editDistance;
                        edlibFreeAlignResult(edResult);
                    }
                    edCigar->insert(edCigar->end(), readLen - readAlnStart_new, 'I');
                    edMD->insert(edMD->end(), readLen - readAlnStart_new, '-');
                    // update score
                    // editScore_clip -= (readLen - readAlnStart_new);
                    tmpSam.cigar = edlibCigar_toString(edCigar);
                    tmpSam.md = edlibMD_toString(edMD, edCigar);
                    tmpSam.posEnd = refAlnStart_new;
                    tmpSam.qEnd = readAlnStart_new;
                    // tmpSam.alnScore = editScore + editScore_clip + readLen;
                    // tmpSam.alnScore = editScore;
                    tmpSam.nmCount = editScore;
                    // tmpSam.clipCount = editScore_clip;
                    // push the split
                    if(numAnchorsSoFar > 1)
                    {
                        map.samList.push_back(tmpSam);
                        // map.totalScore += tmpSam.alnScore - editScore_split;
                        // map.totalScore += tmpSam.alnScore;
                    }
                    // reset
                    edCigar->clear();
                    edMD->clear();
                    editScore = 0;
                    // editScore_clip = 0;

                    ////////////////////// check the middle part of the split (if the reverse complement of the middle part aligns well)
                    if(readAlnStart_new < readAlnEnd_new && refAlnStart_new < refAlnEnd_new)
                    {
                        bwt_str_pac2char(refAlnStart_new, refAlnLen_new, refAlnSeq);
                        edResult = edlibAlign(query+readAlnStart_new, readAlnLen_new, refAlnSeq, refAlnLen_new, edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH));
                        reverseComplement(query+readAlnStart_new, readAlnSeq_rev, readAlnLen_new);
                        edResult_rev = edlibAlign(readAlnSeq_rev, readAlnLen_new, refAlnSeq, refAlnLen_new, edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH));
                        if((1 - ((double)edResult_rev.editDistance / readAlnLen_new)) > (1 - ((double)edResult.editDistance / readAlnLen_new))
                            && (1 - ((double)edResult_rev.editDistance / readAlnLen_new)) > _pf_reverseSim)
                        {
                            tmpSam.flag = (isRev ? 0 : 16); // Note that is opposite of the real direction
                            tmpSam.pos = refAlnStart_new;
                            tmpSam.qStart = readAlnStart_new;
                            tmpSam.posEnd = refAlnEnd_new;
                            tmpSam.qEnd = readAlnEnd_new;
                            edCigar->insert(edCigar->end(), readAlnStart_new, 'I');
                            edMD->insert(edMD->end(), readAlnStart_new, '-');
                            // editScore_clip -= readAlnStart_new;
                            //
                            // edlibGetCigar(edResult_rev.alignment, edResult_rev.alignmentLength, EDLIB_CIGAR_STANDARD, edCigar);
                            edlibCigar_pushback(edResult_rev, edCigar);
                            edlibMD_pushback(readAlnSeq_rev, refAlnSeq, edResult_rev, edMD);
                            editScore -= edResult_rev.editDistance;
                            edCigar->insert(edCigar->end(), readLen - readAlnEnd_new, 'I');
                            edMD->insert(edMD->begin(), readLen - readAlnEnd_new, '-');
                            // editScore_clip -= (readLen - readAlnEnd_new);
                            tmpSam.cigar = edlibCigar_toString(edCigar);
                            tmpSam.md = edlibMD_toString(edMD, edCigar);
                            // tmpSam.alnScore = editScore + editScore_clip + readLen;
                            // tmpSam.alnScore = editScore;
                            tmpSam.nmCount = editScore;
                            // tmpSam.clipCount = editScore_clip;
                            // push the split
                            map.samList.push_back(tmpSam);
                            // map.totalScore += tmpSam.alnScore - editScore_split;
                            // map.totalScore += tmpSam.alnScore;
                            // reset
                            edCigar->clear();
                            edMD->clear();
                            editScore = 0;
                            // editScore_clip = 0;
                        }
                        edlibFreeAlignResult(edResult);
                        edlibFreeAlignResult(edResult_rev);
                    }
                    
                    ////////////////////// make the second part of the split alignment
                    if(readAlnEnd_new < readAlnEnd || refAlnEnd_new < refAlnEnd)
                    {
                        reverseComplement(query+readAlnStart, readAlnSeq, readAlnLen);
                        reverseComplement(refAlnSeq, refAlnSeq_rev, refAlnLen);
                        edResult = edlibAlign(readAlnSeq, (readAlnEnd - readAlnEnd_new), refAlnSeq_rev, (refAlnEnd - refAlnEnd_new), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH));
                        // edlibGetCigarReverse(edResult.alignment, edResult.alignmentLength, EDLIB_CIGAR_STANDARD, edCigar);
                        edlibCigar_pushfront(edResult, edCigar);
                        edlibMD_pushfront(readAlnSeq, refAlnSeq_rev, edResult, edMD);
                        editScore -= edResult.editDistance;
                        edlibFreeAlignResult(edResult);
                    }
                    edCigar->insert(edCigar->begin(), readAlnEnd_new, 'I');
                    edMD->insert(edMD->begin(), readAlnEnd_new, '-');
                    // editScore_clip -= readAlnEnd_new;
                    tmpSam.flag = (isRev ? 16 : 0);
                    tmpSam.pos = refAlnEnd_new;
                    tmpSam.qStart = readAlnEnd_new;
                    numAnchorsSoFar = 0;
                }
                else
                {
                    editScore -= edResult.editDistance;
                    // edlibGetCigar(edResult.alignment, edResult.alignmentLength, EDLIB_CIGAR_STANDARD, edCigar);
                    edlibCigar_pushback(edResult, edCigar);
                    edlibMD_pushback(query+readAlnStart, refAlnSeq, edResult, edMD);
                    edlibFreeAlignResult(edResult);
                }
            }
            else
            {
                editScore -= edResult.editDistance;
                // edlibGetCigar(edResult.alignment, edResult.alignmentLength, EDLIB_CIGAR_STANDARD, edCigar);
                edlibCigar_pushback(edResult, edCigar);
                edlibMD_pushback(query+readAlnStart, refAlnSeq, edResult, edMD);
                edlibFreeAlignResult(edResult);
            }
        }
        else // no need to do smith-waterman, either insertion or deletion!
        {
            if(readAlnLen > 0)
            {
                edCigar->insert(edCigar->end(), readAlnLen, 'I');
                edMD->insert(edMD->end(), readAlnLen, '-');
                // update the total score; edit distance => unit score
                editScore -= readAlnLen;
            }
            else
            {
                edCigar->insert(edCigar->end(), refAlnLen, 'D');
                bwt_str_pac2char(refAlnStart, refAlnLen, refAlnSeq);
                for(j = 0; j < refAlnLen; j++)
                    edMD->push_back(refAlnSeq[j]);
                // update the total score; edit distance => unit score
                editScore -= refAlnLen;
            }
        }
        numAnchorsSoFar++;
    }
    // // last seed
    // fprintf(stderr, "\tseed\t%u\t%u\t%u\t%u\t%u\n", chain.seeds[i].qPos, chain.seeds[i].qPos + chain.seeds[i].len - 1, 
    //  chain.seeds[i].tPos, chain.seeds[i].tPos + chain.seeds[i].len - 1, chain.seeds[i].len);
    // fprintf(stderr, "\tseed\tquery\t%.*s\n", chain.seeds[i].len, query + chain.seeds[i].qPos);
    // // print reference
    // bwt_str_pac2int(chain.seeds[i].tPos, chain.seeds[i].len, refAlnSeq);
    // fprintf(stderr, "\tseed\t  ref\t");
    // for(j=0; j<chain.seeds[i].len; j++)
    //  fprintf(stderr, "%c", "ACGTN"[refAlnSeq[j]]);
    // fprintf(stderr, "\n");

    edCigar->insert(edCigar->end(), chain.seeds[i].len, 'M');
    // edCigar->insert(edCigar->begin(), chain.seeds[i].len, '=');
    edMD->insert(edMD->end(), chain.seeds[i].len, '=');

    // set sam posEnd
    tmpSam.posEnd = chain.seeds[i].tPos + chain.seeds[i].len - 1;
    tmpSam.qEnd = chain.seeds[i].qPos + chain.seeds[i].len - 1;

    // extend after last seed
    readAlnStart = chain.seeds[i].qPos + chain.seeds[i].len;
    readAlnLen = readLen - readAlnStart;
    refAlnLen = readAlnLen + 20; // TODO: how much extra sequence?
    if(readAlnLen > 0)
    {
        if(chain.seeds[i].tPos + chain.seeds[i].len + refAlnLen - 1 <= chrEnd)
        {
            refAlnStart = chain.seeds[i].tPos + chain.seeds[i].len;
            bwt_str_pac2char(refAlnStart, refAlnLen, refAlnSeq);
            
            edResult = edlibAlign(query+readAlnStart, readAlnLen, refAlnSeq, refAlnLen, edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH));
            
            LOG2({
                fprintf(stderr, "\tend\tqLen: %d\talnEdit: %d\talnSim: %f\n", 
                    readAlnLen, edResult.editDistance, (1 - ((float)edResult.editDistance / readAlnLen)));
            });

            if(readAlnLen > _pf_clipLen && (1 - ((float)edResult.editDistance / readAlnLen)) < _pf_clipSim)
            {
                convertChar2int(readAlnSeq_ksw, query+readAlnStart, readAlnLen);
                bwt_str_pac2int(refAlnStart, refAlnLen, refAlnSeq_ksw);

                ksw_extend(readAlnLen, readAlnSeq_ksw, refAlnLen, refAlnSeq_ksw, 5, _pf_kswMatrix_clip, _pf_kswGapOpen_clip, _pf_kswGapExtend_clip, 40, 0, 40, readAlnLen, &qLen_ksw, &tLen_ksw, 0, 0, 0);
                if(qLen_ksw > 0 && qLen_ksw < readAlnLen) // perform a new edlib alignment with new coordinates
                {
                    edlibFreeAlignResult(edResult); // free old results
                    edResult = edlibAlign(query+readAlnStart, qLen_ksw, refAlnSeq, tLen_ksw, edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH));
                    // edlibGetCigar(edResult.alignment, edResult.alignmentLength, EDLIB_CIGAR_STANDARD, edCigar);
                    edlibCigar_pushback(edResult, edCigar);
                    edlibMD_pushback(query+readAlnStart, refAlnSeq, edResult, edMD);
                    // update score
                    editScore -= edResult.editDistance;
                    // fix sam pos
                    tmpSam.posEnd = refAlnStart + edResult.endLocations[0];
                    tmpSam.qEnd = readAlnStart + qLen_ksw;

                    edCigar->insert(edCigar->end(), readAlnLen - qLen_ksw, 'I');
                    edMD->insert(edMD->end(), readAlnLen - qLen_ksw, '-');
                    // update score
                    // editScore_clip -= (readAlnLen - qLen_ksw);
                }
                else // use the already calculated results
                {
                    editScore -= edResult.editDistance;
                    // edlibGetCigar(edResult.alignment, edResult.alignmentLength, EDLIB_CIGAR_STANDARD, edCigar);
                    edlibCigar_pushback(edResult, edCigar);
                    edlibMD_pushback(query+readAlnStart, refAlnSeq, edResult, edMD);
                    // fix sam posEnd
                    tmpSam.posEnd = refAlnStart + edResult.endLocations[0];
                    tmpSam.qEnd = readLen;
                }
            }
            else
            {
                editScore -= edResult.editDistance;
                // edlibGetCigar(edResult.alignment, edResult.alignmentLength, EDLIB_CIGAR_STANDARD, edCigar);
                edlibCigar_pushback(edResult, edCigar);
                edlibMD_pushback(query+readAlnStart, refAlnSeq, edResult, edMD);
                // fix sam posEnd
                tmpSam.posEnd = refAlnStart + edResult.endLocations[0];
                tmpSam.qEnd = readLen;
            }

            edlibFreeAlignResult(edResult);
        }
        else // not enough sequence left on the chromosome to align => soft-clip
        {
            edCigar->insert(edCigar->end(), readAlnLen, 'I');
            edMD->insert(edMD->end(), readAlnLen, '-');
            // update score
            // editScore_clip -= readAlnLen;
        }
    }

    // fprintf(stderr, "score: %d\n", alnScore);
    // fprintf(stderr, "cigar: %s\n", soutCigar.str().c_str());

    // edCigar->push_back(0);  // Null character termination.
    // edCigarStr.assign(edCigar->begin(), edCigar->end());

    // fixCigar(alnCigar, edCigarStr);
    // tmpSam.cigar = edCigarStr;
    tmpSam.cigar = edlibCigar_toString(edCigar);
    tmpSam.md = edlibMD_toString(edMD, edCigar);
    // map.cigar = (char*) malloc((alnCigar.size() + 1) * sizeof(char));
    // strcpy(map.cigar, alnCigar.c_str());
    // tmpSam.alnScore = editScore + editScore_clip + readLen;
    // tmpSam.alnScore = editScore;
    tmpSam.nmCount = editScore;
    // tmpSam.clipCount = editScore_clip;

    map.samList.push_back(tmpSam);
    // map.totalScore += tmpSam.alnScore - editScore_split;
    // map.totalScore += tmpSam.alnScore;

    delete edCigar;
    delete edMD;

    // fprintf(stderr, "%s\n", alnCigar.c_str());
    // exit(0);
}
/**********************************************/
