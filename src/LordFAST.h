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

#ifndef __PACFAST_BWT__
#define __PACFAST_BWT__

#include <string>
#include <vector>
#include "Reads.h"

typedef struct
{
    uint32_t tPos;
    uint32_t qPos : 20;
    uint32_t len  : 12;
} Seed_t;

typedef struct
{
    Seed_t *list;
    uint32_t num; // size of the list
} SeedList;

typedef struct
{
    uint32_t tStart;
    uint32_t tEnd;
    uint8_t isReverse;
    float score;
} Win_t;

typedef struct
{
    Win_t *list;
    uint32_t num; // size of the list
} WinList;

typedef struct
{
    int readIdx;
    uint32_t cnt;
} WinCount_t;

bool compareWin(const Win_t& w1, const Win_t& w2);
bool compareSeed(const Seed_t& s1, const Seed_t& s2);

typedef struct
{
    Seed_t *seeds; // list of seeds contained in the chain
    uint32_t chainLen; // chain length
    float score;
} Chain_t;

typedef struct
{
    Chain_t *list;
    uint32_t num; // size of the list
} ChainList;

// bool compareChain(const Chain_t& c1, const Chain_t& c2);

class Sam_t
{
public:
    char *rName;
    // int32_t rLen;
    uint32_t rStart;
    // uint32_t rEnd;
    uint32_t qStart;
    uint32_t qEnd;
    uint16_t flag;
    int16_t mapQ;
    uint32_t pos;
    uint32_t posEnd;
    int32_t alnScore;
    int32_t nmCount;
    // char *cigar;
    std::string cigar;
    std::string md;
    std::string sa;
};

class SamList_t
{
public:
    std::vector<Sam_t> samList; // a vector of sam elements; useful for split mappings
    int32_t totalScore; // sum of the score of all splits
};

class MapInfo
{
public:
    char *qName;
    char *seq;
    char *seq_rev;
    char *qual;
    char *qual_rev;
    SamList_t *mappings;
};

bool compareSam(const SamList_t& s1, const SamList_t& s2);

void initializeFAST();
void finalizeFAST();
void initFASTChunk(Read *seqList, int seqListSize);
// void finalizeFASTChunk();
void mapSeqMT();

#endif /*__PACFAST_BWT__*/
