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

#include "BWT.h"

#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include "bwa.h"

#include <string>
#include <vector>

bwaidx_t *_fmd_index;
bwtCache_t *_fmd_cacheTable;
int32_t kCache = 12;

extern "C"{
int bwa_index(int argc, char *argv[]);
}

int bwt_count_exact(const bwt_t *bwt, const char *str, int len, bwtint_t *sa_begin, bwtint_t *sa_end);
int64_t bwt_count_exact_cached(const bwt_t *bwt, const char *str, int len, bwtint_t *sa_begin, bwtint_t *sa_end);

// char charMap[4] = {'A', 'C', 'G', 'T'};
// std::vector<std::string> allKmers;

// void getKmer(std::string s, int k)
// {
//  if(k==kCache)
//  {
//    allKmers.push_back(s);
//    // fprintf(stdout, "%s\n", s.c_str());
//    return;
//  }
//  for(int i=0; i<4; i++)
//  {
//    getKmer(s+charMap[i], k+1);
//  }
// }

int bwt_cache_gen(char *ref_path)
{
    if ((_fmd_index = bwa_idx_load(ref_path, BWA_IDX_ALL)) == 0)
        return 1;

    const bwt_t *bwt = _fmd_index->bwt;

    int32_t i, j, k;
    int32_t cs = (int32_t) pow(4, kCache); // size of the total cache table
    int32_t os; // old size
    int32_t ni; // new index

    bwtint_t bk, bl, bok, bol; // bwt variables

    _fmd_cacheTable = (bwtCache_t*) malloc(cs * sizeof(bwtCache_t));
    if(_fmd_cacheTable == NULL)
    {
        return 1;
    }
    _fmd_cacheTable[0].beg = 0;
    _fmd_cacheTable[0].end = bwt->seq_len;

    // fprintf(stdout, "%lld (%llu, %llu)\n", _fmd_cacheTable[0].occ, _fmd_cacheTable[0].beg, _fmd_cacheTable[0].end);

    for(k = 0; k < kCache; k++)
    {
        os = (int32_t) pow(4, k);
        for(i = os - 1; i >= 0; i--)
        {
            if(_fmd_cacheTable[i].beg > _fmd_cacheTable[i].end)
            {
                for(j = 3; j >= 0; j--)
                {
                    // calculate the new index
                    ni = i * 4 + j;
                    // no need to find the new interval
                    _fmd_cacheTable[ni].beg = _fmd_cacheTable[i].beg;
                    _fmd_cacheTable[ni].end = _fmd_cacheTable[i].end;
                }
            }
            else
            {
                bk = _fmd_cacheTable[i].beg;
                bl = _fmd_cacheTable[i].end;
                for(j = 3; j >= 0; j--)
                {
                    // calculate the new index
                    ni = i * 4 + j;
                    // find the new interval
                    bwt_2occ(bwt, bk - 1, bl, j, &bok, &bol);
                    _fmd_cacheTable[ni].beg = bwt->L2[j] + bok + 1;
                    _fmd_cacheTable[ni].end = bwt->L2[j] + bol;
                }
            }
        }
    }

    char *cache_path;
    FILE *fp;
    cache_path = (char*) malloc(strlen(ref_path) + 10);
    strcpy(cache_path, ref_path);
    strcat(cache_path, ".cache");
    fp = fopen(cache_path, "wb");
    if(fp == NULL)
    {
        fprintf(stderr, "[BWT_CACHE_GEN] Could not open file: %s\n", cache_path);
        return 1;
    }

    fwrite(&kCache, sizeof(int32_t), 1, fp);
    fwrite(&cs, sizeof(int32_t), 1, fp);
    fwrite(_fmd_cacheTable, sizeof(bwtCache_t), cs, fp);
    fflush(fp);
    fclose(fp);

    free(cache_path);
    free(_fmd_cacheTable);
    return 0;
}

int bwt_index(char *ref_path)
{
    bwa_verbose = 2;
    // reset optind
    optind = 1;
    // prepare command line arguments
    char *argv_index[3] = {"dummy"};
    argv_index[1] = ref_path;
    if(bwa_index(2, argv_index))
        return 1;
    if(bwt_cache_gen(ref_path))
        return 1;
    return 0;
}

int bwt_cache_load(char *ref_path)
{
    int32_t cs;
    char *cache_path;
    FILE *fp;
    cache_path = (char*) malloc(strlen(ref_path) + 10);
    strcpy(cache_path, ref_path);
    strcat(cache_path, ".cache");
    fp = fopen(cache_path, "rb");
    if(fp == NULL)
    {
        fprintf(stderr, "[BWT_CACHE_LOAD] Could not open file: %s\n", cache_path);
        return 1;
    }

    kCache = 0;
    fread(&kCache, sizeof(int32_t), 1, fp);
    fread(&cs, sizeof(int32_t), 1, fp);
    _fmd_cacheTable = (bwtCache_t*) malloc(cs * sizeof(bwtCache_t));
    if(_fmd_cacheTable == NULL)
    {
        return 1;
    }
    fread(_fmd_cacheTable, sizeof(bwtCache_t), cs, fp);
    fclose(fp);

    free(cache_path);
    return 0;
}

int bwt_load(char *ref_path)
{
    bwa_verbose = 2;
    //
    char *idx_path;
    int len_path;
    FILE *fp;
    len_path = strlen(ref_path);
    idx_path = (char*) malloc(len_path + 10);
    strcpy(idx_path, ref_path);
    strcat(idx_path, ".bwt");
    if((fp = fopen(idx_path, "rb")) == NULL)
    {
        fprintf(stderr, "[BWT_LOAD] Could not locate index file: %s\n", idx_path);
        fprintf(stderr, "[BWT_LOAD] Try to build the index...\n");
        if(bwt_index(ref_path))
            return 1;
    }
    else
    {
        fclose(fp);
    }

    // if ((_fmd_index = bwa_idx_load(ref_path, BWA_IDX_BWT|BWA_IDX_BNS)) == 0)
    if ((_fmd_index = bwa_idx_load(ref_path, BWA_IDX_ALL)) == 0)
        return 1;

    if(bwt_cache_load(ref_path))
        return 1;

    // // sanity check for bwt cache
    // std::string kmer = "";
    // getKmer(kmer, 0);

    // bwtint_t sp, ep, occ;
    // for(int i = 0; i < allKmers.size(); i++)
    // {
    //  occ = bwt_count_exact(_fmd_index->bwt, allKmers[i].c_str(), kCache, &sp, &ep);
    //  fprintf(stdout, "%s %10lld %10llu %10llu\n", allKmers[i].c_str(), occ, sp, ep);
    // }
    // for(int i = 0; i < allKmers.size(); i++)
    // {
    //  occ = bwt_count_exact_cached(_fmd_index->bwt, allKmers[i].c_str(), kCache, &sp, &ep);
    //  fprintf(stderr, "%s %10lld %10llu %10llu\n", allKmers[i].c_str(), occ, sp, ep);
    // }
    // exit(0);
    // // end of sanity check

    free(idx_path);
    return 0;
}

int bwt_count_exact(const bwt_t *bwt, const char *str, int len, bwtint_t *sa_begin, bwtint_t *sa_end)
{
    bwtint_t k, l, ok, ol;
    ubyte_t c;
    int i;
    k = 0;
    l = bwt->seq_len;
    for (i = len - 1; i >= 0; --i)
    {
        c = nst_nt4_table[ str[i] ];
        if (c > 3) return 0; // no match
        bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
        k = bwt->L2[c] + ok + 1;
        l = bwt->L2[c] + ol;
        if (k > l) return 0; // no match
    }
    *sa_begin = k;
    *sa_end = l;
    return l - k + 1;
}

int64_t bwt_count_exact_cached(const bwt_t *bwt, const char *str, int len, bwtint_t *sa_begin, bwtint_t *sa_end)
{
    bwtint_t k, l, ok, ol;
    ubyte_t c;
    int i;
    int32_t idx = 0;
    for (i = len - 1; i >= len - kCache; --i)
    {
        c = nst_nt4_table[ str[i] ];
        if (c > 3) return 0; // no match
        // idx = (idx << 2) | c;
        idx = idx * 4 + c;
    }
    // *sa_begin = _fmd_cacheTable[idx].beg;
    // *sa_end = _fmd_cacheTable[idx].end;
    // return _fmd_cacheTable[idx].occ;

    if(_fmd_cacheTable[idx].beg > _fmd_cacheTable[idx].end) return 0;

    k = _fmd_cacheTable[idx].beg;
    l = _fmd_cacheTable[idx].end;
    for (i = len - kCache - 1; i >= 0; --i)
    {
        c = nst_nt4_table[ str[i] ];
        if (c > 3) return 0; // no match
        bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
        k = bwt->L2[c] + ok + 1;
        l = bwt->L2[c] + ol;
        if (k > l) return 0; // no match
    }
    *sa_begin = k;
    *sa_end = l;
    return l - k + 1;
}

uint64_t bwt_get_refGenBWTLen()
{
    return _fmd_index->bwt->seq_len;
}

uint32_t bwt_get_refGenLen()
{
    return _fmd_index->bns->l_pac;
}

#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)

void getLocs_extend_whole_step(char *qSeq, uint32_t qLen, uint32_t hash_count, SeedList *seedForward, SeedList *seedReverse)
{
    int i, m;
    bwtint_t j;
    bwtint_t sp, ep, occ;
    bwtint_t sp_tmp, ep_tmp, occ_tmp;
    bwtint_t sapos;

    double step = (double)qLen / hash_count;
    double seed_pos = 0;
    uint32_t seed_pos_int = 0;
    uint32_t last_pos = 0;

    uint32_t numForward = 0;
    uint32_t numReverse = 0;
  
    for(i=0; i<hash_count; i++)
    {
        m = MIN_ANCHOR_LEN;
        // occ = bwt_count_exact(_fmd_index->bwt, qSeq + seed_pos_int, m, &sp, &ep);
        // fprintf(stderr, "%.*s %llu %llu %llu\n", m, qSeq + seed_pos_int, occ, sp, ep);
        occ = bwt_count_exact_cached(_fmd_index->bwt, qSeq + seed_pos_int, m, &sp, &ep);
        // fprintf(stderr, "%.*s %llu %llu %llu\n", m, qSeq + seed_pos_int, occ, sp, ep);
        // while((occ_tmp = bwt_count_exact(_fmd_index->bwt, qSeq + seed_pos_int, m+1, &sp_tmp, &ep_tmp)) > 0)
        while((occ_tmp = bwt_count_exact_cached(_fmd_index->bwt, qSeq + seed_pos_int, m+1, &sp_tmp, &ep_tmp)) > 0)
        {
            occ = occ_tmp;
            sp = sp_tmp;
            ep = ep_tmp;
            m++;
        }

        // if there are some locations, the number of locations is less than MAX_NUM_HITS, and the seed is not contained
        if(occ > 0 && occ < MAX_NUM_HITS && (seed_pos_int + m) > last_pos)
        {
            // locate
            for(j = sp; j <= ep; j++)
            {
                sapos = bwt_sa(_fmd_index->bwt, j);
                if(sapos >= _fmd_index->bns->l_pac) // reverse strand
                {
                    sapos = (_fmd_index->bns->l_pac << 1) - sapos - m;

                    seedReverse->list[numReverse].tPos = sapos;
                    seedReverse->list[numReverse].qPos = qLen - seed_pos_int - m;
                    seedReverse->list[numReverse].len = m;
                    numReverse++;

                    // // print sequences
                    // fprintf(stderr, "Reverse: ml:%u ql:%u pos:%u qs:%u sa:%llu ts:%llu %.*s\n", m, qLen, seed_pos_int, qLen - seed_pos_int - m, j, sapos, m, qSeq + qLen - seed_pos_int - m);
                    // int64_t k;
                    // for (k = sapos; k < sapos + m; ++k)
                    //  fprintf(stderr, "%c", "ACGTN"[_get_pac(_fmd_index->pac, k)]);
                    // fprintf(stderr, "\n");
                    // for (k = sapos + m - 1; k >= sapos; --k)
                    //  fprintf(stderr, "%c", "ACGTN"[3 - _get_pac(_fmd_index->pac, k)]);
                    // fprintf(stderr, "\n");
                }
                else // forward strand
                {
                    seedForward->list[numForward].tPos = sapos;
                    seedForward->list[numForward].qPos = seed_pos_int;
                    seedForward->list[numForward].len = m;
                    numForward++;

                    // // print sequences
                    // fprintf(stderr, "Forward: ml:%u ql:%u pos:%u qs:%u sa:%llu ts:%llu %.*s\n", m, qLen, seed_pos_int, seed_pos_int, j, sapos, m, qSeq + seed_pos_int);
                    // int64_t k;
                    // for (k = sapos; k < sapos + m; ++k)
                    //  fprintf(stderr, "%c", "ACGTN"[_get_pac(_fmd_index->pac, k)]);
                    // fprintf(stderr, "\n");
                }
            }

            last_pos = seed_pos_int + m;
        }
        seed_pos += step;
        seed_pos_int = (uint32_t)seed_pos;
    }

    seedForward->num = numForward;
    seedReverse->num = numReverse;
}

int bwt_count_exact_backward(const bwt_t *bwt, const char *str, int ePos, bwtint_t *sa_begin, bwtint_t *sa_end, int *sPos)
{
    bwtint_t k, l;
    bwtint_t k_tmp, l_tmp;
    bwtint_t ok, ol;
    ubyte_t c;
    int i;
    k = 0;
    l = bwt->seq_len;
    for (i = ePos; i >= 0; --i)
    {
        c = nst_nt4_table[ str[i] ];
        if (c > 3) break; // no match
        bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
        k_tmp = bwt->L2[c] + ok + 1;
        l_tmp = bwt->L2[c] + ol;
        if (k_tmp > l_tmp) break; // no match
        k = k_tmp;
        l = l_tmp;
    }
    if(ePos - i < MIN_ANCHOR_LEN) return 0; // no match >= MIN_ANCHOR_LEN
    *sa_begin = k;
    *sa_end = l;
    *sPos = i + 1;
    return l - k + 1;
}

void getLocs_extend_whole_step2(char *qSeq, uint32_t qLen, uint32_t hash_count, SeedList *seedForward, SeedList *seedReverse)
{
    bwt_t *bwt = _fmd_index->bwt;
    int64_t l_pac = _fmd_index->bns->l_pac;

    bwtint_t j;
    bwtint_t occ;
    bwtint_t saLoc;
    bwtint_t k, l;
    int m;

    double step = (double)qLen / hash_count;
    double ePos_frac = qLen - 1;
    int ePos = qLen - 1;
    int last_pos = qLen;
    int sPos;

    uint32_t numForward = 0;
    uint32_t numReverse = 0;
  
    while(ePos >= MIN_ANCHOR_LEN - 1)
    {
        occ = bwt_count_exact_backward(_fmd_index->bwt, qSeq, ePos, &k, &l, &sPos);
        m = ePos - sPos + 1;

        // if there are some locations, the number of locations is less than MAX_NUM_HITS
        if(occ > 0 && occ < MAX_NUM_HITS && sPos < last_pos)
        {
            // locate
            for(j = k; j <= l; j++)
            {
                saLoc = bwt_sa(bwt, j);
                if(saLoc >= l_pac) // reverse strand
                {
                    saLoc = (l_pac << 1) - saLoc - m;

                    seedReverse->list[numReverse].tPos = saLoc;
                    seedReverse->list[numReverse].qPos = qLen - sPos - m;
                    seedReverse->list[numReverse].len = m;
                    numReverse++;

                    // // print sequences
                    // fprintf(stderr, "Reverse: ml:%u ql:%u pos:%u qs:%u sa:%llu ts:%llu %.*s\n", m, qLen, i + 1, qLen - i - 1 - m, j, saLoc, m, qSeq + qLen - i - 1 - m);
                    // int64_t z;
                    // for (z = saLoc; z < saLoc + m; ++z)
                    //  fprintf(stderr, "%c", "ACGTN"[_get_pac(_fmd_index->pac, z)]);
                    // fprintf(stderr, "\n");
                    // for (z = saLoc + m - 1; z >= saLoc; --z)
                    //  fprintf(stderr, "%c", "ACGTN"[3 - _get_pac(_fmd_index->pac, z)]);
                    // fprintf(stderr, "\n");
                }
                else // forward strand
                {
                    seedForward->list[numForward].tPos = saLoc;
                    seedForward->list[numForward].qPos = sPos;
                    seedForward->list[numForward].len = m;
                    numForward++;

                    // // print sequences
                    // fprintf(stderr, "Forward: ml:%u ql:%u pos:%u qs:%u sa:%llu ts:%llu %.*s\n", m, qLen, i + 1, i + 1, j, saLoc, m, qSeq + i + 1);
                    // int64_t z;
                    // for (z = saLoc; z < saLoc + m; z++)
                    //  fprintf(stderr, "%c", "ACGTN"[_get_pac(_fmd_index->pac, z)]);
                    // fprintf(stderr, "\n");
                }
            }
            last_pos = sPos;
        }
        ePos_frac -= step;
        ePos = (int)ePos_frac;
    }

    seedForward->num = numForward;
    seedReverse->num = numReverse;
}

void getLocs_extend_whole_step3(char *qSeq, uint32_t qLen, uint32_t hash_count, SeedList *seedForward, SeedList *seedReverse)
{
    typedef struct{int64_t k; int64_t l; int16_t m;} SaInterval;
    SaInterval allIntv[SEQ_MAX_LENGTH];
    int i;
    for(i = 0; i < SEQ_MAX_LENGTH; i++)
    {
        allIntv[i].k = -1;
        allIntv[i].l = -1;
        allIntv[i].m = 0;
    }
    //
    bwt_t *bwt = _fmd_index->bwt;
    int64_t l_pac = _fmd_index->bns->l_pac;
    bwtint_t k, l;
    bwtint_t ok, ol;
    ubyte_t c;
    int pos;
    for(pos = qLen - 1; pos >=0; pos--)
    {
        k = 0;
        l = bwt->seq_len;
        for(i = pos; i >= 0; i--)
        {
            c = nst_nt4_table[ qSeq[i] ];
            if (c > 3) break; // no match
            bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
            k = bwt->L2[c] + ok + 1;
            l = bwt->L2[c] + ol;
            if (k > l) break; // no match
            if(allIntv[i].k == -1 && allIntv[i].l == -1)
            {
                allIntv[i].k = k;
                allIntv[i].l = l;
                allIntv[i].m = pos - i + 1;
            }
        }
    }

    bwtint_t j;
    bwtint_t sapos;

    double step = (double)qLen / hash_count;
    double seed_pos = 0;
    uint32_t seed_pos_int = 0;
    uint32_t last_pos = 0;

    uint32_t numForward = 0;
    uint32_t numReverse = 0;
  
    for(i = 0; i < hash_count; i++)
    {
        // if there are some locations, the number of locations is less than MAX_NUM_HITS
        if(allIntv[seed_pos_int].m >= MIN_ANCHOR_LEN && allIntv[seed_pos_int].k != -1 && allIntv[seed_pos_int].l != -1 &&
            allIntv[seed_pos_int].l - allIntv[seed_pos_int].k + 1 < MAX_NUM_HITS &&
            (seed_pos_int + allIntv[seed_pos_int].m) > last_pos)
        {
            // locate
            for(j = allIntv[seed_pos_int].k; j <= allIntv[seed_pos_int].l; j++)
            {
                sapos = bwt_sa(bwt, j);
                if(sapos >= l_pac) // reverse strand
                {
                    sapos = (l_pac << 1) - sapos - allIntv[seed_pos_int].m;

                    seedReverse->list[numReverse].tPos = sapos;
                    seedReverse->list[numReverse].qPos = qLen - seed_pos_int - allIntv[seed_pos_int].m;
                    seedReverse->list[numReverse].len = allIntv[seed_pos_int].m;
                    numReverse++;

                    // // print sequences
                    // fprintf(stderr, "Reverse: ml:%u ql:%u pos:%u qs:%u sa:%llu ts:%llu %.*s\n", allIntv[seed_pos_int].m, qLen, seed_pos_int, qLen - seed_pos_int - allIntv[seed_pos_int].m, j, sapos, allIntv[seed_pos_int].m, qSeq + qLen - seed_pos_int - allIntv[seed_pos_int].m);
                }
                else // forward strand
                {
                    seedForward->list[numForward].tPos = sapos;
                    seedForward->list[numForward].qPos = seed_pos_int;
                    seedForward->list[numForward].len = allIntv[seed_pos_int].m;
                    numForward++;

                    // // print sequences
                    // fprintf(stderr, "Forward: ml:%u ql:%u pos:%u qs:%u sa:%llu ts:%llu %.*s\n", allIntv[seed_pos_int].m, qLen, seed_pos_int, seed_pos_int, j, sapos, allIntv[seed_pos_int].m, qSeq + seed_pos_int);
                }
            }
            last_pos = seed_pos_int + allIntv[seed_pos_int].m;
        }
        seed_pos += step;
        seed_pos_int = (uint32_t)seed_pos;
    }

    seedForward->num = numForward;
    seedReverse->num = numReverse;
}

void bwt_str_pac2int(uint32_t beg, uint32_t len, uint8_t *seq)
{
    uint32_t k;
    uint32_t l = 0;
    for (k = beg; k < beg + len; k++)
        seq[l++] = _get_pac(_fmd_index->pac, k);
}

void bwt_str_pac2char(uint32_t beg, uint32_t len, char *seq)
{
    uint32_t k;
    uint32_t l = 0;
    for (k = beg; k < beg + len; k++)
        seq[l++] = "ACGTN"[_get_pac(_fmd_index->pac, k)];
}

// void bwt_get_intv_info(uint64_t beg, uint64_t end, char **chr_name, int32_t *chr_len, uint64_t *chr_beg, uint64_t *chr_end, int *is_rev)
// {
//  uint64_t mid = (beg + end) >> 1;

//  //TODO: what if beg < l_pac && end > l_pac
//  if(beg >= _fmd_index->bns->l_pac) // reverse strand
//  {
//    *chr_beg = (_fmd_index->bns->l_pac << 1) - 1 - end;
//    *chr_end = (_fmd_index->bns->l_pac << 1) - 1 - beg;
//    // for (k = end_f; k > beg_f; --k)
//    //  seq[l++] = 3 - _get_pac(pac, k);
//  }
//  else // forward strand
//  {
//    *chr_beg = beg;
//    *chr_end = end;
//    // for (k = beg; k < end; ++k)
//    //  seq[l++] = _get_pac(pac, k);
//  }

//  int rid = bns_pos2rid(_fmd_index->bns, bns_depos(_fmd_index->bns, mid, is_rev));
//  *chr_beg -= _fmd_index->bns->anns[rid].offset;
//  *chr_end -= _fmd_index->bns->anns[rid].offset;
//  *chr_name = _fmd_index->bns->anns[rid].name;
//  *chr_len = _fmd_index->bns->anns[rid].len;
// }

void bwt_get_intv_info(uint64_t beg, uint64_t end, char **chr_name, int32_t *chr_len, uint32_t *chr_beg, uint32_t *chr_end)
{
    int rid;
    // rid = bns_pos2rid(_fmd_index->bns, beg);
    // fprintf(stderr, "beg: %s\n", _fmd_index->bns->anns[rid].name);
    uint64_t mid = (beg + end) >> 1;

    *chr_beg = beg;
    *chr_end = end;

    rid = bns_pos2rid(_fmd_index->bns, mid);
    *chr_beg -= _fmd_index->bns->anns[rid].offset;
    *chr_end -= _fmd_index->bns->anns[rid].offset;
    *chr_name = _fmd_index->bns->anns[rid].name;
    *chr_len = _fmd_index->bns->anns[rid].len;
}

void bwt_get_chr_boundaries(uint64_t beg, uint64_t end, uint32_t *chr_beg, uint32_t *chr_end)
{
    int rid;
    uint64_t mid = (beg + end) >> 1;
    rid = bns_pos2rid(_fmd_index->bns, mid);

    *chr_beg = _fmd_index->bns->anns[rid].offset;
    *chr_end = _fmd_index->bns->anns[rid].offset + _fmd_index->bns->anns[rid].len - 1;

    // fprintf(stderr, "chrName: %s\n", _fmd_index->bns->anns[rid].name);
    // fprintf(stderr, "chrLen: %s\n", _fmd_index->bns->anns[rid].len);
    // fprintf(stderr, " chrBeg: %s\n", *chr_beg);
    // fprintf(stderr, " chrEnd: %s\n", *chr_end);
}

void printSamHeader(FILE *fp)
{
    int i;
    fprintf(fp, "@HD\tVN:1.5\tSO:unsorted\n");
    for(i=0; i<_fmd_index->bns->n_seqs; i++)
    {
        fprintf(fp, "@SQ\tSN:%s\tLN:%d\n", _fmd_index->bns->anns[i].name, _fmd_index->bns->anns[i].len);
    }
    fprintf(fp, "@PG\tID:pacfast\tPN:pacfast\tVN:%s\tCL:%s\n", PROG_VERSION, opt_commandAll);
}
