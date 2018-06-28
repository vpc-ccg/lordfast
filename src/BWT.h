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

#ifndef __BWT_FMD__
#define __BWT_FMD__

#include "LordFAST.h"

int         bwt_index(char *ref_path);
int         bwt_load(char *ref_path);
uint64_t    bwt_get_refGenBWTLen();
uint32_t    bwt_get_refGenLen();
void        getLocs_extend_whole_step(char *qSeq, uint32_t qLen, uint32_t hash_count, SeedList *seedForward, SeedList *seedReverse);
void        getLocs_extend_whole_step2(char *qSeq, uint32_t qLen, uint32_t hash_count, SeedList *seedForward, SeedList *seedReverse);
void        getLocs_extend_whole_step3(char *qSeq, uint32_t qLen, uint32_t hash_count, SeedList *seedForward, SeedList *seedReverse);
void        bwt_get_intv_info(uint64_t beg, uint64_t end, char **chr_name, int32_t *chr_len, uint32_t *chr_beg, uint32_t *chr_end);
void        bwt_get_chr_boundaries(uint64_t beg, uint64_t end, uint32_t *chr_beg, uint32_t *chr_end);
void        bwt_str_pac2int(uint32_t beg, uint32_t len, uint8_t *seq);
void        bwt_str_pac2char(uint32_t beg, uint32_t len, char *seq);
void        printSamHeader(FILE *fp);

typedef struct
{
    uint64_t beg;
    uint64_t end;
    // uint8_t occ;
} bwtCache_t;

#endif //__BWT_FMD__
