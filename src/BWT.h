/*
 * Author: Ehsan Haghshenas (ehaghshe AT sfu DOT ca)
 */


#ifndef __BWT_FMD__
#define __BWT_FMD__

#include "LordFAST.h"

int      bwt_index(char *ref_path);
int      bwt_load(char *ref_path);
uint64_t bwt_get_refGenBWTLen();
uint32_t bwt_get_refGenLen();
void getLocs_extend_whole_step(char *qSeq, uint32_t qLen, uint32_t hash_count, SeedList *seedForward, SeedList *seedReverse);
void getLocs_extend_whole_step2(char *qSeq, uint32_t qLen, uint32_t hash_count, SeedList *seedForward, SeedList *seedReverse);
void getLocs_extend_whole_step3(char *qSeq, uint32_t qLen, uint32_t hash_count, SeedList *seedForward, SeedList *seedReverse);
void bwt_get_intv_info(uint64_t beg, uint64_t end, char **chr_name, int32_t *chr_len, uint32_t *chr_beg, uint32_t *chr_end);
void bwt_str_pac2int(uint32_t beg, uint32_t len, uint8_t *seq);
void bwt_str_pac2char(uint32_t beg, uint32_t len, char *seq);
void printSamHeader(FILE *fp);

typedef struct
{
	uint64_t beg;
	uint64_t end;
	// uint8_t occ;
} bwtCache_t;

#endif //__BWT_FMD__
