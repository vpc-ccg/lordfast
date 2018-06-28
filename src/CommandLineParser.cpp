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
#include <getopt.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include "Common.h"
#include "CommandLineParser.h"

int                     uniqueMode=1;
int                     indexingMode;
int                     searchingMode;
// int                     bestMappingMode = 1;
int                     affineMode = 0;
int                     seqCompressed = 1;
int                     outCompressed;
int                     progressRep = 0;
int                     nohitDisabled = 1;
int                     noSamHeader = 0;
int                     tabOutput = 0;
// int                     errThreshold = -1;
char                    *seqFile;
char                    fileName[3][1000];
char                    outputMap[1000];
// char                    mappingOutputPath[1000];
char                    outputUnmap[1000];
char                    opt_commandAll[2000];
int                     opt_outputBufferSize = 2000000;
unsigned char           MIN_ANCHOR_LEN = 14;
// unsigned int            SAMPLING_INTERVAL = 2000;
unsigned int            SAMPLING_COUNT = 1000;
unsigned int            MAX_MAP = 10;
unsigned int            CHUNK_SIZE = 2000;
unsigned int            MIN_READ_LEN = 1000;
unsigned int            MAX_NUM_HITS = 1000;
unsigned int            CONTIG_SIZE;
unsigned int            CONTIG_MAX_SIZE;
unsigned int            THREAD_COUNT = 1;
double                  MAX_MEMORY = 4;// GB
int                     THREAD_ID[255];
// chainAlg_t              chainAlg = CHAIN_ALG_CLASP;
chainAlg_t              chainAlg = CHAIN_ALG_DPN2;
double                  chainReward = 9.3;
double                  chainPenalty = 11.4;
double                  gapPenalty = 0.15;

#if (defined(__MACH__) && defined(__APPLE__))
#include <mach-o/getsect.h>
#else
extern char _binary_HELP_start;
extern char _binary_HELP_end;
#endif

void printHelp()
{
#if (defined(__MACH__) && defined(__APPLE__))
    size_t i, sz = getsectbyname("binary", "HELP")->size;
    const uint8_t *c =  (const uint8_t*) getsectbyname("binary", "HELP")->addr;
    for (i = 0; i < sz; i++) 
        putchar(c[i]); 
#else
    char *c;
    for (c = &_binary_HELP_start; c != &_binary_HELP_end; c++)
        putchar(*c);
#endif
}

int parseCommandLine (int argc, char *argv[])
{
    int index, len, o;
    char *fastaFile = NULL;
    char *fastaOutputFile = NULL;
    char *indexFile = NULL;
    std::string optarg_str;

    // mappingOutput = (char*) getMem(FILE_NAME_LENGTH);
    // mappingOutputPath = (char*) getMem(FILE_NAME_LENGTH);
    // unmappedOutput = (char*) getMem(FILE_NAME_LENGTH);
    // strcpy(mappingOutput, "output");
    // strcpy(unmappedOutput, "output.nohit");
    // mappingOutputPath[0] = '\0';

    static struct option longOptions[] = 
    {
        {"index",                   required_argument,  0,                  'I'},
        {"search",                  required_argument,  0,                  'S'},
        {"seq",                     required_argument,  0,                  's'},
        {"out",                     required_argument,  0,                  'o'},
        {"threads",                 required_argument,  0,                  't'},
        {"minAnchorLen",            required_argument,  0,                  'k'},
        {"maxRefHit",               required_argument,  0,                  'm'},
        {"minReadLen",              required_argument,  0,                  'l'},
        {"anchorCount",             required_argument,  0,                  'c'},
        {"numMap",                  required_argument,  0,                  'n'},
        {"chainAlg",                required_argument,  0,                  'A'},
        {"progress",                no_argument,        &progressRep,       1},
        {"noSamHeader",             no_argument,        &noSamHeader,       1},
        {"tabOutput",               no_argument,        &tabOutput,         1},
        {"affine",                  no_argument,        &affineMode,        1},
        {"unMapped",                required_argument,  0,                  'u'},
        {"err",                     required_argument,  0,                  'e'},
        {"chainReward",             required_argument,  0,                  'R'},
        {"chainPenalty",            required_argument,  0,                  'P'},
        {"gapPenalty",              required_argument,  0,                  'G'},
        {"help",                    no_argument,        0,                  'h'},
        {"version",                 no_argument,        0,                  'v'},
        // {"sl",                       required_argument,  0,                  'l'},
        // {"mem",                      required_argument,  0,                  'z'},
        {0,0,0,0}
    };



    while ( (o = getopt_long ( argc, argv, "I:S:s:o:t:k:m:l:c:n:A:u:e:R:P:G:hv", longOptions, &index))!= -1 )
    {
        switch (o)
        {
            case 0:
                fprintf(stderr, "[NOTE] option is set: %s\n", longOptions[index].name);
                break;
            case 'I':
                indexingMode = 1;
                fastaFile = optarg;
                break;
            case 'S':
                searchingMode = 1;
                fastaFile = optarg;
                break;
            case 's':
                seqFile = optarg;
                break;
            case 'o':
                strcpy(outputMap, optarg);
                // stripPath (optarg, mappingOutputPath, mappingOutput);
                // sprintf(unmappedOutput, "%s%s.nohit", mappingOutputPath, mappingOutput );
                // fprintf(stderr, "[CHECK] %s\n", mappingOutputPath);
                // fprintf(stderr, "[CHECK] outputMap: %s\n", outputMap);
                break;
            case 't':
                THREAD_COUNT = atoi(optarg);
                if (THREAD_COUNT == 0 || THREAD_COUNT > sysconf( _SC_NPROCESSORS_ONLN ))
                    THREAD_COUNT = sysconf( _SC_NPROCESSORS_ONLN );
                break;
            case 'k':
                MIN_ANCHOR_LEN = atoi(optarg);
                break;
            case 'm':
                MAX_NUM_HITS = atoi(optarg);
                break;
            case 'l':
                MIN_READ_LEN = atoi(optarg);
                break;
            case 'c':
                SAMPLING_COUNT = atoi(optarg);
                break;
            case 'n':
                MAX_MAP = atoi(optarg);
                break;
            case 'A':
                optarg_str = optarg;
                if(optarg_str == "clasp")
                {
                    chainAlg = CHAIN_ALG_CLASP;
                }
                else if(optarg_str == "dp-n2")
                {
                    chainAlg = CHAIN_ALG_DPN2;
                }
                // else if(optarg_str == "lis-n2")
                // {
                //     chainAlg = CHAIN_ALG_LISN2;
                // }
                // else if(optarg_str == "lis-nlogn")
                // {
                //     chainAlg = CHAIN_ALG_LISNLOGN;
                // }
                else
                {
                    fprintf(stderr, "Unknown argument for chainAlg. Using clasp algorithm!\n");
                    // chainAlg = CHAIN_ALG_CLASP;
                    chainAlg = CHAIN_ALG_DPN2;
                }
                break;
            case 'u':
                strcpy(outputUnmap, optarg);
                // fprintf(stderr, "[CHECK] outputUnmap %s\n", outputUnmap);
                break;
            // case 'e':
            //     errThreshold = atoi(optarg);
            //     break;
            case 'R':
                chainReward = atof(optarg);
                break;
            case 'P':
                chainPenalty = atof(optarg);
                break;
            case 'G':
                gapPenalty = atof(optarg);
                break;
            case 'h':
                printHelp();
                exit(EXIT_SUCCESS);
                break;
            case 'v':
                fprintf(stdout, "lordFAST %s\n", PROG_VERSION);
                exit(EXIT_SUCCESS);
                break;
            // case 'l':
            //  SAMPLING_INTERVAL = atoi(optarg);
            //  break;
            // case 'z':
            //  MAX_MEMORY = atoi(optarg);
            //  break;
            default:
                return 1;
        }

    }

// #ifndef MRSFAST_SSE4
//  if (searchingMode)
//      fprintf(stdout, "==> This version is compiled without any SSE4 optimization <==\n");
// #endif

    // if (bestMappingMode)
    // {
    //  nohitDisabled = 1;
    // }

    if(tabOutput)
    {
        noSamHeader = 1;
    }

    if (indexingMode + searchingMode != 1)
    {
        fprintf(stderr, "[parseCommandLine] ERROR: Indexing / Searching mode should be selected\n");
        return 1;
    }

    if (MIN_ANCHOR_LEN > 20 || MIN_ANCHOR_LEN < 10)
    {
        fprintf(stderr, "[parseCommandLine] ERROR: mininum anchor length should be in [10..20]\n");
        return 1;
    }

    // if(SAMPLING_COUNT > SAMPLING_INTERVAL)
    // {
    //  fprintf(stderr, "[parseCommandLine] ERROR: Seed count cannot be higher than segment length\n");
    //  return 1;
    // }

    // if (MAX_MEMORY < 4)
    //  fprintf(stderr, "ERROR: At least 4 GB of memory is required for running pacFAST\n");

    if ( indexingMode )
    {
        CONTIG_SIZE     = 80000000;
        CONTIG_MAX_SIZE = 120000000;
        // CONTIG_SIZE      = 300000000;
        // CONTIG_MAX_SIZE  = 300000000;

        if (fastaFile == NULL)
        {
            fprintf(stderr, "[parseCommandLine] ERROR: Reference(s) should be indicated for indexing\n");
            return 1;
        }
    }

    // if (maxHits)
    // {
    //  if (maxHits < 0)
    //  {
    //      fprintf(stderr, "[parseCommandLine] ERROR: Number of maximum hits must be greater than 0\n");
    //      return 1;
    //  }

    //  if (bestMappingMode)
    //  {
    //      fprintf(stderr, "[parseCommandLine] ERROR: Maximum number of mappings could not be set in best mapping mode. Maximum mappings input ignored\n");
    //      maxHits = 0;
    //  }
    // }

    if ( searchingMode )
    {
        CONTIG_SIZE     = 300000000;
        CONTIG_MAX_SIZE = 300000000;

        
        if (fastaFile == NULL)
        {
            fprintf(stderr, "[parseCommandLine] ERROR: Index File(s) should be indiciated for searching\n");
            return 1;
        }

        if (seqFile == NULL)
        {
            fprintf(stderr, "[parseCommandLine] ERROR: Please indicate a sequence file for searching.\n");
            return 1;
        }
    }

    int i = 0;

    sprintf(fileName[0], "%s", fastaFile);
    sprintf(fileName[1], "%s.index", fileName[0]);
    sprintf(fileName[2], "%s.packed", fileName[0]);

    if (!indexingMode)
    {
        fprintf(stderr, "[NOTE] number of threads: %d\n", THREAD_COUNT);
        for (i = 0; i < 255; i++)
            THREAD_ID[i] = i;
    }

    for(i = 0; i < argc; i++)
    {
        strcat(opt_commandAll, argv[i]);
        strcat(opt_commandAll, " ");
    }

    initCommon();
    return 0;
}
/**********************************************/
// void finalizeCommandParser()
// {
//  freeMem(mappingOutput, FILE_NAME_LENGTH);
//  freeMem(unmappedOutput, FILE_NAME_LENGTH);
//  freeMem(mappingOutputPath, FILE_NAME_LENGTH);
// }
