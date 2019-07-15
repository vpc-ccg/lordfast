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

int                     indexingMode;
int                     searchingMode;
// int                     affineMode = 0;
int                     noSamHeader = 0;
char                    *seqFile;
char                    *refFile;
char                    outputMap[1000];
char                    opt_commandAll[2000];
int                     opt_outputBufferSize = 2000000;
chainAlg_t              chainAlg = CHAIN_ALG_DPN2;
char                    readGroup[1000];
char                    readGroupId[1000];
double                  chainReward = 9.3;
double                  chainPenalty = 11.4;
double                  gapPenalty = 0.15;
long long               memUsage = 0;

int                     THREAD_COUNT = 1;
int                     THREAD_ID[255];
int                     MIN_ANCHOR_LEN = 14;
int                     SAMPLING_COUNT = 1000;
int                     MAX_MAP = 10;
int                     MIN_READ_LEN = 1000;
int                     MAX_REF_HITS = 1000;

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

void printHelp_short()
{
    fprintf(stderr, "usage: lordfast --index ref.fa\n");
    fprintf(stderr, "       lordfast --search ref.fa --seq reads.fa [options]\n");
    fprintf(stderr, "For more details and command line options run \"lordfast --help\"\n");
}

int set_read_group(char *rg_line)
{
    char *p, *q;
    if(strstr(rg_line, "@RG") != rg_line)
    {
        fprintf(stderr, "[ERROR] SAM read group line does not start with @RG\n");
        return 1;
    }
    if(strstr(rg_line, "\t") != NULL)
    {
        fprintf(stderr, "[ERROR] the read group line contained literal <tab> characters. Please replace with escaped tabs: \\t\n");
        return 1;
    }
    // convert escape characters
    for(p = rg_line, q = readGroup; *p; p++)
    {
        if(*p != '\\')
        {
            *q++ = *p;
        }
        else
        {
            p++;
            if(*p == 't') *q++ = '\t';
            else if(*p == 'n') *q++ = '\n';
            else if(*p == 'r') *q++ = '\r';
            else if(*p == '\\') *q++ = '\\';
        }
    }
    *q = '\0';
    // find the read group id
    if((p = strstr(readGroup, "ID:")) == 0)
    {
        fprintf(stderr, "[ERROR] no ID within the read group line\n");
        return 1;
    }
    // copy the read group id
    for(p = p+3, q = readGroupId; *p && *p != '\n' && *p != '\t'; ) *q++ = *p++;
    return 0;
}

int parseCommandLine (int argc, char *argv[])
{
    if(argc < 2)
    {
        printHelp_short();
        return 1;
    }
    int i;
    int index, ch;

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
        {"chainAlg",                required_argument,  0,                  'a'},
        {"readGroup",               required_argument,  0,                  'R'},
        {"noSamHeader",             no_argument,        &noSamHeader,         1},
        // {"affine",                  no_argument,        &affineMode,        1},
        {"chainReward",             required_argument,  0,                  'r'},
        {"chainPenalty",            required_argument,  0,                  'p'},
        {"gapPenalty",              required_argument,  0,                  'g'},
        {"help",                    no_argument,        0,                  'h'},
        {"version",                 no_argument,        0,                  'v'},
        {0,0,0,0}
    };

    while ( (ch = getopt_long ( argc, argv, "I:S:s:o:t:k:m:l:c:n:a:r:R:P:G:hv", longOptions, &index))!= -1 )
    {
        switch (ch)
        {
            case 0:
                fprintf(stderr, "[NOTE] option %s is set\n", longOptions[index].name);
                break;
            case 'I':
                indexingMode = 1;
                refFile = optarg;
                break;
            case 'S':
                searchingMode = 1;
                refFile = optarg;
                break;
            case 's':
                seqFile = optarg;
                break;
            case 'o':
                strcpy(outputMap, optarg);
                break;
            case 't':
                THREAD_COUNT = atoi(optarg);
                if (THREAD_COUNT <= 0 || THREAD_COUNT > sysconf( _SC_NPROCESSORS_ONLN ))
                    THREAD_COUNT = sysconf( _SC_NPROCESSORS_ONLN );
                break;
            case 'k':
                MIN_ANCHOR_LEN = atoi(optarg);
                break;
            case 'm':
                MAX_REF_HITS = atoi(optarg);
                break;
            case 'l':
                MIN_READ_LEN = atoi(optarg);
                if(MIN_READ_LEN < 100) MIN_READ_LEN = 100;
                break;
            case 'c':
                SAMPLING_COUNT = atoi(optarg);
                break;
            case 'n':
                MAX_MAP = atoi(optarg);
                break;
            case 'a':
                if(strcmp(optarg, "clasp") == 0)
                {
                    chainAlg = CHAIN_ALG_CLASP;
                }
                else if(strcmp(optarg, "dp-n2") == 0)
                {
                    chainAlg = CHAIN_ALG_DPN2;
                }
                else
                {
                    fprintf(stderr, "[WARNING] (parseCommandLine) unknown argument for -A/--chainAlg. Using dynamic programming (dp-n2)!\n");
                    chainAlg = CHAIN_ALG_DPN2;
                }
                break;
            case 'R':
                if(set_read_group(optarg))
                    return 1;
                break;
            case 'r':
                chainReward = atof(optarg);
                break;
            case 'p':
                chainPenalty = atof(optarg);
                break;
            case 'g':
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
            default:
                printHelp_short();
                return 1;
        }

    }

    if (indexingMode + searchingMode != 1)
    {
        fprintf(stderr, "[ERROR] (parseCommandLine) indexing / searching mode should be selected\n");
        printHelp_short();
        return 1;
    }

    if ( indexingMode && refFile == NULL )
    {
        fprintf(stderr, "[ERROR] (parseCommandLine) reference file should be indicated for indexing\n");
        printHelp_short();
        return 1;
    }
    if ( searchingMode )
    {
        if (refFile == NULL)
        {
            fprintf(stderr, "[ERROR] (parseCommandLine) reference file should be indiciated for searching\n");
            printHelp_short();
            return 1;
        }
        if (seqFile == NULL)
        {
            fprintf(stderr, "[ERROR] (parseCommandLine) please indicate a sequence file for searching.\n");
            printHelp_short();
            return 1;
        }
    }

    if (MIN_ANCHOR_LEN < 12 || MIN_ANCHOR_LEN > 20)
    {
        fprintf(stderr, "[ERROR] (parseCommandLine) -k/--minAnchorLen requires an argument in [12..20]\n");
        return 1;
    }
    if(SAMPLING_COUNT <= 0)
    {
        fprintf(stderr, "[ERROR] (parseCommandLine) -c/--anchorCount requires a positive integer argument\n");
        return 1;
    }
    if(MAX_MAP <= 0)
    {
        fprintf(stderr, "[ERROR] (parseCommandLine) -n/--numMap requires a positive integer argument\n");
        return 1;
    }
    if(MAX_REF_HITS <= 0)
    {
        fprintf(stderr, "[ERROR] (parseCommandLine) -m/--maxRefHit requires a positive integer argument\n");
        return 1;
    }

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

    return 0;
}
