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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Common.h"
#include "CommandLineParser.h"
// #include "Output.h"
#include "BWT.h"
// #include "PacFAST-BWT.h"
#include "LordFAST.h"

// unsigned char  seqFastq;
// void     printStat();

int main(int argc, char *argv[])
{
    if (parseCommandLine(argc, argv))
        return EXIT_FAILURE;

    // INDEXING
    if (indexingMode)
    {
        if(bwt_index(fileName[0]))
            return EXIT_FAILURE;
    }
    // SEARCHING
    else
    {
        // queryTest();
        Read *seqList;
        unsigned int seqListSize;
        int totalNumOfReads = 0;
        double totalLoadingTime = 0;
        double totalMappingTime = 0;
        double startTime;
        double loadingTime;
        double mappingTime;
        double lstartTime;
        double tmpTime;
        //  double maxMem=0;
        //  int flag;
        //  int chrIndex;

        // Loading BWT-FM index
        startTime = getTime();
        if(bwt_load(fileName[0]))
            return EXIT_FAILURE;

        // if (!initRead(seqFile, 1500000000))
        // if (!initRead(seqFile, 10000000))
        if (!initRead(seqFile, 500000000))
            return EXIT_FAILURE;

        totalLoadingTime += getTime()-startTime;

        // // Preparing output
        // initOutput(mappingOutput, outCompressed);

        fprintf(stderr, "-----------------------------------------------------------------------------------------------------------\n");
        fprintf(stderr, "| %15s | %15s | %15s | %15s | %15s %15s |\n","Genome Name","Loading Time", "Mapping Time", "Memory Usage(M)","Total Mappings","Mapped reads");
        fprintf(stderr, "-----------------------------------------------------------------------------------------------------------\n");

        mappingTime = 0;
        loadingTime = 0;
        //  flag = 1;

        initializeFAST();

        tmpTime = getTime();
        while (readChunk(&seqList, &seqListSize) || seqListSize > 0)
        {
            initFASTChunk(seqList, seqListSize);
            totalLoadingTime += (getTime() - tmpTime);  // readAllReads + initLoadingHashTable
            totalNumOfReads += seqListSize;

            lstartTime = getTime();     
            mapSeqMT();
            mappingTime += getTime() - lstartTime;
            totalMappingTime += mappingTime;

            finalizeFASTChunk();

            // fprintf(stdout, "| %15s | %15.2f | %15.2f | %15.2f | %15lld %15lld |\n",
            //     chrName, loadingTime, mappingTime, maxMem, mappingCnt , mappingCnt);
            // fflush(stdout);
            
        //    chrIndex = 0;
        //    do
        //    {
        //         flag = loadHashTable ( &tmpTime );        // Reading a fragment
        //         loadingTime += tmpTime;

        //         char *chrName_tmp = getRefGenomeName();
        //         char *chrName = getMem(CONTIG_NAME_SIZE);
        //         strcpy(chrName, chrName_tmp);

        //         initFASTContig();
        //         mapSeq(flag, chrIndex);

        //         if (maxMem < getMemUsage())
        //             maxMem = getMemUsage();

        //         if (flag == 0 || flag == 2)
        //         {
        //             totalMappingTime += mappingTime;
        //             totalLoadingTime += loadingTime;

        //             loadingTime = 0;
        //             mappingTime = 0;
        //             maxMem = 0;
        //             chrIndex++;
        //         }
        //         else if (progressRep)
        //         {
        //             fprintf(stdout, "| %15s | %15.2f | %15.2f | %15.2f | %15lld %15lld |\n",
        //                 chrName, loadingTime, mappingTime, maxMem, mappingCnt , mappingCnt);
        //             fflush(stdout);
        //         }

        //         freeMem(chrName, CONTIG_NAME_SIZE);
        //    } while (flag);

            releaseChunk();
            // tmpTime = getTime();
        }
        totalLoadingTime += (getTime() - tmpTime);    // for the last readAllReads call

        finalizeFAST();
        // // finalizeLoadingHashTable();
        // finalizeReads();
        // finalizeOutput();
        // finalizeCommandParser();

        //  fprintf(stdout, "----------------------------------------------------------------------------------------------------------\n");

        //  fprintf(stdout, "%19s%16.2f%18.2f\n\n", "Total:",totalLoadingTime, totalMappingTime);
        //  fprintf(stdout, "%-30s%10.2f\n","Total Time:", totalMappingTime+totalLoadingTime);
        //  fprintf(stdout, "%-30s%10d\n","Total No. of Reads:", totalNumOfReads);
        //  fprintf(stdout, "%-30s%10lld\n","Total No. of Mappings:", mappingCnt);
        //  // fprintf(stdout, "%-30s%10.0f\n","Avg No. of locations verified:", ceil((float)verificationCnt/totalNumOfReads));
        //  if (memUsage > 0)
        //     fprintf(stdout, "Memory Leak: %lld Bytes\n", memUsage);

    }
    return 0;
}
