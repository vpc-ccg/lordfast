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
#include "Common.h"
#include "CommandLineParser.h"
#include "BWT.h"
#include "LordFAST.h"

int main(int argc, char *argv[])
{
    if (parseCommandLine(argc, argv))
        return EXIT_FAILURE;

    // INDEXING
    if (indexingMode)
    {
        if(bwt_index(refFile))
            return EXIT_FAILURE;
    }
    // SEARCHING
    else
    {
        Read *seqList;
        unsigned int seqListSize;
        uint32_t totalNumOfReads = 0;
        double startCpuTime = getCpuTime();
        double startRealTime = getRealTime();
        double ct;
        double rt;

        // Loading BWT-FM index
        if(bwt_load(refFile))
            return EXIT_FAILURE;

        // if (!initRead(seqFile, 500000000))
        if (!initRead(seqFile, 100000000))
            return EXIT_FAILURE;

        initializeFAST();

        while (readChunk(&seqList, &seqListSize) > 0)
        {
            totalNumOfReads += seqListSize;
            initFASTChunk(seqList, seqListSize);

            ct = getCpuTime();
            rt = getRealTime();
            fprintf(stderr, "\tmapping... ");
            
            mapSeqMT();

            fprintf(stderr, "done in %.2f seconds (%.2f CPU seconds)\n", getRealTime()-rt, getCpuTime()-ct);

            releaseChunk();
        }

        finalizeFAST();
        fprintf(stderr, "[NOTE] processed %u reads in %.2f seconds (%.2f CPU seconds)\n", totalNumOfReads, getRealTime()-startRealTime, getCpuTime()-startCpuTime);
    }
    return 0;
}
