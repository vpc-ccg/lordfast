/*
 * Author: Ehsan Haghshenas (ehaghshe AT sfu DOT ca)
 */

#include "Chain.h"

#include <deque>

int par_maxSeedDist = 500;
const int8_t LogTable256[256] = {
    -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
     4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
     5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
     5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
     6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
     7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
     7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
     7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
     7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
     7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
     7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
     7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
     7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
};

inline int32_t ilog2_32(uint32_t v)
{
	uint32_t t, tt;
	if ((tt = v>>16)) return (t = tt>>8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
	return (t = v>>8) ? 8 + LogTable256[t] : LogTable256[v];
}

inline int dist_read(Seed_t &s1, Seed_t &s2)
{
	return s2.qPos - (s1.qPos + s1.len - 1);
}

inline double score_alpha(int distR, int distT, int seedLen)
{
    int minDist = distR < distT ? distR : distT;
    return seedLen < minDist ? seedLen : minDist;
}

inline double score_betha(int distR, int distT, int seedLen)
{
    uint32_t d = distR > distT ? distR - distT : distT - distR;
    int32_t dLog = d ? ilog2_32(d) : 0;
	// return 0.01 * seedLen * d + 0.5 * dLog;
	return 0.01 * 17 * d + 0.5 * dLog;
}

// void chain_seeds(vector<Seed_t> &seedList, Graph_t &graph, deque<int> &bestChain)
void chain_seeds(Seed_t *fragment_list, uint32_t nFragment, Chain_t &bestChain)
{
    // f(i) = max{max_{j<i}{f(i) + a(i, j) - b(i, j)}, w_i}
    // a(i, j) = min{min{y_i - y_j, x_i - x_j}, w_i}
    // b(i, j) = inf     y_j >= y_i || max{y_i - y_j, x_i - x_j} > maxDist
    // b(i, j) = gap_cost

	int i, j;
	// int n = seedList.size();
	int distR, distT;
	double aScore, bScore;
	double dp[nFragment + 1];
    int prev[nFragment + 1];
    
    double bestScore = -1;
    int bestIndex = -1;

    for(i = 0; i < nFragment; i++)
    {
		dp[i] = fragment_list[i].len;
		prev[i] = -1;
        for(j = i - 1; j >=0 ; j--)
        {
			distR = fragment_list[i].qPos - (fragment_list[j].qPos + fragment_list[j].len - 1);
			if(distR <= 0) continue;
			// if(distR < 0) continue;
			if(distR > par_maxSeedDist) break;
			
            distT = fragment_list[i].tPos - (fragment_list[j].tPos + fragment_list[j].len - 1);
            if(distT < 0)
            	continue;
            aScore = score_alpha(distR, distT, fragment_list[i].len);
            // aScore = fragment_list[i].len + 1;
            bScore = score_betha(distR, distT, fragment_list[i].len);

			if(dp[j] + aScore - bScore > dp[i])
			{
				dp[i] = dp[j] + aScore - bScore;
				prev[i] = j;
			}
        }
        // 
        if(dp[i] > bestScore)
        {
            bestScore = dp[i];
            bestIndex = i;
        }
    }

    // for(i = 0; i < nFragment; i++)
    // {
    //     fprintf(stdout, "### %-4d %-6d >%-6d %-6d %-4d %c score: %.2lf\n", i, seedList[i].readPos, seedList[i].nodeId, seedList[i].nodePos, seedList[i].len, (seedList[i].isRev ? '-' : '+'), dp[i]);
    //     // cout<< i << " " <<  " " << dp[i] << endl;
    // }

    // fprintf(stdout, "best score: %.2lf (%d)\n", bestScore, bestIndex);

    std::deque<int> ChainIndex;
    while(bestIndex != -1)
    {
        ChainIndex.push_front(bestIndex);
        bestIndex = prev[bestIndex];
    }

    bestChain.score = bestScore;
    bestChain.chainLen = ChainIndex.size();
    for(i = 0; i < ChainIndex.size(); i++)
    {
    	bestChain.seeds[i] = fragment_list[ChainIndex[i]];
    }
} 
