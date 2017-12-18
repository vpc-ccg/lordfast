#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/times.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "info.h"
#include "debug.h"
#include "container.h"
#include "manopt.h"
#include "fileio.h"
#include "sltypes.h"
#include "slchain.h"
#include "rangetree.h"

#ifdef __cplusplus
}
#endif

#include <algorithm>

#include "Clasp.h"

unsigned char mute = 0;
extern double maxmem;

int clasp_chain_seed_best(Seed_t *fragment_list, uint32_t nFragment, Chain_t &bestChain)
{
	Container *fragments;
	// claspinfo_t info;
	manopt_optionset optset;
	manopt_arg *unflagged;
	manopt_arg *list;
	Uint i, j, k, begin;
	int num;
	time_t start, end;
	double chtime, intime;
	char chainmode = SOP;
	// char chainmode = LIN;
	// double lambda = 0.5;
	double lambda = 0.15;
	double epsilon = 0;
	double minscore = 0;
	int maxgap = -1;
	int minfrag = 0;
	BOOL outputm = 1;
	BOOL outputc = 1;
	BOOL outputf = 1;

	bestChain.score = -1;

	/* initialization */
	fragments = (Container *) malloc(sizeof(Container));
	bl_containerInit(fragments, 1000, sizeof(slmatch_t));
	/* put the fragments in the container */
	slmatch_t frag;

	for (i=0; i < nFragment; i++)
	{
		bl_slmatchInit(&frag, 0);
		frag.p = fragment_list[i].tPos;
		frag.i = fragment_list[i].qPos;
		frag.q = frag.j = fragment_list[i].len;
		frag.scr = fragment_list[i].len;
		bl_containerAdd(fragments, &frag);
	}

	// for (i = 0; i < bl_containerSize(fragments); i++)
	// {
	// 	slmatch_t *match_tmp = (slmatch_t *) bl_containerGet(fragments, i);
	// 	fprintf(stderr, "%d\t%d\t%ld\t%ld\t%f\n", 
	// 		match_tmp->i, match_tmp->j, 
	// 		match_tmp->p, match_tmp->q,
	// 		match_tmp->scr);
	// }

	// exit(0);

	/* sort fragments */
	qsort(fragments->contspace, bl_containerSize(fragments),
	sizeof(slmatch_t), cmp_slmatch_qsort);
	begin = 0;
	for (i = 1; i <= bl_containerSize(fragments); i++)
	{
		/* 
		 * end of fragments list or different database sequence 
		 * --> process fragment[begin]...fragment[i-1], write output
		 *     and free chains (less memory consumption with large input files)
		 */
		if (i == bl_containerSize(fragments) ||
		((slmatch_t *) bl_containerGet(fragments, begin))->subject !=
		((slmatch_t *) bl_containerGet(fragments, i))->subject)
		{
			//fprintf(info.dev, "%d\t%d\n", begin, i-begin);
			if (chainmode == SOP)
			{
				/* only use chaining without clustering if no ids are specified */
				//bl_slChainSop((slmatch_t *) info.fragments->contspace + begin, i - begin,
				//	      info.epsilon, info.lambda);  
				bl_slClusterSop((slmatch_t *) fragments->contspace + begin, i - begin,
						epsilon, lambda, maxgap);
			}
			else
			{    
				//bl_slChainLin((slmatch_t *) info.fragments->contspace + begin, i - begin,
				//	      info.epsilon, info.lambda);
				bl_slClusterLin((slmatch_t *) fragments->contspace + begin, i - begin,
						epsilon, lambda, maxgap);
			}

			for (j = begin; j < i; j++)
			{
				slmatch_t *match = (slmatch_t *) bl_containerGet(fragments, j);

				/* output matches (if desired) */
				// if (outputm)
				// {
				// 	fprintf(stderr, "M\t");
				// 	fprintf(stderr, "%d\t%d\t%ld\t%ld\t%.3f\n", match->i,
				// 	    match->i + match->j - 1, match->p,
				// 	    match->p + match->q - 1, match->scr);
				// }
				if (match->chain)
				{
					slchain_t *chain = (slchain_t *) match->chain;

					if(chain->scr > bestChain.score)
					{
						bestChain.score = chain->scr;
						bestChain.chainLen = 0;
						// chain fragmenst
						for (k = 0; k < bl_containerSize(chain->matches); k++)
						{
							slmatch_t *frag = *(slmatch_t **) bl_containerGet(chain->matches, k);
							bestChain.seeds[bestChain.chainLen].tPos = frag->p;
							bestChain.seeds[bestChain.chainLen].qPos = frag->i;
							bestChain.seeds[bestChain.chainLen++].len = frag->j; // same as frag->q
						}
					}

					// if (outputc && chain->scr >= minscore &&
					// 	bl_containerSize(chain->matches) >= minfrag)
					// {
					// 	fprintf(stderr, "C\t");
					// 	fprintf(stderr, "%d\t%d\t%ld\t%ld\t%.3f\n", chain->i,
					// 	    chain->i + chain->j - 1, chain->p,
					// 	    chain->p + chain->q - 1, chain->scr);
					// }
					// /* output chains and fragments (if requested) */
					// if (outputf && chain->scr >= minscore &&
					// 	bl_containerSize(chain->matches) >= minfrag)
					// {
					// 	for (k = 0; k < bl_containerSize(chain->matches); k++)
					// 	{
					// 		slmatch_t *frag = *(slmatch_t **)
					// 		bl_containerGet(chain->matches, k);
					// 		fprintf(stderr, "F\t");
					// 		fprintf(stderr, "%d\t%d\t%ld\t%ld\t%.3f\n", frag->i,
					// 		frag->i + frag->j - 1, frag->p, frag->p + frag->q - 1,
					// 		frag->scr);
					// 	}
					// }
					bl_slchainDestruct(chain);
					free(chain);
					match->chain = NULL;
				}
			}
			begin = i;
		}
	}

	if (fragments)
	{
		for (i = 0; i < bl_containerSize(fragments); i++)
		{
			slmatch_t *sl = (slmatch_t *) bl_containerGet(fragments, i);
			if (sl->chain != NULL)
			{
				DBG("still at least one chain not freed before end: %d", i);
				exit(-1);
				bl_slchainDestruct(sl->chain);
				free(sl->chain);
			}
		}
		bl_containerDestruct(fragments, bl_slmatchDestruct);
		free(fragments);
	}

	// for (i = 0; i < (*chainNum); i++)
	// {
	// 	fprintf(stderr, "+++\t%u\t%u\t%u\t%u\t%f\n", chainList[i].qStart, chainList[i].qEnd, 
	// 		chainList[i].tStart, chainList[i].tEnd, chainList[i].score);
	// }

	// if(isReverse)
	// 	std::sort_heap (chainList, chainList+(*chainNum), compare_chain);

	return 1;
}

// int clasp_chain_seed_top(Seed_t *fragment_list, uint32_t nFragment, Chain_t *chainList, uint8_t *chainNum, int maxChain, int isReverse)
// {
// 	Container *fragments;
// 	// claspinfo_t info;
// 	manopt_optionset optset;
// 	manopt_arg *unflagged;
// 	manopt_arg *list;
// 	Uint i, j, k, begin;
// 	int num;
// 	time_t start, end;
// 	double chtime, intime;
// 	char chainmode = SOP;
// 	// char chainmode = LIN;
// 	// double lambda = 0.5;
// 	double lambda = 0.2;
// 	double epsilon = 0;
// 	double minscore = 0;
// 	int maxgap = -1;
// 	int minfrag = 0;
// 	BOOL outputm = 0;
// 	BOOL outputc = 1;
// 	BOOL outputf = 0;

// 	/* initialization */
// 	fragments = (Container *) malloc(sizeof(Container));
// 	bl_containerInit(fragments, 1000, sizeof(slmatch_t));
// 	/* put the fragments in the container */
// 	slmatch_t frag;

// 	for (i=0; i < nFragment; i++)
// 	{
// 		bl_slmatchInit(&frag, 0);
// 		frag.p = fragment_list[i].tPos;
// 		frag.i = fragment_list[i].qPos;
// 		frag.q = frag.j = fragment_list[i].len;
// 		frag.scr = fragment_list[i].len;
// 		bl_containerAdd(fragments, &frag);
// 	}

// 	// for (i = 0; i < bl_containerSize(fragments); i++)
// 	// {
// 	// 	slmatch_t *match_tmp = (slmatch_t *) bl_containerGet(fragments, i);
// 	// 	fprintf(stderr, "%d\t%d\t%ld\t%ld\t%f\n", 
// 	// 		match_tmp->i, match_tmp->j, 
// 	// 		match_tmp->p, match_tmp->q,
// 	// 		match_tmp->scr);
// 	// }

// 	// exit(0);

// 	/* sort fragments */
// 	qsort(fragments->contspace, bl_containerSize(fragments),
// 	sizeof(slmatch_t), cmp_slmatch_qsort);
// 	begin = 0;
// 	for (i = 1; i <= bl_containerSize(fragments); i++)
// 	{
// 		/* 
// 		 * end of fragments list or different database sequence 
// 		 * --> process fragment[begin]...fragment[i-1], write output
// 		 *     and free chains (less memory consumption with large input files)
// 		 */
// 		if (i == bl_containerSize(fragments) ||
// 		((slmatch_t *) bl_containerGet(fragments, begin))->subject !=
// 		((slmatch_t *) bl_containerGet(fragments, i))->subject)
// 		{
// 			//fprintf(info.dev, "%d\t%d\n", begin, i-begin);
// 			if (chainmode == SOP)
// 			{
// 				/* only use chaining without clustering if no ids are specified */
// 				//bl_slChainSop((slmatch_t *) info.fragments->contspace + begin, i - begin,
// 				//	      info.epsilon, info.lambda);  
// 				bl_slClusterSop((slmatch_t *) fragments->contspace + begin, i - begin,
// 						epsilon, lambda, maxgap);
// 			}
// 			else
// 			{    
// 				//bl_slChainLin((slmatch_t *) info.fragments->contspace + begin, i - begin,
// 				//	      info.epsilon, info.lambda);
// 				bl_slClusterLin((slmatch_t *) fragments->contspace + begin, i - begin,
// 						epsilon, lambda, maxgap);
// 			}

// 			for (j = begin; j < i; j++)
// 			{
// 				slmatch_t *match = (slmatch_t *) bl_containerGet(fragments, j);

// 				/* output matches (if desired) */
// 				if (outputm)
// 				{
// 					fprintf(stderr, "M\t");
// 					fprintf(stderr, "%d\t%d\t%ld\t%ld\t%.3f\n", match->i,
// 					    match->i + match->j - 1, match->p,
// 					    match->p + match->q - 1, match->scr);
// 				}
// 				if (match->chain)
// 				{
// 					slchain_t *chain = (slchain_t *) match->chain;
// 					if (outputc && chain->scr >= minscore &&
// 						bl_containerSize(chain->matches) >= minfrag)
// 					{
// 						if((*chainNum) < maxChain) // the list has some space, push the chains 
// 						{
// 							chainList[(*chainNum)].qStart = chain->i;
// 							chainList[(*chainNum)].qEnd = chain->i + chain->j - 1;
// 							chainList[(*chainNum)].tStart = chain->p;
// 							chainList[(*chainNum)].tEnd = chain->p + chain->q - 1;
// 							chainList[(*chainNum)].score = chain->scr;
// 							chainList[(*chainNum)++].strand = isReverse;
// 							std::push_heap(chainList, chainList+(*chainNum), compare_chain);
// 						}
// 						else // the list is full
// 						{
// 							if(chain->scr > chainList[0].score) // if the score is higher than the current smallest score, replace it with that
// 							{
//   								std::pop_heap(chainList, chainList+(*chainNum), compare_chain); // the smallest score chain is at the last position
//   								// update the last chain
// 								chainList[(*chainNum)-1].qStart = chain->i;
// 								chainList[(*chainNum)-1].qEnd = chain->i + chain->j - 1;
// 								chainList[(*chainNum)-1].tStart = chain->p;
// 								chainList[(*chainNum)-1].tEnd = chain->p + chain->q - 1;
// 								chainList[(*chainNum)-1].score = chain->scr;
// 								chainList[(*chainNum)-1].strand = isReverse;
// 								std::push_heap(chainList, chainList+(*chainNum), compare_chain);
// 							}
// 						}

// 						// fprintf(stderr, "C\t");
// 						// fprintf(stderr, "%d\t%d\t%ld\t%ld\t%.3f\n", chain->i,
// 						//     chain->i + chain->j - 1, chain->p,
// 						//     chain->p + chain->q - 1, chain->scr);
// 					}
// 					/* output chains and fragments (if requested) */
// 					if (outputf && chain->scr >= minscore &&
// 						bl_containerSize(chain->matches) >= minfrag)
// 					{
// 						for (k = 0; k < bl_containerSize(chain->matches); k++)
// 						{
// 							slmatch_t *frag = *(slmatch_t **)
// 							bl_containerGet(chain->matches, k);
// 							fprintf(stderr, "F\t");
// 							fprintf(stderr, "%d\t%d\t%ld\t%ld\t%.3f\n", frag->i,
// 							frag->i + frag->j - 1, frag->p, frag->p + frag->q - 1,
// 							frag->scr);
// 						}
// 					}
// 					bl_slchainDestruct(chain);
// 					free(chain);
// 					match->chain = NULL;
// 				}
// 			}
// 			begin = i;
// 		}
// 	}

// 	if (fragments)
// 	{
// 		for (i = 0; i < bl_containerSize(fragments); i++)
// 		{
// 			slmatch_t *sl = (slmatch_t *) bl_containerGet(fragments, i);
// 			if (sl->chain != NULL)
// 			{
// 				DBG("still at least one chain not freed before end: %d", i);
// 				exit(-1);
// 				bl_slchainDestruct(sl->chain);
// 				free(sl->chain);
// 			}
// 		}
// 		bl_containerDestruct(fragments, bl_slmatchDestruct);
// 		free(fragments);
// 	}

// 	// for (i = 0; i < (*chainNum); i++)
// 	// {
// 	// 	fprintf(stderr, "+++\t%u\t%u\t%u\t%u\t%f\n", chainList[i].qStart, chainList[i].qEnd, 
// 	// 		chainList[i].tStart, chainList[i].tEnd, chainList[i].score);
// 	// }

// 	if(isReverse)
// 		std::sort_heap (chainList, chainList+(*chainNum), compare_chain);

// 	return 1;
// }
