#include "align.h"
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <stdlib.h>
#include <math.h>

#define BLOCK_SIZE	128
#define MISMATCH	-1.0f
#define MATCH		1.0f

extern "C" bool setCudaDevice(int device);
extern "C" const char* initScoreMatrixTexture(float score_matrix[128][128]);
extern "C" const char* setMatchSequence(char* matchSeq, int matchLen);
extern "C" const char* setRead(float* matchSeq, int matchLen);
extern "C" const char* initSequences(int seqALen, int blockSize, int numBlocks);
extern "C" const char* setSequences(char* seqA, int start, int bytes);
extern "C" const char* copyScores(float* scores, int start, int bytes);
extern "C" const char* initRow(int blockSize, int numBlocks, int seqALen);
extern "C" const char* initScores(int blockSize, int numBlocks);
extern "C" long getTotalGlobalMemory();
extern "C" int getMaxGridBlocks();
extern "C" void runAlignKernel(int blocks, int blockSize, int seqALen, int seqBLen, float gap);
extern "C" void runSumKernel(int blocks, int blockSize, int smSize, float& sum);

Align::Align()
{
    scores = NULL;
}

bool setGPU(int num)
{
    return setCudaDevice(num);
}

double Align::When()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}

void Align::loadMatchSequence(char* seq, int len)
{
    seqB = seq;
    seqBLen = len;
    setMatchSequence(seq, len);
}

void Align::loadRead(float* seq, int len)
{
	seqBPtr = seq;
	seqBLen = len;
	setRead(seq, len);
}

void Align::loadSequences(char* seqs, int seqLen, int numSeq)
{
    seqA = seqs;
    seqALen = seqLen;
    this->numSeq = numSeq;

    long globalMem = getTotalGlobalMemory() * 0.75f;
    //printf("amount of global memory: %ld\n", globalMem);

    numBlocks = (int)ceil((float)numSeq / (float)BLOCK_SIZE);
    maxBlocks = globalMem / ((sizeof(float) + sizeof(char)) * seqALen + sizeof(float)) / BLOCK_SIZE;

    int deviceMaxGridBlocks = getMaxGridBlocks();
    if ( maxBlocks > deviceMaxGridBlocks )
        maxBlocks = deviceMaxGridBlocks;

    initSequences(seqLen, BLOCK_SIZE, maxBlocks);
}

void Align::setGapPenalty(float penalty)
{
    gap = penalty;
}

float* Align::getScores()
{
    return scores;
}

void Align::initScoreMatrix(ScoringType type)
{
    switch (type) {

    case BLOSUM62:
        break;

    case DEFAULT:
    default:
        score_matrix['A']['A'] = MATCH;
        score_matrix['A']['G'] = MISMATCH;
        score_matrix['A']['C'] = MISMATCH;
        score_matrix['A']['T'] = MISMATCH;
        score_matrix['G']['A'] = MISMATCH;
        score_matrix['G']['G'] = MATCH;
        score_matrix['G']['C'] = MISMATCH;
        score_matrix['G']['T'] = MISMATCH;
        score_matrix['C']['A'] = MISMATCH;
        score_matrix['C']['G'] = MISMATCH;
        score_matrix['C']['C'] = MATCH;
        score_matrix['C']['T'] = MISMATCH;
        score_matrix['T']['A'] = MISMATCH;
        score_matrix['T']['G'] = MISMATCH;
        score_matrix['T']['C'] = MISMATCH;
        score_matrix['T']['T'] = MATCH;

        score_matrix['a']['a'] = MATCH;
        score_matrix['a']['g'] = MISMATCH;
        score_matrix['a']['c'] = MISMATCH;
        score_matrix['a']['t'] = MISMATCH;
        score_matrix['g']['a'] = MISMATCH;
        score_matrix['g']['g'] = MATCH;
        score_matrix['g']['c'] = MISMATCH;
        score_matrix['g']['t'] = MISMATCH;
        score_matrix['c']['a'] = MISMATCH;
        score_matrix['c']['g'] = MISMATCH;
        score_matrix['c']['c'] = MATCH;
        score_matrix['c']['t'] = MISMATCH;
        score_matrix['t']['a'] = MISMATCH;
        score_matrix['t']['g'] = MISMATCH;
        score_matrix['t']['c'] = MISMATCH;
        score_matrix['t']['t'] = MATCH;

        // n's are always mismatches
//        score_matrix['n']['a'] = MISMATCH;
//        score_matrix['n']['c'] = MISMATCH;
//        score_matrix['n']['g'] = MISMATCH;
//        score_matrix['n']['t'] = MISMATCH;
//        score_matrix['N']['a'] = MISMATCH;
//        score_matrix['N']['c'] = MISMATCH;
//        score_matrix['N']['g'] = MISMATCH;
//        score_matrix['N']['t'] = MISMATCH;
        score_matrix['n']['a'] = 0;
        score_matrix['n']['c'] = 0;
        score_matrix['n']['g'] = 0;
        score_matrix['n']['t'] = 0;
        score_matrix['N']['a'] = 0;
        score_matrix['N']['c'] = 0;
        score_matrix['N']['g'] = 0;
        score_matrix['N']['t'] = 0;
        score_matrix['N']['N'] = 0;
        score_matrix['n']['n'] = 0;
        score_matrix['a']['n'] = 0;
        score_matrix['c']['n'] = 0;
        score_matrix['g']['n'] = 0;
        score_matrix['t']['n'] = 0;
        score_matrix['A']['n'] = 0;
        score_matrix['C']['n'] = 0;
        score_matrix['G']['n'] = 0;
        score_matrix['T']['n'] = 0;

        break;

    }
    initScoreMatrixTexture(score_matrix);
}

void Align::computeAlignments()
{
    scores = new float[numSeq];

    initRow(BLOCK_SIZE, maxBlocks, seqALen);
    initScores(BLOCK_SIZE, numBlocks);

    int numIterations = (int)ceil((float)numBlocks / (float)maxBlocks);
//    printf("max blocks: %d, %d\n", maxBlocks, numIterations);

    double t0 = When();

    for ( int i = 0; i < numIterations; ++i ) {
        if ( i == numIterations - 1 ) { //last iteration
            int remainingBlocks = numBlocks - (maxBlocks * i);
            int remainingSeqs = numSeq - (maxBlocks * i * BLOCK_SIZE);

            int nBlocks = (int)ceil((float)remainingSeqs/BLOCK_SIZE);

//            fprintf(stderr,"Setting sequences starting at %d, going for %d\n",maxBlocks * i * BLOCK_SIZE,nBlocks * BLOCK_SIZE * seqALen);
            setSequences(seqA, maxBlocks * i * BLOCK_SIZE, nBlocks * BLOCK_SIZE * seqALen);
            runAlignKernel(remainingBlocks, BLOCK_SIZE, seqALen, seqBLen, gap);
            copyScores(scores, maxBlocks * i * BLOCK_SIZE, sizeof(float) * remainingSeqs);
//            fprintf(stderr,"Calling copyScores for the last time.  First is %f\n",scores[0]);
        }
        else {
//        	fprintf(stderr,"Setting sequences starting at %d, going for %d\n",maxBlocks * i * BLOCK_SIZE,maxBlocks * BLOCK_SIZE * seqALen);
            setSequences(seqA, maxBlocks * i * BLOCK_SIZE, maxBlocks * BLOCK_SIZE * seqALen);
            runAlignKernel(maxBlocks, BLOCK_SIZE, seqALen, seqBLen, gap);
            copyScores(scores, maxBlocks * i * BLOCK_SIZE, sizeof(float) * maxBlocks * BLOCK_SIZE);
//            fprintf(stderr,"Calling copyScores again.  First is %f\n",scores[0]);
        }
    }

    double t = When() - t0;

//    printf("\n ----------------------------------------------------------------------\n ");
//    printf( "Computed in %lf seconds", t);
//    printf("\n ----------------------------------------------------------------------\n");

}

void Align::computeDenominator(float& denom)
{
	int numIterations = (int)ceil((float)numBlocks / (float)maxBlocks);


	for ( int i = 0; i < numIterations; ++i ) {
		if ( i == numIterations - 1 ) { //last iteration
			int remainingBlocks = numBlocks - (maxBlocks * i);
			int remainingSeqs = numSeq - (maxBlocks * i * BLOCK_SIZE);

			int nBlocks = (int)ceil((float)remainingSeqs/BLOCK_SIZE);

//			setSequences(seqA, maxBlocks * i * BLOCK_SIZE, nBlocks * BLOCK_SIZE * seqALen);
//			runAlignKernel(remainingBlocks, BLOCK_SIZE, seqALen, seqBLen, gap);
//			copyScores(scores, maxBlocks * i * BLOCK_SIZE, sizeof(float) * remainingSeqs);
			fprintf(stderr,"Running SUM KERNAL, nBlocks=%d, remainingBlocks=%d, remainingSeqs=%d\n",nBlocks,remainingBlocks,remainingSeqs);
			//runSumKernel(remainingBlocks, BLOCK_SIZE, BLOCK_SIZE, denom);
			//runSumKernel(remainingBlocks, BLOCK_SIZE, BLOCK_SIZE*sizeof(float), denom);

			runSumKernel(remainingBlocks, BLOCK_SIZE, nBlocks*BLOCK_SIZE*sizeof(float), denom);
		}
		else {
//			setSequences(seqA, maxBlocks * i * BLOCK_SIZE, maxBlocks * BLOCK_SIZE * seqALen);
//			runAlignKernel(maxBlocks, BLOCK_SIZE, seqALen, seqBLen, gap);
//			copyScores(scores, maxBlocks * i * BLOCK_SIZE, sizeof(float) * maxBlocks * BLOCK_SIZE);
			runSumKernel(maxBlocks, BLOCK_SIZE, BLOCK_SIZE*sizeof(float), denom);
		}
	}

}


