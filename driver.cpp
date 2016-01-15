#include "align.h"
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <stdlib.h>
#include <cmath>

//#define MAX_SEQS 10000000
#define MAX_SEQS 10

double When() {
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}

int main(int argc, char *argv[])
{
	int NUM_SEQ = MAX_SEQS;
	if(argc > 1)
		NUM_SEQ = atoi(argv[1]);

	fprintf(stderr,"Starting align...\n");

	Align* align = new Align();

	//align->setGPU(0);
	align->initScoreMatrix(DEFAULT);
	align->setGapPenalty(-2.0f);

	//generate test data
	const int SEQ_LEN = 35;
	const int SEQ_LEN_B = 35;

	//const int NUM_SEQ = 135;
	const int BLOCK_SIZE = 128;

	const int num_blocks = (int)ceil(((float)NUM_SEQ)/BLOCK_SIZE);
	const unsigned long malloc_size = SEQ_LEN * num_blocks * BLOCK_SIZE;

	char* seqA = (char*)malloc(malloc_size);
	memset(seqA,'n',malloc_size);
	int randNum;
	srand(12345);
	fprintf(stderr,"Setting strings...");
	int blockNum = -1;
	for ( int i = 0; i < NUM_SEQ ; ++i ) {
		if(i % BLOCK_SIZE == 0)
			blockNum++;

		// Inner loop creates the sequence
		for ( int j = 0; j < SEQ_LEN ; ++j ) {
			randNum = rand() % 4;
			if ( randNum == 0 )
				seqA[(i%BLOCK_SIZE)+BLOCK_SIZE*j + blockNum*BLOCK_SIZE*SEQ_LEN] = 'a';
			else if ( randNum == 1 )
				seqA[(i%BLOCK_SIZE)+BLOCK_SIZE*j + blockNum*BLOCK_SIZE*SEQ_LEN] = 'g';
			else if ( randNum == 2 )
				seqA[(i%BLOCK_SIZE)+BLOCK_SIZE*j + blockNum*BLOCK_SIZE*SEQ_LEN] = 'c';
			else
				seqA[(i%BLOCK_SIZE)+BLOCK_SIZE*j + blockNum*BLOCK_SIZE*SEQ_LEN] = 't';
		}
	}
	fprintf(stderr,"Done!\n");

	/*
    char matchSeq[SEQ_LEN_B];
    for ( int i = 0; i < SEQ_LEN_B; ++i ) {
        matchSeq[i] = seqA[i * BLOCK_SIZE];
        seqA[i * BLOCK_SIZE + 1] = matchSeq[i];
        seqA[i * BLOCK_SIZE + 2] = matchSeq[i];
    }
	 */

	//tggtgagcaatgggcgcagagtactttcgtgagtt



	float* x = (float*)malloc(sizeof(float) * SEQ_LEN_B * 4);
	memset(x,0,sizeof(float)*SEQ_LEN_B*4);
	//	for(int i=0; i<SEQ_LEN_B; i++) {
	//		x[i] = (float*)malloc(sizeof(float) * 4);
	//		x[i][0] = 0.0f;
	//		x[i][1] = 0.0f;
	//		x[i][2] = 0.0f;
	//		x[i][3] = 0.0f;
	//	}
	//float x2[SEQ_LEN_B][4];
	int row = 0;
	//tggtgagcaatgggcgcagagtactttcgtgagtt
	x[row*4+0] = 0.0f; x[row*4+1] = 0.0f; x[row*4+2] = 0.0f; x[row*4+3] = 1.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 0.0f; x[row*4+2] = 1.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 0.0f; x[row*4+2] = 1.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 0.0f; x[row*4+2] = 0.0f; x[row*4+3] = 1.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 0.0f; x[row*4+2] = 1.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 1.0f; x[row*4+1] = 0.0f; x[row*4+2] = 0.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 0.0f; x[row*4+2] = 1.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 1.0f; x[row*4+2] = 0.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 1.0f; x[row*4+1] = 0.0f; x[row*4+2] = 0.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 1.0f; x[row*4+1] = 0.0f; x[row*4+2] = 0.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 0.0f; x[row*4+2] = 0.0f; x[row*4+3] = 1.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 0.0f; x[row*4+2] = 1.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 0.0f; x[row*4+2] = 1.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 0.0f; x[row*4+2] = 1.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 1.0f; x[row*4+2] = 0.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 0.0f; x[row*4+2] = 1.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 1.0f; x[row*4+2] = 0.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 1.0f; x[row*4+1] = 0.0f; x[row*4+2] = 0.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 0.0f; x[row*4+2] = 1.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 1.0f; x[row*4+1] = 0.0f; x[row*4+2] = 0.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 0.0f; x[row*4+2] = 1.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 0.0f; x[row*4+2] = 0.0f; x[row*4+3] = 1.0f; row++;
	x[row*4+0] = 1.0f; x[row*4+1] = 0.0f; x[row*4+2] = 0.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 1.0f; x[row*4+2] = 0.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 0.0f; x[row*4+2] = 0.0f; x[row*4+3] = 1.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 0.0f; x[row*4+2] = 0.0f; x[row*4+3] = 1.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 0.0f; x[row*4+2] = 0.0f; x[row*4+3] = 1.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 1.0f; x[row*4+2] = 0.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 0.0f; x[row*4+2] = 1.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 0.0f; x[row*4+2] = 0.0f; x[row*4+3] = 1.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 0.0f; x[row*4+2] = 1.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 1.0f; x[row*4+1] = 0.0f; x[row*4+2] = 0.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 0.0f; x[row*4+2] = 1.0f; x[row*4+3] = 0.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 0.0f; x[row*4+2] = 0.0f; x[row*4+3] = 1.0f; row++;
	x[row*4+0] = 0.0f; x[row*4+1] = 0.0f; x[row*4+2] = 0.0f; x[row*4+3] = 1.0f; row++;
	for(; row<SEQ_LEN_B; row++) {
		x[row*4+0] = 0.0f;
		x[row*4+1] = 0.0f;
		x[row*4+2] = 1.0f;
		x[row*4+3] = 0.0f;
	}



	char matchSeq[SEQ_LEN_B+1];
	//tggtgagcaatgggcgcagagtactttcgtgagtt
	matchSeq[0] = 't';
	matchSeq[1] = 'g';
	matchSeq[2] = 'g';
	matchSeq[3] = 't';
	matchSeq[4] = 'g';
	matchSeq[5] = 'a';
	matchSeq[6] = 'g';
	matchSeq[7] = 'c';
	matchSeq[8] = 'a';
	matchSeq[9] = 'a';
	matchSeq[10] = 't';
	matchSeq[11] = 'g';
	matchSeq[12] = 'g';
	matchSeq[13] = 'g';
	matchSeq[14] = 'c';
	matchSeq[15] = 'g';
	matchSeq[16] = 'c';
	matchSeq[17] = 'a';
	matchSeq[18] = 'g';
	matchSeq[19] = 'a';
	matchSeq[20] = 'g';
	matchSeq[21] = 't';
	matchSeq[22] = 'a';
	matchSeq[23] = 'c';
	matchSeq[24] = 't';
	matchSeq[25] = 't';
	matchSeq[26] = 't';
	matchSeq[27] = 'c';
	matchSeq[28] = 'g';
	matchSeq[29] = 't';
	matchSeq[30] = 'g';
	matchSeq[31] = 'a';
	matchSeq[32] = 'g';
	matchSeq[33] = 't';
	matchSeq[34] = 't';
	matchSeq[SEQ_LEN_B] = '\0';


	double start = When();
	//initialize data, compute, and gets results
#ifdef SEQ_STYLE
	align->loadMatchSequence((char*)&matchSeq, SEQ_LEN_B);
#else
	align->loadRead(x, SEQ_LEN_B);
#endif

	float* scores;

	// For testing times, let's see how long it takes to match 5M seqs in different sizes
	for(int i=0; i<MAX_SEQS; i+=NUM_SEQ) {
		align->loadSequences(seqA, SEQ_LEN, NUM_SEQ);


		align->computeAlignments();
		scores = align->getScores();
	}

	double end = When();
	fprintf(stderr,"Time=%.4f seconds\n",end-start);

#ifdef DEBUG
	for(int i=0; i<NUM_SEQ && i<BLOCK_SIZE; i++) {
		fprintf(stderr,"Score of %f and sequence\n%s\n",scores[i],matchSeq);

		for(int j=0; j<SEQ_LEN; j++) {
			fprintf(stderr,"%c",seqA[i+j*128]);
		}
		fprintf(stderr,"\n");
	}
#endif

	float denom = 0;
	double startD = When();
	align->computeDenominator(denom);
	fprintf(stderr,"Time to Denom=%f seconds and denominator is %f\n",When()-startD,denom);

	return 0;
}
