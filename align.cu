#include<stdio.h>
#include<string.h>
#include<sys/time.h>
#include<cuda.h>

#define BLOCK_SIZE 128

#ifdef SEQ_STYLE
__constant__ char seqB[64000];
#else
//__constant__ float seqB[128][4];
__constant__ float seqB[1600][4];
#endif

texture<float, 2, cudaReadModeElementType> tex;

//__device__ float sdata[256];

char* d_seqA;
float* d_scores;
float* d_row;
int currentDeviceNum = 0;

__device__ void calc(char& a, float* b, float& leftVal, float& diag, float& above, float& max, float& GAP)
{
#ifdef DEBUG
	printf("%.2f %.2f %.2f %.2f against %c\n",b[0],b[1],b[2],b[3],a);
#endif
	// Gets the alignment score from the "texture" (score matrix) with a and b
	max = tex2D(tex, 'a', a) * b[0];
	max += tex2D(tex, 'c', a) * b[1];
	max += tex2D(tex, 'g', a) * b[2];
	max += tex2D(tex, 't', a) * b[3];


    float score2 = leftVal + GAP;

    max += diag;
    float score3 = above + GAP;

    if ( score2 > max )
    {
        max = score2;
    }
    if ( score3 > max )
    {
        max = score3;
    }
}

__device__ void calc(char& a, char& b, float& leftVal, float& diag, float& above, float& max, float& GAP)
{

	// Gets the alignment score from the "texture" (score matrix) with a and b
	max = tex2D(tex, a, b);

    float score2 = leftVal + GAP;

    max += diag;
    float score3 = above + GAP;

    if ( score2 > max )
    {
        max = score2;
    }
    if ( score3 > max )
    {
        max = score3;
    }

//    printf("Comparing %c with %c, left=%f diag=%f above=%f max=%f gap=%f\n",a,b,leftVal,diag,above,max,GAP);

}

__global__ void align(char* dSeqA, float* dScores, float* dRow, int seqALen, int seqBLen, float GAP)
{
    int tid = threadIdx.x;
    int blockRowStart = (blockDim.x * seqBLen) * blockIdx.x;
    int blockSeqAStart = (blockDim.x * seqALen) * blockIdx.x;
    float max;
#ifdef DEBUG
    if(blockIdx.x != 0)
    	printf("[%d] blockDim.x=%d, blockIdx.x=%d, Aligning seqA=%s\n",tid,blockDim.x,blockIdx.x, dSeqA);
#endif

    float v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15;

    max = GAP;
    int i;
    for ( i = 0; i < seqBLen; ++i ) {
        dRow[blockRowStart + i * blockDim.x + tid] = max;
        max += GAP;
    }

    float leftVal, above, diag;

    char a;

    int remaining = seqBLen;
    float loop = 0.0f;
    float oldLeft = 0;
    int idx;
    while ( remaining > 15 )
    {
        float addGap = (GAP * loop * 16.0f);
        v0 = GAP + addGap;
        v1 = GAP * 2.0f + addGap;
        v2 = GAP * 3.0f + addGap;
        v3 = GAP * 4.0f + addGap;
        v4 = GAP * 5.0f + addGap;
        v5 = GAP * 6.0f + addGap;
        v6 = GAP * 7.0f + addGap;
        v7 = GAP * 8.0f + addGap;
        v8 = GAP * 9.0f + addGap;
        v9 = GAP * 10.0f + addGap;
        v10 = GAP * 11.0f + addGap;
        v11 = GAP * 12.0f + addGap;
        v12 = GAP * 13.0f + addGap;
        v13 = GAP * 14.0f + addGap;
        v14 = GAP * 15.0f + addGap;
        v15 = GAP * 16.0f + addGap;

        for ( int row = 0; row < seqALen; ++row ) {
            idx = (int)(loop * 16.0f);
            above = v0;
            if ( loop > 0 )
            //if ( row > 0 )
                leftVal = dRow[blockRowStart + row * blockDim.x + tid];
            else
                leftVal = (row+1) * GAP;

            a = dSeqA[tid + row * blockDim.x + blockSeqAStart];
#ifdef DEBUG
            if(blockIdx.x != 0)
            	printf("at:%d(%c=%d) ",tid + row * blockDim.x + blockSeqAStart,a,a);
#endif

            calc(a, seqB[idx++], leftVal, oldLeft, above, max, GAP);
            oldLeft = leftVal;
            leftVal = max;

            diag = above;
            v0 = max;

            above = v1;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            leftVal = max;
            diag = above;
            v1 = max;

            above = v2;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            leftVal = max;
            diag = above;
            v2 = max;


            above = v3;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            leftVal = max;
            diag = above;
            v3 = max;

            above = v4;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            leftVal = max;
            diag = above;
            v4 = max;

            above = v5;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            leftVal = max;
            diag = above;
            v5 = max;

            above = v6;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            leftVal = max;
            diag = above;
            v6 = max;

            above = v7;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            leftVal = max;
            diag = above;
            v7 = max;

            above = v8;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            leftVal = max;
            diag = above;
            v8 = max;

            above = v9;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            leftVal = max;
            diag = above;
            v9 = max;

            above = v10;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            leftVal = max;
            diag = above;
            v10 = max;

            above = v11;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            leftVal = max;
            diag = above;
            v11 = max;

            above = v12;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            leftVal = max;
            diag = above;
            v12 = max;

            above = v13;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            leftVal = max;
            diag = above;
            v13 = max;

            above = v14;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            leftVal = max;
            diag = above;
            v14 = max;

            above = v15;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            leftVal = max;
            diag = above;
            v15 = max;

            dRow[blockRowStart + row * blockDim.x + tid] = max;

        }
        remaining -= 16;
        loop += 1.0f;
    }
    if ( remaining == 0 ) {
        dScores[blockDim.x * blockIdx.x + tid] = max;
        return;
    }


    if ( remaining == 1 ) {
        float addGap = (GAP * loop * 16.0f);
        v0 = GAP + addGap;

        idx = (int)(loop * 16.0f);
        for ( int row = 0; row < seqALen; ++row ) {

            above = v0;
            leftVal = dRow[blockRowStart + row * blockDim.x + tid];

            a = dSeqA[tid + row * blockDim.x + blockSeqAStart];

            calc(a, seqB[idx], leftVal, oldLeft, above, max, GAP);
            oldLeft = leftVal;
            v0 = max;
        }

        dScores[blockDim.x * blockIdx.x + tid] = max;
        return;
    }


    else if ( remaining == 2 ) {
        float addGap = (GAP * loop * 16.0f);
        v0 = GAP + addGap;
        v1 = GAP * 2.0f + addGap;

        for ( int row = 0; row < seqALen; ++row ) {
            idx = (int)(loop * 16.0f);
            above = v0;
            leftVal = dRow[blockRowStart + row * blockDim.x + tid];

            a = dSeqA[tid + row * blockDim.x + blockSeqAStart];

            calc(a, seqB[idx++], leftVal, oldLeft, above, max, GAP);
            oldLeft = leftVal;
            v0 = max;

            leftVal = max;
            diag = above;
            above = v1;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v1 = max;
        }
        dScores[blockDim.x * blockIdx.x + tid] = max;
        return;
    }


    else if ( remaining == 3 ) {
        float addGap = (GAP * loop * 16.0f);
        v0 = GAP + addGap;
        v1 = GAP * 2.0f + addGap;
        v2 = GAP * 3.0f + addGap;

        for ( int row = 0; row < seqALen; ++row ) {
            idx = (int)(loop * 16.0f);
            above = v0;
            leftVal = dRow[blockRowStart + row * blockDim.x + tid];

            a = dSeqA[tid + row * blockDim.x + blockSeqAStart];

            calc(a, seqB[idx++], leftVal, oldLeft, above, max, GAP);
            oldLeft = leftVal;
            v0 = max;

            leftVal = max;
            diag = above;
            above = v1;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v1 = max;

            leftVal = max;
            diag = above;
            above = v2;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v2 = max;
        }
        dScores[blockDim.x * blockIdx.x + tid] = max;
        return;

    }

    else if ( remaining == 4 ) {
        float addGap = (GAP * loop * 16.0f);
        v0 = GAP + addGap;
        v1 = GAP * 2.0f + addGap;
        v2 = GAP * 3.0f + addGap;
        v3 = GAP * 4.0f + addGap;

        for ( int row = 0; row < seqALen; ++row ) {
            idx = (int)(loop * 16.0f);
            above = v0;
            leftVal = dRow[blockRowStart + row * blockDim.x + tid];

            a = dSeqA[tid + row * blockDim.x + blockSeqAStart];

            calc(a, seqB[idx++], leftVal, oldLeft, above, max, GAP);
            oldLeft = leftVal;
            v0 = max;

            leftVal = max;
            diag = above;
            above = v1;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v1 = max;

            leftVal = max;
            diag = above;
            above = v2;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v2 = max;

            leftVal = max;
            diag = above;
            above = v3;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v3 = max;
        }
        dScores[blockDim.x * blockIdx.x + tid] = max;
        return;

    }

    else if ( remaining == 5 ) {
        float addGap = (GAP * loop * 16.0f);
        v0 = GAP + addGap;
        v1 = GAP * 2.0f + addGap;
        v2 = GAP * 3.0f + addGap;
        v3 = GAP * 4.0f + addGap;
        v4 = GAP * 5.0f + addGap;

        for ( int row = 0; row < seqALen; ++row ) {
            idx = (int)(loop * 16.0f);
            above = v0;
            leftVal = dRow[blockRowStart + row * blockDim.x + tid];

            a = dSeqA[tid + row * blockDim.x + blockSeqAStart];

            calc(a, seqB[idx++], leftVal, oldLeft, above, max, GAP);
            oldLeft = leftVal;
            v0 = max;

            leftVal = max;
            diag = above;
            above = v1;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v1 = max;

            leftVal = max;
            diag = above;
            above = v2;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v2 = max;

            leftVal = max;
            diag = above;
            above = v3;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v3 = max;

            leftVal = max;
            diag = above;
            above = v4;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v4 = max;
        }
        dScores[blockDim.x * blockIdx.x + tid] = max;
        return;

    }


    else if ( remaining == 6 ) {
        float addGap = (GAP * loop * 16.0f);
        v0 = GAP + addGap;
        v1 = GAP * 2.0f + addGap;
        v2 = GAP * 3.0f + addGap;
        v3 = GAP * 4.0f + addGap;
        v4 = GAP * 5.0f + addGap;
        v5 = GAP * 6.0f + addGap;

        for ( int row = 0; row < seqALen; ++row ) {
            idx = (int)(loop * 16.0f);
            above = v0;
            leftVal = dRow[blockRowStart + row * blockDim.x + tid];

            a = dSeqA[tid + row * blockDim.x + blockSeqAStart];

            calc(a, seqB[idx++], leftVal, oldLeft, above, max, GAP);
            oldLeft = leftVal;
            v0 = max;

            leftVal = max;
            diag = above;
            above = v1;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v1 = max;

            leftVal = max;
            diag = above;
            above = v2;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v2 = max;

            leftVal = max;
            diag = above;
            above = v3;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v3 = max;

            leftVal = max;
            diag = above;
            above = v4;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v4 = max;

            leftVal = max;
            diag = above;
            above = v5;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v5 = max;
        }
        dScores[blockDim.x * blockIdx.x + tid] = max;
        return;

    }


    else if ( remaining == 7 ) {
        float addGap = (GAP * loop * 16.0f);
        v0 = GAP + addGap;
        v1 = GAP * 2.0f + addGap;
        v2 = GAP * 3.0f + addGap;
        v3 = GAP * 4.0f + addGap;
        v4 = GAP * 5.0f + addGap;
        v5 = GAP * 6.0f + addGap;
        v6 = GAP * 7.0f + addGap;

        for ( int row = 0; row < seqALen; ++row ) {
            idx = (int)(loop * 16.0f);
            above = v0;
            leftVal = dRow[blockRowStart + row * blockDim.x + tid];

            a = dSeqA[tid + row * blockDim.x + blockSeqAStart];

            calc(a, seqB[idx++], leftVal, oldLeft, above, max, GAP);
            oldLeft = leftVal;
            v0 = max;

            leftVal = max;
            diag = above;
            above = v1;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v1 = max;

            leftVal = max;
            diag = above;
            above = v2;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v2 = max;

            leftVal = max;
            diag = above;
            above = v3;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v3 = max;

            leftVal = max;
            diag = above;
            above = v4;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v4 = max;

            leftVal = max;
            diag = above;
            above = v5;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v5 = max;

            leftVal = max;
            diag = above;
            above = v6;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v6 = max;
        }
        dScores[blockDim.x * blockIdx.x + tid] = max;
        return;

    }


    else if ( remaining == 8 ) {
        float addGap = (GAP * loop * 16.0f);
        v0 = GAP + addGap;
        v1 = GAP * 2.0f + addGap;
        v2 = GAP * 3.0f + addGap;
        v3 = GAP * 4.0f + addGap;
        v4 = GAP * 5.0f + addGap;
        v5 = GAP * 6.0f + addGap;
        v6 = GAP * 7.0f + addGap;
        v7 = GAP * 8.0f + addGap;

        for ( int row = 0; row < seqALen; ++row ) {
            idx = (int)(loop * 16.0f);
            above = v0;
            leftVal = dRow[blockRowStart + row * blockDim.x + tid];

            a = dSeqA[tid + row * blockDim.x + blockSeqAStart];

            calc(a, seqB[idx++], leftVal, oldLeft, above, max, GAP);
            oldLeft = leftVal;
            v0 = max;

            leftVal = max;
            diag = above;
            above = v1;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v1 = max;

            leftVal = max;
            diag = above;
            above = v2;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v2 = max;

            leftVal = max;
            diag = above;
            above = v3;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v3 = max;

            leftVal = max;
            diag = above;
            above = v4;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v4 = max;

            leftVal = max;
            diag = above;
            above = v5;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v5 = max;

            leftVal = max;
            diag = above;
            above = v6;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v6 = max;

            leftVal = max;
            diag = above;
            above = v7;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v7 = max;
        }
        dScores[blockDim.x * blockIdx.x + tid] = max;
        return;

    }

    else if ( remaining == 9 ) {
        float addGap = (GAP * loop * 16.0f);
        v0 = GAP + addGap;
        v1 = GAP * 2.0f + addGap;
        v2 = GAP * 3.0f + addGap;
        v3 = GAP * 4.0f + addGap;
        v4 = GAP * 5.0f + addGap;
        v5 = GAP * 6.0f + addGap;
        v6 = GAP * 7.0f + addGap;
        v7 = GAP * 8.0f + addGap;
        v8 = GAP * 9.0f + addGap;

        for ( int row = 0; row < seqALen; ++row ) {
            idx = (int)(loop * 16.0f);
            above = v0;
            leftVal = dRow[blockRowStart + row * blockDim.x + tid];

            a = dSeqA[tid + row * blockDim.x + blockSeqAStart];

            calc(a, seqB[idx++], leftVal, oldLeft, above, max, GAP);
            oldLeft = leftVal;
            v0 = max;

            leftVal = max;
            diag = above;
            above = v1;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v1 = max;

            leftVal = max;
            diag = above;
            above = v2;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v2 = max;

            leftVal = max;
            diag = above;
            above = v3;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v3 = max;

            leftVal = max;
            diag = above;
            above = v4;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v4 = max;

            leftVal = max;
            diag = above;
            above = v5;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v5 = max;

            leftVal = max;
            diag = above;
            above = v6;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v6 = max;

            leftVal = max;
            diag = above;
            above = v7;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v7 = max;

            leftVal = max;
            diag = above;
            above = v8;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v8 = max;
        }
        dScores[blockDim.x * blockIdx.x + tid] = max;
        return;

    }

    else if ( remaining == 10 ) {
        float addGap = (GAP * loop * 16.0f);
        v0 = GAP + addGap;
        v1 = GAP * 2.0f + addGap;
        v2 = GAP * 3.0f + addGap;
        v3 = GAP * 4.0f + addGap;
        v4 = GAP * 5.0f + addGap;
        v5 = GAP * 6.0f + addGap;
        v6 = GAP * 7.0f + addGap;
        v7 = GAP * 8.0f + addGap;
        v8 = GAP * 9.0f + addGap;
        v9 = GAP * 10.0f + addGap;

        for ( int row = 0; row < seqALen; ++row ) {
            idx = (int)(loop * 16.0f);
            above = v0;
            leftVal = dRow[blockRowStart + row * blockDim.x + tid];

            a = dSeqA[tid + row * blockDim.x + blockSeqAStart];

            calc(a, seqB[idx++], leftVal, oldLeft, above, max, GAP);
            oldLeft = leftVal;
            v0 = max;

            leftVal = max;
            diag = above;
            above = v1;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v1 = max;

            leftVal = max;
            diag = above;
            above = v2;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v2 = max;

            leftVal = max;
            diag = above;
            above = v3;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v3 = max;

            leftVal = max;
            diag = above;
            above = v4;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v4 = max;

            leftVal = max;
            diag = above;
            above = v5;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v5 = max;

            leftVal = max;
            diag = above;
            above = v6;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v6 = max;

            leftVal = max;
            diag = above;
            above = v7;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v7 = max;

            leftVal = max;
            diag = above;
            above = v8;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v8 = max;

            leftVal = max;
            diag = above;
            above = v9;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v9 = max;
        }
        dScores[blockDim.x * blockIdx.x + tid] = max;
        return;

    }


    else if ( remaining == 11 ) {
        float addGap = (GAP * loop * 16.0f);
        v0 = GAP + addGap;
        v1 = GAP * 2.0f + addGap;
        v2 = GAP * 3.0f + addGap;
        v3 = GAP * 4.0f + addGap;
        v4 = GAP * 5.0f + addGap;
        v5 = GAP * 6.0f + addGap;
        v6 = GAP * 7.0f + addGap;
        v7 = GAP * 8.0f + addGap;
        v8 = GAP * 9.0f + addGap;
        v9 = GAP * 10.0f + addGap;
        v10 = GAP * 11.0f + addGap;

        for ( int row = 0; row < seqALen; ++row ) {
            idx = (int)(loop * 16.0f);
            above = v0;
            leftVal = dRow[blockRowStart + row * blockDim.x + tid];

            a = dSeqA[tid + row * blockDim.x + blockSeqAStart];

            calc(a, seqB[idx++], leftVal, oldLeft, above, max, GAP);
            oldLeft = leftVal;
            v0 = max;

            leftVal = max;
            diag = above;
            above = v1;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v1 = max;

            leftVal = max;
            diag = above;
            above = v2;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v2 = max;

            leftVal = max;
            diag = above;
            above = v3;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v3 = max;

            leftVal = max;
            diag = above;
            above = v4;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v4 = max;

            leftVal = max;
            diag = above;
            above = v5;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v5 = max;

            leftVal = max;
            diag = above;
            above = v6;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v6 = max;

            leftVal = max;
            diag = above;
            above = v7;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v7 = max;

            leftVal = max;
            diag = above;
            above = v8;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v8 = max;

            leftVal = max;
            diag = above;
            above = v9;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v9 = max;

            leftVal = max;
            diag = above;
            above = v10;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v10 = max;
        }
        dScores[blockDim.x * blockIdx.x + tid] = max;
        return;

    }


    else if ( remaining == 12 ) {
        float addGap = (GAP * loop * 16.0f);
        v0 = GAP + addGap;
        v1 = GAP * 2.0f + addGap;
        v2 = GAP * 3.0f + addGap;
        v3 = GAP * 4.0f + addGap;
        v4 = GAP * 5.0f + addGap;
        v5 = GAP * 6.0f + addGap;
        v6 = GAP * 7.0f + addGap;
        v7 = GAP * 8.0f + addGap;
        v8 = GAP * 9.0f + addGap;
        v9 = GAP * 10.0f + addGap;
        v10 = GAP * 11.0f + addGap;
        v11 = GAP * 12.0f + addGap;

        for ( int row = 0; row < seqALen; ++row ) {
            idx = (int)(loop * 16.0f);
            above = v0;
            leftVal = dRow[blockRowStart + row * blockDim.x + tid];

            a = dSeqA[tid + row * blockDim.x + blockSeqAStart];

            calc(a, seqB[idx++], leftVal, oldLeft, above, max, GAP);
            oldLeft = leftVal;
            v0 = max;

            leftVal = max;
            diag = above;
            above = v1;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v1 = max;

            leftVal = max;
            diag = above;
            above = v2;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v2 = max;

            leftVal = max;
            diag = above;
            above = v3;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v3 = max;

            leftVal = max;
            diag = above;
            above = v4;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v4 = max;

            leftVal = max;
            diag = above;
            above = v5;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v5 = max;

            leftVal = max;
            diag = above;
            above = v6;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v6 = max;

            leftVal = max;
            diag = above;
            above = v7;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v7 = max;

            leftVal = max;
            diag = above;
            above = v8;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v8 = max;

            leftVal = max;
            diag = above;
            above = v9;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v9 = max;

            leftVal = max;
            diag = above;
            above = v10;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v10 = max;

            leftVal = max;
            diag = above;
            above = v11;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v11 = max;
        }
        dScores[blockDim.x * blockIdx.x + tid] = max;
        return;

    }


    else if ( remaining == 13 ) {
        float addGap = (GAP * loop * 16.0f);
        v0 = GAP + addGap;
        v1 = GAP * 2.0f + addGap;
        v2 = GAP * 3.0f + addGap;
        v3 = GAP * 4.0f + addGap;
        v4 = GAP * 5.0f + addGap;
        v5 = GAP * 6.0f + addGap;
        v6 = GAP * 7.0f + addGap;
        v7 = GAP * 8.0f + addGap;
        v8 = GAP * 9.0f + addGap;
        v9 = GAP * 10.0f + addGap;
        v10 = GAP * 11.0f + addGap;
        v11 = GAP * 12.0f + addGap;
        v12 = GAP * 13.0f + addGap;

        for ( int row = 0; row < seqALen; ++row ) {
            idx = (int)(loop * 16.0f);
            above = v0;
            leftVal = dRow[blockRowStart + row * blockDim.x + tid];

            a = dSeqA[tid + row * blockDim.x + blockSeqAStart];

            calc(a, seqB[idx++], leftVal, oldLeft, above, max, GAP);
            oldLeft = leftVal;
            v0 = max;

            leftVal = max;
            diag = above;
            above = v1;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v1 = max;

            leftVal = max;
            diag = above;
            above = v2;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v2 = max;

            leftVal = max;
            diag = above;
            above = v3;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v3 = max;

            leftVal = max;
            diag = above;
            above = v4;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v4 = max;

            leftVal = max;
            diag = above;
            above = v5;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v5 = max;

            leftVal = max;
            diag = above;
            above = v6;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v6 = max;

            leftVal = max;
            diag = above;
            above = v7;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v7 = max;

            leftVal = max;
            diag = above;
            above = v8;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v8 = max;

            leftVal = max;
            diag = above;
            above = v9;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v9 = max;

            leftVal = max;
            diag = above;
            above = v10;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v10 = max;

            leftVal = max;
            diag = above;
            above = v11;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v11 = max;

            leftVal = max;
            diag = above;
            above = v12;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v12 = max;
        }
        dScores[blockDim.x * blockIdx.x + tid] = max;
        return;

    }


    else if ( remaining == 14 ) {
        float addGap = (GAP * loop * 16.0f);
        v0 = GAP + addGap;
        v1 = GAP * 2.0f + addGap;
        v2 = GAP * 3.0f + addGap;
        v3 = GAP * 4.0f + addGap;
        v4 = GAP * 5.0f + addGap;
        v5 = GAP * 6.0f + addGap;
        v6 = GAP * 7.0f + addGap;
        v7 = GAP * 8.0f + addGap;
        v8 = GAP * 9.0f + addGap;
        v9 = GAP * 10.0f + addGap;
        v10 = GAP * 11.0f + addGap;
        v11 = GAP * 12.0f + addGap;
        v12 = GAP * 13.0f + addGap;
        v13 = GAP * 14.0f + addGap;

        for ( int row = 0; row < seqALen; ++row ) {
            idx = (int)(loop * 16.0f);
            above = v0;
            leftVal = dRow[blockRowStart + row * blockDim.x + tid];

            a = dSeqA[tid + row * blockDim.x + blockSeqAStart];

            calc(a, seqB[idx++], leftVal, oldLeft, above, max, GAP);
            oldLeft = leftVal;
            v0 = max;

            leftVal = max;
            diag = above;
            above = v1;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v1 = max;

            leftVal = max;
            diag = above;
            above = v2;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v2 = max;

            leftVal = max;
            diag = above;
            above = v3;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v3 = max;

            leftVal = max;
            diag = above;
            above = v4;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v4 = max;

            leftVal = max;
            diag = above;
            above = v5;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v5 = max;

            leftVal = max;
            diag = above;
            above = v6;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v6 = max;

            leftVal = max;
            diag = above;
            above = v7;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v7 = max;

            leftVal = max;
            diag = above;
            above = v8;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v8 = max;

            leftVal = max;
            diag = above;
            above = v9;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v9 = max;

            leftVal = max;
            diag = above;
            above = v10;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v10 = max;

            leftVal = max;
            diag = above;
            above = v11;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v11 = max;

            leftVal = max;
            diag = above;
            above = v12;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v12 = max;

            leftVal = max;
            diag = above;
            above = v13;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v13 = max;
        }
        dScores[blockDim.x * blockIdx.x + tid] = max;
        return;

    }


    else if ( remaining == 15 ) {
        float addGap = (GAP * loop * 16.0f);
        v0 = GAP + addGap;
        v1 = GAP * 2.0f + addGap;
        v2 = GAP * 3.0f + addGap;
        v3 = GAP * 4.0f + addGap;
        v4 = GAP * 5.0f + addGap;
        v5 = GAP * 6.0f + addGap;
        v6 = GAP * 7.0f + addGap;
        v7 = GAP * 8.0f + addGap;
        v8 = GAP * 9.0f + addGap;
        v9 = GAP * 10.0f + addGap;
        v10 = GAP * 11.0f + addGap;
        v11 = GAP * 12.0f + addGap;
        v12 = GAP * 13.0f + addGap;
        v13 = GAP * 14.0f + addGap;
        v14 = GAP * 15.0f + addGap;

        for ( int row = 0; row < seqALen; ++row ) {
            idx = (int)(loop * 16.0f);
            above = v0;
            leftVal = dRow[blockRowStart + row * blockDim.x + tid];

            a = dSeqA[tid + row * blockDim.x + blockSeqAStart];

            calc(a, seqB[idx++], leftVal, oldLeft, above, max, GAP);
            oldLeft = leftVal;
            v0 = max;

            leftVal = max;
            diag = above;
            above = v1;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v1 = max;

            leftVal = max;
            diag = above;
            above = v2;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v2 = max;

            leftVal = max;
            diag = above;
            above = v3;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v3 = max;

            leftVal = max;
            diag = above;
            above = v4;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v4 = max;

            leftVal = max;
            diag = above;
            above = v5;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v5 = max;

            leftVal = max;
            diag = above;
            above = v6;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v6 = max;

            leftVal = max;
            diag = above;
            above = v7;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v7 = max;

            leftVal = max;
            diag = above;
            above = v8;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v8 = max;

            leftVal = max;
            diag = above;
            above = v9;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v9 = max;

            leftVal = max;
            diag = above;
            above = v10;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v10 = max;

            leftVal = max;
            diag = above;
            above = v11;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v11 = max;

            leftVal = max;
            diag = above;
            above = v12;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v12 = max;

            leftVal = max;
            diag = above;
            above = v13;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v13 = max;

            leftVal = max;
            diag = above;
            above = v14;
            calc(a, seqB[idx++], leftVal, diag, above, max, GAP);
            v14 = max;
        }
        dScores[blockDim.x * blockIdx.x + tid] = max;
        return;

    }


}

/**
 * This is the reduceSum function as described in the Project Requirements.
 * Will be needed to sum scores over each alignment
 */
#ifdef __DEVICE_EMULATION__
#define EMUSYNC __syncthreads()
#else
#define EMUSYNC
#endif
__global__ void reduceSum1(float *g_idata, float *g_odata)
{
	extern __shared__ float sdata[];

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;

    unsigned int i = blockIdx.x*(blockDim.x*2) + threadIdx.x;
    sdata[tid] = g_idata[i] + g_idata[i+blockDim.x];
#ifdef __DEVICE_EMULATION__
    fprintf(stderr,"[%d] trying to sync...%f\n",tid,sdata[tid]);
#endif
    __syncthreads();

    // do reduction in shared mem
    for(unsigned int s=blockDim.x/2; s>32; s>>=1)
    {
        if (tid < s)
        {
            sdata[tid] += sdata[tid + s];
        }
        __syncthreads();
    }

#ifndef __DEVICE_EMULATION__
    if (tid < 32)
#endif
    {
        sdata[tid] += sdata[tid + 32]; EMUSYNC;
        sdata[tid] += sdata[tid + 16]; EMUSYNC;
        sdata[tid] += sdata[tid +  8]; EMUSYNC;
        sdata[tid] += sdata[tid +  4]; EMUSYNC;
        sdata[tid] += sdata[tid +  2]; EMUSYNC;
        sdata[tid] += sdata[tid +  1]; EMUSYNC;
#ifdef __DEVICE_EMULATION__
        fprintf(stderr,"[%d] data is %f\n",tid,sdata[tid]);
#endif
    }

    // write result for this block to global mem
    if (tid == 0) *g_odata += sdata[0];
}

extern "C" bool setCudaDevice(int device)
{
    if ( cudaSetDevice(device) != cudaSuccess )
        return false;
    currentDeviceNum = device;
    return true;
}

extern "C" const char* initScoreMatrixTexture(float score_matrix[128][128])
{

    // Allocate CUDA array in device memory
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
    cudaArray* cuArray;
    cudaMallocArray(&cuArray, &channelDesc, 128, 128);
    cudaMemcpyToArray(cuArray, 0, 0, score_matrix, 128*128*sizeof(float), cudaMemcpyHostToDevice);

    // Set texture parameters
    tex.addressMode[0] = cudaAddressModeWrap;
    tex.addressMode[1] = cudaAddressModeWrap;
    tex.filterMode	= cudaFilterModePoint;
    tex.normalized	= false;

    // Bind the array to the texture
    cudaBindTextureToArray(tex, cuArray, channelDesc);
    if ( cudaGetLastError() != cudaSuccess )
  	return cudaGetErrorString(cudaGetLastError());
    return NULL;
}

extern "C" const char* setMatchSequence(char* matchSeq, int matchLen)
{
    //copy match sequence to constant memory
    cudaMemcpyToSymbol(seqB, matchSeq, matchLen, 0, cudaMemcpyHostToDevice);
    if ( cudaGetLastError() != cudaSuccess )
  	return cudaGetErrorString(cudaGetLastError());
    return NULL;
}

extern "C" const char* setRead(float* matchSeq, int matchLen)
{
    //copy match sequence to constant memory
										// Need to copy one for each character (4)
    cudaMemcpyToSymbol((char*)seqB, (char*)matchSeq, matchLen*4*sizeof(float), 0, cudaMemcpyHostToDevice);

    if ( cudaGetLastError() != cudaSuccess )
    	return cudaGetErrorString(cudaGetLastError());

    return NULL;
}

extern "C" const char* initSequences(int seqALen, int blockSize, int numBlocks)
{
    cudaMalloc( (void**)&d_seqA, numBlocks * blockSize * seqALen );
    if ( cudaGetLastError() != cudaSuccess )
  	return cudaGetErrorString(cudaGetLastError());
    return NULL;
}

extern "C" const char* setSequences(char* seqA, int start, int bytes)
{
    cudaMemcpy( (void*)d_seqA, (void*)&(seqA[start]), bytes, cudaMemcpyHostToDevice );
    if ( cudaGetLastError() != cudaSuccess )
  	return cudaGetErrorString(cudaGetLastError());
    return NULL;
}

extern "C" const char* copyScores(float* scores, int start, int bytes)
{
    cudaMemcpy( (void*)&(scores[start]), (void*)d_scores, bytes, cudaMemcpyDeviceToHost );
    if ( cudaGetLastError() != cudaSuccess )
    	return cudaGetErrorString(cudaGetLastError());
    return NULL;
}

extern "C" const char* initRow(int blockSize, int numBlocks, int seqALen)
{
    cudaMalloc( (void**)&d_row, sizeof(float) * numBlocks * blockSize * seqALen );
    if ( cudaGetLastError() != cudaSuccess )
    	return cudaGetErrorString(cudaGetLastError());
    return NULL;
}

extern "C" const char* initScores(int blockSize, int numBlocks)
{
    cudaMalloc( (void**)&d_scores, sizeof(float) * numBlocks * blockSize );
    if ( cudaGetLastError() != cudaSuccess )
    	return cudaGetErrorString(cudaGetLastError());
    return NULL;
}

extern "C" long getTotalGlobalMemory()
{
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, currentDeviceNum);
    return deviceProp.totalGlobalMem;
}

extern "C" int getMaxGridBlocks()
{
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, currentDeviceNum);
    return deviceProp.maxGridSize[0];
}

extern "C" void runAlignKernel(int blocks, int blockSize, int seqALen, int seqBLen, float gap)
{
    align<<<blocks, blockSize>>>(d_seqA, d_scores, d_row, seqALen, seqBLen, gap);
}

extern "C" void runSumKernel(int blocks, int blockSize, int smSize, float& sum)
{
	fprintf(stderr,"Making %d blocks with a block size of %d and a shared memory size of %d(%d)\n",blocks,blockSize,smSize,sizeof(float));
	reduceSum1<<<blocks, blockSize, smSize>>>(d_scores, &sum);
	//reduceSum1<<<blocks, blockSize>>>(d_scores, &sum);
}


