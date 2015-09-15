/* const_include.h
 *
 * Copyright (C) 2009-2014 Nathan Clement
 *
 * This file is part of GNUMAP.
 *
 * GNUMAP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GNUMAP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNUMAP.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef CONST_INCLUDE_H
#define CONST_INCLUDE_H

#define gVERSION "3.1.0 BETA"

#include <cmath>
#include <string>
#include <sys/time.h>

typedef unsigned long long BIT64_t;

//If it's wanted to have a base-pair resolution, un-comment these lines and comment the next two.
// NOTE:  CANNOT CHANGE THIS TYPE WITHOUT RUINING MPI
typedef unsigned int GEN_TYPE;
#define GEN_PACK	8

extern unsigned int gGEN_SIZE;

//#define MAX_MISMATCH 100

// This flag will create a lot of verbose output.
//#define DEBUG

#define MAX_LN_SZ	1024
#define MAX_NAME_SZ	1024
#define MAX_CIGAR_SZ	1024
#define MAX_MER_SIZE	32

#define NUM_CHAR_VALS	5
#ifdef _INDEL
#define NUM_SNP_VALS	NUM_CHAR_VALS+2
#else
#define NUM_SNP_VALS	NUM_CHAR_VALS
#endif
#define A_POS 0
#define C_POS 1
#define G_POS 2
#define T_POS 3
#define N_POS 4
#define DEL_POS 5 /* read gap, genome insertion */
#define INS_POS 6 /* genome gap, read insertion */

#define READS_PER_PROC	1024*2

#define MUTEX_LOCK(t) pthread_mutex_lock(t)
#define MUTEX_UNLOCK(t) pthread_mutex_unlock(t)

#if defined( INTDISC )
#define FLTFROM255(t, c) ((float)(c)/255.0f*(t))
#define CHARFROMFLT(t, f) (unsigned char)(round((f)/(t)*255))

inline void adjust255(float prevV, float newV, unsigned char *curvals) {
	int i;
	if(prevV != 0 && prevV!=newV && newV>1)
		i=0;
	for(i=0; i<NUM_SNP_VALS; i++) {
		curvals[i] = CHARFROMFLT(newV, FLTFROM255(prevV, curvals[i]));
	}
}
#endif


/******** Needed for Probcons *********/
//int VERSION = 1;
const int NumInsertStates = 2;
/**************************************/

extern unsigned char g_bs_CONVERSION[256];
extern unsigned char g_gen_CONVERSION[256];

extern unsigned char gINT2BASE[16];
extern int gVERBOSE;
extern int gMER_SIZE;
extern unsigned int gBUFFER_SIZE;
extern const unsigned int gREAD_BUFFER;
extern unsigned int gHASH_SIZE;
extern unsigned int gMAX_HASH_SIZE;
extern bool gUNIQUE;

extern bool gMATCH_POS_STRAND;
extern bool gMATCH_NEG_STRAND;

extern float gALIGN_SCORE;
extern float gALIGN_SCORES[256][4];
extern float gPHMM_ALIGN_SCORES[256][4];

extern int gMAX_GAP;
extern float gGAP;
extern float gADJUST;
extern float gMATCH;
extern float gTRANSITION;
extern float gTRANSVERSION;
extern float gGAP;


extern int gGEN_SKIP;
extern unsigned int gBISULFITE;
extern unsigned int gBISULFITE2;
extern unsigned int gATOG;
extern unsigned int gSNP;
extern float gSNP_PVAL;
extern bool gSNP_MONOP;
extern bool gPRINT_VCF;

extern int gBIN_SIZE;

extern char* g_adaptor;
extern float gMIN_ADAPTOR_DIFF;
extern unsigned int gMIN_CHOPPED_BASES;

extern bool gILLUMINA;

extern bool gSAVE;
extern bool gREAD;
extern char* gSAVE_FN;
extern char* gREAD_FN;

extern int iproc, nproc;
extern bool gMPI_LARGEMEM;
extern bool gSAM2GMP; //added by CJ
extern unsigned int gEDGEALIGN; //added by CJ

enum FileType {
	CHROM=0,
	SEQUENCE=1
};

enum SeqType {
	INT,
	PRB,
	FASTA,
	FASTQ
};

struct Read {
	char* name;
	Read() { name = 0; length=0; }
	
	Read(float** p, int l) {
		pwm = p; length = l;
	}

	unsigned int length;
	float** pwm;
	std::string seq;
	std::string fq;
};

struct GenomeLocation {
	GenomeLocation() {
		amount = 0.0;
		packed_char = 0;
	}

	GenomeLocation(GEN_TYPE pc, float amt)
		: amount(amt), packed_char(pc) {}

	float amount;
	GEN_TYPE packed_char;
};

#define READ_FAILED		-1
#define READ_TOO_SHORT	-2
#define READ_TOO_POOR	-3
#define READ_TOO_MANY	999999
#define SAME_DIFF		0.00001

#define POS_STRAND	0
#define NEG_STRAND	1

struct TopReadOutput {
	char READ_NAME[MAX_NAME_SZ];
	char CHR_NAME[MAX_NAME_SZ];
	unsigned long CHR_POS;
	int strand;
	int MAPQ;
	float A_SCORE;
	float POST_PROB;
	int SIM_MATCHES;
	char CIGAR[MAX_CIGAR_SZ];
	unsigned int readIndex;
	std::string consensus;
	std::string qual;
};

/**********************************************
 * These next are for Command-Line Parsing
 **********************************************/
#define INVALID_ARG -1
#define NO_MATCHING_ARG -2
#define UNKNOWN_ARG -3
#define PARSE_ERROR -10

struct CMDLine {
	int errnum;
	int errpos;
};
/**********************************************/

// Needed mostly for testing, but also other stuff
/* Return the current time in seconds, using a float precision number. */
inline double When() {
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}

#endif //CONST_INCLUDE_H
