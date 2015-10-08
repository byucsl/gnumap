/* const_define.h
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

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include "const_include.h"

unsigned int gBISULFITE = 0;
unsigned int gBISULFITE2 = 0;
unsigned int gATOG = 0; 
unsigned int gSNP = 0;
float gSNP_PVAL = 0.001;
bool gSNP_MONOP = false;
bool gPRINT_VCF = false;
unsigned char gINT2BASE[16];
unsigned char g_bs_CONVERSION[256];
unsigned char g_gen_CONVERSION[256];

int gVERBOSE=1;
#define DEF_MER_SIZE	10
int gMER_SIZE=0;
//unsigned int gMAX_MATCHES = 100;
unsigned int gMAX_MATCHES = 1000;
bool gUNIQUE = false;
const unsigned int gREAD_BUFFER = 1024 * 5;
unsigned int gBUFFER_SIZE = 1024*512;
unsigned int gMEM_BUFFER_SIZE = 67108864; //1024*1024*64;
//unsigned int gMEM_BUFFER_SIZE = 536870912; //1024*1024*512;
//unsigned int gMEM_BUFFER_SIZE = 1048576; //1024*1024;
//unsigned int gMEM_BUFFER_SIZE = 1024*1024*64; //1024*1024;
//unsigned int gBUFFER_SIZE = 1024;
unsigned int gHASH_SIZE = 1048576;	//gHASH_SIZE is set as 4^10
//unsigned int gMAX_KMER_SIZE = 1000;
unsigned int gMAX_KMER_SIZE = 0;
float gALIGN_SCORE = 0.9;
bool perc = true;
float gCUTOFF_SCORE = 0.0;

float gADJUST = 0.25;
float gMATCH = 3;
float gTRANSITION=-2;
float gTRANSVERSION=-3;
float gGAP = -4;
unsigned int gEDGEALIGN=3;

/*
float gTRANSITION=-1;
float gTRANSVERSION=-1;
float gGAP = -1;
*/

// These are the defaults
float gPHMM_match=0.98;
float gPHMM_syn=0.01;
float gPHMM_nsyn=0.005;
//float gPHMM_match=0.97;
//float gPHMM_syn=0.01;
//float gPHMM_nsyn=0.01;

int gMAX_GAP = 3;
char* g_adaptor = NULL;
int gPRINT_FULL=0;
float gALIGN_SCORES[256][4];
float gPHMM_ALIGN_SCORES[256][4];
int gSEQ_LENGTH=-1;

bool gFAST = false;	// Do we perform a fast alignment or a more thorough one?
bool gMATCH_POS_STRAND = true;
bool gMATCH_NEG_STRAND = true;
bool gSAM2GMP = true; // added by CJ on 04/23/2012, an option for printing gmp file in an output

float gMIN_ADAPTOR_DIFF = 0.85;	//used to define the minimum difference from the given g_adaptor sequence to the consensus
unsigned int gMIN_CHOPPED_BASES = 4;	//if there's less than 4 bases that match the g_adaptor, we're not going to chop them.
int gGEN_SKIP = 0;	// The number of bases we skip as we hash the genome
unsigned int gGEN_SIZE = 8;	// The number of bases each GEN segment contains (per each unsigned int)
unsigned int gJUMP_SIZE = 0;	// The number of bases we jump in a read (per hash)
unsigned int gMIN_JUMP_MATCHES = 2;

int gSGREX = false;

// This is for making MPI work
bool gMPI_LARGEMEM = false;

bool gILLUMINA = false;

//I'm not sure what BIN_SIZE means...It's not used anymore--only appearance is in Genome.cpp:365
//int BIN_SIZE = 10;
int gBIN_SIZE = 1;
 
// gSAVE and gREAD are global variables used to indicate whether the hashed genome should be
// read and saved or created from fasta files.
bool gSAVE=false;
bool gREAD=false;
char* gSAVE_FN;
char* gREAD_FN;

int bad_files =0;

void printReadPWM(Read* read);

int iproc=0,nproc=1;

void InitProg() {
	memset(g_gen_CONVERSION,4,256);
    g_gen_CONVERSION[(int)'a'] = g_gen_CONVERSION[(int)'A'] = 0;
    g_gen_CONVERSION[(int)'c'] = g_gen_CONVERSION[(int)'C'] = 1;
    g_gen_CONVERSION[(int)'g'] = g_gen_CONVERSION[(int)'G'] = 2;
    g_gen_CONVERSION[(int)'t'] = g_gen_CONVERSION[(int)'T'] = 3;
    g_gen_CONVERSION[(int)'n'] = g_gen_CONVERSION[(int)'N'] = 4;
    g_gen_CONVERSION[(int)'\n'] = g_gen_CONVERSION[(int)'\r'] = 5;
    // Strange white-space characters
    g_gen_CONVERSION[10] = g_gen_CONVERSION[11] = g_gen_CONVERSION[12] = g_gen_CONVERSION[13] = 5;
    g_gen_CONVERSION[(int)'\0'] = 6;
    g_gen_CONVERSION[(int)'>'] = 7;

    memset(g_bs_CONVERSION,4,256);
    g_bs_CONVERSION[(unsigned int)'a'] = g_bs_CONVERSION[(unsigned int)'A'] = 0;
    g_bs_CONVERSION[(unsigned int)'c'] = g_bs_CONVERSION[(unsigned int)'C'] = 1;
    g_bs_CONVERSION[(unsigned int)'g'] = g_bs_CONVERSION[(unsigned int)'G'] = 2;
    g_bs_CONVERSION[(unsigned int)'t'] = g_bs_CONVERSION[(unsigned int)'T'] = 3;
    g_bs_CONVERSION[(unsigned int)'n'] = g_bs_CONVERSION[(unsigned int)'N'] = 4;
    g_bs_CONVERSION[(unsigned int)'\n'] = g_bs_CONVERSION[(unsigned int)'\r'] = 5;
    // Strange white-space characters
    g_bs_CONVERSION[10] = g_bs_CONVERSION[11] = g_bs_CONVERSION[12] = g_bs_CONVERSION[13] = 5;
    g_bs_CONVERSION[(unsigned int)'\0'] = 6;
    g_bs_CONVERSION[(unsigned int)'>'] = 7;

    memset(gINT2BASE,'?',16);
    gINT2BASE[0] = 'a';
    gINT2BASE[1] = 'c';
    gINT2BASE[2] = 'g';
    gINT2BASE[3] = 't';
    gINT2BASE[4] = 'n';	
}
