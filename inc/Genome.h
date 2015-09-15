/* Genome.h
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

#ifndef GENOME_H
#define GENOME_H

#ifdef	MPI_RUN
#include <mpi.h>
#endif

#include <iostream>
#include <sstream>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <ctime>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <math.h>
#include <pthread.h>
#include <execinfo.h>
#include <time.h>

#include "gsl/gsl_cdf.h"
#include "gvector.h"
#ifndef GENOME_STL
#include "hash_map.h"
#endif
#include "bin_seq.h"
#include "const_include.h"
#include "Reader.h"

#include "centers.h"

/* 
 * The Genome class is used to store the Genome that is loaded in a decent size.
 */

// There are 4 for each nucleotide, plus one for gaps and one for insertions
#define MAX_NAME_LEN	200
#define MAX_FN_LEN		1024

struct HashLocation {
	HashLocation() {
		size = 0;
	}

	unsigned long* hash_arr;
	unsigned long size;
};

// An ugly hack, but it made it object-oriented.
#if defined(GENOME_SUFFIX)
#define MAP_USE_STRING
	typedef string MAP_TYPE;
typedef map<MAP_TYPE, HashLocation> map_type;
#else // GENOME_SUFFIX
#if defined(GENOME_STL)
// The following are for GENOME_STL, but needed here
//#define MAP_USE_STRING
#ifdef MAP_USE_STRING
	typedef string MAP_TYPE;
#else //MAP_USE_STRING
	typedef unsigned long MAP_TYPE;
#endif //MAP_USE_STRING

typedef map<MAP_TYPE, HashLocation> map_type;

#else // GENOME_STL

typedef hash_map<unsigned int, HashLocation> map_type;

#endif // GENOME_STL
#endif // GENOME_STL

class Genome {
	friend class bin_seq;
	
	public:
		/*
		 * Constructor.  Creates new Genome with filename.
		 * fn - the filename of the file to be parsed into the genome.
		 */
		Genome(const char* fn);
		Genome();
		
		virtual ~Genome();
		
		void use(const char* fn, unsigned long, unsigned long);
		/**
		 * Wrapper function that just calls the previous one with 0,0
		 */
		void use(const char* fn);
		
		/*
		 * Returns a string from begin to size.  Forward strand of DNA
		 * @par begin - the beginning of the string to return
		 * @par size - the size of the string to return.
		 */
		string GetString(const unsigned long begin, const unsigned int size);
		string GetString(const unsigned long begin, const unsigned int size, 
					const unsigned int hash_start, const bool extend_past);
		char GetChar(const unsigned long begin);
		
		/*
		 * GetMatches will return the vector that represents the given hash.
		 *
		 * @return the vector of the given hash.
		 */
		virtual HashLocation* GetMatches(unsigned int hash) = 0;
		virtual HashLocation* GetMatches(string &) = 0;
		
		/**
		 * LoadGenome is responsible for loading the genome into memory--whether it be by reading from
		 * a binary file or from a set of fasta files, it should be the same.

		 */
		void LoadGenome();
		/**
		 * This function was pulled out of LoadGenome so it could also be used in Store Genome.
		 * It simply allocates space for the 'read' and 'genome' arrays, as necessary.
		 */
		void make_extra_arrays();

		/*
		 * hash_and_store() will cause the genome to be hashed--from the file fn given above.
		 * it will then store the genome in a smaller binary array.
		 */
		virtual void hash_and_store() = 0;
		void inc_counter(int,unsigned int&, unsigned long&, unsigned int);
		void ReplaceSpaceWithUnderscore(string &str);

		
		/**
		 * count is used to determine the genome size to split for using MPI
		 */
		unsigned long count();
		
		/**
		 * StoreGenome is used by the sam2sgr program to store the genome (doesn't need to hash it)
		 */
		virtual void StoreGenome() { StoreGenome(true); }
		virtual void StoreGenome(bool make_extra);
		
		/*! 
		 * AddScore will add the given score, amt, to the genomic position, pos.
		 *
		 * @param pos - the position on the genome.
		 * @param amt - the amount to be given to this specific position on the genome.
		 */
		float GetScore(const unsigned long &pos);
		void AddScore(const unsigned long &pos, const float &amt);
		void AddScore(const unsigned long &pos, const float &amt, float a_add, float b_add, float gen_a, float gen_b);
		void AddScore(const unsigned long &pos, const float &amt, float a_add, float c_add, float g_add, float t_add,
					float gen_a, float gen_c, float gen_g, float gen_t);

		void AddSeqScore(unsigned long, const float*, const float);
		void AddSeqScore(unsigned long, const float, unsigned int which);
		void AddGenScore(unsigned long, const float, char which);
		
		/*! 
		 * GetPos is used mostly in printing out the matches.  It will find the position
		 * on the genome and the corresponding chromosome name.
		 *
		 * @param pos - the position to find the local chromosomal placement from.
		 *
		 * @return a string representing the chromosomal name and local position.
		 */
		string GetPos(const unsigned long &pos, const int strand);
		string GetPos(const pair<unsigned long,int> p_pos);
		
		/*!
		 * GetPosPair is used by the ScoredSeq class to return a pair with the sequence name and position--
		 * instead of the single string from GetPos
		 *
		 * @param pos The GNUMAP position of the sequence
		 */
		pair<string,unsigned long> GetPosPair(const unsigned long pos);
		
		/**
		 * GetAbsolutePosition will do the opposite from GetPos, meaning it will take the
		 * chromosome name and local position and produce the absolute position used by GNUMAP
		 */
		unsigned long GetAbsolutePosition(const char* chr, const unsigned long &pos);

		/**
		 * GetStringAtPosition will go one step further from GetAbsolutePosition and will
		 * actually produce the string at that position
		 */
		string GetStringRelative(const char* chr, const unsigned long &pos, const unsigned int length);

		
#ifdef SET_POS
		/**
		 * Will be used to set the thread ID for this genomic position.
		 * If this location has already been searched by my thread, it will return a -1 and
		 * not set the location.
		 * 
		 * The thread value at pos will always be set to myID after this function is called.
		 */
		//efficiency is not proved, so we're tossing it out
		bool getSetThreadID(long unsigned int pos, unsigned int myID);
		void unsetThreadID(long unsigned int pos, unsigned int myID);
#endif

		/**
		 * This is needed for the MPI_Reduce call on this genome.  We need to pass back the
		 * pointer to the GenomeLocation object
		 */
		float* GetGenomeAmtPtr() { return amount_genome; }
#ifdef DISCRETIZE
		center_d* GetGenomeAllotPtr() { return read_allot; }
#elif defined( INTDISC )
		unsigned char* GetGenomeAllotPtr() { return reads; }
#else
		float* GetGenomeAPtr() { return reads[A_POS]; }
		float* GetGenomeCPtr() { return reads[C_POS]; }
		float* GetGenomeGPtr() { return reads[G_POS]; }
		float* GetGenomeTPtr() { return reads[T_POS]; }
		float* GetGenomeNPtr() { return reads[N_POS]; }
#ifdef _INDEL
		float* GetGenomeIPtr() { return reads[INS_POS]; }
		float* GetGenomeDPtr() { return reads[DEL_POS]; }
#endif
#endif // DISCRETIZE
		
		unsigned long size();
		
		/*!
		 * PrintFinal will print out the entire genome in the SGR format used for viewing in the
		 * Integrated Genome Browser (IGB).  Also prints in SGREX if desired.
		 *
		 * Made a modification, so now it can also print in append mode.  Previous function
		 * is a stub call to the current function.
		 *
		 * @param fn - the filename to read to (will have .sgr appended to the end)
		 * @param append - should I append or write a new file?
		 */
		void PrintFinalSNP(const char* fn, bool append);
		void PrintSNPCall(const unsigned long, FILE*);
		void PrintFinalBisulfite(const char* fn, bool append);
		void PrintFinalVCF(const char* fn, bool append);
		void PrintFinalSGR(const char* fn, bool append);
		void PrintFinal(const char* fn, bool append);
		void PrintFinal(const char* fn);
		void AppendFinal(const char* fn);


		inline unsigned int max_pos(float[], unsigned int);

		/*
		 * is_snp determines if this position is a SNP or not.  It takes as parameters
		 * the character in the genome and the float array represeting the coverage for
		 * all nucleotides at this position.  Currently, there are only two characteristics
		 * that make it a SNP:  It must exceed a pre-defined MIN_COVERAGE amount, and it
		 * must exceed a MIN_RATIO.  
		 * The ratio is the value of the max nucleotide divided by the second highest 
		 * nucleotide
		 * 
		 * @param to_check The character at this position in the genome
		 * @param chars An array representing the nucloetides at all positions
		 * @return -1 if it is not a SNP, else the snp position if it is
		 */
		//inline int is_snp(char, float[5]);
		/**
		 * @param chars An array representing the nucloetides at all positions
		 * @param snp_pos The position for the snp(?)
		 * @return the p-value of the maximum position being the correct nucleotide
		 */
		double is_snp(float chars[5], int &, int&, bool&);
		inline double LRT(float chars[5], int &snp_pos);
		inline double dipLRT(float chars[5], int &, int&, bool&);

		
		//! Mostly used in testing...
		virtual void print_to_file(map_type& gh, char* ofn) = 0;

	protected:
		void Init(const char* fn);
		//unsigned long int store(vector<GenomeLocation> *g, unsigned char* gen, 
		//void store(gvector<GenomeLocation> *g, unsigned char* gen, 
		//						unsigned long int, unsigned long &);
		void store(gvector<GEN_TYPE> *pg, unsigned char* gen, 
								unsigned long int, unsigned long &);
		
		/*!
		 * fix_hash will remove all the high-density matches in the hash.  Anything over the
		 * threshold (MAX_HASH_SIZE) will be removed.
		 */
		virtual void fix_hash() = 0;

		inline int base2int(char c);
		inline char int2base(int i);
				
		string reverse(string &s);

		// Reads the genome file from memory
		unsigned int readGen(char* fn);
		virtual void readHash(FILE* readF, bool ldebug) = 0;
		virtual void readGenome(FILE* readF, bool ldebug);
		// Saves the genome file to memory
		unsigned int saveGen(char* fn);
		virtual void saveHash(FILE* saveF, bool ldebug) = 0;
		virtual void saveGenome(FILE* saveF, bool ldebug);
		
		GEN_TYPE offset;
		unsigned int bin_offset;
		Reader* reader;

		//genome is the compressed genome
		//GenomeLocation* genome;
		GEN_TYPE* packed_genome;
		float* amount_genome;
		
		// Should be a char so we can have 8 bits
		char* gs_positions;
		unsigned long gs_positions_size;
		
		unsigned long genome_size;

#if defined( DISCRETIZE )
		// All we need is a character array for the reads as we'll discretize them
		center_d* read_allot;
#elif defined( INTDISC )
		// "Parts per 255" for each character
		// one-dimensional array of size genome_size*NUM_SNP_VALS
		unsigned char* reads;
#else
		float* reads[NUM_SNP_VALS];
#endif

		//gen_piece is used to store a portion of the genome as we read it in.
		//unsigned char* gen_piece;
		unsigned long num_hash_elements;
		vector<pair<string,unsigned long> > names;
		
		unsigned long my_start;	// The position we should start reading from
		unsigned long my_end;
		bool read_all;
		bool include_hash;	// So we can just store the genome for sam2sgr
};

#endif
