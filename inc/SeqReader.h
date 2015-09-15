/* SeqReader.h
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

#ifndef _SEQREADER_H
#define _SEQREADER_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "const_include.h"
#include "bin_seq.h"
#include "UnitTest.h"
#include "Exception.h"

#include <execinfo.h>


using namespace std;

const unsigned int READ_BUFFER = 1024;

struct SeqFile {
	SeqType type;
	string file_name;

	SeqFile() {
		file_name = "";
	}

	bool operator ()(SeqFile &other) {
		if( file_name == other.file_name)
			return true;
		return false;
	}
};

class SeqReader {
	public:
	
		/*
		 * The Constructor with a given file name.  Used to set everything up, determine the number
		 * of sequence files, and create the set of file strings to use.
		 * 
		 * @param fn - the string received from the SeqDriver program.
		 */
		SeqReader(const char* fn);
		SeqReader();
		SeqReader(const SeqReader &);

		SeqReader& operator =(const SeqReader&);
		
		void Construct(const char* fn);
		
		~SeqReader();

		/*!
		 * use will enable you to use this file after the SeqDriver has already
		 * been initialized.  
		 *
		 * @param fn - the file(s) to use as sequences
		 *				(can be either _int or _prb)
		 */
		void use(const char* fn);
		void use(const string &str);
		
		/**
		 * I suddenly realized that by just "eating" sequences in the SeqManager we
		 * were creating an incredible memory leak.  So, instead of actually reading them
		 * in, we'll just increment the file size here.
		 *
		 * @param fn The file to use
		 * @param nBurn The number of sequences to burn
		 */
		void use(const char* fn, unsigned int nBurn);
		void use(const string &str, unsigned int nBurn);

		
		string GetFilename() {
			return this_filename;
		}

		/*
		 * GetNextSequence will return the next available PWM for the sequence.
		 * If it has run out of sequences, it will read more in from memory.
		 *
		 * @return the PWM representing the next sequence.
		 */
		//vector<vector<double> > GetNextSequence();
		//vector<vector<double> > PeekNextSequence();
		Read* GetNextSequence();
		Read* PeekNextSequence();

		void PrintSeq(double &seq);		
		string GetConsensus(vector<vector<double> > &read);
		string GetConsensus(Read* read);
		static Read* getReadFromStrings(char* name, string &consensus, string &qual);
		
		unsigned int GetNumSeqs();
		unsigned int GetNumSeqs(const string&);

		
		/*!
		 * Used with the SeqDriver in determining how many files to use...
		 */
		unsigned int GetNumSeqFiles();
		vector<string>& GetSeqFileVector();
		
		static bool Test(ostream &os, unsigned int &);

	private:
		void Init(const char* fn);
		void ReadBatch();
		void clearReads(Read* readArr[]);
		
		inline bool burnPRBINT(unsigned int nBurn);
		inline bool burnFASTA(unsigned int nBurn);
		inline bool burnFASTQ(unsigned int Burn);
		
		bool get_more_int();
		bool get_more_prb();
		bool get_more_fasta();
		bool get_more_fastq();

		static double Q2Prb_ill(double Q);
		static double Q2Prb_std(double Q);
		static inline float Compare(string, string, int);
		void FixReads(vector<vector<double> > &read);
		unsigned int FixReads2(string seq_read);
		//void FixReads(Read &pwm);
		
		void find_type();
		SeqType find_type(const string&);

		static char max_char(vector<double> &chr);
		//char max_char(float chr[4]);
	
		//vector<vector<vector<double> > > reads;
		//Read reads[READ_BUFFER];	//an array of pointers to reads
		Read* reads[READ_BUFFER];	//an array of pointers to reads
		unsigned int num_reads;
		unsigned int seq_counter;

		bool file_finished;
		unsigned long offset;	//The offset in the individual sequence file.
		unsigned int seq_num;	//To keep track which sequence we are on--global for the entire file
		//const char* filename;
		//SeqType type;
		string this_filename;
		SeqType this_type;
		vector<string> seq_files;	//A list of the sequence files to read.

};

#endif
