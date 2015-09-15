/* GenomeMem.h
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

#ifndef GENOME_MEM_H
#define GENOME_MEM_H

#include <iostream>
#include <sstream>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <ctime>

#include "UnitTest.h"
#include "gvector.h"
#include "hash_map.h"
#include "bin_seq.h"
#include "const_include.h"
#include "Reader.h"
#include "Genome.h"


/* 
 * The GenomeMem class is used to store the Genome that is loaded in a decent size.
 */

class GenomeMem : public Genome {
	
	public:
		/*
		 * Constructor.  Creates new Genome with filename.
		 * fn - the filename of the file to be parsed into the genome.
		 */
		GenomeMem(const char* fn);
		GenomeMem();
		
		~GenomeMem();
		
		/*
		 * GetMatches will return the vector that represents the given hash.
		 *
		 * @return the vector of the given hash.
		 */
		virtual HashLocation* GetMatches(unsigned int hash);
		virtual HashLocation* GetMatches(string &);

		/*
		 * hash_and_store() will cause the genome to be hashed--from the file fn given above.
		 * it will then store the genome in a smaller binary array.
		 */
		virtual void hash_and_store();	

		//! Mostly used in testing...
		void print_to_file(hash_map<unsigned int, HashLocation> & gh, char* ofn);

		static bool Test(ostream &os, unsigned int&);

	protected:
		
		/*!
		 * fix_hash will remove all the high-density matches in the hash.  Anything over the
		 * threshold (MAX_HASH_SIZE) will be removed.
		 */
		void fix_hash();
		void fix_hash(hash_map<unsigned int, vector<unsigned long> > *,
							 gvector<GEN_TYPE> *);

		// Reads the genome file from memory
		virtual void readHash(FILE* readF, bool ldebug);
		// Saves the genome file to memory
		virtual void saveHash(FILE* saveF, bool ldebug);
		
		/**
		 * convertToVector will convert the hashed-and-stored genome to the non-STL version
		 *
		 * @post gh and g will be deleted
		 */
		void convertToVector(hash_map<unsigned int, vector<unsigned long> > *gh,
		                     //gvector<GenomeLocation> *g);
		                     gvector<GEN_TYPE> *pg);

		//gen_piece is used to store a portion of the genome as we read it in.
		//unsigned char* gen_piece;
		hash_map<unsigned int, HashLocation> gen_hash;
		unsigned long num_hash_elements;
		hash_map<unsigned int, vector<unsigned long> > gen_hash_vec;
		
};

#endif
