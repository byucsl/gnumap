/* Genome.cpp
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

#include "Genome.h"

/*
 * Constructor.  Creates new Genome with filename.
 * fn - the filename of the file to be parsed into the genome.
 */
Genome::Genome(const char* fn) :
	offset(0), bin_offset(0), packed_genome(NULL), amount_genome(NULL), gs_positions(NULL), genome_size(0) {

	Init(fn);
	reader = new Reader();
}

Genome::Genome() :
	offset(0), bin_offset(0), reader(NULL), packed_genome(NULL), amount_genome(NULL), gs_positions(NULL), genome_size(0) {

}

Genome::~Genome() {
	//fprintf(stderr,"[-/%d] Genome destructor...\n",iproc);
	delete reader;
	
	if(packed_genome)
    {
		free(packed_genome);
    }

	if(amount_genome)
    {
		free(amount_genome);
    }
#ifdef SET_POS
	if(gs_positions)
    {
		free(gs_positions);
    }
#endif
	
	if(gSNP)
    {
#if defined( DISCRETIZE )
		delete[] read_allot;
#elif defined( INTDISC )
		// To maintain cache consistency, store them as reads[pos][char]
		delete[] reads;
		reads = 0;
#else

        for(unsigned int i=0; i<NUM_SNP_VALS; i++)
        {
            delete[] reads[i];
            reads[i] = 0;
        }
#endif // DISCRETIZE/INTDISC
    }

    if(gBISULFITE)
    {
#if defined( DISCRETIZE )
		delete[] read_allot;
#elif defined( INTDISC )
		delete[] reads;
		reads = 0;
#else
        for(unsigned int i=0; i<NUM_SNP_VALS; i++)
        {
            delete[] reads[i];
            reads[i] = 0;
        }
#endif // DISCRETIZE/INTDISC
    }

    if(gATOG)
    {
    	// We can have the SNP option and the ATOG option at the same time
#if defined( DISCRETIZE )
		delete[] read_allot;
#elif defined( INTDISC )
		if(reads)
        {
			delete[] reads;
			reads = 0;
		}
#else
        for(unsigned int i=0; i<NUM_SNP_VALS; i++)
        {
        	if(reads[i])
            {
				delete[] reads[i];
				reads[i] = 0;
			}
        }
#endif // DISCRETIZE/INTDISC
    }
} 

void Genome::Init(const char* fn)
{

	// Don't need to init names if we're going to do so in the read file
	if(!gREAD)
    {
		if(names.size()) //if we've already set up the "names" vector, just return.
        {
			return;
        }

		char files[strlen(fn)];
		strcpy(files,fn);
		
		const char* delims = ", \"\n";
		char* token = strtok(files,delims);
		if(gVERBOSE > 1)
        {
			cout << "Found:" << endl;
        }

		while(token != NULL)
        {
			if(gVERBOSE)
            {
				cout << "\t" << token << endl;
            }
			string temp_str = token;
			
			// Get rid of spaces that make it impossible to handle
			ReplaceSpaceWithUnderscore(temp_str);
			
			names.push_back(pair<string,unsigned long>(temp_str,0));
			
			token = strtok(NULL,delims);
		}

		//Push a dummy string on the end so we can know when to end...
		names.push_back(pair<string,unsigned long>("end",0));
	}

	//clear the *_spots files
#if defined( DISCRETIZE )
	read_allot = 0;
#elif defined( INTDISC )
	reads = 0;
#else
	for(int i=0; i<NUM_SNP_VALS; i++)
    {
		reads[i] = 0;
	}
#endif // DISCRETIZE/INTDISC

}

void Genome::use(const char* fn) {
	use(fn,0,0);
}

void Genome::use(const char* fn, unsigned long s, unsigned long e) {
	my_start = s;
	my_end = e;
	if(my_start >= my_end)
		// read everything
		read_all = true;
	else
		read_all = false;
	
	//fprintf(stderr,"Reading Genome from %lu to %lu, read_all=%d\n",s,e,read_all);
	
	if(reader != NULL)
    {
		cout << "non-null reader" << endl;
		delete reader;
	}

	Init(fn);
	reader = new Reader();
}

void Genome::make_extra_arrays() {
	//if we are doing the bisulfite sequencing, we want to record what nucleotide
	//comes at each position
	if(gSNP) {
		printf("Finding SNPs\n");
#if defined( DISCRETIZE )
		read_allot = new center_d[genome_size];
		if(!read_allot)
			throw new Exception("Could not allocate space for SNP calling");
		memset(read_allot,0,genome_size*sizeof(char));
#elif defined( INTDISC )
		reads = new unsigned char[genome_size*NUM_SNP_VALS];
		if(!reads)
			throw new Exception("Could not allocate space for SNP calling");
		memset(reads,0,genome_size*NUM_SNP_VALS*sizeof(char));
#else
		for(unsigned int i=0; i<NUM_SNP_VALS; i++) {
			reads[i] = new float[genome_size];
			if(!reads[i])
				throw new Exception("Could not allocate space for SNP calling");
			memset(reads[i],0,genome_size*sizeof(float));
		}
#endif
	
	}
	if(gBISULFITE) {
		// Might already have been done under SNP calling
		printf("Using bisulfite conversion techniques\n");
#if defined( DISCRETIZE )	
		if(!read_allot) {
			read_allot = new center_d[genome_size];
			if(!read_allot)
				throw new Exception("Could not allocate space for BISULFITE calling");
			memset(read_allot,0,genome_size*sizeof(char));
		}
#elif defined( INTDISC )
		if(!reads) {
			reads = new unsigned char[genome_size*NUM_SNP_VALS];
			if(!reads)
				throw new Exception("Could not allocate space for Bisulfite calling");
			memset(reads,0,genome_size*NUM_SNP_VALS*sizeof(char));
		}
#else
		for(unsigned int i=0; i<NUM_SNP_VALS; i++) {
			// Might already have been done under SNP calling
			if(!reads[i])
				reads[i] = new float[genome_size];
			if(!reads[i])
				throw new Exception("Could not allocate space for Bisulfite calling");
			memset(reads[i],0,genome_size*sizeof(float));
		}
#endif
		
	}
	if(gATOG) {
		printf("Using AtoG conversion techniques\n");
		// Might already have been done under SNP calling
#if defined( DISCRETIZE )
		if(!read_allot) {
			read_allot = new center_d[genome_size];
			if(!read_allot)
				throw new Exception("Could not allocate space for A to G calling");
			memset(read_allot,0,genome_size*sizeof(char));
		}
#elif defined( INTDISC )
		if(!reads) {
			reads = new unsigned char[genome_size*NUM_SNP_VALS];
			if(!reads)
				throw new Exception("Could not allocate space for Bisulfite calling");
			memset(reads,0,genome_size*NUM_SNP_VALS*sizeof(char));
		}
#else
		for(unsigned int i=0; i<NUM_SNP_VALS; i++) {
			if(!reads[i])
				reads[i] = new float[genome_size];
				
			if(!reads[i])
				throw new Exception("Could not allocate space for A to G calling");
			memset(reads[i],0,genome_size*sizeof(float));
		}
#endif // DISCRETIZE/INTDISC
	}
}

void Genome::ReplaceSpaceWithUnderscore(string &str) {
	for(unsigned int i=0; i<str.length(); i++) {
		if(str[i] == ' ' || str[i] == 10 || str[i] == 11 || str[i] == 12 || str[i] == 13)
			str[i] = '_';
	}

}


#define READ_EXTRA	1024

/*
 * store() will store the genome into a smaller array of BIT64_t portions.
 * Each character is represented by 4 bits (because we require the 4 nucleotides
 * and an additional N character), and can be accessed directly.  
 *
 * One of the problems we have now is the last few characters in gen.  They will
 * appear the second time at the beginning of the genome yet again, throwing off the 
 * numbering system. 
 *
 * @par gen - the piece of the genome to store.
 * @par num - the number of characters to store.
 * @par num_bases - the base counter.  And tells when to start recording
 */
//void Genome::store(gvector<GenomeLocation> *g, unsigned char* gen, 
//						unsigned long int num, unsigned long &num_bases) {
void Genome::store(gvector<GEN_TYPE> *pg, unsigned char* gen, 
						unsigned long int num, unsigned long &num_bases) {
	// DONE!!
	if(num_bases > my_end+READ_EXTRA && !read_all)
		return;
		
	bool start_yet = num_bases >= my_start ? true : false;

	//fprintf(stderr,"READ_ALL: %d start_yet: %d\n",read_all,start_yet);

	unsigned long gen_pos=0;
//	unsigned long int num_bases=0;
	
	while(gen_pos<num) {
		GEN_TYPE temp = 0;
		
		if(!read_all && num_bases > my_end+READ_EXTRA)
			break;

		//if( read_all || (start_yet && (num_bases <= my_end+READ_EXTRA)) ) {
			if(bin_offset != 0) {
				//fprintf(stderr,"Is this the error?  Offset is %d and my_start is %ul\n",offset,my_start);
				//temp = (*g)[offset-my_start/gGEN_SIZE].packed_char;
				//temp = (*g)[offset].packed_char;
				temp = (*pg)[offset];
				pg->pop_back();
			}
		//}
			
		for(; (bin_offset<GEN_PACK) && gen_pos<num; gen_pos++) {	//can fit 16 chars into a BIT64_t 
			//int int_to_add = base2int(gen[gen_pos]);
			int int_to_add = g_gen_CONVERSION[(int)gen[gen_pos]];
			
			
			//If it's a valid char.
			if(int_to_add < 5) {
				num_bases++;
				if(num_bases > my_start)
					start_yet = true;
				if(!read_all && num_bases > my_end+READ_EXTRA) {
					// remember this location, then return
					//g->push_back(GenomeLocation(temp,0));
					pg->push_back(temp);
					return;	// DONE!
				}

				if(start_yet) {
					temp <<= 4;
					temp |= (int_to_add & 0xf);
					//fprintf(stderr,"Adding %c: %x(%d), now %0.8x\n",gen[gen_pos],int_to_add & 0xf,int_to_add,temp);
					++bin_offset;

				}
									
			}
			//If we come accross the flag for the end of the stream--should never happen...
			else if(int_to_add == 6) { 
				gen_pos = num;
				//cout << "Found end." << endl;
				break;
			}
			else if(int_to_add == 7) {	//come accross a new chromosome name.
				string temp_str = (char*)(gen+gen_pos);
				//while(gen[gen_pos] != '\n' && gen[gen_pos] != '\0' && gen_pos<num )
				while( gen[gen_pos] != '\n' && gen[gen_pos] != '\0' )
					gen_pos++;	//increment the counter until we've found the end of the name
			}
			else ;//do nothing--it's an invalid character and not an end-of-line
		}
				
		if(bin_offset != 0)	{ //if we've made any changes to this value...
			try {
				//g->push_back(pair<GEN_TYPE,float>(temp,0));
				if(read_all || (start_yet && num_bases <= my_end+READ_EXTRA)) {
					//fprintf(stderr,"Pushing %0.8x onto genome\n",temp);
					//g->push_back(GenomeLocation(temp,0));
					pg->push_back(temp);
				}
			}
			catch(bad_alloc& e) {
				printf("Found ERROR!!!:  %lu elements out of max of %lu\n",pg->size(),pg->max_size());
				//printf("Total space is %lu\n",g->size()*sizeof(GenomeLocation));
				printf("Total space is %lu\n",pg->size()*sizeof(GEN_TYPE));
				throw e;
			}
		}
		else ;	//no changes.  Coming out of the previous for-loop
		
		
		//means we've filled it fully--so increment the offset.
		if(bin_offset == GEN_PACK) {
			if(start_yet) {
				offset++;
				bin_offset = 0;
			}
		}


	}
	
	//fprintf(stdout,"Just read %ld bases out of %lu assigned\n",num_bases, num);
	return;
}

char Genome::GetChar(const unsigned long pos) {

	unsigned long pos_in_gen = (pos-my_start)/GEN_PACK;
	assert(pos_in_gen <= genome_size/GEN_PACK);

	GEN_TYPE t_temp = packed_genome[pos_in_gen];
	unsigned long bit_offset = (pos-my_start)%GEN_PACK;
	return gINT2BASE[0xf & (t_temp >> (4*(GEN_PACK-bit_offset-1)))];
}

float Genome::GetScore(const unsigned long &pos) {

	unsigned long pos_temp = (pos-my_start)/gGEN_SIZE;
	//return genome[pos_temp].amount;
	return amount_genome[pos_temp];
}

/**
 * INVARIANT:  must add to the genome previous to this call
 *			thus, amount_genome[pos] will never be 0
 */
void Genome::AddSeqScore(unsigned long pos, const float* amt, const float scale) {
	//fprintf(stderr,"float* - pos: %lu, myStart: %lu\n",pos,my_start);
	unsigned int loc = (pos-my_start)/gGEN_SIZE;
	
/*
	if(amt[A_POS] * scale > 0.0)
		fprintf(stderr,"For position %lu, total A is %f = %f*%f\n",pos,amt[A_POS]*scale,amt[A_POS],scale);
	if(amt[C_POS] * scale > 0.0)
		fprintf(stderr,"For position %lu, total C is %f = %f*%f\n",pos,amt[C_POS]*scale,amt[C_POS],scale);
	if(amt[G_POS] * scale > 0.0)
		fprintf(stderr,"For position %lu, total G is %f = %f*%f\n",pos,amt[G_POS]*scale,amt[G_POS],scale);
	if(amt[T_POS] * scale > 0.0)
		fprintf(stderr,"For position %lu, total T is %f = %f*%f\n",pos,amt[T_POS]*scale,amt[T_POS],scale);
	if(amt[N_POS] * scale > 0.0)
		fprintf(stderr,"For position %lu, total N is %f = %f*%f\n",pos,amt[N_POS]*scale,amt[N_POS],scale);
//*/

#if defined( DISCRETIZE )
	int i;
	center_t *map_c = centers[read_allot[loc]];
	center_t gen_amt[DIM_CENTERS];
	// Convert to floating point representation
	for(i=0; i<DIM_CENTERS; i++)
		// we've added in the value previous to this call.  Subtract it out first
		gen_amt[i] = map_c[i]*(amount_genome[loc]-scale);
	
	gen_amt[A_POS] += amt[A_POS]*scale;
	gen_amt[C_POS] += amt[C_POS]*scale;
	gen_amt[G_POS] += amt[G_POS]*scale;
	gen_amt[T_POS] += amt[T_POS]*scale;
	gen_amt[N_POS] += amt[N_POS]*scale;


	// Convert back, including the added value
	for(i=0; i<DIM_CENTERS; i++)
		gen_amt[i] /= amount_genome[loc];
#ifdef FIND_EXHAUSTIVE
	read_allot[loc] = find_center_exhaustive(gen_amt);
#else
	read_allot[loc] = find_center(gen_amt);
#endif

#elif defined( INTDISC )
	// What do we need to do here?  Use the "parts per 255" again
	// Need to scale the variable "scale", then we can add it to the position
	// then round it to the closest 255
	float loc_amt = amount_genome[loc];
	float prev_amt = loc_amt-scale;
	unsigned int reads_loc = loc*NUM_SNP_VALS;
	unsigned char newC;
	
	//fprintf(stderr,"prev: %d %d %d %d %d\n",reads[reads_loc+A_POS],reads[reads_loc+C_POS],reads[reads_loc+G_POS],reads[reads_loc+T_POS],reads[reads_loc+N_POS]);
	/*
	float fnA = FLTFROM255(prev_amt, reads[reads_loc+A_POS]) + amt[A_POS]*scale;
	float fnC = FLTFROM255(prev_amt, reads[reads_loc+C_POS]) + amt[C_POS]*scale;
	float fnG = FLTFROM255(prev_amt, reads[reads_loc+G_POS]) + amt[G_POS]*scale;
	float fnT = FLTFROM255(prev_amt, reads[reads_loc+T_POS]) + amt[T_POS]*scale;
	float fnN = FLTFROM255(prev_amt, reads[reads_loc+N_POS]) + amt[N_POS]*scale;
	fprintf(stderr,"prev_floats(%f): %f %f %f %f %f\n",prev_amt,fnA-amt[A_POS]*scale,fnC-amt[C_POS]*scale,fnG-amt[G_POS]*scale,fnT-amt[T_POS]*scale,fnN-amt[N_POS]*scale);
	fprintf(stderr,"floats(%f): %f %f %f %f %f\n",loc_amt,fnA,fnC,fnG,fnT,fnN);
	*/

	for(unsigned int i=0; i<NUM_SNP_VALS; i++) {
		newC = CHARFROMFLT(loc_amt, FLTFROM255(prev_amt, reads[reads_loc+i]) + amt[i]*scale);
		reads[reads_loc+i] = newC;
	}

	//fprintf(stderr,"post: %hu %hu %hu %hu %hu\n",reads[reads_loc+A_POS],reads[reads_loc+C_POS],reads[reads_loc+G_POS],reads[reads_loc+T_POS],reads[reads_loc+N_POS]);
#else
	reads[A_POS][loc] += amt[A_POS] * scale;
	reads[C_POS][loc] += amt[C_POS] * scale;
	reads[G_POS][loc] += amt[G_POS] * scale;
	reads[T_POS][loc] += amt[T_POS] * scale;

	if(NUM_SNP_VALS > N_POS)
		reads[N_POS][loc] += amt[N_POS] * scale;
#ifdef _INDEL
	if(NUM_SNP_VALS > DEL_POS)
		reads[DEL_POS][loc] += amt[DEL_POS] * scale;
	if(NUM_SNP_VALS > INS_POS)
		reads[INS_POS][loc] += amt[INS_POS] * scale;
#endif // _INDEL
#endif // DISCRETIZE/INTDISC
}

/* 
 * Invariant:  We must call "AddScore" before calling this function
 */
void Genome::AddSeqScore(unsigned long pos, const float amt, unsigned int which) {
	unsigned int loc = (pos-my_start)/gGEN_SIZE;
	//fprintf(stderr,"pos: %lu, myStart: %lu, loc:%u \n",pos,my_start,loc);

#if defined( DISCRETIZE )
	// If these haven't been instantiated, don't make any changes
	if(!read_allot)
		return;

	int i;
	//fprintf(stderr,"At loc %u, read_allot is %d\n",loc,read_allot[loc]);
	center_t *map_c = centers[read_allot[loc]];
	center_t gen_amt[DIM_CENTERS];
	for(i=0; i<DIM_CENTERS; i++)
		gen_amt[i] = map_c[i]*(amount_genome[loc]-amt);
	
	gen_amt[which] += amt;
	for(i=0; i<DIM_CENTERS; i++)
		gen_amt[i] /= amount_genome[loc];

	read_allot[loc] = find_center(gen_amt);
#elif defined( INTDISC )
	if(!reads)
		return;

	float loc_amt = amount_genome[loc];
	float prev_amt = loc_amt-amt;
	unsigned int reads_loc = loc*NUM_SNP_VALS;
	unsigned char newC;
	
	for(unsigned int i=0; i<NUM_SNP_VALS; i++) {
		if(i == which) {
			newC = CHARFROMFLT(loc_amt, FLTFROM255(prev_amt, reads[reads_loc+which]) + amt);
			reads[reads_loc+which] = newC;
		}
		else {
			newC = CHARFROMFLT(loc_amt, FLTFROM255(prev_amt, reads[reads_loc+i]));
			reads[reads_loc+i] = newC;
		}
	}
#else
	if(which < NUM_CHAR_VALS && reads[which])
		reads[which][loc] += amt;
#endif	// DISCRETIZE/INTDISC

}

string Genome::GetPos(const pair<unsigned long,int> p_pos) {
	return GetPos(p_pos.first,p_pos.second);
}
string Genome::GetPos(const unsigned long &pos, const int strand) {
	int chr_num = 0;

	while(names[chr_num].second <= pos && ((unsigned int)chr_num < names.size()) )
		chr_num++;
		
	ostringstream to_return;
	
	chr_num--;
	to_return << names[chr_num].first << ":";
	if(strand == POS_STRAND)
		to_return << "+";
	else
		to_return << "-";
	to_return << (pos - names[chr_num].second);
	
	return to_return.str();
}

pair<string,unsigned long> Genome::GetPosPair(const unsigned long pos) {
	
	int chr_num = 0;

	while(names[chr_num].second <= pos && ((unsigned int)chr_num < names.size()) )
		chr_num++;
		
	chr_num--;

	string chr_name = names[chr_num].first;
	unsigned long chr_pos = (pos - names[chr_num].second);
	
	return pair<string,unsigned long>(chr_name,chr_pos);
}

unsigned long Genome::GetAbsolutePosition(const char* chr, const unsigned long &pos) {
	unsigned long absPos = 0;
	bool found = false;
	
	//fprintf(stderr,"Size of names is %u\n",names.size());
	for(unsigned int i=0; i<names.size(); i++) {
		if(names[i].first.compare(chr) == 0) {
			absPos = names[i].second+pos;
			found = true;
		}
	}
	
	if(!found) {
		fprintf(stderr,"Chromosome is >%s<\n",chr);
		string excp = "Could not find chromosome in Genome file: ";
		excp += chr;
		throw new Exception(excp);
	}

	return absPos;
}

string Genome::GetStringRelative(const char* chr, const unsigned long &pos, const unsigned int length) {
	unsigned long abs_pos = GetAbsolutePosition(chr,pos);
	return GetString(abs_pos,length);
}


#ifdef SET_POS
bool Genome::getSetThreadID(long unsigned int pos, unsigned int tid) {
	if(pos < my_start)
		return false;

	//fprintf(stderr,"pos=%lu, genome_size=%lu, my_start=%lu, final_pos=%lu\n",pos,genome_size,my_start,pos-my_start);
	pos -= my_start;

	assert(pos < genome_size);

	if(gs_positions[pos/GEN_PACK] & (1<<(tid%8)))
		return true;

	unsigned int new_tid = 1<<(tid%8);
	gs_positions[pos/GEN_PACK] |= new_tid;
	return false;
}

void Genome::unsetThreadID(long unsigned int pos, unsigned int tid) {
	if(pos < my_start)
		return;

	pos -= my_start;


	assert(pos < genome_size);
	// We'll clear all the threads at this position
	gs_positions[pos/GEN_PACK] = 0;
}
#endif

unsigned int Genome::max_pos(float array[], unsigned int array_size) {
	float* pos = max_element(array, array+array_size); \
	unsigned int max_pos = pos - &array[0]; \
	return max_pos; \
}

const unsigned int MIN_COVERAGE = 5;
const float MIN_RATIO = 1.5;
const float MAX_PVAL = 0.01;
const float MAX_RATIO = 0.0;

double Genome::LRT(float chars[5], int &snp_pos) {
	snp_pos = max_pos(chars, 5);

	double sum = chars[0] + chars[1] + chars[2] + chars[3] + chars[4];

	double ratio = pow(.2,sum) 
					/ ( pow(chars[snp_pos]/sum,(double)chars[snp_pos]) 
					    * pow((sum-chars[snp_pos])/sum/4,sum-chars[snp_pos]) );
	
	//double pval = 1-pchisq(-2*log(ratio),1);
	double pval = 1-gsl_cdf_chisq_P(-2*log(ratio),1);
	
	return pval;
	
}

#define MONO_DIP_RATIO	3.0f
#define DIFF_UNIF_PRIOR	0.2

double Genome::dipLRT(float chars[NUM_CHAR_VALS], int &snp_pos1, int &snp_pos2, bool &dip) {

	/*******************************************
	 * If the ratio is too great, we're just going 
	 * to define it as a monoploid SNP
	 *******************************************/
	snp_pos1 = max_pos(chars, NUM_CHAR_VALS);
		
	double ratio1,ratio2;
	double sum = 0;
	for(unsigned int i=0; i<NUM_CHAR_VALS; i++) {
		sum += chars[i];
	}

	ratio1 = pow(.2,sum) 
					/ (pow(chars[snp_pos1]/sum,(double)chars[snp_pos1]) 
					  * pow(((sum-chars[snp_pos1])/sum)/(NUM_CHAR_VALS-1),sum-chars[snp_pos1]));
					
	double pval1 = 1-gsl_cdf_chisq_P(-2*log(ratio1),1);
	
	// Find the second highest SNP
	float second[NUM_CHAR_VALS];
	for(unsigned int i=0; i<NUM_CHAR_VALS; i++) {
		second[i] = chars[i];
	}
	second[snp_pos1] = 0;
	snp_pos2 = max_pos(second, NUM_CHAR_VALS);
	
	double pval2;
	if(chars[snp_pos1]/chars[snp_pos2] > MONO_DIP_RATIO || snp_pos1 == snp_pos2) {
		pval2 = MAX_PVAL;
		ratio2 = MAX_RATIO;
		snp_pos2 = -1;
	}
	else {
		// Add in the uniform prior
		for(unsigned int i=0; i<NUM_CHAR_VALS; i++) {
			chars[i] += DIFF_UNIF_PRIOR;
			sum += DIFF_UNIF_PRIOR;
		}
		
		// Calculate a new ratio and pval with the uniform prior
		ratio1 = pow(.2,sum) 
					/ (pow(chars[snp_pos1]/sum,(double)chars[snp_pos1]) 
					  * pow(((sum-chars[snp_pos1])/sum)/(NUM_CHAR_VALS-1),sum-chars[snp_pos1]));
		pval1 = 1-gsl_cdf_chisq_P(-2*log(ratio1),1);

		// Calculate the diploid ratio and prior
		//rat = (.2)^sum/((x[m]/sum)^(x[m])*(x[m2]/sum)^(x[m2])
		//		* (1-(x[m]+x[m2])/sum)^(sum-x[m]-x[m2]))
		double numerator = pow(.2,sum);

		double denominator = ( pow(chars[snp_pos1]/sum, (double)chars[snp_pos1])
				    * pow(chars[snp_pos2]/sum, (double)chars[snp_pos2])
				    * pow(((sum - chars[snp_pos1] + chars[snp_pos2])/sum)/3, 
				    	  sum - chars[snp_pos1] - chars[snp_pos2]));	
		ratio2 = numerator/denominator;

#ifdef DEBUG
		fprintf(stderr,"x is %f %f %f %f %f %f %f\n",chars[0],chars[1],chars[2],chars[3],chars[4],chars[5],chars[6]);
		fprintf(stderr,"ratio2 is %f (.2^%f/denom=%f/%f)\n",ratio2,sum,numerator,denominator);
		fprintf(stderr,"Result from gsl is %e\n",gsl_cdf_chisq_P(-2*ratio2, 2));
#endif


		//pval2 = 1-pchisq(-2*log(ratio2), 2);
		pval2 = 1-gsl_cdf_chisq_P(-2*log(ratio2), 2);

		// Remove the uniform prior
		for(unsigned int i=0; i<NUM_CHAR_VALS; i++) {
			chars[i] -= DIFF_UNIF_PRIOR;
			sum -= DIFF_UNIF_PRIOR;
		}

	}
	
#ifdef DEBUG
	fprintf(stderr,"pval1 is %f and pval2 is %f\n",pval1,pval2);
	fprintf(stderr,"Ratio is %f/%f=%f\n",chars[snp_pos1],chars[snp_pos2],chars[snp_pos1]/chars[snp_pos2]);
#endif
	
	// They're too high for a chi-square distribution, use the ratio
	if(pval2 == 0 && pval1 == 0) {
		if(ratio2 < ratio1 && (chars[snp_pos1]/chars[snp_pos2] < MONO_DIP_RATIO)) {
			dip = true;
#ifdef DEBUG
			fprintf(stderr,"Result: Diploid.  Pval is %f\n\n",pval2);
#endif
		}
		else {
			dip = false;
#ifdef DEBUG
			fprintf(stderr,"Result: Monoploid.  Pval is %f\n\n",pval2);
#endif
		}

		return 0;
	}

	if(pval2 < pval1 && (chars[snp_pos1]/chars[snp_pos2] < MONO_DIP_RATIO)) {
	//if(ratio2 < ratio1 && (chars[snp_pos1]/chars[snp_pos2] < MONO_DIP_RATIO)) {
#ifdef DEBUG
		fprintf(stderr,"Result: Diploid.  Pval is %f\n\n",pval2);
		fprintf(stderr,"Y pval2: %e<pval1: %e, chars[1]=%f,chars[2]=%f\n",pval2,pval1,chars[snp_pos1],chars[snp_pos2]);
#endif
		dip = true;
		return pval2;
	}
	else {
#ifdef DEBUG
		fprintf(stderr,"Result: Monoploid.  Pval is %f\n\n",pval2);
#endif
		dip = false;
		return pval1;
	}
}

double Genome::is_snp(float chars[NUM_CHAR_VALS], int &snp_pos1, int &snp_pos2, bool &dip) {
	
	/*
	//add a diffuse/uniform prior to chars
	for(int i=0; i<5; i++) {
		chars[i] += DIFF_UNIF_PRIOR;
	}
	*/

	double pval;
	if(gSNP_MONOP) {
		dip = false;
		pval = LRT(chars,snp_pos1);
	}
	else
		pval = dipLRT(chars,snp_pos1,snp_pos2,dip);

	/*
	// remove the prior from chars
	for(int i=0; i<5; i++) {
		chars[i] -= DIFF_UNIF_PRIOR;
	}
	*/

	return pval;

}

/**
 * is_indel will determine if this is an indel or not. 
 */
//inline int Genome::is_indel

void Genome::AppendFinal(const char* fn) {
	PrintFinal(fn, true);
}
void Genome::PrintFinal(const char* fn) {
	PrintFinal(fn, false);
}

void Genome::PrintFinal(const char* fn, bool append) {
	if(gPRINT_VCF) {
		// Not sure how to handle VCF and Bisulfite.  We'll let them deal with it.
		PrintFinalVCF(fn, append);
		return;
	}

	if(gSNP) {
		PrintFinalSNP(fn, append);
	}
	// This one is still sgrex, but not SNP calling
	else if(gBISULFITE || gATOG) {	// if any of the reads have been instantiated...
		PrintFinalBisulfite(fn, append);
	}
	else {
		PrintFinalSGR(fn, append);
	}
}

#define MIN_PRINT	0.001

void Genome::PrintFinalSNP(const char* fn, bool append) {
	// Instead of printing out a file for each, we want to print
	// out one file (dubbed the '.sgrex' file for 'sgr-extra') to be a tab-deliminated file
	// for the number of a's, c's, g's, t's, and gaps at each location in the genome.
	char gmp_fn[strlen(fn)+10];
	strcpy(gmp_fn,fn);
	strcat(gmp_fn,".gmp");
	fprintf(stderr, "SNP File is: %s\n",gmp_fn);
	
	FILE* gmp_file;
	if(append)
		gmp_file = fopen(gmp_fn,"a");
	else
		gmp_file = fopen(gmp_fn,"w");
		
	if(!gmp_file)
		perror("Could not open .gmp file");
		
	// don't output the Header
	//fprintf(sgrex_file,"chr\tpos\tsgr_val\tread_a\tread_c\tread_g\tread_t\tread_n\tread_gap\tsnp\n");
	//fprintf(sgrex_file,"#chr\tpos\tsgr_val\tread_a\tread_c\tread_g\tread_t\tread_n\tsnp\n");
	
	if(gVERBOSE > 1) 
		fprintf(stderr,"NAMES: %ld\n",names.size());
	
	unsigned long count=0;
	for(unsigned int j=0; j<names.size()-1; j++) {
	//for(unsigned int j=0,count=0; names[j].first != "end"; j++) {
		if(gVERBOSE > 1) 
			fprintf(stderr,"\t[%d]: %s  %lu\n",j+1,names[j+1].first.c_str(),names[j+1].second);
		for(; count<names[j+1].second; count+=gGEN_SIZE) {
			if(count < my_start && !read_all)
				continue;
			if(count > my_end && !read_all)
				break;
				
			unsigned long locus = (count-my_start)/gGEN_SIZE;

			//if(genome[locus].amount > (float)0.0) {	//if something matched here...
			if(amount_genome[locus] > MIN_PRINT) {	//if something matched here...
			
								 //  chr pos  sgr -- want 1-based indexing
				fprintf(gmp_file, "%s\t%ld\t%.5f", names[j].first.c_str(), (count-names[j].second)+1,
							//genome[locus].amount);
							amount_genome[locus]);

				// Print out each of the read positions
#if defined( DISCRETIZE )
				unsigned int dim;
				center_t *map_c = centers[read_allot[locus]];
				center_t read_amt[DIM_CENTERS];
				for(dim=0; dim<DIM_CENTERS; dim++) {
					read_amt[dim] = map_c[dim]*amount_genome[locus];
				}
				for(dim=0; dim<DIM_CENTERS; dim++) {
					fprintf(gmp_file,"\t%.5f",read_amt[dim]);
				}
#elif defined( INTDISC )
				for(unsigned int read_count=0; read_count<NUM_SNP_VALS; read_count++) {
					fprintf(gmp_file,"\t%.5f",reads[locus*NUM_SNP_VALS+read_count]/255
									*amount_genome[locus]);
				}
#else
				for(unsigned int read_count=0; read_count<NUM_SNP_VALS; read_count++) {
					fprintf(gmp_file,"\t%.5f",reads[read_count][locus]);
				}
				
#endif // DISCRETIZE/INTDISC

				PrintSNPCall(count, gmp_file);
				
				// Carriage return
				fprintf(gmp_file,"\n");
			}
		}
	}
	
	// Close the file
	fclose(gmp_file);
}

void Genome::PrintSNPCall(const unsigned long count, FILE* gmp_file) {
	// Get the position here
	char at_pos = GetChar(count);
	unsigned long locus = (count-my_start)/gGEN_SIZE;
	// Determine if it's a snp
#if defined( DISCRETIZE )
	center_t *map_c = centers[read_allot[locus]];
	float char_arr[DIM_CENTERS];
	for(unsigned int i=0; i<DIM_CENTERS; i++)
		char_arr[i] = map_c[i]*amount_genome[locus];
#elif defined( INTDISC )
	float char_arr[NUM_SNP_VALS];
	for(unsigned int i=0; i<NUM_SNP_VALS; i++) {
		char_arr[i] = reads[locus*NUM_SNP_VALS+i]/255
								*amount_genome[locus];
	}
#else
#ifdef _INDEL										
	float char_arr[NUM_SNP_VALS] =
						{	reads[0][locus],
							reads[1][locus],
							reads[2][locus],
							reads[3][locus],
							reads[4][locus],
							reads[5][locus],
							reads[6][locus]};
#else
	float char_arr[NUM_SNP_VALS] =
						{	reads[0][locus],
							reads[1][locus],
							reads[2][locus],
							reads[3][locus],
							reads[4][locus]};
#endif // _INDEL
#endif // DISCRETIZE/INTDISC

	// Sanity check
	float sum = 0;
	for(unsigned int i=0; i<NUM_SNP_VALS; i++) {
		sum += char_arr[i];
	}
	if(gVERBOSE > 1) {
		if(amount_genome[(count-my_start)/gGEN_SIZE] - sum > 0.01) {
			fprintf(stderr,"ERROR at %lu:  %f vs sum of %f\n",(count-my_start)/gGEN_SIZE,amount_genome[(count-my_start)/gGEN_SIZE],sum);
		}
	}

	int snp_pos1, snp_pos2;
	bool dip = false;
	double p_val = is_snp(char_arr, snp_pos1, snp_pos2, dip);


	// If the max pos is not what's in the genome
	// 	or it's a diploid snp
	if(snp_pos1 != g_gen_CONVERSION[(int)at_pos] || dip) {
		if( p_val < gSNP_PVAL) {
			if(dip) {
				fprintf(gmp_file,"\tY:%c->%c/%c p_val=%.2e", 
					at_pos, gINT2BASE[snp_pos1], gINT2BASE[snp_pos2], p_val);
			}
			else {
				fprintf(gmp_file,"\tY:%c->%c p_val=%.2e", 
					at_pos, gINT2BASE[snp_pos1], p_val);
			}
		}
		else
			if(dip) {
				fprintf(gmp_file,"\tN:%c->%c/%c p_val=%.2e", 
					at_pos, gINT2BASE[snp_pos1], gINT2BASE[snp_pos2], p_val);
			}
			else {
				fprintf(gmp_file,"\tN:%c->%c p_val=%.2e", 
					at_pos, gINT2BASE[snp_pos1], p_val);
			}
	}
	else
		//fprintf(gmp_file,"\tN:%c p_val=%.2e", at_pos, p_val);
		fprintf(gmp_file,"\tN");

}

void Genome::PrintFinalBisulfite(const char* fn, bool append) {
	// We got a new update from Evan.  Instead of printing out a file for each, we want to print
	// out one file (dubbed the '.gmp' file--previously '.sgrex') to be a tab-deliminated file
	// for the number of a's, c's, g's, t's, and gaps at each location in the genome.
	char gmp_fn[strlen(fn)+10];
	strcpy(gmp_fn,fn);
	strcat(gmp_fn,".gmp");
	if(gBISULFITE)
		printf("Bisulfite File is: %s\n",gmp_fn);
	if(gATOG)
		printf("RNA edit File is: %s\n",gmp_fn);

	FILE* gmp_file;
	if(append)
		gmp_file = fopen(gmp_fn,"a");
	else
		gmp_file = fopen(gmp_fn,"w");
		
	if(!gmp_file)
		perror("Could not open .gmp file");
	// don't output the Header
	//fprintf(gmp_file,"chr\tpos\tsgr_val\tread_a\tread_c\tread_g\tread_t\tread_n\tread_gap\n");
	//fprintf(gmp_file,"#chr\tpos\tsgr_val\tread_a\tread_c\tread_g\tread_t\tread_n\n");
	
	if(gVERBOSE > 1) 	
		fprintf(stderr,"Size of NAMES: %ld, size of gGEN_SIZE:%d\n",names.size(),gGEN_SIZE);
	
	double total_stuff = 0;
	
	unsigned long count=0;
	for(unsigned int j=0; j<names.size()-1; j++) {
		if(gVERBOSE > 1) 
			fprintf(stderr,"\t[%d/%lu]: %s\n",j,names.size(),names[j].first.c_str());
		for(; count<names[j+1].second; count+=gGEN_SIZE) {
			if(count < my_start && !read_all)
				continue;
			if(count > my_end && !read_all)
				break;

			unsigned long locus = (count-my_start)/gGEN_SIZE;

			char at_pos = GetChar(count);
			int gen_count_spot;	// Later, we'll ensure there's actually info here
			
			// Check if it's bisulfite or RNA editing we're doing
			if(gBISULFITE) {
				if(gMATCH_POS_STRAND && !gBISULFITE2) {
					// only print if there is a 'c' nucleotide here
					if(at_pos != 'c')
						continue;
					gen_count_spot = 1;
				}
				else {
					if(at_pos != 'g')
						continue;
					gen_count_spot = 2;
				}
			}
			else {	// gATOG
				if(gMATCH_POS_STRAND) {
					if(at_pos != 'a')
						continue;
					gen_count_spot = 0;
				}
				else {	// gMATCH_NEG_STRAND
					if(at_pos != 't')
						continue;
					gen_count_spot = 3;
				}
			}
			
			total_stuff += amount_genome[locus];
			
			//if(genome[locus].amount > (float)0.0) {	//if something matched here...
			if(amount_genome[locus] > 0.0f) {	//if something matched here...
			//if(reads[gen_count_spot][locus] > MIN_PRINT) {	//if something matched here...
								 //  chr pos  sgr -- want 1-based indexing
				fprintf(gmp_file, "%s\t%ld\t%f", names[j].first.c_str(), (count-names[j].second)+1,
							//genome[locus].amount);
							amount_genome[locus]);

				// Print out each of the read positions
#if defined( DISCRETIZE )
				unsigned int dim;
				center_t *map_c = centers[read_allot[locus]];
				center_t read_amt[DIM_CENTERS];
				for(dim=0; dim<DIM_CENTERS; dim++) {
					read_amt[dim] = map_c[dim]*amount_genome[locus];
				}
				for(dim=0; dim<DIM_CENTERS; dim++) {
					fprintf(gmp_file,"\t%.5f",read_amt[dim]);
				}
#elif defined( INTDISC )
				for(unsigned int read_count=0; read_count<NUM_CHAR_VALS; read_count++) {
					fprintf(gmp_file,"\t%.5f",
								reads[locus*NUM_SNP_VALS+read_count]/255 
								*amount_genome[locus]);
				}
#else
				// Only print out acgtn, not indel
				for(unsigned int read_count=0; read_count<NUM_CHAR_VALS; read_count++) {
					fprintf(gmp_file,"\t%.5f",reads[read_count][locus]);
				}

#endif // DISCRETIZE/INTDISC

				// Just for testing right now
				//PrintSNPCall(count,gmp_file);
				
				// Carriage return
				fprintf(gmp_file,"\n");
			}
		}
	}
	//fprintf(stderr,"Total stuff=%f\n",total_stuff);
	
	// Close the file
	fclose(gmp_file);
}

void Genome::PrintFinalVCF(const char* fn, bool append) {
	char* vcf_fn = new char[strlen(fn)+5];
	strcpy(vcf_fn, fn);
	strcat(vcf_fn, ".vcf");
	fprintf(stderr, "VCF file is %s\n", vcf_fn);
	
	FILE* vcf_file;
	if(append) {
		vcf_file = fopen(vcf_fn, "a");
	}
	else {
		// If we're not appending, we need to write the header lines
		vcf_file = fopen(vcf_fn, "w");
		fprintf(stderr, "Just opened vcf file %s.  Successful?%c\n", vcf_fn, vcf_file?'Y':'N');
		fprintf(vcf_file, "##fileformat=VCFv4.0\n");
		time_t rawtime;
		time ( &rawtime );
		fprintf(vcf_file, "##fileDate=%s\n", ctime (&rawtime));
		fprintf(vcf_file, "##source=GNUMAP v"gVERSION"\n");
		fprintf(vcf_file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
		fprintf(stderr, "Finished writing header to VCF file\n");
	}

	if(!vcf_file)
		perror("Could not open .vcf file!");

	unsigned long snp_count=0;
	unsigned long count=0;
	for(unsigned int j=0; j<names.size(); j++) {
		fprintf(stderr, "Printing out for name %s (my_start=%lu,my_end=%lu,read_all=%c)\n", names[j].first.c_str(), my_start,my_end,read_all?'y':'n');
		for(; count<names[j+1].second; count+=gGEN_SIZE) {
			if(count < my_start && !read_all)
				continue;
			if(count > my_end && !read_all) {
				fprintf(stderr, "\t...skipping\n");
				break;
			}

			unsigned long locus = (count-my_start)/gGEN_SIZE;
			if(amount_genome[locus] <= 0.0f)
				continue;

			char at_pos = GetChar(count);
#if defined( DISCRETIZE )
			center_t *map_c = centers[read_allot[locus]];
			float char_arr[DIM_CENTERS];
			for(unsigned int i=0; i<DIM_CENTERS; i++)
				char_arr[i] = (float)map_c[i]*amount_genome[locus];
#elif defined( INTDISC )
			
			float char_arr[NUM_SNP_VALS];
			for(unsigned int i=0; i<NUM_SNP_VALS; i++) {
				char_arr[i] = (float)reads[locus*NUM_SNP_VALS+i]/255
										*amount_genome[locus];
			}
#else
#ifdef _INDEL										
			float char_arr[NUM_SNP_VALS] =
						{	reads[0][locus],
							reads[1][locus],
							reads[2][locus],
							reads[3][locus],
							reads[4][locus],
							reads[5][locus],
							reads[6][locus]};
#else
			float char_arr[NUM_SNP_VALS] =
						{	reads[0][locus],
							reads[1][locus],
							reads[2][locus],
							reads[3][locus],
							reads[4][locus]};
#endif // _INDEL
#endif // DISCRETIZE/INTDISC

			int snp_pos1, snp_pos2;
			bool dip = false;
			double p_val = is_snp(char_arr, snp_pos1, snp_pos2, dip);


			// If the max pos is not what's in the genome
			// 	or it's a diploid snp
			if(snp_pos1 != g_gen_CONVERSION[(int)at_pos] || dip) {
				if( p_val < gSNP_PVAL) {
					if(dip) {
						fprintf(vcf_file,"%s\t%u\tsnp%u\t%c\t%c%c\t.\t.\tDiploid;pval=%.5f;coverage=%.5f;ratio=%.2f\n",
								names[j].first.c_str(), (count-names[j].second)+1, snp_count++, 
								at_pos, gINT2BASE[snp_pos1],gINT2BASE[snp_pos2], 
								p_val, amount_genome[locus], char_arr[snp_pos2]/char_arr[snp_pos1]);
					}
					else {
						fprintf(vcf_file,"%s\t%u\tsnp%u\t%c\t%c\t.\t.\tMonoploid;pval=%.5f;coverage=%.5f\n",
								names[j].first.c_str(), (count-names[j].second)+1, snp_count++, 
								at_pos, gINT2BASE[snp_pos1],
								p_val, amount_genome[locus]);
					}
				}
			}


		}
	}
	
}

void Genome::PrintFinalSGR(const char* fn, bool append) {
	char* filename = new char[strlen(fn)+5];
	sprintf(filename,"%s.sgr",fn);
	
	FILE* sgr_file;
	if(append)
		sgr_file = fopen(filename,"a");
	else
		sgr_file = fopen(filename,"w");
	
	//fprintf(stderr,"Printing for %u genomes to %lu\n",names.size(),names[names.size()-1].second);
	unsigned long count=0;
	for(unsigned int i=0; i<names.size()-1; i++) {
		for(; count<names[i+1].second; count+=gGEN_SIZE) {
			if(count < my_start && !read_all)
				continue;
			if(count > my_end && !read_all)
				break;

			//if(genome[(count-my_start)/gGEN_SIZE].amount> 0)
			if(amount_genome[(count-my_start)/gGEN_SIZE] > MIN_PRINT) {
				//fprintf(stderr,"count=%u, start=%lu, end=%lu/%u=%lu, size=%lu, at: %f\n",count,my_start,my_end,gGEN_SIZE,my_end/gGEN_SIZE,genome_size,amount_genome[(count-my_start)/gGEN_SIZE]);
				
									// chr pos  sgr -- want 1-based indexing
					fprintf(sgr_file, "%s\t%ld\t%.5f\n",names[i].first.c_str(),(count-names[i].second)+1,
					amount_genome[(count-my_start)/gGEN_SIZE]);
			}
		}
	}
	
	fclose(sgr_file);
	printf("Finished printing to %s\n",filename);
	delete[] filename;
}

unsigned long Genome::size() {
	return genome_size;
}

