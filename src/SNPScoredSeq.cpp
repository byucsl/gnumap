/* SNPScoredSeq.cpp
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

#include "SNPScoredSeq.h"

#ifndef SNP_NW
void SNPScoredSeq::score(double denom, Genome &gen, unsigned int len, Read &read, 
						pthread_mutex_t &lock) {

	if(!positions.size())	//if there were no matches, just return.
		return;

	double total_score = log_align_score / denom;

	// we only need to compute this alignment string once because, while there may be many 
	// locations for each sequence, we're ensuring that the sequences are the exact same.

	bin_seq bs;

	string gen_string = sequence;
	string read_consensus;
	float** hmm_aligned;

	if(firstStrand == NEG_STRAND) {
		//fprintf(stderr,"Is on the negative strand! %s\n",__FILE__);
		// Copy it so we don't need to reverse it again
		float** rev_pwm = reverse_comp_cpy(read.pwm,read.length);
		Read rc(rev_pwm,read.length);
		rc.name = read.name;
		read_consensus = GetConsensus(rc);
		
		hmm_aligned = bs.pairHMM(rc, read_consensus, gen_string);
		
		// Clean up the copy
		for(unsigned int i=0; i<rc.length; i++) {
			delete[] rev_pwm[i];
		}
		delete[] rev_pwm;
	}
	else {
		read_consensus = GetConsensus(read);
		hmm_aligned = bs.pairHMM(read, read_consensus, gen_string);
	}

	// Need to create the reverse compliment as well
	float** hmm_aligned_rev = reverse_comp_cpy_phmm(hmm_aligned, read.length);

	//fprintf(stderr,"Aligning (size=%u) %s vs %s\n",positions.size(),read_consensus.c_str(),gen_string.c_str());
	
	set<pair<unsigned long,int> >::iterator sit;
	
	MUTEX_LOCK(&lock);
	for(sit=positions.begin(); sit!=positions.end(); sit++) {
		// If it's the same strand as the first, we'll just add it to the genome normally.
		// Otherwise, we need to take the reverse-complement of the phmm sequence and add 
		// it to the genome that way
		if(sit->second == firstStrand) {
			for(unsigned int i=0; i<gen_string.size(); i++) {
				gen.AddScore( (*sit).first+i,total_score );
				//pair<string,unsigned long> chrspot = gen.GetPosPair((*sit).first);

				gen.AddSeqScore((*sit).first+i, hmm_aligned[i], total_score);
			}
		}
		else {
			for(unsigned int i=0; i<gen_string.size(); i++) {
				gen.AddScore( (*sit).first+i,total_score );
				//pair<string,unsigned long> chrspot = gen.GetPosPair((*sit).first);

				gen.AddSeqScore((*sit).first+i, hmm_aligned_rev[i], total_score);
			}
		}
	}
	MUTEX_UNLOCK(&lock);
	
	// Free HMM arrays
	freeHMMSequence(hmm_aligned,read.length);
	freeHMMSequence(hmm_aligned_rev,read.length);
}

#else
/* 
 *This is the version of SCORE that doesn't use a pHMM.  Just to test.
 */
void SNPScoredSeq::score(double denom, Genome &gen, unsigned int len, Read &read, pthread_mutex_t &lock) {

	if(!positions.size())	//if there were no matches, just return.
		return;

	double total_score = log_align_score / denom;
	//set<pair<unsigned int,bool> >::iterator sit;

	//we only need to compute this alignment string once because, while there may be many locations for each
	//sequence, we're ensuring that the sequences are the exact same.

	bin_seq bs;	// global one
	unsigned int seq_size = sequence.size();

	string read_consensus = GetConsensus(read);

	string gen_string = gen.GetString((*positions.begin()).first, seq_size);
	string aligned = bs.get_align_score_w_traceback(read, read_consensus, gen_string).first;

	set<pair<unsigned long,int> >::iterator sit;

	MUTEX_LOCK(&lock);
	for(sit=positions.begin(); sit!=positions.end(); sit++) {
		if(sit->second == firstStrand) {
		for(unsigned int i=0; i<aligned.size(); i++) {
			gen.AddScore( (*sit).first+i,total_score );

			gen.AddSeqScore((*sit).first+i,total_score,g_gen_CONVERSION[(unsigned int)aligned[i]]);
		}
		}
		else {
		for(unsigned int i=0; i<aligned.size(); i++) {
			gen.AddScore( (*sit).first+i,total_score );

			gen.AddSeqScore((*sit).first+i,total_score,g_gen_CONVERSION[(unsigned int)aligned[i]]);
		}

		}
	}
	MUTEX_UNLOCK(&lock);
}
#endif
