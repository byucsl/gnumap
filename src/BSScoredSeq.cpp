/* BSScoredSeq.cpp
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

#include "BSScoredSeq.h"

void BSScoredSeq::score(double denom, Genome &gen, unsigned int len, Read &read, pthread_mutex_t &lock) {

	if(!positions.size())	//if there were no matches, just return.
		return;

	double total_score = log_align_score / denom;
	//set<pair<unsigned int,bool> >::iterator sit;

	//we only need to compute this alignment string once because, while there may be many locations for each
	//sequence, we're ensuring that the sequences are the exact same.

	bin_seq bs;	// global one

	string read_consensus;
	string aligned, rev_comp;
	string gen_string = sequence;
	
	// If it's on the revcomp strand
	if(firstStrand == NEG_STRAND) {
		//fprintf(stderr,"Is on the negative strand! %s\n",__FILE__);
		// Copy it so we don't need to reverse it again
		float** rev_pwm = reverse_comp_cpy(read.pwm,read.length);
		Read rc(rev_pwm,read.length);
		rc.name = read.name;
		read_consensus = GetConsensus(rc);
		
        // TODO: I think that this is broken too, it is using the outdated 
        // bs.get_align_score_w_traceback method, does this need to be fixed?
		aligned = bs.get_align_score_w_traceback(rc, read_consensus, gen_string).first;
		
		
		// Clean up the copy
		for(unsigned int i=0; i<rc.length; i++) {
			delete[] rev_pwm[i];
		}
		delete[] rev_pwm;
	}
	else {
		read_consensus = GetConsensus(read);
		aligned = bs.get_align_score_w_traceback(read, read_consensus, gen_string).first;
	}

	rev_comp = reverse_comp(aligned);
	
	//fprintf(stderr,"Aligning %s with %s, result %s\n",gen_string.c_str(),read_consensus.c_str(),aligned.c_str());

	set<pair<unsigned long,int> >::iterator sit;
	
	MUTEX_LOCK(&lock);
	for(sit=positions.begin(); sit!=positions.end(); sit++) {
		if(sit->second == firstStrand) {
			for(unsigned int i=0; i<aligned.size(); i++) {
				gen.AddScore( (*sit).first+i,total_score );
				gen.AddSeqScore( ( *sit ).first + i, total_score, g_gen_CONVERSION[ ( int ) aligned[ i ] ] );
			}
		}
		else {
			for(unsigned int i=0; i<rev_comp.size(); i++) {
				gen.AddScore( (*sit).first+i,total_score );
				gen.AddSeqScore((*sit).first+i,total_score,g_gen_CONVERSION[(int)rev_comp[i]]);
			}
		}
	}
	MUTEX_UNLOCK(&lock);
}

inline string BSScoredSeq::GetConsensus(const Read &read) const {
	string seq = "";

	for(unsigned int i=0; i<read.length; i++) {
		seq += max_char(read.pwm[i]);
	}

	return seq;
}

/*
 * Because this is used to produce the consensus sequence, we need ambiguity characters...
 */
//inline char max_char(vector<double> &chr) {
inline char BSScoredSeq::max_char(const float* chr) const {
	if((chr[0] == chr[1]) && (chr[0] == chr[2]) && (chr[0] == chr[3]))
		return 'n';
	if(chr[0] >= chr[1]) {
		if(chr[0] >= chr[2]) {
			if(chr[0] >= chr[3])
				return 'a';
			else
				return 't';
			}
		else {
			if(chr[2] >= chr[3])
				return 'g';
			else
				return 't';
		}
	}
	else {
		if(chr[1] >= chr[2]) {
			if(chr[1] >= chr[3])
				return 'c';
			else
				return 't';
		}
		else {
			if(chr[2] >= chr[3])
				return 'g';
			else return 't';
		}
	}

	return 'n';
}

