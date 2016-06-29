/* ScoredSeq.h
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

#ifndef SCORED_SEQ_H
#define SCORED_SEQ_H

#include <string>
#include <set>
#include <cmath>
#include <pthread.h>
#include "const_include.h"
#include "SequenceOperations.h"
#include "Genome.h"

class ScoredSeq {
	protected:
		static const bool usePairHMM = false;
		
		void freeHMMSequence(float** phmmSeq, unsigned int length) {
			for(unsigned int i=0; i<length; i++) {
				delete[] phmmSeq[i];
			}
			delete[] phmmSeq;
		}
		
		
		/* Return the current time in seconds, using a float precision number. */
		double When() {
			struct timeval tp;
			gettimeofday(&tp, NULL);
			return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
		}

		// We want to know if the read matched at this location is the positive or 
		// negative strand.  If it is, all the matching genomic locations are backwards.
		int firstStrand;


		inline string GetConsensus(const Read &read) const {
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
		inline char max_char(const float* chr) const {
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

	public:
		string sequence;			// the genomic string this sequence matched to
		double align_score;			// the total (raw) alignment score of this sequence
		double log_align_score;		// logged alignment score
		
		//the int in this pair represents which strand:
		//	0 (POS_STRAND)_ = positive strand (+)
		//	1 (NEG_STRAND) = complementary strand (-)
		set<pair<unsigned long,int> > positions;

		ScoredSeq() : firstStrand(NEG_STRAND), sequence("<empty>"), align_score(-1),
			log_align_score(exp(align_score)) {
			positions.clear();
		}
	
		ScoredSeq(string seq, double as, unsigned long pos, int strand) :
			align_score(as) {

			if(strand == NEG_STRAND)  {
				// already been done...
				//sequence = reverse_comp(seq);
				sequence = seq;
			}
			else {
				//fprintf(stderr,"Strand is positive\n");
				sequence = seq;
			}

			log_align_score = exp(align_score);
			positions.insert(pair<unsigned long,int>(pos,strand));
			firstStrand = strand;
		}

		ScoredSeq(const ScoredSeq &other) {
			Init(other);
		}
		
		virtual ScoredSeq& operator =(const ScoredSeq &other) {
			//if( !(this == other) )
			Init(other);
			return *this;
		}
		
		virtual ~ScoredSeq() { }
		
		inline void Init(const ScoredSeq &other) {
			positions.clear();
			sequence = other.sequence;
			align_score = other.align_score;
			log_align_score = other.log_align_score;
			//positions = other.positions;
			for(set<pair<unsigned long,int> >::iterator sit=other.positions.begin(); 
					sit != other.positions.end(); ++sit)
				positions.insert(pair<unsigned long,int>((*sit).first,(*sit).second));
			
			firstStrand = other.firstStrand;
		}
		
		inline int compare(ScoredSeq *a, ScoredSeq *b) {
			return a->sequence.compare(b->sequence);
			/*
			if(a == b)
				return 0;
			else {
				if(a > b)
					return 1;
				else
					return -1;
			}
			*/
		}
		inline bool operator() (ScoredSeq &a, ScoredSeq &b) {
			//printf("ADDR: Comparing %s with %s\n",a.sequence.c_str(),b.sequence.c_str());
			return a.sequence.compare(b.sequence);
		}
		
		//The return value is 1 if the set doesn't contain an element whose sort key matches the parameter key. 
		//0 if the set does contain an element with a matching key.
		inline bool operator() (const ScoredSeq *a, const ScoredSeq *b) const {
			//printf("Comparing %s with %s\n",a->sequence.c_str(),b->sequence.c_str());
			return (a->sequence != b->sequence);
		}
		inline bool operator() (const double &a, const float &b) {
			return a != b;
		}
		inline bool operator == (ScoredSeq &other) {
			if( (sequence == other.sequence) )
				return true;
			return false;
		}
		inline bool operator == (ScoredSeq *other) {
			if( (sequence == other->sequence) )
				return true;
			return false;
		}
		inline bool operator< (const ScoredSeq &other) const {
			//printf("< adr Comparing %s with %s\n",sequence.c_str(),other.sequence.c_str());
			//return align_score < other.align_score;
			return (sequence < other.sequence);
		}
		inline bool operator< (const ScoredSeq* other) const {
			//printf("< ptr Comparing %s with %s\n",sequence.c_str(),other->sequence.c_str());
			//return align_score < other->align_score;
			return (sequence < other->sequence);
		}
		inline bool operator> (const ScoredSeq &other) const {
			//printf("> adr Comparing %s with %s\n",sequence.c_str(),other.sequence.c_str());
			double my_total = log_align_score / positions.size();
			double other_total = other.log_align_score / other.positions.size();
			return my_total > other_total;
		}
		inline bool operator> (const ScoredSeq *other) const {
			//printf("> ptr Comparing %s with %s\n",sequence.c_str(),other->sequence.c_str());
			double my_total = log_align_score / positions.size();
			double other_total = other->log_align_score / other->positions.size();
			return my_total > other_total;
		}
		
		inline bool is_greater(const ScoredSeq &other) {
			//printf("Logs: %f=%e %f=%e\n",other.align_score,exp(other.align_score),align_score,exp(align_score));
			double my_total = log_align_score;
			double other_total = other.log_align_score;
			return my_total > other_total;
		}
		
		double get_score() {
			return align_score;
		}
		double get_log_align_score() {
			return log_align_score;
		}
		
		bool add_spot(unsigned long pos, int strand) {
			//cout << "inserting: " << pos << endl;
			if( positions.insert( pair<unsigned long,int>(pos,strand) ).second ) {
#ifdef DEBUG
				cerr << "inserted. \tsize: " << positions.size() << endl;
#endif
				return true;
			}
			else {
#ifdef DEBUG
				cerr << "already in set. \tsize: " << positions.size() << endl;
#endif
				return false;
			}
		}
		
		/*
		 * score will add the log-odds score (adjusted by denom) to each
		 * position on the Genome gen.  Additionally, this can be passed
		 * a flag to specify the number of positions after the start position
		 * to add to (for adding to different lengths of sequences).
		 */		 
		virtual void score(double denom, Genome &gen, unsigned int len, Read &read, pthread_mutex_t &lock) = 0;
		 
		/**
		 * Since the SAM format took just a little while for me to understand, I'll explain it
		 * here.  According to specs, each sequence is as follows:
		 * 
		 * <QNAME> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL> <
		 * 
		 * - For GNUMAP, the QNAME will be used if it is a fastq or fasta string--otherwise
		 *   the query name will be "seq" followed by the sequence number in the file.
		 * - The only bits in FLAG that are important are 0x0010 (where 1 represents an alignment
		 *   to the forward strand and 0 represents reverse)
		 * - RNAME is the chromosome mapped to
		 * - POS is the position on the chromosome
		 * - MAPQ is the mapping quality, where MAPQ = 10 * log_10(1-p(x)) where p(x) is GNUMAP's
		 *   posterior probability of mapping to that specific location
		 * - CIGAR is the alignment differences, where I=>Insertion, D=>Deletion and M=>Mismatch
		 *   on the genomic strand
		 * - MRNM will be '*' for all sequences (it's the mate-pair name)
		 * - MPOS is ignored.  Always '0'
		 * - ISIZE is ignored.  Always '0'
		 * - SEQ is the sequence
		 * - QUAL is the Phred-based quality of the sequences
		 * - OPT is a variable optional field of the format TAG:VTYPE:VALUE
		 * 		for GNUMAP, we'll use the format XA:f:<float> to represent a TAG with name
		 *		XA, the type is a float and the value is the alignment score.
		 * 
		 * @param denom The denominator for this particular read
		 * @param read The GNUMAP-style Read for this sequence
		 * @param consensus The consensus sequence for this Read
		 * @param gen The Genome that has been used for this program
		 * 
		 * @return A string in the SAM format
		 */
		vector<TopReadOutput> get_SAM(double denom, Read &read, string consensus, Genome &gen, unsigned int rIndex) {
			bin_seq bs;

			vector<TopReadOutput> ret_vec;
			if(positions.size() == 0)
				return ret_vec;
				
			double total_score = log_align_score / denom;
			
			// Mapping quality is defined as -10 * log_10(1-p)
			int MAPQ;
			if(total_score == 1)	// it was the only match.  Highest probability of .999
				MAPQ = 30;
			else
				MAPQ = (int)round(-10 * log(1-total_score)/log(10));
			if(MAPQ > 30)
				MAPQ = 30;	// only allow a max of 30
			
			string QUAL = str2qual(read);
			//fprintf(stderr,"QUAL for %s is %s\n",read.name,QUAL.c_str());

			string gen_string = sequence;
			pair<string,string> aligned;
			//if((*positions.begin()).second == NEG_STRAND)
			//	gen_string = reverse_comp(sequence);
			if(firstStrand == NEG_STRAND) {
				float** rev_pwm = reverse_comp_cpy(read.pwm,read.length);
				Read rc(rev_pwm,read.length);
				rc.name = read.name;
				string read_consensus = GetConsensus(rc);
				
                if( gNW )
                {
				    aligned = bs.get_align_score_w_traceback(rc,read_consensus,gen_string);
                }
                else
                {
                    aligned.first = read_consensus;
                    aligned.second = string();
                }

				//fprintf(stderr,"[%d] Aligning %s with %s, result %s\n",firstStrand,gen_string.c_str(),read_consensus.c_str(),aligned.first.c_str());
		
				// Clean up the copy
				for(unsigned int i=0; i<rc.length; i++) {
					delete[] rev_pwm[i];
				}
				delete[] rev_pwm;
			}
			else {
                if( gNW )
                {
                    aligned = bs.get_align_score_w_traceback(read,consensus,gen_string);
                }
                else
                {
                    aligned.first = consensus;
                    aligned.second = string();
                }

				//fprintf(stderr,"[%d] Aligning %s with %s, result %s\n",firstStrand,gen_string.c_str(),consensus.c_str(),aligned.first.c_str());
			}


			string CIGAR = aligned.second;
			if(CIGAR.length() == 0)
            {
                if( gNW )
                {
				    CIGAR = "*";
                }
                else
                {
                    CIGAR = std::to_string( consensus.size() ) + "M";
                }
            }
			else
            {
				fix_CIGAR_for_deletions(CIGAR);
			}

			// Print all the locations
			set<pair<unsigned long, int> >::iterator sit;
			for(sit = positions.begin(); sit != positions.end(); sit++) {
				pair<string,unsigned long> seq_pos = gen.GetPosPair((*sit).first);

				TopReadOutput tro;
				strncpy(tro.READ_NAME,read.name,MAX_NAME_SZ-1);
				tro.READ_NAME[MAX_NAME_SZ-1] = '\0';
				strncpy(tro.CHR_NAME,seq_pos.first.c_str(),MAX_NAME_SZ-1);
				tro.CHR_NAME[MAX_NAME_SZ-1] = '\0';
				// We're off by one (zero base), so increment it
				tro.CHR_POS = seq_pos.second+1;
				tro.strand = (*sit).second;
				tro.MAPQ = MAPQ;
				strncpy(tro.CIGAR,CIGAR.c_str(),MAX_CIGAR_SZ-1);
				tro.CIGAR[MAX_CIGAR_SZ-1] = '\0';
				tro.readIndex = rIndex;
				
				tro.consensus = consensus;
				tro.qual = QUAL;
				tro.A_SCORE = align_score;
				tro.SIM_MATCHES = positions.size();
				tro.POST_PROB = total_score;
				
				ret_vec.push_back(tro);
				//fprintf(stderr,"[%d/%d] Sequence >%s< aligns to %s:%lu with MAPQ=%d\n",
						//iproc,nproc,tro.READ_NAME,tro.CHR_NAME,tro.CHR_POS,tro.MAPQ);
			}

			return ret_vec;
		}

};

#endif
