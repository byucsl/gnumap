/* SNPScoredSeq.h
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

#ifndef SNP_SCORED_SEQ_H
#define SNP_SCORED_SEQ_H

#include "ScoredSeq.h"

class SNPScoredSeq : public ScoredSeq {
	
	public:
		SNPScoredSeq() :ScoredSeq() {}
	
		SNPScoredSeq(string seq, double as, unsigned long pos, bool us) 
			: ScoredSeq(seq,as,pos,us) {}
		
		SNPScoredSeq(const ScoredSeq &other) : ScoredSeq(other) {}
		
		SNPScoredSeq& operator =(const ScoredSeq &other) {
			Init(other);
			return *this;
		}
		
		virtual ~SNPScoredSeq() {}
	
		/*
		 * score will add the log-odds score (adjusted by denom) to each
		 * position on the Genome gen.  Additionally, this can be passed
		 * a flag to specify the number of positions after the start position
		 * to add to (for adding to different lengths of sequences).
		 */
		virtual void score(double, Genome &, unsigned int, Read &, pthread_mutex_t&);
		
	private:
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


};

#endif
