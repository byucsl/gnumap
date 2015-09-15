/* BSScoredSeq.h
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

#ifndef BS_SCORED_SEQ_H
#define BS_SCORED_SEQ_H

#include "ScoredSeq.h"

class BSScoredSeq : public ScoredSeq {

	public:
		BSScoredSeq() :ScoredSeq() {}
	
		BSScoredSeq(string seq, double as, unsigned long pos, bool us) 
			: ScoredSeq(seq,as,pos,us) {}
		
		BSScoredSeq(const BSScoredSeq &other) : ScoredSeq(other) {}
		
		virtual ~BSScoredSeq() {}
		
		/*
		 * score will add the log-odds score (adjusted by denom) to each
		 * position on the Genome gen.  Additionally, this can be passed
		 * a flag to specify the number of positions after the start position
		 * to add to (for adding to different lengths of sequences).
		 */
		virtual void score(double, Genome &, unsigned int, Read &, pthread_mutex_t&);
		
	private:
		inline string GetConsensus(const Read &read) const ;
		
		/*
		 * Because this is used to produce the consensus sequence, we need ambiguity characters...
		 */
		inline char max_char(const float* chr) const;
};

#endif

