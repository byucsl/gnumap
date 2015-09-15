/* BWAGenome.h
 *
 * Copyright 2009-2014 Jamison Dance.
 * Created by Jamison Dance on 12/4/09.
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
 * Subclass of GenomeMem to interface with the bwa code. It allows
 * you to align a sequence against a bwa index and get a list of locations
 * in the genome back.
 */

#ifndef BWAGENOME_H
#define BWAGENOME_H

#include "GenomeMem.h"
#include "UnitTest.h"
#include "Reader.h"
#include "bwtaln.h"

//Just a wrapper to call methods on the bwa code
class BWAGenome : public GemomeMem {
public:
	BWAGenome();
	
	BWAGenome(const char * bwaIndexLocation);
	
	~BWAGenome();
	string GetString(const unsigned long begin, const unsigned int size);
	HashLocation* GetMatches(unsigned int hash);
private:
	bwt_t *bwt[2];
};
#endif
