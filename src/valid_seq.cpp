/* valid_seq.cpp
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

#include "valid_seq.h"

// map<string,set<range> > valid_seqs;
// int seq_size;
// int range_size;

valid_seq::valid_seq() {
	total_seqs = 0;
}

valid_seq::valid_seq(char* fn) :
	total_seqs(0) {
	store(fn);
	cout << "total seqs: " << total_seqs << endl;
}

void valid_seq::store(const char* fn) {
	total_seqs = 0;
	ifstream infile(fn);
	string line="";
	const char* DELIMS = "\t\n ";
	
	if(!infile || infile.bad())
		throw "Bad file";
	
	while(line == "") 		// Ignore any leading blank lines
		getline(infile,line);	// The first line contains important info
	
	//cout << "Line: " << line << endl;
	
	char* temp = new char[100];
	strcpy(temp,line.c_str());
	
	char* result = strtok(temp,DELIMS);
	seq_size = atoi(result);
	
	result = strtok(NULL,DELIMS);
	assert(result != NULL);
	range_size = atoi(result);
	assert(range_size >= seq_size);
	
	int shift = seq_size/2 - range_size/2;
	
	while(getline(infile,line)) {
		
		strcpy(temp,line.c_str());
		result = strtok(temp,DELIMS);
		
		if(result == NULL)
			break;
		
		string chr_name = result;

		result = strtok(NULL,DELIMS);
		unsigned int pos = atol(result);
		
		// We want to find the middle of the sequence, then set the range as the surrounding pieces.
		unsigned int start = pos + shift;
		range chr_range(start, start+range_size);
		
		if( (valid_seqs[chr_name].insert(chr_range)).second)
			total_seqs++;
	}

	// cout << "total seqs: " << total_seqs << endl;

	delete[] temp;
}

bool valid_seq::is_valid(string &chr_pos, int size) {
	//cout << "chr_pos: " << chr_pos << endl;
	
	const char* DELIMS = " :";
	
	char* temp = new char[100];
	strcpy(temp,chr_pos.c_str());
	
	char* result = strtok(temp,DELIMS);
	string chr = result;
	
	result = strtok(NULL,DELIMS);
	assert(result != NULL);
	unsigned int pos = atol(result);
	
	range seq_range(pos,pos+size);
	
	if(valid_seqs[chr].find(seq_range) != valid_seqs[chr].end())
		return true;
	
	delete[] temp;
	return false;
	
}

unsigned int valid_seq::get_next_rand(unsigned int &sequence) {
	range to_return;
		
	// The num given in this function is only a counter from DataGen.  Because of this,
	// we need to mod it with the total number of sequences.
	unsigned int seq_num = sequence % total_seqs;
	unsigned int counter = 0;
	map<string,set<range> >::iterator mit = valid_seqs.begin();

	for(unsigned int i=0; i<valid_seqs.size(); i++) {
		if( (counter + (*mit).second.size()) > seq_num ) {
			set<range>::iterator sit = (*mit).second.begin();
			
			//Just increment through until we have the right range.
			for(unsigned int j=0; j < seq_num - counter; j++, ++sit);
			
			to_return = *sit;
			break;
		}
		else
			++mit;
	}
	
	float scale = rand()/float(RAND_MAX);	// in order to be truly RANDOM
	unsigned int spot = (unsigned int)(scale * (range_size - seq_size));
											// this allows for range_size-seq_size to be 0
	return to_return.begin + spot;
}

