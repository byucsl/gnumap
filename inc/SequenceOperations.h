/* SequenceOperations.h
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

#ifndef SEQUENCE_OPERATIONS_H
#define SEQUENCE_OPERATIONS_H

#include <string>
#include <algorithm>
#include <cstring>
#include "const_include.h"

using namespace std;

inline void fix_CIGAR_for_deletions(string &cigar) {
	int i;
    if(cigar[cigar.size()-1] == 'D') {
        // find the first character before the number of deletions in the CIGAR string
        for(i=cigar.size()-2; i>=0; i--) {
            if(!isdigit(cigar[i]))
                break;
        }
        cigar = cigar.substr(0,i+1);
    }
}

inline string trim_end_string(string &s, unsigned int start2trim0) {
  string temp = "";
  //int len = s.length();
  
  for(unsigned int i=0; i<start2trim0; i++)
  {
    temp+=s[i];
  }
  
  return temp;
}

inline string reverse_comp(string &s) {
	string temp = "";
	int len = s.length();
	
	for(int i=1; i<=len; i++) {
		switch(s[len-i]) {
			case 'a':
				temp += 't';
				break;
			case 'A':
				temp += 'T';
				break;
			case 't':
				temp += 'a';
				break;
			case 'T':
				temp += 'A';
				break;
			case 'c':
				temp += 'g';
				break;
			case 'C':
				temp += 'G';
				break;
			case 'g':
				temp += 'c';
				break;
			case 'G':
				temp += 'C';
				break;
			case '-':
				temp += '-';
				break;
			default:
				temp += 'n';
				break;
		}
	}
	
	return temp;
}

inline string reverse_qual(string &s) {
	string temp = "";
	int len = s.length();
	
	for(int i=1; i<=len; i++) {
		temp += s[len-i];
	}
	
	return temp;
}

inline string reverse_CIGAR(char* s) {
	string temp_CIGAR;
	string number = "";
	for(unsigned int i=0; i<strlen(s); i++) {
		if(48 <= (int)s[i] && (int)s[i] <= 58) {
			number += s[i];
		}
		else {	// found a character.  Reset number and prepend to CIGAR
			temp_CIGAR = number+s[i] + temp_CIGAR;
			number = "";
		}
	}

	return temp_CIGAR;
}

inline void reverse_comp(float** pwm, unsigned int length) {
	float temp_pwm[length][4];


	//for(r_vit = v.rbegin(); r_vit != v.rend(); ++r_vit) {
	for(unsigned int i = 0; i < length; ++i) {
		//need to switch places with the characters: a-t, g-c
		//also need to rev-com it.
		temp_pwm[length-1-i][A_POS] = pwm[i][T_POS];
		temp_pwm[length-1-i][C_POS] = pwm[i][G_POS];
		temp_pwm[length-1-i][G_POS] = pwm[i][C_POS];
		temp_pwm[length-1-i][T_POS] = pwm[i][A_POS];
	}
		
	for(unsigned int i=0; i<length; i++) {
		for(unsigned int j=0; j<4; j++) {
			pwm[i][j] = temp_pwm[i][j];
		}
	}
	
	//cout << "temp's size: " << temp.size() << endl;
	//cout << "originals size: " << v.size() << endl; 
}

inline float** reverse_comp_cpy(float** pwm, unsigned int length) {
	float** temp_pwm = new float*[length];

	for(unsigned int i=0; i<length; i++) {
		temp_pwm[length-1-i] = new float[4];
		temp_pwm[length-1-i][A_POS] = pwm[i][T_POS];
		temp_pwm[length-1-i][C_POS] = pwm[i][G_POS];
		temp_pwm[length-1-i][G_POS] = pwm[i][C_POS];
		temp_pwm[length-1-i][T_POS] = pwm[i][A_POS];
	}

	return temp_pwm;
}


inline float** reverse_comp_cpy_phmm(float** pwm, unsigned int length) {
	float** temp_pwm = new float*[length];

	for(unsigned int i=0; i<length; i++) {
		temp_pwm[length-1-i] = new float[NUM_SNP_VALS];
		temp_pwm[length-1-i][A_POS] = pwm[i][T_POS];
		temp_pwm[length-1-i][C_POS] = pwm[i][G_POS];
		temp_pwm[length-1-i][G_POS] = pwm[i][C_POS];
		temp_pwm[length-1-i][T_POS] = pwm[i][A_POS];
		temp_pwm[length-1-i][N_POS] = pwm[i][N_POS];
#ifdef _INDEL
		temp_pwm[length-1-i][INS_POS] = pwm[i][INS_POS];
		temp_pwm[length-1-i][DEL_POS] = pwm[i][DEL_POS];
#endif
	}

	return temp_pwm;
}

inline unsigned int max_flt(float array[], unsigned int array_size) {
	float* pos = max_element(array, array+array_size);
	unsigned int max_pos = pos - &array[0];
	return max_pos;
}

// Need the MAX_PRB for fasta files where there's a 100% call at each base
// Otherwise, it won't print out any characters
#define MAX_PRB	0.9999

inline string str2qual(Read &r) {
	if(r.fq.size() > 0)
		return r.fq;

	string qual = "";

	unsigned int i;
	for(i=0; i<r.length; i++) {
		int max_pos = max_flt(r.pwm[i],4);
		if(r.pwm[i][max_pos] > MAX_PRB)
			//if(gILLUMINA) 	qual += (char)((-10 * log((1-MAX_PRB)/MAX_PRB)/log(10.))+33);
			//else 			
			qual += (char)((-10 * log(1-MAX_PRB)/log(10.))+33);
			
		else {
			//fprintf(stderr,"max_pos of %d and resulting character of %c\n",
			//		max_pos, (char)((-10 * log(1-r.pwm[i][max_pos])/log(10))+33));
			//if(gILLUMINA) 	qual += (char)((-10 * log((1-r.pwm[i][max_pos])/r.pwm[i][max_pos])/log(10.))+33);
			//else			
			qual += (char)((-10 * log(1-r.pwm[i][max_pos])/log(10))+33);
		}
	}

	return qual;
}

inline void delete_read(Read* temp_read) {
	//done with the read.  Delete it.
	for(unsigned int j=0; j<temp_read->length; j++) {
		delete[] temp_read->pwm[j];
	}
	delete[] temp_read->pwm;
	if(temp_read->name)
		delete[] temp_read->name;
	
	delete temp_read;	
}

#endif
