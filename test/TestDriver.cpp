/* TestDriver.cpp
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

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "bin_seq.h"
#include "Reader.h"
#include "GenomeMem.h"
#include "SeqReader.h"
#include "const_include.h"
#include "const_define.h"
#include "centers.h"

using namespace std;

const char* genome_file;
const char* output_file;
const char* pos_matrix=NULL;
unsigned int output_type=0;
		/* output_type is: 
		 * 0 = all
		 * 1 = fasta
		 * 2 = fq
		 * 3 = fa&fq
		 * 4 = prb&int
		 */
bool t_fasta=true;
bool t_fq=true;
bool t_int=true;
bool t_prb=true;
unsigned int num_gen=1000;
int mean=0;
int sd = 0;
bool RAND=0;
bool EQ_SPACE=0;

#include "a_matrices.c"

int main() {
	// Use a default MER_SIZE here
	gMER_SIZE = 10;
	gGEN_SKIP = 1;
	
	// Need to do setup work for the rest of the files.
	memset(g_gen_CONVERSION,4,256);
	g_gen_CONVERSION[(int)'a'] = g_gen_CONVERSION[(int)'A'] = 0;
	g_gen_CONVERSION[(int)'c'] = g_gen_CONVERSION[(int)'C'] = 1;
	g_gen_CONVERSION[(int)'g'] = g_gen_CONVERSION[(int)'G'] = 2;
	g_gen_CONVERSION[(int)'t'] = g_gen_CONVERSION[(int)'T'] = 3;
	g_gen_CONVERSION[(int)'n'] = g_gen_CONVERSION[(int)'N'] = 4;
	g_gen_CONVERSION[(int)'\n'] = g_gen_CONVERSION[(int)'\r'] = 5;
	g_gen_CONVERSION[(int)'\0'] = 6;
	g_gen_CONVERSION[(int)'>'] = 7;

	memset(g_bs_CONVERSION,4,256);
	g_bs_CONVERSION[(int)'a'] = g_bs_CONVERSION[(int)'A'] = 0;
	g_bs_CONVERSION[(int)'c'] = g_bs_CONVERSION[(int)'C'] = 1;
	g_bs_CONVERSION[(int)'g'] = g_bs_CONVERSION[(int)'G'] = 2;
	g_bs_CONVERSION[(int)'t'] = g_bs_CONVERSION[(int)'T'] = 3;
	g_bs_CONVERSION[(int)'n'] = g_bs_CONVERSION[(int)'N'] = 4;
	g_bs_CONVERSION[(int)'\n'] = g_bs_CONVERSION[(int)'\r'] = 5;
	g_bs_CONVERSION[(int)'\0'] = 6;
	g_bs_CONVERSION[(int)'>'] = 7;

	memset(gINT2BASE,'?',16);
	gINT2BASE[0] = 'a';
	gINT2BASE[1] = 'c';
	gINT2BASE[2] = 'g';
	gINT2BASE[3] = 't';
	gINT2BASE[4] = 'n';

	// Move this code out of the main body
	setup_alignment_matrices();

#ifdef DISCRETIZE
	fprintf(stderr,"Setting up centers matrices\n");
	// Set up the matrices used to discretize the read values
	init_centers();
	init_lookup();
#endif

	bool all_passed = true;
	unsigned int warnings = 0;

	float array[] = {0.96, 4.1, 0.24, 0.33, 1.45};
	int array_size = sizeof(array) / sizeof(array[0]);
	float* pos = max_element(array, array+array_size);
	unsigned int diff = pos - &array[0];
	printf("Diff: %u\n",diff);
		
/*	size_t find_num = 1234;
	cout << "mod 17: " << find_num % 17 << endl;
	cout << "mod 16: " << find_num % 16 << endl;
	
	int test = 0x01001000;
	cout << "test is: >" << test << "<\n";
	test = test & 0x11111111;
	cout << "test is: >" << test << "<\n";
	
	vector<char> test_v;
	test_v.push_back('1');
	test_v.push_back('2');
	test_v.push_back('3');
	test_v.push_back('4');
	test_v.push_back('5');

	char test_a[5] = {'1','2','3','4','5'};
	char test_axa[5][5];

	cout << "int: " << sizeof(unsigned int) << endl;
	cout << "char: " << sizeof(char) << endl;
	cout << "unsigned char: " << sizeof(unsigned char) << endl;
	cout << "short: " << sizeof(unsigned short) << endl;
	cout << "size_t: " << sizeof(size_t) << endl;
	cout << "wchar_t: " << sizeof(wchar_t) << endl;
	cout << "vector: " << sizeof(test_v) << endl;
	cout << "array: " << sizeof(test_a) << endl;
	cout << "2-dim array: " << sizeof(test_axa) << endl;
	cout << "bool: " << sizeof(bool) << endl;
	cout << "unsigned: " << sizeof(unsigned) << endl;
	cout << "unsigned long: " << sizeof(unsigned long) << endl;
	cout << "unsigned long long: " << sizeof(unsigned long long) << endl;
	cout << "signed long long: " << sizeof(signed long long) << endl;
	cout << "double: " << sizeof(double) << endl;

	signed long long sh = 0;
	for(int i=0; i<45; i++) {
		cout << i << ": " << sh << endl;
		sh = sh << 1;
		sh += 1;
	}
	size_t st = ~0;
	cout << "largest size_t: " << st << endl;
	BIT64_t b64 = ~0;
	cout << "largest 64-bit: " << b64 << endl;
	unsigned int uint = ~0;
	cout << "largest unsigned int: " << uint << endl;
	unsigned long ulong = ~0;
	cout << "largest unsigned long: " << ulong << endl;
	signed long long sll = ~0;
	cout << "largest signed long long: " << sll << endl;
	unsigned short sshort = ~0;
	cout << "largest short: " << sshort+100000 << endl;
	unsigned int nint = 1;
	nint *= -1;
	cout << "negative uing: " << nint << endl;

	char c = (char)-128;
	cout << "negative character: " << c << endl; 
	cout << "back to int: " << (int)c << endl;
	test_v.push_back('6');
	test_v.push_back('7');
	test_v.push_back('8');
	test_v.push_back('9');
	test_v.push_back('0');
	test_v.push_back('a');
	test_v.push_back('b');
	test_v.push_back('c');
	
	char test_b[13] = {'1','2','3','4','5','6','7','8','9','0','a','b','c'};
	cout << "vector: " << sizeof(test_v) << endl;
	cout << "array: " << sizeof(test_b) << endl;
	cout << test_v.size() << endl;
	
	string str = "actgccc";
	bin_seq bs;
	cout << "bin_seq: " << sizeof(bs) << endl;
	
	
	vector<char>* s_vector = new vector<char>[50];
	
	s_vector[5].push_back('c');
	cout << "in the vector array: " << s_vector[5][0] << endl;



	double test1 = 1.7995674832123e-09;
	float test2 = test1;
	cout << "Test1: " << test1 << " test2: " << test2 << endl;
	*/


	cout << "\n******************Testing******************" << endl;
///*
	cout << "Testing bin_seq:" << endl;
	if(bin_seq::Test(cerr, warnings))
		cout << "Success" << endl << endl;
	else {
		cout << "Failure" << endl << endl;
		all_passed = false;
	}
//*/	
/*
	cout << "Testing Reader:" << endl;
	if(Reader::Test(cout,warnings))
		cout << "Success" << endl << endl;
	else {
		cout << "Failure" << endl << endl;
		all_passed = false;
	}
//*/
///*
	cout << "Testing Genome:" << endl;
	if(GenomeMem::Test(cerr, warnings))
		cout << "Success" << endl << endl;
	else {
		cout << "Failure" << endl << endl;
		all_passed = false;
	}
//*/
///*
	cout << "Testing SeqReader" << endl;
	try {
	if(SeqReader::Test(cerr, warnings))
		cout << "Success" << endl << endl;
	else {
		cout << "Failure" << endl << endl;
		all_passed = false;
	}
	}
	catch(Exception* e) {
		fprintf(stderr,"ERROR: %s\n",e->GetMessage());
	}
//*/	
	
	if(warnings)
		cout << "Number of Warnings: " << warnings << endl;
	if(all_passed)
		cout << "*******************************************\n    All tests have passed successfully." << endl;
	else
		cout << "*******************************************\n              Failed Tests." << endl;
		
	cout << endl;
	
	return 0;
}
