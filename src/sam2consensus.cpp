/* sam2consensus.cpp
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

#include <pthread.h>
#include <iostream>
#include <string>
#include <fstream>

#include "Exception.h"
#include "const_include.h"
#include "const_define.h"
#include "ScoredSeq.h"
#include "NormalScoredSeq.h"
#include "SNPScoredSeq.h"
#include "BSScoredSeq.h"
#include "SeqReader.h"
#ifdef GENOME_STL
#include "GenomeSTL.h"
typedef GenomeSTL GENOME_t;
#else
#include "GenomeMem.h"
typedef GenomeMem GENOME_t;
#endif


using namespace std;

#define MAX_NAME_SZ		1024
#define MAX_CIGAR_SZ	1024
#define gREAD_BUFFER	1024


const char* sam_file;
const char* genome_file;
const char* output_file = "sam2consensus.out";
const char* pos_matrix = NULL;
ofstream of;
Genome *gGen;

bool has_printed_error = false;
int gSTART_POS = 1; // The starting offset for the genome.
					// Many programs use 1 instead

/*****************************************************
 * These next functions are used in parsing command-
 * line arguments
 *****************************************************/
int ParseCmdLine(const int argc, const char* argv[]);
void GetParseError(ostream &os, const char* argv[]);
int set_arg_ext(const char* param, const char* assign, int &count);
int set_arg(const char* param, const char* assign, int &count);
CMDLine this_cmd;
/*****************************************************/


pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;



void usage(int has_error, char* errmessage) {
	if(has_error) cout << endl << errmessage << endl ;

	cout << "Usage: sam2consensus [options] <file_to_parse>\n"
		 << "  -g, --genome=STRING          Genome .fa file(s)\n"
		 << "  -o, --output=STRING          Output file\n"
		 << "  -v, --verbose=INT            Verbose (default=0)\n"
		 << "\n"
		 << "Help options:\n"
		 << "  -?, --help                   Show this help message\n";
	exit(1);
}

vector<string> splitFileNames(const char* fn) {
	vector<string> names;
	
	//Separate the different files.
	char files[strlen(fn)+1];
	strcpy(files,fn);
	
	const char* delims = " \'\"\n,";
	char* token = strtok(files,delims);
	if(gVERBOSE > 1)
		cout << "Found:" << endl;

	while(token != NULL) {
		if(gVERBOSE > 1)
			cout << "\t" << token << endl;
		string temp_str = token;
		names.push_back(temp_str);
		
		token = strtok(NULL,delims);
	}
	
	return names;
}


#define MAX_QUAL_CHAR	'!'

/**
 * Will fill the vector records with gREAD_BUFFER TopReadOutput objects
 *
 * @return false if there is nothing more to read (file is finished)
 */
bool fillVector(ifstream &in, TopReadOutput records[gREAD_BUFFER]) {
	//fprintf(stderr,"Filling vector again...\n");
	unsigned int readsRead = 0;
	unsigned int linesRead = 0;
	
	while(!in.eof() && readsRead < gREAD_BUFFER) {
		string line;
		getline(in, line);
		linesRead++;
		if(line.size() < 1) {
			fprintf(stderr,"Short line at line %u (%s).  Continuing\n",linesRead,line.c_str());
			continue;
		}
		if(line[0] == '#') {	// means it a comment
			fprintf(stderr,"Ignoring comment %s\n",line.substr(1).c_str());
			continue;
		}
		//extended comment as well
		if(line[0] == '<' && line[1] == '<') {
			string comment_tag = line.substr(2);
			fprintf(stderr,"Found comment line with tag %s\n",comment_tag.c_str());
			while(line.compare(comment_tag) != 0) {
				getline(in,line);
				linesRead++;
			}
			fprintf(stderr,"Finished comment line\n");
			continue;
		}
		records[readsRead].qual = "";
		
		//fprintf(stderr,"Found line %s\n",line.c_str());
			
		char c_line[line.size()];
		strcpy(c_line,line.c_str());
		char delims[] = " \t\n";
		char *result = 0;
		
		/*
		 * Have a SAM line of the format:
		 * <QNAME> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL>
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
		 * - MRNM will be '=' for all sequences (it's the mate-pair name)
		 * - MPOS is ignored.  Always '0'
		 * - ISIZE is ignored.  Always '0'
		 * - SEQ is the sequence
		 * - QUAL is the Phred-based quality of the sequences
		 */
		unsigned int temp_ui;
		
		result = strtok(c_line,delims);
		int i;
		// Continue until the end of the tokenizer
		for(i=0; result; i++) {
			switch(i) {
				case 0: // QNAME
					strncpy(records[readsRead].READ_NAME,result,MAX_NAME_SZ);
				break;
				case 1: // FLAG
					temp_ui = (unsigned int)atoi(result);
					if(temp_ui & 0x0010)
						records[readsRead].strand = true;
					else
						records[readsRead].strand = false;
				break;
				case 2:	// RNAME
					strncpy(records[readsRead].CHR_NAME,result,MAX_NAME_SZ);
				break;
				case 3: // POS
					// We need to subtract off the start position for 1-based programs
					records[readsRead].CHR_POS = (unsigned long)atol(result)-gSTART_POS;
				break;
				case 4: // MAPQ
					records[readsRead].MAPQ = atoi(result);
				break;
				case 5: // CIGAR
					strncpy(records[readsRead].CIGAR,result,MAX_CIGAR_SZ);
				break;
				case 6: // MRNM
					// do nothing
				break;
				case 7:	// MPOS
					// do nothing
				break;
				case 8: // ISIZE
					// do nothing
				break;
				case 9: // SEQ
					records[readsRead].consensus = result;
				break;
				case 10: // QUAL
					records[readsRead].qual = result;
				break;
				
				default: // Just ignore the rest of the line
					break;
					/*
					fprintf(stderr,"Error at position %d with token %s in line: %s\n",i,result,line.c_str());
					throw new Exception("SAM file not formatted correctly");
					*/
			}
			
			// Get the next one
			result = strtok(NULL,delims);
			
		} // end for
		// Make sure you've filled it properly
		if(i < 11) {
			if(records[readsRead].qual.size() == 0 && i == 10) {
				if(!has_printed_error) {
					fprintf(stderr,"ERROR:  Missing qualities?\n");
					fprintf(stderr,"\tCan still recover, but there might be something wrong with your file.\n");
					fprintf(stderr,"\tPress Ctrl+c to abort\n");
					has_printed_error = true;
				}
				
				// Make up the qualities as if it's a fasta file
				records[readsRead].qual = "";
				for(unsigned int i=0; i<records[readsRead].consensus.size(); i++) {
					records[readsRead].qual += MAX_QUAL_CHAR;
				}
			}
			else {
				fprintf(stderr,"Error (didn't have 11 elements, instead %d) in line: %s\n",i,line.c_str());
				throw new Exception("SAM file not formatted correctly");
			}
		}
		readsRead++;
	}
	
	//fprintf(stderr,"Stored %u reads to vector, and eof is %d\n",readsRead,in.eof());
	if(readsRead < gREAD_BUFFER)
		records[readsRead].MAPQ = -1;
	
	if(in.eof())
		return false;
	return true;

}

unsigned int GetConsensusFromSAM(const char* fn) {
	fprintf(stderr,"Reading from SAM File %s\n",fn);

	unsigned int total_records = 0;

	ifstream in;
	in.open(fn);
	
	if(in.fail())
		throw new Exception("The previous operation failed.  Check your SAM file\n");
	if(in.bad())
		throw new Exception("Stream state is undefined; the stream can no longer be used.  Check your SAM file\n");
	
	//vector<TopReadOutput> records;
	TopReadOutput records[gREAD_BUFFER];
	bool finished = false;
		
	while(!finished) {
		finished = !fillVector(in, records);
		
		for(unsigned int i=0; i<gREAD_BUFFER; i++) {
			// Means there's not another record
			if(records[i].MAPQ < 0)
				break;
				
			unsigned int len = records[i].consensus.size();
			string consensus = gGen->GetStringRelative(records[i].CHR_NAME, records[i].CHR_POS, len);
			of << records[i].READ_NAME << "\t" << consensus << "\n";
			
			total_records++;

			// Reset this position
			records[i].MAPQ = -1;

		} // end for
		if(gVERBOSE)
			fprintf(stderr,"Found %u records thus far\n",total_records);

	} //end while(!finished)
	in.close();

	
	return total_records;
}

/**
 * Reads all the SAM files that are contained in 'fn', parsing through the string
 *
 * @return The total number of SAM records
 */
unsigned int BatchConsensus(const char* fn) {
	unsigned int total_records = 0;
	
	vector<string> files = splitFileNames(fn);
	
	if(gVERBOSE)
		fprintf(stderr,"Reading %lu SAM files\n",files.size());
	
	for(unsigned int i=0; i<files.size(); i++) {
		total_records += GetConsensusFromSAM(files[i].c_str());
		if(gVERBOSE)
			fprintf(stderr,"Found %u SAM records so far...\n",total_records);
	}
	
	return total_records;
}

/* readPWM is used for getting a PWM matrix representing the match (m) and mismatch (mm)
 * penalties.
 * The input is a file, which should be a 4x4 matrix of the format:
 * 				 A   C   G   T
 * 			A    m  mm  mm  mm
 *			C   mm   m  mm  mm
 * 			G   mm  mm   m  mm
 * 			T   mm  mm  mm   m
 * For ease in reading, the letters defining the rows and columns are disregarded.
 *
 *
 * Note:  readPWM will modify the global gALIGN_SCORES variable.
 */
void readPWM(const char* fn) {
	float temp_ALIGN_SCORES[5][4];

	int THIS_BUFFER = 100;
	
	ifstream in;
	in.open(fn);
	char temp_chr[THIS_BUFFER];
	bool contains_labels=false;
	float a,c,g,t;
	char label[THIS_BUFFER];
	char alabel[THIS_BUFFER];
	char clabel[THIS_BUFFER];
	int count=0;
	int num_lines=5;

	// Get all the lines of input
	while(in.getline(temp_chr,THIS_BUFFER) && count < num_lines) {
		if(sscanf(temp_chr,"%s %s %s %s",alabel,clabel,label,label) == 4 
				&& (tolower(*alabel) == 'a')
				&& (tolower(*clabel) == 'c') ) {
			if(gVERBOSE > 1)
				printf("Matched here: %s\n",temp_chr);
			in.getline(temp_chr,THIS_BUFFER);
			contains_labels=true;
		}
		if(gVERBOSE > 1)
			cout << "Line: " << temp_chr << endl;

		if(contains_labels) {
			if(sscanf(temp_chr,"%s %f %f %f %f",label,&a,&c,&g,&t) != 5) {
				char* error = new char[THIS_BUFFER];
				strcat(error,"Error in Score File: ");
				strcat(error,temp_chr);
				throw(error);
			}
		}
		else {
			if(sscanf(temp_chr,"%f %f %f %f",&a,&c,&g,&t) != 4) {
				char* error = new char[THIS_BUFFER];
				strcat(error,"Error in Score File: ");
				strcat(error,temp_chr);
				throw(error);
			}
		}
		
		temp_ALIGN_SCORES[count][0] = (float)a;
		temp_ALIGN_SCORES[count][1] = (float)c;
		temp_ALIGN_SCORES[count][2] = (float)g;
		temp_ALIGN_SCORES[count][3] = (float)t;

		count++;
	}
	if(count < num_lines) {
		throw new Exception("Error in Score File:  Not enough lines");
	}

	cout << "Count: " << count << "\tNum lines: " << num_lines << endl;

#ifdef DEBUG
	//for testing purposes...
	if(gVERBOSE > 1) {
		cout << "Matrix: " << endl;
		cout << "\t0\t1\t2\t3" << endl;
		for(int i=0; i<5; i++) {
			cout << i << "\t";
			for(int j=0; j<4; j++) 
				cout << temp_ALIGN_SCORES[i][j] << "\t";
			cout << endl;
		}
		cout << endl;
		throw new Exception("Passed.");
	}
#endif

	for(unsigned int i=0; i<4; i++) {
		gALIGN_SCORES[(unsigned int)'a'][i] = temp_ALIGN_SCORES[0][i];
		gALIGN_SCORES[(unsigned int)'c'][i] = temp_ALIGN_SCORES[1][i];
		gALIGN_SCORES[(unsigned int)'g'][i] = temp_ALIGN_SCORES[2][i];
		gALIGN_SCORES[(unsigned int)'t'][i] = temp_ALIGN_SCORES[3][i];
		gALIGN_SCORES[(unsigned int)'n'][i] = temp_ALIGN_SCORES[4][i];
	}
}

/**
 * getPWM will return the current substitution matrix in a string that can be
 * printed to the console
 */
string getPWM() {
	ostringstream pwm;
	pwm << "\t\tA\tC\tG\tT\n";
	pwm << "\tA";
	for(unsigned int i=0; i<4; i++) {
		pwm << "\t" << gALIGN_SCORES[(unsigned int)'a'][i];
	}
	pwm << "\n\tC";
	for(unsigned int i=0; i<4; i++) {
		pwm << "\t" << gALIGN_SCORES[(unsigned int)'c'][i];
	}
	pwm << "\n\tG";
	for(unsigned int i=0; i<4; i++) {
		pwm << "\t" << gALIGN_SCORES[(unsigned int)'g'][i];
	}
	pwm << "\n\tT";
	for(unsigned int i=0; i<4; i++) {
		pwm << "\t" << gALIGN_SCORES[(unsigned int)'t'][i];
	}
	pwm << "\n\tN";
	for(unsigned int i=0; i<4; i++) {
		pwm << "\t" << gALIGN_SCORES[(unsigned int)'n'][i];
	}
	pwm << endl;
	
	return pwm.str();
}

#include "a_matrices.c"


int main(const int argc, const char* argv[]) {
	cerr << "This is get_batch_consensus, Version 1.0, for public and private use." << endl;

	double time_prog_start = When();
	
	// Move this into const_define.h
	InitProg();

	//set up alignment matrix
	// Move this code out of the main body
	setup_alignment_matrices();


	cout << endl << "Command Line Arguments:  ";
	for(int i=0; i<argc; i++)
		cout << argv[i] << ' ';
	cout << endl;
	
	gVERBOSE = 0;

	if(argc == 1) 
		usage(0,"");

	if(argc < 3) {
		usage(1, "Need at least a genome and a SAM file to find consensus.");
		exit(1);
	}


	// Parse the command line
	int rc = ParseCmdLine(argc,argv);
	// an error occurred during option processing
	if(rc != 0) {
		GetParseError(cerr, argv);
		return -1;
	}

	ostringstream params;
	params << "Parameters: " << endl;
	params << "\tVerbose: " << gVERBOSE << endl;
	params << "\tGenome file(s): " << genome_file << endl;
	params << "\tOutput file: " << output_file << endl;
	params << "\tSAM file: " << sam_file << endl;

	if(gVERBOSE) {
		cerr << params.str() << endl;
	}

	/************************************************/
	/* Finished initialization, now do the stuff!	*/
	/************************************************/
	

	of.open(output_file);
	if(of.bad() || of.fail()) {
		cerr << "ERROR:\n\t";
		perror("Error in reading output file");
		return -1;
	}
	
	unsigned int total_records = 0;
	try {
		gGen = new GENOME_t;
		gGen->use(genome_file);
		gGen->StoreGenome(false);
		
		total_records = BatchConsensus(sam_file);
		fprintf(stderr,"Total Read: %u\n",total_records);
		
		fprintf(stderr,"\nFinished! Output written to %s\n",output_file);

		delete gGen;
	}
	catch(Exception *e) {
		fprintf(stderr,"ERROR: %s\n",e->GetMessage());
	}
	
	cerr << "\nFinished.\n"
		 << "\tTotal time: " << When()-time_prog_start << " seconds\n"
		 << "\tFound " << total_records << " records\n";
}


/************************************************************
 * COMMAND-LINE STUFF                                       *
 ************************************************************/

int ParseCmdLine(const int argc, const char* argv[]) {
	
	bool found_other = false;

	for(int i=1; i<argc; i++) {
		if(*argv[i] == '-') { //check for normal
			int set_ret;

			if(*((argv[i])+1) == '-')  { //check for extended
				//printf("Found extended\n");
				set_ret = set_arg_ext(argv[i],argv[i+1],i);
			}
			else {
				set_ret = set_arg(argv[i],argv[i+1],i);
				i++;	//increment past the next expression
			}
			
			if(set_ret != 0) {
				this_cmd.errnum = set_ret;
				this_cmd.errpos = i;
				return -1;
			}

		}
		else {
			//if we've already found the the other cmd arg, it doesn't work. 
			if(found_other) {
				this_cmd.errnum = INVALID_ARG;
				this_cmd.errpos = i;
				return INVALID_ARG;
			}
			else {	//set the other arg to this pointer
				sam_file = argv[i];
				found_other = true;
			}
		}
	}

	return 0;

}

void GetParseError(ostream &os, const char* argv[]) {
	if(this_cmd.errnum == PARSE_ERROR) {
		os << "Irregular Parameter in: " << argv[this_cmd.errpos] << endl;
	}
	else if(this_cmd.errnum == NO_MATCHING_ARG) {
		os << "No matching arg in: " << argv[this_cmd.errpos] << endl;
	}
	else if(this_cmd.errnum == INVALID_ARG) {
		os << "Please specify only a single SAM file (found additional at >"
		   << argv[this_cmd.errpos] << "<)" << endl;
	}
	else {
		os << "Unknown error " << this_cmd.errnum << " in: " << argv[this_cmd.errpos] << endl;
	}

}


int set_arg(const char* param, const char* assign, int &count) {
	switch(*(param+1)) {	//switch on the first character
		case 'g':
			genome_file = assign;
			break;
		case 'o':
			output_file = assign;
			break;
		case 'v':
			if(sscanf(assign,"%d",&gVERBOSE) < 1)
				return PARSE_ERROR;
			break;
		case '?':
			usage(0,"");
		default:
			return PARSE_ERROR;
	}

	return 0;

}

enum {
	PAR_ORIG=256,
	PAR_GENOME,
	PAR_OUTPUT,
	PAR_VERBOSE,
};

int set_arg_ext(const char* param, const char* assign, int &count) {
	int which = 0;
	int adjust=0;

	if(strncmp(param+2,"genome=",7) == 0) {	//check to see if it's the genome
		//which = 'g';
		which = PAR_GENOME;
		adjust = 8;	//the size of "genome" plus 2
	}
	else if(strncmp(param+2,"output=",7) == 0) {
		//which = 'o';
		which = PAR_OUTPUT;
		adjust = 8;
	}
	else if(strncmp(param+2,"verbose=",8) == 0) {
		//which = 'v';
		which = PAR_VERBOSE;
		adjust = 9;
	}
	else 
		return NO_MATCHING_ARG;

	param += ++adjust;

	switch(which) {	//switch on the first character
		//case 'g':
		case PAR_GENOME:
			genome_file = param;
			break;
		//case 'o':
		case PAR_OUTPUT:
			output_file = param;
			break;
		//case 'v':
		case PAR_VERBOSE:
			if(sscanf(param,"%d",&gVERBOSE) < 1)
				return PARSE_ERROR;
			break;
		default:
			return NO_MATCHING_ARG;
	}

	return 0;
}
