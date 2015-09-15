#include <pthread.h>
#include <omp.h>
#include <iostream>
#include <string>
#include <fstream>
#include <signal.h>

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
#else
#include "GenomeMem.h"
#endif
#include "SequenceOperations.h"

#ifdef GENOME_STL
typedef GenomeSTL GENOME_t;
#else
typedef GenomeMem GENOME_t;
#endif
using namespace std;

#define MAX_NAME_SZ		1024
#define MAX_CIGAR_SZ	1024
#define gREAD_BUFFER	1024*512


const char* sam_file;
const char* genome_file;
const char* output_file = "sam2sgr";
const char* pos_matrix = NULL;
int gNUM_THREADS=1;
Genome *gGen;

bool has_printed_error = false;
int gSTART_POS = 1; // The starting offset for the genome.
		// Most programs use 1

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

#include "a_matrices.c"

void usage(int has_error, char* errmessage) {
	if(has_error) cout << endl << errmessage << endl ;

	cout << "Usage: sam2sgr [options] <file_to_parse>\n"
		 << "  -g, --genome=STRING          Genome .fa file(s)\n"
		 << "  -o, --output=STRING          Output file\n"
		 << "  -v, --verbose=INT            Verbose (default=0)\n"
		 << "  -c, --num_proc=INT           Number of processors to use\n"
		 << "  -s, --start_pos=INT          The start position for the mapping \n"
		 << "                               program (GNUMAP's default is currently\n"
		 << "                               set to 1, which is default)\n"
		 << "  -S, --subst_file=STRING      Position-Weight Matrix file for DNA\n"
		 << "                               Substitutions\n"
		 << "  -G, --gap_penalty=DOUBLE     Gap Penalty for Substitution Matrix\n"
		 << "  -M, --max_gap=INT            Maximum Number of Gaps to use in Alignment\n"
		 << "  -0, --print_full             Print locations for the entire sequence, not\n"
		 << "                               just for the beginning.\n"
		 << "  -b, --bs_seq                 Flag to turn on the C to T conversion, used in\n"
		 << "                               bisulfite sequence analysis\n"
		 << "  -d, --a_to_g                 Flag that allows for A to G conversion\n"
		 << "  --bin_size=INT               The resolution for GNUMAP (default: 8)\n"
		 << "  --snp                        Turn on SNP mapping (will output a .sgrex file)\n"
		 << "  --snp_pval=DOUBLE            P-Value cutoff for calling SNPs\n"
		 << "                               (default: 1e-4)\n"
		 << "  --snp_monop                  Flag that turns on monoploid SNP calling\n"
		 << "                               (default: diploid SNP calling)\n"
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
bool fillVector(ifstream &in, TopReadOutput* records) {
	//fprintf(stderr,"Filling vector again...\n");
	unsigned int readsRead = 0;
	unsigned int linesRead = 0;
	
	while(!in.eof() && readsRead < gREAD_BUFFER) {
		string line;
		getline(in, line);
		linesRead++;
		if(line.size() < 1) {
			//fprintf(stderr,"Short line at line %u (%s).  Continuing\n",linesRead,line.c_str());
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
						records[readsRead].strand = POS_STRAND;	// SAM requires to be written in pos
					else
						records[readsRead].strand = POS_STRAND;
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
					/* Don't revCOMP because SAM requires forward strand
					if(records[readsRead].strand == NEG_STRAND) {
						string consensus = result;
						records[readsRead].consensus = reverse_comp(consensus);
					}*/
				break;
				case 10: // QUAL
					records[readsRead].qual = result;
					/*
					if(records[readsRead].strand == NEG_STRAND) {
						string qual = result;
						records[readsRead].qual = reverse_qual(qual);
					}*/
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
				fprintf(stderr,"==ERROR  (didn't have 11 elements, instead %d) in line: %s\n",i,line.c_str());
				//throw new Exception("SAM file not formatted correctly");
				fprintf(stderr,"==  Skipping record.\n");
				records[readsRead].MAPQ = -1;
				continue;
			}
		}

		// If there's a star, it's not a real match
		if(strcmp(records[readsRead].CHR_NAME,"*") == 0) {
			records[readsRead].MAPQ = -1;
		}
		else
			readsRead++;
	}
	
	//fprintf(stderr,"Stored %u reads to vector, and eof is %d\n",readsRead,in.eof());
	for(; readsRead < gREAD_BUFFER; readsRead++)
		records[readsRead].MAPQ = -1;
	
	if(in.eof())
		return false;
	return true;

}

unsigned int AddSamtoGenome(const char* fn) {
	fprintf(stderr,"\nReading from SAM File %s\n",fn);

	unsigned int total_records = 0;
	double file_start_time = When();

	ifstream in;
	in.open(fn);
	
	if(in.fail())
		throw new Exception("The previous operation failed.  Check your SAM file\n");
	if(in.bad())
		throw new Exception("Stream state is undefined; the stream can no longer be used.  Check your SAM file\n");
	
	//vector<TopReadOutput> records;
	TopReadOutput* records = new TopReadOutput[gREAD_BUFFER];
	bool finished = false;
#pragma omp parallel shared(finished,lock,records) reduction(+:total_records) default(shared)
{
	while(!finished) {
#pragma omp  single
{
		finished = !fillVector(in, records);
		if(gVERBOSE && total_records != 0)
			fprintf(stderr,"Read %u SAM, %.2fs total\n",total_records,When()-file_start_time);
}
		
		
		ScoredSeq* ss = 0;
		int i;
#pragma omp for private(ss,i,gALIGN_SCORES,gPHMM_ALIGN_SCORES)
		for(i=0; i<(int)gREAD_BUFFER; i++) {
			// Means there's not another record
			if(records[i].MAPQ < 0) {
				continue;
			}
			//fprintf(stderr," -- %d\n",i);
				
			
			// Needs the sequence, alignment score, position and strand:
			// ScoredSeq(string seq, double as, unsigned long pos, bool us)
			double alignment_score = 1-pow(10.0,(double)(records[i].MAPQ)/-10.0);
			// Expecting a log-aligned score here.  We'll just make do with what we have
			alignment_score = log(alignment_score);

			if(alignment_score > 1) {
				fprintf(stderr,"Alignment score is %f from MAPQ of %d\n",alignment_score,records[i].MAPQ);
			}
			unsigned int len = gPRINT_FULL ? records[i].consensus.size() : 1;

			
			Read* r = NULL;

			// Need the absolute genome position
			unsigned long abs_gen_pos;
			try {
				r = SeqReader::getReadFromStrings(records[i].READ_NAME, 
										records[i].consensus, records[i].qual);
				abs_gen_pos = gGen->GetAbsolutePosition(records[i].CHR_NAME,
																		records[i].CHR_POS);
			}
			catch(Exception *e) {
				fprintf(stderr,"ERRROR with sequence %s: %s\n",records[i].READ_NAME,e->GetMessage());
				fprintf(stderr," -- Continuing --\n");
				//exit(0);
				continue;
			}
			
			// Take the revc out here so we don't need to do it inside the ScoredSeq constructor
			//string ssStr = records[i].strand == POS_STRAND ?
			//		records[i].consensus : reverse_comp(records[i].consensus);
			string ssStr = records[i].consensus;
					
			if(gSNP && !gBISULFITE)
				ss = new SNPScoredSeq(ssStr, alignment_score, 
						 				abs_gen_pos, records[i].strand);
			else if(gBISULFITE || gATOG)
				//ss = new BSScoredSeq(ssStr, alignment_score, 
				ss = new BSScoredSeq(ssStr, alignment_score, 
										records[i].CHR_POS, records[i].strand);
			else
				ss = new NormalScoredSeq(ssStr, alignment_score, 
										records[i].CHR_POS, records[i].strand);
						
						
			if( (gBISULFITE || gATOG) && !pos_matrix ) {
				// Change the scoring matrices, then change them back
				if(records[i].strand == POS_STRAND) {
					if(gBISULFITE) {
						gALIGN_SCORES[(int)'c'][3] = gMATCH;
						gPHMM_ALIGN_SCORES[(int)'c'][3] = gPHMM_match;
					}
					if(gATOG) {
						gALIGN_SCORES[(int)'a'][2] = gMATCH;
						gPHMM_ALIGN_SCORES[(int)'a'][2] = gPHMM_match;
					}
				}
				else {
					if(gBISULFITE) {
						gALIGN_SCORES[(int)'g'][0] = gMATCH;
						gPHMM_ALIGN_SCORES[(int)'g'][0] = gPHMM_match;
					}
					if(gATOG) {
						gALIGN_SCORES[(int)'t'][1] = gMATCH;
						gPHMM_ALIGN_SCORES[(int)'t'][1] = gPHMM_match;
					}
				}
				ss->score(1.0, *gGen, len, *r, lock);		

				reset_alignment_matrices();
			}
			else {
				// When we score it, just use a denominator of 1 so it doesn't correct for
				// anything.
				ss->score(1.0, *gGen, len, *r, lock);
			}
			//fprintf(stderr,"Scoring at position %lu with alignment score of %f(%d)\n",	
				//records[i].CHR_POS,alignment_score,records[i].MAPQ);
			
			// Free up the record
			delete ss;
			// Free up the Read
			delete_read(r);
			
			total_records++;
			records[i].MAPQ = -1;
		} // end for

	} //end while(!finished)
}
	in.close();

	delete[] records;
	return total_records;
}

/**
 * Reads all the SAM files that are contained in 'fn', parsing through the string
 *
 * @return The total number of SAM records
 */
unsigned int ReadAllSam(const char* fn) {
	unsigned int total_records = 0;
	
	vector<string> files = splitFileNames(fn);
	
	if(gVERBOSE)
		fprintf(stderr,"Reading %lu SAM files\n",files.size());
	
	for(unsigned int i=0; i<files.size(); i++) {
		total_records += AddSamtoGenome(files[i].c_str());
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
bool HAS_QUIT = false;

void sig_handler(int signum) {
	fprintf(stderr,"Caught signal: %d\n",signum);
	if(signum == SIGINT && !HAS_QUIT) {
		fprintf(stderr,"\n\n----- SIGINT -----\n");
		fprintf(stderr,"Just received a signal to quit the program ");
		fprintf(stderr,"(likely through Ctrl-C)\n");
		fprintf(stderr,"\nIn the event this was an error, GNUMAP is printing genome output.\n");
		fprintf(stderr,"Press Ctrl-C again or issue SIGINT again to kill the program\n");
		fprintf(stderr,"------------------\n");

		HAS_QUIT = true;

		// Now we tell the program to finish printing output
		if(gBISULFITE || gATOG)	// Print every base
			gGen->PrintFinalSNP(output_file, false);
		else // Just do what's normal
			gGen->PrintFinal(output_file);

		return;
	}
	else if(signum == SIGINT) {
		fprintf(stderr,"\n\n----- SIGINT -----\nQuitting...\n---------------\n");
		exit(0);
	}

	// Print the backtrace
	void* bt[10];
	size_t size = backtrace(bt,10);
	backtrace_symbols_fd(bt, size, 2);

	exit(0);
}

int main(const int argc, const char* argv[]) {
	cerr << "This is GNUMAP's SAM2SGR, Version 2.0, for public and private use." << endl;

	//sigset_t sigmask;
	//sigemptyset(&sigmask);
	struct sigaction stop_action;
	stop_action.sa_handler = sig_handler;
	sigaction(SIGINT, &stop_action, NULL);

	double time_prog_start = When();
	
	InitProg();

	//set up alignment matrix
	// Move this code out of the main body
	setup_alignment_matrices();


	cout << endl << "Command Line Arguments:  ";
	for(int i=0; i<argc; i++)
		cout << argv[i] << ' ';
	cout << endl;
	
	gVERBOSE = 0;

	if(argc == 1) //means only the bin/gnumap parameter
		usage(0,"");

	if(argc < 3) {
		usage(1, "Need at least a genome and a SAM file to convert.");
		exit(1);
	}


	if(argc == 1) //means only the bin/gnumap parameter
		usage(0,"");

	if(argc < 3) {
		usage(1, "Need at least a genome and a SAM file to convert.");
		exit(1);
	}

	// Parse the command line
	int rc = ParseCmdLine(argc,argv);
	// an error occurred during option processing
	if(rc != 0) {
		GetParseError(cerr, argv);
		return -1;
	}

	if(pos_matrix != NULL) {
		try {
			readPWM(pos_matrix);
		}
		catch(const char* err) {
			cerr << "ERROR: \n\t" << err << endl;
			return -1;
		}
		catch(Exception *e) {
			cerr << "Error: \n\t" << e->GetMessage() << endl;
			delete e;
			return -1;
		}
	}
	
	ostringstream params;
	params << "Parameters: " << endl;
	params << "\tVerbose: " << gVERBOSE << endl;
	params << "\tGenome file(s): " << genome_file << endl;
	params << "\tOutput file: " << output_file << endl;
	params << "\tSAM file: " << sam_file << endl;
	params << "\tNumber of threads: " << gNUM_THREADS << endl;
	if(pos_matrix != NULL) {
		params << "\tUsing User-Defined Alignment Scores: " << endl;
		params << getPWM();
	}	
	else 
		params << "\tUsing Default Alignment Scores" << endl;
	params << "\tGap score: " << gGAP << endl;
	params << "\tMaximum Gaps: " << gMAX_GAP << endl;
	if(gSNP) {
		params << "\tEmploying SNP calling" << endl;
		params << "\t  SNP p-value is "<<gSNP_PVAL<<endl;
		if(gSNP_MONOP)
			params << "\t  Only allowing for monoploid SNPs"<<endl;
	}
	if(gBISULFITE) {
		gPRINT_FULL = true;
		params << "\tUsing BISULFITE conversion" << endl;
		g_bs_CONVERSION[(int)'c'] = g_bs_CONVERSION[(int)'C'] = 3;	//change it so every 'c' will be *hashed* as a 't'
		//gALIGN_SCORES[3][1] = MATCH;
		gALIGN_SCORES[(int)'c'][3] = gMATCH;
	}
	if(gATOG) {
		gPRINT_FULL = true;
		params << "\tUsing A to G conversion" << endl;
		g_bs_CONVERSION[(int)'a'] = g_bs_CONVERSION[(int)'A'] = 2;	//change it so every 'a' will be *hashed* as a 'g'
		//gALIGN_SCORES[2][0] = MATCH;
		gALIGN_SCORES[(int)'a'][2] = gMATCH;
	}
	if(gGEN_SIZE != 8) {
		if(gGEN_SIZE > 8)
			usage(1,"Invalid bin size (must be 8 or less)\n");
		params << "\tUsing irregular bin size of " << gGEN_SIZE << endl;
	}

#ifdef DEBUG
	printf("Alignment scores:\n");
	printf("\tA\tC\tG\tT\n");
	printf("A\t%f\t%f\t%f\t%f\n",gALIGN_SCORES[(int)'a'][0],gALIGN_SCORES[(int)'a'][1],
			gALIGN_SCORES[(int)'a'][2],gALIGN_SCORES[(int)'a'][3]);
	printf("C\t%f\t%f\t%f\t%f\n",gALIGN_SCORES[(int)'c'][0],gALIGN_SCORES[(int)'c'][1],
			gALIGN_SCORES[(int)'c'][2],gALIGN_SCORES[(int)'c'][3]);
	printf("G\t%f\t%f\t%f\t%f\n",gALIGN_SCORES[(int)'g'][0],gALIGN_SCORES[(int)'g'][1],
			gALIGN_SCORES[(int)'g'][2],gALIGN_SCORES[(int)'g'][3]);
	printf("T\t%f\t%f\t%f\t%f\n",gALIGN_SCORES[(int)'t'][0],gALIGN_SCORES[(int)'t'][1],
			gALIGN_SCORES[(int)'t'][2],gALIGN_SCORES[(int)'t'][3]);
	printf("N\t%f\t%f\t%f\t%f\n",gALIGN_SCORES[(int)'n'][0],gALIGN_SCORES[(int)'n'][1],
			gALIGN_SCORES[(int)'n'][2],gALIGN_SCORES[(int)'n'][3]);
#endif

	if(gVERBOSE) {
		cerr << params.str() << endl;
	}

	/************************************************/
	/* Finished initialization, now do the stuff!	*/
	/************************************************/
	
	omp_set_num_threads(gNUM_THREADS);

	
	//genome_file = argv[1];
	//sam_file = argv[2];
	//output_file = "testout.sgrexbak";
	unsigned int total_records = 0;
	try {
		gGen = new GENOME_t;
		gGen->use(genome_file);
		gGen->StoreGenome();
		
		total_records = ReadAllSam(sam_file);
		fprintf(stderr,"\nTotal Read: %u\n",total_records);
		
		fprintf(stderr,"\nFinished!  Printing final .sgr/.gmp file\n");
		if(gBISULFITE || gATOG)	// Print every base
			gGen->PrintFinalSNP(output_file, false);
		else // Just do what's normal
			gGen->PrintFinal(output_file);
			
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
	double temp_dbl;

	switch(*(param+1)) {	//switch on the first character
		case 'g':
			genome_file = assign;
			break;
		case 'o':
			output_file = assign;
			break;
		case 'c':
			if(sscanf(assign,"%d",&gNUM_THREADS) < 1)
				return PARSE_ERROR;
			break;
		case 'v':
			if(sscanf(assign,"%d",&gVERBOSE) < 1)
				return PARSE_ERROR;
			break;
		case 's':
			if(sscanf(assign,"%d",&gSTART_POS) < 1)
				return PARSE_ERROR;
			break;
		case 'h':
			if(sscanf(assign,"%d",&gMAX_HASH_SIZE) < 1)
				return PARSE_ERROR;
			break;
		case 'S':
			pos_matrix = assign;
			break;
		case 'G':
			if(sscanf(assign,"%lf",&temp_dbl) < 1)
				return PARSE_ERROR;
			gGAP = (float)temp_dbl;
			break;
		case 'M':
			if(sscanf(assign,"%d",&gMAX_GAP) < 1)
				return PARSE_ERROR;
			break;
		case '0':
			count--;
			gPRINT_FULL = 1;
			break;
		case 'b':
			count--;
			
			// set gSNP so it prints the full output
			gSNP = true;
			
			gBISULFITE = 1;
			gGEN_SIZE = 1;
			gPRINT_FULL = 1;
			break;
		case 'd':
			count--;
			
			// set gSNP so it prints the full output
			gSNP = true;
			
			gATOG = 1;
			gGEN_SIZE = 1;
			gPRINT_FULL = 1;
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
	PAR_NPROC,
	PAR_VERBOSE,
	PAR_START_POS,
	PAR_SUBST_FILE,
	PAR_GAP_PENALTY,
	PAR_ADAPTOR,
	PAR_PRINT_FULL,
	PAR_BS_SEQ,
	PAR_A_TO_G,
	PAR_GEN_SIZE,
	PAR_SNP,
	PAR_SNP_PVAL,
	PAR_SNP_MONOP
};

int set_arg_ext(const char* param, const char* assign, int &count) {
	double temp_dbl;
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
	else if(strncmp(param+2,"num_proc=",9) == 0) {
		which = PAR_NPROC;
		adjust = 10;
	}
	else if(strncmp(param+2,"verbose=",8) == 0) {
		//which = 'v';
		which = PAR_VERBOSE;
		adjust = 9;
	}
	else if(strncmp(param+2,"start_pos=",10) == 0) {
		which = PAR_START_POS;
		adjust = 11;
	}
	else if(strncmp(param+2,"subst_file=",9) == 0) {
		//which = 'S';
		which = PAR_SUBST_FILE;
		adjust = 10;
	}
	else if(strncmp(param+2,"gap_penalty=",12) == 0) {
		//which = 'G';
		which = PAR_GAP_PENALTY;
		adjust = 13;
	}
	else if(strncmp(param+2,"print_full",10) == 0) {
		//which = '0';
		which = PAR_PRINT_FULL;
		adjust = 12;
	}
	else if(strncmp(param+2,"bs_seq",6) == 0) {
		//which = 'b';
		which = PAR_BS_SEQ;
		adjust = 8;
	}
	else if(strncmp(param+2,"a_to_g",6) == 0) {
		//which = 'd';
		which = PAR_A_TO_G;
		adjust = 8;
	}
	else if(strncmp(param+2,"bin_size=",9) == 0) {
		//which = 3;
		which = PAR_GEN_SIZE;
		adjust = 10;
	}
	else if(strncmp(param+2,"snp",3) == 0) {
		which = PAR_SNP;
		adjust = 5;
	}
	else if(strncmp(param+2,"snp_pval",8) == 0) {
		which = PAR_SNP_PVAL;
		adjust = 10;
	}
	else if(strncmp(param+2,"snp_monop",9) == 0) {
		which = PAR_SNP_MONOP;
		adjust = 11;
	}
	else if(strcmp(param+2,"help") == 0) {
		usage(0,""); //usage will exit immediately
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
		case PAR_NPROC:
			if(sscanf(param,"%d",&gNUM_THREADS) < 1)
				return PARSE_ERROR;
			break;
		//case 'v':
		case PAR_VERBOSE:
			if(sscanf(param,"%d",&gVERBOSE) < 1)
				return PARSE_ERROR;
			break;
		//case 's':
		case PAR_START_POS:
			if(sscanf(param,"%d",&gSTART_POS) < 1)
				return PARSE_ERROR;
			break;
		//case 'S':
		case PAR_SUBST_FILE:
			pos_matrix = param;
			break;
		//case 'G':
		case PAR_GAP_PENALTY:
			if(sscanf(param,"%lf",&temp_dbl) < 1)
				return PARSE_ERROR;
			gGAP = (float)temp_dbl;
			break;
		//case 'A':
		case PAR_ADAPTOR:
			g_adaptor = new char[strlen(assign)];
			for(unsigned int i=0; i<strlen(assign); i++)
				g_adaptor[i] = (char)param[i];
			break;
		//case '0':
		case PAR_PRINT_FULL:
			gPRINT_FULL = 1;
			break;
		//case 'b':
		case PAR_BS_SEQ:			
			gBISULFITE = 1;
			gGEN_SIZE = 1;
			gPRINT_FULL = 1;
			break;
		//case 'd':
		case PAR_A_TO_G:						
			gATOG = 1;
			gGEN_SIZE = 1;
			gPRINT_FULL = 1;
			break;
		case PAR_GEN_SIZE:
			if(sscanf(param,"%d",&gGEN_SIZE) < 1)
				return PARSE_ERROR;
			break;
		case PAR_SNP:
			gSNP = true;
			gGEN_SIZE = 1;
			gPRINT_FULL = 1;
			break;
		case PAR_SNP_PVAL:
			if(sscanf(param,"%f",&gSNP_PVAL) < 1)
				return PARSE_ERROR;
			break;
		case PAR_SNP_MONOP:
			gSNP_MONOP = 1;
			break;
		default:
			return NO_MATCHING_ARG;
	}

	return 0;

}
