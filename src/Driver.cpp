/* Driver.cpp
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

#ifdef	MPI_RUN
#include <mpi.h>
#endif

#ifdef OMP_RUN
#include <omp.h>
#endif

#include <signal.h>
#include <vector>
#include <ctime>
#include <sys/time.h>
#include <sys/resource.h>
#include <cmath>
#include <set>
#include <map>
#include <pthread.h>
#include <unistd.h>

#include <algorithm>
#include <exception>
#include "const_include.h"
#include "const_define.h"	// Defines many things in this file
#include "GenomeBwt.h"
#include "bin_seq.h"
#include "SeqManager.h"
#include "ScoredSeq.h"
#include "BSScoredSeq.h"
#include "NormalScoredSeq.h"
#include "SNPScoredSeq.h"
//#include "HMMScoredSeq.h"
#include "SequenceOperations.h"
#include "centers.h"
#include "Exception.h"
#include "CudaDriver.h"

#ifdef DEBUG_NW
unsigned int num_nw[8] = {0,0,0,0,0,0,0,0};
#endif

typedef set<ScoredSeq*,ScoredSeq > seq_set;
typedef map<string,ScoredSeq*> seq_map;


const char* seq_file = 0;
const char* genome_file = 0;
const char* output_file = "gnumap_out";
ofstream of;
const char* pos_matrix = NULL;
bool gPRINT_ALL_SAM = false;

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

/*****************************************************
 * These next few are all for multi-threading
 *****************************************************/
unsigned int gNUM_THREADS = 1;
volatile bool setup_complete = false;	// global signal for complete setup
volatile bool finished = false;			// global signal for finishing everything
unsigned int *finished_arr = 0;					// Array used to tell when threads are done

double gTimeBegin;

pthread_attr_t attr;
pthread_mutex_t lock;
pthread_mutex_t nwlock;
pthread_mutex_t total_read_lock;
pthread_mutex_t cond_lock;
#ifdef DEBUG_TIME
double cond_lock_time = 0;
double read_lock_time = 0;
double write_lock_time = 0;
#endif

pthread_mutex_t comm_barrier_lock;
pthread_cond_t comm_barrier_cond = PTHREAD_COND_INITIALIZER;
volatile unsigned int cond_count = 0;
volatile unsigned int cond_thread_num;
//pthread_mutex_t read_lock;

// The writer thread
pthread_t writer;
pthread_mutex_t write_lock;

void* parallel_thread_run(void* t_opts);
struct thread_rets* omp_thread_run(struct thread_opts*);
void* mpi_thread_run(void* t_opts);
void getMoreReads(unsigned int &begin, unsigned int &end);
volatile unsigned int gReadsDone = 0;
volatile unsigned int iter_num = 0;
unsigned int nReadPos;
unsigned int nReadPerProc;

/*
const unsigned int nReads = 1024*8;
Read* gReadArray[nReads];
seq_map gReadLocs[nReads];
double gReadDenominator[nReads];
double gTopReadScore[nReads];
*/
unsigned int nReads;
Read** gReadArray;	// allocate dynamically
seq_map* gReadLocs;	// allocate dynamically
double* gReadDenominator;
double* gTopReadScore;

vector<TopReadOutput> gTopReadOutput;

#ifdef MPI_RUN
MPI::Status status;
#ifdef DISCRETIZE
/* center_e is relly a character (so it holds 128 bits) */
void MPISumCenters( const void *in, void *inout, int len, const MPI::Datatype &dptr ) {
	int i;
	center_d *a = (center_d*)in;
	center_d *b = (center_d*)inout;

//void MPISumCenters( const center_d *in, center_d *intou, int len, const MPI::Datatype &dptr ) {
//	int i;
	center_d c;

	for(i=0; i< len; i++) {
		// The sum has been pre-calculated and is just a lookup in this table
		c = sum_mat[*a][*b];
		*b = c;
		a++; b++;
	}
}
#endif // DISCRETIZE
#endif // MPI_RUN

typedef GenomeBwt GENOME_t;
GENOME_t gGen;
string cl_args;

CudaDriver* cudaDriver;

SeqManager* gSM;

void usage(int has_error, const char * errmessage) {
	if(has_error) cout << endl << errmessage << endl ;

	cout << "Usage: gnumap [options] <file_to_parse>\n"
		 << "  -g, --genome=STRING          Genome .fa file(s)\n"
		 << "  -o, --output=STRING          Output file\n"
		 << "  -v, --verbose=INT            Verbose (default=0)\n"
		 << "  -c, --num_proc=INT           Number of processors on which to run\n"
		 << "  -B, --buffer=INT             Buffer size\n"
		 << "  -?, --help                   Show this help message\n"
		 << "\n"
		 << "Options for Alignment\n"
		 << "  -a, --align_score=DOUBLE     Limit for sequence alignment (default: 90%)\n"
		 << "  -r, --raw                    Use raw score when determining alignment cutoff\n"
		 << "  -G, --gap_penalty=DOUBLE     Gap Penalty for Substitution Matrix\n"
		 << "  -A, --adaptor=STRING         Adaptor sequence to remove from sequences\n"
		 << "  -M, --max_gap=INT            Maximum Number of Gaps to use in Alignment\n"
		 << "  -S, --subst_file=STRING      Position-Weight Matrix file for DNA\n"
		 << "                               Substitutions\n"
		 << "  -T, --max_match=INT          Maximum number of matches for a given\n"
		 << "                               sequence (default: "<<gMAX_MATCHES<<")\n"
		 << "  -u, --unique_only            Only match sequences that map uniquely\n"
		 << "                               (default: not included)\n"
		 << "  -j, --jump=INT               The number of bases to jump in the sequence indexing\n"
		 << "                               (default: mer_size)\n"
		 << "  -k, --num_seed=INT           The total number of seed hits that must match to a\n"
		 << "                               location before it is considered for alignment\n"
		 << "                               (default: "<<gMIN_JUMP_MATCHES<<")\n"
		 << "  --up_strand                  Will only search the positive strand of the genome\n"
		 << "                               for matching location (will not look for reverse\n"
		 << "                               compliment match to the genome.\n"
		 << "  --down_strand                Will only search the negaitve strand (opposite of\n"
		 << "                               --up_strand command)\n"
		 << "  --no_nw                      This will disable the Needleman-Wunsch alignments and\n"
		 << "                               only use hit count as the basis for alignment. Score is\n"
		 << "                               calculated by summing the number of hits for a position.\n"
		 << "  --gpu                        This will enable the GPU to assist in the alignment.\n"
		 << "\n"
		 << "Options for Read Quality\n"
		 << "  -q, --read_quality=DOUBLE    Read quality cutoff:  won't align reads if they are\n"
		 << "                               below this cutoff value\n"
		 << "                               (default=0.0)\n"
		 << "  --illumina                   Defines the fastq file as Illumina file (otherwise\n"
		 << "                               does nothing)\n"
		 << "\n"
		 << "Options for creating index and genome\n"
		 << "  -m, --mer_size=INT           Mer size (default="<<DEF_MER_SIZE<<")\n"
		 << "  -s, --gen_skip=INT           Number of bases to skip when the genome is aligned\n"
		 << "  --bin_size=INT               The resolution for GNUMAP (default: 8)\n"
         << "  -h, --max_kmer=INT           Kmers that occur above this threshold are not used.\n"
		 << "\n"
		 << "Options for Printing:\n"
		 << "  -0, --print_full             Print locations for the entire sequence, not\n"
		 << "                               just for the beginning.\n"
		 << "  --print_all_sam              Include all possible SAM records in output\n"
		 << "  --vcf                        Prints VCF format instead of gmp (default: false)\n"
		 << "\n"
		 << "Options for extra GNUMAP runs\n"
		 << "  --snp                        Turn on SNP mapping (will output a .sgrex file)\n"
		 << "  --snp_pval=DOUBLE            P-Value cutoff for calling SNPs\n"
		 << "                               (default: "<<gSNP_PVAL<<")\n"
		 << "  --snp_monop                  Flag that turns on monoploid SNP calling\n"
		 << "                               (default: diploid SNP calling)\n"
		 << "  -b, --bs_seq                 Flag to turn on the C to T conversion, used in\n"
		 << "                               bisulfite sequence analysis\n"
		 << "  --b2                         Flag to turn on bisulfite sequencing, searching the reverse\n"
		 << "                               strand of -b\n"
		 << "  -d, --a_to_g                 Flag that allows for A to G conversion\n"
		 << "  --fast                       Perform a fast alignment (at some reduction\n"
		 << "                               of accuracy)\n"
		 << "                               (default: none)\n"
		 << "  --assembler                  Turns on flags for assembler output (currently in ALPHA)\n"
		 << "\n"
		 << "For MPI usage:\n"
		 << "  --MPI_largemem               If the run requires a large amount of memory, this\n"
		 << "                               flag will spread it accross several nodes.\n"
		 << "                               (default:  not included)\n"
		 << "\n";
	exit(1);
}

/*
struct GenomeLocation {
	GenomeLocation() {
		amount = 0.0;
		packed_char = 0;
	}

	GenomeLocation(GEN_TYPE pc, float amt)
		: amount(amt), packed_char(pc) {}

	float amount;
	GEN_TYPE packed_char;
};
*/

/************************************
 * MPI Structures and DataTypes 	*
 ************************************/
 #ifdef MPI_RUN
 /**
  * This is an MPI function for performing a reduce with user-defined data types.
  *
  * @param in One of the vectors of length /len/
  * @param inout The other vector to sum, which will contain the final result as well.
  * 		Also length /len/.
  * @param len The length of /in/ and /inout/
  * @param dptr The MPI Datatype of this function (not used?)
  */
 void sumGen( GenomeLocation *in, GenomeLocation *inout, int *len, MPI_Datatype *dptr ) {
 	int i;
 	
 	for(i=0; i<*len; i++) {
 		inout->amount = in->amount+inout->amount;
		
 		++in; ++inout;
 	}
 }
 
 #endif

/**
 * We have a loss of precision with some log additions. So we'll just add in the log space
 * and go from there.
 */
/*
inline float LogSum(float a, float b) {
	if(a > b && (a != 0)
		return a + log(1+exp(b-a));
	else
		return b + log(1 + exp(a-b));
}
*/
/**
 * This is the same as the above function, only it requires the input to be the powers
 * of e.  This way, we don't need to take the log of anything.
 */
/*inline float ExpLogSum(float expA, float expB) {
	
}
*/

/**
 * Because this is used to produce the consensus sequence, we need ambiguity characters...
 */
//inline char max_char(vector<double> &chr) {
inline char max_char(float* chr) {
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
 
//inline string GetConsensus(vector<vector<double> > &read) {
inline string GetConsensus(Read &read) {
	if(read.seq.size() > 0) {
		return read.seq;
	}

	string seq = "";

	for(unsigned int i=0; i<read.length; i++) {
		seq += max_char(read.pwm[i]);
	}

	return seq;
}

//I don't think this function is even used...
//vector<vector<double> > subvector(vector<vector<double> > &vec, int start, int length) {
float** subvector(float** &vec, int start, int length) {
	//vector<vector<double> > to_return;
	float** to_return = new float*[length];

	for(int i=start,j=0; i<start+length; i++,j++) {
		to_return[j] = new float[4];
		for(int k=0; k<4; k++) {
			to_return[j][k] = vec[i][k];
		}
		delete[] vec[i];
	}
	
	delete[] vec;	//we want to free the memory for vec because we're creating a new
					//one, right??
	
	return to_return;
}

struct thread_opts
{
	thread_opts(GENOME_t &g, SeqManager &s, int ti) :
		gen(g), sm(s), thread_id(ti) {}
	
	GENOME_t& gen;
	SeqManager &sm;
	//vector<vector<double> >* seqs;
	unsigned int num_seqs;
	int thread_id;
};

struct thread_rets {
	ostringstream* of;
	unsigned int good_seqs;
	unsigned int bad_seqs;
};


void clearMapAt(unsigned int index) {
	if(gReadLocs[index].empty()) {
		return;
	}
		
	seq_map::iterator sit;
	
	for(sit = gReadLocs[index].begin(); sit != gReadLocs[index].end(); ++sit) { 
		if((*sit).second)
			delete (*sit).second;
	}
	gReadLocs[index].clear();
}

#ifdef SET_POS
void clearPositions(GENOME_t &gen, vector<unsigned long> &positions, int thread_id) {
	for(unsigned int i=0; i<positions.size(); i++) {
		gen.unsetThreadID(positions[i],thread_id);
	}
}
#endif

#include "align_seq2_raw.cpp"

//string print_top_matches(GENOME_t &gen, vector<vector<double> > &search, string &consensus) {
//string print_top_matches(GENOME_t &gen, Read &search, string &consensus, 
//			unsigned int &num_not_matched, unsigned int &num_matched) {
void set_top_matches( GENOME_t &gen, unsigned int rIndex, string &consensus, 
			unsigned int &num_not_matched, unsigned int &num_matched, unsigned int &thread_id )
{

	bin_seq bs;
	double min_align_score = 0;
	double top_align_score = 0;
	double denominator = 0;
	

	Read* search = gReadArray[ rIndex ];
	seq_map* unique = &gReadLocs[ rIndex ];

	// Too short to align
	if( search->length < ( unsigned int ) gMER_SIZE )
    {
#ifdef DEBUG
		fprintf( stderr, "%s Too short to align\n", consensus.c_str() );
#endif
		num_not_matched++;
		gReadDenominator[ rIndex ] = 0;
		gTopReadScore[ rIndex ] = READ_TOO_SHORT;
		clearMapAt( rIndex );
		return;
	}


	// We need to do the length-1 because it will go to the base at that position.
    // @masakistan: it sets the maximum alignment score by doing a needleman-wunsch
    // with an exact match which can then be used to calculate the minimum alinment score needed
    // we don't need this if we're not doing needleman-wunsch
    if( gNW )
    {
        double max_align_score = bs.get_align_score( *search, consensus, 0, search->length - 1 );

        // Too low of quality to align
        if( max_align_score < gCUTOFF_SCORE )
        { 
#ifdef DEBUG
            fprintf( stderr, "%s:%s didn't meet cutoff: %f vs %f\n", consensus.c_str(), search->fq.c_str(),
                                max_align_score, gCUTOFF_SCORE );
#endif
            num_not_matched++;
            gReadDenominator[ rIndex ] = 0;
            gTopReadScore[ rIndex ] = READ_TOO_POOR;
            clearMapAt( rIndex );
            return;
        }

#ifdef DEBUG
        if(gVERBOSE > 2)
        {
            cerr << "\nNext Read ("<< search->name << ") aligned with " << max_align_score << ", less than " << gCUTOFF_SCORE << "?\n";
            cerr << consensus << "\n";
        }
#endif

        if( perc )
        {
            min_align_score = gALIGN_SCORE * max_align_score;
        }
        else
        {
            min_align_score = gALIGN_SCORE;
        }
    }
    else
    {
        // TODO: figure out what this value should actually be set to
        min_align_score = gMIN_JUMP_MATCHES;
    }

	// Match the positive strand
	if( gMATCH_POS_STRAND )
    {
        //cerr << "checking positive strand" << endl;
		//fprintf(stderr, "Aligning to UP_Strand\n" );
		bool aligned = align_sequence( gen, *unique, *search, consensus, 
									min_align_score, denominator, top_align_score,
									POS_STRAND, thread_id );
		if( !aligned )
        {

#ifdef DEBUG
			fprintf( stderr, "%s didn't align forward: too many!\n", search->name );
#endif
			clearMapAt( rIndex );
			
			num_matched++;
			
			gReadDenominator[ rIndex ] = 0;
			gTopReadScore[ rIndex ] = READ_TOO_MANY;
			return;
		}
	}

	/************************************************************/
	/* We need to do the reverse of this as well...				*/
	/************************************************************/
	// If we only want to do the positive strand, skip the reverse compliment stuff
	if( gMATCH_NEG_STRAND )
    {
        //cerr << "checking negative strand" << endl;
		//fprintf(stderr, "Aligning to DOWN_Strand\n" );
		//TODO: trim consensus before flipping and then send a query!
		//consensus = trim_end_string(consensus,search->length);
		string rev_consensus = reverse_comp( consensus );
		// Copy it so we don't need to reverse it again
		float** rev_pwm = reverse_comp_cpy( search->pwm, search->length );
		Read rc( rev_pwm, search->length );
		rc.name = search->name;

		//printReadPWM(&search);

        
#ifdef DEBUG
		if(gVERBOSE > 2)
        {
			cerr << "\tPerforming reverse compliment alignment with sequence " << consensus << " vs " << rev_consensus << endl;
		}
#endif
		

		bool aligned = align_sequence(gen, *unique, rc, rev_consensus, 
								min_align_score, denominator, top_align_score,
								NEG_STRAND, thread_id);

		// Clean up the copy
		for( unsigned int i = 0; i < rc.length; i++ )
        {
			delete[] rev_pwm[ i ];
		}
		delete[] rev_pwm;

#ifdef DEBUG
		if(gVERBOSE > 2)
        {
			cerr << "\tFound " << gReadLocs[ rIndex ].size() << " matches for this sequence\n";
		}
#endif
		if( !aligned )
        {
			// Too many sequences
#ifdef DEBUG
			fprintf( stderr, "%s didn't align backward: too many!\n", search->name );
#endif
			clearMapAt( rIndex );
			num_matched++;
			
			gReadDenominator[ rIndex ] = 0;
			gTopReadScore[ rIndex ] = READ_TOO_MANY;
			
			return;
		}
		
	}

#ifdef DEBUG
	fprintf( stderr, "**[%d/%d] Sequence [%s] has %lu matching locations\n", iproc, nproc, consensus.c_str(), unique->size() );
#endif

	if( gReadLocs[ rIndex ].size() == 0 )
    {
#ifdef DEBUG
        fprintf( stderr, "unmatched!\n" );
#endif
		num_not_matched++;
		gReadDenominator[ rIndex ] = 0;
		gTopReadScore[ rIndex ] = 0;
		return;
	}

	//fprintf(stderr,"Index %d has %u elements\n",rIndex,gReadLocs[rIndex].size());
	num_matched++;
	gReadDenominator[ rIndex ] = denominator;
#ifdef DEBUG
	fprintf( stderr, "Index %d has denominator of %f\n", rIndex, denominator );
#endif
	gTopReadScore[ rIndex ] = top_align_score;
	return;
}

void create_match_output(GENOME_t &gen, unsigned int rIndex, string &consensus) {

#ifdef DEBUG
	fprintf( stderr,"[%d/%d] [%s] has %lu matches\n", iproc, nproc, gReadArray[rIndex]->name, gReadLocs[rIndex].size() );
#endif

	if(gReadLocs[rIndex].size() == 0 )
    {
#ifdef DEBUG
		if(gTopReadScore[rIndex] != READ_TOO_POOR)
        {
			fprintf(stderr,"[%d/%d] Sequence [%s] does not have any matches\n", iproc, nproc, gReadArray[rIndex]->name);
        }
#endif
		return;
	}
	
	//to_return << endl;
	//to_return << "Denominator: " << denominator << endl;
	//seq_set::iterator sit;
	seq_map::iterator sit;

	ScoredSeq *max = new NormalScoredSeq;	// But empty...

	unsigned int len = gPRINT_FULL ? consensus.length() : 1;
	
	for(sit = gReadLocs[rIndex].begin(); sit != gReadLocs[rIndex].end(); ++sit) {

		(*sit).second->score(gReadDenominator[rIndex],gGen,len,*gReadArray[rIndex],lock);
/*
		vector<TopReadOutput> tro = (*sit).second->get_SAM(gReadDenominator[rIndex],*gReadArray[rIndex],consensus,gen,rIndex);
		for(unsigned int i=0; i<tro.size(); i++) {
			if(strcmp(tro[i].CHR_NAME, "chr2" ) == 0 )
				fprintf(stderr, "Found one here\n" );
		}
*/		
		// Here's the flag for printing a SAM record at every position
		if(gPRINT_ALL_SAM)
        {
			vector<TopReadOutput> out_vec = 
					//max->get_SAM(gReadDenominator[rIndex],*gReadArray[rIndex],consensus,gen,rIndex);
					(*sit).second->get_SAM(gReadDenominator[rIndex],*gReadArray[rIndex],consensus,gen,rIndex);
#ifdef DEBUG
			fprintf(stderr,"[%d/%d] \tFor sequence, %ld positions:\n",iproc,nproc,out_vec.size());
#endif
			MUTEX_LOCK(&write_lock);
			vector<TopReadOutput>::iterator vit = out_vec.begin();
			for(; vit != out_vec.end(); vit++) {
#ifdef DEBUG
				fprintf(stderr,"[%d/%d] Pushing back %s (%f)...size: %lu\n",iproc,nproc,vit->CHR_NAME,vit->POST_PROB, gTopReadOutput.size());
#endif
				gTopReadOutput.push_back(*vit);
			}
			MUTEX_UNLOCK(&write_lock);
		}
		
		
		// Free up the record
		if( (*(*sit).second).is_greater(*max)) {
//		if( (*(*sit).second).is_greater(*max, gReadDenominator[rIndex])) {
			//fprintf(stderr,"cur(%f) is greater than max(%f)\n",(*(*sit).second).get_score(),max->get_score());
			ScoredSeq *temp = max;
			max = (*sit).second;
			delete temp;
		}
		else 
			delete (*sit).second;
	}	
	//fprintf(stderr,"MAX IS %f and size is %d and DENOM is %f\n",max->get_score(), gReadLocs[rIndex].size(), gReadDenominator[rIndex]);
	
	// We just want to clear any spurious elements out of it
	gReadLocs[rIndex].clear();

	// We've already printed the max, so we'll go ahead and quit now
	if(gPRINT_ALL_SAM) {
		delete max;
		return;
	}
	
	// Only do this if I have one of the top reads...
	// we've got rounding error here (apparently)...
	if(max->get_score() > gTopReadScore[rIndex]-SAME_DIFF) {
		if(	gTopReadScore[rIndex] != READ_TOO_SHORT &&
			gTopReadScore[rIndex] != READ_TOO_POOR &&
			gTopReadScore[rIndex] != READ_TOO_MANY) {
			
			vector<TopReadOutput> out_vec = 
					max->get_SAM(gReadDenominator[rIndex],*gReadArray[rIndex],consensus,gen,rIndex);
#ifdef DEBUG_TIME
			double wait_begin = When();
#endif
			MUTEX_LOCK(&write_lock);
#ifdef DEBUG_TIME
			double wait_end = When();
			write_lock_time += wait_end - wait_begin;
#endif
			
			vector<TopReadOutput>::iterator vit = out_vec.begin();
			for(; vit != out_vec.end(); vit++) {
				//fprintf(stderr, "Pushing back...\n" );
				gTopReadOutput.push_back(*vit);
			}
			MUTEX_UNLOCK(&write_lock);
		}
		// We want to add this to the vector if this isn't the largemem option.  If it is, we'll
		// only add if this is the writer node
		else if(!gMPI_LARGEMEM || (iproc == 0 ) ) {
			TopReadOutput tro;
			tro.consensus = consensus;
			tro.qual = str2qual(*gReadArray[rIndex]);
			strncpy(tro.READ_NAME,gReadArray[rIndex]->name,MAX_NAME_SZ);
			//tro.MAPQ = gTopReadScore[rIndex];	// Flag for a bad read
			// Do we need to have the specific failure process?
			tro.MAPQ = READ_FAILED;

			//fprintf(stderr, "Setting it to READ_FAILED here\n" );

#ifdef DEBUG_TIME			
			double wait_begin = When();
#endif
			MUTEX_LOCK(&write_lock);
#ifdef DEBUG_TIME
			double wait_end = When();
			write_lock_time += wait_end - wait_begin;
#endif
			
			gTopReadOutput.push_back(tro);
			MUTEX_UNLOCK(&write_lock);
		}
	}
	//else
	//	fprintf(stderr,"Failed here: %f vs %f\n",max->get_score(),gTopReadScore[rIndex]);
	
	if(max->get_score() > gTopReadScore[rIndex]) {
		fprintf(stderr,"ERROR:  max score greater than TopReadScore: %f vs %f...\n\n",max->get_score(),gTopReadScore[rIndex]);
	}
	
	delete max;
	
}

/* GetPWM is used for getting a PWM matrix representing the match (m) and mismatch (mm)
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
 * Note:  GetPWM will modify the global gALIGN_SCORES variable.
 */
void readPWM(const char* fn) {
	gADJUST = 1;	// don't adjust the scores

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
				&& (tolower(*clabel) == 'c') )
        {
			if(gVERBOSE > 1)
            {
				printf("Matched here: %s\n",temp_chr);
            }
			in.getline(temp_chr,THIS_BUFFER);
			contains_labels=true;
		}
		if(gVERBOSE > 1)
        {
			cout << "Line: " << temp_chr << endl;
        }

		if(contains_labels)
        {
			if(sscanf(temp_chr,"%s %f %f %f %f",label,&a,&c,&g,&t) != 5) {
				char* error = new char[THIS_BUFFER];
				strcat(error, "Error in Score File: " );
				strcat(error,temp_chr);
				throw(error);
			}
		}
		else {
			if(sscanf(temp_chr,"%f %f %f %f",&a,&c,&g,&t) != 4) {
				char* error = new char[THIS_BUFFER];
				strcat(error, "Error in Score File: " );
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
	if(gVERBOSE > 1)
    {
		cout << "Matrix: " << endl;
		cout << "\t0\t1\t2\t3" << endl;
		for(int i=0; i<5; i++) {
			cout << i << "\t";
			for(int j=0; j<4; j++) 
				cout << temp_ALIGN_SCORES[i][j] << "\t";
			cout << endl;
		}
		cout << endl;
		//throw new Exception("Passed.");
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

string PWMtoString() {
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

void printReadPWM(Read* read) {
	printf("Name: %s\n",read->name);
	printf("Length: %d\n",read->length);
	for(unsigned int i=0; i<read->length; i++) {
		printf("%f %f %f %f\n",read->pwm[i][0],read->pwm[i][1],read->pwm[i][2],read->pwm[i][3]);
	}
}


void sig_handler(int signum) {
	fprintf(stderr,"Process-%d just caught signal %d\n",iproc, signum);

	// Print the backtrace
	void* bt[10];
	size_t size = backtrace(bt, 10 );
	backtrace_symbols_fd(bt, size, 2);

	exit(0);
}

/* 
 * We'll pull out this code from the main body into a separate function
 */
void setup_matrices() {
	memset(g_gen_CONVERSION,4, 256 );
	g_gen_CONVERSION[(int)'a'] = g_gen_CONVERSION[(int)'A'] = 0;
	g_gen_CONVERSION[(int)'c'] = g_gen_CONVERSION[(int)'C'] = 1;
	g_gen_CONVERSION[(int)'g'] = g_gen_CONVERSION[(int)'G'] = 2;
	g_gen_CONVERSION[(int)'t'] = g_gen_CONVERSION[(int)'T'] = 3;
	g_gen_CONVERSION[(int)'n'] = g_gen_CONVERSION[(int)'N'] = 4;
	g_gen_CONVERSION[(int)'\n'] = g_gen_CONVERSION[(int)'\r'] = 5;
	// Strange white-space characters
	g_gen_CONVERSION[10] = g_gen_CONVERSION[11] = g_gen_CONVERSION[12] = g_gen_CONVERSION[13] = 5;
	g_gen_CONVERSION[(int)'\0'] = 6;
	g_gen_CONVERSION[(int)'>'] = 7;

	memset(g_bs_CONVERSION,4, 256 );
	g_bs_CONVERSION[(unsigned int)'a'] = g_bs_CONVERSION[(unsigned int)'A'] = 0;
	g_bs_CONVERSION[(unsigned int)'c'] = g_bs_CONVERSION[(unsigned int)'C'] = 1;
	g_bs_CONVERSION[(unsigned int)'g'] = g_bs_CONVERSION[(unsigned int)'G'] = 2;
	g_bs_CONVERSION[(unsigned int)'t'] = g_bs_CONVERSION[(unsigned int)'T'] = 3;
	g_bs_CONVERSION[(unsigned int)'n'] = g_bs_CONVERSION[(unsigned int)'N'] = 4;
	g_bs_CONVERSION[(unsigned int)'\n'] = g_bs_CONVERSION[(unsigned int)'\r'] = 5;
	// Strange white-space characters
	g_bs_CONVERSION[10] = g_bs_CONVERSION[11] = g_bs_CONVERSION[12] = g_bs_CONVERSION[13] = 5;
	g_bs_CONVERSION[(unsigned int)'\0'] = 6;
	g_bs_CONVERSION[(unsigned int)'>'] = 7;

	memset(gINT2BASE,'?', 16 );
	gINT2BASE[0] = 'a';
	gINT2BASE[1] = 'c';
	gINT2BASE[2] = 'g';
	gINT2BASE[3] = 't';
	gINT2BASE[4] = 'n';
	gINT2BASE[5] = 'I';
	gINT2BASE[6] = 'D';


#ifdef DEBUG
    if( gVERBOSE > 1 )
    {
        printf("int\tbs\tgen\n");
        for(int i=0; i<256; i++) {
            printf("%d\t%c\t%d\t%d\n",i,(char)i,g_bs_CONVERSION[i],g_gen_CONVERSION[i]);
        }
        printf("int2base\n");
        for(int i=0; i<16; i++) {
            printf("%d\t%c\n",i,gINT2BASE[i]);
        }
    }
#endif

}

#include "a_matrices.c"

/**
 * Instead of having each of these be a constant, allocate them according to the number
 * of threads we are requesting.
 */
void alloc_nReads() {
	nReads = READS_PER_PROC * gNUM_THREADS;

	gReadArray = new Read*[nReads];
	gReadLocs = new seq_map[nReads];
	gReadDenominator = new double[nReads];
	gTopReadScore = new double[nReads];

	memset(gReadArray,0,nReads*sizeof(Read*));
	memset(gReadDenominator,0,nReads*sizeof(double));
	memset(gTopReadScore,0,nReads*sizeof(double));
}

/** 
 * Clean up after ourselves
 */
void clean_nReads() {
	delete[] gReadArray;
	delete[] gReadLocs;
	delete[] gReadDenominator;
	delete[] gTopReadScore;
}

int main(const int argc, const char* argv[]) {

#ifdef MPI_RUN
	// If we don't have MPI_THREAD_SERIALIZED, it won't allow us to run with threads
	int provided = MPI::Init_thread((int&)argc, (char**&)argv, MPI_THREAD_SERIALIZED);
	//MPI::Init_thread((int&)argc, (char**&)argv, MPI_THREAD_SERIALIZED);
	///*
	fprintf(stderr,"Provided is %d and MPI_THREAD_SERIALIZED is %d\n",provided,MPI_THREAD_SERIALIZED);
	fprintf(stderr,"Single: %d\n",MPI_THREAD_SINGLE);
	fprintf(stderr,"FUNNELED: %d\n",MPI_THREAD_FUNNELED);
	fprintf(stderr,"SERIALIZED: %d\n",MPI_THREAD_SERIALIZED);
	fprintf(stderr,"MULTIPLE: %d\n",MPI_THREAD_MULTIPLE);
	//*/
	nproc = MPI::COMM_WORLD.Get_size();
	iproc = MPI::COMM_WORLD.Get_rank();
	//fprintf(stderr,"iproc:%d nproc:%d\n",iproc,nproc);
#endif

	double prog_start_time = When();

	// We want to install a sig handler because this is throwing a bus error somewhere...
	struct sigaction bus_action;
	bus_action.sa_handler = sig_handler;
	sigaction(SIGBUS, &bus_action, NULL);

	//Print command line args
	cerr << "This is GNUMAP, Version "gVERSION", for public and private use." << endl;
	cerr << "# Using BWT/Suffix Array.\n";
	
	// move this code out of the main body.
	setup_matrices();

#ifdef DISCRETIZE
	// Set up the matrices used to discretize the read values
	init_centers();
	init_lookup();
	init_sums();
	//print_stats(stderr);
#endif

    cl_args = "";
	cout << endl << "Command Line Arguments:  ";
	for(int i=0; i<argc; i++)
    {
		cout << argv[i] << ' ';
        cl_args += argv[ i ];
        cl_args += " ";
    }
	cout << endl;
	
	gVERBOSE = 0;

	if(argc == 1) //means only the bin/gnumap parameter
    {
		usage(0, "" );
    }

	if(argc < 2)
    {
		usage(1, "Need at least a genome and a file to map.");
		exit(1);
	}

	int rc = ParseCmdLine(argc,argv);
	// an error occurred during option processing
	if(rc != 0)
    {
		GetParseError(cerr, argv);
		return -1;
	}


	// Move this code out of the main body
	setup_alignment_matrices();


	/* Manage errors that might occur */	
	// If they don't include a fasta file for sequences
	if( (seq_file == NULL) )
		usage(1, "Specify a single file e.g., sequences.fa\n");
	if(!genome_file)
		usage(1, "Specify a genome to map to with the -g flag\n");

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

	// Allocate space for dynammically-allocated arrays
	alloc_nReads();
	
	ostringstream params;
	params << "Parameters: " << endl;
	params << "\tVerbose: " << gVERBOSE << endl;
	params << "\tGenome file(s): " << genome_file << endl;
	params << "\tOutput file: " << output_file << endl;
	params << "\tMax matches: " << gMAX_MATCHES << endl;
    if(gPRINT_ALL_SAM)
    {
		params << "\t\tPrinting all SAM records\n";
    }

	params << "\tSequence file(s): " << seq_file << endl;
	
    if(!perc)
    {
		params << "\tAlign score: " << gALIGN_SCORE << endl;
    }
	else
    {
		params << "\tAlign percentage: " << gALIGN_SCORE*100 << "%" << endl;
    }
	
    if(gCUTOFF_SCORE != 0)
    {
		if(gCUTOFF_SCORE < 0)
        {
			usage(1, "-q: Invalid cutoff score.  Must be >= 0");
        }

		params << "\tUsing cutoff score of " << gCUTOFF_SCORE << endl;
	}

    if( gNW )
    {
        params << "\tUsing Needleman-Wunsch alignments" << endl;
    }
    else
    {
        params << "\tUsing kmer hit counts as proxy for full sequence alignment" << endl;
    }

#ifdef MPI_RUN
	// You must compile the MPI library to work with multiple threads otherwise
	// it will actually be SLOWER with multiple threads than a single one (it will
	// default to MPI_THREAD_FUNNELED or MPI_THREAD_SINGLE)
	if(provided != MPI_THREAD_SERIALIZED && provided != MPI_THREAD_MULTIPLE) {
		cerr << "\nWarning:  MPI libraries do not support multiple threads. "
			 << "Turning this option off\n\n";
		gNUM_THREADS = 1;
	}
#endif
	params << "\tNumber of threads: " << gNUM_THREADS << endl;
	if(gFAST)
    {
		if(!gMER_SIZE)	//if the user didn't specify one
        {
			gMER_SIZE = 14;
        }

		gJUMP_SIZE = gMER_SIZE;
		params << "\tUsing FAST alignment mode" << endl;
	}

	if(gGEN_SKIP != 0)
    {
		if(gGEN_SKIP < 0)
        {
			usage(1, "-s: Invalid Genome Skip size\n");
        }

		params << "\tSkipping " << gGEN_SKIP << " bases during indexing step" << endl;
	}
	
    gGEN_SKIP++;	// We're mod'ing in the function, so we want to mod by GEN_SKIP+1
    
    if(gMAX_KMER_SIZE > 0)
    {
		params << "\tLargest Kmer Size: " << gMAX_KMER_SIZE << endl;
    }
	
    if(!gMER_SIZE)	//if the user didn't specify one
    {
		gMER_SIZE = DEF_MER_SIZE;
    }
	
    params << "\tMer size: " << gMER_SIZE << endl;
	
    if(gJUMP_SIZE)
    {
		if(gJUMP_SIZE < 1)
        {
			usage(1, "-j/--jump: Invalid jump size\n");
        }
	}
	else
    {
		gJUMP_SIZE = gMER_SIZE/2;
	}
	
    params << "\tUsing jump size of " << gJUMP_SIZE << endl;
	
    if(gMIN_JUMP_MATCHES < 1)
    {
		usage(1, "-k/--num_seed: Invalid matching seed number\n");
    }

	params << "\tUsing min seed matches for each sequence of " << gMIN_JUMP_MATCHES << endl;

	if(pos_matrix != NULL)
    {
		params << "\tUsing User-Defined Alignment Scores: " << endl;
		params << PWMtoString();
	}	
	else
    {
		params << "\tUsing Default Alignment Scores" << endl;
    }

	params << "\tGap score: " << gGAP << endl;
	params << "\tMaximum Gaps: " << gMAX_GAP << endl;
	
	params << "\tUsing GPU: ";
	if(gpu) {
		params << "True" << endl;
	}
	else {
		params << "False" << endl;
	}
	
    if(g_adaptor)
    {
		params << "\tUsing Adaptor Sequence: " << g_adaptor << endl;
	}

	if(gSNP)
    {
		params << "\tEmploying SNP calling" << endl;
		params << "\t  SNP p-value is "<<gSNP_PVAL<<endl;
		if(gSNP_MONOP)
        {
			params << "\t  Only allowing for monoploid SNPs"<<endl;
        }
	}

	if(gBISULFITE)
    {
		gPRINT_FULL = true;
		params << "\tUsing BISULFITE conversion" << endl;
		//gALIGN_SCORES[3][1] = MATCH;

		// Bisulfite conversion only works with one strand, so don't let them have both
		/*
		 * We don't want to match to only one strand:  All 4, BSW,BSWR,BSC,BSCR
		if(gMATCH_POS_STRAND && gMATCH_NEG_STRAND)
			usage(1, "Please specify a strand to match with --up_strand or --down_strand");
			// exit
		*/

		if(gMATCH_POS_STRAND && !gBISULFITE2)
        {
			g_bs_CONVERSION[(int)'c'] = g_bs_CONVERSION[(int)'C'] = 3;	//change it so every 'c' will be *indexed* as a 't'
			// If they've supplied a matrix, we don't want to mess with this
			if(pos_matrix == NULL)
            {
				gALIGN_SCORES[(int)'c'][3] = gMATCH;
            }
		}
		else // we don't want to match to the positive strand, by default only mapping to one strand
        {
			// when doing the opposite strand, we want to allow 'G->A' (the rev comp)
			g_bs_CONVERSION[(int)'g'] = g_bs_CONVERSION[(int)'G'] = 0;	//change it so every 'g' will be *indexed* as a 'a'
			// If they've supplied a matrix, we don't want to mess with this
			if(pos_matrix == NULL)
            {
				gALIGN_SCORES[(int)'g'][0] = gMATCH;
            }

		}
				
	}

	if(gATOG)
    {
		gPRINT_FULL = true;
		params << "\tUsing A to G conversion" << endl;
	
		/*
		 * We don't want to match to only one strand:  All 4, REW,REWR,REC,RECR
		// RNA editing only works with one strand, so don't let them have both
		if(gMATCH_POS_STRAND && gMATCH_NEG_STRAND)
			usage(1, "Please specify a strand to match with --up_strand or --down_strand");
			// exit
		 */

		if(gMATCH_POS_STRAND)
        {	
			g_bs_CONVERSION[(int)'a'] = g_bs_CONVERSION[(int)'A'] = 2;	//change it so every 'a' will be *indexed* as a 'g'
			// If they've supplied a matrix, we don't want to mess with this
			if(pos_matrix == NULL)
            {
				gALIGN_SCORES[(int)'a'][2] = gMATCH;
            }
		}
		else // only if we don't want to match to both strands...
        {
			// when doing the opposite strand, we want to allow 'T->C' (the rev comp)
			g_bs_CONVERSION[(int)'t'] = g_bs_CONVERSION[(int)'T'] = 1;	//change it so every 't' will be *indexed* as a 'c'
			// If they've supplied a matrix, we don't want to mess with this
			if(pos_matrix == NULL)
            {
				gALIGN_SCORES[(int)'t'][1] = gMATCH;
            }
		}

	}

	if(gGEN_SIZE != 8)
    {
		if(gGEN_SIZE < 1)
        {
			usage(1, "Invalid bin size (must be more than 0)\n" );
        }
		if(gGEN_SIZE > 8)
        {
			usage(1, "Invalid bin size (must be 8 or less)\n" );
        }
		params << "\tUsing irregular bin size of " << gGEN_SIZE << endl;
	}

	if(!gMATCH_NEG_STRAND)
    {
		params << "\tOnly matching to positive strand of genome\n";
	}

	if(!gMATCH_POS_STRAND)
    {
		params << "\tOnly matching to negative strand of genome\n";
	}

#ifdef DEBUG
    if( gVERBOSE > 0 )
    {
        printf("Alignment scores:\n");
        printf("\tA\tC\tG\tT\n");
        printf("A\t%f\t%f\t%f\t%f\n",1.0/gADJUST*gALIGN_SCORES[(int)'a'][0],1.0/gADJUST*gALIGN_SCORES[(int)'a'][1],
                1.0/gADJUST*gALIGN_SCORES[(int)'a'][2],1.0/gADJUST*gALIGN_SCORES[(int)'a'][3]);
        printf("C\t%f\t%f\t%f\t%f\n",1.0/gADJUST*gALIGN_SCORES[(int)'c'][0],1.0/gADJUST*gALIGN_SCORES[(int)'c'][1],
                1.0/gADJUST*gALIGN_SCORES[(int)'c'][2],1.0/gADJUST*gALIGN_SCORES[(int)'c'][3]);
        printf("G\t%f\t%f\t%f\t%f\n",1.0/gADJUST*gALIGN_SCORES[(int)'g'][0],1.0/gADJUST*gALIGN_SCORES[(int)'g'][1],
                1.0/gADJUST*gALIGN_SCORES[(int)'g'][2],1.0/gADJUST*gALIGN_SCORES[(int)'g'][3]);
        printf("T\t%f\t%f\t%f\t%f\n",1.0/gADJUST*gALIGN_SCORES[(int)'t'][0],1.0/gADJUST*gALIGN_SCORES[(int)'t'][1],
                1.0/gADJUST*gALIGN_SCORES[(int)'t'][2],1.0/gADJUST*gALIGN_SCORES[(int)'t'][3]);
        printf("N\t%f\t%f\t%f\t%f\n",1.0/gADJUST*gALIGN_SCORES[(int)'n'][0],1.0/gADJUST*gALIGN_SCORES[(int)'n'][1],
                1.0/gADJUST*gALIGN_SCORES[(int)'n'][2],1.0/gADJUST*gALIGN_SCORES[(int)'n'][3]);
    }
#endif

	//if(gVERBOSE)
    {
		cerr << params.str() << endl;
	}
	
	// This works for MPI and for everything else
	if(nproc > 1)
    {
		// rename the output file to have the processor number after it
		char* output_ext = new char[strlen(output_file)+10];
		sprintf(output_ext,"%s_%d",output_file,iproc);
		output_file = output_ext;

	}

	char* sam_file = new char[strlen(output_file)+10];
	sprintf(sam_file,"%s.sam",output_file);

	// Everyone opens up their output file
	of.open(sam_file);
	// Check to make sure the output file/directory exists
	if(of.bad() || of.fail())
    {
		cerr << "ERROR:\n\t";
		perror("Error in reading output file");
		return -1;
	}
	
	// Don't write anything to the output file
	/*
	// only write the params if you're the first node		
	if(iproc == 0 ) {
		of << "<<HEADER\n";
		of << params.str();
		of << "\nHEADER\n";
	}
	*/

	if(gVERBOSE)
    {
		cerr << "Indexing the genome." << endl;
    }
	// Use the global genome: gGen
	//GENOME_t gen; 
	unsigned long my_start=0, my_end=0;
	
    	if( gpu ) {
		if(!CudaDriver::initDevice()) {
			// There are no CUDA devices found, so exit with an error
			fprintf(stderr, "GNUMAP will run without CUDA devices.");
			gpu = false;
		}
		else {
			cudaDriver = new CudaDriver();
			cudaDriver->getAvailableMem();
		}
	}

	try {

		//unsigned long gen_size=0;
		

#ifdef MPI_RUN
		// don't need to do this if we're only reading in a file
		if(!gREAD) {
			if(gMPI_LARGEMEM && iproc == 0 ) {
				GENOME_t notherGen;
				// we don't care about these parameters.  Just read the whole thing
				notherGen.use(genome_file);
				gen_size = notherGen.count();
				//fprintf(stderr,"Size of genome here: %lu\n",gen_size);
			}

			//fprintf(stderr,"Size of genome after here: %lu\n",gen_size);
			
			my_start = 0;
			my_end = gen_size;
			// tell each node what part they need to get
			//MPI_Broadcast((char*)&my_start,sizeof(my_start),MPI_CHAR,
			MPI::COMM_WORLD.Bcast(&gen_size, 1, MPI_UNSIGNED_LONG, 0);

			// If we want to have large memory MPI, we should have declared the flag.
			if(gMPI_LARGEMEM) {
				my_start = gen_size/nproc * iproc;
				my_end = my_start + gen_size/nproc;
			}
			//return 0;
			
			fprintf(stderr,"[%d/%d] gen_size=%lu, my_start=%lu, my_end=%lu\n",iproc,nproc,gen_size,my_start,my_end);
		}
#endif

		gGen.use(genome_file, my_start, my_end);
		gGen.LoadGenome();

	}
	catch(const char* err) {
		cerr << "ERROR: \n\t" << err << endl;
		return -1;
	}
	catch(bad_alloc&) {
		perror("ERROR HERE");
		cerr << "ERROR: \n\tUnable to allocate enough memory.\n\t(Pick a smaller mer size?)" << endl;
		cerr << "\t\tTime: " << When()-prog_start_time << endl;
		return -1;
	}
	catch(exception &e) {
		cerr << "ERROR: \n\t" << e.what() << endl;
		return -1;
	}
	catch(Exception *e) {
		fprintf(stderr,"ERROR: \n\t%s\n",e->GetMessage());
	}
	catch(...) {
		cerr << "ERROR: \n\tUnknown error." << endl;
		return -1;
	}

	if(gVERBOSE)
    {
		fprintf(stderr,"\nTime to index: %.0f seconds\n",When()-prog_start_time);
    }
	
	// Don't want to print any header information in this file
	//if(iproc == 0 )
		// Output SAM by default
		//of << "# <QNAME>\t<FLAG>\t<RNAME>\t<POS>\t<MAPQ>\t<CIGAR>\t<MRNM>\t<MPOS>\t<ISIZE>\t<SEQ>\t<QUAL>\t<OPT>\n";

	gSM = new SeqManager(gReadArray, seq_file, nReads, gNUM_THREADS);

	/*
	if(!iproc)
	{
		int i=0;
		char hostname[256];
		gethostname(hostname, sizeof(hostname));
		fprintf(stderr,"PID %d on %s ready for attach\n",getpid(),hostname);
		while(0 == i)
			sleep(5);
	}
	*/

	/************************************************************************/
	/* Here is where we'll begin the process of matching and writing to the	*/
	/* output file.															*/
	/************************************************************************/	
	unsigned int seqs_matched = 0;
	unsigned int seqs_not_matched = 0;
	cond_thread_num = gNUM_THREADS;
	setup_complete = true;

	/************************************
	 * Do the threading					*
	 ************************************/
	
	pthread_t pthread_hand[gNUM_THREADS];
	void *pthread_ret[gNUM_THREADS];
	pthread_mutex_init(&lock,NULL);
	pthread_mutex_init(&nwlock,NULL);
	pthread_mutex_init(&write_lock,NULL);
	pthread_mutex_init(&cond_lock,NULL);
	pthread_mutex_init(&total_read_lock,NULL);
	pthread_mutex_init(&comm_barrier_lock,NULL);
	
	
	//Here is where we set the threading information.
	pthread_setconcurrency(gNUM_THREADS);
	pthread_attr_init(&attr);
	pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

	//An array of pointers so we can delete them when we're done.
	thread_opts* t_o_ptr[gNUM_THREADS];
	//Create an array to signify when we're done
	finished_arr = new unsigned int[gNUM_THREADS];
	for(unsigned int i=0; i<gNUM_THREADS; i++) {
		finished_arr[i] = 0;		
	}
	
#ifdef OMP_RUN
    cerr << "open mp" << endl;
	cond_thread_num = 1;
	omp_set_num_threads(gNUM_THREADS);
	thread_opts* t_o = new thread_opts(gGen,*gSM, 0 );
	thread_rets* ret = omp_thread_run(t_o);

	unsigned int seqs_matched = 0, seqs_not_matched = 0;
	seqs_matched = ret->good_seqs;
	seqs_not_matched = ret->bad_seqs;

#else	// If OMP_RUN not defined
	for(unsigned int i=0; i<gNUM_THREADS; i++) {

		//The parameters to pass to the thread
		thread_opts* t_o = new thread_opts(gGen,*gSM,i);
		t_o_ptr[i] = t_o;

#ifdef MPI_RUN
		if(gMPI_LARGEMEM) 
        {
            //cerr << "mpi large mem" << endl;
			pthread_create(&pthread_hand[i], NULL, mpi_thread_run, (void*)t_o);
        }
		else
        {
            //cerr << "parallel thread" << endl;
			pthread_create(&pthread_hand[i], NULL, parallel_thread_run, (void*)t_o);
        }

#else //MPI_RUN not defined
        //cerr << "parallel thread 2" << endl;
		pthread_create(&pthread_hand[i], NULL, parallel_thread_run, (void*)t_o);
#endif //end MPI_RUN
	}
	
	
	/********************************************************/
	/* wait for all the threads to join back together... 	*/
	/********************************************************/
	for (unsigned int i=0; i<gNUM_THREADS; i++)
    {
		pthread_join(pthread_hand[i], &pthread_ret[i]);
	}
	
	/********************************************************/
	/* clean up the memory leaks stuff					 	*/
	/********************************************************/
	for(unsigned int i=0; i<gNUM_THREADS; i++)
    {
		//delete[] seq_ptr[i];
		delete t_o_ptr[i];
	}
	
	for(unsigned int i=0; i<gNUM_THREADS; i++)
    {
		seqs_matched += ((thread_rets*)pthread_ret[i])->good_seqs;
		seqs_not_matched += ((thread_rets*)pthread_ret[i])->bad_seqs;

		delete (thread_rets*)pthread_ret[i];
	}
#endif //end OMP_RUN

#if defined(DEBUG_NW)
    unsigned int total_nw = 0;
    unsigned int max_nw = 0;
    unsigned int min_nw = num_nw[ 0 ];

    for( int i = 0; i < nproc; i++ )
    {
        total_nw += num_nw[ i ];

        if( num_nw[ i ] > max_nw )
        {
            max_nw = num_nw[ i ];
        }
        if( num_nw[ i ] < min_nw )
        {
            min_nw = num_nw[ i ];
        }
    }

    fprintf( stderr, "Total NW: %u, Max NW: %u, Min NW: %u, DIFF: %u\n", total_nw, max_nw, min_nw, max_nw - min_nw );
#endif

	
	if(gVERBOSE)
    {
        cerr << endl;
    }
		
	
	if(gVERBOSE)
    {
		cerr << "\n[" << iproc << "/-] Time since start: " << When()-prog_start_time << endl;
    }
	if(gVERBOSE)
    {
		cerr << "\n[" << iproc << "/-] Printing output." << endl;
    }

#ifdef MPI_RUN
	// If it's not the largemem (where we each print out our genome), we need to reduce
	if(!gMPI_LARGEMEM && nproc > 1)
    {		
		fprintf(stderr,"[-/%d] Sending genomes...",iproc);
		float* gen_ptr = gGen.GetGenomeAmtPtr();
		unsigned int gen_size = (unsigned int)gGen.size()/gGEN_SIZE;

		/****************************************************************
		 * We only want to create a portion of the genome, bit by bit
		 ****************************************************************/
		//float* sendGenome = (float*)malloc(sizeof(float)*gen_size);
		size_t commGenSize = gen_size < gMEM_BUFFER_SIZE ? gen_size : gMEM_BUFFER_SIZE;
#ifdef INTDISC
		// We need a larger space for sending genomes from the integer discretization method
		float* sendGenome = (float*)calloc(sizeof(float),commGenSize*NUM_SNP_VALS);
		//fprintf(stderr, "%d - Successfully allocated genomes\n", iproc);
#else //INTDISC
		float* sendGenome = (float*)calloc(sizeof(float),commGenSize);
#endif //INTDISC
		if(!sendGenome) {
			fprintf(stderr,"ERROR!!! Could not create genome of size %lu bytes!!!\n",sizeof(float)*commGenSize);
			assert(sendGenome && "sendGenome could not be initialized");
		}
//#ifdef DEBUG
		else {
#ifdef INTDISC
			fprintf(stderr,"Successfully created genome of size %lu bytes\n",sizeof(float)*commGenSize*NUM_SNP_VALS);
#else //INTDISC
			fprintf(stderr,"Successfully created genome of size %lu bytes\n",sizeof(float)*commGenSize);
#endif //INTDISC
		}
//#endif //DEBUG

#ifdef DISCRETIZE
		// Define an MPI datatype to pass
		MPI::Datatype MPIctype = MPI::UNSIGNED_CHAR.Create_contiguous(1);
		MPIctype.Commit();

		// Create the operation for summing centers
		MPI::Op MPImyOp;
		MPImyOp.Init( MPISumCenters, true );
#endif //DISCRETIZE

		//for(size_t i=0; i<gGen.size(); i+=gMEM_BUFFER_SIZE) {
		for(size_t i=0; i<gen_size; i+=gMEM_BUFFER_SIZE) {
			size_t remaining = gen_size - i;
			size_t spots_used = remaining < gMEM_BUFFER_SIZE ?
				remaining : gMEM_BUFFER_SIZE;
#ifdef DEBUG			
			fprintf(stderr,"[%d/%d] For this round, starting at %lu and going for %lu spots, and commGenSize is %lu\n",
					iproc,nproc,i,spots_used,commGenSize*NUM_SNP_VALS);
			fprintf(stderr,"[%d/%d] There are %lu spots remaining, buffer of %u\n",
					iproc,nproc,remaining,gMEM_BUFFER_SIZE);
#endif //DEBUG
			//fprintf(stderr, "[%d] BEF: %f\n", iproc, gen_ptr[500035]);
			memcpy(sendGenome,gen_ptr+i,(size_t)(spots_used*sizeof(float)));
			MPI::COMM_WORLD.Allreduce(sendGenome, gen_ptr+i, spots_used, MPI::FLOAT, MPI_SUM);
			//MPI::COMM_WORLD.Reduce(sendGenome, gen_ptr+i, spots_used, MPI::FLOAT, MPI_SUM, 0);
			//fprintf(stderr, "[%d] AFT: %f\n", iproc, gen_ptr[500035]);

			if(gSNP || gBISULFITE || gATOG) {
#ifdef DISCRETIZE
				if(gVERBOSE > 0)
                {
					fprintf(stderr, "Getting character allocations...\n" );
                }
				center_d* read_ptr = gGen.GetGenomeAllotPtr();
				
				// We already have summed our amount pointer, so adjust our weights and then sum them
				//for(unsigned int j=0; j<gen_size; j++) {
				for(unsigned int j=0; j<spots_used; j++) {
					//fprintf(stderr,"[%d: %u] ",iproc,j);
					//char prev = read_ptr[j+i];
					// The previous values exist in sendGenome
					read_ptr[j+i] = adjust_center(sendGenome[j], gen_ptr[j+i], read_ptr[j+i]);
				}
			
				// To reduce, use the space already allocated in sendGenome
				center_d *sendReads = (center_d*)sendGenome;
				memcpy(sendReads, read_ptr+i, spots_used*sizeof(center_d));
				//MPI::COMM_WORLD.Reduce((void*)sendReads, (void*)read_ptr, gen_size, MPIctype, MPImyOp, 0);
				MPI::COMM_WORLD.Reduce((void*)sendReads, (void*)(read_ptr+i), spots_used, MPIctype, MPImyOp, 0);

#elif defined( INTDISC )
				unsigned char* read_ptr = gGen.GetGenomeAllotPtr();
				
				// We already have summed our amount pointer, so adjust our weights and then sum them
				//for(unsigned int j=0; j<gen_size; j++)
				for(unsigned int j=0; j<spots_used; j++)
					adjust255(sendGenome[j], gen_ptr[i+j], read_ptr + (i+j)*NUM_SNP_VALS);
				
				unsigned char* sendReads = (unsigned char*)sendGenome;
				// sizeof(float)=4*sizeof(char), so we can fit all four chars in this one array
				memcpy(sendReads, read_ptr+i*NUM_SNP_VALS, (spots_used*NUM_SNP_VALS)*sizeof(unsigned char));
				//fprintf(stderr, "[%d] BEF: %u %u %u %u %u\n", iproc, read_ptr[500035*NUM_SNP_VALS], read_ptr[500035*NUM_SNP_VALS+1], read_ptr[500035*NUM_SNP_VALS+2], read_ptr[500035*NUM_SNP_VALS+3], read_ptr[500035*NUM_SNP_VALS+4]);
				MPI::COMM_WORLD.Reduce(sendReads, read_ptr+i*NUM_SNP_VALS, spots_used*NUM_SNP_VALS, MPI::UNSIGNED_CHAR, MPI_SUM, 0);
				//fprintf(stderr, "[%d] AFT: %u %u %u %u %u\n", iproc, read_ptr[500035*NUM_SNP_VALS], read_ptr[500035*NUM_SNP_VALS+1], read_ptr[500035*NUM_SNP_VALS+2], read_ptr[500035*NUM_SNP_VALS+3], read_ptr[500035*NUM_SNP_VALS+4]);
#else  //not defined INTDISC or DISCRETIZED
				if(gVERBOSE > 0)
                {
					fprintf(stderr, "A..." );
                }
				// Do the read A's now
				gen_ptr = gGen.GetGenomeAPtr();
				if(gen_ptr)
                {
					memcpy(sendGenome,&gen_ptr[i],spots_used*sizeof(float));
					MPI::COMM_WORLD.Reduce(sendGenome, &gen_ptr[i], spots_used, MPI::FLOAT, MPI_SUM, 0);
				}
				if(gVERBOSE > 0)
                {
					fprintf(stderr, "C..." );
                }
				// C's
				gen_ptr = gGen.GetGenomeCPtr();
				if(gen_ptr)
                {
					memcpy(sendGenome,&gen_ptr[i],spots_used*sizeof(float));
					MPI::COMM_WORLD.Reduce(sendGenome, &gen_ptr[i], spots_used, MPI::FLOAT, MPI_SUM, 0);
				}
				if(gVERBOSE > 0)
                {
					fprintf(stderr, "G..." );
                }
				// G's
				gen_ptr = gGen.GetGenomeGPtr();
				if(gen_ptr)
                {
					memcpy(sendGenome,&gen_ptr[i],spots_used*sizeof(float));
					MPI::COMM_WORLD.Reduce(sendGenome, &gen_ptr[i], spots_used, MPI::FLOAT, MPI_SUM, 0);
				}
				if(gVERBOSE > 0)
                {
					fprintf(stderr, "T..." );
                }
				// T's
				gen_ptr = gGen.GetGenomeTPtr();
				if(gen_ptr)
                {
					memcpy(sendGenome,&gen_ptr[i],spots_used*sizeof(float));
					MPI::COMM_WORLD.Reduce(sendGenome, &gen_ptr[i], spots_used, MPI::FLOAT, MPI_SUM, 0);
				}
				if(gVERBOSE > 0)
                {
					fprintf(stderr, "N..." );
                }
				// N's
				gen_ptr = gGen.GetGenomeNPtr();
				if(gen_ptr)
                {
					memcpy(sendGenome,&gen_ptr[i],spots_used*sizeof(float));
					MPI::COMM_WORLD.Reduce(sendGenome, &gen_ptr[i], spots_used, MPI::FLOAT, MPI_SUM, 0);
				}
#ifdef _INDEL
				if(gVERBOSE > 0)
                {
					fprintf(stderr, "Insertion..." );
                }
				// Insertions's
				gen_ptr = gGen.GetGenomeIPtr();
				if(gen_ptr)
                {
					memcpy(sendGenome,&gen_ptr[i],spots_used*sizeof(float));
					MPI::COMM_WORLD.Reduce(sendGenome, &gen_ptr[i], spots_used, MPI::FLOAT, MPI_SUM, 0);
				}
				if(gVERBOSE > 0)
                {
					fprintf(stderr, "Deletion..." );
                }
				// Deletions's
				gen_ptr = gGen.GetGenomeDPtr();
				if(gen_ptr)
                {
					memcpy(sendGenome,&gen_ptr[i],spots_used*sizeof(float));
					MPI::COMM_WORLD.Reduce(sendGenome, &gen_ptr[i], spots_used, MPI::FLOAT, MPI_SUM, 0);
				}
#endif // _INDEL
#endif // DISCRETIZE/INTDISC/NONE

			}
		}
		
		fprintf(stderr,"\n[-/%d] Finished!  Printing final .sgr/.gmp file\n",iproc);
		// Then print out our new genome
		if ((gSAM2GMP) && (iproc == 0 )) //added by CJ to select an option of creating gmp file, 04/23/2012
        {
			gGen.PrintFinal(output_file);
        }

#ifdef DISCRETIZE
		MPImyOp.Free();
		MPIctype.Free();
#endif
			
		free(sendGenome);
	}
	else
		// Even if we're doing MPI, have everyone print to their (separate) output files
		if (gSAM2GMP) //added by CJ to select an option of creating gmp file, 04/23/2012
        {
			gGen.PrintFinal(output_file);
        }
		
#else //MPI_RUN
	if (gSAM2GMP) //added by CJ to select an option of creating gmp file, 04/23/2012
    {
		gGen.PrintFinal(output_file);
    }
#endif //MPI_RUN

#ifdef MPI_RUN
	//fprintf(stderr,"[-/%d] Reducing seqs_matched from [%u] and not matched from [%u]...\n",iproc,seqs_matched,seqs_not_matched);
	unsigned int bseqs_matched = seqs_matched;
	unsigned int bseqs_not_matched = seqs_not_matched;
	if(gMPI_LARGEMEM)
    {
		// This won't be accuracte, but we'll just get the best matches and return them
		MPI::COMM_WORLD.Reduce(&bseqs_matched, &seqs_matched, 1, MPI_UNSIGNED, MPI_MAX, 0);
		MPI::COMM_WORLD.Reduce(&bseqs_not_matched, &seqs_not_matched, 1, MPI_UNSIGNED, MPI_MIN, 0);
	}
	else
    {
		// Here, we DO want to get a sum over all the matches
		MPI::COMM_WORLD.Reduce(&bseqs_matched, &seqs_matched, 1, MPI_UNSIGNED, MPI_SUM, 0);
		MPI::COMM_WORLD.Reduce(&bseqs_not_matched, &seqs_not_matched, 1, MPI_UNSIGNED, MPI_SUM, 0);
	}
	//fprintf(stderr,"[-/%d] After reduction, seqs_matched=[%u] and not matched=[%u]\n",iproc,seqs_matched,seqs_not_matched);
	
#endif //MPI_RUN

	if(iproc == 0 )
    {
		double time_prog_end = When();
		
		ostringstream stats;
		stats 	<< "\n#Finished.\n";
		stats	<< "#\tTotal Time: " << time_prog_end-prog_start_time << " seconds.\n";
#ifdef DEBUG_TIME
		stats   << "#\tTime spent waiting: " << cond_lock_time << " + " << write_lock_time << "=" 
											<< cond_lock_time+write_lock_time << " seconds.\n";
#endif //DEBUG_TIME
		stats	<< "#\tFound " << seqs_matched+seqs_not_matched << " sequences.\n";
		stats	<< "#\tSequences matched: " << seqs_matched << "\n";
		stats	<< "#\tSequences not matched: " << seqs_not_matched << "\n";
		stats	<< "#\tOutput written to " << sam_file << "\n";
		
		cerr << stats.str();
		// Don't print stats to output file
		//of << stats.str();		
	}

	of.close();
#ifdef DEBUG
	fprintf(stderr, "[%d/-] Freeing memory...\n", iproc);
#endif //DEBUG
	if(!gMPI_LARGEMEM)
    {
		delete gSM;
    }
	delete[] finished_arr;
	delete[] sam_file;

	if(gpu) {
		delete cudaDriver;
	}

	// Deallocate space requested for these arrays
#ifdef DEBUG
	fprintf(stderr, "[%d/-] Cleaning reads...\n", iproc);
#endif
	clean_nReads();
	
#ifdef DEBUG		
	fprintf(stderr,"[%d/-] End of program...waiting...\n",iproc);
#endif
#ifdef MPI_RUN
	MPI::Finalize();
#endif

	return 0;
}


void comm_cond_wait() {
	// We only want to wait in here if we need to communicate anything.  Otherwise, just continue
	if(!gMPI_LARGEMEM)
    {
		return;
    }
		
#ifdef MPI_RUN
#ifdef DEBUG_TIME
	double begin_time = When();
#endif

	MUTEX_LOCK(&cond_lock);

#ifdef DEBUG_TIME
	double end_time = When();
	cond_lock_time += end_time-begin_time;
#endif

	if(++cond_count < cond_thread_num)
    {
		//fprintf(stderr,"Thread %d/%d waiting...\n",cond_count,cond_thread_num);
		
#ifdef DEBUG_TIME		
		begin_time = When();
#endif

		pthread_cond_wait(&comm_barrier_cond, &cond_lock);

#ifdef DEBUG_TIME
		end_time = When();
		cond_lock_time += end_time-begin_time;
#endif
		
		MUTEX_UNLOCK(&cond_lock);	
	}
	else
    {
		//fprintf(stderr,"[%d/%d] Thread %d/%d made it in!\n",iproc,nproc,cond_count,cond_thread_num);
		// Do your stuff

		// MPI communicate everything you need to		
#ifdef DEBUG_TIME		
		double start_time = When();
		// set the global variable to be my start time
		gTimeBegin = start_time;
#endif
		
		// Do a reduce-all with a SUM on the denominators
		double sendDenominator[ nReads ];
		memcpy( sendDenominator, gReadDenominator, nReads * sizeof( double ) );

		MPI::COMM_WORLD.Allreduce(sendDenominator, gReadDenominator, nReads, MPI::DOUBLE, MPI_SUM);

#ifdef DEBUG_TIME
		double end_time = When();
#endif

#ifdef DEBUG_NW
		unsigned int total_nw = 0;
		unsigned int max_nw = 0;
		unsigned int min_nw = num_nw[ 0 ];

		for( int i = 0; i < nproc; i++ )
        {
			total_nw += num_nw[ i ];

			if( num_nw[ i ] > max_nw )
            {
				max_nw = num_nw[ i ];
            }
			if( num_nw[ i ] < min_nw )
            {
				min_nw = num_nw[ i ];
            }
		}

		fprintf( stderr,"[%d] Total time for Denoms is %f.  Total NW: %u, Max NW: %u, Min NW: %u, DIFF: %u\n", iproc, end_time - start_time, total_nw, max_nw, min_nw, max_nw - min_nw );
#endif		
#ifdef DEBUG_TIME
		start_time = When();
#endif
		
		// Do a reduce-all with a MAX on the top reads
		double sendReads[nReads];
		memcpy(sendReads, gTopReadScore, nReads*sizeof(double));

		MPI::COMM_WORLD.Allreduce( sendReads, gTopReadScore, nReads, MPI::DOUBLE, MPI_MAX );
#ifdef DEBUG_TIME
		end_time = When();
		fprintf( stderr, "[%d] Total time for TopReads is %f\n", iproc, end_time - start_time );
#endif

		cond_count = 0;
		MUTEX_UNLOCK( &cond_lock );
		pthread_cond_broadcast( &comm_barrier_cond );
	}
#endif
}

void write_cond_wait() {
#ifdef DEBUG_TIME
	double begin_time = When();
#endif
	MUTEX_LOCK(&cond_lock);
#ifdef DEBUG_TIME
	double end_time = When();
	cond_lock_time += end_time-begin_time;
#endif
	if(++cond_count < cond_thread_num)
    {
		//fprintf(stderr,"Thread %d/%d waiting...\n",cond_count,cond_thread_num);
#ifdef DEBUG_TIME
		begin_time = When();
#endif
		pthread_cond_wait(&comm_barrier_cond, &cond_lock);
#ifdef DEBUG_TIME
		end_time = When();
		cond_lock_time += end_time-begin_time;
#endif
		
		MUTEX_UNLOCK(&cond_lock);
	}
	else
    {
		//fprintf(stderr,"Thread %d/%d made it in!\n",cond_count,cond_thread_num);

#ifdef MPI_RUN
		if(gMPI_LARGEMEM) {

			//fprintf(stderr,"\nSize of vector: %u\n",gTopReadOutput.size());
			unsigned int total_outputed = 0;
			
			// Each machine prints out their own reads
			
#ifdef DEBUG_TIME
			begin_time = When();
#endif
			// Lock on the vector so others don't use it
			MUTEX_LOCK(&write_lock);
#ifdef DEBUG_TIME
			end_time = When();
			write_lock_time += end_time-begin_time;
#endif

            // print out all of the reads to a sam file
			vector<TopReadOutput>::iterator vit;
			for(vit = gTopReadOutput.begin(); vit != gTopReadOutput.end(); vit++) {
				if((*vit).MAPQ >= 0) {
					total_outputed++;
					
					of << (*vit).READ_NAME << "\t";
					if((*vit).strand == POS_STRAND)
						of << 0x0000 << "\t";
					else 
						of << 0x0010 << "\t";
					of << (*vit).CHR_NAME << "\t" << (*vit).CHR_POS << "\t";
					of << (*vit).MAPQ << "\t";
					if((*vit).strand == POS_STRAND)
						of << (*vit).CIGAR << "\t";
					else
						of << reverse_CIGAR((*vit).CIGAR) << "\t";
					of << "*\t0\t0\t";
					if((*vit).strand == POS_STRAND) {
						of << (*vit).consensus << "\t";
						of << (*vit).qual << "\t";	
					}
					else {	// need to do the reverse compliment of the sequence and qual
						of << reverse_comp((*vit).consensus) << "\t";
						of << reverse_qual((*vit).qual) << "\t";
					}
					// We're adjusting the alignment score internally.  The math works, but it's not
					// intuitive unless we adjust the score back when we print it out
					of << "XA:f:" << (float)((*vit).A_SCORE)*(1.0/gADJUST) << "\t";
					of << "XP:f:" << (float)((*vit).POST_PROB) << "\t";
					of << "X0:i:" << (*vit).SIM_MATCHES << "\n";
				}
				else { // It didn't really match
					of << (*vit).READ_NAME << "\t";
					of << 0x0200 << "\t*\t0\t0\t*\t=\t0\t0\t";
					of << (*vit).consensus << "\t";
					of << (*vit).qual << "\n";
				}
			}

			// Clear the output vector
			gTopReadOutput.clear();
			MUTEX_UNLOCK(&write_lock);
				
		
		} // endif(gMPI_LARGEMEM)
		else {
			// Write it all to the file 
#ifdef DEBUG_TIME
			begin_time = When();
#endif
			MUTEX_LOCK(&write_lock);
#ifdef DEBUG_TIME
			end_time = When();
			write_lock_time += end_time-begin_time;
#endif
			vector<TopReadOutput>::iterator vit;
			for(vit = gTopReadOutput.begin(); vit != gTopReadOutput.end(); vit++) {
				//string consensus = (*vit).consensus;
				//string QUAL = (*vit).qual;
				
				of << (*vit).READ_NAME << "\t";
				if((*vit).strand == POS_STRAND)
					of << 0x0000 << "\t";
				else
					of << 0x0010 << "\t";
				of << (*vit).CHR_NAME << "\t" << (*vit).CHR_POS << "\t";
				of << (*vit).MAPQ << "\t";
				if((*vit).strand == POS_STRAND)
					of << (*vit).CIGAR << "\t";
				else
					of << reverse_CIGAR((*vit).CIGAR) << "\t";
				of << "*\t0\t0\t";
				//of << consensus << "\t";
				//of << QUAL.c_str() << "\t";
				if((*vit).strand == POS_STRAND) {
					of << (*vit).consensus << "\t";
					of << (*vit).qual << "\t";	
					//of << consensus << "\t";
					//of << QUAL.c_str() << "\t";	
				}
				else {	// need to do the reverse compliment of the sequence and qual
					of << reverse_comp((*vit).consensus) << "\t";
					of << reverse_qual((*vit).qual) << "\t";
					//of << reverse_comp(consensus) << "\t";
					//of << reverse_qual(QUAL) << "\t";
				}
				// We're adjusting the alignment score internally.  The math works, but it's not
				// intuitive unless we adjust the score back when we print it out
				of << "XA:f:" << (float)((*vit).A_SCORE)*(1.0/gADJUST) << "\t";
				of << "XP:f:" << (float)((*vit).POST_PROB) << "\t";
				of << "X0:i:" << (*vit).SIM_MATCHES << "\n";
			}
			// Clear the output vector
			gTopReadOutput.clear();
			MUTEX_UNLOCK(&write_lock);
		}
#endif
		
		cond_count = 0;
		MUTEX_UNLOCK(&cond_lock);
		pthread_cond_broadcast(&comm_barrier_cond);
	}// endelse(++cond_count < cond_thread_num)
}


void single_write_cond_wait(int thread_no) {
	
	//fprintf(stderr,"[%d] Thread %d made it in!\n",thread_no,cond_count);
#ifdef DEBUG_TIME
	double wait_begin = When();
#endif
	MUTEX_LOCK(&write_lock);
#ifdef DEBUG_TIME
	double wait_end = When();
	write_lock_time += wait_end - wait_begin;
#endif
	// Write it all to the file 
	vector<TopReadOutput>::iterator vit;
	for(vit = gTopReadOutput.begin(); vit != gTopReadOutput.end(); vit++) {
		// Means it's a valid match
		if((*vit).MAPQ >= 0) {

			of << (*vit).READ_NAME << "\t";
			if((*vit).strand == POS_STRAND)
				of << 0x0000 << "\t";
			else
				of << 0x0010 << "\t";
			of << (*vit).CHR_NAME << "\t" << (*vit).CHR_POS << "\t";
			of << (*vit).MAPQ << "\t";
			if((*vit).strand == POS_STRAND)
				of << (*vit).CIGAR << "\t";
			else
				of << reverse_CIGAR((*vit).CIGAR) << "\t";
			of << "*\t0\t0\t";
			if((*vit).strand == POS_STRAND) {
				of << (*vit).consensus << "\t";
				of << (*vit).qual << "\t";	
			}
			else {	// need to do the reverse compliment of the sequence and qual
				of << reverse_comp((*vit).consensus) << "\t";
				of << reverse_qual((*vit).qual) << "\t";
			}
			//of << (*vit).consensus << "\t";
			//of << (*vit).qual << "\t";
			// We're adjusting the alignment score internally.  The math works, but it's not
			// intuitive unless we adjust the score back when we print it out
			of << "XA:f:" << (float)((*vit).A_SCORE)*(1.0/gADJUST) << "\t";
			of << "XP:f:" << (float)((*vit).POST_PROB) << "\t";
			of << "X0:i:" << (*vit).SIM_MATCHES << "\n";
		}
		else { // It didn't really match
			of << (*vit).READ_NAME << "\t";
			of << 0x0200 << "\t*\t0\t0\t*\t=\t0\t0\t";
			of << (*vit).consensus << "\t";
			of << (*vit).qual << "\n";
		}
	}
	// Clear the output vector
	gTopReadOutput.clear();
	
	MUTEX_UNLOCK(&write_lock);
}


// We just want to make sure we're clean before we delete all the reads...
void clean_cond_wait(bool my_finished) {
#ifdef DEBUG_TIME
	double begin_time = When();
#endif
	MUTEX_LOCK(&cond_lock);
#ifdef DEBUG_TIME
	double end_time = When();
	cond_lock_time += end_time-begin_time;
#endif

	finished += my_finished;

	if(++cond_count < cond_thread_num) {
#ifdef DEBUG_TIME
		begin_time = When();
#endif
		pthread_cond_wait(&comm_barrier_cond, &cond_lock);
#ifdef DEBUG_TIME
		end_time = When();
		cond_lock_time += end_time-begin_time;
#endif
		
		MUTEX_UNLOCK(&cond_lock);
		
	}
	else {
		//fprintf(stderr, "Exiting the output thingy\n" );
		// Now that we're done, reset the global SeqManager's edit point
		gSM->resetCounter();
		gReadsDone = 0;
		iter_num++;
		
		cond_count = 0;
		pthread_cond_broadcast(&comm_barrier_cond);
		MUTEX_UNLOCK(&cond_lock);
	}
}

void single_clean_cond_wait(bool my_finished, int thread_no, bool verbose) {
	
	if(!my_finished && thread_no)
		return;

#ifdef DEBUG_TIME
	double begin_time = When();
#endif
	MUTEX_LOCK(&cond_lock);
#ifdef DEBUG_TIME
	double end_time = When();
	cond_lock_time += end_time-begin_time;
#endif
	
	finished_arr[thread_no] = my_finished;

	// Only have one thread check
	if(thread_no == 0 ) {
		unsigned int total = 0;
		// We want to make sure every thread has checked in here to verify that they're done.
		for(unsigned int i=0; i<gNUM_THREADS; i++) {
			total += finished_arr[i];
		}
		if(total == gNUM_THREADS)
        {
			finished = true;
        }
			
		iter_num++;
		if(verbose)
        {
			gSM->resetCounter();
        }
			
	}
	MUTEX_UNLOCK(&cond_lock);
}

/**
 * Thread runner for non-MPI required runs
 * This is NOT MPI_largemem.
 */
void* parallel_thread_run(void* t_opts) {
	//cout << "New Thread\n";
	//fprintf(stderr, "CALLING PARALLEL_THREAD_RUN\n" );
    
	unsigned int thread_id = ((thread_opts*)t_opts)->thread_id;

    //cout << "Calling thread no: " << thread_id << endl;

	unsigned int numRPT = READS_PER_PROC;

	unsigned int read_begin = numRPT*thread_id;
	unsigned int read_end = read_begin + numRPT;

	unsigned int good_seqs = 0;
	unsigned int bad_seqs = 0;

	unsigned int seqs_processed = 0;
	
	bool my_finished = false;

    // print the header lines
    if( thread_id == 0 )
    {
        vector< pair< string,unsigned long> >::iterator header_it;

        for( header_it = gGen.GetNames()->begin(); header_it != gGen.GetNames()->end() - 1; header_it++ )
        {
            of << "@SQ\tSN:" << ( *header_it ).first << "\tLN:" << ( ( * ( header_it + 1 ) ).second - ( *header_it ).second ) << endl;
        }
        
        of << "@PG\tID:gnumap\tPN:gnumap\tVN:" << gVERSION << "\tCL:" << cl_args << endl;
    }

	// Fill the read buffer
	while( !my_finished )
    {

		//fprintf(stdout,"[%d/%d at %d] Thread continuing, finshed=%d, mine=%d\n",
		//		iproc,thread_id,iter_num,finished,my_finished);

		//fprintf(stderr,"[%d/%d] Getting more reads from %u to %u\n",thread_id,iproc,read_begin,read_end);
		my_finished = !gSM->getMoreReads( read_begin, read_end, false );
		//fprintf(stderr,"[%d/%d] Got more reads successfully? %s from %d to %d\n",
		//		thread_id,iproc,!my_finished ? "Y" : "**N**",read_begin,read_end);


		for( unsigned int k = read_begin; k < read_end; k++ )
        {
			//vector<vector<double> > temp_vect = ((thread_opts*)t_opts)->seqs[k];
			Read* temp_read = gReadArray[k];
			if( !temp_read )
            {
				break;
            }
			
            //cerr << thread_id << ":\t" << temp_read->name << endl;
			string consensus = GetConsensus( *temp_read );
			set_top_matches( ( ( thread_opts* ) t_opts )->gen, k, consensus, bad_seqs, good_seqs, thread_id );
			seqs_processed++;
			
		}  //end for loop
		
		// We'll break these two loops up so they can be doing different work instead of
		// getting caught up on all the mutex locks.
		for( unsigned int k = read_begin; k < read_end; k++ )
        {
			if( !gReadArray[ k ] )
            {
				break;
            }
			
			string consensus = GetConsensus( *( gReadArray[ k ] ) );
			create_match_output( ( ( thread_opts* )t_opts )->gen, k, consensus );

			delete_read( gReadArray[ k ] );
			gReadArray[ k ] = 0;
		}
		
		// This just sends one processor to print everything out
		if(thread_id == 0 )
        {
			single_write_cond_wait(thread_id);
        }
			
		single_clean_cond_wait(my_finished, thread_id, true);
	} // end while(!my_finished)
		
	// Make sure we wait for everyone to finish here
	if(thread_id == 0 )
    {
		while(!finished)
        {
			single_clean_cond_wait(my_finished, thread_id, false);
        }
    }

	//fprintf(stderr,"[%d/%d] getting ready to enter single_write_cond_wait...\n",iproc,thread_id);
	// Print them out here, just to make sure
	if(thread_id == 0 )
    {
		single_write_cond_wait(thread_id);
    }

	//fprintf(stderr,"[%d/%d] Finished after processing %d/%u reads with %u good and %u bad\n",iproc,thread_id,seqs_processed,gSM->getTotalSeqCount(),good_seqs,bad_seqs);

	struct thread_rets* ret = new thread_rets;
	ret->good_seqs = good_seqs;
	ret->bad_seqs = bad_seqs;
	
	return (void*)ret;
}

/**
 * This is just a block of code, not really a function call.  We just want to put everything in here...
 * 
 * By including the flag -DOMP_RUN, this function will use openmp to split up the reads.  Just for testing
 * to see if I could do better (and I can)
 */
#ifdef OMP_RUN
thread_rets* omp_thread_run(thread_opts* t_opts) {
	fprintf(stderr, "CALLING OMP_THREAD_RUN\n" );

	unsigned int i,j, nReadsRead, good_seqs=0, bad_seqs=0;
#ifdef DEBUG_TIME
	double total_time=0, comm_time=0;
#endif

#pragma omp parallel shared(t_opts,gReadArray,gReadDenominator,gTopReadScore)
{
	while(!finished) {

#ifdef DEBUG_NW
		fprintf(stderr,"[%d/%d] Num nw's: %u\n",iproc,omp_get_thread_num(),num_nw[omp_get_thread_num()]);
		num_nw[omp_get_thread_num()] = 0;
#endif

#pragma omp single
{
		t_opts->sm.fillReadArray(gReadArray,nReads,nReadsRead);
		fprintf(stderr,"[-/%d] Num Reads to read: %u and found %u.  Finished is %d\n",iproc,nReads,nReadsRead,finished);

}

		// Use omp to parallelize this code
#pragma omp for private(i) schedule(dynamic, 128 ) reduction(+:good_seqs) reduction(+:bad_seqs)
		for(i=0; i<nReadsRead; i++) {
			//fprintf(stderr,"[%d:%d] Setting match output for %d\n",iproc,omp_get_thread_num(),i);
			string consensus = GetConsensus(*(gReadArray[i]));
			unsigned int thread_no = omp_get_thread_num();
			set_top_matches(t_opts->gen, i, consensus, bad_seqs, good_seqs, thread_no);
		}

// wait for all the threads to sync
#pragma omp single
{
#ifdef DEBUG_TIME
		double start_comm = When();
#endif
		comm_cond_wait();
#ifdef DEBUG_TIME
		comm_time += When()-start_comm;
		fprintf(stderr,"[%d] Made it past the comm_cond_wait: %f\n",iproc,When()-start_comm);
#endif
}

#pragma omp for private(i) schedule(dynamic, 128 )
		for(i=0; i<nReadsRead; i++) {
			//fprintf(stderr,"[%d:%d] Creating match output for %d\n",iproc,omp_get_thread_num(),i);
			string consensus = GetConsensus(*(gReadArray[i]));
			create_match_output(t_opts->gen, i, consensus);
		}

#pragma omp single
{
#ifdef DEBUG_TIME
		double start = When();
#endif
		write_cond_wait();
#ifdef DEBUG_TIME
		total_time += When()-start;
		fprintf(stderr,"[%d] Made it past the write_cond_wait: %f\n",iproc,When()-start);
#endif
}

#pragma omp for private(i,j) schedule(static)
		for(i=0; i<nReadsRead; i++) {
			if(!gReadArray[i])
				continue;

			for(j=0; j<gReadArray[i]->length; j++)
				delete[] gReadArray[i]->pwm[j];

			delete[] gReadArray[i]->pwm;
			if(gReadArray[i]->name)
				delete[] gReadArray[i]->name;

			delete gReadArray[i];
			gReadArray[i] = 0;
			gReadDenominator[i] = 0;
			gTopReadScore[i] = 0;
		}

#pragma omp single
{
		gSM->printProgress();
#ifdef DEBUG_TIME
		fprintf(stderr,"[-/%d] Seconds used in write_comm: %f\tcomm_cond: %f\n",iproc,total_time,comm_time);
#endif
		finished = gSM->isFinished();
}


	}
}
	
#ifdef DEBUG_TIME
	fprintf(stderr,"[-/%d] Total time used in write_comm: %f\tcomm_cond: %f\n",iproc,total_time,comm_time);
#endif
	
	struct thread_rets* ret = new thread_rets;
	ret->good_seqs = good_seqs;
	ret->bad_seqs = bad_seqs;
	return ret;
}

#endif

/**
 * This function will be used by each thread when the MPI_largmem flag is used.
 */
void* mpi_thread_run(void* t_opts) {
	fprintf(stderr, "CALLING MPI_THREAD_RUN\n" );

	while(!setup_complete);	// Busy-Wait for the setup to complete
	
	unsigned int thread_id = ((thread_opts*)t_opts)->thread_id;
	unsigned int i,j,begin=0,end=0;
	unsigned int bad_seqs=0, good_seqs=0;
	
	bool my_finished = false;
		
	// main loop body
	while(!finished) {

#ifdef DEBUG_NW
		// For testing, reset the number of NW's we're doing
		num_nw[thread_id] = 0;
#endif
		 
		// We'll be looping until we're finished with everything.  
		//Just read a portion at a time
		MUTEX_LOCK(&total_read_lock);
		while(gReadsDone < nReads && !my_finished) {
		
			// Get your portion of the reads here
			// It's locking, so we don't need to worry about multi-entrances
			//my_finished = !((thread_opts*)t_opts)->sm.getMoreReads(begin, end, true);
			my_finished = !gSM->getMoreReads(begin, end, true);

			gReadsDone += end-begin;
			MUTEX_UNLOCK(&total_read_lock);


			//fprintf(stderr,"[%u/%u] Thread %u has seqs %u-%u\n",iproc,thread_id,thread_id,begin,end);
			// Do the first part of the alignment
			unsigned int count = 0;
			for(i=begin; i<end; i++, count++) {
	
				string consensus = GetConsensus(*(gReadArray[i]));
				set_top_matches(((thread_opts*)t_opts)->gen, i, consensus, bad_seqs, good_seqs, thread_id);
	
				// For some reason, the percentage complete doesn't work...
				//if(gVERBOSE)
				//	fprintf(stderr,"\r%d%% complete",(int)(((double)i)/gSM->getNumSeqsPProc())*100);
	
			}
			
			MUTEX_LOCK(&total_read_lock);
		} // End of reading portions
		MUTEX_UNLOCK(&total_read_lock);
		

		//fprintf(stderr,"[%d/%d] Thread %d about to enter cond_wait\n",iproc,nproc,((thread_opts*)t_opts)->thread_id);
		// Wait for each thread to finish
		comm_cond_wait();

#ifdef DEBUG_NW
		double b_time = When();
		double e_time = When();
		fprintf( stderr, "%u - [%d/%d] Total time for comm_wait: %f, since start: %f, difference: %f, num_nw: %u\n",
				iter_num, iproc, thread_id, e_time-b_time, gTimeBegin - b_time, e_time - gTimeBegin, num_nw[ thread_id ] );
		//fprintf(stderr,"Thread %d just exited cond_wait\n",((thread_opts*)t_opts)->thread_id);
#elif defined(DEBUG_TIME)
		fprintf(stderr,"%u - [%d/%d] Total time for comm_wait: %f, since start: %f, difference: %f\n",
				iter_num,iproc,thread_id, e_time-b_time, gTimeBegin-b_time, e_time-gTimeBegin);
#endif

		// Now just set the begin and end to be an equivalent portion of reads
		unsigned int readsPerProc = nReads/gNUM_THREADS;
		begin = readsPerProc*thread_id;
		end = readsPerProc*(thread_id+1);
		end = (end > nReads) ? nReads : end;
		
		for(i=begin; i<end; i++) {
				
			// We'll just continue if it's empty
			if(!gReadArray[i])
				continue;

			string consensus = GetConsensus(*(gReadArray[i]));
			create_match_output(((thread_opts*)t_opts)->gen, i, consensus);
		}

		//fprintf(stderr,"Thread %d about to enter write cond_wait\n",((thread_opts*)t_opts)->thread_id);
		// Wait for each thread to finish, then write them.
#ifdef DEBUG_TIME
		b_time = When();
#endif
		write_cond_wait();
#ifdef DEBUG_TIME
		e_time = When();
		fprintf(stderr,"[%d/%d] Total time for write_wait: %f\n", iproc,thread_id, e_time-b_time);
#endif
		
		// Delete the reads you've used
		for(i=begin; i<end; i++) {
			if(!gReadArray[i])
				continue;
			
			for(j=0; j<gReadArray[i]->length; j++) {
				delete[] gReadArray[i]->pwm[j];
			}
			delete[] gReadArray[i]->pwm;
			if(gReadArray[i]->name)
				delete[] gReadArray[i]->name;

			delete gReadArray[i];
			gReadArray[i] = 0;
			gReadDenominator[i] = 0;
			gTopReadScore[i] = 0;
			
			//clearMapAt(i);
		}
		
		// Even if we're not using MPI, we should wait here.
		clean_cond_wait(my_finished);
		//fprintf(stderr,"Thread %d finished, starting over\n",((thread_opts*)t_opts)->thread_id);
	}
	
	struct thread_rets* ret = new thread_rets;
	ret->good_seqs = good_seqs;
	ret->bad_seqs = bad_seqs;
	
	return (void*)ret;
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
				//printf("Found extended: %s\n",argv[i]);
				set_ret = set_arg_ext(argv[i],argv[i+1],i);
			}
			else if(strlen(argv[i]) > 2) {
				set_ret = PARSE_ERROR;
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
				seq_file = argv[i];
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
		os << "Please specify only a single sequence file, or multiple files in \"quotations\" ("
		   << argv[this_cmd.errpos] << ")" << endl;
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
		case 'l':
			if(sscanf(assign,"%d",&gSEQ_LENGTH) < 1)
				return PARSE_ERROR;
			break;
		case 'o':
			output_file = assign;
			break;
		case 'a':
			if(sscanf(assign,"%lf",&temp_dbl)< 1)
				return PARSE_ERROR;
			gALIGN_SCORE = (float)temp_dbl;
			break;
		case 'p':
			count--;
			perc = true; // but it's already true.  Just ignore it
			break;
		case 'r':
			count--;
			perc = false;
			break;
		case 'q':
			if(sscanf(assign,"%lf",&temp_dbl) < 1)
				return PARSE_ERROR;
			gCUTOFF_SCORE = (float)temp_dbl;
			break;
		case 'v':
			if(sscanf(assign,"%d",&gVERBOSE) < 1)
				return PARSE_ERROR;
			break;
		case 'c':
			if(sscanf(assign,"%d",&gNUM_THREADS) < 1)
				return PARSE_ERROR;
			break;
		case 'm': 
			if(sscanf(assign,"%d",&gMER_SIZE) < 1)
				return PARSE_ERROR;
			break;
		case 'B':
			if(sscanf(assign,"%d",&gBUFFER_SIZE) < 1)
				return PARSE_ERROR;
			break;
		case 'u':
			gUNIQUE = true;
			break;
		case 'T':
			if(sscanf(assign,"%d",&gMAX_MATCHES) < 1)
				return PARSE_ERROR;
			break;
		case 'h':
			if(sscanf(assign,"%d",&gMAX_KMER_SIZE) < 1)
				return PARSE_ERROR;
			break;
		case 'S':
			pos_matrix = assign;
			break;
		case 'G':
			//if(sscanf(assign,"%lf",&gGAP) < 1)
			if(sscanf(assign,"%lf",&temp_dbl) < 1)
				return PARSE_ERROR;
			gGAP = (float)temp_dbl;
			break;
		case 'M':
			if(sscanf(assign,"%d",&gMAX_GAP) < 1)
				return PARSE_ERROR;
			break;
		case 'A':
			g_adaptor = new char[strlen(assign)+1];
			unsigned int i;
			for(i=0; i<strlen(assign); i++)
				g_adaptor[i] = (char)assign[i];
			// Null-terminate the string
			g_adaptor[i] = '\0';
			//g_adaptor = assign;
			break;
		case '0':
			count--;
			gPRINT_FULL = 1;
			break;
		case 'b':
			count--;
			gBISULFITE = 1;
			gGEN_SIZE = 1;
			gPRINT_FULL = 1;
			break;
		case 'd':
			count--;
			gATOG = 1;
			gGEN_SIZE = 1;
			gPRINT_FULL = 1;
			break;
		case 's':
			if(sscanf(assign,"%d",&gGEN_SKIP) < 1)
				return PARSE_ERROR;
			break;
		case 'j':
			if(sscanf(assign,"%d",&gJUMP_SIZE) < 1)
				return PARSE_ERROR;
			break;
		case 'k':
			if(sscanf(assign,"%d",&gMIN_JUMP_MATCHES) < 1)
				return PARSE_ERROR;
			break;
		case '?':
			usage(0, "" );
		default:
			return PARSE_ERROR;
	}

	return 0;

}

enum {
	PAR_ORIG=256,
	PAR_GENOME,
	PAR_LENGTH,
	PAR_OUTPUT,
	PAR_ALIGN_SCORE,
	PAR_PERCENT, // include backward compatability
	PAR_RAW,
	PAR_CUTOFF_SCORE,
	PAR_VERBOSE,
	PAR_NUM_PROC,
	PAR_MER_SIZE,
	PAR_BUFFER,
	PAR_MAX_MATCH,
	PAR_UNIQUE,
	PAR_MAX_KMER,
	PAR_SUBST_FILE,
	PAR_GAP_PENALTY,
	PAR_MAX_GAP,
	PAR_ADAPTOR,
	PAR_PRINT_FULL,
	PAR_PRINT_ALL_SAM,
	PAR_VCF,
	PAR_BS_SEQ,
	PAR_A_TO_G,
	PAR_FAST,
	PAR_GEN_SKIP,
	PAR_READ,
	PAR_SAVE,
	PAR_GEN_SIZE,
	PAR_JUMP,
	PAR_N_MATCH,
	PAR_SNP,
	PAR_ILL,
	PAR_LARGEMEM,
	PAR_UPSTRAND,
	PAR_DOWNSTRAND,
	PAR_SNP_PVAL,
	PAR_SNP_MONOP,
	PAR_NO_GMP,
	PAR_ASSEMBLER,
	PAR_NO_NW,
	PAR_GPU
};

int set_arg_ext(const char* param, const char* assign, int &count) {
	double temp_dbl;
	int which = 0;
	int adjust=0;

    // TODO: fix the 'adjust' parameter as it appears to be inconsistent
	if( strncmp(param+2, "genome=", 7 ) == 0 ) {	//check to see if it's the genome
		//which = 'g';
		which = PAR_GENOME;
		adjust = 8;	//the size of "genome" plus 2
	}
	else if( strncmp(param+2, "length=", 7 ) == 0 ) {
		//which = 'l';
		which = PAR_LENGTH;
		adjust = 8;
	}
	else if( strncmp(param+2, "output=", 7 ) == 0 ) {
		//which = 'o';
		which = PAR_OUTPUT;
		adjust = 8;
	}
	else if( strncmp( param + 2, "align_score=", 12 ) == 0 ) {
		//which = 'a';
		which = PAR_ALIGN_SCORE;
		adjust = 13;
	}
	// only use strcmp when there won't be a second argument
	else if( strcmp( param + 2, "percent") == 0 ) {
		//which = 'p';
		which = PAR_PERCENT;
		adjust = 9;
	}
	else if( strcmp( param + 2, "no_nw" ) == 0 )
	{
		which = PAR_NO_NW;
		adjust = 7;
	}
	else if( strcmp( param + 2, "gpu" ) == 0 ) {
		which = PAR_GPU;
		adjust = 5;
	}
	else if( strcmp( param + 2, "raw") == 0 ) {
		//which = 'p';
		which = PAR_RAW;
		adjust = 9;
	}
	else if( strncmp( param + 2, "read_quality=", 13 ) == 0 ) {
		which = PAR_CUTOFF_SCORE;
		adjust = 15;
	}
	else if( strncmp( param + 2, "verbose=", 8 ) == 0 ) {
		//which = 'v';
		which = PAR_VERBOSE;
		adjust = 9;
	}
	else if( strncmp( param + 2, "num_proc=", 9 ) == 0 ) {
		//which = 'c';
		which = PAR_NUM_PROC;
		adjust = 10;
	}
	else if( strncmp( param + 2, "mer_size=", 9 ) == 0 ) {
		//which = 'm';
		which = PAR_MER_SIZE;
		adjust = 10;
	}
	else if( strncmp( param + 2, "buffer=", 7 ) == 0 ) {
    cerr << "gsm done!" << endl;
		//which = 'B';
		which = PAR_BUFFER;
		adjust = 8;
	}
	else if( strncmp( param + 2, "max_match=", 9 ) == 0 ) {
		//which = 'T';
		which = PAR_MAX_MATCH;
		adjust = 10;
	}
	else if( strncmp( param + 2, "subst_file=", 9 ) == 0 ) {
		//which = 'S';
		which = PAR_SUBST_FILE;
		adjust = 10;
	}
	else if( strncmp( param + 2, "gap_penalty=", 12 ) == 0 ) {
		//which = 'G';
		which = PAR_GAP_PENALTY;
		adjust = 13;
	}
	else if( strncmp( param + 2, "max_gap=", 8 ) == 0 ) {
		//which = 'M';
		which = PAR_MAX_GAP;
		adjust = 9;
	}
	else if( strncmp( param + 2, "adaptor=", 8 ) == 0 ) {
		//which = 'A';
		which = PAR_ADAPTOR;
		adjust = 9;
	}
	else if( strcmp( param + 2, "print_full" ) == 0 ) {
		//which = '0';
		which = PAR_PRINT_FULL;
		adjust = 12;
	}
	else if( strcmp( param + 2, "print_all_sam" ) == 0 ) {
		which = PAR_PRINT_ALL_SAM;
		adjust = 13;
	}
	else if( strcmp( param + 2, "vcf" ) == 0 ) {
		which = PAR_VCF;
		adjust = 3;
	}
	else if( strcmp( param + 2, "bs_seq" ) == 0 ) {
		//which = 'b';
		which = PAR_BS_SEQ;
		adjust = 8;
	}
	else if( strcmp( param + 2, "b2" ) == 0 ) {
		gBISULFITE2 = 1;
		which = PAR_BS_SEQ;
		adjust = 4;
	}
	else if( strcmp( param + 2, "a_to_g") == 0 ) {
		//which = 'd';
		which = PAR_A_TO_G;
		adjust = 8;
	}
	else if( strcmp( param + 2, "fast") == 0 ) {
		//which = 'f';
		which = PAR_FAST;
		adjust = 6;
	}
	else if( strncmp( param + 2,"gen_skip=", 9 ) == 0 ) {
		//which = 's';
		which = PAR_GEN_SKIP;
		adjust = 10;
	}
	else if( strncmp( param + 2,"bin_size=", 9 ) == 0 ) {
		//which = 3;
		which = PAR_GEN_SIZE;
		adjust = 10;
	}
	else if( strncmp( param + 2,"jump=", 5 ) == 0 ) {
		which = PAR_JUMP;
		adjust = 6;
	}
	else if( strncmp( param + 2,"num_seed=", 10 ) == 0 ) {
		which = PAR_N_MATCH;
		adjust = 11;
	}
	else if( strcmp( param + 2, "snp" ) == 0 ) {
		which = PAR_SNP;
		adjust = 5;
	}
	else if( strncmp( param + 2,"snp_pval=", 9 ) == 0 ) {
		which = PAR_SNP_PVAL;
		adjust = 11;
	}
	else if( strcmp( param + 2, "snp_monop" ) == 0 ) {
		which = PAR_SNP_MONOP;
		adjust = 11;
	}
	else if( strcmp( param + 2, "illumina" ) == 0 ) {
		which = PAR_ILL;
		adjust = 10;
	}
	else if( strcmp( param + 2, "MPI_largemem" ) == 0 ) {
		which = PAR_LARGEMEM;
		adjust = 14;
	}
	else if( strcmp( param + 2, "up_strand" ) == 0 ) {
		which = PAR_UPSTRAND;
		adjust = 11;
	}
	else if( strcmp( param + 2, "down_strand" ) == 0 ) {
		which = PAR_DOWNSTRAND;
		adjust = 13;
	}
	else if( strcmp( param + 2, "no_gmp" ) == 0 ) {
		which = PAR_NO_GMP;
		adjust = 8;
	}
	else if( strcmp( param + 2, "assembler" ) == 0 ) {
		which = PAR_ASSEMBLER;
		adjust = 11;
	}
	else if( strcmp( param + 2, "help" ) == 0 ) {
		usage(0, "" ); //usage will exit immediately
	}
	else 
		return NO_MATCHING_ARG;

	param += ++adjust;

	switch(which) {	//switch on the first character
		//case 'g':
		case PAR_GENOME:
			genome_file = param;
			break;
		//case 'l':
		case PAR_LENGTH:
			if(sscanf(param,"%d",&gSEQ_LENGTH) < 1)
				return PARSE_ERROR;
			break;
		//case 'o':
		case PAR_OUTPUT:
			output_file = param;
			break;
		//case 'a':
		case PAR_ALIGN_SCORE:
			//if(sscanf(param,"%lf",&ALIGN_SCORE) < 1)
			if(sscanf(param,"%lf",&temp_dbl) < 1)
				return PARSE_ERROR;
			//	return (PARSE_ERROR * (int)(*param));
			gALIGN_SCORE = (float)temp_dbl;
			break;
		//case 'p':
		case PAR_PERCENT:
			//count--;
			perc = true;
			break;
		//case 'r':
		case PAR_RAW:
			//count--;
			perc = false;
			break;
		// case 'q':
		case PAR_NO_NW:
			gNW = false;
			break;
		case PAR_GPU:
			gpu = true;
			break;
		case PAR_CUTOFF_SCORE:
			if(sscanf(param,"%lf",&temp_dbl) < 1)
				return PARSE_ERROR;
			gCUTOFF_SCORE = (float)temp_dbl;
			break;
		//case 'v':
		case PAR_VERBOSE:
			if(sscanf(param,"%d",&gVERBOSE) < 1)
				return PARSE_ERROR;
			break;
		//case 'c':
		case PAR_NUM_PROC:
			if(sscanf(param,"%d",&gNUM_THREADS) < 1)
				return PARSE_ERROR;
			break;
		//case 'm': 
		case PAR_MER_SIZE:
			if(sscanf(param,"%d",&gMER_SIZE) < 1)
				return PARSE_ERROR;
			break;
		//case 'B':
		case PAR_BUFFER:
			if(sscanf(param,"%d",&gBUFFER_SIZE) < 1)
				return PARSE_ERROR;
			break;
		//case 'T':
		case PAR_MAX_MATCH:
			if(sscanf(param,"%d",&gMAX_MATCHES) < 1)
				return PARSE_ERROR;
			break;
		//case 'h':
		case PAR_MAX_KMER:
			if(sscanf(param,"%d",&gMAX_KMER_SIZE) < 1)
				return PARSE_ERROR;
			break;
		//case 'S':
		case PAR_SUBST_FILE:
			pos_matrix = param;
			break;
		//case 'G':
		case PAR_GAP_PENALTY:
			//if(sscanf(param,"%lf",&gGAP) < 1)
			if(sscanf(param,"%lf",&temp_dbl) < 1)
				return PARSE_ERROR;
			gGAP = (float)temp_dbl;
			break;
		//case 'M':
		case PAR_MAX_GAP:
			if(sscanf(param,"%d",&gMAX_GAP) < 1)
				return PARSE_ERROR;
			break;
		//case 'A':
		case PAR_ADAPTOR:
			g_adaptor = new char[strlen(assign)];
			for(unsigned int i=0; i<strlen(assign); i++)
				g_adaptor[i] = tolower((char)param[i]);
			//g_adaptor = (char*)param;
			break;
		//case '0':
		case PAR_PRINT_FULL:
			//count--;	//don't decrement in the extended version
			gPRINT_FULL = 1;
			break;
		case PAR_PRINT_ALL_SAM:
			gPRINT_ALL_SAM = true;
			break;
		case PAR_VCF:
			gPRINT_VCF = true;
			break;
		//case 'f':
		case PAR_FAST:
			//count--;
			gFAST = true;
			break;
		//case 's':
		case PAR_GEN_SKIP:
			if(sscanf(param,"%d",&gGEN_SKIP) < 1)
				return PARSE_ERROR;
			break;
		case PAR_GEN_SIZE:
			if(sscanf(param,"%d",&gGEN_SIZE) < 1)
				return PARSE_ERROR;
			break;
		case PAR_JUMP: 
			if(sscanf(param,"%d",&gJUMP_SIZE) < 1)
				return PARSE_ERROR;
			break;
		case PAR_N_MATCH:
			if(sscanf(param,"%d",&gMIN_JUMP_MATCHES) < 1)
				return PARSE_ERROR;
			break;
		case PAR_SNP_PVAL:
			if(sscanf(param,"%f",&gSNP_PVAL) < 1)
				return PARSE_ERROR;
			break;
		case PAR_SNP_MONOP:
			gSNP_MONOP = 1;
			// let it fall through to PAR_SNP
			// We don't want to force them to specify both flags
		case PAR_SNP:
			gSNP = true;
			gGEN_SIZE = 1;
			gPRINT_FULL = 1;
			break;
		//case 'b':
		case PAR_BS_SEQ:
			//count--;
			gBISULFITE = true;
			gGEN_SIZE = 1;
			gPRINT_FULL = 1;
			break;
		//case 'd':
		case PAR_A_TO_G:
			//count--;
			gATOG = true;
			gGEN_SIZE = 1;
			gPRINT_FULL = 1;
			break;
		case PAR_ILL:
			gILLUMINA = true;
			break;
		case PAR_LARGEMEM:
			gMPI_LARGEMEM = true;
			break;
		// Don't match to the negative strand
		case PAR_UPSTRAND:
			gMATCH_NEG_STRAND = false;
			break;
		case PAR_DOWNSTRAND:
			gMATCH_POS_STRAND = false;
			break;
		case PAR_NO_GMP:
			gSAM2GMP = false;
			break;
		case PAR_ASSEMBLER:
			break;
		default:
			return NO_MATCHING_ARG;
	}

	return 0;

}
