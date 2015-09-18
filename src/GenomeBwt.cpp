#include "GenomeBwt.h"

/*
 * Constructor.  Creates new Genome with filename.
 * fn - the filename of the file to be parsed into the genome.
 */
GenomeBwt::GenomeBwt(const char* fn) : Genome(fn)
{
    std::cout << "constructor 1" << std::endl;
	Init(fn);
	reader = new Reader();
}

GenomeBwt::GenomeBwt() : Genome()
{
    std::cout << "constructor 2" << std::endl;
}

GenomeBwt::~GenomeBwt() {
	//fprintf(stderr,"[-/%d] Genome destructor...\n",iproc);
	//TODO: fix destructor for bwt/sa
} 

char * GenomeBwt::bwa_idx_infer_prefix(const char *hint)
{
    char *prefix;
    int l_hint;
    FILE *fp;
    l_hint = strlen(hint);
    //std::cout << "hint:\t" << hint << std::endl;
    //std::cout << "hint len:\t" << l_hint << std::endt;
    prefix = ( char* ) malloc(l_hint + 3 + 4 + 1 + 7 );
    strcpy(prefix, hint);
    //std::cout << "size:\t" << sizeof( prefix ) << std::endl;
    //std::cout << prefix << std::endl;
    strcpy(prefix + l_hint, ".64.gnumap.bwt");
    //std::cout << "new prefix:\t" << prefix << "\t" << strlen( prefix ) << std::endl;
    if ((fp = fopen(prefix, "rb")) != 0) {
        fclose(fp);
        prefix[l_hint + 3] = 0;
        return prefix;
    } else {
        strcpy(prefix + l_hint, ".gnumap.bwt");
        if ((fp = fopen(prefix, "rb")) == 0) {
            free(prefix);
            return 0;
        } else {
            fclose(fp);
            prefix[l_hint] = 0;
            return prefix;
        }
    }
}

bwt_t* GenomeBwt::bwa_idx_load_bwt(const char *hint)
{
    int bwa_verbose = 0;
    char *tmp, *prefix;
    bwt_t *bwt;
    prefix = bwa_idx_infer_prefix(hint);

    //std::cout << "inferred!\t" << prefix << std::endl;

    if (prefix == 0) {
        if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to locate the index files\n", __func__);
        return 0;
    }
    tmp = ( char* ) calloc(strlen(prefix) + 5 + 7, 1);
    strcat(strcpy(tmp, prefix), ".gnumap.bwt"); // FM-index
    bwt = bwt_restore_bwt(tmp);
    strcat(strcpy(tmp, prefix), ".gnumap.sa");  // partial suffix array (SA)
    bwt_restore_sa(tmp, bwt);
    free(tmp); free(prefix);
    return bwt;
}

bwaidx_t * GenomeBwt::bwa_idx_load_from_disk(const char *hint, int which)
{
    int bwa_verbose = 0;
    bwaidx_t *idx;
    char *prefix;
    prefix = bwa_idx_infer_prefix(hint);
    if (prefix == 0) {
        if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to locate the index files\n", __func__);
        return 0;
    }
    idx = ( bwaidx_t* ) calloc(1, sizeof(bwaidx_t));
    if (which & BWA_IDX_BWT) idx->bwt = bwa_idx_load_bwt(hint);
    if (which & BWA_IDX_BNS) {
    int i, c;
    idx->bns = bns_restore(prefix);
    for (i = c = 0; i < idx->bns->n_seqs; ++i)
        if (idx->bns->anns[i].is_alt) ++c;
            if (bwa_verbose >= 3)
                fprintf(stderr, "[M::%s] read %d ALT contigs\n", __func__, c);
        if (which & BWA_IDX_PAC) {
            idx->pac = ( uint8_t* ) calloc(idx->bns->l_pac/4+1, 1);
            err_fread_noeof(idx->pac, 1, idx->bns->l_pac/4+1, idx->bns->fp_pac); // concatenated 2-bit encoded sequence
            err_fclose(idx->bns->fp_pac);
            idx->bns->fp_pac = 0;
        }
    }
    free(prefix);
    return idx;
}

void GenomeBwt::Init(const char* fn) {

    ref_genome_fn = ( char* ) malloc( sizeof( char* ) * ( strlen( fn ) + 20 ) );
    strcpy( ref_genome_fn, fn );
    
	//clear the *_spots files
#if defined( DISCRETIZE )
	read_allot = 0;
#elif defined( INTDISC )
	reads = 0;
#else
	for(int i=0; i<NUM_SNP_VALS; i++) {
		reads[i] = 0;
	}
#endif // DISCRETIZE/INTDISC
}

void GenomeBwt::use(const char* fn) {
    std::cout << "genome::use" << std::endl;
	use(fn,0,0);
}

void GenomeBwt::use(const char* fn, unsigned long s, unsigned long e) {
	my_start = s;
	my_end = e;
	if(my_start >= my_end)
		// read everything
		read_all = true;
	else
		read_all = false;
	
	//fprintf(stderr,"Reading Genome from %lu to %lu, read_all=%d\n",s,e,read_all);
	
	if(reader != NULL) {
		cout << "non-null reader" << endl;
		delete reader;
	}
	Init(fn);
	reader = new Reader();
}

void GenomeBwt::print_to_file(hash_map<unsigned int, HashLocation> & gh, char* ofn) { 
	printf("Printing to file: %s\n",ofn);
	int SEQ_LENGTH = gMER_SIZE; //This is an archaic #def'd element that needed to be removed.
	ofstream out;
	out.open(ofn,ios_base::out);
	
	bin_seq bs;
	
	unsigned int i;
		
	for(i=0; i<gHASH_SIZE; ++i) {
		if(gh[i].size == 0)
			continue;
		//fprintf(stderr,"[%d] arry position %u (%s) has size %lu\n",
		//		iproc,i,bs.hash2str(i,gMER_SIZE).c_str(),gh[i].size);
		out << i << "\t" << gh[i].size;
		unsigned int j;
		out << '\t' << bs.hash2str(i,gMER_SIZE);
		for(j=0; j<gh[i].size; j++) {
			//fprintf(stderr,"\t[%d] gh[%u].hash_arr[%u]\n",iproc,i,j);
			out << '\t' << gh[i].hash_arr[j] << "[" << GetPos(gh[i].hash_arr[j],true) << "] :" 
				<< GetString(gh[i].hash_arr[j],SEQ_LENGTH);
		}
		out << endl;
		
	}

	out << "Finished\n" << endl;
	
	/*
	string gen_str = GetString(0,genome_size);
	cout << gen_str << endl;
	
	for(unsigned int i=0; i<names.size(); i++) {
		cout << names[i].first << "\t" << names[i].second << endl;
	}
	*/
	out.close();
}

unsigned int GenomeBwt::readGen(char* fn) {
    //std::cout << "readgen:\t" << fn << std::endl;
	index = bwa_idx_load_from_disk( fn, BWA_IDX_ALL );

    if( index == NULL )
    {
        //perror( "failed to load genome index!\n" );
        //throw "Could not load genome index";
        return 0;
    }
    
	return 1;
}

void GenomeBwt::make_extra_arrays() {
	//if we are doing the bisulfite sequencing, we want to record what nucleotide
	//comes at each position
	if(gSNP) {
		printf("Finding SNPs\n");
#if defined( DISCRETIZE )
		read_allot = new center_d[genome_size];
		if(!read_allot)
			throw new Exception("Could not allocate space for SNP calling");
		memset(read_allot,0,genome_size*sizeof(char));
#elif defined( INTDISC )
		reads = new unsigned char[genome_size*NUM_SNP_VALS];
		if(!reads)
			throw new Exception("Could not allocate space for SNP calling");
		memset(reads,0,genome_size*NUM_SNP_VALS*sizeof(char));
#else
		for(unsigned int i=0; i<NUM_SNP_VALS; i++) {
			reads[i] = new float[genome_size];
			if(!reads[i])
				throw new Exception("Could not allocate space for SNP calling");
			memset(reads[i],0,genome_size*sizeof(float));
		}
#endif
	
	}
	if(gBISULFITE) {
		// Might already have been done under SNP calling
		printf("Using bisulfite conversion techniques\n");
#if defined( DISCRETIZE )	
		if(!read_allot) {
			read_allot = new center_d[genome_size];
			if(!read_allot)
				throw new Exception("Could not allocate space for BISULFITE calling");
			memset(read_allot,0,genome_size*sizeof(char));
		}
#elif defined( INTDISC )
		if(!reads) {
			reads = new unsigned char[genome_size*NUM_SNP_VALS];
			if(!reads)
				throw new Exception("Could not allocate space for Bisulfite calling");
			memset(reads,0,genome_size*NUM_SNP_VALS*sizeof(char));
		}
#else
		for(unsigned int i=0; i<NUM_SNP_VALS; i++) {
			// Might already have been done under SNP calling
			if(!reads[i])
				reads[i] = new float[genome_size];
			if(!reads[i])
				throw new Exception("Could not allocate space for Bisulfite calling");
			memset(reads[i],0,genome_size*sizeof(float));
		}
#endif
		
	}
	if(gATOG) {
		printf("Using AtoG conversion techniques\n");
		// Might already have been done under SNP calling
#if defined( DISCRETIZE )
		if(!read_allot) {
			read_allot = new center_d[genome_size];
			if(!read_allot)
				throw new Exception("Could not allocate space for A to G calling");
			memset(read_allot,0,genome_size*sizeof(char));
		}
#elif defined( INTDISC )
		if(!reads) {
			reads = new unsigned char[genome_size*NUM_SNP_VALS];
			if(!reads)
				throw new Exception("Could not allocate space for Bisulfite calling");
			memset(reads,0,genome_size*NUM_SNP_VALS*sizeof(char));
		}
#else
		for(unsigned int i=0; i<NUM_SNP_VALS; i++) {
			if(!reads[i])
				reads[i] = new float[genome_size];
				
			if(!reads[i])
				throw new Exception("Could not allocate space for A to G calling");
			memset(reads[i],0,genome_size*sizeof(float));
		}
#endif // DISCRETIZE/INTDISC
	}
}

void GenomeBwt::LoadGenome()
{
    //std::cout << "load genome" << std::endl;
    unsigned int read_status = readGen( ref_genome_fn );
    
    if( read_status )
    {
        std::cout << "Genome reference already exists and loaded successfully!" << std::endl;
    }
    else
    {
        std::cout << "Could not find reference genome index! Building now." << std::endl;
        index_and_store();
        readGen( ref_genome_fn );
    }


    for( int chr_id = 0; chr_id < index->bns->n_seqs; chr_id++ )
    {
        names.push_back( pair< string, unsigned long >( index->bns->anns[ chr_id ].name, index->bns->anns[ chr_id ].offset ) );
    }

    names.push_back( pair< string, unsigned long >( "end", index->bns->l_pac ) );

#ifdef DEBUG
		for( unsigned int i = 0; i < names.size(); i++ )
        {
			printf( "Name[%d]: %s\tBegin: %ld\n", i, names[ i ].first.c_str(), names[ i ].second );
		}
		cout.flush();
#endif

	make_extra_arrays();

    // only save if the save flag is set and you haven't read it from disk
    // no need to save if you're already able to read it from disk.    
	if( gSAVE && !gREAD )
    {
		fprintf( stderr, "\tSaving to file: %s\n", gSAVE_FN );
		saveGen( gSAVE_FN );
	}
	
    {
	    amount_genome = ( float* ) malloc( sizeof( float ) * index->bns->l_pac / gGEN_SIZE );
	    gs_positions = ( char* ) malloc( index->bns->l_pac );
	    fprintf( stderr,"Success!\n" );
	    gs_positions_size = index->bns->l_pac;
    }
}

void GenomeBwt::index_and_store()
{
    //std::cout << "index_and_store" << std::endl;
    //std::cout << "trying to index:\t" << ref_genome_fn << std::endl;

    //extern int bwa_index(int argc, char *argv[]);

    char* indexing_args[] = { "index", ref_genome_fn };
    bwa_index( 2, indexing_args );
}

void GenomeBwt::ReplaceSpaceWithUnderscore(string &str)
{
	for(unsigned int i=0; i<str.length(); i++) {
		if(str[i] == ' ' || str[i] == 10 || str[i] == 11 || str[i] == 12 || str[i] == 13)
			str[i] = '_';
	}
}


void GenomeBwt::StoreGenome( bool make_extra )
{
    StoreGenome();
}

#define READ_EXTRA	1000

void GenomeBwt::StoreGenome() {
    //std::cout << "store genome" << std::endl;
	// We'll set it in here so we don't need to take care of it elsewhere
	gMER_SIZE = MAX_MER_SIZE;
	
	// Don't hash anything--not needed
	include_hash = false;
	index_and_store();
	
	// Still need to create extra arrays.
	make_extra_arrays();
}

/*
 * This function will just count the genome--it won't store anything.
 */
unsigned long GenomeBwt::count()
{
    return index->bns->l_pac;
}

unsigned int GenomeBwt::saveGen(char* fn)
{
    //deprecated, no need to save anymore
    std::cout << "GenomeBwt::saveGen Deprecated" << std::endl;
    return 0;
}

/*
 * Returns a string from begin to begin+size.
 * the first position in the genome is position 0.
 * @par begin - the beginning of the string to return
 * @par size - the size of the string to return.
 */
string GenomeBwt::GetString(const unsigned long begin, const unsigned int size)
{
#ifdef DEBUG
	fprintf(stderr,"[%d] Getting string at %lu",iproc,begin);
#endif

    uint64_t l_pac = index->bns->l_pac;
    uint8_t* rseq;
    int64_t rlen;

    // check to see if the sequence falls on a chromosome boundary
    int rid = bns_intv2rid( index->bns, begin, begin + size );
    if( rid < 0 )
    {
        return "";
    }

    rseq = bns_get_seq(index->bns->l_pac, index->pac, begin, begin + size, &rlen);
    
	string ref_seq( rlen, 'n' );
    for( int i = 0; i < rlen; ++i )
    {
        ref_seq[ i ] = "acgtn"[ ( int ) rseq[ i ] ];
    }
    
#ifdef DEBUG
    fprintf(stderr,"\n");
#endif

	return ref_seq;
}

inline char GenomeBwt::GetChar(const unsigned long pos) {

	/*unsigned long pos_in_gen = (pos-my_start)/GEN_PACK;
	assert(pos_in_gen <= genome_size/GEN_PACK);

	GEN_TYPE t_temp = packed_genome[pos_in_gen];
	unsigned long bit_offset = (pos-my_start)%GEN_PACK;

	return gINT2BASE[0xf & (t_temp >> (4*(GEN_PACK-bit_offset-1)))];*/
    
    int64_t rlen;
    return "acgtn"[ bns_get_seq( index->bns->l_pac, index->pac, pos, pos + 1, &rlen )[ 0 ] ];
}

uint64_t GenomeBwt::get_sa_coord( uint64_t sa_pos )
{
    //std::cout << "\t" << index->bwt->sa[ sa_pos ] << std::endl;
    //std::cout << "\t" << bwt_sa( index->bwt, sa_pos ) << std::endl;
    return bwt_sa( index->bwt, sa_pos );
}

void GenomeBwt::get_sa_int( string& seq, uint64_t* in_start, uint64_t* in_end )
{
    unsigned char* c_seq = new unsigned char[ seq.size() ];

    for (int i = 0; i < seq.size(); ++i) // convert to 2-bit encoding if we have not done so
    {
        //std::cout << "\ti:\t" << i << "\t" << seq[ i ] << "\t" << ( int ) seq[ i ] << std::endl;
        c_seq[i] = seq[i] < 4? seq[ i ] : nst_nt4_table[ ( int ) seq[ i ] ];
    }
   
    //printf( "%d\n", c_seq ); 
    //std::cout << "c seq:\t" << ( long ) c_seq << std::endl;
    //std::cout << "t seq:\t" << seq << std::endl;

    //printf( "sa size:\t%d\n", index->bwt->n_sa );
    //printf( "bwt size:\t%d\n", index->bwt->bwt_size );
    //printf( "seq size:\t%d\n", index->bwt->seq_len );
    bwtint_t start, end;

    int search_result = bwt_match_exact( index->bwt, seq.size(), c_seq, &start, &end );
    //std::cout << "search result:\t" << search_result << std::endl;

    if( search_result > 0 )
    {
        *in_start = start;
        *in_end = end;
    }
    else
    {
        *in_start = -1;
        *in_end = -1;
    }

    //std::cout << "start:\t" << start << "\nend:\t" << end << std::endl;
}

/*HashLocation* GenomeBwt::GetMatches(string &str) {
	pair<int, unsigned long> p_hash = bin_seq::get_hash(str);

    std::cout << "GetMatches" << std::endl;
    std::cout << "\tfind:\t" << str << std::endl;

	if(!p_hash.first)
		return NULL;
		
	//return &gen_hash[p_hash.second];
    return NULL;
}

HashLocation* GenomeBwt::GetMatches(unsigned int hash) {

	//return &gen_hash[hash];
    return NULL;
}*/

float GenomeBwt::GetScore(const unsigned long &pos) {

	unsigned long pos_temp = (pos-my_start)/gGEN_SIZE;
	//return amount_genome[pos_temp];
	return amount_genome[pos];
}

void GenomeBwt::AddScore(const unsigned long &pos, const float &amt) {
//	if(amt > 1)
//		fprintf(stderr,"Genome Error!! %d\n",__LINE__);
		
	//unsigned long pos_temp = pos-my_start;
	//amount_genome[pos_temp/gGEN_SIZE] += amt;
	amount_genome[ pos / gGEN_SIZE ] += amt;
}

/**
 * INVARIANT:  must add to the genome previous to this call
 *			thus, amount_genome[pos] will never be 0
 */
void GenomeBwt::AddSeqScore(unsigned long pos, const float* amt, const float scale) {
	//fprintf(stderr,"pos: %lu, myStart: %lu\n",pos,my_start);
	unsigned int loc = (pos-my_start)/gGEN_SIZE;
	
/*
	if(amt[A_POS] * scale > 0.0)
		fprintf(stderr,"For position %lu, total A is %f = %f*%f\n",pos,amt[A_POS]*scale,amt[A_POS],scale);
	if(amt[C_POS] * scale > 0.0)
		fprintf(stderr,"For position %lu, total C is %f = %f*%f\n",pos,amt[C_POS]*scale,amt[C_POS],scale);
	if(amt[G_POS] * scale > 0.0)
		fprintf(stderr,"For position %lu, total G is %f = %f*%f\n",pos,amt[G_POS]*scale,amt[G_POS],scale);
	if(amt[T_POS] * scale > 0.0)
		fprintf(stderr,"For position %lu, total T is %f = %f*%f\n",pos,amt[T_POS]*scale,amt[T_POS],scale);
	if(amt[N_POS] * scale > 0.0)
		fprintf(stderr,"For position %lu, total N is %f = %f*%f\n",pos,amt[N_POS]*scale,amt[N_POS],scale);
//*/

#if defined( DISCRETIZE )
	int i;
	center_t *map_c = centers[read_allot[loc]];
	center_t gen_amt[DIM_CENTERS];
	// Convert to floating point representation
	for(i=0; i<DIM_CENTERS; i++)
		// we've added in the value previous to this call.  Subtract it out first
		gen_amt[i] = map_c[i]*(amount_genome[loc]-scale);
	
	gen_amt[A_POS] += amt[A_POS]*scale;
	gen_amt[C_POS] += amt[C_POS]*scale;
	gen_amt[G_POS] += amt[G_POS]*scale;
	gen_amt[T_POS] += amt[T_POS]*scale;
	gen_amt[N_POS] += amt[N_POS]*scale;


	// Convert back, including the added value
	for(i=0; i<DIM_CENTERS; i++)
		gen_amt[i] /= amount_genome[loc];
	read_allot[loc] = find_center(gen_amt);

#elif defined( INTDISC )
	// What do we need to do here?  Use the "parts per 255" again
	// Need to scale the variable "scale", then we can add it to the position
	// then round it to the closest 255
	float loc_amt = amount_genome[loc];
	float prev_amt = loc_amt-scale;
	unsigned int reads_loc = loc*NUM_SNP_VALS;
	unsigned char newC;
	
	//fprintf(stderr,"prev: %d %d %d %d %d\n",reads[reads_loc+A_POS],reads[reads_loc+C_POS],reads[reads_loc+G_POS],reads[reads_loc+T_POS],reads[reads_loc+N_POS]);
	/*
	float fnA = FLTFROM255(prev_amt, reads[reads_loc+A_POS]) + amt[A_POS]*scale;
	float fnC = FLTFROM255(prev_amt, reads[reads_loc+C_POS]) + amt[C_POS]*scale;
	float fnG = FLTFROM255(prev_amt, reads[reads_loc+G_POS]) + amt[G_POS]*scale;
	float fnT = FLTFROM255(prev_amt, reads[reads_loc+T_POS]) + amt[T_POS]*scale;
	float fnN = FLTFROM255(prev_amt, reads[reads_loc+N_POS]) + amt[N_POS]*scale;
	fprintf(stderr,"prev_floats(%f): %f %f %f %f %f\n",prev_amt,fnA-amt[A_POS]*scale,fnC-amt[C_POS]*scale,fnG-amt[G_POS]*scale,fnT-amt[T_POS]*scale,fnN-amt[N_POS]*scale);
	fprintf(stderr,"floats(%f): %f %f %f %f %f\n",loc_amt,fnA,fnC,fnG,fnT,fnN);
	*/

	for(unsigned int i=0; i<NUM_SNP_VALS; i++) {
		newC = CHARFROMFLT(loc_amt, FLTFROM255(prev_amt, reads[reads_loc+i]) + amt[i]*scale);
		reads[reads_loc+i] = newC;
	}

	//fprintf(stderr,"post: %hu %hu %hu %hu %hu\n",reads[reads_loc+A_POS],reads[reads_loc+C_POS],reads[reads_loc+G_POS],reads[reads_loc+T_POS],reads[reads_loc+N_POS]);
#else
	reads[A_POS][loc] += amt[A_POS] * scale;
	reads[C_POS][loc] += amt[C_POS] * scale;
	reads[G_POS][loc] += amt[G_POS] * scale;
	reads[T_POS][loc] += amt[T_POS] * scale;

	if(NUM_SNP_VALS > N_POS)
		reads[N_POS][loc] += amt[N_POS] * scale;
#ifdef _INDEL
	if(NUM_SNP_VALS > DEL_POS)
		reads[DEL_POS][loc] += amt[DEL_POS] * scale;
	if(NUM_SNP_VALS > INS_POS)
		reads[INS_POS][loc] += amt[INS_POS] * scale;
#endif // _INDEL
#endif // DISCRETIZE/INTDISC
}

/* 
 * Invariant:  We must call "AddScore" before calling this function
 */
void GenomeBwt::AddSeqScore(unsigned long pos, const float amt, unsigned int which) {
	unsigned int loc = (pos-my_start)/gGEN_SIZE;
	//fprintf(stderr,"pos: %lu, myStart: %lu, loc:%u \n",pos,my_start,loc);

#if defined( DISCRETIZE )
	// If these haven't been instantiated, don't make any changes
	if(!read_allot)
		return;

	int i;
	//fprintf(stderr,"At loc %u, read_allot is %d\n",loc,read_allot[loc]);
	center_t *map_c = centers[read_allot[loc]];
	center_t gen_amt[DIM_CENTERS];
	for(i=0; i<DIM_CENTERS; i++)
		gen_amt[i] = map_c[i]*(amount_genome[loc]-amt);
	
	gen_amt[which] += amt;
	for(i=0; i<DIM_CENTERS; i++)
		gen_amt[i] /= amount_genome[loc];

	read_allot[loc] = find_center(gen_amt);
#elif defined( INTDISC )
	if(!reads)
		return;

	float loc_amt = amount_genome[loc];
	float prev_amt = loc_amt-amt;
	unsigned int reads_loc = loc*NUM_SNP_VALS;
	unsigned char newC;
	
	for(unsigned int i=0; i<NUM_SNP_VALS; i++) {
		if(i == which) {
			newC = CHARFROMFLT(loc_amt, FLTFROM255(prev_amt, reads[reads_loc+which]) + amt);
			reads[reads_loc+which] = newC;
		}
		else {
			newC = CHARFROMFLT(loc_amt, FLTFROM255(prev_amt, reads[reads_loc+i]));
			reads[reads_loc+i] = newC;
		}
	}
#else
	if(reads[which])
		reads[which][loc] += amt;
#endif	// DISCRETIZE/INTDISC

}

string GenomeBwt::GetPos(const pair<unsigned long,int> p_pos) {
	return GetPos(p_pos.first,p_pos.second);
}
string GenomeBwt::GetPos(const unsigned long &pos, const int strand) {
	int chr_num = 0;

	while(names[chr_num].second <= pos && ((unsigned int)chr_num < names.size()) )
		chr_num++;
		
	ostringstream to_return;
	
	chr_num--;
	to_return << names[chr_num].first << ":";
	if(strand == POS_STRAND)
		to_return << "+";
	else
		to_return << "-";
	to_return << (pos - names[chr_num].second);
	
	return to_return.str();
}

pair<string,unsigned long> GenomeBwt::GetPosPair(const unsigned long pos)
{
    int chr_id = bns_pos2rid( index->bns, pos );
    int chr_base_pos = pos - index->bns->anns[ chr_id ].offset;

	return pair< string, unsigned long >( index->bns->anns[ chr_id ].name, chr_base_pos );
}

unsigned long GenomeBwt::GetAbsolutePosition(const char* chr, const unsigned long &pos) {
	unsigned long absPos = 0;
	bool found = false;
	
	//fprintf(stderr,"Size of names is %u\n",names.size());
	for(unsigned int i=0; i<names.size(); i++) {
		if(names[i].first.compare(chr) == 0) {
			absPos = names[i].second+pos;
			found = true;
		}
	}
	
	if(!found) {
		fprintf(stderr,"Chromosome is >%s<\n",chr);
		string excp = "Could not find chromosome in Genome file: ";
		excp += chr;
		throw new Exception(excp);
	}

	return absPos;
}

string GenomeBwt::GetStringRelative(const char* chr, const unsigned long &pos, const unsigned int length) {
	unsigned long abs_pos = GetAbsolutePosition(chr,pos);
	return GetString(abs_pos,length);
}


#ifdef SET_POS
bool GenomeBwt::getSetThreadID(long unsigned int pos, unsigned int tid) {
	if(pos < my_start)
		return false;

	//fprintf(stderr,"pos=%lu, genome_size=%lu, my_start=%lu, final_pos=%lu\n",pos,genome_size,my_start,pos-my_start);
	pos -= my_start;

	//assert(pos < genome_size);
	assert( pos < index->bns->l_pac );

	if( gs_positions[ pos / GEN_PACK ] & ( 1 << ( tid % 8 ) ) )
    {
		return true;
    }

	unsigned int new_tid = 1<<(tid%8);
	gs_positions[pos/GEN_PACK] |= new_tid;
	return false;
}

void GenomeBwt::unsetThreadID(long unsigned int pos, unsigned int tid) {
	if(pos < my_start)
		return;

	pos -= my_start;


	//assert(pos < genome_size);
	assert( pos < index->bns->l_pac );
	// We'll clear all the threads at this position
	gs_positions[pos/GEN_PACK] = 0;
}
#endif

inline bool GenomeBwt::is_max_pos(char to_check, float a, float c, float g, float t, float n) {
	switch(to_check) {
		case 'a':
			return !(a >= c && a >= g && a >= t && a >= n);
		case 'c':
			return !(c >= a && c >= g && c >= t && c >= n);
		case 'g':
			return !(g >= c && g >= a && g >= t && g >= n);
		case 't':
			return !(t >= c && t >= g && t >= a && t >= n);
		case 'n':
			return !(n >= c && n >= g && n >= t && n >= a);
	}

	return false;
}

inline unsigned int GenomeBwt::max_pos(float array[], unsigned int array_size) {
	float* pos = max_element(array, array+array_size);
	unsigned int max_pos = pos - &array[0];
	return max_pos;
}

const unsigned int MIN_COVERAGE = 5;
const float MIN_RATIO = 1.5;
const float MAX_PVAL = 0.01;
const float MAX_RATIO = 0.0;

inline double GenomeBwt::LRT(float chars[5], int &snp_pos) {
	snp_pos = max_pos(chars, 5);

	double sum = chars[0] + chars[1] + chars[2] + chars[3] + chars[4];

	double ratio = pow(.2,sum) 
					/ ( pow(chars[snp_pos]/sum,(double)chars[snp_pos]) 
					    * pow((sum-chars[snp_pos])/sum/4,sum-chars[snp_pos]) );
	
	//double pval = 1-pchisq(-2*log(ratio),1);
	double pval = 1-gsl_cdf_chisq_P(-2*log(ratio),1);
	
	return pval;
	
}

#define MONO_DIP_RATIO	3.0f
#define DIFF_UNIF_PRIOR	0.2

inline double GenomeBwt::dipLRT(float chars[NUM_CHAR_VALS], int &snp_pos1, int &snp_pos2, bool &dip) {

	/*******************************************
	 * If the ratio is too great, we're just going 
	 * to define it as a monoploid SNP
	 *******************************************/
	snp_pos1 = max_pos(chars, NUM_CHAR_VALS);
		
	double ratio1,ratio2;
	double sum = 0;
	for(unsigned int i=0; i<NUM_CHAR_VALS; i++) {
		sum += chars[i];
	}

	ratio1 = pow(.2,sum) 
					/ (pow(chars[snp_pos1]/sum,(double)chars[snp_pos1]) 
					  * pow(((sum-chars[snp_pos1])/sum)/(NUM_CHAR_VALS-1),sum-chars[snp_pos1]));
					
	double pval1 = 1-gsl_cdf_chisq_P(-2*log(ratio1),1);
	
	// Find the second highest SNP
	float second[NUM_CHAR_VALS];
	for(unsigned int i=0; i<NUM_CHAR_VALS; i++) {
		second[i] = chars[i];
	}
	second[snp_pos1] = 0;
	snp_pos2 = max_pos(second, NUM_CHAR_VALS);
	
	double pval2;
	if(chars[snp_pos1]/chars[snp_pos2] > MONO_DIP_RATIO || snp_pos1 == snp_pos2) {
		pval2 = MAX_PVAL;
		ratio2 = MAX_RATIO;
		snp_pos2 = -1;
	}
	else {
		// Add in the uniform prior
		for(unsigned int i=0; i<NUM_CHAR_VALS; i++) {
			chars[i] += DIFF_UNIF_PRIOR;
			sum += DIFF_UNIF_PRIOR;
		}
		
		// Calculate a new ratio and pval with the uniform prior
		ratio1 = pow(.2,sum) 
					/ (pow(chars[snp_pos1]/sum,(double)chars[snp_pos1]) 
					  * pow(((sum-chars[snp_pos1])/sum)/(NUM_CHAR_VALS-1),sum-chars[snp_pos1]));
		pval1 = 1-gsl_cdf_chisq_P(-2*log(ratio1),1);

		// Calculate the diploid ratio and prior
		//rat = (.2)^sum/((x[m]/sum)^(x[m])*(x[m2]/sum)^(x[m2])
		//		* (1-(x[m]+x[m2])/sum)^(sum-x[m]-x[m2]))
		double numerator = pow(.2,sum);

		double denominator = ( pow(chars[snp_pos1]/sum, (double)chars[snp_pos1])
				    * pow(chars[snp_pos2]/sum, (double)chars[snp_pos2])
				    * pow(((sum - chars[snp_pos1] + chars[snp_pos2])/sum)/3, 
				    	  sum - chars[snp_pos1] - chars[snp_pos2]));	
		ratio2 = numerator/denominator;

#ifdef DEBUG
		fprintf(stderr,"x is %f %f %f %f %f %f %f\n",chars[0],chars[1],chars[2],chars[3],chars[4],chars[5],chars[6]);
		fprintf(stderr,"ratio2 is %f (.2^%f/denom=%f/%f)\n",ratio2,sum,numerator,denominator);
		fprintf(stderr,"Result from gsl is %e\n",gsl_cdf_chisq_P(-2*ratio2, 2));
#endif


		//pval2 = 1-pchisq(-2*log(ratio2), 2);
		pval2 = 1-gsl_cdf_chisq_P(-2*log(ratio2), 2);

		// Remove the uniform prior
		for(unsigned int i=0; i<NUM_CHAR_VALS; i++) {
			chars[i] -= DIFF_UNIF_PRIOR;
			sum -= DIFF_UNIF_PRIOR;
		}

	}
	
#ifdef DEBUG
	fprintf(stderr,"pval1 is %f and pval2 is %f\n",pval1,pval2);
	fprintf(stderr,"Ratio is %f/%f=%f\n",chars[snp_pos1],chars[snp_pos2],chars[snp_pos1]/chars[snp_pos2]);
#endif
	
	// They're too high for a chi-square distribution, use the ratio
	if(pval2 == 0 && pval1 == 0) {
		if(ratio2 < ratio1 && (chars[snp_pos1]/chars[snp_pos2] < MONO_DIP_RATIO)) {
			dip = true;
#ifdef DEBUG
			fprintf(stderr,"Result: Diploid.  Pval is %f\n\n",pval2);
#endif
		}
		else {
			dip = false;
#ifdef DEBUG
			fprintf(stderr,"Result: Monoploid.  Pval is %f\n\n",pval2);
#endif
		}

		return 0;
	}

	if(pval2 < pval1 && (chars[snp_pos1]/chars[snp_pos2] < MONO_DIP_RATIO)) {
	//if(ratio2 < ratio1 && (chars[snp_pos1]/chars[snp_pos2] < MONO_DIP_RATIO)) {
#ifdef DEBUG
		fprintf(stderr,"Result: Diploid.  Pval is %f\n\n",pval2);
		fprintf(stderr,"Y pval2: %e<pval1: %e, chars[1]=%f,chars[2]=%f\n",pval2,pval1,chars[snp_pos1],chars[snp_pos2]);
#endif
		dip = true;
		return pval2;
	}
	else {
#ifdef DEBUG
		fprintf(stderr,"Result: Monoploid.  Pval is %f\n\n",pval2);
#endif
		dip = false;
		return pval1;
	}
}

inline double GenomeBwt::is_snp(float chars[NUM_CHAR_VALS], int &snp_pos1, int &snp_pos2, bool &dip) {
	
	/*
	//add a diffuse/uniform prior to chars
	for(int i=0; i<5; i++) {
		chars[i] += DIFF_UNIF_PRIOR;
	}
	*/

	double pval;
	if(gSNP_MONOP) {
		dip = false;
		pval = LRT(chars,snp_pos1);
	}
	else
		pval = dipLRT(chars,snp_pos1,snp_pos2,dip);

	/*
	// remove the prior from chars
	for(int i=0; i<5; i++) {
		chars[i] -= DIFF_UNIF_PRIOR;
	}
	*/

	return pval;

}

/**
 * is_indel will determine if this is an indel or not. 
 */
//inline int GenomeBwt::is_indel

void GenomeBwt::AppendFinal(const char* fn) {
	PrintFinal(fn, true);
}
void GenomeBwt::PrintFinal(const char* fn) {
	PrintFinal(fn, false);
}

void GenomeBwt::PrintFinal(const char* fn, bool append) {
	if(gSNP) {
		PrintFinalSNP(fn, append);
	}
	// This one is still sgrex, but not SNP calling
	else if(gBISULFITE || gATOG) {	// if any of the reads have been instantiated...
		PrintFinalBisulfite(fn, append);
	}
	else {
		PrintFinalSGR(fn, append);
	}
}

#define MIN_PRINT	0.001

void GenomeBwt::PrintFinalSNP(const char* fn, bool append) {
	// Instead of printing out a file for each, we want to print
	// out one file (dubbed the '.sgrex' file for 'sgr-extra') to be a tab-deliminated file
	// for the number of a's, c's, g's, t's, and gaps at each location in the genome.
	char gmp_fn[strlen(fn)+10];
	strcpy(gmp_fn,fn);
	strcat(gmp_fn,".gmp");
	printf("SNP File is: %s\n",gmp_fn);
	
	FILE* gmp_file;
	if(append)
		gmp_file = fopen(gmp_fn,"a");
	else
		gmp_file = fopen(gmp_fn,"w");
		
	if(!gmp_file)
		perror("Could not open .gmp file");
		
	// don't output the Header
	//fprintf(sgrex_file,"chr\tpos\tsgr_val\tread_a\tread_c\tread_g\tread_t\tread_n\tread_gap\tsnp\n");
	//fprintf(sgrex_file,"#chr\tpos\tsgr_val\tread_a\tread_c\tread_g\tread_t\tread_n\tsnp\n");
	
	if(gVERBOSE > 1) 
		fprintf(stderr,"NAMES: %ld\n",names.size());
	
	unsigned long count=0;
	for(unsigned int j=0; j<names.size()-1; j++) {
	//for(unsigned int j=0,count=0; names[j].first != "end"; j++) {
		if(gVERBOSE > 1) 
			fprintf(stderr,"\t[%d]: %s  %lu\n",j+1,names[j+1].first.c_str(),names[j+1].second);
		for(; count<names[j+1].second; count+=gGEN_SIZE) {
			if(count < my_start && !read_all)
				continue;
			if(count > my_end && !read_all)
				break;
				
			unsigned long locus = (count-my_start)/gGEN_SIZE;

			//if(genome[locus].amount > (float)0.0) {	//if something matched here...
			if(amount_genome[locus] > MIN_PRINT) {	//if something matched here...
			
								 //  chr pos  sgr -- want 1-based indexing
				fprintf(gmp_file, "%s\t%ld\t%.5f", names[j].first.c_str(), (count-names[j].second)+1,
							//genome[locus].amount);
							amount_genome[locus]);

				// Print out each of the read positions
#if defined( DISCRETIZE )
				unsigned int dim;
				center_t *map_c = centers[read_allot[locus]];
				center_t read_amt[DIM_CENTERS];
				for(dim=0; dim<DIM_CENTERS; dim++) {
					read_amt[dim] = map_c[dim]*amount_genome[locus];
				}
				for(dim=0; dim<DIM_CENTERS; dim++) {
					fprintf(gmp_file,"\t%.5f",read_amt[dim]);
				}
#elif defined( INTDISC )
				for(unsigned int read_count=0; read_count<NUM_SNP_VALS; read_count++) {
					fprintf(gmp_file,"\t%.5f",reads[locus*NUM_SNP_VALS+read_count]/255
									*amount_genome[locus]);
				}
#else
				for(unsigned int read_count=0; read_count<NUM_SNP_VALS; read_count++) {
					fprintf(gmp_file,"\t%.5f",reads[read_count][locus]);
				}
				
#endif // DISCRETIZE/INTDISC

				PrintSNPCall(count, gmp_file);
				
				// Carriage return
				fprintf(gmp_file,"\n");
			}
		}
	}
	
	// Close the file
	fclose(gmp_file);
}

void GenomeBwt::PrintSNPCall(const unsigned long count, FILE* gmp_file) {
	// Get the position here
	char at_pos = GetChar(count);
	unsigned long locus = (count-my_start)/gGEN_SIZE;
	// Determine if it's a snp
#if defined( DISCRETIZE )
	center_t *map_c = centers[read_allot[locus]];
	float char_arr[DIM_CENTERS];
	for(unsigned int i=0; i<DIM_CENTERS; i++)
		char_arr[i] = map_c[i]*amount_genome[locus];
#elif defined( INTDISC )
	float char_arr[NUM_SNP_VALS];
	for(unsigned int i=0; i<NUM_SNP_VALS; i++) {
		char_arr[i] = reads[locus*NUM_SNP_VALS+i]/255
								*amount_genome[locus];
	}
#else
#ifdef _INDEL										
	float char_arr[NUM_SNP_VALS] =
						{	reads[0][locus],
							reads[1][locus],
							reads[2][locus],
							reads[3][locus],
							reads[4][locus],
							reads[5][locus],
							reads[6][locus]};
#else
	float char_arr[NUM_SNP_VALS] =
						{	reads[0][locus],
							reads[1][locus],
							reads[2][locus],
							reads[3][locus],
							reads[4][locus]};
#endif // _INDEL
#endif // DISCRETIZE/INTDISC

	// Sanity check
	float sum = 0;
	for(unsigned int i=0; i<NUM_SNP_VALS; i++) {
		sum += char_arr[i];
	}
	if(gVERBOSE > 1) {
		if(amount_genome[(count-my_start)/gGEN_SIZE] - sum > 0.01) {
			fprintf(stderr,"ERROR at %lu:  %f vs sum of %f\n",(count-my_start)/gGEN_SIZE,amount_genome[(count-my_start)/gGEN_SIZE],sum);
		}
	}

	int snp_pos1, snp_pos2;
	bool dip = false;
	double p_val = is_snp(char_arr, snp_pos1, snp_pos2, dip);


	// If the max pos is not what's in the genome
	// 	or it's a diploid snp
	if(snp_pos1 != g_gen_CONVERSION[(int)at_pos] || dip) {
		if( p_val < gSNP_PVAL) {
			if(dip) {
				fprintf(gmp_file,"\tY:%c->%c/%c p_val=%.2e", 
					at_pos, gINT2BASE[snp_pos1], gINT2BASE[snp_pos2], p_val);
			}
			else {
				fprintf(gmp_file,"\tY:%c->%c p_val=%.2e", 
					at_pos, gINT2BASE[snp_pos1], p_val);
			}
		}
		else
			if(dip) {
				fprintf(gmp_file,"\tN:%c->%c/%c p_val=%.2e", 
					at_pos, gINT2BASE[snp_pos1], gINT2BASE[snp_pos2], p_val);
			}
			else {
				fprintf(gmp_file,"\tN:%c->%c p_val=%.2e", 
					at_pos, gINT2BASE[snp_pos1], p_val);
			}
	}
	else
		//fprintf(gmp_file,"\tN:%c p_val=%.2e", at_pos, p_val);
		fprintf(gmp_file,"\tN");

}

void GenomeBwt::PrintFinalBisulfite(const char* fn, bool append) {
	// We got a new update from Evan.  Instead of printing out a file for each, we want to print
	// out one file (dubbed the '.gmp' file--previously '.sgrex') to be a tab-deliminated file
	// for the number of a's, c's, g's, t's, and gaps at each location in the genome.
	char gmp_fn[strlen(fn)+10];
	strcpy(gmp_fn,fn);
	strcat(gmp_fn,".gmp");
	if(gBISULFITE)
		printf("Bisulfite File is: %s\n",gmp_fn);
	if(gATOG)
		printf("RNA edit File is: %s\n",gmp_fn);

	FILE* gmp_file;
	if(append)
		gmp_file = fopen(gmp_fn,"a");
	else
		gmp_file = fopen(gmp_fn,"w");
		
	if(!gmp_file)
		perror("Could not open .gmp file");
	// don't output the Header
	//fprintf(gmp_file,"chr\tpos\tsgr_val\tread_a\tread_c\tread_g\tread_t\tread_n\tread_gap\n");
	//fprintf(gmp_file,"#chr\tpos\tsgr_val\tread_a\tread_c\tread_g\tread_t\tread_n\n");
	
	if(gVERBOSE > 1) 	
		fprintf(stderr,"Size of NAMES: %ld\n",names.size());
	
	double total_stuff = 0;
	
	unsigned long count=0;
	for(unsigned int j=0; j<names.size()-1; j++) {
		if(gVERBOSE > 1) 
			fprintf(stderr,"\t[%d/%lu]: %s\n",j,names.size(),names[j].first.c_str());
		for(; count<names[j+1].second; count+=gGEN_SIZE) {
			if(count < my_start && !read_all)
				continue;
			if(count > my_end && !read_all)
				break;

			unsigned long locus = (count-my_start)/gGEN_SIZE;

			char at_pos = GetChar(count);
			int gen_count_spot;	// Later, we'll ensure there's actually info here
			
			// Check if it's bisulfite or RNA editing we're doing
			if(gBISULFITE) {
				if(gMATCH_POS_STRAND && !gBISULFITE2) {
					// only print if there is a 'c' nucleotide here
					if(at_pos != 'c')
						continue;
					gen_count_spot = 1;
				}
				else {
					if(at_pos != 'g')
						continue;
					gen_count_spot = 2;
				}
			}
			else {	// gATOG
				if(gMATCH_POS_STRAND) {
					if(at_pos != 'a')
						continue;
					gen_count_spot = 0;
				}
				else {	// gMATCH_NEG_STRAND
					if(at_pos != 't')
						continue;
					gen_count_spot = 3;
				}
			}
			
			total_stuff += amount_genome[locus];
			
			//if(genome[locus].amount > (float)0.0) {	//if something matched here...
			if(amount_genome[locus] > 0.0f) {	//if something matched here...
			//if(reads[gen_count_spot][locus] > MIN_PRINT) {	//if something matched here...
								 //  chr pos  sgr -- want 1-based indexing
				fprintf(gmp_file, "%s\t%ld\t%f", names[j].first.c_str(), (count-names[j].second)+1,
							//genome[locus].amount);
							amount_genome[locus]);

				// Print out each of the read positions
#if defined( DISCRETIZE )
				unsigned int dim;
				center_t *map_c = centers[read_allot[locus]];
				center_t read_amt[DIM_CENTERS];
				for(dim=0; dim<DIM_CENTERS; dim++) {
					read_amt[dim] = map_c[dim]*amount_genome[locus];
				}
				for(dim=0; dim<DIM_CENTERS; dim++) {
					fprintf(gmp_file,"\t%.5f",read_amt[dim]);
				}
#elif defined( INTDISC )
				for(unsigned int read_count=0; read_count<NUM_CHAR_VALS; read_count++) {
					fprintf(gmp_file,"\t%.5f",
								reads[locus*NUM_SNP_VALS+read_count]/255 
								*amount_genome[locus]);
				}
#else
				// Only print out acgtn, not indel
				for(unsigned int read_count=0; read_count<NUM_CHAR_VALS; read_count++) {
					fprintf(gmp_file,"\t%.5f",reads[read_count][locus]);
				}

#endif // DISCRETIZE/INTDISC

				// Just for testing right now
				//PrintSNPCall(count,gmp_file);
				
				// Carriage return
				fprintf(gmp_file,"\n");
			}
		}
	}
	//fprintf(stderr,"Total stuff=%f\n",total_stuff);
	
	// Close the file
	fclose(gmp_file);
}

void GenomeBwt::PrintFinalSGR(const char* fn, bool append) {
	char* filename = new char[strlen(fn)+5];
	sprintf(filename,"%s.sgr",fn);
	
	FILE* sgr_file;
	if(append)
		sgr_file = fopen(filename,"a");
	else
		sgr_file = fopen(filename,"w");
	
	//fprintf(stderr,"Printing for %u genomes to %lu\n",names.size(),names[names.size()-1].second);
	unsigned long count = 0;
    /*for( unsigned int i = 0; i < index->bns->n_seqs - 1; i++)
    {
		for(; count < names[ i + 1 ].second; count += gGEN_SIZE )
        {
			if(count < my_start && !read_all)
            {
				continue;
            }
			if(count > my_end && !read_all)
            {
				break;
            }

			//if(genome[(count-my_start)/gGEN_SIZE].amount> 0)
			if( amount_genome[ ( count - my_start ) / gGEN_SIZE ] > MIN_PRINT )
            {
                // chr pos  sgr -- want 1-based indexing
                fprintf(sgr_file, "%s\t%ld\t%.5f\n",names[i].first.c_str(),(count-names[i].second)+1,
                amount_genome[(count-my_start)/gGEN_SIZE]);
			}
		}
	}*/

	for(unsigned int i=0; i<names.size()-1; i++)
    {
		for(; count < names[ i + 1 ].second; count += gGEN_SIZE )
        {
			if(count < my_start && !read_all)
            {
				continue;
            }
			if(count > my_end && !read_all)
            {
				break;
            }

			//if(genome[(count-my_start)/gGEN_SIZE].amount> 0)
			if( amount_genome[ ( count - my_start ) / gGEN_SIZE ] > MIN_PRINT )
            {
                // chr pos  sgr -- want 1-based indexing
                fprintf(sgr_file, "%s\t%ld\t%.5f\n",names[i].first.c_str(),(count-names[i].second)+1,
                amount_genome[(count-my_start)/gGEN_SIZE]);
			}
		}
	}
	
	fclose(sgr_file);
	printf("Finished printing to %s\n",filename);
	delete[] filename;
}

unsigned long GenomeBwt::size() {
	//return genome_size;
    return index->bns->l_pac;
}

bool GenomeBwt::Test(ostream &os, unsigned int &warnings) {
	bool success = true;
	GenomeBwt gensnp;
	
	// Test the max_pos function
	float a[7] = {0, 1, 2, 3, 4, 5, 6};
	TEST(gensnp.max_pos(a,7) == 6);
	a[5] = 7;
	TEST(gensnp.max_pos(a,7) == 5);
	
	
	float gamma_val = 4/2.0;
	fprintf(stderr,"Gamma(%d/2=%f)=%f, should be 1.0\n",4,gamma_val,tgamma(gamma_val));
	fprintf(stderr,"pchisq(%d,%d)=%f, should be 0.9544997\n",4,1,gsl_cdf_chisq_P(4,1));
	//TEST(dchisq(4,1) == 0.026995);
	TEST(tgamma(gamma_val) == 1.0);

	int i=0,j;
	bool b;
	float x[] = {0,0.000001, 0.008909,0, 0, 0, 0, 0};
	float pval = gensnp.is_snp(x,i,j,b);
	fprintf(stderr,"%f, %f, %f, %f, %f\n",x[0],x[1],x[2],x[3],x[4]);
	fprintf(stderr,"PVAL: %f at %d and second is %d (PVAL should be 0.7359 and second should be -1)\n",pval,i,b?j:-1);
	x[0] = 0;
	x[1] = 0.00172;
	x[2] = 0;
	x[3] = 7.44906;
	x[4] = 0;
	pval = gensnp.is_snp(x,i,j,b);
	fprintf(stderr,"%f, %f, %f, %f, %f\n",x[0],x[1],x[2],x[3],x[4]);
	fprintf(stderr,"PVAL: %.1e at %d and second is %d (PVAL should be 9.2e-06 and second should be -1)\n",pval,i,b?j:-1);
	float monop[] = {0.00084, 5.76128, 0.99737, 0.00077, 0.00000, 0, 0};
	float monop2[] = {0.00000, 1.19790, 0.67151, 9.79130, 0.00000, 0, 0};
	float dip[] = {9.32711, 0.57861, 25.21611, 0.00566, 0.00000, 0, 0};
	float dip2[] = {0.00003, 0.00799, 31.70140, 8.83931, 0.00000, 0, 0};
	float dip3[] = {2.21139, 205.58350, 216.35576, 3.84918, 0.00000, 0, 0};
	pval = gensnp.is_snp(monop,i,j,b);
	fprintf(stderr,"%f, %f, %f, %f, %f\n",monop[0],monop[1],monop[2],monop[3],monop[4]);
	fprintf(stderr,"PVAL: %.2e at %d and second is %d (should be 6.64e-04, -1)\n",pval,i,b?j:-1);
	TEST(!b);
	pval = gensnp.is_snp(dip,i,j,b);
	fprintf(stderr,"%f, %f, %f, %f, %f\n",dip[0],dip[1],dip[2],dip[3],dip[4]);
	fprintf(stderr,"PVAL: %.2e at %d and second is %d (should be 9.99e-16, 0)\n",pval,i,b?j:-1);
	TEST(b && j==0);
	pval = gensnp.is_snp(monop2,i,j,b);
	fprintf(stderr,"%f, %f, %f, %f, %f\n",monop2[0],monop2[1],monop2[2],monop2[3],monop2[4]);
	fprintf(stderr,"PVAL: %.2e at %d and second is %d (should be 6.59e-06, -1)\n",pval,i,b?j:-1);
	TEST(!b);
	pval = gensnp.is_snp(dip2,i,j,b);
	fprintf(stderr,"%f, %f, %f, %f, %f\n",dip2[0],dip2[1],dip2[2],dip2[3],dip2[4]);
	fprintf(stderr,"PVAL: %.2e at %d and second is %d (second should probably be 3)\n",pval,i,b?j:-1);	
	TEST_W(b && j==3);
	pval = gensnp.is_snp(dip3,i,j,b);
	fprintf(stderr,"%f, %f, %f, %f, %f\n",dip3[0],dip3[1],dip3[2],dip3[3],dip3[4]);
	fprintf(stderr,"PVAL: %.2e at %d and second is %d (second should probably be 1)\n",pval,i,b?j:-1);
	TEST_W(b && j==1);
	float noSNP[] = {0.07155, 194.00175, 201.30124, 11.62545, 0.00000};
	pval = gensnp.is_snp(noSNP,i,j,b);
	fprintf(stderr,"%f, %f, %f, %f, %f\n",noSNP[0],noSNP[1],noSNP[2],noSNP[3],noSNP[4]);
	fprintf(stderr,"PVAL: %.2e at %d and second is %d (second should probably be 1)\n",pval,i,b?j:-1);
	TEST_W(b && j==1);
	//return success;
	 
	unsigned char* genome = new unsigned char[1000010];
	ifstream in("test/short.fa");
	string name;
	getline(in,name);
	in.read((char*)genome,1000000);
	
	//&genome.front() = "actgctgccggtN\ncccgtcggcaacgtgNNNNNgaa\n\0";
	
	try {
		os << "**TESTING HASH_AND_STORE**\n";

		unsigned int prev_MER_SIZE = gMER_SIZE;
		gMER_SIZE = 10;
		gHASH_SIZE = (unsigned int) pow(4.0,gMER_SIZE);
		int prev_snp = gSNP;
		int prev_gensize = gGEN_SIZE;
		gSNP = 1;
		gGEN_SIZE = 1;
		GenomeBwt genYumei;
		char fn[] = "test/yumei_test.fa";
		genYumei.use(fn);
		genYumei.LoadGenome();
		os << "Beginning string: " << genYumei.GetString(0,1000) << endl;
		string str_str = genYumei.GetString(0,20);
		string char_str = "";
		for(unsigned int i=0; i<20; i++) {
			char_str += genYumei.GetChar(i);
		}
		TEST(str_str == char_str);

		genYumei.AddScore(0,1);
		float toAdd[5] = {0.0, 0.2, 0.2, 0.6, 0.0};
		float total[5] = {0.0, 0.2, 0.2, 0.6, 0.0};
		
		os << "Adding 1 with {0.0, 0.2, 0.2, 0.6, 0.0}" << endl;
		genYumei.AddSeqScore(0,toAdd,1);
#ifdef DISCRETIZE
		center_t *cent = centers[genYumei.GetGenomeAllotPtr()[0]];
		os << "\tallotcenter: " << (int)genYumei.GetGenomeAllotPtr()[0] << " vs " << find_center(total)
		   << " vs exhaustive: " << find_center_exhaustive(total) << "\n";
		os << "\tfloat representation is: "<<cent[0]<<","<<cent[1]<<","<<cent[2]<<","<<cent[3]<<","<<cent[4]<< "\n";
		os << "\t              should be: "<<total[0]<<","<<total[1]<<","<<total[2]<<","<<total[3]<<","<<total[4]<<"\n";
#elif defined( INTDISC )
		os << "a: " << (float)genYumei.GetGenomeAllotPtr()[0+A_POS]/255.0 << "\n";
		os << "c: " << (float)genYumei.GetGenomeAllotPtr()[0+C_POS]/255.0 << "\n";
		os << "g: " << (float)genYumei.GetGenomeAllotPtr()[0+G_POS]/255.0 << "\n";
		os << "t: " << (float)genYumei.GetGenomeAllotPtr()[0+T_POS]/255.0 << "\n";
		os << "n: " << (float)genYumei.GetGenomeAllotPtr()[0+N_POS]/255.0 << "\n";
#else
		os << "a: " << genYumei.GetGenomeAPtr()[0] << "\n";
		os << "c: " << genYumei.GetGenomeCPtr()[0] << "\n";
		os << "g: " << genYumei.GetGenomeGPtr()[0] << "\n";
		os << "t: " << genYumei.GetGenomeTPtr()[0] << "\n";
		os << "n: " << genYumei.GetGenomeNPtr()[0] << "\n";
#endif
		float toAdd2[5] = {0.4, 0.2, 0.1, 0.3, 0.0};
		total[0] = .2/1.5;
		total[1] = .3/1.5;
		total[2] = .25/1.5;
		total[3] = .75/1.5;
		total[4] = 0;
		os << "Adding 2 with .5* {0.4, 0.2, 0.1, 0.3, 0.0}" << endl;
		genYumei.AddScore(0,.5);
		genYumei.AddSeqScore(0,toAdd2,.5);

#if defined( DISCRETIZE )
		cent = centers[genYumei.GetGenomeAllotPtr()[0]];
		os << "\tallotcenter: " << (int)genYumei.GetGenomeAllotPtr()[0] << " vs " << find_center(total)
		   << " vs exhaustive: " << find_center_exhaustive(total) << "\n";
		os << "\tfloat representation is: "<<cent[0]<<","<<cent[1]<<","<<cent[2]<<","<<cent[3]<<","<<cent[4]<< "\n";
		os << "\t              should be: "<<total[0]<<","<<total[1]<<","<<total[2]<<","<<total[3]<<","<<total[4]<<"\n";
#elif defined( INTDISC )
		os << "a: " << (float)genYumei.GetGenomeAllotPtr()[0+A_POS]/255.0 << " vs " << total[A_POS] << "\n";
		os << "c: " << (float)genYumei.GetGenomeAllotPtr()[0+C_POS]/255.0 << " vs " << total[C_POS] << "\n";
		os << "g: " << (float)genYumei.GetGenomeAllotPtr()[0+G_POS]/255.0 << " vs " << total[G_POS] << "\n";
		os << "t: " << (float)genYumei.GetGenomeAllotPtr()[0+T_POS]/255.0 << " vs " << total[T_POS] << "\n";
		os << "n: " << (float)genYumei.GetGenomeAllotPtr()[0+N_POS]/255.0 << " vs " << total[N_POS] << "\n";
		// We'll test to make sure there's not less than a 1% error
		TEST_W((float)genYumei.GetGenomeAllotPtr()[0+A_POS]/255.0 - (0.0+.5*.4)/1.5 < 0.01);
		TEST_W((float)genYumei.GetGenomeAllotPtr()[0+C_POS]/255.0 - (0.2+.5*.2)/1.5 < 0.01);
		TEST_W((float)genYumei.GetGenomeAllotPtr()[0+G_POS]/255.0 - (0.2+.5*.1)/1.5 < 0.01);
		TEST_W((float)genYumei.GetGenomeAllotPtr()[0+T_POS]/255.0 - (0.6+.5*.3)/1.5 < 0.01);
		TEST_W((float)genYumei.GetGenomeAllotPtr()[0+N_POS]/255.0 - (0.0+.5*.0)/1.5 < 0.01);
#else
		os << "a: " << genYumei.GetGenomeAPtr()[0] << "\n";
		os << "c: " << genYumei.GetGenomeCPtr()[0] << "\n";
		os << "g: " << genYumei.GetGenomeGPtr()[0] << "\n";
		os << "t: " << genYumei.GetGenomeTPtr()[0] << "\n";
		os << "n: " << genYumei.GetGenomeNPtr()[0] << "\n";
#endif
		
		os << "Adding 3 with .7* {0.7, 0.0, 0.0, 0.0, 0.0}" << endl;
		genYumei.AddScore(0,.7);
		genYumei.AddSeqScore(0,.7,A_POS);
		
		total[0] = .9/2.2;
		total[1] = .3/2.2;
		total[2] = .25/2.2;
		total[3] = .75/2.2;
		total[4] = 0;

#if defined( DISCRETIZE )
		cent = centers[genYumei.GetGenomeAllotPtr()[0]];
		os << "\tallotcenter: " << (int)genYumei.GetGenomeAllotPtr()[0] << " vs " << find_center(total)
		   << " vs exhaustive: " << find_center_exhaustive(total) << "\n";
		os << "\tfloat representation is: "<<cent[0]<<","<<cent[1]<<","<<cent[2]<<","<<cent[3]<<","<<cent[4]<< "\n";
		os << "\t              should be: "<<total[0]<<","<<total[1]<<","<<total[2]<<","<<total[3]<<","<<total[4]<<"\n";
#elif defined( INTDISC )
		os << "a: " << (float)genYumei.GetGenomeAllotPtr()[0+A_POS]/255.0 << " vs " << (0.0+.5*.4+.7)/2.2 << "\n";
		os << "c: " << (float)genYumei.GetGenomeAllotPtr()[0+C_POS]/255.0 << " vs " << (0.2+.5*.2)/2.2 << "\n";
		os << "g: " << (float)genYumei.GetGenomeAllotPtr()[0+G_POS]/255.0 << " vs " << (0.2+.5*.1)/2.2 << "\n";
		os << "t: " << (float)((int)(genYumei.GetGenomeAllotPtr()[0+T_POS]))/255.0 << " vs " << (0.6+.5*.3)/2.2 << "\n";
		os << "n: " << (float)genYumei.GetGenomeAllotPtr()[0+N_POS]/255.0 << " vs " << (0.0+.5*.0)/2.2 << "\n";
		// We'll test to make sure there's not less than a 1% error
		TEST_W((float)genYumei.GetGenomeAllotPtr()[0+A_POS]/255.0 - (0.0+.5*.4+.7)/2.2 < 0.01);
		TEST_W((float)genYumei.GetGenomeAllotPtr()[0+C_POS]/255.0 - (0.2+.5*.2)/2.2 < 0.01);
		TEST_W((float)genYumei.GetGenomeAllotPtr()[0+G_POS]/255.0 - (0.2+.5*.1)/2.2 < 0.01);
		TEST_W((float)genYumei.GetGenomeAllotPtr()[0+T_POS]/255.0 - (0.6+.5*.3)/2.2 < 0.01);
		TEST_W((float)genYumei.GetGenomeAllotPtr()[0+N_POS]/255.0 - (0.0+.5*.0)/2.2 < 0.01);
#else
		os << "a: " << genYumei.GetGenomeAPtr()[0] << "\n";
		os << "c: " << genYumei.GetGenomeCPtr()[0] << "\n";
		os << "g: " << genYumei.GetGenomeGPtr()[0] << "\n";
		os << "t: " << genYumei.GetGenomeTPtr()[0] << "\n";
		os << "n: " << genYumei.GetGenomeNPtr()[0] << "\n";
#endif

		os << endl;

		gSNP = prev_snp;
		gGEN_SIZE = prev_gensize;
		gMER_SIZE = prev_MER_SIZE;
		gHASH_SIZE = (unsigned int) pow(4.0,gMER_SIZE);


		GenomeBwt gen45;
		gREAD = true;
		gREAD_FN = "/data/public/genomes/human/gnumap/gnumap_human_s2_h100k_m11";
		gen45.use("/data/public/genomes/human/gnumap/gnumap_human_s2_h100k_m11");
		gen45.LoadGenome();
		os << "Getting string at the end of chr4: " << 1877852952-200 << " " << gen45.GetPos(1877852952-200,true) << endl;
		os << gen45.GetString(1877852952-200,200) << endl;
		os << "Getting string at the end of chrX: " << 2866253375-200 << " " << gen45.GetPos(2866253375-200,true) << endl;
		os << gen45.GetString(2866253375-200,200) << endl;
		gen45.AddScore(2866253375,4.76);
		os << "Getting value at 2866253375: " << gen45.GetScore(2866253375) << endl;
		os << "                    Should be 4.76" << endl;
		gREAD = false;
		
	}
	catch(const char* erro) {
		os << "\t**Error! Did you include the proper test files?" << endl;
		TEST_W(0);
	}
	catch (...) {
		os << "Unknown Error!" << endl;
	}

	delete[] genome;
	return success;
}
