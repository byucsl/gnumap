/* Genome.cpp
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

#include "Genome.h"

/*
 * Constructor.  Creates new Genome with filename.
 * fn - the filename of the file to be parsed into the genome.
 */
Genome::Genome(const char* fn) :
	offset(0), bin_offset(0), packed_genome(NULL), amount_genome(NULL), gs_positions(NULL), genome_size(0) {

	Init(fn);
	reader = new Reader();
}

Genome::Genome() :
	offset(0), bin_offset(0), reader(NULL), packed_genome(NULL), amount_genome(NULL), gs_positions(NULL), genome_size(0) {

}

Genome::~Genome() {
	//fprintf(stderr,"[-/%d] Genome destructor...\n",iproc);
	delete reader;
	
	if(packed_genome)
		free(packed_genome);
	if(amount_genome)
		free(amount_genome);
#ifdef SET_POS
	if(gs_positions)
		free(gs_positions);
#endif
	
	if(gSNP)
#if defined( DISCRETIZE )
		delete[] read_allot;
#elif defined( INTDISC )
		// To maintain cache consistency, store them as reads[pos][char]
		delete[] reads;
		reads = 0;
#else
        for(unsigned int i=0; i<NUM_SNP_VALS; i++) {
            delete[] reads[i];
            reads[i] = 0;
        }
#endif // DISCRETIZE/INTDISC

    if(gBISULFITE) {
#if defined( DISCRETIZE )
		delete[] read_allot;
#elif defined( INTDISC )
		delete[] reads;
		reads = 0;
#else
        for(unsigned int i=0; i<NUM_SNP_VALS; i++) {
            delete[] reads[i];
            reads[i] = 0;
        }
#endif // DISCRETIZE/INTDISC
    }

    if(gATOG) {
    	// We can have the SNP option and the ATOG option at the same time
#if defined( DISCRETIZE )
		delete[] read_allot;
#elif defined( INTDISC )
		if(reads) {
			delete[] reads;
			reads = 0;
		}
#else
        for(unsigned int i=0; i<NUM_SNP_VALS; i++)
        	if(reads[i]) {
				delete[] reads[i];
				reads[i] = 0;
			}
#endif // DISCRETIZE/INTDISC
    }
	
} 

void Genome::Init(const char* fn) {

	// Don't need to init names if we're going to do so in the read file
	if(!gREAD) 	{
		if(names.size()) //if we've already set up the "names" vector, just return.
			return;

		char files[strlen(fn)];
		strcpy(files,fn);
		
		const char* delims = ", \"\n";
		char* token = strtok(files,delims);
		if(gVERBOSE > 1)
			cout << "Found:" << endl;

		while(token != NULL) {
			if(gVERBOSE)
				cout << "\t" << token << endl;
			string temp_str = token;
			
			// Get rid of spaces that make it impossible to handle
			ReplaceSpaceWithUnderscore(temp_str);
			
			names.push_back(pair<string,unsigned long>(temp_str,0));
			
			token = strtok(NULL,delims);
		}

		//Push a dummy string on the end so we can know when to end...
		names.push_back(pair<string,unsigned long>("end",0));
	}

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

void Genome::use(const char* fn) {
	use(fn,0,0);
}

void Genome::use(const char* fn, unsigned long s, unsigned long e) {
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

void Genome::readGenome(FILE* readF, bool ldebug) {
	int ret;

	// Read in the genome.
	if(ldebug) fprintf(stderr,">> Reading genome from file\n");
	unsigned long temp_genome_size = 0;

	// The size of the genome
	fread(&temp_genome_size, sizeof(unsigned long), 1, readF);
	if(ldebug) fprintf(stderr,"\tgenome_size: %lu\n",temp_genome_size);

	//genome = new GenomeLocation[temp_genome_size];
	//genome = (GenomeLocation*)malloc(sizeof(GenomeLocation) * temp_genome_size);
	packed_genome = (GEN_TYPE*)calloc(sizeof(GEN_TYPE), temp_genome_size/GEN_PACK);
	amount_genome = (float*)calloc(sizeof(float), temp_genome_size/gGEN_SIZE);
	memset(amount_genome, 0, sizeof(float) * temp_genome_size/gGEN_SIZE);
#ifdef SET_POS
	gs_positions = (char*)calloc(sizeof(char), temp_genome_size);
	memset(gs_positions, 0, sizeof(char)*temp_genome_size);
	gs_positions_size = temp_genome_size;
#endif

	if(ldebug) fprintf(stderr,">> Finished reading sizes, now reading actual genome\n");
	ret = fread(packed_genome,sizeof(GEN_TYPE),temp_genome_size/GEN_PACK,readF);
	/*
	for(unsigned long i=0; i<temp_genome_size; i++) {
		GEN_TYPE genome_char = 0;
		fread(&(genome_char),sizeof(GEN_TYPE),1,readF);
		//genome[i].packed_char = genome_char;
		//genome[i].amount = 0;
		packed_genome[i] = genome_char;
		amount_genome[i] = 0;
	}
	*/
	genome_size = temp_genome_size;

	if(ldebug) fprintf(stderr,">> Reading chromosome names from file...\n");
	int num_chrs = 0;
	fscanf(readF,"\n%d\n",&num_chrs);
	//fprintf(stderr,"\n%d\n",num_chrs);
	
	/* Read the chromosome names from a file */
	for(int i=0; i<num_chrs; i++) {
		char name[MAX_NAME_SZ];
		unsigned long chr_pos;
		fscanf(readF,"%s %ld\n",name,&chr_pos);
		//fprintf(stderr,"Name: %s Size: %ld\n",name,chr_pos);
		names.push_back(pair<string,unsigned long>(name,chr_pos));
	}

#ifdef DEBUG	
	for(unsigned int i=0; i<names.size(); i++) {
		fprintf(stderr,"names[%u] (%lu): %s\n",i,names[i].second,names[i].first.c_str());
	}
#endif

}


unsigned int Genome::readGen(char* fn) {
	unsigned int ret;
	bool ldebug = false;

	/**********************/
	// Open the file
	/**********************/
	FILE* readF;

	readF = fopen(fn,"r");
	
	if(!readF) {
		printf("Trying to open: >%s<\n",fn);
		perror("File wouldn't open");
		throw new Exception("Could not open Genome read file");
	}

	int version = 0;

	/**********************/
	// Check the version number
	/**********************/
	if(fscanf(readF,"Version 1.%d\n",&version) == 1) {
		if(ldebug) fprintf(stderr,">> VERSION 1.%d\n",version);
	}
	else {
		if(ldebug) fprintf(stderr,">> VERSION 1.0\n");
	}

	/**********************/
	// Check the number of machines is acceptable
	/**********************/
	if(ldebug) fprintf(stderr,">> Reading number of machines...");
	int temp_nproc,temp_iproc;
	ret = fscanf(readF,"%*s %d %d\n",&temp_iproc,&temp_nproc);
	if(ret != 2)	// ERROR
		assert(!"ERROR:  Improper formatting for saved genome on first line");

	// If we're all using the same genome, we don't need to ensure this property
	fprintf(stderr,"%d/%d\n",iproc,nproc);
	if(gMPI_LARGEMEM) {
		if(temp_nproc != nproc)
			assert(!"ERROR:  Invalid number of processors with MPI when reading genome");
		if(temp_iproc != iproc)
			assert(!"ERROR:  Invalid processor number with MPI when reading genome");
	}
	else {	// want to make sure there's only one node, otherwise we'll only be reading a piece of the genome
		if(temp_nproc != 1 && temp_nproc != 0)
			assert(!"ERROR:  Invalid number of processors with MPI when reading genome");
	}
	
	// Read the genome start and end
	ret = fscanf(readF,"%*s %lu %lu\n",&my_start,&my_end);
	if(ret != 2)
		assert(!"ERROR:  Could not determine genome start or end\n");
	if(ldebug) fprintf(stderr,"my_start: %lu my_end: %lu\n",my_start,my_end);
	// Read a bool for read_all
	int temp_read;
	ret = fscanf(readF,"%*s %d\n",&temp_read);
	read_all = temp_read;
	if(ret != 1)
		assert(!"ERROR: Could not identify read_all flag");
	if(ldebug) fprintf(stderr,">> read_all? %c %d\n",read_all ? 'Y':'N',read_all);

	// Read in the hash (refactored to allow for expansion)
	readHash(readF, ldebug);
		
	// Read in the genome (refactored to allow for expansion)
	readGenome(readF, ldebug);

	fclose(readF);

	return 0;
}

void Genome::make_extra_arrays() {
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

void Genome::LoadGenome() {
	if(gREAD) {
		char read_fn[strlen(gREAD_FN)+10];
		if(gMPI_LARGEMEM && nproc > 1) {
			sprintf(read_fn, "%s_%d", gREAD_FN, iproc);
		}
		else {
			sprintf(read_fn, "%s", gREAD_FN);
		}
		fprintf(stderr,"\tReading from file: %s and iproc:%d nproc:%d\n",read_fn, iproc, nproc);// read the file in from memory
		readGen(read_fn);
		fprintf(stderr,"\tFinished reading file %s.\n",read_fn);
	}
	else {
		include_hash = true;
		hash_and_store();
	}

		
#ifdef DEBUG
		for(unsigned int i=0; i<names.size(); i++) {
			printf("Name[%d]: %s\tBegin: %ld\n",i,names[i].first.c_str(),names[i].second);
		}
		cout.flush();
#endif
	//commmented by CJ on 04/23/2012 to separate a postprocess, sam2gmp
	if (gSAM2GMP)
		make_extra_arrays();
	
	if(gSAVE) {
		char save_fn[strlen(gSAVE_FN)+10];
		if(gMPI_LARGEMEM && nproc > 1)
			sprintf(save_fn, "%s_%d", gSAVE_FN, iproc);
		else
			sprintf(save_fn, "%s", gSAVE_FN);

		fprintf(stderr,"\tSaving to file: %s\n",save_fn);
		saveGen(save_fn);
	}
	
	// Some print-to-hash stuff
	/*
	char fn[20];
	strcpy(fn,"printed_hash");
#ifdef MPI_RUN
	int iproc = MPI::COMM_WORLD.Get_rank();
	sprintf(fn,"%s_%d","printed_hash",iproc);
#endif
	print_to_file(gen_hash,fn);
	//*/
}

void Genome::ReplaceSpaceWithUnderscore(string &str) {
	for(unsigned int i=0; i<str.length(); i++) {
		if(str[i] == ' ' || str[i] == 10 || str[i] == 11 || str[i] == 12 || str[i] == 13)
			str[i] = '_';
	}

}


#define READ_EXTRA	1024

/*
 * store() will store the genome into a smaller array of BIT64_t portions.
 * Each character is represented by 4 bits (because we require the 4 nucleotides
 * and an additional N character), and can be accessed directly.  
 *
 * One of the problems we have now is the last few characters in gen.  They will
 * appear the second time at the beginning of the genome yet again, throwing off the 
 * numbering system. 
 *
 * @par gen - the piece of the genome to store.
 * @par num - the number of characters to store.
 * @par num_bases - the base counter.  And tells when to start recording
 */
//void Genome::store(gvector<GenomeLocation> *g, unsigned char* gen, 
//						unsigned long int num, unsigned long &num_bases) {
void Genome::store(gvector<GEN_TYPE> *pg, unsigned char* gen, 
						unsigned long int num, unsigned long &num_bases) {
	// DONE!!
	if(num_bases > my_end+READ_EXTRA && !read_all)
		return;
		
	bool start_yet = num_bases >= my_start ? true : false;

	//fprintf(stderr,"READ_ALL: %d start_yet: %d\n",read_all,start_yet);

	unsigned long gen_pos=0;
//	unsigned long int num_bases=0;
	
	while(gen_pos<num) {
		GEN_TYPE temp = 0;
		
		if(!read_all && num_bases > my_end+READ_EXTRA)
			break;

		//if( read_all || (start_yet && (num_bases <= my_end+READ_EXTRA)) ) {
			if(bin_offset != 0) {
				//fprintf(stderr,"Is this the error?  Offset is %d and my_start is %ul\n",offset,my_start);
				//temp = (*g)[offset-my_start/gGEN_SIZE].packed_char;
				//temp = (*g)[offset].packed_char;
				temp = (*pg)[offset];
				pg->pop_back();
			}
		//}
			
		for(; (bin_offset<GEN_PACK) && gen_pos<num; gen_pos++) {	//can fit 16 chars into a BIT64_t 
			//int int_to_add = base2int(gen[gen_pos]);
			int int_to_add = g_gen_CONVERSION[(int)gen[gen_pos]];
			
			
			//If it's a valid char.
			if(int_to_add < 5) {
				num_bases++;
				if(num_bases > my_start)
					start_yet = true;
				if(!read_all && num_bases > my_end+READ_EXTRA) {
					// remember this location, then return
					//g->push_back(GenomeLocation(temp,0));
					pg->push_back(temp);
					return;	// DONE!
				}

				if(start_yet) {
					temp <<= 4;
					temp |= (int_to_add & 0xf);
					//fprintf(stderr,"Adding %c: %x(%d), now %0.8x\n",gen[gen_pos],int_to_add & 0xf,int_to_add,temp);
					++bin_offset;

				}
									
			}
			//If we come accross the flag for the end of the stream--should never happen...
			else if(int_to_add == 6) { 
				gen_pos = num;
				//cout << "Found end." << endl;
				break;
			}
			else if(int_to_add == 7) {	//come accross a new chromosome name.
				string temp_str = (char*)(gen+gen_pos);
				//while(gen[gen_pos] != '\n' && gen[gen_pos] != '\0' && gen_pos<num )
				while( gen[gen_pos] != '\n' && gen[gen_pos] != '\0' )
					gen_pos++;	//increment the counter until we've found the end of the name
			}
			else ;//do nothing--it's an invalid character and not an end-of-line
		}
				
		if(bin_offset != 0)	{ //if we've made any changes to this value...
			try {
				//g->push_back(pair<GEN_TYPE,float>(temp,0));
				if(read_all || (start_yet && num_bases <= my_end+READ_EXTRA)) {
					//fprintf(stderr,"Pushing %0.8x onto genome\n",temp);
					//g->push_back(GenomeLocation(temp,0));
					pg->push_back(temp);
				}
			}
			catch(bad_alloc& e) {
				printf("Found ERROR!!!:  %lu elements out of max of %lu\n",pg->size(),pg->max_size());
				//printf("Total space is %lu\n",g->size()*sizeof(GenomeLocation));
				printf("Total space is %lu\n",pg->size()*sizeof(GEN_TYPE));
				throw e;
			}
		}
		else ;	//no changes.  Coming out of the previous for-loop
		
		
		//means we've filled it fully--so increment the offset.
		if(bin_offset == GEN_PACK) {
			if(start_yet) {
				offset++;
				bin_offset = 0;
			}
		}


	}
	
	//fprintf(stdout,"Just read %ld bases out of %lu assigned\n",num_bases, num);
	return;
}

void Genome::StoreGenome(bool make_extra) {
	// We'll set it in here so we don't need to take care of it elsewhere
	gMER_SIZE = MAX_MER_SIZE;
	
	// Don't hash anything--not needed
	include_hash = false;
	hash_and_store();
	
	// Still need to create extra arrays.
	if(make_extra) make_extra_arrays();
}

/*
 * This function will just count the genome--it won't store anything.
 */
unsigned long Genome::count() {
	
	// Set some values that are necessary for the 'store()' function
	offset = 0;
	bin_offset = 0;

	Reader reader_vec;
	unsigned char gen_piece[gBUFFER_SIZE+1];

	if(gVERBOSE)
		printf("Counting Genome to Reserve Size...\n");

	
	// If we've already calculated the genome size (read from a file), just return the size
	if(genome_size != 0) {
		fprintf(stderr,"Already determined genome size (probably because of a --read_file argument\n");
		return genome_size*(nproc > 1 ? nproc : 1);
	}

	bin_seq bs;	

	unsigned long counter = 0;
	unsigned long int chr_begin_pos = 0;
	unsigned int j;

	for(j=0; j<names.size()-1; j++) {
		reader_vec.use(names[j].first);
		string chr_name = reader_vec.GetName();
		ReplaceSpaceWithUnderscore(chr_name);
		
		//unsigned long chr_begin_pos = counter;
		chr_begin_pos = counter;
						
		bool keep_reading = true;
		unsigned int i;

		while(keep_reading) { 
	
			keep_reading = reader_vec.read(gen_piece);
	
			for(i=0; i < gBUFFER_SIZE-gMER_SIZE; i++) {
					
				pair<int ,unsigned int> hash = bs.get_hash(i,gen_piece,gBUFFER_SIZE-gMER_SIZE);
						
				//if hash.first is -2, it means it's a new 'chromosome'.
				if(hash.first == -2) {
					//printf("New Genome at counter: %lu and i: %u\n",counter,i);
					
					//Find the name of this genome.
					//at position i, we should have the '>'
					char gen_name[MAX_NAME_LEN+5];	//the name can only be 100 chars long
					unsigned int gen_pos=0; 
					
					//we want to read the name of the next chromosome that's been found
					while(gen_piece[++i] != '\n' && i < gBUFFER_SIZE) { 
						if(gen_pos < MAX_NAME_LEN)
							gen_name[gen_pos] = gen_piece[i];
						gen_pos++;
					}
					//Three things could happen here: 
					// 1. The end of the name (\n) extends beyone the end of the buffer
					//		(it didn't finish reading the name)
					// 2. The name extends into the last MER_SIZE bases, but doesn't 
					//		continue past the end of the buffer
					
					//*** CASE 1 ***//
					if(i >= gBUFFER_SIZE) {	//it didn't finish reading the name
						//fprintf(stderr,"Special Case #1(%u)\n",i);
						//fprintf(stdout,"Special Case #1\n");
						char* line = reader_vec.ReadLine();
						string temp = line;
						delete[] line;
						
						// Set up the 'chromosome' name
						unsigned int chr_name_pos = 0;
						while(gen_pos < MAX_NAME_LEN && chr_name_pos < temp.size())	//don't overflow the gen_name
							gen_name[gen_pos++] = temp[chr_name_pos++];
						gen_name[gen_pos] = '\0';
					}
					
					//*** CASE 2 ***//
					// Here, there's not enough characters on the end of the sequence to read
					// a new hash, so it throws off the numbering and reading(?).
					// By default, i < gBUFFER_SIZE
					else if(i + gMER_SIZE >= gBUFFER_SIZE) {
						// Shift the reader's offset by the number of bases remaining in the sequence
						// so it will read them again.
						reader_vec.ShiftOffset(gMER_SIZE - (gBUFFER_SIZE - i));
						
					}

					gen_name[gen_pos > MAX_NAME_LEN ? MAX_NAME_LEN : gen_pos] = '\0'; 
					
					string str_gen_name = gen_name;
					ReplaceSpaceWithUnderscore(str_gen_name);

					vector<pair<string,unsigned long> >::iterator vit = names.begin();
					//Find the position of this genome name
					for(unsigned int name_pos=0; name_pos < j; name_pos++)
						vit++;
					vit++;	// the vector inserts to the one "previous", but we want to insert
							// after the current name.

					names.insert(vit,pair<string,unsigned long>(str_gen_name,0));

					//take care of the former name.
					//1. set the name in the genome to have the correct position
					//2. set the chr_begin_pos var to be the number of chars we've stored so far.
					names[j] = pair<string,unsigned long>(chr_name,chr_begin_pos);
					chr_name = str_gen_name;
				
					counter += hash.second;	//increment the counter the number of bases we've skipped

					chr_begin_pos = counter;
					j++;
				}
				

				//if hash.first is -1, it means it's found the end-of-line char.
				else if(hash.first == -1 || i >= gBUFFER_SIZE-gMER_SIZE) {
					break;
				}
				else {	// can't be -2 or -1, so it must be nonnegative.
					counter++;
				}
			}
		}

		names[j] = pair<string,unsigned long>(chr_name,chr_begin_pos);
		
		//there are cases where the end of the line throws off our numbering
		//We'll fix it here.
		unsigned long int end = i;
		int count=0;
		for(unsigned int p=end; gen_piece[p] != '\0'; p++) {
			//if(base2int(gen_piece[p]) <= 4) {
			if(g_gen_CONVERSION[(int)gen_piece[p]] <= 4) {
				counter++;
				count++;
			}
		}
	}
	names[j] = pair<string,unsigned long>("end",counter);
	
//	fprintf(stderr,"Size of genome: %lu vs %lu\n",counter,chr_begin_pos);

		
#ifdef DEBUG 
		for(unsigned int i=0; i<names.size(); i++) {
			printf("Name[%d]: %s\tBegin: %lu\n",i,names[i].first.c_str(),names[i].second);
		}
		fflush(stdout);
		printf("Last Hash number: %lu\n",counter);
#endif

	if(gVERBOSE)
		printf("...Finished!\n\n");
	
	return counter;

}


void Genome::saveGenome(FILE* saveF, bool ldebug) {
	/***********************************/
	/* Write the genome out to a file. */
	/***********************************/
	if(ldebug) fprintf(stderr,">> Writing genome out to file...\n");
	if(ldebug) fprintf(stderr,"\tgenome_size: %lu\n",genome_size);
	unsigned long tempgenome_size = genome_size;
	fwrite(&tempgenome_size, sizeof(unsigned long), 1, saveF);
	if(ldebug) fprintf(stderr,">> Finished writing sizes, now writing actual genome\n");
	fwrite(packed_genome,sizeof(GEN_TYPE),genome_size/GEN_PACK,saveF);

	if(ldebug) fprintf(stderr,">> Writing chromosome names to file...\n");
	fprintf(saveF,"\n%d\n",(int)names.size());
	
	/* Write the chromosome names to a file */
	if(ldebug) fprintf(stderr,">> Writing actual chromosomes to file...\n");
	for(unsigned int i=0; i<names.size(); i++) {
		fprintf(saveF,"%s %ld\n",names[i].first.c_str(),names[i].second);
	}
}

unsigned int Genome::saveGen(char* fn) {
	/**********************/
	// Open the file
	/**********************/
	FILE* saveF;
	bool ldebug = true;

	saveF = fopen(fn,"w");
	if(!saveF)
		throw "Could not open Genome save file";

	/**********************/
	// Write the version number
	/**********************/
	fprintf(saveF,"Version 1.1\n");

	/**********************/
	// Write the machine information
	/**********************/
	if(ldebug)  {
		fprintf(stderr,">> Writing machine information:\n");
		fprintf(stderr,"\t%s %d %d\n","iproc,nproc: ",iproc,nproc);
		fprintf(stderr,"\t%s %lu %lu\n","start,end: ",my_start,my_end);
		fprintf(stderr,"\t%s %d\n","read_all?",read_all?1:0);
	}
	if(gMPI_LARGEMEM)
		fprintf(saveF,"%s %d %d\n","nproc,iproc:",iproc,nproc);
	else	// If we're not usig the MPI_LARGEMEM option, we're always writing with one node
		fprintf(saveF,"%s %d %d\n","nproc,iproc:",1,0);
	fprintf(saveF,"%s %lu %lu\n","start,end:",my_start,my_end);
	fprintf(saveF,"%s %d\n","read_all:",read_all?1:0);

	// Write the hash to a file
	saveHash(saveF, ldebug);

	// Write the genome to a file
	saveGenome(saveF, ldebug);

	// Close the file
	fclose(saveF);

	return 0;
}

/*
 * Returns a string from begin to begin+size.
 * the first position in the genome is position 0.
 * @par begin - the beginning of the string to return
 * @par size - the size of the string to return.
 */
string Genome::GetString(const unsigned long begin, const unsigned int size) { 
	return GetString(begin, size, begin, false);
}
string Genome::GetString(const unsigned long begin, const unsigned int s, 
		const unsigned int hash_start, bool extend_past) { 
#ifdef DEBUG
	fprintf(stderr,"[%d] Getting string at %lu",iproc,begin);
#endif
	unsigned int size = s;
 
	string to_return = "";

	/***********************Check for Bounds******************************/
	//We don't want to get a string beyond the end of the chromosome...
	int chr_num=0;
	//											it's already been pre-incremented
	while(names[chr_num].second <= begin && ((unsigned int)chr_num <= names.size()) ) chr_num++;
	//while(names[++chr_num].second <= begin && ((unsigned int)chr_num <= names.size()) );
	
	//means it's not valid--will extend beyond the chromosome file
	if( (names[chr_num].second < begin+size) )  {
		if(extend_past)
			size = names[chr_num].second - begin;
		else
			//printf("Extending Beyond the Chromosome\n");
			return "";
	}
	
	// Make sure we're not trying to read before the genome began in a multi-part run
	if(my_start > begin) {
		return "";
	}
	/*********************************************************************/
	
	unsigned long pos_in_gen = (begin-my_start)/GEN_PACK; 
	unsigned int temp_size = size;

	// Check over-bounds error
	//fprintf(stderr,"[%d] with begin=%lu, Assertion checks %lu <= %lu\n",
	//		iproc,begin,pos_in_gen,genome_size/GEN_PACK);
	assert(pos_in_gen <= genome_size/GEN_PACK);

	//GEN_TYPE t_temp = genome[pos_in_gen].packed_char;
	GEN_TYPE t_temp = packed_genome[pos_in_gen];

	//cerr << "\tbegin: " << begin << endl;
	//cerr << "\tmy_start: " << my_start << endl;
	//cerr << "\tbegin-my_start: " << begin-my_start << endl;
	//cerr << "\tpos_in_gen: " << pos_in_gen << endl;
	//cerr << "\tt_temp: " << t_temp << endl;
		
	//if(gVERBOSE > 3)
		//fprintf(stderr,"[%d/-] Getting position: %lu, max size is %lu\n",iproc,pos_in_gen,genome_size/GEN_PACK);
		//cerr << "getting position: " << t_temp << " max size is: " << genome_size/GEN_PACK << endl;
	
	unsigned int count=0;

	while(count < size) {
		//fprintf(stderr,"\tt_temp at %lu/%lu: %.8x %.8x\n",pos_in_gen+1,genome_size/GEN_PACK,t_temp,packed_genome[pos_in_gen+1]);

		string temp = "";
	
		if(count == 0) {	//get the last few characters off the BIT64_t
			unsigned long bit_offset = ((begin-my_start)%GEN_PACK);

			for(unsigned int j=GEN_PACK-bit_offset; j>0; j--) {
				//shift 0xf left 4 bases for each j.
				//then we don't need to reverse the string.
				//printf("Shifting by: %d = %lld\n",4*(j-1),(t_temp>>(4*(j-1))) & 0xf);
				temp += gINT2BASE[0xf & (t_temp >> (4*(j-1)))];
#ifdef DEBUG
                fprintf(stderr,"%d",0xf & (t_temp >> (4*(j-1))));
#endif
				//t_temp >>= 4;
				//fprintf(stderr,"Round %d, temp is %s\n",count,temp.c_str());
				count++;
			}
			//temp = reverse(temp);
			if(temp.size() > size)
				return temp.substr(0,size);
		}
		else {	//we now want the first offset characters off the BIT64_t.
			for(unsigned int j=GEN_PACK; j>0; j--) {
				
				//fprintf(stderr,"Shifting by: %d = %lld\n",4*(j-1),(t_temp>>(4*(j-1))) & 0xf);
				//temp += int2base(t_temp & 0xf);
				//temp += int2base((t_temp >> (4*(j-1))) & 0xf);
				temp += gINT2BASE[0xf & (t_temp >> (4*(j-1)))];
				//t_temp >>= 4;
#ifdef DEBUG
                fprintf(stderr,"%d",0xf & (t_temp >> (4*(j-1))));
#endif
			}
			//temp = reverse(temp);
			temp = temp.substr(0,temp_size);
			count += temp.size();
			//fprintf(stderr,"Round %d, temp is %s and return is %s\n",count,temp.c_str(),to_return.c_str());
		}
		to_return += temp;
		
		temp_size = size - count;
		
		++pos_in_gen;
		if((pos_in_gen)  >= genome_size/GEN_PACK) {
				t_temp = (GEN_TYPE)0x44444444;
		}
		else
			//t_temp = genome[pos_in_gen].packed_char;
			t_temp = packed_genome[pos_in_gen];
	}
#ifdef DEBUG
    fprintf(stderr,"\n");
#endif
	
	return to_return;
}

char Genome::GetChar(const unsigned long pos) {

	unsigned long pos_in_gen = (pos-my_start)/GEN_PACK;
	assert(pos_in_gen <= genome_size/GEN_PACK);

	GEN_TYPE t_temp = packed_genome[pos_in_gen];
	unsigned long bit_offset = (pos-my_start)%GEN_PACK;
	return gINT2BASE[0xf & (t_temp >> (4*(GEN_PACK-bit_offset-1)))];
}

float Genome::GetScore(const unsigned long &pos) {

	unsigned long pos_temp = (pos-my_start)/gGEN_SIZE;
	//return genome[pos_temp].amount;
	return amount_genome[pos_temp];
}

void Genome::AddScore(const unsigned long &pos, const float &amt) {
//	if(amt > 1)
//		fprintf(stderr,"Genome Error!! %d\n",__LINE__);
		
	unsigned long pos_temp = pos-my_start;
	//fprintf(stderr,"pos: %lu, myStart: %lu pos_temp: %lu\n",pos,my_start,pos_temp);
	//if(pos_temp > genome_size)
	//	fprintf(stderr,"Adding %f to %lu with max of %lu\n",amt,pos_temp/gGEN_SIZE,genome_size);
	amount_genome[pos_temp/gGEN_SIZE] += amt;
}

/**
 * INVARIANT:  must add to the genome previous to this call
 *			thus, amount_genome[pos] will never be 0
 */
void Genome::AddSeqScore(unsigned long pos, const float* amt, const float scale) {
	//fprintf(stderr,"float* - pos: %lu, myStart: %lu\n",pos,my_start);
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
#ifdef FIND_EXHAUSTIVE
	read_allot[loc] = find_center_exhaustive(gen_amt);
#else
	read_allot[loc] = find_center(gen_amt);
#endif

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
void Genome::AddSeqScore(unsigned long pos, const float amt, unsigned int which) {
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
	if(which < NUM_CHAR_VALS && reads[which])
		reads[which][loc] += amt;
#endif	// DISCRETIZE/INTDISC

}

string Genome::GetPos(const pair<unsigned long,int> p_pos) {
	return GetPos(p_pos.first,p_pos.second);
}
string Genome::GetPos(const unsigned long &pos, const int strand) {
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

pair<string,unsigned long> Genome::GetPosPair(const unsigned long pos) {
	
	int chr_num = 0;

	while(names[chr_num].second <= pos && ((unsigned int)chr_num < names.size()) )
		chr_num++;
		
	chr_num--;

	string chr_name = names[chr_num].first;
	unsigned long chr_pos = (pos - names[chr_num].second);
	
	return pair<string,unsigned long>(chr_name,chr_pos);
}

unsigned long Genome::GetAbsolutePosition(const char* chr, const unsigned long &pos) {
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

string Genome::GetStringRelative(const char* chr, const unsigned long &pos, const unsigned int length) {
	unsigned long abs_pos = GetAbsolutePosition(chr,pos);
	return GetString(abs_pos,length);
}


#ifdef SET_POS
bool Genome::getSetThreadID(long unsigned int pos, unsigned int tid) {
	if(pos < my_start)
		return false;

	//fprintf(stderr,"pos=%lu, genome_size=%lu, my_start=%lu, final_pos=%lu\n",pos,genome_size,my_start,pos-my_start);
	pos -= my_start;

	assert(pos < genome_size);

	if(gs_positions[pos/GEN_PACK] & (1<<(tid%8)))
		return true;

	unsigned int new_tid = 1<<(tid%8);
	gs_positions[pos/GEN_PACK] |= new_tid;
	return false;
}

void Genome::unsetThreadID(long unsigned int pos, unsigned int tid) {
	if(pos < my_start)
		return;

	pos -= my_start;


	assert(pos < genome_size);
	// We'll clear all the threads at this position
	gs_positions[pos/GEN_PACK] = 0;
}
#endif

unsigned int Genome::max_pos(float array[], unsigned int array_size) {
	float* pos = max_element(array, array+array_size); \
	unsigned int max_pos = pos - &array[0]; \
	return max_pos; \
}

const unsigned int MIN_COVERAGE = 5;
const float MIN_RATIO = 1.5;
const float MAX_PVAL = 0.01;
const float MAX_RATIO = 0.0;

double Genome::LRT(float chars[5], int &snp_pos) {
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

double Genome::dipLRT(float chars[NUM_CHAR_VALS], int &snp_pos1, int &snp_pos2, bool &dip) {

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

double Genome::is_snp(float chars[NUM_CHAR_VALS], int &snp_pos1, int &snp_pos2, bool &dip) {
	
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
//inline int Genome::is_indel

void Genome::AppendFinal(const char* fn) {
	PrintFinal(fn, true);
}
void Genome::PrintFinal(const char* fn) {
	PrintFinal(fn, false);
}

void Genome::PrintFinal(const char* fn, bool append) {
	if(gPRINT_VCF) {
		// Not sure how to handle VCF and Bisulfite.  We'll let them deal with it.
		PrintFinalVCF(fn, append);
		return;
	}

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

void Genome::PrintFinalSNP(const char* fn, bool append) {
	// Instead of printing out a file for each, we want to print
	// out one file (dubbed the '.sgrex' file for 'sgr-extra') to be a tab-deliminated file
	// for the number of a's, c's, g's, t's, and gaps at each location in the genome.
	char gmp_fn[strlen(fn)+10];
	strcpy(gmp_fn,fn);
	strcat(gmp_fn,".gmp");
	fprintf(stderr, "SNP File is: %s\n",gmp_fn);
	
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

void Genome::PrintSNPCall(const unsigned long count, FILE* gmp_file) {
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

void Genome::PrintFinalBisulfite(const char* fn, bool append) {
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
		fprintf(stderr,"Size of NAMES: %ld, size of gGEN_SIZE:%d\n",names.size(),gGEN_SIZE);
	
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

void Genome::PrintFinalVCF(const char* fn, bool append) {
	char* vcf_fn = new char[strlen(fn)+5];
	strcpy(vcf_fn, fn);
	strcat(vcf_fn, ".vcf");
	fprintf(stderr, "VCF file is %s\n", vcf_fn);
	
	FILE* vcf_file;
	if(append) {
		vcf_file = fopen(vcf_fn, "a");
	}
	else {
		// If we're not appending, we need to write the header lines
		vcf_file = fopen(vcf_fn, "w");
		fprintf(stderr, "Just opened vcf file %s.  Successful?%c\n", vcf_fn, vcf_file?'Y':'N');
		fprintf(vcf_file, "##fileformat=VCFv4.0\n");
		time_t rawtime;
		time ( &rawtime );
		fprintf(vcf_file, "##fileDate=%s\n", ctime (&rawtime));
		fprintf(vcf_file, "##source=GNUMAP v"gVERSION"\n");
		fprintf(vcf_file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
		fprintf(stderr, "Finished writing header to VCF file\n");
	}

	if(!vcf_file)
		perror("Could not open .vcf file!");

	unsigned long snp_count=0;
	unsigned long count=0;
	for(unsigned int j=0; j<names.size(); j++) {
		fprintf(stderr, "Printing out for name %s (my_start=%lu,my_end=%lu,read_all=%c)\n", names[j].first.c_str(), my_start,my_end,read_all?'y':'n');
		for(; count<names[j+1].second; count+=gGEN_SIZE) {
			if(count < my_start && !read_all)
				continue;
			if(count > my_end && !read_all) {
				fprintf(stderr, "\t...skipping\n");
				break;
			}

			unsigned long locus = (count-my_start)/gGEN_SIZE;
			if(amount_genome[locus] <= 0.0f)
				continue;

			char at_pos = GetChar(count);
#if defined( DISCRETIZE )
			center_t *map_c = centers[read_allot[locus]];
			float char_arr[DIM_CENTERS];
			for(unsigned int i=0; i<DIM_CENTERS; i++)
				char_arr[i] = (float)map_c[i]*amount_genome[locus];
#elif defined( INTDISC )
			
			float char_arr[NUM_SNP_VALS];
			for(unsigned int i=0; i<NUM_SNP_VALS; i++) {
				char_arr[i] = (float)reads[locus*NUM_SNP_VALS+i]/255
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

			int snp_pos1, snp_pos2;
			bool dip = false;
			double p_val = is_snp(char_arr, snp_pos1, snp_pos2, dip);


			// If the max pos is not what's in the genome
			// 	or it's a diploid snp
			if(snp_pos1 != g_gen_CONVERSION[(int)at_pos] || dip) {
				if( p_val < gSNP_PVAL) {
					if(dip) {
						fprintf(vcf_file,"%s\t%u\tsnp%u\t%c\t%c%c\t.\t.\tDiploid;pval=%.5f;coverage=%.5f;ratio=%.2f\n",
								names[j].first.c_str(), (count-names[j].second)+1, snp_count++, 
								at_pos, gINT2BASE[snp_pos1],gINT2BASE[snp_pos2], 
								p_val, amount_genome[locus], char_arr[snp_pos2]/char_arr[snp_pos1]);
					}
					else {
						fprintf(vcf_file,"%s\t%u\tsnp%u\t%c\t%c\t.\t.\tMonoploid;pval=%.5f;coverage=%.5f\n",
								names[j].first.c_str(), (count-names[j].second)+1, snp_count++, 
								at_pos, gINT2BASE[snp_pos1],
								p_val, amount_genome[locus]);
					}
				}
			}


		}
	}
	
}

void Genome::PrintFinalSGR(const char* fn, bool append) {
	char* filename = new char[strlen(fn)+5];
	sprintf(filename,"%s.sgr",fn);
	
	FILE* sgr_file;
	if(append)
		sgr_file = fopen(filename,"a");
	else
		sgr_file = fopen(filename,"w");
	
	//fprintf(stderr,"Printing for %u genomes to %lu\n",names.size(),names[names.size()-1].second);
	unsigned long count=0;
	for(unsigned int i=0; i<names.size()-1; i++) {
		for(; count<names[i+1].second; count+=gGEN_SIZE) {
			if(count < my_start && !read_all)
				continue;
			if(count > my_end && !read_all)
				break;

			//if(genome[(count-my_start)/gGEN_SIZE].amount> 0)
			if(amount_genome[(count-my_start)/gGEN_SIZE] > MIN_PRINT) {
				//fprintf(stderr,"count=%u, start=%lu, end=%lu/%u=%lu, size=%lu, at: %f\n",count,my_start,my_end,gGEN_SIZE,my_end/gGEN_SIZE,genome_size,amount_genome[(count-my_start)/gGEN_SIZE]);
				
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

unsigned long Genome::size() {
	return genome_size;
}

