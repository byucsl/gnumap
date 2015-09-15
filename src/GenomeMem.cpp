/* GenomeMem.cpp
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

#include "GenomeMem.h"

/*
 * Constructor.  Creates new Genome with filename.
 * fn - the filename of the file to be parsed into the genome.
 */
GenomeMem::GenomeMem(const char* fn) : Genome(fn) {
}

GenomeMem::GenomeMem() : Genome() {

}

GenomeMem::~GenomeMem() {
    if(gen_hash.has_data_ptr())
    	gen_hash.free_data_ptr();
} 


void GenomeMem::print_to_file(hash_map<unsigned int, HashLocation> & gh, char* ofn) { 
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

void GenomeMem::readHash(FILE* readF, bool ldebug) {
	if(ldebug) fprintf(stderr,">> Reading hash from file...\n");
	unsigned long temp_gen_hash_size;
	//fscanf(readF,"%ld\n",&gen_hash_size);
	fread(&temp_gen_hash_size, sizeof(unsigned long), 1, readF);
	gHASH_SIZE = temp_gen_hash_size;
	// fprintf(stderr,"gen_hash_size: %d\n",gHASH_SIZE);
	gen_hash.set_size(gHASH_SIZE);

	// Read in the MER_SIZE (this value will change) 
	fread(&gMER_SIZE, sizeof(int), 1, readF);
	// gMER_SIZE = 10;
	if(ldebug) fprintf(stderr,"\tMer Size: %d\n",gMER_SIZE);

	// Read in the number of hash elements (so we can allocate the array once)
	fread(&num_hash_elements, sizeof(unsigned long), 1, readF);
	if(ldebug) fprintf(stderr,"\tNum Hash Elements: %lu\n",num_hash_elements);
	unsigned long* location_array = 
			(unsigned long*)calloc(sizeof(unsigned long), num_hash_elements);
	assert(location_array);
	// Read the entire array from memory
	fread(location_array, sizeof(unsigned long), num_hash_elements, readF);

	// Save the data pointer
	gen_hash.set_data_ptr((void*)location_array);
	
	// Read in the hash
	unsigned long count = 0;
	while(true) {
		unsigned int index,hash_size;
		//fscanf(readF,"%d %d\n",&index,&hash_size);
		fread(&index, sizeof(unsigned int), 1, readF);
		fread(&hash_size, sizeof(unsigned int), 1, readF);
		//fprintf(stderr,"%d %d\n",index,hash_size);
		if(!hash_size) {	//This is a marker for the end of the hash.
			if(ldebug) fprintf(stderr,">> Reading end marker\n");
			break;
		}

		gen_hash[index].size = hash_size;
		// malloc a new array for this location
		//unsigned long* temp_arr = new unsigned long[hash_size];
		gen_hash[index].hash_arr = location_array+count;
		count += hash_size;
		/*
		for(unsigned int j=0; j<hash_size; j++) {
			unsigned long hash_spot;
			fread(&hash_spot,sizeof(unsigned long),1,readF);
			temp_arr[j] = hash_spot;
			//gen_hash[index].push_back(hash_spot);
		}
		*/
	}

	// fprintf(stderr,"Finished reading in hash.  Reading in the genome now\n");
}


void GenomeMem::hash_and_store() {
	
	// Records what position we want to start saving at
	bool start_yet = my_start == 0 ? true : false;
	
	// Set some values that are necessary for the 'store()' function
	offset = 0;
	bin_offset = 0;

	hash_map<unsigned int, vector<unsigned long> > * gen_hash_vec
		= new hash_map<unsigned int, vector<unsigned long> >;
	//gvector<GenomeLocation> * genome_vec
	//	= new gvector<GenomeLocation>;
	gvector<GEN_TYPE> *genome_vec
		= new gvector<GEN_TYPE>;

	Reader reader_vec;
	unsigned char gen_piece[gBUFFER_SIZE+1];

	printf("Hashing Genome...\n");
	if(include_hash)
		gen_hash_vec->set_size(gHASH_SIZE);
	
	bin_seq bs;	

	// should be 1-based offset, so we'll start the counter at 1 as well.
	//unsigned long counter = 0;
	unsigned long counter = 0;
	unsigned long int chars_stored = 0;
	unsigned long int chr_begin_pos = 0;
	unsigned int j;

	for(j=0; j<names.size()-1; j++) {
		reader_vec.use(names[j].first);
		string chr_name = reader_vec.GetName();
		//fprintf(stderr,"Name here is: %s\n",chr_name.c_str());
		ReplaceSpaceWithUnderscore(chr_name);
		
		//unsigned long chr_begin_pos = counter;
		chr_begin_pos = chars_stored;
						
		unsigned long file_size = reader_vec.GetFileSize();
		//cerr << "File Size: " << file_size << endl;
		bool keep_reading = true;
		unsigned int i;

		while(keep_reading) { 
	
			keep_reading = reader_vec.read(gen_piece);
	
			//If we're not at the end of the file, then we can throw it the whole chunk,
			//without the extra gMER_SIZE bases at the end.
			//When we come to the end of the file, we can give it everything.
			//We'll also return the number of characters stored.
			//fprintf(stderr,"Buffer size is: %u and mer size is: %u\n",gBUFFER_SIZE,gMER_SIZE);
			if(chars_stored != counter) {
				if(gVERBOSE)
					fprintf(stderr,"Possible Error at %s::%d: NE: chars_stored: %lu counter: %lu\n",
							__FILE__,__LINE__,chars_stored,counter);
					
				chars_stored = counter;
			}
			//else {
			//	fprintf(stdout,"EQ: chars_stored: %u counter: %u\n",chars_stored,counter);
			//}
			
			if(keep_reading) {
				try {
					//fprintf(stderr, "Keep trying...\n");
					// We'll have it increment our counter for us
					store(genome_vec,gen_piece,gBUFFER_SIZE-gMER_SIZE,chars_stored);
				}
				catch(bad_alloc& e) {
					printf("Caught bad alloc around store!\n");
					throw e;
				}
				//printf("keeping reading\n");
			}
			else {
				try {
					fprintf(stderr, "Not keep trying\n");
					store(genome_vec,gen_piece,gBUFFER_SIZE,chars_stored);
				}
				catch(bad_alloc& e) {
					printf("Caught bad alloc around store!\n");
					throw e;
				}
			}
			
	 
			for(i=0; i < gBUFFER_SIZE-gMER_SIZE; i++) {
					
				pair<int ,unsigned int> hash = bs.get_hash(i,gen_piece,gBUFFER_SIZE-gMER_SIZE);
				//fprintf(stderr, "hash elem: %u\n", hash.first);

				//if hash.first is -1, it means it's found the end-of-line char.
				if(hash.first == -1 || i >= gBUFFER_SIZE-gMER_SIZE) {
					//if(hash.first == -1)
						//printf("Found the eof char\n");
					break;
				}

				else if(hash.first >= 0) {	// can't be -2 or -1, so it must be nonnegative.
					if(hash.first != 0)	{	//don't want to store N's.
						if(include_hash && (counter % gGEN_SKIP == 0)) { //This allows us a smaller hash size in total.
							try {
								if(start_yet) {
								//if(start_yet && 
								//		(gMAX_HASH_SIZE > 0 &&(*gen_hash_vec)[hash.second].size() > gMAX_HASH_SIZE))
									(*gen_hash_vec)[hash.second].push_back(counter);
								}
							}
							catch(bad_alloc& e) {
								printf("Caught bad alloc around push_back!\n");
								throw e;
							}
							//fprintf(stderr,"Hash of %u (%d) is %s\n",hash.second,hash.first,bs.hash2str(hash.second,gMER_SIZE).c_str());
						}
						
						// print output
						if( gVERBOSE && ((counter % (file_size/50) ) == 0)) {
							cout << ".";
							cout.flush();
						}

					}

					counter++;

					// Check to see if we've come across the starting point
					if(counter >= my_start)
						start_yet = true;
					if(counter > my_end && !read_all) {
						//fprintf(stderr,"Exiting HERE2! (at %lu vs set:%lu)\n",counter,my_end);
						keep_reading = false;
						break;
					}
					//printf("<%u>",counter);
				}
				
				//if hash.first is -2, it means it's a new 'chromosome'.
				else if(hash.first == -2) {
					//printf("New Genome at counter: %lu and i: %u\n",counter,i);
					
					//Find the name of this genome.
					//at position i, we should have the '>'
					char gen_name[MAX_NAME_LEN+5];	//the name can only be 100 chars long
					unsigned int gen_pos=0; 
					unsigned int name_pos = i;
					
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
						
						// also need to store any potentially missed characters
						// We want to include the name so we can find it...
						// But if the name portion is short, we want to at least the last
						// gMER_SIZE bases
						if(name_pos + gMER_SIZE < gBUFFER_SIZE)  {
							// Don't need to store or increment the counter
						}
						else {
							// We need to store the last characters before the name
							//chars_stored += store(genome_vec,gen_piece+gBUFFER_SIZE-gMER_SIZE,gMER_SIZE);
							//unsigned int temp_num = store(genome_vec,gen_piece+gBUFFER_SIZE-gMER_SIZE,gMER_SIZE);
							// Get from where the name began to at least MER_SIZE chars
							//unsigned int temp_num = store(genome_vec,gen_piece+name_pos,gMER_SIZE);
							unsigned int to_read = gBUFFER_SIZE - gMER_SIZE;
							store(genome_vec,gen_piece+to_read,gMER_SIZE, chars_stored);

						}
					
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
						
						//fprintf(stdout,"MAde it to special case #3: %lu (%u)--%u %s\n",
						//	chars_stored,i,gBUFFER_SIZE-i-1,gen_piece+gBUFFER_SIZE-(gBUFFER_SIZE-i-1));
						//fprintf(stdout,"Incrementing offset by %d and hash.second is %d\n",gMER_SIZE - (gBUFFER_SIZE - i),hash.second);
					}
					else {
						;//printf("Nothing special :( i: %u vs OTHER: %u diff: %u\n",i,gBUFFER_SIZE-gMER_SIZE,gBUFFER_SIZE-gMER_SIZE-i);
					}

					gen_name[gen_pos > MAX_NAME_LEN ? MAX_NAME_LEN : gen_pos] = '\0'; 
					
					string str_gen_name = gen_name;
					ReplaceSpaceWithUnderscore(str_gen_name);
					if(gVERBOSE > 1) 
						printf("Found new genome at pos %lu: %s\n",chr_begin_pos,gen_name);

					vector<pair<string,unsigned long> >::iterator vit = names.begin();
					//Find the position of this genome name
					for(unsigned int name_pos=0; name_pos < j; name_pos++)
						vit++;
					vit++;	// the vector inserts to the one "previous", but we want to insert
							// after the current name.

					//printf("Inserting: %s, %d at %d before %s, %d ",str_gen_name.c_str(),0,j,(*vit).first.c_str(),(*vit).second);
					//printf("and after %s, %d",(*(vit-1)).first.c_str(),(*(vit-1)).second);
					names.insert(vit,pair<string,unsigned long>(str_gen_name,0));

					//take care of the former name.
					//1. set the name in the genome to have the correct position
					//2. set the chr_begin_pos var to be the number of chars we've stored so far.
					names[j] = pair<string,unsigned long>(chr_name,chr_begin_pos);
					chr_name = str_gen_name;
				
					counter += hash.second;	//increment the counter the number of bases we've skipped
					
					chr_begin_pos = counter;
					j++;
					
					// Check to see if we've come across the starting point
					if(counter >= my_start)
						start_yet = true;
					// We want to return here, right?  Because we don't care about the rest of the genome...
					if(counter > my_end && !read_all) {
						//fprintf(stderr,"Exiting HERE!\n");
						keep_reading = false;
						// Increment j because that's what will happen at the end of the loop
						j++;
						break;
					}
				}
				
			}
		}
		if(gVERBOSE)
			cout << endl;

		//fprintf(stderr,"Making name at %d (pos=%lu) to %s\n",j,chr_begin_pos,chr_name.c_str());
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

		// Check to see if we've come across the starting point
		if(counter >= my_start)
			start_yet = true;
		if(counter > my_end && !read_all) {
			//fprintf(stderr,"Exiting HERE3!\n");
			keep_reading = false;
			j++;
			break;
		}

	}

	chr_begin_pos = chars_stored;

	//fprintf(stderr,"Making name at %d (pos=%lu) to %s\n",j,chr_begin_pos,"end");
	names[j] = pair<string,unsigned long>("end",chr_begin_pos);

	fprintf(stderr,"Size of genome: %lu\n",chr_begin_pos);

		
#ifdef DEBUG
		for(unsigned int i=0; i<names.size(); i++) {
			printf("Name[%d]: %s\tBegin: %lu\n",i,names[i].first.c_str(),names[i].second);
		}
		fflush(stdout);
		printf("Last Hash number: %lu\n",chars_stored);
#endif

	fix_hash(gen_hash_vec,genome_vec);

	// Need to convert it to the real, vector genome.
	convertToVector(gen_hash_vec, genome_vec);

	// fprintf(stderr,"End of hash_and_store\n");
}


//void GenomeMem::convertToVector(hash_map<unsigned int, vector<unsigned long> > *gh, 
//							 gvector<GenomeLocation> *g) {
void GenomeMem::convertToVector(hash_map<unsigned int, vector<unsigned long> > *gh, 
							 gvector<GEN_TYPE> *pg) {

	fprintf(stderr,"[%d/-] Converting to Vector...\n",iproc);
		
	if(include_hash) {
		fprintf(stderr,"\t[%d/-] Trying to create hash of size %lu\n",iproc,gh->size());
		
		unsigned long count = 0;
		unsigned long max_size = 0;
		unsigned long min_size = ~0;
	
		// Set up the hashed genome
		gen_hash.set_size(gh->size());
	
		
		// This doesn't work if we create a bunch of fragmented hash locations, so we'll just 
		// malloc one large array and set each location to be an index into that array.
		for(unsigned int i=0; i<gh->size(); i++) {
			count += (*gh)[i].size();
		}

		// When we write to the file, we want to remember this information.
		// So we've made it a private class variable
		num_hash_elements = count;
		//fprintf(stderr,"This should be 8: %lu (size of size_t)\n",sizeof(size_t));
		fprintf(stderr,"Trying to create entire array of size %lu (%lu total positions of size %lu)\n",
				sizeof(unsigned long) * count,count,sizeof(unsigned long));
		unsigned long* location_array = 
				(unsigned long*)calloc(sizeof(unsigned long), count);
		//fprintf(stderr,"asserting location array...\n");
		assert(location_array);
		
		//fprintf(stderr,"Saving data ptr...\n");
		// Save the data pointer
		gen_hash.set_data_ptr((void*)location_array);
		//fprintf(stderr,"Saved Data Ptr successfully!\n");

		count = 0;
		for(unsigned int i=0; i<gh->size(); i++) {
			gen_hash[i].size = (*gh)[i].size();
		
			// These next three lines are for debug purposes
			if(gen_hash[i].size > max_size)
				max_size = gen_hash[i].size;
	
			if(gen_hash[i].size < min_size && gen_hash[i].size != 0)
				min_size = gen_hash[i].size;
			
			//gen_hash[i].hash_arr = new unsigned long[(*gh)[i].size()];
			//gen_hash[i].hash_arr = 
			//		(unsigned long*)malloc(sizeof(unsigned long) * (*gh)[i].size());
			// Set the array to be a location in the array
			gen_hash[i].hash_arr = &location_array[count];
			for(unsigned int j=0; j<(*gh)[i].size(); j++) {
				gen_hash[i].hash_arr[j] = (*gh)[i][j];
			}
			count += gen_hash[i].size;
	
			
		}

		// Free this space created by the genome vector
		delete gh;
	
		fprintf(stderr,"\t[%d/-] Finished create hash.\n",iproc);
		fprintf(stderr,"\t[%d/-] Stats: Total hashes is %lu, Longest hash is %lu, shortest is %lu, and average is %f\n",
				iproc,count,max_size,min_size,((double)count)/gen_hash.size());
	}
	
	// Set up the genome
	if(packed_genome)
		free(packed_genome);
	if(amount_genome)
		free(amount_genome);
	
	//fprintf(stderr,"Creating new genome with location of %lu * %u = %lu\n",g->size(),gGEN_SIZE,g->size() * gGEN_SIZE);
	fprintf(stderr,"\t[%d/-] Trying to create a new genome with a size of %lu...",
			iproc,pg->size() * GEN_PACK);
	
	genome_size = pg->size() * GEN_PACK;

	//packed_genome = (GEN_TYPE*)malloc(sizeof(GEN_TYPE) * pg->size());
	//memcpy(packed_genome,pg->getArrPtr(),sizeof(GEN_TYPE) * pg->size());
	// Use the gvector's internal memory location.  Must delete when we're done.
	packed_genome = pg->saveArrPtr();
#ifdef DEBUG	
	fprintf(stderr,"Saved Pointer...");
	fprintf(stderr,"malloc'ing %ld...",sizeof(float) * genome_size/gGEN_SIZE);
#endif

	amount_genome = (float*)calloc(sizeof(float), genome_size/gGEN_SIZE);
	fprintf(stderr,"Success!\n");
	
#ifdef SET_POS
	if(gs_positions)
		free(gs_positions);
	fprintf(stderr,"\t[%d/-] Trying to malloc %lu elements for positions array...",
			iproc,pg->size());
	gs_positions = (char*)malloc(pg->size());
	fprintf(stderr,"Success!\n");
	gs_positions_size = pg->size();
#endif

	
#ifdef SET_POS
	if(!packed_genome || !amount_genome || !gs_positions)
		throw(new Exception("Could not malloc genome"));
#else
	if(!packed_genome || !amount_genome)
		throw(new Exception("Could not malloc genome"));
#endif

	memset(amount_genome, 0, sizeof(float)*genome_size/gGEN_SIZE);
#ifdef SET_POS
	memset(gs_positions, 0, pg->size());
#endif
	
	// Clear up the memory we've used
	delete pg;
	
	fprintf(stderr,"[%d/-] Finished Vector Conversion\n",iproc);
	
}



/********************************/
/* Write the hash out to a file */
/********************************/
void GenomeMem::saveHash(FILE* saveF, bool ldebug) {
	if(ldebug) fprintf(stderr,">> Writing hash out to file...\n");
	if(ldebug) fprintf(stderr,"\tMer Size: %d\n\tNum Hash Elements: %lu\n",gMER_SIZE,num_hash_elements);
	//fprintf(saveF,"%ld\n",gen_hash.size());
	unsigned long tempsize = gen_hash.size();
	fwrite(&tempsize, sizeof(unsigned long), 1, saveF);
	fwrite(&gMER_SIZE, sizeof(int), 1, saveF);
	fwrite(&num_hash_elements, sizeof(unsigned long), 1, saveF);
	fwrite(gen_hash.get_data_ptr(), sizeof(unsigned long), num_hash_elements, saveF);

	//fprintf(stderr,"%ld\n",gen_hash.size());
	for(unsigned int i=0; i<gen_hash.size(); i++) {
		if(gen_hash[i].size > 0) {
			unsigned int hash_size = gen_hash[i].size;
			fwrite(&i, sizeof(unsigned int), 1, saveF);
			fwrite(&hash_size, sizeof(unsigned int), 1, saveF);
			
		}
	}
	// Write out an end marker
	if(ldebug) fprintf(stderr,">> Writing end marker\n");
	//fprintf(saveF,"0 0\n");
	unsigned int tempmarker = 0;
	fwrite(&tempmarker, sizeof(unsigned int), 1, saveF);
	fwrite(&tempmarker, sizeof(unsigned int), 1, saveF);
}
	
void GenomeMem::fix_hash() {
	//Means there is no user-defined max-hash.
	if(gMAX_HASH_SIZE <= 0)
		return;

	for(unsigned int i=0; i<gHASH_SIZE; i++)
		if(gen_hash[i].size > gMAX_HASH_SIZE) {
#ifdef DEBUG
			if(gVERBOSE > 1)
				printf("Clearing hash at pos %d with size %d\n",i,gen_hash[i].size);
#endif
			// Don't delete it here because we've allocated a huge array elsewhere
			//delete[] gen_hash[i].hash_arr;
			gen_hash[i].hash_arr = 0;
			gen_hash[i].size = 0;
		}
}

void GenomeMem::fix_hash(hash_map<unsigned int, vector<unsigned long> > *gh, 
							 gvector<GEN_TYPE> *pg) {
							 
	//Means there is no user-defined max-hash.
	if(gMAX_HASH_SIZE <= 0)
		return;

	for(unsigned int i=0; i<gHASH_SIZE; i++)
		if((*gh)[i].size() > gMAX_HASH_SIZE) {
#ifdef DEBUG
			if(gVERBOSE > 1)
				printf("Clearing hash at pos %d with size %d\n",i,(*gh)[i].size);
#endif
			(*gh)[i].clear();
		}
}

HashLocation* GenomeMem::GetMatches(string &str) {
	pair<int, unsigned long> p_hash = bin_seq::get_hash(str);
	if(!p_hash.first)
		return NULL;
		
	return &gen_hash[p_hash.second];
}

HashLocation* GenomeMem::GetMatches(unsigned int hash) {

	return &gen_hash[hash];
}

bool GenomeMem::Test(ostream &os, unsigned int &warnings) {
	bool success = true;
	GenomeMem gensnp;
	
	// Test the max_pos function
	float a[7] = {0, 1, 2, 3, 4, 5, 6};
	//TEST(gensnp.max_pos(a,7) == 6);
	a[5] = 7;
	//TEST(gensnp.max_pos(a,7) == 5);
	
	
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
	float noSNP2[] = {0.52999, 1.57500, 7.02047, 2.49373, 0.00000};
	pval = gensnp.is_snp(noSNP2,i,j,b);
	fprintf(stderr,"%f, %f, %f, %f, %f\n",noSNP2[0],noSNP2[1],noSNP2[2],noSNP2[3],noSNP2[4]);
	fprintf(stderr,"PVAL: %.2e at %d and second is %d (second should probably be none)\n",pval,i,b?j:-1);
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
		GenomeMem genYumei;
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


		GenomeMem gen45;
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
