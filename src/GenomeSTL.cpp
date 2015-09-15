/* GenomeSTL.cpp
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

#include "GenomeSTL.h"

GenomeSTL::GenomeSTL() : Genome(), hash_data(0) {
#ifdef MAP_USE_STRING
	fprintf(stderr, "# Using STL genome type: string\n");
#else
	fprintf(stderr, "# Using STL genome type: unsigned int\n");
#endif
}

GenomeSTL::GenomeSTL(const char* fn) : Genome(fn), hash_data(0) {
#ifdef MAP_USE_STRING
	fprintf(stderr, "Using STL genome type: string\n");
#else
	fprintf(stderr, "Using STL genome type: unsigned int\n");
#endif
}

GenomeSTL::~GenomeSTL() {
	if(hash_data)
		delete[] hash_data;
}

void GenomeSTL::print_to_file(map<MAP_TYPE, HashLocation> & gh, char* ofn) { 
	printf("Printing to file: %s\n",ofn);
	int SEQ_LENGTH = gMER_SIZE; //This is an archaic #def'd element that needed to be removed.
	ofstream out;
	out.open(ofn,ios_base::out);
	
	bin_seq bs;
	
	unsigned int i=0;
		
	map<MAP_TYPE, HashLocation>::const_iterator git;
	for(git=gh.begin(); git!=gh.end(); ++git) {
		if(git->second.size == 0)
			continue;
		//fprintf(stderr,"[%d] arry position %u (%s) has size %lu\n",
		//		iproc,i,bs.hash2str(i,gMER_SIZE).c_str(),gh[i].size);
		out << i++ << "\t" << git->second.size;
		unsigned int j;
		for(j=0; j<git->second.size; j++) {
			//fprintf(stderr,"\t[%d] gh[%u].hash_arr[%u]\n",iproc,i,j);
			out << '\t' << git->second.hash_arr[j] << "[" << GetPos(git->second.hash_arr[j],true) << "] :" 
				<< GetString(git->second.hash_arr[j],SEQ_LENGTH);
		}
		out << endl;
		
	}

	out << "Finished\n" << endl;
	
	out.close();
}

void GenomeSTL::readHash(FILE* readF, bool ldebug) {
	
	if(ldebug) fprintf(stderr,">> Reading hash from file...\n");
	unsigned long temp_gen_hash_size;
	//fscanf(readF,"%ld\n",&gen_hash_size);
	fread(&temp_gen_hash_size, sizeof(unsigned long), 1, readF);
	gHASH_SIZE = temp_gen_hash_size;
	// fprintf(stderr,"gen_hash_size: %d\n",gHASH_SIZE);

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

	// Read in the hash
	unsigned long count = 0;
	unsigned int hash_size;
#ifdef MAP_USE_STRING
	char mer[gMER_SIZE+1];
#else
	MAP_TYPE mer;
#endif
	while(true) {
		//fscanf(readF,"%d %d\n",&index,&hash_size);
#ifdef MAP_USE_STRING
		// Need to not just read a single char
		fread(&mer, sizeof(char), gMER_SIZE, readF);
		mer[gMER_SIZE] = 0;
#else
		fread(&mer, sizeof(MAP_TYPE), 1, readF);
#endif
		fread(&hash_size, sizeof(unsigned int), 1, readF);
		//fprintf(stderr,"%d %d\n",index,hash_size);
		if(!hash_size) {	//This is a marker for the end of the hash.
			if(ldebug) fprintf(stderr,">> Reading end marker\n");
			break;
		}

		gen_hash[mer].size = hash_size;
		// malloc a new array for this location
		//unsigned long* temp_arr = new unsigned long[hash_size];
		gen_hash[mer].hash_arr = location_array+count;
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

void GenomeSTL::hash_and_store() {
	
	// Records what position we want to start saving at
	bool start_yet = my_start == 0 ? true : false;
	
	// Set some values that are necessary for the 'store()' function
	offset = 0;
	bin_offset = 0;

	gvector<GEN_TYPE> *genome_vec
		= new gvector<GEN_TYPE>;
	map<MAP_TYPE, vector<unsigned long> > *gen_hash_vec
		= new map<MAP_TYPE, vector<unsigned long> >;

	Reader reader_vec;
	unsigned char gen_piece[gBUFFER_SIZE+1];

	printf("Hashing Genome...\n");
	
	bin_seq bs;	

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
		//cout << "File Size: " << file_size << endl;
		bool keep_reading = true;
		unsigned int i;

		while(keep_reading) { 
	
			keep_reading = reader_vec.read(gen_piece);
	
			//If we're not at the end of the file, then we can throw it the whole chunk,
			//without the extra gMER_SIZE bases at the end.
			//When we come to the end of the file, we can give it everything.
			//We'll also return the number of characters stored.
			//printf("Buffer size is: %u and mer size is: %u\n",gBUFFER_SIZE,gMER_SIZE);
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
					store(genome_vec,gen_piece,gBUFFER_SIZE,chars_stored);
				}
				catch(bad_alloc& e) {
					printf("Caught bad alloc around store!\n");
					throw e;
				}
			}
			
	 
			for(i=0; i < gBUFFER_SIZE-gMER_SIZE; i++) {
					
#ifdef MAP_USE_STRING
				pair<int,MAP_TYPE> hash = bs.get_hash_str(i,gen_piece,gBUFFER_SIZE-gMER_SIZE);
#else
				pair<int,MAP_TYPE> hash = bs.get_hash(i,gen_piece,gBUFFER_SIZE-gMER_SIZE);
#endif

				if(hash.first >= 0) {	//if it's zero, that's a result of an N character.
					if(hash.first != 0)	{	//don't want to store N's.
						if(include_hash && (counter % gGEN_SKIP == 0)) { //This allows us a smaller hash size in total.
							try {
								if(start_yet)
									(*gen_hash_vec)[hash.second].push_back(counter);
							}
							catch(bad_alloc& e) {
								printf("Caught bad alloc around push_back!\n");
								throw e;
							}
							//printf("Hash of %s (%d) is %s\n",hash.second.c_str(),hash.first,hash.second.c_str());
						}
						
						// print output
						if( gVERBOSE && ((counter % (file_size/50) ) == 0)) {
							cout << ".";
							cout.flush();
						}

					}
				}
						
				//if hash.first is -2, it means it's a new 'chromosome'.
				if(hash.first == -2) {
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
						//fprintf(stdout,"Incrementing offset by %d and hash.second is %s\n",gMER_SIZE - (gBUFFER_SIZE - i),hash.second.c_str());
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

#ifdef MAP_USE_STRING
					counter += hash.second.size();	//increment the counter the number of bases we've skipped
#else
					counter += hash.second;	//increment the counter the number of bases we've skipped
#endif
					
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
				

				//if hash.first is -1, it means it's found the end-of-line char.
				else if(hash.first == -1 || i >= gBUFFER_SIZE-gMER_SIZE) {
					//if(hash.first == -1)
						//printf("Found the eof char\n");
					break;
				}
				else {	// can't be -2 or -1, so it must be nonnegative.
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

void GenomeSTL::convertToVector(map<MAP_TYPE, vector<unsigned long> > *gh, 
							 gvector<GEN_TYPE> *pg) {

	fprintf(stderr,"[%d/-] Converting to Vector...\n",iproc);
		
	if(include_hash) {
		fprintf(stderr,"\t[%d/-] Trying to create hash of size %lu\n",iproc,gh->size());
		
		unsigned long count = 0;
		unsigned long max_size = 0;
		unsigned long min_size = ~0;
	
		
		// This doesn't work if we create a bunch of fragmented hash locations, so we'll just 
		// malloc one large array and set each location to be an index into that array.
		map<MAP_TYPE, vector<unsigned long> >::iterator mit;
		for(mit=gh->begin(); mit!=gh->end(); ++mit) {
			count += mit->second.size();
		}

		// When we write to the file, we want to remember this information.
		// So we've made it a private class variable
		num_hash_elements = count;
		//fprintf(stderr,"This should be 8: %lu (size of size_t)\n",sizeof(size_t));
		fprintf(stderr,"Trying to create entire array of size %lu (%lu total positions of size %lu)\n",
				sizeof(unsigned long) * count,count,sizeof(unsigned long));
		LOC_ARR_TYPE* location_array = 
				(LOC_ARR_TYPE*)calloc(sizeof(unsigned long), count);
		//fprintf(stderr,"asserting location array...\n");
		assert(location_array);
		
		//fprintf(stderr,"Saving data ptr...\n");
		// Save the data pointer
		hash_data = location_array;

		count = 0;
		for(mit=gh->begin(); mit!=gh->end(); ++mit) {
			MAP_TYPE s = mit->first;
			gen_hash[s].size = mit->second.size();
		
			// These next three lines are for debug purposes
			if(gen_hash[s].size > max_size)
				max_size = gen_hash[s].size;
	
			if(gen_hash[s].size < min_size && gen_hash[s].size != 0)
				min_size = gen_hash[s].size;
			
			gen_hash[s].hash_arr = &location_array[count];
			for(unsigned int j=0; j<mit->second.size(); j++) {
				gen_hash[s].hash_arr[j] = mit->second[j];
			}
			count += gen_hash[s].size;
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
void GenomeSTL::saveHash(FILE* saveF, bool ldebug) {
	if(ldebug) fprintf(stderr,">> Writing hash out to file...\n");
	if(ldebug) fprintf(stderr,"\tMer Size: %d\n\tNum Hash Elements: %lu\n",gMER_SIZE,num_hash_elements);
	//fprintf(saveF,"%ld\n",gen_hash.size());
	unsigned long tempsize = gen_hash.size();
	fwrite(&tempsize, sizeof(unsigned long), 1, saveF);
	fwrite(&gMER_SIZE, sizeof(int), 1, saveF);
	fwrite(&num_hash_elements, sizeof(unsigned long), 1, saveF);
	fwrite(hash_data, sizeof(unsigned long), num_hash_elements, saveF);

	//fprintf(stderr,"%ld\n",gen_hash.size());
	unsigned int hash_size;
	map<MAP_TYPE, HashLocation>::iterator mit;
	for(mit=gen_hash.begin(); mit!=gen_hash.end(); ++mit) {
		if(mit->second.size > 0) {
			hash_size = mit->second.size;
#ifdef MAP_USE_STRING
			fwrite(mit->first.c_str(), sizeof(char), gMER_SIZE, saveF);
#else
			fwrite(&(mit->first), sizeof(MAP_TYPE), 1, saveF);
#endif
			fwrite(&hash_size, sizeof(unsigned int), 1, saveF);
		}
	}
	// Write out an end marker
	if(ldebug) fprintf(stderr,">> Writing end marker\n");
	//fprintf(saveF,"0 0\n");
	unsigned int tempmarker = 0;
#ifdef MAP_USE_STRING
	char tempMer[gMER_SIZE];
	fwrite(tempMer, sizeof(char), gMER_SIZE, saveF);
#else
	unsigned long temp_ul = 0;
	fwrite(&temp_ul, sizeof(MAP_TYPE), 1, saveF);
#endif
	fwrite(&tempmarker, sizeof(unsigned int), 1, saveF);
}

void GenomeSTL::fix_hash() {
	//Means there is no user-defined max-hash.
	if(gMAX_HASH_SIZE <= 0)
		return;

	map<MAP_TYPE, HashLocation>::iterator mit;
	for(mit=gen_hash.begin(); mit!=gen_hash.end(); ++mit) {
		if(mit->second.size > gMAX_HASH_SIZE) {
#ifdef DEBUG
			if(gVERBOSE > 1)
				printf("Clearing hash at pos %d with size %d\n",i,mit->second.size);
#endif
			// Don't delete it here because we've allocated a huge array elsewhere
			//delete[] gen_hash[i].hash_arr;
			mit->second.hash_arr = 0;
			mit->second.size = 0;
		}
	}
}

void GenomeSTL::fix_hash(map<MAP_TYPE, vector<unsigned long> > *gh, 
							 gvector<GEN_TYPE> *pg) {
							 
	//Means there is no user-defined max-hash.
	if(gMAX_HASH_SIZE <= 0)
		return;

	map<MAP_TYPE, vector<unsigned long> >::iterator mit;
	for(mit=gh->begin(); mit!=gh->end(); ++mit)
		if(mit->second.size() > gMAX_HASH_SIZE) {
#ifdef DEBUG
			if(gVERBOSE > 1)
				printf("Clearing hash at pos %d with size %d\n",i,mit->second.size());
#endif
			mit->second.clear();
		}
}

HashLocation* GenomeSTL::GetMatches(string &str) {
#ifdef MAP_USE_STRING
	if(gen_hash.count(str))
		return &gen_hash[str];
	return NULL;
#else
	pair<int, unsigned long> p_hash = bin_seq::get_hash(str);
	if(!p_hash.first)
		return NULL;
		
	if(gen_hash.count(p_hash.second))
		return &gen_hash[p_hash.second];
	return NULL;
#endif
}

HashLocation* GenomeSTL::GetMatches(unsigned int hash) {
#ifdef MAP_USE_STRING
	string s = bin_seq::hash2str(hash, gMER_SIZE);
	return &gen_hash[s];
#else
	return &gen_hash[hash];
#endif
}

