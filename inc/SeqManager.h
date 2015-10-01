/* SeqManager.h
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

#ifndef _SEQMANAGER_H
#define _SEQMANAGER_H

#include <vector>
#include <string>
#include <pthread.h>
#include <cassert>
#include "const_include.h"
#include "SeqReader.h"

class SeqManager {

	public:
		SeqManager(Read** reads, const char* readsStr, const unsigned int nR, const unsigned int nP) 
			: g_nReads(nR), g_nProc(nP) {

			Init(reads, readsStr, nR, nP);
			
		}
		
		~SeqManager() {
#ifdef DEBUG_TIME
			fprintf(stderr,"Total wait time is %f\n",wait_time);
#endif
			if(sr)
            {
				delete sr;
            }
		}
		
		void Init(Read** reads, const char* readsStr, const unsigned int nR, const unsigned int nP)
        {
			num_files = 0;
			finished = false;

			g_readPtr = reads;

			sr = new SeqReader(readsStr);
		
			//if(gVERBOSE) 
				fprintf(stderr,"Matching %d file(s):\n",sr->GetNumSeqFiles());
		
			readFiles = sr->GetSeqFileVector();
			rit = readFiles.begin();

			getNextFile();

			totalSeqs = 0;
			readCount = 0;
			fileReadCount = 0;
			//readPerFile = 0;
			readsPerProc = READS_PER_PROC;
			//readsPerProc = rpp;
#ifdef DEBUG
			fprintf(stderr,"Reads per processor: %d (%d/%d perhaps?)\n",readsPerProc,g_nReads,g_nProc);
#else
			fprintf(stderr,"Reads per processor: %d\n",readsPerProc);
#endif

			// init the read lock
			pthread_mutex_init(&read_lock,NULL);
			wait_time = 0;
		}
		
		/**
		 * Resets the counter--adds to the beginning of the global read array
		 */
		void resetCounter()
        {
			if(gVERBOSE)
            {
				printProgress();
            }
			readCount = 0;
		}
		
		unsigned int getTotalSeqCount()
        {
			return totalSeqs;
		}

		void printProgress()
        {
			fprintf(stderr,"[%d/0] %.0f%% reads complete\n",iproc,((float)fileReadCount)/readPerFile * 100);
		}
		
		unsigned int getCurrNumReads()
        {
			return readPerFile;
		}

		/**
		 * Here, we'll read all the reads we can into an array passed in.  This allows the mother
		 * processor to decide how to divy them out.
		 * 
		 * Used by OpenMP
		 * 
		 * @param readArr The read pointer array to which we should add the reads
		 * @param size The size of the read pointer array (readArr)
		 * @param num_reads Will be the number of reads we've added to the array
		 * @return true if we added any reads, false if we're entirely done adding reads
		 */	
		bool fillReadArray(Read** readArr, const unsigned int size, unsigned int &num_reads)
        {
			unsigned int i;
			Read* temp;

#ifdef DEBUG_TIME
			double begin_time = When();
#endif
			MUTEX_LOCK(&read_lock);
#ifdef DEBUG_TIME
			double end_time = When();
#endif

			num_reads = 0;

			if(finished) {
				// Unlock the mutex first, then return
				MUTEX_UNLOCK(&read_lock);
				return false;
			}

#ifdef DEBUG_TIME
			wait_time += end_time-begin_time;
			//fprintf(stderr,"Wait time for fillReadArray is %f, at %f\n",end_time-begin_time,wait_time);
#endif

			for( i = 0; i < size; i++ )
            {
				try
                {
					while( ( temp = sr->GetNextSequence() ) == NULL )
                    {
						if( !getNextFile() )
                        {
							//fprintf(stderr,"CANNOT READ MORE FROM FILE...RETURNING\n");
							MUTEX_UNLOCK(&read_lock);
							return false;
						}
					}
				}
				catch( Exception* e )
                {
					fprintf(stderr,"ERROR(%s:%u):  %s\n",__FILE__,__LINE__,e->GetMessage());
					if(!getNextFile())
                    {
						MUTEX_UNLOCK(&read_lock);
						return false;
					}
					else
                    {
						MUTEX_UNLOCK(&read_lock);
						return true;
					}
				}
				

				// Set the pointer
				readArr[i] = temp;
				num_reads++;
			}

			fileReadCount += num_reads;
			
			MUTEX_UNLOCK(&read_lock);

			return true;	//only return false if we couldn't get the next file
		}
	
		
		/**
		 * Responsible for dishing out reads, as needed
		 * Parameter of 'set' added to be backwards compatable
		 * 
		 * @par begin The beginning location of where the reads will be placed
		 * @par end The ending location of where the reads will be placed
		 * @par set true iff we should set 'begin', false if we should use the value of 'begin'
		 */
		bool getMoreReads(unsigned int &begin, unsigned int &end, bool set)
        {
			// Ensure there's only one thread in here

#ifdef DEBUG_TIME
			double begin_time = When();
#endif
			MUTEX_LOCK(&read_lock);
#ifdef DEBUG_TIME
			double end_time = When();
#endif

			if(finished) {
				//fprintf(stderr,"[-/%d] CANNOT READ MORE--RETURNING\n",iproc);
				if(set) {
					begin = 0; end = 0;
				}
				MUTEX_UNLOCK(&read_lock);
				return false;
			}
			
#ifdef DEBUG_TIME
			wait_time += end_time-begin_time;
			//fprintf(stderr,"Wait time for fillReadArray is %f, at %f\n",end_time-begin_time,wait_time);
#endif
			
				
			Read* temp;
			unsigned int i;
			unsigned int num_reads=0;
			// if we don't set, we'll decrease the number to read at a time
			unsigned int num_to_read = set ? readsPerProc/4 : end-begin;
			bool ret_val = true;
			
			if(set)
            {
				begin = readCount;
            }

			//fprintf(stderr,"[-/%d] Starting with begin=%u,end=%u and fileProgress=%d/%d\n",iproc,begin,end,fileReadCount,readPerFile);
			
			//for(i=0; i< (set ? readsPerProc : end-begin) && fileReadCount+i<readPerFile; i++) {
			for(i=0; i< num_to_read && begin+i<g_nReads; i++)
            {
				// If we've read too many sequences from this file, 
				// if there's no more reads, get the next sequence from the next file
				if(fileReadCount+i >= readPerFile || ((temp = sr->GetNextSequence()) == NULL))
                {
					//fprintf(stderr,"No more reads.  Trying next file\n");
					// get the next file, then return
					totalSeqs += num_reads;
					if(!getNextFile())
                    {
						ret_val = false;
					}

					if(set)
                    {
						end = begin+num_reads;
						readCount = end;
					}

					MUTEX_UNLOCK(&read_lock);
					return ret_val;
				}

/*
				while((temp = sr->GetNextSequence()) == NULL) {
					fprintf(stderr,"******[-/%d] temp is NULL (i=%d) and should be off by %u\n",g_iproc,i,num_reads);
					if(!getNextFile()) {
						//fprintf(stderr,"[-/%d] CANNOT READ MORE FROM FILE...RETURNING (%u-%u) and finished is %d\n",g_iproc,begin,end,finished);
						//fprintf(stderr,"[-/%d] i:%u, total=%u\n",g_iproc,i,fileReadCount);
						ret_val = false;
						goto AFTER_FOR;
					}
				}
*/
				//fprintf(stderr,"[-/%d] temp is NOT NULL!!!\n");

				if(!temp->name || !temp->pwm)
                {
					fprintf(stderr,"ERROR:  NO NAME!!!\n");
                }

				// Set it at the global pointer
				g_readPtr[begin+i] = temp;
				num_reads++;
			}
			
			if(set)
            {
				end = begin+num_reads;
				readCount = end;
			}

			//fprintf(stderr,"Read %d reads, start is %d and end is %d\n",num_reads,begin,end);
			totalSeqs += num_reads;
			fileReadCount += num_reads;
			
			//fprintf(stderr,"Read from %u to %u and has read %u(%u) seqs\n",begin,set ? readsPerProc : end-begin,num_reads,totalSeqs);
			// Release the lock
			MUTEX_UNLOCK(&read_lock);
			return ret_val;
		}
		
		bool isFinished()
        {
			return finished;
		}
		
	private:
	
		bool getNextFile()
        {
			bool invalid = true;
			while( invalid )
            {
				if( rit == readFiles.end() )
                {
					//if(rit == readFiles.begin())
					//fprintf(stderr,"*********Read from %u files (out of %u expected)\n",num_files,readFiles.size());
					//fprintf(stderr,"Found %u sequences\n",totalSeqs);
					finished = true;
					return false;
				}

				try
                {
					unsigned int nBurn = 0;
					unsigned int total_reads = sr->GetNumSeqs(*rit);
				
					// Need to burn sequences to be at the right pointer for each node
					if( nproc > 1 && !gMPI_LARGEMEM )
                    {
						// use +1 because we want to get all of them
						unsigned int each = total_reads / nproc + 1;
						unsigned int begin = each * iproc;
						
						nBurn = begin;
						readPerFile = each;
					}					
					else
                    {
						readPerFile = total_reads;
					}

					//fprintf(stderr,"[-/%d] Should burn %u reads and have %u reads total. total_reads:%u and read:%u\n",
					//		iproc,nBurn,readPerFile,total_reads,totalSeqs);

					// instruct the SeqReader to use this file
					sr->use( *rit, nBurn );
					
					// increment the vector iterator
					rit++;
					num_files++;

					invalid = false;
					fileReadCount = 0;
				}
				catch( const char* err )
                {
					fprintf(stderr,"ERROR(%s:%u): \n\t%s\n",__FILE__,__LINE__,err);
					fprintf(stderr,"\tIn file: %s\n",(*rit).c_str());
					rit++;
					continue;
				}
				catch( Exception *e )
                {
					fprintf(stderr,"ERROR(%s:%u): \n\t%s\n",__FILE__,__LINE__,e->GetMessage());
					fprintf(stderr,"\tIn file %s\n",(*rit).c_str());
					delete e;
					rit++;
					continue;
				}
				catch(...)
                {
					fprintf(stderr,"ERROR(%s:%u): \n\tImproper sequence length or file.\n",__FILE__,__LINE__);
					fprintf(stderr,"\tIn file %s\n",(*rit).c_str());
					rit++;
					continue;
				}
			}
			
			if(gVERBOSE)
            {
				fprintf(stderr,"\n[-/%d] Matching %u sequences of: %s\n",iproc, readPerFile,sr->GetFilename().c_str());
            }
						
			return true;
		}
		
		/* Return the current time in seconds, using a float precision number. */
		double When() {
			struct timeval tp;
			gettimeofday(&tp, NULL);
			return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
		}

		Read** g_readPtr;			// Read pointer from Driver
		SeqReader *sr;				// SeqReader to get Reads from
		vector<string> readFiles;	// The vector from the parsed read files
		vector<string>::iterator rit;	// Vector Iterator for the read files
		unsigned int readPerFile;	// Number of reads per file
		unsigned int readCount;		// Location in the global array for the Read to go
		unsigned int fileReadCount;	// For the percentage of reads that are read
		unsigned int totalSeqs;
		unsigned int g_nReads;		// Number of reads in the global array
		unsigned int g_nProc;		// Number of processors in the program
		unsigned int readsPerProc;	// Number of reads per processor=g_nReads/g_nProc
		
		bool finished;				// Flag when we're finished
		
		pthread_mutex_t read_lock;	// Lock to ensure no duplication

		unsigned int num_files;
		double wait_time;
};

#endif // _SEQ_MANAGER_H
