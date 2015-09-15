/* gvector.h
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

#ifndef G_VECTOR_H
#define G_VECTOR_H

#include <cassert>
#include <cstdlib>
#include <execinfo.h>
#include "const_include.h"
#include "Exception.h"

#define GROW_SIZE	2.0

//#define FREE(p) {fprintf(stderr,"Freeing mem here\n");free(p);}
#define FREE(p) {free(p);}

template <typename TYPE>
class gvector {
	public:
		gvector<TYPE>():
			_cur_pos(0), _max_size(DEFAULT_SIZE), retain(true) {
			_vector = NULL;
			Init(DEFAULT_SIZE);
		}
		
		gvector<TYPE>(unsigned long size) :
			_cur_pos(0), _max_size(size), retain(true) {
			//_vector = new MAP_TYPE[size];
			_vector = NULL;
			Init(size);
		}

		~gvector<TYPE>(){
			if(_vector && retain)
				FREE(_vector);
		}

		TYPE& operator[](unsigned long v) {
			// Just throw an error.  You shouldn't overpass the bounds
			if(v > _max_size) {
				void* bt[10];
				size_t size = backtrace(bt,10);
				fprintf(stderr, "Error: overstepped bounds: %lu vs %lu\n",v,_max_size);
				backtrace_symbols_fd(bt, size, 2);
			}
			assert(v <= _max_size);

			return _vector[v];
		}

		void set_size(unsigned long new_size) {
			resize(new_size);
		}
		
		unsigned long size() {
			return _cur_pos;
		}

		void clear() {
			Init(DEFAULT_SIZE);
		}

		/**
		 * The push_back() function appends val to the end of the vector.
		 */
		void push_back(const TYPE &e) {

			if(_cur_pos+1 > _max_size)
				resize();

			_vector[_cur_pos++] = e;
		}

		/**
		 * The pop_back method removes the last element of the vector
		 */
		void pop_back() {
			_cur_pos--;
		}

		void resize(unsigned long new_size) {
			if(gVERBOSE > 1)
				fprintf(stderr,"Resizing Vector to %lu elements\n",new_size);
			if(gVERBOSE > 1)
				fprintf(stderr,"Hot swap...");
			
			// Allocate a new array
			TYPE* temp_vec = (TYPE*)malloc(sizeof(TYPE) * new_size);

			// Means it didn't malloc--it's too big
			if(!temp_vec) {
				//unsigned long inc_size = MAX_LONG;
				fprintf(stderr,"ERROR: Cannot handle a vector of size %lu\n",new_size);
				throw new Exception("Could not allocate array");
			}

			//TYPE* hash_temp = new TYPE[new_size];
			// Copy from the old one
			if(_max_size < new_size)
				memcpy(temp_vec,_vector,sizeof(TYPE) * _max_size);
			else
				memcpy(temp_vec,_vector,sizeof(TYPE) * new_size);
				
			// Delete the old one
			if(_vector && retain)
				FREE(_vector);
			
			// Swap the two
			_vector = temp_vec;
			_max_size = new_size;

			retain = true;
			if(gVERBOSE > 1)
				fprintf(stderr,"Done!\n");
		}

		unsigned long max_size() {
			return _max_size;
		}

		/**
		 * This function will return a pointer to the memory--so it doesn't need to
		 * be malloc'd again.  This should save on memory so we don't need to duplicate it.
		 * 
		 * We need to make sure that after we've return it, we don't delete it so it doesn't
		 * cause memory issues.  So set retain to false and we won't delete when we're done.
		 *
		 * @return a pointer to the internal memory
		 */
		TYPE* saveArrPtr() {
			retain = false;
			return _vector;
		}

		TYPE* getArrPtr() {
			return _vector;
		}

	private:

		/**
		 * The resize method will resize it based on pre-defined conditions
		 */
		inline void resize() {
			// If the max size is below 2G, we'll double it
			if(_max_size < ONE_G) {
				resize((unsigned long)(_max_size * GROW_SIZE));
			}
			// Otherwise, we'll only increment by 1G
			else {
				resize(_max_size+ONE_G);
			}
		}

		void Init(unsigned long size) {
			
			_cur_pos = 0;
			_max_size = size;
			if(_vector && retain)
				FREE(_vector);

			_vector = (TYPE*)malloc(sizeof(TYPE) * size);
			retain = true;
		}
	
		TYPE* _vector;
		unsigned long _cur_pos;
		unsigned long _max_size;

		static const unsigned int ONE_K=1024;
		static const unsigned int ONE_M=ONE_K*1024;
		static const unsigned long ONE_G=ONE_M*1024;

		static const unsigned long MAX_DOUBLE_SIZE=1024*1024*1024;
		static const unsigned long MAX_DOUBLE_INC=1024*1024*1024;
		static const unsigned int MAX_LONG=~0;

		static const unsigned long DEFAULT_SIZE=ONE_M;

		bool retain;	// If true, don't delete memory when you destruct
};


#endif //G_VECTOR_H
