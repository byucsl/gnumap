/* hash_map.h
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

#ifndef HASH_MAP_H
#define HASH_MAP_H

#include "const_include.h"

template <typename KEY_TYPE, class MAP_TYPE>
class hash_map {
	public:
		hash_map<KEY_TYPE, MAP_TYPE>():
			num_elements(0), max_size(DEFAULT_SIZE) {
			//_hash_map = new MAP_TYPE[max_size];
			_hash_map = NULL;
			Init(DEFAULT_SIZE);
		}
		
		hash_map<KEY_TYPE, MAP_TYPE>(size_t size) :
			num_elements(0), max_size(size) {
			//_hash_map = new MAP_TYPE[size];
			_hash_map = NULL;
			Init(size);
		}

		~hash_map<KEY_TYPE, MAP_TYPE>(){
			//fprintf(stderr,"Hash Map Destructor\n");
			delete[] _hash_map;
			//free(_hash_map);
		}

		/*
		MAP_TYPE& operator[](KEY_TYPE &hash_spot) {
			if(static_cast<unsigned long>(hash_spot) > max_size)
				resize(max_size * 10);
		
			return _hash_map[static_cast<unsigned long>(hash_spot)];
		}
		*/
		
		MAP_TYPE& operator[](size_t hash_spot) {
			if(hash_spot > max_size)
				resize(max_size * 10);
		
			return _hash_map[hash_spot];
		}

		void set_size(size_t new_size) {
			resize(new_size);
		}
		
		size_t size() {
			return max_size;
		}

		void clear() {
			Init(DEFAULT_SIZE);
		}
		
		/* These two functions are intended to help with un-fragmented memory locations.
		 * They don't do any memory allocation themselves (I won't delete this pointer),
		 * but I'll save it.
		 *
		 * @precondition dpt MUST have been allocated with malloc, NOT new
		 *
		 * @postcondition dpt will be saved until free_data_ptr is called
		 */
		void set_data_ptr(void* dpt) {
			_data_ptr = dpt;
		}
		void* get_data_ptr() {
			return _data_ptr;
		}
		bool has_data_ptr() {
			return _data_ptr ? true : false;
		}
		void free_data_ptr() {
			free(_data_ptr);
		}

	private:
		void resize(size_t new_size) {
			if(gVERBOSE > 1) 
				fprintf(stderr,"HASH_MAP: Calling resize(%lu) with former at %lu and max at %lu\n",new_size,num_elements,max_size);
			MAP_TYPE* hash_temp = new MAP_TYPE[new_size];
			if(!hash_temp)
				throw "Unable to allocate hash.  Try a smaller mer size";

			//MAP_TYPE* hash_temp = (MAP_TYPE*)malloc(sizeof(MAP_TYPE)*new_size);
			//for(unsigned long i=0; i<max_size && i<new_size; i++)
			// If we haven't set an element, don't try and set it in the new array
			for(unsigned long i=0; i<num_elements && i<new_size; i++)
				hash_temp[i] = _hash_map[i];
				
			if(_hash_map)
				delete[] _hash_map;
				//free(_hash_map);
			_hash_map = hash_temp;
			max_size = new_size;
		}

		void Init(size_t size) {
			num_elements = 0;
			max_size = size;
			if(_hash_map)
				delete[] _hash_map;
				//free(_hash_map);

			_hash_map = new MAP_TYPE[size];
			//_hash_map = (MAP_TYPE*)malloc(sizeof(MAP_TYPE)*size);

			_data_ptr = 0;
		}
	
		MAP_TYPE* _hash_map;
		size_t num_elements;
		size_t max_size;
		
		static const size_t DEFAULT_SIZE=1024;
		
		void* _data_ptr;
};


#endif //HASH_MAP_H
