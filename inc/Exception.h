/* Exception.h
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

#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <string>
#include <cstring>

using namespace std;

struct Exception {
	public:
		Exception(const string msg) {
			errstr = new char[msg.size() + 1];
			strcpy(errstr,msg.c_str());
		}

		Exception(const char* msg) {
			errstr = new char[strlen(msg)];
			strcpy(errstr,msg);
		}

		Exception(const char* msg, const char* FN, const unsigned int LN) {
			errstr = new char[strlen(msg) + strlen(FN) + 25];
			sprintf(errstr,"!!%s:%u %s",FN,LN,msg);
		}

		~Exception() {
			if(errstr)
				delete[] errstr;
		}

		char* GetMessage() {
			return errstr;
		}

		void SetMessage(const char* msg) {
			errstr = new char[strlen(msg)];
			strcpy(errstr,msg);
		}

	private:
		char* errstr;
};

#endif

