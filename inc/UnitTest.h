/* UnitTest.h
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

#ifndef UNIT_TEST_H
#define UNIT_TEST_H

#define TEST(cond) \
	if (!(cond)) { \
		success = false; \
		os << "!!  Test FAILED at [" << __FILE__ << ", " << __LINE__ << "]  !!" << endl; \
	} \
	else { \
		os << "==  Test Passed at [" << __FILE__ << ", " << __LINE__ << "]  ==" << endl; \
	}

#define TEST_W(cond) \
	if (!(cond)) { \
		warnings++;	\
		os << "!!  Test WARNING at [" << __FILE__ << ", " << __LINE__ << "]  !!" << endl; \
	} \
	else { \
		os << "==  Test Passed at [" << __FILE__ << ", " << __LINE__ << "]  ==" << endl; \
	}

#endif

