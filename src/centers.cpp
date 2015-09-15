/* centers.cpp
 *
 * Copyright (C) 2009-2014 Scott Lloyd
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

/*
$Id: $

Description: Setup and find profile space centers

$Log: $
*/

#include <stdio.h>
#include <stdlib.h> /* qsort, malloc, abort */
#include <string.h> /* memcpy */
#include <math.h> /* modff, sqrtf */

#include "centers.h"

// #undef DIM_CENTERS
// #define DIM_CENTERS 3
// #undef __SSE2__

//#define METHOD_G 1
//#define METHOD_H 1
//#define METHOD_I 1
#define METHOD_J 1
//#define METHOD_K 1

#define MNAME(x) #x
#define MVALUE(x) MNAME(x)

#define RAGGED_ARRAY 1

#define EPSILON 0.0001f

#define PRINT_RATIONALS 0

#define lengthof(a) (sizeof(a)/sizeof(a[0]))

size_t ccount, max_centers;
center_t **centers;
unsigned int** sum_mat;

static center_t *new_center(void);
static int compare_centers(const void *a, const void *b);
static void unique_centers(void);

#ifdef METHOD_G
#define METHOD G
/* Generates all points for LEVELS in DIM_CENTERS dimensions. */
/* All points are normalized, sorted, and then duplicates are removed. */
//#define LEVELS 1 // @ DIM 3 -> N_CENTERS =   7
//#define LEVELS 2 // @ DIM 3 -> N_CENTERS =  19
//#define LEVELS 3 // @ DIM 3 -> N_CENTERS =  49
//#define LEVELS 2 // @ DIM 5 -> N_CENTERS = 211

void init_centers(void)
{
	int i, sum;
	int coord[DIM_CENTERS];
	center_t *pos;
	//int icount = 0;

	sum = 0;
	for (i = 0; i < DIM_CENTERS; i++) coord[i] = 0;
	for (;;) {
		/* increment */
		for (i = 0; i < DIM_CENTERS; i++) {
			++sum;
			if (++coord[i] <= LEVELS) break;
			sum -= coord[i];
			coord[i] = 0;
		}
		if (i == DIM_CENTERS) break;
		//icount++;
		//for (i = DIM_CENTERS-1; i >= 0; i--) printf(" %d", coord[i]); printf("\n");

		/* allocate space */
		pos = new_center();

		/* calculate point in profile space from lattice coord */
		for (i = 0; i < DIM_CENTERS; i++) {
			pos[(DIM_CENTERS-1)-i] = (center_t)coord[i]/sum;
		}
	}
	//printf("icount:%d\n", icount);
	qsort(centers, N_CENTERS, sizeof(center_t *), compare_centers);
	unique_centers();
}
#endif

#ifdef METHOD_H
#define METHOD H
/* Generates points regularly spaced on a hyperplane */
//#define LEVELS 4 // @ DIM 3 -> N_CENTERS = 15
//#define LEVELS 5 // @ DIM 3 -> N_CENTERS = 21
//#define LEVELS 6 // @ DIM 3 -> N_CENTERS = 28
//#define LEVELS 4 // @ DIM 5 -> N_CENTERS = 70
//#define LEVELS 5 // @ DIM 5 -> N_CENTERS = 126

void init_centers(void)
{
	int i, sum;
	int coord[DIM_CENTERS];
	center_t *pos;
	//int icount = 0;

	sum = 0;
	for (i = 0; i < DIM_CENTERS; i++) coord[i] = 0;
	for (;;) {
		/* increment */
		for (i = 0; i < DIM_CENTERS; i++) {
			++sum;
			if (++coord[i] <= LEVELS && sum <= LEVELS) break;
			sum -= coord[i];
			coord[i] = 0;
		}
		if (i == DIM_CENTERS) break;
		if (sum < LEVELS) {
			coord[0] += LEVELS-sum;
			sum = LEVELS;
		}
		//icount++;
		//for (i = DIM_CENTERS-1; i >= 0; i--) printf(" %d", coord[i]); printf("\n");

		/* allocate space */
		pos = new_center();

		/* calculate point in profile space from lattice coord */
		for (i = 0; i < DIM_CENTERS; i++) {
			pos[(DIM_CENTERS-1)-i] = (center_t)coord[i]/LEVELS;
		}
	}
	//printf("icount:%d\n", icount);
}
#endif

#ifdef METHOD_I
#define METHOD I
/* Generates points along dimensions at 1/2 the distance to zero. */
/* LEVELS must be a power of 2 */
//#define LEVELS  4 // @ DIM 3 -> N_CENTERS =   9
//#define LEVELS  8 // @ DIM 3 -> N_CENTERS =   9
//#define LEVELS 16 // @ DIM 3 -> N_CENTERS =   9
//#define LEVELS 16 // @ DIM 5 -> N_CENTERS = 185

void init_centers(void)
{
	int i, sum;
	int coord[DIM_CENTERS];
	center_t *pos;
	//int icount = 0;

	sum = 0;
	for (i = 0; i < DIM_CENTERS; i++) coord[i] = 0;
	for (;;) {
		/* increment */
		for (i = 0; i < DIM_CENTERS; i++) {
			int last = coord[i];
			coord[i] = last ? last << 1 : 1;
			sum += coord[i] - last;
			if (coord[i] <= LEVELS && sum <= LEVELS) break;
			sum -= coord[i];
			coord[i] = 0;
		}
		if (i == DIM_CENTERS) break;
		//icount++;
		if (sum != LEVELS) continue;
		//for (i = DIM_CENTERS-1; i >= 0; i--) printf(" %d", coord[i]); printf("\n");

		/* allocate space */
		pos = new_center();

		/* calculate point in profile space from lattice coord */
		for (i = 0; i < DIM_CENTERS; i++) {
			pos[(DIM_CENTERS-1)-i] = (center_t)coord[i]/LEVELS;
		}
	}
	//printf("icount:%d\n", icount);
}
#endif

#ifdef METHOD_J
#define METHOD J
/* Generates points in a lattice first */
/* and then finds the closest point on a hyperplane */
/*
DIM 5 :
	LEVELS 4 :
		+0-0 -> N_CENTERS = 125, filter A = 113, filter B =  77 (note symmetry +4)
		+0-1 -> N_CENTERS =  59, filter A =  51, filter B =  35
		+3-1 -> N_CENTERS =  39, filter A =  31, filter B =  23
	LEVELS 5 :
		+0-0 -> N_CENTERS = 251, filter A = 219, filter B = 155
		+0-1 -> N_CENTERS = 129, filter A = 117, filter B =  81 (note symmetry -4)
		+3-1 -> N_CENTERS =  74, filter A =  70, filter B =  50
*/
#define BL 0 /* BAND low offset */
#define BH 1 /* BAND high offset */
#define elim(x,y) (pos[x]+pos[y] > 0.75 && fabsf(pos[x]-pos[y]) < 0.30) /* filter A */
//#define elim(x,y) (pos[x]+pos[y] > 0.50 && fabsf(pos[x]-pos[y]) < 0.30) /* filter B */

void init_centers(void)
{
	int i, sum;
	int coord[DIM_CENTERS];
	center_t *pos;
	//int icount = 0;

	sum = 0;
	for (i = 0; i < DIM_CENTERS; i++) coord[i] = 0;
	for (;;) {
		/* increment */
		for (i = 0; i < DIM_CENTERS; i++) {
			++sum;
			if (++coord[i] <= LEVELS && sum <= LEVELS) break;
			sum -= coord[i];
			coord[i] = 0;
		}
		if (i == DIM_CENTERS) break;
		if (sum < (LEVELS-BAND)) {
			coord[0] += (LEVELS-BAND)-sum;
			sum = (LEVELS-BAND);
		}
		/* at this point, sum is between: LEVELS-BAND <= sum <= LEVELS */

		/* narrow the band */
		#if BL
		if (sum < LEVELS-BAND+BL) continue;
		#endif
		#if BH
		if (sum > LEVELS-BH) continue;
		#endif
		/* at this point, sum is between: LEVELS-BAND+BL <= sum <= LEVELS-BH */

		//icount++;
		//for (i = DIM_CENTERS-1; i >= 0; i--) printf(" %d", coord[i]); printf(" = %d\n", sum);

		/* allocate space */
		pos = new_center();

		/* calculate point in profile space from lattice coord */
		for (i = 0; i < DIM_CENTERS; i++) {
			pos[(DIM_CENTERS-1)-i] = (center_t)coord[i]/LEVELS;
		}
		/* find closest point on hyperplane */
		if (sum < LEVELS) {
			center_t t = 1;
			for (i = 0; i < DIM_CENTERS; i++) t -= pos[i];
			t /= DIM_CENTERS;
			for (i = 0; i < DIM_CENTERS; i++) pos[i] += t;
		}
		// This seems to help even with 60 centers (BAliBASE RV65)
		#if defined(elim) && DIM_CENTERS >= 4
		if (
			elim(0,1) || elim(0,3) ||
			elim(1,0) || elim(1,2) ||
			elim(2,1) || elim(2,3) ||
			elim(3,0) || elim(3,2)
			) {
			free(pos); ccount--;
		}
		#endif
	}
	// add a few more centers
	// leaving out 'all gaps' seems to help (BAliBASE RV65)
	#if BH
	for (i = 0; i < DIM_CENTERS-1; i++) {
		pos = new_center();
		memset(pos, 0, sizeof(center_t)*DIM_CENTERS);
		pos[i] = 1;
		// pos = new_center();
		// memset(pos, 0, sizeof(center_t)*DIM_CENTERS);
		// pos[i] = (center_t)1.0/RLEVELS;
		// pos[DIM_CENTERS-1] = (center_t)1.0 - pos[i];
	}
	#endif
	#if 0
	for (i = 0; i < DIM_CENTERS-1; i++) {
		for (j = RLEVELS >> 1; j > 0; j >>= 1) {
			center_t t = 1;
			pos = new_center();
			memset(pos, 0, sizeof(center_t)*DIM_CENTERS);
			pos[i] = (center_t)j/RLEVELS;
			t -= pos[i];
			t /= DIM_CENTERS;
			for (k = 0; k < DIM_CENTERS; k++) pos[k] += t;
		}
	}
	qsort(centers, N_CENTERS, sizeof(center_t *), compare_centers);
	unique_centers();
	#endif
	//printf("icount:%d\n", icount);
}
#endif

#ifdef METHOD_K
#define METHOD K

#if 0 /* 5 */
#if DIM_CENTERS != 5
#error DIM_CENTERS must equal 5
#endif
static center_t _centers[][5] = {
	{1.00, 0.00, 0.00, 0.00, 0.00},
	{0.00, 1.00, 0.00, 0.00, 0.00},
	{0.00, 0.00, 1.00, 0.00, 0.00},
	{0.00, 0.00, 0.00, 1.00, 0.00},
	{0.00, 0.00, 0.00, 0.00, 1.00},
};
#endif

#if 1 /* 31 */
#if DIM_CENTERS != 5
#error DIM_CENTERS must equal 5
#endif
static center_t _centers[][5] = {
	/*      A         C         G         T         - */
	{0.050000, 0.050000, 0.050000, 0.050000, 0.800000},
	{0.050000, 0.050000, 0.050000, 0.300000, 0.550000},
	{0.050000, 0.050000, 0.050000, 0.550000, 0.300000},
	{0.050000, 0.050000, 0.050000, 0.800000, 0.050000},
	{0.050000, 0.050000, 0.300000, 0.050000, 0.550000},
	{0.050000, 0.050000, 0.300000, 0.300000, 0.300000},
//	{0.050000, 0.050000, 0.300000, 0.550000, 0.050000},
	{0.050000, 0.050000, 0.550000, 0.050000, 0.300000},
//	{0.050000, 0.050000, 0.550000, 0.300000, 0.050000},
	{0.050000, 0.050000, 0.800000, 0.050000, 0.050000},
	{0.050000, 0.300000, 0.050000, 0.050000, 0.550000},
	{0.050000, 0.300000, 0.050000, 0.300000, 0.300000},
	{0.050000, 0.300000, 0.050000, 0.550000, 0.050000},
	{0.050000, 0.300000, 0.300000, 0.050000, 0.300000},
	{0.050000, 0.300000, 0.300000, 0.300000, 0.050000},
//	{0.050000, 0.300000, 0.550000, 0.050000, 0.050000},
	{0.050000, 0.550000, 0.050000, 0.050000, 0.300000},
	{0.050000, 0.550000, 0.050000, 0.300000, 0.050000},
//	{0.050000, 0.550000, 0.300000, 0.050000, 0.050000},
	{0.050000, 0.800000, 0.050000, 0.050000, 0.050000},
	{0.300000, 0.050000, 0.050000, 0.050000, 0.550000},
	{0.300000, 0.050000, 0.050000, 0.300000, 0.300000},
//	{0.300000, 0.050000, 0.050000, 0.550000, 0.050000},
	{0.300000, 0.050000, 0.300000, 0.050000, 0.300000},
	{0.300000, 0.050000, 0.300000, 0.300000, 0.050000},
	{0.300000, 0.050000, 0.550000, 0.050000, 0.050000},
	{0.300000, 0.300000, 0.050000, 0.050000, 0.300000},
	{0.300000, 0.300000, 0.050000, 0.300000, 0.050000},
	{0.300000, 0.300000, 0.300000, 0.050000, 0.050000},
//	{0.300000, 0.550000, 0.050000, 0.050000, 0.050000},
	{0.550000, 0.050000, 0.050000, 0.050000, 0.300000},
//	{0.550000, 0.050000, 0.050000, 0.300000, 0.050000},
	{0.550000, 0.050000, 0.300000, 0.050000, 0.050000},
//	{0.550000, 0.300000, 0.050000, 0.050000, 0.050000},
	{0.800000, 0.050000, 0.050000, 0.050000, 0.050000},

	{1.000000, 0.000000, 0.000000, 0.000000, 0.000000},
	{0.000000, 1.000000, 0.000000, 0.000000, 0.000000},
	{0.000000, 0.000000, 1.000000, 0.000000, 0.000000},
	{0.000000, 0.000000, 0.000000, 1.000000, 0.000000},
//	{0.000000, 0.000000, 0.000000, 0.000000, 1.000000},
};
#endif

void init_centers(void)
{
	int i;

	if ((centers = (center_t **)malloc(sizeof(center_t *)*lengthof(_centers))) == NULL) abort();
	for (i = 0; i < lengthof(_centers); i++) centers[i] = _centers[i];
	ccount = lengthof(_centers);
}
#endif

static center_t *new_center(void)
{
	if (max_centers == 0) {
		max_centers = 128;
		if ((centers = (center_t **)malloc(sizeof(center_t *)*max_centers)) == NULL) abort();
	} else if (ccount == max_centers) {
		max_centers <<= 1;
		if ((centers = (center_t **)realloc(centers, sizeof(center_t *)*max_centers)) == NULL) abort();
	}
	if ((centers[ccount] = (center_t *)malloc(sizeof(center_t)*DIM_CENTERS)) == NULL) abort();
	return centers[ccount++];
}

static int compare_centers(const void *a, const void *b)
{
	center_t *pa = *(center_t **)a;
	center_t *pb = *(center_t **)b;
	int i;

	for (i = 0; i < DIM_CENTERS; i++) {
		center_t delta = pa[i] - pb[i];
		if (delta < -EPSILON) return -1;
		else if (delta > EPSILON) return 1;
	}
	return 0;
}

static void unique_centers(void)
{
	unsigned int i, j;

	for (i = 0; i < N_CENTERS-1; i++) {
		if (compare_centers(&centers[i], &centers[i+1])) continue;
		for (j = i+1; j < N_CENTERS-1; j++) memcpy(centers[j], centers[j+1], sizeof(center_t)*DIM_CENTERS);
		free(centers[j]);
		N_CENTERS--;
		i--;
	}
}

/* 
 * Adjust the center based on previous and new values
 */
center_d adjust_center(float prevV, float newV, center_d d) {
	int i;
	// No value. Don't adjust
	if(prevV == 0)
		return 0;

	center_t c[DIM_CENTERS];
	center_t *map = centers[d];
	float running_sum = 0;
	for(i=0; i<DIM_CENTERS; i++) {
		c[i] = map[i]*prevV/newV;
		running_sum += c[i];
	}
	// Needs to be normalized
	for(i=0; i<DIM_CENTERS; i++) {
		c[i] /= running_sum;
	}
#ifdef FIND_EXHAUSTIVE
	return find_center_exhaustive(c);
#else
	return find_center(c);
#endif
}

#if PRINT_RATIONALS
static void to_rational(center_t value, int *a, int *b)
{
	int i;
	center_t ipart;

	for (i = 1; i < 1024; i++) {
		if (modff(value*i, &ipart) < EPSILON) break;
	}
	*a = (int)ipart;
	*b = i;
}
#endif

void print_centers(FILE *stream)
{
	unsigned int i, j;

	fprintf(stream, "centers:%lu\n", N_CENTERS);
	for (i = 0; i < N_CENTERS; i++) {
		fprintf(stream, "%3d", i);
		for (j = 0; j < DIM_CENTERS; j++) fprintf(stream, " %f", centers[i][j]);
		// for (j = 0; j < DIM_CENTERS; j++) fprintf(stream, " %08X", *((int *)(&centers[i][j])));
		// fprintf(stream, "\t{");
		// for (j = 0; j < DIM_CENTERS; j++) fprintf(stream, "%f%s", centers[i][j], j < DIM_CENTERS-1 ? ", " : "");
		// fprintf(stream, "},");
		#if PRINT_RATIONALS
		for (j = 0; j < DIM_CENTERS; j++) {
			int a, b;
			to_rational(centers[i][j], &a, &b);
			fprintf(stream, " %d/%d", a, b);
		}
		#endif
		fprintf(stream, "\n");
	}
}

// use: c += dx*dx; or c += fabsf(dx);
#define distance2(a, b, c) \
{ \
	register int jj; \
	c = 0; \
	for (jj = 0; jj < DIM_CENTERS; jj++) { \
		register center_t dx = (b)[jj]-(a)[jj]; \
		c += dx*dx; \
	} \
}

// use: return sqrtf(temp); or return temp; */
center_t distance(center_t *pos1, center_t *pos2)
{
	center_t temp;

	distance2(pos1, pos2, temp);
	return sqrtf(temp);
}

/* find nearest partition center by exhaustive search */

int find_center_exhaustive(center_t *pos)
{
	unsigned int ii, index = 0;
	center_t dist = DIM_CENTERS; /* farthest distance */

	for (ii = 0; ii < N_CENTERS; ii++) {
		register center_t temp;
		distance2(pos, centers[ii], temp);
		if (temp < dist) {
			index = ii;
			dist = temp;
		}
	}
	return index;
}

/* NOTES:
 * Finding the closest point on the hyperplane to build the reduce_lut
   does not seem to change the choice of center very much.
 * Using an additive vs. a squared distance metric to find the closest
   center does not seem to change the choice of center very much.
 * Increasing RLEVELS (sample points in the table) does seem to reduce
   the differences between a table-lookup search and an exhaustive search.
 * The sum of coords range (BAND) from quantization is DIM_CENTERS-1.
 * Using a reduce lookup table with power of 2 dimensions is about the
   same speed as with [RLEVELS+1].
 * The ragged array reduce lookup table is somewhat slower, but the
   total memory usage is better.
 * The Intel Core2 Duo E6600 has 2x32K of L1 cache, 4M of shared L2 cache,
   Supplemental SSE3 (SSSE3), and EM64T (FPGA system).
 * The Intel Core2 Duo T5500 has 2x32K of L1 cache, 2M of shared L2 cache,
   Supplemental SSE3 (SSSE3), and EM64T (laptop).
 * See http://en.wikipedia.org/wiki/List_of_Intel_Core_2_microprocessors
 * See http://en.wikipedia.org/wiki/SSSE3
 */

static size_t rcount;
static size_t reduce_size;

#if RAGGED_ARRAY == 1
#if DIM_CENTERS == 3
unsigned char ***reduce_lut;
#elif DIM_CENTERS == 4
unsigned char ****reduce_lut;
#elif DIM_CENTERS == 5
unsigned char *****reduce_lut asm ("reduce_lut");
#endif

static void *dup(void *aptr, size_t nbytes)
{
	void *tmp = malloc(nbytes);

	if (tmp == NULL) abort();
	reduce_size += nbytes;
	return memcpy(tmp, aptr, nbytes);
}

void init_lookup(void)
{
	int i, sum;
	int coord[DIM_CENTERS];
	center_t pos[DIM_CENTERS];
	void *arrn[DIM_CENTERS][RLEVELS+1];
	unsigned char arr0[RLEVELS+1];

	sum = 0;
	for (i = 0; i < DIM_CENTERS; i++) coord[i] = 0;
	for (;;) {
		/* increment */
		for (i = 0; i < DIM_CENTERS; i++) {
			++sum;
			if (++coord[i] <= RLEVELS && sum <= RLEVELS) break;
			switch (i) {
			case 0:
				arrn[1][coord[1]] = dup(arr0, coord[0]*sizeof(unsigned char));
				break;
			case DIM_CENTERS-1:
#if DIM_CENTERS == 3
				reduce_lut = (unsigned char***)
#elif DIM_CENTERS == 4
				reduce_lut = (unsigned char****)
#elif DIM_CENTERS == 5
				reduce_lut = (unsigned char*****)
#endif
						dup(arrn[i], coord[i]*sizeof(void *));
				break;
			default:
				arrn[i+1][coord[i+1]] = dup(arrn[i], coord[i]*sizeof(void *));
			}
			sum -= coord[i];
			coord[i] = 0;
		}
		if (i == DIM_CENTERS) break;
		if (sum < (RLEVELS-BAND)) {
			coord[0] += (RLEVELS-BAND)-sum;
			sum = (RLEVELS-BAND);
		}
		rcount++;
		//for (i = DIM_CENTERS-1; i >= 0; i--) printf(" %d", coord[i]); printf("\n");

		for (i = 0; i < DIM_CENTERS; i++) pos[(DIM_CENTERS-1)-i] = (center_t)coord[i]/RLEVELS;
		arr0[coord[0]] = find_center_exhaustive(pos);
	}
}
#endif

#if RAGGED_ARRAY != 1
#if DIM_CENTERS == 3
unsigned char reduce_lut[RLEVELS+1][RLEVELS+1][RLEVELS+1];
#elif DIM_CENTERS == 4
unsigned char reduce_lut[RLEVELS+1][RLEVELS+1][RLEVELS+1][RLEVELS+1];
#elif DIM_CENTERS == 5
unsigned char reduce_lut[RLEVELS+1][RLEVELS+1][RLEVELS+1][RLEVELS+1][RLEVELS+1] asm ("reduce_lut");
#endif

void init_lookup(void)
{
	int i, sum;
	int coord[DIM_CENTERS];

	sum = 0;
	for (i = 0; i < DIM_CENTERS; i++) coord[i] = 0;
	for (;;) {
		/* increment */
		for (i = 0; i < DIM_CENTERS; i++) {
			++sum;
			if (++coord[i] <= RLEVELS && sum <= RLEVELS) break;
			sum -= coord[i];
			coord[i] = 0;
		}
		if (i == DIM_CENTERS) break;
		if (sum < (RLEVELS-BAND)) {
			coord[0] += (RLEVELS-BAND)-sum;
			sum = (RLEVELS-BAND);
		}
		rcount++;
		//for (i = DIM_CENTERS-1; i >= 0; i--) printf(" %d", coord[i]); printf("\n");

		{
			int i, index;
			center_t pos[DIM_CENTERS];

			for (i = 0; i < DIM_CENTERS; i++) pos[(DIM_CENTERS-1)-i] = (center_t)coord[i]/RLEVELS;
			index = find_center_exhaustive(pos);
			#if 0
			{
				int index2;
				center_t t, sum2 = 0;
				center_t pos2[DIM_CENTERS];
				/* find closest point to hyperplane first */
				t = 1;
				for (i = 0; i < DIM_CENTERS; i++) t -= pos[i];
				t /= DIM_CENTERS;
				for (i = 0; i < DIM_CENTERS; i++) sum2 += pos2[i] = t + pos[i];
				index2 = find_center_exhaustive(pos2);
				if (index != index2) {
					static int diff_count;
					int jj;
					diff_count++;
					printf(" -- dist:%d t:%f indexA:%d indexB:%d diff_count:%d\n", sum-RLEVELS, t, index, index2, diff_count);
					printf("posA : ");
					for (jj = 0; jj < DIM_CENTERS; jj++) {
						fprintf(stdout, " %f", pos[jj]);
					}
					fputc('\n', stdout);
					printf("posB : ");
					for (jj = 0; jj < DIM_CENTERS; jj++) {
						fprintf(stdout, " %f", pos2[jj]);
					}
					fprintf(stdout, " = %f", sum2);
					fputc('\n', stdout);
					printf("centerA : ");
					for (jj = 0; jj < DIM_CENTERS; jj++) {
						fprintf(stdout, " %f", centers[index][jj]);
					}
					fputc('\n', stdout);
					printf("centerB : ");
					for (jj = 0; jj < DIM_CENTERS; jj++) {
						fprintf(stdout, " %f", centers[index2][jj]);
					}
					fputc('\n', stdout);
				}
				index = index2;
			}
			#endif
			/* Least significant coord needs to go on the right for efficiency */
			#if DIM_CENTERS == 3
			reduce_lut[coord[2]][coord[1]][coord[0]] = index;
			#elif DIM_CENTERS == 4
			reduce_lut[coord[3]][coord[2]][coord[1]][coord[0]] = index;
			#elif DIM_CENTERS == 5
			reduce_lut[coord[4]][coord[3]][coord[2]][coord[1]][coord[0]] = index;
			#endif
		}
	}
	reduce_size = sizeof(reduce_lut);
}
#endif

#if __GNUC__ && __SSE2__ && RAGGED_ARRAY && DIM_CENTERS == 5 && 0
// 66.6 Mpos/s with RLEVELS = 16 and pos unaligned
#define FIND a
#if __x86_64
#define deref "mov\t(%1,%2,8), %1\n\t"
#else
#define deref "mov\t(%1,%2,4), %1\n\t"
#endif
typedef float v4sf __attribute__ ((vector_size (16)));
const v4sf rlev = {(float)RLEVELS, (float)RLEVELS, (float)RLEVELS, (float)RLEVELS};

int find_center(center_t *pos)
{
	int ret;
	unsigned char *ptr;
	long off;

	__asm__(
		"movups	(%3), %%xmm0\n\t"
		"movss	16(%3), %%xmm1\n\t"
		"mulps	%4, %%xmm0\n\t"
		"mulss	%4, %%xmm1\n\t"

		"mov	%5, %1\n\t"
#if __x86_64
		"xor	%2, %2\n\t"
#endif
		"cvttss2si	%%xmm0, %%edx\n\t"
		"psrldq	$4, %%xmm0\n\t"
		deref
		"cvttss2si	%%xmm0, %%edx\n\t"
		"psrldq	$4, %%xmm0\n\t"
		deref
		"cvttss2si	%%xmm0, %%edx\n\t"
		"psrldq	$4, %%xmm0\n\t"
		deref
		"cvttss2si	%%xmm0, %%edx\n\t"
		deref
		"cvttss2si	%%xmm1, %%edx\n\t"
		"movzbl	(%1,%2), %0"

		: "=r"(ret), "=r"(ptr), "=d"(off)
		: "r"(pos), "x"(rlev), "m"(reduce_lut)
		: "xmm0", "xmm1"
	);

	return(ret);
}

#elif __GNUC__ && __SSE2__ && !RAGGED_ARRAY && DIM_CENTERS == 5 && 0
//  23.1 Mpos/s with RLEVELS = 64 and pos aligned to 16 bytes
//  80.6 Mpos/s with RLEVELS = 32 and pos aligned to 16 bytes
//  90.0 Mpos/s with RLEVELS = 16 and pos unaligned
// 118.6 Mpos/s with RLEVELS = 16 and pos aligned to 16 bytes
#define FIND b
#define R1 (RLEVELS+1)
typedef float v4sf __attribute__ ((vector_size (16)));
typedef   int v4si __attribute__ ((vector_size (16)));
const v4sf rlev = {(float)RLEVELS, (float)RLEVELS, (float)RLEVELS, (float)RLEVELS};
const v4si idx1 = {R1*R1*R1*R1, 0, R1*R1, 0};
const v4si idx2 = {R1*R1*R1, 0, R1, 0};

int find_center(center_t *pos)
{
	int ret;

	__asm__(
		"movups	(%1), %%xmm0\n\t"
		"movss	16(%1), %%xmm2\n\t"
		"mulps	%2, %%xmm0\n\t"
		"mulss	%2, %%xmm2\n\t"

		"cvttps2dq	%%xmm0, %%xmm0\n\t"
		"pshufd	$0xF5, %%xmm0, %%xmm1\n\t" /* 11,11,01,01 = 0xF5 */
		"pmuludq	%3, %%xmm0\n\t"
		"pmuludq	%4, %%xmm1\n\t"
		"paddq	%%xmm1, %%xmm0\n\t"
		"pshufd	$0x4E, %%xmm0, %%xmm1\n\t" /* 01,00,11,10 = 0x4E */
		"paddq	%%xmm1, %%xmm0\n\t"
		"movd	%%xmm0,	%%edi\n\t"

		"cvttss2si	%%xmm2, %%edx\n\t"
		"movzbl	reduce_lut(%%edi,%%edx), %0"

		: "=r"(ret)
		: "r"(pos), "x"(rlev), "x"(idx1), "x"(idx2)
		: "xmm0", "xmm1", "xmm2", "edx", "edi"
	);

	return(ret);
}

#else
// 66.9 Mpos/s with RLEVELS = 16 - ragged array
// 76.8 Mpos/s with RLEVELS = 16 - static array
#define FIND c

int find_center(center_t *pos)
{
	int i, j, k, l, m;

	//return find_center_exhaustive(pos);
	i = (int)(pos[0]*RLEVELS);
	j = (int)(pos[1]*RLEVELS);
	k = (int)(pos[2]*RLEVELS);
	#if DIM_CENTERS > 3
	l = (int)(pos[3]*RLEVELS);
	#else
	l = 0;
	#endif
	#if DIM_CENTERS > 4
	m = (int)(pos[4]*RLEVELS);
	#else
	m = 0;
	#endif
	#if 0
	{
		int qsum = i+j+k+l+m;
		if (qsum < RLEVELS-BAND || qsum > RLEVELS) {
			printf(" -- %d %d %d %d %d = %d\n", i, j, k, l, m, qsum);
		}
	}
	#endif
	#if DIM_CENTERS == 3
	return reduce_lut[i][j][k];
	#elif DIM_CENTERS == 4
	return reduce_lut[i][j][k][l];
	#elif DIM_CENTERS == 5
	return reduce_lut[i][j][k][l][m];
	#endif
}
#endif

void print_stats(FILE *stream)
{
	fprintf(stream, "method: %s%s\n", MVALUE(METHOD), MVALUE(FIND));
	fprintf(stream, "centers: %lu\n", ccount);
	fprintf(stream, "RLEVELS: %u\n", RLEVELS);
	fprintf(stream, "LEVELS: %u\n", LEVELS);
	fprintf(stream, "DIM_CENTERS: %d\n", DIM_CENTERS);
	fprintf(stream, "rcount: %lu\n", rcount);
	fprintf(stream, "size reduce_lut: %lu\n", reduce_size);
	fprintf(stream, "N_CENTERS: %lu\n", N_CENTERS);

	/*
	int i,j;
	for(i=0; i<ccount; i++) {
		fprintf(stream, "%d",i);
		for(j=0; j<DIM_CENTERS; j++) {
			fprintf(stream, "\t%f",__LINE__);
		}
		fprintf(stream, "\n");
	}
	*/
}

/* for summing two centers */
void init_sums() {
	unsigned int i,j,k;
	/* Init and clear sum_mat matrix */
	sum_mat = (unsigned int**)malloc(sizeof(unsigned int*)*N_CENTERS);
	for(i=0; i<N_CENTERS; i++) {
		sum_mat[i] = (unsigned int*)malloc(sizeof(unsigned int)*N_CENTERS);
		memset(sum_mat[i], 0, sizeof(unsigned int)*N_CENTERS);
	}

	/* Calculate the centers as a weighted average, and then the closest point */
	for(i=0; i<N_CENTERS; i++) {
		sum_mat[i][i] = i;
		for(j=i+1; j<N_CENTERS; j++) {
			center_t* ta = centers[i];
			center_t* tb = centers[j];
			center_t tboth[DIM_CENTERS];
			
			/* Sum both of them up, keeping track of the running total */
			float total = 0;
			for(k=0; k<DIM_CENTERS; k++) {
				tboth[k] = ta[k]+tb[k];
				total+=tboth[k];
			}   
			/* Normalize so it doesn't crash the find_center function call */ 
			for(k=0; k<DIM_CENTERS; k++)
				tboth[k] /= total;

#ifdef FIND_EXHAUSTIVE
            sum_mat[i][j] = sum_mat[j][i] = find_center_exhaustive(tboth);
#else
            sum_mat[i][j] = sum_mat[j][i] = find_center(tboth);
#endif

        }
    }

}

#ifdef TEST_CENTERS

#include <ctype.h> /* isdigit */

#include "ticks.h" /* timing macros: tget, tsec */


/* align a number to the next multiple of "size" where "size" is a power of 2 */
#define ALIGN2(n, size) ( ((size_t)(n) + ((size_t)(size)-1)) & ~((size_t)(size)-1) )

static void *malloc1D(size_t n, size_t s)
{
	return malloc(n*s);
}

static void *malloc2D(size_t m, size_t n, size_t s)
{
	size_t tabsize = ALIGN2(m*sizeof(void *),16);
	size_t stride = ALIGN2(n*s,16);
	void *tmp = (void *)ALIGN2(malloc(tabsize + m*stride + 16),16);
	void *arr;
	int i;

	if (tmp == NULL) return tmp;
	arr = tmp + tabsize;
	for (i = 0; i < m; i++) {
		((void **)tmp)[i] = arr;
		arr += stride;
	}
	return tmp;
}

static void rand_point(center_t *pos)
{
	int i;
	center_t sum = 0;

	for (i = 0; i < DIM_CENTERS; i++) sum += pos[i] = (center_t)rand()/RAND_MAX;
	for (i = 0; i < DIM_CENTERS; i++) pos[i] /= sum;
}

#ifndef VERSION
#define VERSION ?.?
#endif

#define OFLAG 0x01
#define SFLAG 0x02
#define VFLAG 0x1000

#define DEFAULT_ITR 1000
#define DEFAULT_LEN (1 << 14)

int flags; /* argument flags */
size_t iarg = DEFAULT_ITR; /* number of iterations */
size_t larg = DEFAULT_LEN; /* length of profile */

int main(int argc, char *argv[])
{
	int nok = 0;
	char *s;

	while (--argc > 0 && (*++argv)[0] == '-')
		for (s = argv[0]+1; *s; s++)
			switch (*s) {
			case 'i':
				if (isdigit(s[1])) iarg = atoi(s+1);
				else nok = 1;
				s += strlen(s+1);
				break;
			case 'l':
				if (isdigit(s[1])) larg = atoi(s+1);
				else nok = 1;
				s += strlen(s+1);
				break;
			case 'o':
				flags |= OFLAG;
				break;
			case 's':
				flags |= SFLAG;
				iarg = 1;
				larg = 100;
				break;
			case 'v':
				flags |= VFLAG;
				break;
			default:
				nok = 1;
				fprintf(stderr, " -- not an option: %c\n", *s);
				break;
			}

	if (flags & VFLAG) fprintf(stderr, "Version: %s\n", MVALUE(VERSION));
	if (nok || (argc > 0 && *argv[0] == '?')) {
		fprintf(stderr, "Usage: centers -osv -i<int> -l<int>\n");
		fprintf(stderr, "  -i  iteration count, default: %d\n", DEFAULT_ITR);
		fprintf(stderr, "  -l  length of profile, default: %d\n", DEFAULT_LEN);
		fprintf(stderr, "  -o  output VTK point cloud\n");
		fprintf(stderr, "  -s  show differences between methods\n");
		fprintf(stderr, "  -v  version\n");
		exit(EXIT_FAILURE);
	}

	{
		int i, j, index1, index2;
		double qerror1 = 0, qerror2 = 0;
		double init_time;
		double reduce_time;
		double reduce_length;
		tick_t start, finish;
		center_t **cprof; /* continuous profile */
		unsigned char *dprof; /* discrete profile */

		tget(start);
		init_centers();
		init_lookup();
		tget(finish);
		init_time = tsec(finish, start);
		print_stats(stderr);

		if (flags & OFLAG) {
			int i;
			fprintf(stdout, "\n# vtk DataFile Version 3.0\n");
			fprintf(stdout, "Point Cloud\n");
			fprintf(stdout, "ASCII\n");
			fprintf(stdout, "DATASET POLYDATA\n");
			fprintf(stdout, "POINTS %lu float\n", N_CENTERS);
			print_centers(stdout);
			fprintf(stdout, "VERTICES 1 %lu\n", N_CENTERS+1);
			fprintf(stdout, "%lu", N_CENTERS);
			for (i = 0; i < N_CENTERS; i++) fprintf(stdout, " %d", i);
			fprintf(stdout, "\n");
		}

		srand(10000);
		cprof = malloc2D(larg, DIM_CENTERS, sizeof(center_t));
		if (cprof == NULL) {
			fprintf(stderr, " -- error: can't allocate continuous profile memory\n");
			exit(EXIT_FAILURE);
		}
		dprof = malloc1D(larg, sizeof(unsigned char));
		if (dprof == NULL) {
			fprintf(stderr, " -- error: can't allocate discrete profile memory\n");
			exit(EXIT_FAILURE);
		}
		for (i = 0; i < larg; i++) rand_point(cprof[i]);

		tget(start);
		for (j = 0; j < iarg; j++) {
			for (i = 0; i < larg; i++) dprof[i] = find_center(cprof[i]);
		}
		tget(finish);
		reduce_time = tsec(finish, start);
		reduce_length = larg * iarg;

		for (i = 0; i < larg; i++) {
			index1 = find_center_exhaustive(cprof[i]);
			index2 = dprof[i];
			qerror1 += distance(cprof[i], centers[index1]);
			qerror2 += distance(cprof[i], centers[index2]);

			if (flags & SFLAG && index1 != index2) {
				static int diff_count;
				int jj;
				diff_count++;
				fprintf(stdout, " -- indexA:%d indexB:%d diff_count:%d\n", index1, index2, diff_count);
				fprintf(stdout, "original: ");
				for (jj = 0; jj < DIM_CENTERS; jj++) {
					fprintf(stdout, " %f", cprof[i][jj]);
				}
				fputc('\n', stdout);
				fprintf(stdout, "centerA : ");
				for (jj = 0; jj < DIM_CENTERS; jj++) {
					fprintf(stdout, " %f", centers[index1][jj]);
				}
				fputc('\n', stdout);
				fprintf(stdout, "centerB : ");
				for (jj = 0; jj < DIM_CENTERS; jj++) {
					fprintf(stdout, " %f", centers[index2][jj]);
				}
				fputc('\n', stdout);
			}
		}

		int jj;
		for(jj=0; jj<DIM_CENTERS; jj++) {
			fprintf(stdout, " %f", centers[index2][jj]);
		}
		fprintf(stdout, "\n");

		fprintf(stderr, "Method 1 quant. error (avg): %f\n", qerror1/larg);
		fprintf(stderr, "Method 2 quant. error (avg): %f\n", qerror2/larg);
		fprintf(stderr, "Init.  Time      (s): %f\n", init_time);
		fprintf(stderr, "Reduce Length  (pos): %f\n", reduce_length);
		fprintf(stderr, "Reduce Time      (s): %f\n", reduce_time);
		fprintf(stderr, "Reduce Perf. (pos/s): %f\n", reduce_length/reduce_time);
	}
	return(EXIT_SUCCESS); /* or EXIT_FAILURE */
}
#endif

/* TODO:
 * Try other gcc -m options to improve performance.
 * Try finding lookup points only on the hyperplane instead of in a band
   near the hyperplane. This should reduce the lookup table size and
   improve cache performance but will add more operations.
 * Try a search for the profile center using a kd-tree or other structure.
 * Reorder symbols assigned to the dimensions (A, C, G, T, -)? -- Not now
 */

 /* NOTES: (old)
 * Only initialize the hyperplane in the hyperspace matrix.
   - use a sparse matrix method?
   - use a new hyperplane indexing method?
 * Encode the profiles using fractions instead of floats?
 * May be able to get away with only indexing 4 dimensions since
   the 5th is dependent on the others.
 * Use a ragged array for the reduction look-up table. Generate the
   array as a static constant in C source.
 * Will quantizing to RLEV always index to the hyperplane or will
   a guard index be needed on the ends of a range?
 * Try levels for each dimension scaled by log base2
   (e.g. 1/16, 1/8, 1/4, 1/2).
 * Using smaller quanta (1/4) for RLEV than LEVEL seems to help
   alignment quility. If RLEV and LEVEL are the same there shouldn't be
   much difference. See if rounding or the distance metric is the difference.
 * Leverage the hyperplane organization in a gradiant search for the closest
   center in profile space.
 */
