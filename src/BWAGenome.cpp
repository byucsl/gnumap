/*
 * BWAGenome.cpp
 *
 * Copyright (C) 2009-2014 Jamison Dance
 * Created by Jamison Dance on 12/4/09.
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

#include "BWAGenome.h"

//Load the BWT into memory
BWAGenome::BWAGenome(const char * bwaIndexLocation) {
	char *str = (char*)calloc(strlen(bwaIndexLocation) + 10, 1);
	strcpy(str, prefix); strcat(str, ".bwt");  bwt[0] = bwt_restore_bwt(str);
	strcpy(str, prefix); strcat(str, ".rbwt"); bwt[1] = bwt_restore_bwt(str);
	free(str);
}

BWAGenome::~BWAGenome() {
	// remove the bwt from memory.
	bwt_destroy(bwt[0]);
	bwt_destroy(bwt[1]);
}

/**
 * Return a string of length "size" begining at position "begin" in the genome.
 */
string BWAGenome::GetString(const unsigned long begin, const unsigned int size) {
	string foo = "";
	return foo;
}
/**
 * Return a HashLocation with an array of all the locations on the genome where
 * this sequence is found.
 */
HashLocation* BWAGenome::GetMatches(unsigned int hash, int length) {
	
	/*string bin_seq::hash2str(const unsigned int h, int length) {*/
	char sequenceString[length];
	unsigned int temp = hash;
	
	for(int i = 1; i <= length; i++) {
		sequenceString[length-i] = gINT2BASE[(temp % 4)];
		temp >>= 2;
	}
	
	sequenceString[length] = '\0';
	
	//}*/
	//basic algo:
	//put hash into a seq struct
	//call bwa_cal_sa_reg_gap? to get suffix array coordinates calculated.
	//turn them into locations on the genome by calling bwa_cal_pac_pos?
	//put them in a HashLocation
	
	//have to take reverse complement of the sequences b/c of BWA weirdness.
	gap_opt_t *opt;
	opt = gap_init_opt();
	
	//set up the sequence
	bwa_seq_t p;
	kseq_t sequence;
	sequence.f = NULL;
	//s is the string. set it equal to the sequence.
	sequence.seq.s = sequenceString;
	sequence.seq.l = length;
	sequence.seq.m = 0; //JD dunno what m is for.
	/*
	 typedef struct __kstring_t {
	 size_t l, m;
	 char *s;
	 } kstring_t;
	 
	 typedef struct {
	 kstring_t name, comment, seq, qual;
	 int last_char;
	 kstream_t *f;
	 } kseq_t;
	 
	 //JD struct to store the sequences
	 typedef struct {
	 char *name;
	 ubyte_t *seq, *rseq, *qual;
	 uint32_t len:20, strand:1, type:2, dummy:1, extra_flag:8;
	 uint32_t n_mm:8, n_gapo:8, n_gape:8, mapQ:8;
	 int score;
	 int clip_len;
	 // alignments in SA coordinates
	 int n_aln;
	 bwt_aln1_t *aln;
	 // alignment information
	 bwtint_t sa, pos; //pos is the position in the genome!
	 uint64_t c1:28, c2:28, seQ:8; // number of top1 and top2 hits; single-end mapQ
	 int n_cigar;
	 uint16_t *cigar;
	 // for multi-threading only
	 int tid;
	 // NM and MD tags
	 uint32_t full_len:20, nm:12;
	 char *md;
	 } bwa_seq_t;
	 */
	int l = length;
	int n_seqs , i;
	long n_trimmed = 0
	long n_tot = 0;
	p->tid = -1; // no assigned to a thread
	p->qual = 0;
	p->full_len = p->clip_len = p->len = l;
	n_tot += p->full_len;
	p->seq = (ubyte_t*)calloc(p->len, 1);
	for (i = 0; i != p->full_len; ++i) {
		p->seq[i] = nst_nt4_table[(int)seq->seq.s[i]];
	}
	
	//JD we might not care about the quality. I don't know what this is.
	if (seq->qual.l) { // copy quality
		p->qual = (ubyte_t*)strdup((char*)seq->qual.s);
		if (trim_qual >= 1) {
			n_trimmed += bwa_trim_read(opt->trim_qual, p);
		}
	}
	
	p->rseq = (ubyte_t*)calloc(p->full_len, 1);
	memcpy(p->rseq, p->seq, p->len);
	seq_reverse(p->len, p->seq, 0); // *IMPORTANT*: will be reversed back in bwa_refine_gapped()
	seq_reverse(p->len, p->rseq, opt->mode & BWA_MODE_COMPREAD);
	p->name = strdup((const char*)seq->name.s);
	{ // trim /[12]$
		int t = strlen(p->name);
		if (t > 2 && p->name[t-2] == '/' && (p->name[t-1] == '1' || p->name[t-1] == '2')) {
			p->name[t-2] = '\0';
		}
	}
	
	
	
	//opt->max_gape = opte;
	opt->mode &= ~BWA_MODE_GAPE;
	
	bwa_aln_core(gBWA_INDEX, seq_file, opt);
}
