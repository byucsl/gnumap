/* SeqReader.cpp
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

#include "SeqReader.h"

SeqReader::SeqReader() :
	seq_counter(0), offset(0), seq_num(0), this_filename("") {

	for(unsigned int i=0; i<READ_BUFFER; i++) {
		reads[i] = 0;	//clear it all
		//reads[i].pwm = 0;
	}

}

SeqReader::SeqReader(const char* fn) :
	seq_counter(0), offset(0), seq_num(0), this_filename("") {
	Construct(fn);	
}

SeqReader::SeqReader(const SeqReader &other) {
	this_filename = other.this_filename;
	this_type = other.this_type;
	seq_files = other.seq_files;
	seq_num = other.seq_num;
}

SeqReader& SeqReader::operator =(const SeqReader &other) {
	this_filename = other.this_filename;
	this_type = other.this_type;
	seq_files = other.seq_files;
	seq_num = other.seq_num;
	return *this;
}


void SeqReader::Construct(const char* fn) {

	//fprintf(stderr,"CALLING CLEARREAD (%d).  DID WE INTEND TO?\n",__LINE__);
	clearReads(reads);
	
	//Separate the different files.
	char files[strlen(fn)+1];
	strcpy(files,fn);
	
	const char* delims = " \"\n,";
	char* token = strtok(files,delims);
	if(gVERBOSE > 1)
		cout << "Found:" << endl;

	while(token != NULL) {
		if(gVERBOSE > 1)
			cout << "\t" << token << endl;
		string temp_str = token;
		seq_files.push_back(temp_str);
		
		token = strtok(NULL,delims);
	}

}
		
SeqReader::~SeqReader() {
	//in.close();
	
	/*
	// Print the backtrace
	void* bt[10];
	size_t size = backtrace(bt,10);
	backtrace_symbols_fd(bt, size, 2);
	
	fprintf(stderr,"CALLING CLEARREAD (%d).  DID WE INTEND TO?\n",__LINE__);
	*/
	clearReads(reads);
}

void SeqReader::clearReads(Read* readArr[]) {
	//fprintf(stderr,"Calling the ClearRead function\n");
	if(!readArr)
		return;
	
	// changing our mind.  We don't want to delete everything here.
	// we want the SeqDriver to delete them after he's done with them.
	/*Read* tempR = readArr;
	while(tempR->pwm) {
		tempR = 0;
		for(unsigned int i=0; i<tempR->length; i++) {
			delete[] tempR->pwm[i];
		}
		delete[] tempR->pwm;
		delete tempR;	//don't delete this reference.  Keep the array around.
		tempR++;
	}*/
	for(unsigned int i=0; i<READ_BUFFER; i++) {
		readArr[i] = 0;
	}
	num_reads = 0;
}

void SeqReader::Init(const char* fn) {
	file_finished = false;
	offset=0;
	seq_counter=0;
	//reads.clear();
	
	if(reads[0]) {
		//fprintf(stderr,"CALLING CLEARREAD (%d).  DID WE INTEND TO?\n",__LINE__);
		clearReads(reads);
	}

	find_type();	
}

void SeqReader::ReadBatch() {

	if(this_type == PRB)
		get_more_prb();
	else if(this_type == INT)
		get_more_int();
	else if(this_type == FASTQ)
		get_more_fastq();
	else if(this_type == FASTA)
		get_more_fasta();

}

/*
 * There are three file types we're interested it:
 * PRB: \d+[ ]\d+[ ]\d+[ ]\d+(\t\d+[ ]+\d+[ ]+\d+[ ]+\d+)+
 * INT: ((\d*\.?\d+[ ]+){4}\t)+ (no tab on the end of the file
 * FASTQ (from maq.sourceforge.net/fastq.shtml): 
 *   <fastq>   := <block>+
 *   <block>   := @<seqname>\n<seq>\n[+][<seqname>]?\n<qual>\n
 *   <seqname> := [A-Za-z0-9_.:-]+
 *   <seq>     := [A-Za-z\n\.~]+
 *   <qual>    := [!-~\n]+
 
 *   And <qual> is the quality of the highest base, where the Sanger format is calculated
 *   like such:
 *   Q = (ASCI(QUAL[x]) - 33)
 *   And Solexa does:
 *   Q = (ASCI(QUAL[x]) - 64)
 
 *   Then use Q2Prb to find the probability
 */ 
void SeqReader::find_type() {
	ifstream infile(this_filename.c_str());
	
	if(infile.bad() || infile.fail())	
		throw new Exception("Invalid File",__FILE__,__LINE__);
	
	// Look for FASTQ format
	char first;
	if((first = infile.get()) == '@') {
		this_type = FASTQ;
		infile.putback(first);
		return;
	}
	// Look for FASTA
	else if(first == '>') {
		this_type = FASTA;
		infile.putback(first);
		return;
	}

	// Put back the character you've pulled off
	infile.putback(first);
	
	string temp;
	getline(infile,temp, '\t');
	
	//The _int.txt file will have 4 tab-delimintated numbers before any spaces.
	//The _prb.txt file (created by Illumina's pipeline) will have 4 numbers--
	//  separated by spaces--before any tabs.
	string::size_type loc = temp.find(' ',0);
	if(loc != string::npos)
		this_type = PRB;
	else
		this_type = INT;
}

SeqType SeqReader::find_type(const string &str) {
	ifstream infile(str.c_str());
	
	
	if(infile.bad() || infile.fail()) {
		fprintf(stderr,"Bad file is %s\n",str.c_str());
		throw new Exception("Invalid File",__FILE__,__LINE__);
	}
	
	// Look for FASTQ format
	char first;
	if((first = infile.get()) == '@') {
		return FASTQ;
	}
	// Look for FASTA
	else if(first == '>') {
		return FASTA;
	}

	// Put back the character you've pulled off
	infile.putback(first);
	
	string temp;
	getline(infile,temp, '\t');
	
	//The _int.txt file will have 4 tab-delimintated numbers before any spaces.
	//The _prb.txt file (created by Illumina's pipeline) will have 4 numbers--
	//  separated by spaces--before any tabs.
	string::size_type loc = temp.find(' ',0);
	if(loc != string::npos)
		return PRB;
	else
		return INT;
}

void SeqReader::use(const string &str) {
	this_filename = str;
	Init(this_filename.c_str());
	ReadBatch();
}
void SeqReader::use(const char* fn) {
	
	this_filename = fn;
	Init(fn);
	ReadBatch();
}

void SeqReader::use(const string &str, unsigned int nBurn) {
	use(str.c_str(),nBurn);
}
void SeqReader::use(const char* fn, unsigned int nBurn) {
	this_filename = fn;
	Init(fn);
	
	if(nBurn == 0) {
		ReadBatch();
		return;
	}
			
	// Now we need to burn some reads
	switch(this_type) {
	case PRB:
	case INT:
		burnPRBINT(nBurn);
		break;
	case FASTA:
		burnFASTA(nBurn);
		break;
	case FASTQ:
		burnFASTQ(nBurn);
		break;
	default:
		break;
	}
	
	seq_num += nBurn;
	ReadBatch();
}

inline bool SeqReader::burnPRBINT(unsigned int nBurn) {
	ifstream inf(this_filename.c_str());
	//inf.clear();
	
	unsigned int line_count=0;
	
	char temp[gBUFFER_SIZE];
	
	while(line_count < nBurn) {
		inf.getline(temp, gBUFFER_SIZE);
		offset += inf.gcount();
		line_count++;

		if(inf.eof()) {
			file_finished = true;
			fprintf(stderr,"***************ERROR!  Closed before we could burn. Only got %u\n",line_count);
			inf.close();
			return false;
		}
	}

#ifdef DEBUG
	unsigned int extra = 0;
	while(!inf.eof()) {
		inf.getline(temp,gBUFFER_SIZE);
		extra++;
	}
	fprintf(stderr,"Burned %u reads and could still read %u\n",line_count,extra);
#endif

	inf.close();
	return true;
}

inline bool SeqReader::burnFASTA(unsigned int nBurn) {
	//fprintf(stderr,"[%d/-] Burning FASTA!!\n",iproc);
	ifstream inf(this_filename.c_str());
	//inf.clear();
	
	unsigned int read_count=0;
	
	char temp[gBUFFER_SIZE];
	
	while(read_count < nBurn) {
		inf.getline(temp, gBUFFER_SIZE);
		if(temp[0] == '>') {
			read_count++;
			// If this puts us over the edge, don't record the last line
			if(read_count >= nBurn)
				break;
		}
		offset = inf.tellg();

		if(inf.eof()) {
			file_finished = true;
			inf.close();
			//fprintf(stderr,"Just burned %u reads\n",read_count);
			return false;
		}
	}
	//fprintf(stderr,"Just burned %u reads, not EOF, temp is %s and offset is %u\n",read_count,temp,offset);
	inf.close();
	
	return true;
}

inline bool SeqReader::burnFASTQ(unsigned int nBurn) {
	ifstream inf(this_filename.c_str());
	
	unsigned int read_count=0;
	
	char temp[gBUFFER_SIZE];
	
	while(read_count < nBurn) {
		// Read in 4 lines
		inf.getline(temp, gBUFFER_SIZE);
		offset += inf.gcount();
		inf.getline(temp, gBUFFER_SIZE);
		offset += inf.gcount();
		inf.getline(temp, gBUFFER_SIZE);
		offset += inf.gcount();
		inf.getline(temp, gBUFFER_SIZE);
		offset += inf.gcount();
		
		read_count++;

		if(inf.eof()) {
			file_finished = true;
			inf.close();
			return false;
		}
	}

	inf.close();
	return true;
}

Read* SeqReader::GetNextSequence() {
	bool has_more=true;

	
	//if(seq_counter >= reads.size()) {
	if(seq_counter >= num_reads) {
		if(this_type == PRB)
			has_more = get_more_prb();
		else if(this_type == INT)
			has_more = get_more_int();
		else if(this_type == FASTQ)
			has_more = get_more_fastq();
		else if(this_type == FASTA)
			has_more = get_more_fasta();
		else
			throw new Exception("Could not determine sequence file type");
	}
	
	if(!has_more)
		return NULL;

	return reads[seq_counter++]; 
}

//vector<vector<double> > SeqReader::PeekNextSequence() {
Read* SeqReader::PeekNextSequence() {

	//if(seq_counter > reads.size()) {
	if(seq_counter > num_reads) {
		if(this_type == PRB)
			get_more_prb();
		else
			get_more_int();
	}
	
	if(seq_counter >= num_reads)
		return NULL;
	
	return reads[seq_counter];
}

bool SeqReader::get_more_prb() {

	//fprintf(stderr,"Getting more reads, seq_num is at %u and num_reads is at %u\n",seq_num,num_reads);
	
	//fprintf(stderr,"CALLING CLEARREAD (%d).  DID WE INTEND TO?\n",__LINE__);
	clearReads(reads);
	//reads.clear();

	if(file_finished)
		return false;

	seq_counter = 0;

	ifstream in;
	in.open(this_filename.c_str());
	in.seekg(offset);
	
	if(in.eof()) {
		file_finished = true;
		return false;
	}

	char* token;
	const char* delims = "\n\r\t ";

#ifdef DEBUG
	if(gVERBOSE > 1)
		cout << "Reading more reads" << endl;
#endif
	
	char *temp_chr = new char[gBUFFER_SIZE];
	//READ_BUFFER is the amount of reads we'll read into memory.
	for(unsigned int to_read = 0; to_read < READ_BUFFER; to_read++, seq_num++) {

#ifdef DEBUG
		if(gVERBOSE > 1 && (to_read % 100 == 0 ) )	{
			cout << ".";
			cout.flush();
		}
#endif

		vector< vector<double> > one_read;
		
		if(in.eof()) {
			delete[] temp_chr;
			file_finished = true;
			in.close();
			if(to_read > 0)
				return true;
			return false;
		}

		in.getline(temp_chr,gBUFFER_SIZE);
		offset += in.gcount();
		if(in.gcount() == 0) {
			//perror("Maybe an invalid read...");
			//printf("It's happening at %d out of %d...\n",to_read,READ_BUFFER);
			continue;
		}

		token = strtok(temp_chr,delims);
		int first = 1;

		while(token != NULL) {
			vector<double> one_char(4);
	
			char* a_check;
			if(first) {
				a_check = token;
				first = 0;
			}
			else
				 a_check = strtok(NULL,delims);

			char* c_check = strtok(NULL,delims);
			char* g_check = strtok(NULL,delims);
			char* t_check = strtok(NULL,delims);

			//If they're all NULL, then we're done
			if(!a_check && !c_check && !g_check && !t_check)
				break;

			//this would only have one not being the right length.
			//that should throw an error--they all should be nothing.
			if( (!a_check) || (!c_check) || (!g_check) || (!t_check) ) {
				if(gVERBOSE > 1) {
					fprintf(stderr,"to_read(%d) A: %s C: %s G: %s T: %s\n",to_read,a_check,c_check,g_check,t_check);
				}
				string errmsg = "Bad Sequence Probability File at ";
				errmsg += a_check;
				errmsg += ",";
				errmsg += c_check;
				errmsg += ",";
				errmsg += g_check;
				errmsg += ",";
				errmsg += t_check;
				throw new Exception(errmsg);
			}
				
			double a = atof(a_check);
			double c = atof(c_check);
			double g = atof(g_check);
			double t = atof(t_check);
			
			a = Q2Prb_ill(a);
			c = Q2Prb_ill(c);
			g = Q2Prb_ill(g);
			t = Q2Prb_ill(t);

			// Normalize them so they all add up to 1
			double sum = a+c+g+t;
			a /= sum;
			c /= sum;
			g /= sum;
			t /= sum;

#ifdef DEBUG
			if(gVERBOSE > 2)
				cout << "a: " << a << "\tc: "  << c << "\tg: " << g << "\tt: "  << t << endl;
				cout << "\ttotal: " << a+c+g+t << endl;
#endif

			one_char[0] = a;
			one_char[1] = c;
			one_char[2] = g;
			one_char[3] = t;

			one_read.push_back(one_char);
			
		}
		
		//  if there was an g_adaptor sequence defined, remove it from the end
		//  of this read.
		if(g_adaptor)
			FixReads(one_read);
			
		//Read* temp_read = new Read;
		reads[num_reads] = new Read;
		reads[num_reads]->length = one_read.size();
		reads[num_reads]->pwm = new float*[one_read.size()];

		for(unsigned int i=0; i<one_read.size(); i++) {
			reads[num_reads]->pwm[i] = new float[4];
			for(unsigned int j=0; j<4; j++) {
				reads[num_reads]->pwm[i][j] = (float)one_read[i][j];
			}
		}

		reads[num_reads]->name = new char[20];
		sprintf(reads[num_reads]->name,"seq%u",seq_num);

#ifdef DEBUG
		//for(unsigned int q=0; q<reads[num_reads]->length; q++) {
		//	printf("%f %f %f %f\n",reads[num_reads]->pwm[q][0],reads[num_reads]->pwm[q][1],
		//						   reads[num_reads]->pwm[q][2],reads[num_reads]->pwm[q][3]);
		//}
		string consensus = GetConsensus(reads[num_reads]); 
		printf("[%d] Length: %d\tconsensus: %s\n",num_reads,reads[num_reads]->length,consensus.c_str());
		//cout << "consensus: " << consensus << "\t";
		printf("\n\n");
#endif
		
		num_reads++;
		//&reads[num_reads++] = temp_read;
		//reads.push_back(one_read);
	}
	
	in.close();
	
	delete[] temp_chr;
	return true;

}

inline double SeqReader::Q2Prb_ill(double Q) {
	//double answer = 1.0 - 1.0/ (pow( 10.0,(Q/10.0) ) + 1.0);
	double answer = 1.0 - 1.0/ (pow( 10.0,(Q/10.0) )); //illumina+64 v1.5
	return answer > 1.0 ? 1.0 : answer;
}
inline double SeqReader::Q2Prb_std(double Q) {
	double answer = 1-exp((-Q/10.0)*log(10.0));
	
	return answer > 1.0 ? 1.0 : answer;
}

bool SeqReader::get_more_int() {
	clearReads(reads);
	//reads.clear();
	
	if(file_finished)
		return false;
		
	seq_counter = 0;
	ifstream in;
	in.open(this_filename.c_str());
	in.seekg(offset);
	
	if(in.eof()) {
		file_finished = true;
		return false;
	}

	char* token;
	const char* delims = "\t ";

	if(gVERBOSE > 1)
		cout << "Reading more reads" << endl;
	
	char *temp_chr = new char[gBUFFER_SIZE];
	//READ_BUFFER is the amount of reads we'll read into memory.
	for(unsigned int to_read = 0; to_read < READ_BUFFER; to_read++) {
#ifdef DEBUG
		if(gVERBOSE > 1 && (to_read % 100 == 0 ) )	{
			cout << ".";
			cout.flush();
		}
#endif
		vector< vector<double> > one_read;

		//string temp;
		//getline(in,temp);

		in.getline(temp_chr,gBUFFER_SIZE);
		offset += in.gcount();
		if(in.eof()) {
			delete[] temp_chr;
			file_finished = true;
			if(to_read > 0)
				return true;
			return false;
		}

#ifdef DEBUG
		if(gVERBOSE > 2)
			cout << temp_chr << endl;
#endif
			
		token = strtok(temp_chr,delims);

		//There are 4 params that don't matter to us...
		for(int i=0; i<4; i++)
			token = strtok(NULL,delims);

		int first = 1;

		while(token != NULL) {
			vector<double> one_char(4);

			char* a_check;
			if(first) {
				a_check = token;
				first = 0;
			}
			else
				 a_check = strtok(NULL,delims);

			char* c_check = strtok(NULL,delims);
			char* g_check = strtok(NULL,delims);
			char* t_check = strtok(NULL,delims);
			
			if(!a_check && !c_check && !g_check && !t_check)
				break;

			//this would only have one not being the right length.
			//that should throw an error--they all should be nothing.
			if( (!a_check) || (!c_check) || (!g_check) || (!t_check) )
				throw new Exception("Bad Sequence Intensity File");
			
			double a = atof(a_check);
			double c = atof(c_check);
			double g = atof(g_check);
			double t = atof(t_check);
			
			//The next few lines deal with a negative sequence probability,
			//so that the probability of a base is never > 1 
			//
			double shift = 0;
			if(a < shift)
				shift = a;
			if(c < shift)
				shift = c;
			if(g < shift)
				shift = g;
			if(t < shift)
				shift = t;
				
			if(shift < 0) {
				shift = shift * -1;
				a += shift;
				c += shift;
				g += shift;
				t += shift;
			}
			
			double total = a + c + g + t;
			
			if(total != 0) {
				a /= total;
				c /= total;
				g /= total;
				t /= total;
			}
#ifdef DEBUG
			if(gVERBOSE > 2)
				cout << "a: " << a << "\tc: "  << c << "\tg: " << g << "\tt: "  << t << endl;
#endif
			one_char[0] = a;
			one_char[1] = c;
			one_char[2] = g;
			one_char[3] = t;

			one_read.push_back(one_char);
			//one_read[i] = one_char;
			//delete[] begin;
		}

		//  if there was an g_adaptor sequence defined, remove it from the end
		//  of this read.
		if(g_adaptor)
			FixReads(one_read);

		//Read* temp_read = new Read;
		reads[num_reads] = new Read;
		reads[num_reads]->length = one_read.size();
		reads[num_reads]->pwm = new float*[one_read.size()];

		for(unsigned int i=0; i<one_read.size(); i++) {
			reads[num_reads]->pwm[i] = new float[4];
			for(unsigned int j=0; j<4; j++) {
				reads[num_reads]->pwm[i][j] = (float)one_read[i][j];
			}
		}

		reads[num_reads]->name = new char[10];
		sprintf(reads[num_reads]->name,"seq%u",seq_num);
		
		
#ifdef DEBUG
		for(unsigned int q=0; q<reads[num_reads]->length; q++) {
			printf("%f %f %f %f\n",reads[num_reads]->pwm[q][0],reads[num_reads]->pwm[q][1],
								   reads[num_reads]->pwm[q][2],reads[num_reads]->pwm[q][3]);
		}
		string consensus = GetConsensus(reads[num_reads]); 
		printf("[%d] Length: %d\tconsensus: %s\n",num_reads,reads[num_reads]->length,consensus.c_str());
		//cout << "consensus: " << consensus << "\t";
		printf("\n\n");
#endif

		num_reads++;
		//&reads[num_reads++] = temp_read;
		//reads.push_back(one_read);
	}
	
	in.close();
	
	delete[] temp_chr;

	return true;
}

/************************************/
/*           FASTA					*/
/************************************/
bool SeqReader::get_more_fasta() {
	clearReads(reads);
	//reads.clear();

	if(file_finished)
		return false;

	seq_counter = 0;

	ifstream in;
	in.open(this_filename.c_str());
	in.seekg(offset);
	
	if(in.eof()) {
		file_finished = true;
		return false;
	}

	if(gVERBOSE > 1)
		cerr << "Reading more reads FASTA style, with offset=" <<offset<< endl;
	
	string name, sequence, temp_sequence;

	//READ_BUFFER is the amount of reads we'll read into memory.
	for(unsigned int to_read = 0; to_read < READ_BUFFER; to_read++) {
#ifdef DEBUG
		if(gVERBOSE > 1 && (to_read % 100 == 0 ) )	{
			cout << ".";
			cout.flush();
		}
#endif
		if(in.eof()) {
			file_finished = true;
			in.close();

			if(to_read > 0) {
				return true;
			}
			return false;
		}

		getline(in,name);
		// Check to see if we're at the end of the file now
		if(in.eof()) {
			file_finished = true;
			
			in.close();
			if(to_read > 0)
				return true;
			return false;
		}
		getline(in,sequence);
		while(in.peek() != '>' && !in.eof()) {
			getline(in,temp_sequence);
			sequence += temp_sequence;
		}
#ifdef DEBUG
		if(iproc) fprintf(stderr,"[%d/-] At %u, Name of %s and sequence with %s\n",
				iproc,to_read,name.c_str(),sequence.c_str());
#endif
		vector< vector<double> > one_read;

		for(unsigned int i=0; i<sequence.size(); i++) {
			vector<double> one_char(4);

			// Set them all to zero
			one_char[A_POS] = 0.0;
			one_char[C_POS] = 0.0;
			one_char[G_POS] = 0.0;
			one_char[T_POS] = 0.0;

			switch(tolower(sequence[i])) {
				case 'a':
					one_char[A_POS] = 1.0;
					break;
				case 'c':
					one_char[C_POS] = 1.0;
					break;
				case 'g':
					one_char[G_POS] = 1.0;
					break;
				case 't':
					one_char[T_POS] = 1.0;
					break;

				/**********************************************************/
				// Here are the ambiguity characters.  We'll allow this
				/**********************************************************/
				case 'r': // G or A for puRine
					one_char[A_POS] = 0.5;
					one_char[G_POS] = 0.5;
					break;
				case 'y': // T or C for pYrimidine
					one_char[C_POS] = 0.5;
					one_char[T_POS] = 0.5;
					break;
				case 'k': // G or T for Ketone
					one_char[G_POS] = 0.5;
					one_char[T_POS] = 0.5;
					break;
				case 'm': // A or C for aMino group
					one_char[A_POS] = 0.5;
					one_char[C_POS] = 0.5;
					break;
				case 's': // G or C for Strong interaction
					one_char[C_POS] = 0.5;
					one_char[G_POS] = 0.5;
					break;
				case 'w': // A or T for Weak interaction
					one_char[A_POS] = 0.5;
					one_char[T_POS] = 0.5;
					break;
				case 'b': // C, G, T (not A) B comes after a
					one_char[C_POS] = 1.0/3.0;
					one_char[G_POS] = 1.0/3.0;
					one_char[T_POS] = 1.0/3.0;
					break;
				case 'd': // A, G, T (not C) D comes after c
					one_char[A_POS] = 1.0/3.0;
					one_char[G_POS] = 1.0/3.0;
					one_char[T_POS] = 1.0/3.0;
					break;
				case 'h': // A, C, T (not G) H comes after g
					one_char[A_POS] = 1.0/3.0;
					one_char[C_POS] = 1.0/3.0;
					one_char[T_POS] = 1.0/3.0;
					break;
				case 'v': // a, c, g (not t or u) V comes after u
					one_char[A_POS] = 1.0/3.0;
					one_char[C_POS] = 1.0/3.0;
					one_char[G_POS] = 1.0/3.0;
					break;
				case 'n':
					one_char[A_POS] = 0.25;
					one_char[C_POS] = 0.25;
					one_char[G_POS] = 0.25;
					one_char[T_POS] = 0.25;
					break;
				
				// Whitespace is just ignored
				case 10:
				case 11:
				case 12:
				case 13: //line feed
				case ' ':
					continue;

				default:
					fprintf(stderr,"[%d/-] Error at fasta sequence (#%u) with name=>%s< and sequence=>%s<\n",
						iproc,num_reads,name.c_str(),sequence.c_str());
					fprintf(stderr,"Bad char is (%d)%c at pos %d\n",(int)sequence[i],sequence[i],i);
					char* errstr = new char[22];
					strcpy(errstr,"Invalid Sequence at ");
					strcat(errstr,&sequence[i]);
					throw errstr;
					break;	// redundant
			}

#ifdef DEBUG
			if(gVERBOSE > 2)
				cout << "a: " << one_char[A_POS] << "\tc: "  << one_char[A_POS] << "\tg: " << one_char[A_POS] << "\tt: "  << one_char[A_POS] << endl;
#endif

			one_read.push_back(one_char);
		}

		//  if there was an g_adaptor sequence defined, remove it from the end
		//  of this read.
		if(g_adaptor)
			FixReads(one_read);

		reads[num_reads] = new Read;
		reads[num_reads]->length = one_read.size();
		reads[num_reads]->pwm = new float*[one_read.size()];
		reads[num_reads]->name = new char[name.size()+1];
		strcpy(reads[num_reads]->name,name.c_str()+1);	//+1 for the '>' at the beginning
		reads[num_reads]->seq = sequence;

		for(unsigned int i=0; i<one_read.size(); i++) {
			reads[num_reads]->pwm[i] = new float[4];
			for(unsigned int j=0; j<4; j++) {
				reads[num_reads]->pwm[i][j] = (float)one_read[i][j];
			}
		}

#ifdef DEBUG
		for(unsigned int q=0; q<reads[num_reads]->length; q++) {
			printf("%f %f %f %f\n",reads[num_reads]->pwm[q][0],reads[num_reads]->pwm[q][1],
								   reads[num_reads]->pwm[q][2],reads[num_reads]->pwm[q][3]);
		}
		string consensus = GetConsensus(reads[num_reads]); 
		printf("[%d] Length: %d\tconsensus: %s\n",num_reads,reads[num_reads]->length,consensus.c_str());
		//cout << "consensus: " << consensus << "\t";
		printf("\n\n");
#endif

		num_reads++;
	}
	
	offset = in.tellg();
	in.close();

	return true;
}

/************************************/
/*           FASTQ					*/
/************************************/
bool SeqReader::get_more_fastq() {
	unsigned int J;
	if(gVERBOSE > 1)
		fprintf(stderr,"[%d/-] Reading more reads FASTQ style",iproc);

	clearReads(reads);
	
	if(file_finished)
		return false;
		
	//reads.clear();
	seq_counter = 0;

	ifstream in;
	in.open(this_filename.c_str());
	in.seekg(offset);
	
	if(in.eof()) {
		file_finished = true;
		return false;
	}

	string name, sequence, name_second, fastq;
	
	//READ_BUFFER is the amount of reads we'll read into memory.
	for(unsigned int to_read = 0; to_read < READ_BUFFER; to_read++) {
#ifdef DEBUG
		if(gVERBOSE > 1 && (to_read % 100 == 0 ) )	{
			cout << ".";
			cout.flush();
		}
#endif
		vector< vector<double> > one_read;

		if(in.eof()) {
			file_finished = true;
			in.close();

			if(to_read > 0)
				return true;
			return false;
		}

		getline(in,name);
		
		// Make sure there aren't any blank lines
		while(name=="" && !in.eof()) {
			//fprintf(stderr,"[%d/-] reado: %u Name: >%s<, eof: %c\n",iproc,to_read,name.c_str(),in.eof() ? 'Y':'N');
			getline(in,name);
		}

		// Check to see if we're at the end of the read
		if(in.eof()) {
			file_finished = true;
			in.close();
			if(to_read > 0)
				return true;
			return false;
		}
		getline(in,sequence);
		getline(in,name_second);
		getline(in,fastq);

		// Let's do some checking to see if we can recover from mistakes in fastq files
		while(name[0] != '@' || name_second[0] != '+' || sequence.length() > fastq.length()) {
			fprintf(stderr, "--Found error...atteming to resolve...\n");
			//ERROR
			// First, check to see if the lines are formatted correctly (with the '@' and '+')
			if(name[0] != '@' || name_second[0] != '+') {
				fprintf(stderr,"--ERROR at sequence %s (name[s] formatted incorrectly)",name.c_str());
				fprintf(stderr,"\n\tTrying to recover...\n");
				while( (name[0] != '@' || name_second[0] != '+') && !in.eof() ) {
					// Shift them all, see if we recover.
					name = sequence;
					sequence = name_second;
					name_second = fastq;
					getline(in,fastq);

					// Check to see if we're at the end of the read
					if(in.eof()) {
						file_finished = true;
						in.close();
						if(to_read > 0)
							return true;
						return false;
					}
				}
				if(in.eof()) {
					file_finished = true;
					in.close();
					if(to_read > 0)
						return true;
					return false;
				}

			}

			// Next, check to make sure the sequences are of the right length
			if(sequence.length() > fastq.length()) {
				fprintf(stderr,"--ERROR at sequence %s (length of fastq and sequence not equal)",name.c_str());
				fprintf(stderr,"\n\tTrying to recover...\n");
				// Just read in 4 more.  We're assuming an error.
				getline(in,name);
				getline(in,sequence);
				getline(in,name_second);
				getline(in,fastq);
				
				// Check to see if we're at the end of the read
				if(in.eof()) {
					file_finished = true;
					in.close();
					if(to_read > 0)
						return true;
					return false;
				}
			}

		}
		
		//  if there was an g_adaptor sequence defined, remove it from the end
		//  of this read.
		J=sequence.size();
		if(g_adaptor)
			J=FixReads2(sequence);
		
		for(int i=0; i<J; i++) {
			vector<double> one_char(4);

			int Q = (int)fastq[i];
			double max_prb;
			if(gILLUMINA) {
				// The former ILLUMINA version used to do this.  Not anymore.
				Q -= 64;
				//Q -= 33;
				max_prb = Q2Prb_ill((double)Q);
			}
			else {
				Q -= 33;	// Standard method
				max_prb = Q2Prb_std((double)Q);
			}
				
			double other_prb = (1 - max_prb) / 3;

			// Check for Errors
			if(max_prb < 0) {
				// Can recover if we've subtracted 64 instead of 33
				if(gILLUMINA) {
					fprintf(stderr,"WARNING:  --illumina flag caused undesirable side-effects.\n\tTurning it off\n");
					gILLUMINA = 0;
					//reset this sequence
					one_read.clear();
					i = -1;
					continue;
				}
				fprintf(stderr,"Found Q of %d and original was %d\n",Q,fastq[i]);
				char* errstr = new char[100];
				strcpy(errstr,"Invalid Fastq Character? at sequence: ");
				strcat(errstr,name.c_str());
				errstr[99] = '\0';
				fprintf(stderr,"ERROR STRING: %s",errstr);
				throw errstr;
			}

			//fprintf(stderr,"Total sum: %f+3*%f=%f\n",max_prb,other_prb,max_prb+3*other_prb);

			switch(tolower(sequence[i])) {
				case 'a':
					one_char[A_POS] = max_prb;
					one_char[C_POS] = other_prb;
					one_char[G_POS] = other_prb;
					one_char[T_POS] = other_prb;
					break;
				case 'c':
					one_char[A_POS] = other_prb;
					one_char[C_POS] = max_prb;
					one_char[G_POS] = other_prb;
					one_char[T_POS] = other_prb;
					break;
				case 'g':
					one_char[A_POS] = other_prb;
					one_char[C_POS] = other_prb;
					one_char[G_POS] = max_prb;
					one_char[T_POS] = other_prb;
					break;
				case 't':
					one_char[A_POS] = other_prb;
					one_char[C_POS] = other_prb;
					one_char[G_POS] = other_prb;
					one_char[T_POS] = max_prb;
					break;
				// What do we do for an 'N'?
				// Just ignore the prb and give them all equal
				case 'n':
					one_char[A_POS] = other_prb;
					one_char[C_POS] = other_prb;
					one_char[G_POS] = other_prb;
					one_char[T_POS] = other_prb;
					break;
				default:
					//cout << "ERROR:  Could not determine fastq character! (" << tolower(sequence[i]) << "\n";
					one_char[A_POS] = other_prb;
					one_char[C_POS] = other_prb;
					one_char[G_POS] = other_prb;
					one_char[T_POS] = other_prb;
					break;
			}

#ifdef DEBUG
			if(gVERBOSE > 2)
				cout << sequence[i] << "(" << fastq[i] << ") has\ta: " << one_char[A_POS] << "\tc: "  << one_char[C_POS] << "\tg: " << one_char[G_POS] << "\tt: "  << one_char[T_POS] << endl;
#endif
			one_read.push_back(one_char);
			
			
		}

		
/*
#ifdef DEBUG
		if(gVERBOSE > 3) {
			string consensus = GetConsensus(one_read); 
			cout << "consensus: " << consensus << "\t";
			bin_seq bs;
			cout << "max align: " << bs.get_align_score(one_read,consensus) << endl;
			
		}
#endif
*/
		
		//Read* temp_read = new Read;
		reads[num_reads] = new Read;
		reads[num_reads]->length = J;
		reads[num_reads]->pwm = new float*[J];
		reads[num_reads]->name = new char[name.size()+1];
		strcpy(reads[num_reads]->name,name.c_str()+1);
		reads[num_reads]->seq = sequence;
		reads[num_reads]->fq  = fastq;
		
		for(unsigned int i=0; i<J; i++) {
			reads[num_reads]->pwm[i] = new float[4];
			for(unsigned int j=0; j<4; j++) {
				reads[num_reads]->pwm[i][j] = (float)one_read[i][j];
			}
		}

#ifdef DEBUG
		//for(unsigned int q=0; q<reads[num_reads]->length; q++) {
		//	printf("%f %f %f %f\n",reads[num_reads]->pwm[q][0],reads[num_reads]->pwm[q][1],
		//						   reads[num_reads]->pwm[q][2],reads[num_reads]->pwm[q][3]);
		//}
		string consensus = GetConsensus(reads[num_reads]); 
		printf("[%d] Length: %d\tconsensus: %s\n",num_reads,reads[num_reads]->length,consensus.c_str());
		//cout << "consensus: " << consensus << "\t";
		printf("\n\n");
#endif

		num_reads++;
	}
	
	offset = in.tellg();
	in.close();

	return true;
}

inline float SeqReader::Compare(string consense, string g_adaptor, int pos) { 
	unsigned int num_same=0;

	unsigned int j,i;
	for(i=pos,j=0; j<g_adaptor.length() && i<consense.length(); i++,j++) {
		if(consense[i] == g_adaptor[j])
			num_same++;
	}

	float perc_identity = (float)num_same / j;
	return perc_identity;
}

void SeqReader::FixReads(vector<vector<double> > &m_read) { 
//void SeqReader::FixReads(Read &m_read) { 
	string s_adaptor = g_adaptor;
	string consense = GetConsensus(m_read);
	unsigned int c_len = consense.length();
	unsigned int a_len = s_adaptor.length();

	//if(c_len-a_length < HASH_SIZE)
	//	return;	//don't adapt it if it'll be shorter than the hash size

	//Chop Adaptor length here.
	// we only are worried about the last *i* characters of the sequence
	//int i;
	//for(i=a_len-1; i>0; i--) {
	//	if(consense.substr(c_len-i) == s_adaptor.substr(0,i))
	//		break;
	//}
	unsigned int i;
	float diff=0;
	for(i=0; i<c_len-gMIN_CHOPPED_BASES && i<a_len; i++) {
		diff = Compare(consense,s_adaptor,i);
		if(diff >= gMIN_ADAPTOR_DIFF)
			break;
	}

#ifdef DEBUG
	if(gVERBOSE > 1)
		printf("Chopping read %s at pos %d (chopped %d bases) with diff %f\n",consense.c_str(),i,c_len-i,diff);
#endif
	// if there was a portion of the g_adaptor that matched the read,
	// pop it off the read.
	if( (c_len-i) > gMIN_CHOPPED_BASES) {
		//printf("Fixing read by %d\n",c_len-i); 
		for(unsigned int j=c_len; j>i; j--) {
			m_read.pop_back();
			//delete[] m_read[i];
			//m_read->length--;
		}
	}
	//don't return anything.  Modify the *read* reference.

#ifdef DEBUG
	if(gVERBOSE > 1) {
		printf("              %*s\n",i+(int)strlen(g_adaptor),g_adaptor);
		printf("              %s\n",GetConsensus(m_read).c_str());
	}
#endif

}

unsigned int SeqReader::FixReads2(string seq_read) { 
//void SeqReader::FixReads(Read &m_read) { 
	string s_adaptor = g_adaptor;
	//string consense = GetConsensus(m_read);
	unsigned int c_len = seq_read.length();
	//unsigned int a_len = s_adaptor.length();

	unsigned int i;
	float diff=0;
	for(i=0; i<c_len-gMIN_CHOPPED_BASES; i++) {
		diff = Compare(seq_read,s_adaptor,i);
		if(diff >= gMIN_ADAPTOR_DIFF)
			break;
	}

return i;
}

string SeqReader::GetConsensus(Read* read) {
	vector<vector<double> > pwm;
	for(unsigned int i=0; i<read->length; i++) {
		vector<double> one_chr;
		for(unsigned int j=0; j<4; j++) {
			one_chr.push_back(read->pwm[i][j]);
		}
		pwm.push_back(one_chr);
	}
	
	return GetConsensus(pwm);
}

string SeqReader::GetConsensus(vector<vector<double> > &read) {
//string SeqReader::GetConsensus(Read &pwm) {
	string seq = "";

	for(unsigned int i=0; i<read.size(); i++) {
		seq += max_char(read[i]);
	//for(unsigned int i=0; i<read.length; i++) {
	//	seq += max_char(read.pwm[i]);
	}

	return seq;
}

Read* SeqReader::getReadFromStrings(char* name, string &consensus, string &qual) {
	Read* n = new Read;
	n->name = new char[strlen(name)+1];
	strcpy(n->name,name);
	n->length = qual.size();
	n->pwm = new float*[qual.size()];
	
	for(unsigned int i=0; i<consensus.size(); i++) {
		n->pwm[i] = new float[4];
	
		int Q = (int)qual[i];
		Q -= 33;
		double max_prb = Q;

		// Check for Errors
		if(max_prb < 0) {
			fprintf(stderr,"Found Q of %d and original was %d\n",Q,qual[i]);
			char* errstr = new char[100];
			strcpy(errstr,"Invalid Fastq Character? at sequence: ");
			strcat(errstr,name);
			errstr[99] = '\0';
			fprintf(stderr,"ERROR STRING: %s",errstr);
			throw new Exception(errstr);
		}

		max_prb = Q2Prb_ill(max_prb);
		double other_prb = (1 - max_prb) / 3;

		switch(tolower(consensus[i])) {
			case 'a':
				n->pwm[i][A_POS] = max_prb;
				n->pwm[i][C_POS] = other_prb;
				n->pwm[i][G_POS] = other_prb;
				n->pwm[i][T_POS] = other_prb;
				break;
			case 'c':
				n->pwm[i][A_POS] = other_prb;
				n->pwm[i][C_POS] = max_prb;
				n->pwm[i][G_POS] = other_prb;
				n->pwm[i][T_POS] = other_prb;
				break;
			case 'g':
				n->pwm[i][A_POS] = other_prb;
				n->pwm[i][C_POS] = other_prb;
				n->pwm[i][G_POS] = max_prb;
				n->pwm[i][T_POS] = other_prb;
				break;
			case 't':
				n->pwm[i][A_POS] = other_prb;
				n->pwm[i][C_POS] = other_prb;
				n->pwm[i][G_POS] = other_prb;
				n->pwm[i][T_POS] = max_prb;
				break;
			// What do we do for an 'N'?
			// Just ignore the prb and give them all equal
			case 'n':
				n->pwm[i][A_POS] = .25;
				n->pwm[i][C_POS] = .25;
				n->pwm[i][G_POS] = .25;
				n->pwm[i][T_POS] = .25;
				break;
			default:
				n->pwm[i][A_POS] = .25;
				n->pwm[i][C_POS] = .25;
				n->pwm[i][G_POS] = .25;
				n->pwm[i][T_POS] = .25;
				break;
		}
	
	}
	
	return n;
}

/*
 * Because this is used to produce the consensus sequence, we need ambiguity characters...
 */
inline char SeqReader::max_char(vector<double> &chr) {
//inline char SeqReader::max_char(float chr[4]) {
	if((chr[0] == chr[1]) && (chr[0] == chr[2]) && (chr[0] == chr[3]))
		return 'n';
	if(chr[0] >= chr[1]) {
		if(chr[0] >= chr[2]) {
			if(chr[0] >= chr[3])
				return 'a';
			else
				return 't';
			}
		else {
			if(chr[2] >= chr[3])
				return 'g';
			else
				return 't';
		}
	}
	else {
		if(chr[1] >= chr[2]) {
			if(chr[1] >= chr[3])
				return 'c';
			else
				return 't';
		}
		else {
			if(chr[2] >= chr[3])
				return 'g';
			else return 't';
		}
	}

	return 'n';
}

unsigned int SeqReader::GetNumSeqs() {
	ifstream inf(this_filename.c_str());
	//inf.clear();
	
	unsigned int line_count=0;
	
	char temp[gBUFFER_SIZE];
	
	switch(this_type) {
	case PRB:
	case INT:
		while(inf) {
			inf.getline(temp, gBUFFER_SIZE);
			if(inf.eof())
				break;
			line_count++;
		}
		break;
	case FASTQ:
		while(inf) {
			inf.getline(temp, gBUFFER_SIZE);
			if(inf.eof())
				break;
			line_count++;
		}
		line_count /= 4;

		break;

	case FASTA:
		while(inf) {
			inf.getline(temp,gBUFFER_SIZE);
			if(inf.eof())
				break;
			if(temp[0] == '>')
				line_count++;
		}
		break;
	default:
		line_count = 0;
		break;
	}
	
	inf.close();
		
	return line_count;
}

unsigned int SeqReader::GetNumSeqs(const string& str) {
	SeqType this_type = find_type(str);

	ifstream inf(str.c_str());
	inf.clear();
	
	unsigned int line_count=0;
	
	char temp[gBUFFER_SIZE];
	
	switch(this_type) {
	case PRB:
	case INT:
		while(inf) {
			inf.getline(temp, gBUFFER_SIZE);
			if(inf.eof())
				break;
			line_count++;
		}
		break;
	case FASTQ:
		while(inf) {
			inf.getline(temp, gBUFFER_SIZE);
			if(inf.eof())
				break;
			line_count++;
		}
		line_count /= 4;

		break;

	case FASTA:
		while(inf) {
			inf.getline(temp,gBUFFER_SIZE);
			if(inf.eof())
				break;
			if(temp[0] == '>')
				line_count++;
		}
		break;
	default:
		fprintf(stderr,"ERROR: DEFAULT!\n");
		line_count = 0;
		break;
	}
	
	inf.close();
		
	return line_count;
}

unsigned int SeqReader::GetNumSeqFiles() {
	return seq_files.size();
}

vector<string>& SeqReader::GetSeqFileVector() {
	return seq_files;
}
		
bool SeqReader::Test(ostream &os, unsigned int &warnings) {
	bool success = true;

	try {
		SeqReader sr;
		sr.use("test/shortall_prb.txt");

		unsigned int i=0;
		//while(sr.GetNextSequence().size() > 0)
		while(sr.GetNextSequence())
			i++;
		cout << "i=" << i << endl;
		TEST(i == sr.GetNumSeqs());
		
		sr.use("test/shortall.fa");
		i=0;
		//while(sr.GetNextSequence().size() > 0)
		while(sr.GetNextSequence())
			i++;
		cout << "I: " << i << endl;
		TEST(i == sr.GetNumSeqs());
		
		sr.use("test/shortall_seqs.fq");
		i=0;
		cout << sr.GetNumSeqs() << " vs ";
		while(sr.GetNextSequence())
			i++;
		cout << "i=" << i << endl;
		TEST(i == sr.GetNumSeqs());
	}
	catch(const char* error) {
		os << "\t********Error! Did you include the proper test files?" << endl;
	}
	catch(...) {
		os << "\t********Error! Did you include the proper test files?" << endl;
	}
	return success;
}

