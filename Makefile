# Makefile
#
# Copyright (C) 2009-2014 Nathan Clement
#
# This file is part of GNUMAP.
#
# GNUMAP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GNUMAP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNUMAP.  If not, see <http://www.gnu.org/licenses/>.

PLAIN_EXE_NAME = bin/gnumap-plain
MPI_EXE_NAME = bin/gnumap
TEST_BIN_FILE = bin/gnutest


EXE_OBJ_FILES = obj/Driver.o $(OBJ_FILES) \
		obj/NormalScoredSeq.o obj/BSScoredSeq.o obj/SNPScoredSeq.o
OBJ_FILES = obj/centers.o obj/bin_seq.o obj/Reader.o obj/SeqReader.o obj/Genome.o $(BWT_OBJ_FILES)
		
INC_FILES = inc/const_include.h inc/const_define.h inc/Exception.h inc/SeqManager.h inc/gvector.h \
			inc/ScoredSeq.h

CONV_EXE_NAME = bin/sam2sgr

CONV_OBJ_FILES = obj/sam2sgr.o obj/Genome.o obj/Reader.o obj/bin_seq.o obj/SeqReader.o \
		obj/NormalScoredSeq.o obj/BSScoredSeq.o obj/SNPScoredSeq.o $(BWT_OBJ_FILES)

CONV_INC_FILES = inc/const_include.h inc/const_define.h

BWT_OBJ_FILES = obj/GenomeBwt.o obj/bwt.o obj/utils.o obj/bntseq.o obj/bwtindex.o \
					 obj/is.o obj/bwt_gen.o obj/QSufSort.o

BATCH_CONS_EXE = bin/sam2consensus
BATCH_CONS_OBJ_FILES = obj/sam2consensus.o obj/Genome.o obj/Reader.o obj/bin_seq.o $(BWT_OBJ_FILES)
BATCH_CONS_INC_FILES = inc/const_include.h inc/const_define.h

GSL_LIB_FILE = lib/lib/libgsl.a
GSL_LIB_DIR = lib/gsl-1.9/
GSL_ZIP_FILE = lib/gsl-1.9.tar.gz
		
TEST_OBJ_FILES = obj/TestDriver.o $(OBJ_FILES)

DEBUG_FLAGS = -m64 -Wall -g -rdynamic

EXTRA_FLAGS =

OPT_FLAGS = -m64 -O3
ERR_FLAGS = -Wall

FLAGS = $(OPT_FLAGS) -DGENOME_BWT -std=c++0x
FLAGS = -DDEBUG_NW -DDEBUG_TIME $(DEBUG_FLAGS) -std=c++0x
#FLAGS = $(PROFILE_FLAGS)
#FLAGS = $(DEBUG_FLAGS) -DDEBUG
#FLAGS = $(DEBUG_FLAGS) -DSET_POS
#FLAGS = $(OPT_FLAGS) -DSET_POS -D_INDEL
#FLAGS = $(OPT_FLAGS) -Wall -DSET_POS $(PROFILE_FLAGS)
#FLAGS = $(OPT_FLAGS) -Wall -DSET_POS $(OMP_FLAGS)

PROFILE_FLAGS = -m64 -pg -O3
#FLAGS= $(PROFILE_FLAGS)

#Here's how to compile with mpic++
MPIXX = mpic++ -DMPI_RUN
MPICC = mpicc -DMPI_RUN
# Intel's MPI compiler
#MPIXX = mpiCC
OMP_FLAGS=-DOMP_RUN -fopenmp


MPI_EXIST = $(shell mpic++ -v 2>&1 | grep -o version)
ifeq ($(strip $(MPI_EXIST)),)
    GXX=g++
	GCC=gcc
    BUILDTARGET=plain
	BUILDEXE=$(PLAIN_EXE_NAME)
else
    GXX=$(MPIXX)
	GCC=$(MPICC)
    BUILDTARGET=mpi
	BUILDEXE=$(MPI_EXE_NAME)
endif

GXX=g++
#GXX=icpc
BUILDTARGET=plain
BUILDEXE=$(PLAIN_EXE_NAME)

INC = -Iinc/ -I$(GSL_LIB_DIR)
# For some reason, I thought we needed to do dynamic linking.  Running a few tests,
# it doesn't seem we do after all, so we'll just pull this out, but leave it in just in case.
LIB = -lz -lm -dynamic -lpthread -Llib/lib -Wl,-Bstatic -lgsl -lgslcblas -Wl,-Bdynamic
#LIB =  -Llib/lib -lgsl -lgslcblas -dynamic -lpthread

prog : $(BUILDTARGET)

all : $(BUILDTARGET) conv batch-consensus

mpi : $(GSL_LIB_FILE) $(MPI_EXE_NAME) 
	-@ echo ""; echo "MPI Build Successful"; echo ""
	
plain : $(GSL_LIB_FILE) $(PLAIN_EXE_NAME) 
	-@ echo ""; echo "Successful"; echo ""

conv : $(CONV_EXE_NAME)
	-@ echo ""; echo "bin/sam2sgr Build Successful"; echo ""

batch-consensus : $(BATCH_CONS_EXE)
	-@ echo ""; echo "bin/sam2consensus Build Successful"; echo ""

example : $(BUILDEXE)
	$(BUILDEXE) -g examples/Cel_gen.fa -o gnumap.output -a .9 -v 1 \
	examples/example_sequences_prb.txt --illumina

example-threaded : $(BUILDEXE)
	$(BUILDEXE) -g examples/Cel_gen.fa -o gnumap.output -a .9 -v 1 -c 8 \
	examples/example_sequences_prb.txt

example-snp : $(BUILDEXE)
	$(BUILDEXE) -g examples/Cel_gen.fa -o example_snp.output -a .9 -v 1 \
	--snp -j 5 examples/example_sequences_prb.txt

test : test_bin 
	./$(TEST_BIN_FILE)

test_bin : $(TEST_OBJ_FILES)
	$(GXX) $(FLAGS) -o $(TEST_BIN_FILE) $(TEST_OBJ_FILES) $(INC) $(LIB)
	
# -pg flag to compile AND link 
# run program to generate gmon.out 
# run gprof to print reports 
#  gprof -p bin/gnumap gmon.out # print flat profile 
#  gprof -q bin/gnumap gmon.out # print call graph 
profile : $(BUILDTARGET) example-threaded
	gprof -p $(BUILDEXE) gmon.out > gprof-p.txt
	gprof -q $(BUILDEXE) gmon.out > gprof-q.txt

$(PLAIN_EXE_NAME) : $(EXE_OBJ_FILES) $(INC_FILES) inc/ScoredSeq.h inc/SeqManager.h
	$(GXX) $(FLAGS) -o $(PLAIN_EXE_NAME) $(EXE_OBJ_FILES) $(INC) $(LIB) $(EXTRA_FLAGS)

bin/gnumap-plain-debug : $(EXE_OBJ_FILES) $(INC_FILES) inc/ScoredSeq.h inc/SeqManager.h
	$(GXX) $(FLAGS) -o bin/gnumap-plain-debug $(EXE_OBJ_FILES) $(INC) $(LIB) $(EXTRA_FLAGS)

$(MPI_EXE_NAME) : $(EXE_OBJ_FILES) $(INC_FILES) inc/ScoredSeq.h inc/SeqManager.h
	$(MPIXX) $(FLAGS) -o $(MPI_EXE_NAME) $(EXE_OBJ_FILES) $(INC) $(LIB) $(EXTRA_FLAGS)
	
obj/Driver.o : src/Driver.cpp $(INC_FILES) inc/a_matrices.c inc/SequenceOperations.h \
				inc/align_seq2_raw.cpp
	$(GXX) $(FLAGS) -o $@ -c $< $(INC)

obj/TestDriver.o : test/TestDriver.cpp
	$(GXX) $(FLAGS) -o obj/TestDriver.o -c test/TestDriver.cpp $(INC)

$(CONV_EXE_NAME) : $(CONV_OBJ_FILES) $(CONV_INC_FILES)
	$(GXX) $(FLAGS) -o $(CONV_EXE_NAME) $(CONV_OBJ_FILES) $(INC) $(LIB) $(EXTRA_FLAGS) -fopenmp
obj/sam2sgr.o : src/sam2sgr.cpp inc/a_matrices.c
	$(GXX) $(FLAGS) -o $@ -c $< $(INC) -fopenmp

$(BATCH_CONS_EXE) : $(BATCH_CONS_OBJ_FILES) $(BATCH_CONS_INC_FILES)
	$(GXX) $(FLAGS) -o $(BATCH_CONS_EXE) $(BATCH_CONS_OBJ_FILES) $(INC) $(LIB) $(EXTRA_FLAGS)
obj/sam2consensus.o : src/sam2consensus.cpp inc/a_matrices.c
	$(GXX) $(FLAGS) -o $@ -c $< $(INC)

obj/GenomeMem.o : src/GenomeMem.cpp inc/GenomeMem.h
	$(GXX) $(FLAGS) -o $@ -c $< $(INC)

obj/utils.o : src/utils.c inc/utils.h
	$(GXX) $(FLAGS) -o $@ -c $< $(INC)

obj/bwt.o : src/bwt.c inc/bwt.h
	$(GXX) $(FLAGS) -o $@ -c $< $(INC)

obj/bntseq.o : src/bntseq.c inc/bntseq.h
	$(GXX) $(FLAGS) -o $@ -c $< $(INC)

obj/GenomeBwt.o : src/GenomeBwt.cpp inc/GenomeBwt.h
	$(GXX) $(FLAGS) -o $@ -c $< $(INC)

obj/bwtindex.o : src/bwtindex.c
	$(GCC) $(FLAGS) -o $@ -c $< $(INC)

obj/is.o : src/is.c
	$(GCC) $(FLAGS) -o $@ -c $< $(INC)

obj/bwt_gen.o : src/bwt_gen.c
	$(GCC) $(FLAGS) -o $@ -c $< $(INC)

obj/QSufSort.o : src/QSufSort.c inc/QSufSort.h
	$(GCC) $(FLAGS) -o $@ -c $< $(INC)

obj/%.o : src/%.cpp inc/%.h
	$(GXX) $(FLAGS) -o $@ -c $< $(INC)

#	valgrind --tool=memcheck --leak-check=full --show-reachable=yes --suppressions=test/string.supp $(EXE_NAME) \
#	-g /data/public/hci/mouse/chr1_both.seq -o test/val_test \
#	-a .7 -p -c 1 -v 1 test/multiple_files_prb.txt
#	-a .7 -p -c 1 -v 1 test/very_short_prb.txt

val : bin
	valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes --suppressions=mpi.supp -v \
	$(BUILDEXE) \
	-g "test/shortall.fa" -o val_test \
	-a .9 -p -v 1 test/shortall_prb.txt
#	valgrind --tool=memcheck --leak-check=full --show-reachable=yes --suppressions=test/string.supp -v \
#	$(EXE_NAME) \
#	-g ~/SeqTrack/data/evan/snomiRNA_2006_noneg.seq -o test/val_test \
#	-a .7 -p -c 8 -v 1 test/multiple_files_prb.txt

val-snp : bin
	valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes -v \
	$(BUILDEXE) \
	-g "test/shortall.fa" -o val_test --snp \
	-a .9 -p -v 1 test/shortall_prb.txt

# Compile the gsl library if we need to
$(GSL_LIB_FILE) :
	if [ ! -d $(GSL_LIB_DIR) ]; then \
		tar -xzf $(GSL_ZIP_FILE) -C lib/; \
	fi
	cd $(GSL_LIB_DIR) && \
	if [ ! -e Makefile ]; then \
		echo "Couldn't find $(GSL_LIB_DIR)Makefile"; \
		./configure --prefix="$(CURDIR)/lib/"; \
	fi && \
	make && make install \

clean : 
	-rm -f obj/*
	-rm -f $(BUILDEXE)
	-rm -f $(TEST_BIN_FILE)
	-rm -f $(BATCH_CONS_EXE)
	-rm -f $(CONV_EXE_NAME)
	
deep-clean: clean
	-rm -rf lib/lib
	-rm -rf lib/gsl-1.9
	-rm -rf lib/bin/
	-rm -rf lib/share/
	-rm -rf lib/include/
