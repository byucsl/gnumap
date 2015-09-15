#! /bin/bash

echo g++ -DMPI_RUN -m64 -O3 \
	-o bin/gnumap obj/Driver.o obj/centers.o obj/bin_seq.o obj/Reader.o obj/SeqReader.o obj/GenomeMem.o obj/NormalScoredSeq.o obj/BSScoredSeq.o obj/SNPScoredSeq.o \
	-Iinc/ -Ilib/gsl-1.14/ -Llib/lib -lgsl -lgslcblas \
	-static -pthread -I/bluescr/nclement/include -L/bluescr/nclement/lib \
		-lmpi_cxx -lmpi -lrdmacm \
		-Wl,--whole-archive -libverbs -Wl,--no-whole-archive \
		-lrt -Wl,--export-dynamic -lnsl -lutil -lm -ltorque -ldl
g++ -DMPI_RUN -m64 -O3 \
	-o bin/gnumap obj/Driver.o obj/centers.o obj/bin_seq.o obj/Reader.o obj/SeqReader.o obj/GenomeMem.o obj/NormalScoredSeq.o obj/BSScoredSeq.o obj/SNPScoredSeq.o \
	-Iinc/ -Ilib/gsl-1.14/ -Llib/lib -lgsl -lgslcblas \
	-static -pthread -I/bluescr/nclement/include -L/bluescr/nclement/lib \
		-lmpi_cxx -lmpi -lrdmacm \
		-Wl,--whole-archive -libverbs -Wl,--no-whole-archive \
		-lrt -Wl,--export-dynamic -lnsl -lutil -lm -ltorque -ldl
