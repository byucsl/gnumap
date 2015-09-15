#! /bin/bash

make clean
make -j 4 -f Makefile.stl
make clean
make -j 4 -f Makefile.stl.str
make clean
make -j 4 all
