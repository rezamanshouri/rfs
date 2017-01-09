CC=g++
CFLAGS=-O2 -g -std=c++11
OMPFLAGS=-fopenmp
C64FLAGS=$(CFLAGS)
LFLAGS=
DEBUGFLAGS=-g -O0
PROFILEFLAGS=-pg
OBJS=my_rfs_supertree
all: $(OBJS)

my_rfs_supertree: my_rfs_supertree.cpp *.h
	$(CC) $(CFLAGS) -o my_rfs_supertree my_rfs_supertree.cpp


.PHONY: test
.PHONY: debug
.PHONY: hyb
.PHONY: profile

debug:
	$(CC) $(LFLAGS) $(DEBUGFLAGS) -o my_rfs_supertree my_rfs_supertree.cpp
profile:
	$(CC) $(LFLAGS) $(DEBUGFLAGS) $(PROFILEFLAGS) -o my_rfs_supertree my_rfs_supertree.cpp
omp:
	$(CC) $(CFLAGS) $(OMPFLAGS) -o my_rfs_supertree-omp my_rfs_supertree.cpp
omp-debug:
	$(CC) $(LFLAGS) $(DEBUGFLAGS) $(OMPFLAGS) -o my_rfs_supertree-omp my_rfs_supertree.cpp
