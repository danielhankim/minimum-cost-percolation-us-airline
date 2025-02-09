# options for icc
#    axW, xW: Xeon,Pentium 4, Pentium M
#    axi, xi: Pentium pro, pentium II
#    axK, xK: pentium III
#    tpp7: Pentium 4, Xeon
#    tpp6: Pentium II, III, Pro
#    mcpu=pentiumpro: Pentium Pro, Pentium II, Pentium III
#    mcpu=pentium4: Pentium 4 (default)
#CC=gcc
OT=-O3 -tpp7 -axW -xW -unroll0
OT=-ffast-math -g -O3
OP=$(OT)
LOP=-lm -o
#LOP=-lm -o

# options for mpi
#CC=hcc
#CC=mpicc
#OT=-O3 -tpp6 -axi -xi -unroll0
#OP=$(OT)
#LOP=-lm -lmpi -o

# Options for gcc
#CC=gcc #-march=native
#OT=-O3 -funroll-all-loops -fexpensive-optimizations 
#OP=$(OT)
#LOP=-lm -lgsl -lgslcblas -o

# For debugging with x86
#CC=gcc
#OT=-g -Wall -rdynamic
#OP= $(OT)
#LOP=-lm -o
# For profiling with x86
#CC=gcc
#OT=-pg -Wall
#OP=$(OT)
#LOP=-lm -o

# For debugging on SUNOS
#CC=gcc 
#OT=-g -Wall 
#LOP=-lm -o

# For debugging with Efence
#CC=gcc -march=i686
#OT=-g -Wall 
#OP=$(OT)
#LOP=-lm -lefence -o
##LOP=-L/afs/nd.edu/user26/syook/util/local/lib -lm -lefence -o

# For SUN
#CC=cc
#OT=-fast -xunroll=0 xtarget=ultra2 -xarch=v8 -xO3
####OT=-xunroll=0 -xO3
#OP=$(OT)
#LOP=-lm -o

CHDR=./
CSRC=./
MSRC=./

MAIN=min_cost_percolation
VER=
#SIZE=L128

EXE=$(MAIN)$(VER)$(SIZE)
TAG=$(EXE).out

OBJ=$(MAIN).o mt19937-64.o basic_functions.o binary_search_tree.o percolation_functions.o heap.o

$(TAG)	: $(OBJ)
	$(CC) $(OBJ) $(LOP) $(TAG)
	rm *.o

$(MAIN).o :
	$(CC) $(OP) -I$(CHDR) -c $(MSRC)$(MAIN).c

mt19937-64.o	:
	$(CC) $(OP) -I$(CHDR) -c $(CSRC)mt19937-64.c

basic_functions.o	:
	$(CC) $(OP) -I$(CHDR) -c $(CSRC)basic_functions.c

binary_search_tree.o	:
	$(CC) $(OP) -I$(CHDR) -c $(CSRC)binary_search_tree.c

percolation_functions.o	:
	$(CC) $(OP) -I$(CHDR) -c $(CSRC)percolation_functions.c

heap.o	:
	$(CC) $(OP) -I$(CHDR) -c $(CSRC)heap.c

clear	:
	rm data/*; rm index.dat; cp work/words.lst ./

clean	:
	rm *.o
