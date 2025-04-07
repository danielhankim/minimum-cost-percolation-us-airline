OT=-ffast-math -O3
OP=$(OT)
LOP=-lm -o

CHDR=./
CSRC=./
MSRC=./

MAIN=min_cost_percolation

EXE=$(MAIN)
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
