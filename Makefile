#CC=clang-3.8
CC=gcc
CFLAGS=-fopenmp --std=gnu99 -Wall -Werror -Ofast -g
#CFLAGS=-fopenmp --std=gnu99 -Wall -Werror -O0 -g -DDEBUG
#CFLAGS=-fopenmp --std=gnu99 -Wall -Werror -O0 -g
GIT_COMMIT=$(shell ./git_status.sh)
#TODO: add "-dirty" if the working dir has local changes

LIBS=-lrt -lm
DEPS=sim.h network.h tech.h
OBJ=main.o sim.o network.o tech.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) -DGIT_COMMIT=\"$(GIT_COMMIT)\"

sim: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: all
all: sim

.PHONY: clean
clean:
	rm *.o
	rm -f sim
