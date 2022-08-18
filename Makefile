#CC=clang-3.8
CC=gcc
CFLAGS=--std=gnu99 -Wall -pedantic -Werror -g
GIT_COMMIT=$(shell ./git_status.sh)

LIBS=-lrt -lm
DEPS=sim.h print.h command.h network.h arch.h
OBJ=main.o sim.o command.o network.o arch.o
OBJDIR=obj

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) -DGIT_COMMIT=\"$(GIT_COMMIT)\"

sim: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: all
all: CFLAGS += -Ofast
all: sim

.PHONY: debug
debug: CFLAGS += -DDEBUG -no-pie -pg -O0
debug: sim

.PHONY: clean
clean:
	rm *.o
	rm -f sim
