#CC=clang-3.8
CC=gcc
CFLAGS=--std=gnu99 -Wall -pedantic -Werror -g -pg
RELFLAGS=-Ofast
DEBUGFLAGS=-DDEBUG -no-pie -pg -O0
GIT_COMMIT=$(shell ./git_status.sh)

LIBS=-lrt -lm
DEPS=sim.h print.h command.h network.h arch.h
OBJ=main.o sim.o command.o network.o arch.o
DEBUGDIR=debug
RELDIR=release

RELOBJ=$(addprefix $(RELDIR)/,$(OBJ))
RELEXE=$(RELDIR)/sim
DEBUGOBJ=$(addprefix $(DEBUGDIR)/,$(OBJ))
DEBUGEXE=$(DEBUGDIR)/sim

.PHONY: all sim release debug clean prep

all: prep sim

sim: release
	cp $(RELEXE) sim

release: $(RELEXE)

debug: prep $(DEBUGEXE)

$(RELDIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(RELFLAGS) -DGIT_COMMIT=\"$(GIT_COMMIT)\"

$(RELEXE): $(RELOBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(RELFLAGS) $(LIBS)

$(DEBUGDIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(DEBUGFLAGS) -DGIT_COMMIT=\"$(GIT_COMMIT)\"

$(DEBUGEXE): $(DEBUGOBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(DEBUGFLAG) $(LIBS)

clean:
	rm -f $(RELEXE)
	rm -f $(DEBUGEXE)
	rm -f $(RELOBJ)
	rm -f $(DEBUGOBJ)

prep:
	@mkdir -p $(RELDIR) $(DEBUGDIR)
