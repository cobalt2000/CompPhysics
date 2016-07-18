SHELL = /bin/sh
RM = rm -f
TAR = tar zcvf
CC = mpic++
CFLAGS = -O3 -I$(HOME)/sprng5/include
CLIBS = -L$(HOME)/sprng5/lib -lsprng

all:	simple Makefile

simple:	simple.cpp
	$(CC) $(CFLAGS) -o simple simple.cpp $(CLIBS)

run:	simple
	mpiexec -np 1 ./simple

clean:
	$(RM) simple
	$(RM) *~

tar:	clean
	cd .. ; $(TAR) Simple.tgz simple/

