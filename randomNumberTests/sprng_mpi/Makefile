SHELL = /bin/sh
RM = rm -f
TAR = tar zcvf
CC = mpic++
CFLAGS = -O3 -I$(HOME)/sprng5/include
CLIBS = -L$(HOME)/sprng5/lib -lsprng

all:	sprng_mpi Makefile

sprng_mpi:	sprng_mpi.cpp
	$(CC) $(CFLAGS) -o sprng_mpi sprng_mpi.cpp $(CLIBS)

run:	sprng_mpi
	mpiexec -np 2 ./sprng_mpi

clean:
	$(RM) sprng_mpi
	$(RM) *~

tar:	clean
	cd .. ; $(TAR) Sprng_mpi.tgz sprng_mpi/

