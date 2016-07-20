SHELL = /bin/sh
RM = rm -f
TAR = tar zcvf
CPP = mpic++
CC = mpicc
CPPFLAGS = -O3 -I$(HOME)/sprng5/include
CPPLIBS = -L$(HOME)/sprng5/lib -lsprng
CFLAGS = -llapack -lblas -lm
CLIBS = -L/usr/local/linalg/lib/

all:	RandomEigenProblem Makefile

RandomEigenProblem:	RandomEigenProblem.o eigenVec.o
	$(CPP) $(CFLAGS) -o RandomEigenProblem RandomEigenProblem.o eigenVec.o $(CLIBS)

RandomEigenProblem.o:	RandomEigenProblem.cpp eigenVec.o
	$(CPP) $(CPPFLAGS) -c RandomEigenProblem.cpp eigenVec.o $(CPPLIBS)

eigenVec.o:	eigenVec.c eigenVec.h
	$(CC) $(CFLAGS) -c eigenVec.c $(CLIBS)

run:	all
	mpiexec -np 2 ./sprng_mpi

clean:
	$(RM) RandomEigenProblem
	$(RM) *~

tar:	clean
	cd .. ; $(TAR) Sprng_mpi.tgz sprng_mpi/

