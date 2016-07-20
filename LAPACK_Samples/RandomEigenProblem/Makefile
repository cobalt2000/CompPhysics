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
	$(CPP) $(CFLAGS) -o RandomEigenProblem RandomEigenProblem.o eigenVec.o $(CLIBS) $(CPPLIBS)

RandomEigenProblem.o:	RandomEigenProblem.cpp eigenVec.o
	$(CPP) $(CPPFLAGS) -c RandomEigenProblem.cpp  $(CPPLIBS)

eigenVec.o:	eigenVec.c eigenVec.h
	$(CC) $(CFLAGS) -c eigenVec.c $(CLIBS)

run:	all
	mpiexec -np 2 ./RandomEigenProblem

clean:
	$(RM) RandomEigenProblem eigenVec.o RandomEigenProblem.o
	$(RM) output*
	$(RM) *~

tar:	clean
	cd .. ; $(TAR) Sprng_mpi.tgz sprng_mpi/

