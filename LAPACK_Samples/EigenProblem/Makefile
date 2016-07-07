SHELL = /bin/sh
RM = rm -f
TAR = tar zcvf
CC = gcc -llapack -lblas -lm -L/usr/local/linalg/lib/
CFLAGS = -O

all:	eigenDriver Makefile

eigenDriver:	eigenDriver.o eigenVec.o

eigenDriver.o:	eigenDriver.c eigenVec.o

eigenVec.o:	eigenVec.c eigenVec.h

run:	eigenDriver
	./eigenDriver

clean:
	$(RM) eigenDriver eigenDriver.o eigenVec.o 
	$(RM) *~

tar:	clean
	cd .. ; $(TAR) EigenProblem.tgz EigenProblem/

