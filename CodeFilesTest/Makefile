SHELL = /bin/sh
RM = rm -f
TAR = tar zcvf
CC = gcc -llapack -lblas -lm -L/usr/local/linalg/lib/
CFLAGS = -O

all:	DriverTest Makefile

eigenDriver:	DriverTest.o HiggsLoops.o

eigenDriver.o:	DriverTest.c HiggsLoops.o

eigenVec.o:	HiggsLoops.c HiggsLoops.h

run:	DriverTest
	./DriverTest

clean:
	$(RM) DriverTest DriverTest.o HiggsLoops.o 
	$(RM) *~

tar:	clean
	cd .. ; $(TAR) HiggsLoops.tgz HiggsLoops/

