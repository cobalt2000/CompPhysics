SHELL = /bin/sh
RM = rm -f
TAR = tar zcvf
CC = gcc 
CFLAGS = -O

all:	DriverTest Makefile

DriverTest:	DriverTest.o HiggsLoops.o

DriverTest.o:	DriverTest.c HiggsLoops.o

HiggsLoops.o:	HiggsLoops.c HiggsLoops.h

run:	DriverTest
	./DriverTest

clean:
	$(RM) DriverTest DriverTest.o HiggsLoops.o 
	$(RM) *~

tar:	clean
	cd .. ; $(TAR) HiggsLoops.tgz HiggsLoops/

