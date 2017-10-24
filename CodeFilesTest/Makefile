SHELL = /bin/sh
RM = rm -f
TAR = tar zcvf
CC = gcc -lm 
CFLAGS = -O3

all:	DriverTest Makefile

DriverTest:	DriverTest.o HiggsLoops.o

DriverTest.o:	DriverTest.c HiggsLoops.o

HiggsLoops.o:	HiggsLoops.c HiggsLoops.h

b2sgDriver:	b2sgDriver.o b_to_sgamma.o

b2sgDriver.o:	b2sgDriver.c b_to_sgamma.o

b_to_sgamma.o:	b_to_sgamma.c b_to_sgamma.h

run:	DriverTest
	./DriverTest

b2sgrun:	b2sgDriver Makefile
		./b2sgDriver

clean:
	$(RM) DriverTest DriverTest.o HiggsLoops.o 
	$(RM) *~

tar:	clean
	cd .. ; $(TAR) HiggsLoops.tgz HiggsLoops/

