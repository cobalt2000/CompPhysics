SHELL = /bin/sh
RM = rm -f
TAR = tar zcvf
CC = gcc -lm 
CFLAGS = -O

all:	b2sgDriver Makefile

b2sgDriver:	b2sgDriver.o Particle_struct.o b_to_sgamma.o

b2sgDriver.o:	b2sgDriver.c

Particle_struct.o:	Particle_struct.c Particle_struct.h

b_to_sgamma.o:	b_to_sgamma.c b_to_sgamma.h

run:	DriverTest
	./DriverTest

clean:
	$(RM) DriverTest DriverTest.o HiggsLoops.o 
	$(RM) *~

tar:	clean
	cd .. ; $(TAR) HiggsLoops.tgz HiggsLoops/

