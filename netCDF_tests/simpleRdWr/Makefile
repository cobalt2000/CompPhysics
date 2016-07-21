SHELL = /bin/sh
RM = rm -f
TAR = tar zcvf
CC = gcc 
CFLAGS = -O -I/usr/local/netcdf/include
LINK = -L/usr/local/netcdf/lib -lnetcdf

all:	Makefile simple_nc4_wr simple_nc4_rd

simple_nc4_wr:	simple_nc4_wr.c
	$(CC) $(CFLAGS) simple_nc4_wr.c $(LINK) -o write.exe

simple_nc4_rd:	simple_nc4_rd.c
	$(CC) $(CFLAGS) simple_nc4_rd.c $(LINK) -o read.exe

run:	main
	./main

clean:
	$(RM) read.exe write.exe simple_nc4.nc
	$(RM) *~

tar:	clean
	cd .. ; $(TAR) NETcdf4_sample.tgz simple_nc4_rd
