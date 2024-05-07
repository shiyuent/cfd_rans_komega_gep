CC=             gcc
CFLAGS=         -g
COBJECTS=       program.o tinyexpr.o nrutil.o tridag.o polcoe.o polint.o
LIBS=           -lm

kom.x:     $(COBJECTS) 
	cc -o kom.x  $(COBJECTS)  $(LIBS) 

$(COBJECTS):    nrutil.h 

clean:
	/bin/rm *.o
	/bin/rm kom.x
