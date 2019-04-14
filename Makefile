CC= gcc
CFLAGS= -g -Wall -O2
OBJS= main.o parse_bam.o randl_list.o exp_fit.o
INC= -Lhtslib -lhts -lpthread -lz -lm -I./gsl/ 
LIB= ./gsl/.libs/libgsl.a ./gsl/cblas/.libs/libgslcblas.a

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@ $(INC)

all: ugat

ugat: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(INC) $(LIB)

submodules:
	cd htslib; make libhts.a; cd ..
	cd gsl; ./autogen.sh; ./configure; make; make libgsl.la; cd ..

clean:
	rm -f $(OBJS)

# DO NOT DELETE THIS LINE -- make depend depends on it.

main.o: parse_bam.h
parse_bam.o: randl_list.h parse_bam.h exp_fit.h
exp_fit.o: exp_fit.h
