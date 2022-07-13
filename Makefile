.PHONY: all clean rebuild run docs

CFLAGS ?= -Ofast -march=native -Wall -fopenmp -lm -lgomp
CC     ?= gcc

OBJS = random.o v3d.o rk4.o #SI_units.o

all: v1

run: v1.o v1
	./v1

clean:
	rm -rf ./*.o

rebuild: clean all

docs:
	doxygen docs.doxy    

#SI_units.o: SI_units.c
#	$(CC) SI_units.c $(CFLAGS) -c -o $@ 
 
random.o: random.c
	$(CC) random.c $(CFLAGS) -c -o $@ 
 
v3d.o: v3d.c
	$(CC) v3d.c $(CFLAGS) -c -o $@

rk4.o: rk4.c
	$(CC) rk4.c $(CFLAGS) -c -o $@

v1.o: 
	$(CC) v1.c $(CFLAGS) -c -o $@ 

v1: $(OBJS) v1.o 
	$(CC) -o $@ $(OBJS) v1.o $(CFLAGS) 
