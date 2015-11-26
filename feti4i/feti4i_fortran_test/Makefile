FC := mpif90
CC := mpicc
MKDIR := mkdir

build: test

feti4i_mod.o: ../feti4i_mod.F90
	${FC} -c -fPIC $^

feti4i.mod: 	../feti4i_mod.F90 feti4i_mod.o
	@true

feti4i_dummy.o: feti4i_dummy.c
	${CC} -c -fPIC -g $^

lib/libfeti4i.so: feti4i_dummy.o
	${MKDIR} lib
	${LINK.F} -o lib/libfeti4i.so feti4i_dummy.o -shared -fPIC

test.o: test.F90 feti4i.mod 
	${FC} -c -fPIC -g test.F90 -o $@ 

test: test.o lib/libfeti4i.so 
	${LINK.F} -o test test.o -Llib -Wl,-rpath=lib -lfeti4i
  
clean:
	${RM} test.o feti4i_dummy.o feti4i_mod.o feti4i.mod
	${RM} -r lib

.PHONY: build clean

