#LINKFLAGS_FOR=-O -g -march=native
FOR = gfortran
LINKFLAGS_FOR=-O3
F90  = ${FOR} ${LINKFLAGS_FOR} -c
OBJS =  vectors.o	\
	main.o		\
        output.o	\
	cell.o
#
all:
	${F90} -c vectors.f95
	${F90} -c cell.f95
	${F90} -c output.f95
	${F90} -c main.f95
	${FOR} ${OBJS} -o White_Rabbit
	${F90} cif2pdb.f95 -o Cheshire_Cat
	${F90} histograms.f95 -o Queen_of_Hearts
	${F90} dft.f95 -o The_Hatter
	g++ allocate.c funciones_aqui.c funciones_par.c flux.cpp -lm -lgsl -lgslcblas  -o Alicia
clean:;         @rm -f $(OBJS) *.mod ; rm -f White_Rabbit Cheshire_Cat Queen_of_Hearts The_Hatter Alicia

