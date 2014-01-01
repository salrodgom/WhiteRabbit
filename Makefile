LINKFLAGS_FOR=-O3 -g -march=native
COMP_FOR = gfortran
all:
	${COMP_FOR} ${LINKFLAGS_FOR} -c vectors.f95
	${COMP_FOR} ${LINKFLAGS_FOR} -c main.f95
	${COMP_FOR} ${LINKFLAGS_FOR} -o White_Rabbit main.o vectors.o
	${COMP_FOR} ${LINKFLAGS_FOR} cif2pdb.f95 -o Cheshire_Cat
	${COMP_FOR} ${LINKFLAGS_FOR} histograms.f95 -o Queen_of_Hearts
	${COMP_FOR} ${LINKFLAGS_FOR} dft.f95 -o The_Hatter
	g++ allocate.c funciones_aqui.c funciones_par.c flux.cpp -lm -lgsl -lgslcblas  -o Alicia
