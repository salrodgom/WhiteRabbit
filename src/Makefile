LINKFLAGS_FOR=-O2 -march=native
COMP_FOR = gfortran
OBJS =  vectors.o	\
	main.o		\
	cif2pdb.o	\
	histograms.o	\
        output.o	\
	dft.o
all:
	${COMP_FOR} ${LINKFLAGS_FOR} -c vectors.f90
	${COMP_FOR} ${LINKFLAGS_FOR} -c cell.f90
	${COMP_FOR} ${LINKFLAGS_FOR} -c output.f90
	${COMP_FOR} ${LINKFLAGS_FOR} -c main.f90
	${COMP_FOR} ${LINKFLAGS_FOR} -o White_Rabbit main.o vectors.o cell.o output.o
	${COMP_FOR} ${LINKFLAGS_FOR} cif2pdb.f90 -o Cheshire_Cat
	${COMP_FOR} ${LINKFLAGS_FOR} histograms.f90 -o Queen_of_Hearts
	${COMP_FOR} ${LINKFLAGS_FOR} dft.f90 -o The_Hatter
	g++ allocate.c funciones_aqui.c funciones_par.c flux.cpp -lm -lgsl -lgslcblas  -o Alicia
clean:;         @rm -f *.mod *.o White_Rabbit Cheshire_Cat Queen_of_Hearts The_Hatter Alicia
