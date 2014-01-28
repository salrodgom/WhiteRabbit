#LINKFLAGS_FOR=-O -g -march=native
LINKFLAGS_FOR=-O3
COMP_FOR = gfortran
OBJS =  vectors.o	\
	main.o		\
	cif2pdb.o	\
	histograms.o	\
        output.o	\
	dft.o
all:
	${COMP_FOR} ${LINKFLAGS_FOR} -c vectors.f95
	${COMP_FOR} ${LINKFLAGS_FOR} -c cell.f95
	${COMP_FOR} ${LINKFLAGS_FOR} -c output.f95
	${COMP_FOR} ${LINKFLAGS_FOR} -c main.f95
	${COMP_FOR} ${LINKFLAGS_FOR} -o White_Rabbit main.o vectors.o cell.o output.o
	${COMP_FOR} ${LINKFLAGS_FOR} cif2pdb.f95 -o Cheshire_Cat
	${COMP_FOR} ${LINKFLAGS_FOR} histograms.f95 -o Queen_of_Hearts
	${COMP_FOR} ${LINKFLAGS_FOR} dft.f95 -o The_Hatter
	g++ allocate.c funciones_aqui.c funciones_par.c flux.cpp -lm -lgsl -lgslcblas  -o Alicia
clean:;         @rm -f $(OBJS)
