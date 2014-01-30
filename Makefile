#LINKFLAGS_FOR=-O -g -march=native -fbounds-check
F90  =  gfortran
LINKFLAGS_FOR= -O3
COMP_FOR  = ${F90} ${LINKFLAGS_FOR}
COMP_FORc = ${COMP_FOR} -c
OBJS =  vectors.o	\
	cell.o		\
	main.o		\
        output.o
COBJS = ${OBJS}		\
	cif2pdb.o 	\
	histograms.o	\
	dft.o
histograms:
	${F90} ${LINKFLAGS_FOR} histograms.f95 -o Queen_of_Hearts
loops:
	g++ allocate.c funciones_aqui.c funciones_par.c flux.cpp -lm -lgsl -lgslcblas -o Alicia
loops_dbg:
	g++ allocate.c funciones_aqui.c funciones_par.c flux.cpp -lm -lgsl -lgslcblas -fbounds-check -o Alicia
all:
	${COMP_FORc} vectors.f95
	${COMP_FORc} cell.f95
	${COMP_FORc} output.f95
	${COMP_FORc} main.f95
	${F90} ${OBJS} -o White_Rabbit
	${COMP_FOR} cif2pdb.f95 -o Cheshire_Cat
	${COMP_FOR} histograms.f95 -o Queen_of_Hearts
	${COMP_FOR} dft.f95 -o The_Hatter
	g++ allocate.c funciones_aqui.c funciones_par.c flux.cpp -lm -lgsl -lgslcblas  -o Alicia
clean:; rm -f $(COBJS) *.mod ; rm -f Alicia Queen_of_Hearts Cheshire_Cat White_Rabbit The_Hatter
