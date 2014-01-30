#LINKFLAGS_FOR=-O -g -march=native -fbounds-check
F90  = gfortran
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
	${F90} ${LINKFLAGS_FOR} histograms.f90 -o Queen_of_Hearts
loops:
	g++ allocate.c funciones_aqui.c funciones_par.c flux.cpp -lm -lgsl -lgslcblas -o Alicia
loops_dbg:
	g++ allocate.c funciones_aqui.c funciones_par.c flux.cpp -lm -lgsl -lgslcblas -fbounds-check -o Alicia
all:
	${COMP_FORc} vectors.f90
	${COMP_FORc} cell.f90
	${COMP_FORc} output.f90
	${COMP_FORc} main.f90
	${F90} ${OBJS} -o White_Rabbit
	${COMP_FOR} cif2pdb.f90 -o Cheshire_Cat
	${COMP_FOR} histograms.f90 -o Queen_of_Hearts
	${COMP_FOR} dft.f90 -o The_Hatter
	g++ allocate.c funciones_aqui.c funciones_par.c flux.cpp -lm -lgsl -lgslcblas  -o Alicia
clean:; rm -f $(COBJS) *.mod ; rm -f Alicia Queen_of_Hearts Cheshire_Cat White_Rabbit The_Hatter
