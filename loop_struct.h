#ifndef _LOOP_STRUCT_H
#define _LOOP_STRUCT_H

struct loop_stat { //Esto genera dos vectores de punteros
  
  int nbins;
  int contmax;
  double A; //el inicio del histograma
  
  int *redes; //numero de redes en cada bin
  double **Lp;//almacena la suma de loops en cada bin
  double **L2p;//almacena el numero de loops al cuadrado
  double **sigma_Lp;
  
  
} ;

#endif

