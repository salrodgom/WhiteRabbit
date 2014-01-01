#ifndef _ALLOCATE_H
#define _ALLOCATE_H

nodo * leer_red_new(int *N,char *nombre);
nodo * red_random_new(int N,double mu,gsl_rng *r,int directed,int weighted);
nodo * red_config(int N,double mu,double sigma,gsl_rng *r,int directed,int weighted,char *type);
nodo * red_config2(int N,double mu,double sigma,gsl_rng *r,int directed);

nodo * simetrizar(nodo *M,int N);
nodo * simetrizar_2(nodo *M,int N);
nodo * randomizar(nodo *M,int N,gsl_rng *r,int times);
nodo * copiar_nodos(nodo * my_original,int N);
nodo * leer_matriz(int *Na,int *Np,char *nombre); //TODO lee "matrices"" no dirigidas!!!
int ** crear_matriz_new(int N,nodo * my_node);
void p_random(double *p,int NSTATES,gsl_rng *r); //genera una lsita de numeros que suman 1

edge * get_edges(nodo *my_node,int N,int *Nlinks,int contmax);
int * get_sources(nodo *M,int N,gsl_rng *r,int *Nsources);
int * get_drains(nodo *M,int N,gsl_rng *r,int *Ndrains);

loop_stat set_loop_stat(int nbins,int contmax);

void liberar_nodos(nodo *my_node,int N);
void liberar_edges(edge *my_edge);
void liberar_int(int *vector);
void liberar_loop_stat(loop_stat my_L);

#endif
