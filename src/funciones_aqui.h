#ifndef _FUNCIONES_AQUI_H
#define _FUNCIONES_AQUI_H

void get_hierarchy(nodo *my_node,int N,int Nsources,int Ndrains);
nodo * crop_network(nodo *my_node,int N,edge *my_link,int Nlinks,long double Nloops_rand,long double Nloops_exp,int contmax);
void resta(long double Nloopsfin,long double Nloops,int cont,edge * my_edge,size_t *imp_index,nodo * my_node,int N,int Nedges,int contmax);
void recalcular_links(edge *my_edge,nodo *my_node,int N,int nodo_in,int nodo_out,int Nedges,int contmax);
void loops_crop(int NODO,nodo *M,edge *Links,int cont,int final,int *rastro,int *estado,int contmax,int Nedges);//saca los loops
void reordenar_importance(edge *my_edge,size_t *imp_index,int Nedges);//devuelve el vector de indices de importancia cambiado
void quitar_link(nodo *my_node,int nodo_in,int nodo_out);//elimina un link de una red, solo lo deja inaccesible, pero no cambia la dimension
nodo * dirigir_red(nodo *my_node,int N,double lambda,gsl_rng *r);
int get_net_number(char *nombre);
double get_double_links_prob(nodo *my_node,int N);

double calc_coherence(nodo *my_node,int N,int Nsources);

#endif
