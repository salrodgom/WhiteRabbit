#ifndef _FUNCIONES_H
#define _FUNCIONES_H
//cabeceras de las funciones

#define F77NAME(a) a ## _

void cosa (void);
//void red_random(nodo * M,int N,double mu,gsl_rng* r,int weighted,int autoloops);
void p_random(double *p,int NSTATES,gsl_rng *r);
void red_regular(nodo *M,int N,int K,gsl_rng *r);
void leer_red(int *N,nodo **M,char *nombre);
void print_red(nodo *M,int N);

void loops(int NODO,nodo *M,int cont,int signo,int final,int * rastro,int *estado,int contmax);
void loops_record(int NODO,nodo *M,edge *Links,int cont,int signo,int final,int * rastro,int * estado,int contmax);
void quitar_nodo(int NODO,nodo * M,int N);
void quitar_nodo_sim(int NODO,nodo * M,int N);
void estrada_index(int N,int **A);
void loops_estrada(int N,int **A);

double power_law(int N,float gamma,gsl_rng *r);

// TOPOLOGIAS----------------------------------------------
void simetrizar(nodo *M,nodo *Msim,int N);
void randomizar(nodo *M,nodo *Mrand ,int N,gsl_rng *r);
void randomize_over(nodo *M,int N,gsl_rng *r,int times);
void ordenar_red(nodo *my_node,int N,int *my_source,int Nsources);
void rewire_directional(nodo *my_node,int N,gsl_rng *r,int *OK_label);
//---------------------------------------------------------
//DAR TAMAÑO-----------------------------------------------
//void crear_matriz(int N,nodo *M,int **A);
//void get_edges(nodo *M,int N,edge **my_edge);
//int * get_sources(nodo *M,int N,gsl_rng *r,int *Nsources);
//int * get_drains(nodo *M,int N,gsl_rng *r,int *Ndrains);
//HIERARCHY------------------------------------------------
void calc_hierarchy(nodo *my_node,int N,int Nbb);
void lu_solve(gsl_matrix *m_a, gsl_vector *v_b, gsl_vector *v_x);
void error_S(const char* reason,const char* file,int line,int gsl_errno);
//---------------------------------------------------------
//HISTOGRAMS-----------------------------------------------
void add_to_histogram(double x, double A, double DELTA_X, int *histogram);
void make_histogram(double *x,int *histo_x,int N,double A,double B,double DX);
void print_histogram(int nbins, double A, double DELTA_X, int events, int *histogram);
void print_histogram_double(int nbins, double A, double DELTA_X, int events, double *histogram);
//---------------------------------------------------------
//--------------------------------------------------------
void mutar(nodo *M2,nodo *M1,int N,gsl_rng *r);
void get_network(int n, nodo *M,double *w);
unsigned long int fact(int n,int stop);

void erase_important_links(nodo *M,int N,edge **my_edge,int *how_many);
void erase_determined_links(nodo *M,int N,edge **my_edge,int who);
void how_many_to_remove(edge *my_edge,long double Nloops_rand, long double Nloops_exp,int contmax,int *who,int *how_many,int Nedges);

extern "C"{
void F77NAME(dnaupd)(int *ido, char *bmat, int *n, char *which,
		     int *nev, double *tol, double *resid,
		     int *ncv, double *V, int *ldv,
		     int *iparam, int *ipntr, double *workd,
		     double *workl, int *lworkl, int *info);	
void F77NAME(dneupd)(int *rvec, char *HowMny, int *select,
		     double *dr,double *di, double *Z, int *ldz,
		     double *sigmar,double *sigmai,double *workev,char *bmat, int *n,
		     char *which, int *nev, double *tol,
		     double *resid, int *ncv, double *V,
		     int *ldv, int *iparam, int *ipntr,
		     double *workd, double *workl,
		     int *lworkl, int *info);
}

// ----------GRAFICOS--------------------------//
//----------------------------------------------//
void colorea(float s,int nodo,char *color);
void define_tipo(int *sources,int Nsources,int* drains,int Ndrains,int nodo,char *color);
void imp_grafo_levels(char *nombre,int N, nodo * my_node,int mode,gsl_rng *r);
void imp_grafo_components(char *nombre,int N, nodo * my_node,int mode,gsl_rng *r,int *sources,int *drains,int Nsources, int Ndrains);


#endif
