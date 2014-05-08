#ifndef _EIGEN_H
#define _EIGEN_H
//EIGENVALUES----------------------------------------------
double get_Lmax(nodo *M,int N,double gamma,double *w);
void espectral(int N,nodo *M,double *lambdamax,double gamma);
void espectral_arpack(int n,nodo *M,double *lambdamax,double gamma,double *w);
void gr_eigenvalues(nodo *M, int *nev, int rvec,double *eval, double *evec,int ncv, int maxitr,int n,double gamma,double *w);
void adj_multipication(nodo *M,int n,double *x, double *y,double gamma);
#endif
