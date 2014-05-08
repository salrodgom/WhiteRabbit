#include "general.h"

//----------------------------------------------------------------------------------------
// Calculo de Av con gsl, calcula todos los AV
void espectral(int N,nodo *M,double *lambdamax,double gamma){
float **A;
A=new float *[N];
printf("espectral\n");
for (int i=0;i<N;i++) {
	A[i]=new float [N];
	for (int j=0;j<N;j++){
		A[i][j]=0;
	}
}
for (int i=0;i<N;i++){
	for(int j=0;j<M[i].kin;j++){
		A[M[i].in_nodos[j]][i]=M[i].in_sig[j];// asigna un peso random entre 0 y 2, por ahora sin signo
		//printf("%f\n",M[i].in_sig[j]);
		//A[i][M[i].in_nodos[j]]=1;
	}
/*	for (int j=0;j<M[i].kout;j++){*/
/*		A[M[i].out_nodos[j]][i]=1;// asigna un peso random entre 0 y 2, por ahora sin signo*/
/*		//A[M[i].out_nodos[j]][i]=1;*/
/*	}*/
}
for (int i=0;i<N;i++) A[i][i] = A[i][i]-gamma;
/*printf("espectral 1 matriz\n");*/
/*for(int i=0;i<N;i++){*/
/*	for (int j=0;j<N;j++){*/
/*		printf("%f ",A[i][j]);*/
/*	}*/
/*	printf("\n");*/
/*}*/
/*printf("\n");*/

//A[1][4]=3.;
//A[4][3]=2.;
//A[3][1]=2.;

//Ahora vamos a diagonalizar esta matriz semipositiva
double *evalR,*evalC;
evalR=new double [N];
evalC=new double [N];
//guardamos sitio para la matriz gsl
gsl_matrix * m = gsl_matrix_calloc (N, N);
//copiamos la matriz a diagonalizar
for (int i = 0; i < N; i++) for (int j = 0; j < N; j++)
	 gsl_matrix_set (m, i, j, 1.*A[i][j] );
printf("matriz gsl:\n");
/*for (int i=0;i<N;i++) {*/
/*	for (int j=0;j<N;j++) {*/
/*		printf("#%2f ",gsl_matrix_get(m,i,j));*/
/*	}*/
/*	printf("#\n");*/
/*}*/
/*printf("#\n");*/
 
//CALCULO :
gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (N);
gsl_vector_complex *eval = gsl_vector_complex_alloc (N);
gsl_matrix_complex *evec = gsl_matrix_complex_alloc (N, N); 
gsl_eigen_nonsymmv(m,eval,evec,w);
gsl_eigen_nonsymmv_free (w);

//gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC); //ESTO ORDENA LOS AV	 
	
//printf("#AV:\n");
for (int i=0;i<N;i++) {
	evalR[i]=GSL_REAL(gsl_vector_complex_get(eval,i));
	evalC[i]=GSL_IMAG(gsl_vector_complex_get(eval,i));
	//printf("av[%d]=%f+%fi; ",i,evalR[i],evalC[i]);
}
//printf("# \n");
*lambdamax=-100.; 
for (int i=0;i<N;i++) {
	if (evalR[i] > *lambdamax) {
		*lambdamax=evalR[i];
	}
}

//printf("#lambdamax gsl=%f\n",*lambdamax);

/*for (int i=0;i<N;i++) {*/
/*	for (int j=0;j<N;j++) {*/
/*		printf("#%.2f ",gsl_matrix_get(m,i,j));*/
/*	}*/
/*	printf("#\n");*/
/*}*/
/*printf("#\n");*/

/*printf("#Avector:\n");*/
/*for (int i=0;i<1;i++){*/
/*	printf("#v[%d]: ",i);*/
/*	for (int j=0;j<N;j++){*/
/*		printf("#%.3f, ",GSL_REAL(gsl_matrix_complex_get(evec,i,j)));*/
/*	}*/
/*	printf("#\n");*/
/*}*/

gsl_matrix_free(m);
gsl_vector_complex_free(eval);
gsl_matrix_complex_free(evec);
delete []evalR;
delete []evalC;
for (int i=0;i<N;i++) delete []A[i];
delete []A;

//printf("#lambdamax=%.3f\n",*lambdamax);
return;
}
//------------------------------------------------------------------------//
//--------------- AUTOVALORES ARPACK --------------------------------------------
double get_Lmax(nodo *M,int N,double gamma,double *w) //devuelve el autovalor con mayor parte real -> inestabilidad
{
	double Lambda;
	
	espectral_arpack(N,M,&Lambda,gamma,w);
	
	//printf("#LAmbda arpack return=%f\n",Lambda);
	
	//espectral(N,M,&Lambda,gamma);
	
	//printf("LAmbda gsl return=%f\n",Lambda);
	
	return Lambda;
}//-------------------------------------------------------------------
void espectral_arpack(int n,nodo *M,double *lambdamax,double gamma,double *w) //llamada a subrutina arpack par calcular los n autovalores
{

	int nev=1; //el numero de AV a sacar
	int rvec=0; //vamos a sacar AVectores
	double *eval, *evec;
	eval = (double*)malloc(nev*sizeof(double));//Solo quiero el mayor
	evec = (double*)malloc(n*(nev+1)*sizeof(double));//tb salen solo un autovetor
	
	gr_eigenvalues(M,&nev,rvec,eval,evec,n,2*n,n,gamma,w); //aquí está tratando todavía con el grafo,M

	*lambdamax = eval[0];
	//printf("#Lmax=%f\n",*lambdamax);
	free(eval);
	free(evec);
	
	return;

}
//--------------------------------------------------------------------
void gr_eigenvalues(nodo *M, int *nev, int rvec,double *eval, double *evec,int ncv, int maxitr,int n,double gamma, double *w) 
{	
	int ido = 0; //TODO
  	int ldv = n; //TODO
 	int iparam[11];//TODO
  	int ipntr[14];//TODO aqui pone 11 pero creo q tendria q tener 14
  	int lworkl;//TODO tamaño del workspace
 	int info = 0;//TODO 
  	char bmat= 'I'; //TODO
  	char which[10]="LR"; //TODO ademas tiene un switch para poderlo cambiar 
  	double tol = 1.0e-4;//TODO
  	double *resid;//TODO
  	resid=(double *) malloc (n*sizeof(double));
  	double *v;//TODO
  	//double workd[3*(n)];//TODO
  	double *workd;
  	workd=(double *) malloc (3*n*sizeof(double));
  	double *workl;//TODO
  	int *select; //TODO este es para el post processing
  	int ierr;//TODO este lo de

	double sigma;//post processing
	
//printf("#arpack matrix:\n");
/*for(int i=0;i<n;i++){*/
/*	for (int j=0;j<n;j++){*/
/*	   if (j != i) {*/
/*		double presente=0.;*/
/*		for (int k=0;k<M[i].kout;k++){*/
/*			if (j== M[i].out_nodos[k]) presente=M[i].out_sig[k];*/
/*		}*/
/*		if (presente!= 0) printf("#%f ",presente);*/
/*		else printf("#%f ",0.);*/
/*	   }*/
/*	   else {*/
/*		double presente=0.;*/
/*		for (int k=0;k<M[i].kout;k++){*/
/*			if (j== M[i].out_nodos[k]) presente=M[i].out_sig[k];*/
/*		}*/
/*		if (presente!= 0) printf("#%f ",presente - gamma);*/
/*		else printf("#%f ",0. - gamma);*/
/*	   }*/
/*	*/
/*	}*/
/*	printf("\n");*/
/*	*/
/*}*/

  //---------ESTOS SON LOS PARAMETROS PARA EL ARPACK----------------------------//
  	iparam[0] = 1;   /* ishfts */ //TODO
 	iparam[2] = maxitr; /* maxitr */ //TODO maxitr=N
  	iparam[6] = 1;   /* mode   */ //TODO
  	//which = "LR"; //es el que necesito por ahora
	//          'LM' -> want the NEV eigenvalues of largest magnitude.
	//          'SM' -> want the NEV eigenvalues of smallest magnitude.
	//          'LR' -> want the NEV eigenvalues of largest real part.
	//          'SR' -> want the NEV eigenvalues of smallest real part.
	//          'LI' -> want the NEV eigenvalues of largest imaginary part.
	//          'SI' -> want the NEV eigenvalues of smallest imaginary part.
  	if (ncv > n) ncv = n; //Este valor tiene q estar acotado asi
  	if (ncv <= *nev) ncv = (*nev)+1;
  	
  	//lworkl = ncv*(ncv+8); //Esto si el problema es simetrico
  	lworkl=3*(ncv*ncv)+6*ncv+1;
  	workl = (double*)malloc(lworkl*sizeof(double));//TODO le damos tamaño a workl con lworkl
  	v = (double*)malloc(n*ncv*sizeof(double));
  	int vez=0;
  	do {
		//TODO para usar este cambiar which a "LA"  	
/*  		F77NAME(dsaupd)(&ido, &bmat, &n, which, nev, &tol, resid, &ncv,*/   //Para problemas simetricos
/*		    v, &ldv, iparam, ipntr, workd, workl,*/
/*		    &lworkl, &info);*/

		F77NAME(dnaupd)(&ido, &bmat, &n, which, nev, &tol, resid, &ncv,
		    v, &ldv, iparam, ipntr, workd, workl,
		    &lworkl, &info);
		    
		if ((ido == -1) || (ido == 1)) adj_multipication(M,n,workd+ipntr[0]-1, workd+ipntr[1]-1,gamma); //esta es la multiplicacion q hay q meter
		
		else break;
	
	//printf("vez %d ido=%d\n",vez, ido);
	vez += 1;
		
  	} while (1);
  	
  	if (info < 0) {
   		printf("dsaupd error %d: %d\n", info, iparam[4]);
   		get_network(n,M,w);
    		free(workl);
    		free(v);
    		exit (1);
  	}
  	
  	*nev = iparam[4];  /* number of converged eigenvalues */
  	 select = (int*)malloc(ncv*sizeof(int));
  	 
/*  	 F77NAME(dseupd)(&rvec, "All", select, eval, v, &ldv, &sigma, */ //para el problema simetrico
/*		  &bmat, &n, which, nev, &tol, resid, &ncv, v, &ldv, */
/*		  iparam, ipntr, workd, workl, &lworkl, &ierr);*/

	//definimos las variables necesarias para el nosimetrico
	int ldz;
   	double sigmar,sigmai;//da un poco igual
   	double *DR,*DI; //aqui metemos la parte real e imaginaria de los AV
  	double *z,*workev; //espacios de trabajo y output
  
 	//DR = (double*)malloc((*nev+1)*sizeof(double));
  	DI = (double*)malloc((*nev+1)*sizeof(double));
  	z = (double*)malloc(n*(*nev+1)*sizeof(double));
  	workev = (double*)malloc(3*ncv*sizeof(double));

	F77NAME(dneupd)(&rvec, "A", select, eval, DI, z, &ldz, &sigmar,&sigmai,
	    workev, &bmat, &n, which, nev, &tol, resid, &ncv, v, &ldv, 
            iparam, ipntr, workd, workl, &lworkl, &ierr);
  	
  	free(workl);
  	free(select);
  	free(workev);
  	
  	if (ierr != 0) { 
    		printf("dseupd error %d\n", ierr);
    		free(v);
    		exit (1);
  	}
  	
  	
  	if (rvec) {
   		memcpy(evec, v, (*nev)*n*sizeof(double));
   	}
  	
  	//for(int i=0;i<*nev;i++) printf("av[%d]=%f\n",i,eval[i]);
  	
  	
  	free(v);
  	free(z);
  	free(DI);
  	
  	free(workd);
  	free(resid);
	return;
}//-------------------------------------------------------------------------
void adj_multipication(nodo *M,int n,double *x, double *y,double gamma)
{
	for(int i=0;i<n;i++){
	y[i]=0.;
		for (int j=0;j<M[i].kout;j++){
			y[i] += x[M[i].out_nodos[j]]*M[i].out_sig[j];
		}
	y[i] -= x[i]*gamma;
	
	}


	return;
}
//------------------------------------------------------------------------

 
