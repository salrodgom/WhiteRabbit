#include "general.h"
#include <gsl/gsl_permutation.h>
//estas variables hace falta porque se meten como externas en la hoja de funciones
int **ciclos;
int cont_ciclos;
double *L,*Lp,*Ln;
double *L_estr;
int error;
//---------------------------------------------------------------------------//

int main (int argc, char *argv[]){

//// ------ GENERADOR DE NUMEROS ALEATORIOS ------//
const gsl_rng_type * T;
gsl_rng * r;
gsl_rng_env_setup();
if (!getenv("GSL_RNG_SEED")) gsl_rng_default_seed = time(0);
T = gsl_rng_default;
r = gsl_rng_alloc(T);
////----------------------------------------------//

//propiedades de las redes--------------
int N,Nlinks,contmax;
double mu,p,gamma,gamma_o,gamma_rand,lambda_rand,lambda_o;
double delta_0,delta,delta2,sigma_delta;
double gamma2_rand,sigma_gammarand;
double lambda; //Parametro de direccion
long double Nloops_exp,Nloops_rand;
char nombre[30];
int number;
nodo *M;
nodo *Mrand;
nodo *Msim;


int *sources,*drains;
int Nsources=0; int Ndrains=0;

//Calculo paralelo -----------------	
//int iCPU = omp_get_num_procs();
//printf("#el n procesadores es %d\n",iCPU);
int iCPU=1;
int num_proc = iCPU;

//////////////////////////////////////////////////////////////////////
int modo=2; //1 comparacion con teorico , 2 busqueda de alfa en reales, 3 - Busqueda de alfa con redes con estructura
//////////////////////////////////////////////////////////////////////

// LEO LA ENTRADA
if(argc < 2 ){
	        printf("ERROR: debes dar el numero de nodos y el nombre\n");
		        return 1;
}
	N=atoi(argv[1]);
	//mu=atof(argv[2]);
	p=atof(argv[2]);
	mu=N*p;
	 //introduzco la conectividad
	//gamma=atof(argv[3]);//introdusco gamma
	contmax=atoi(argv[4]); //la longitud maxima de los loops a buscar
	strcpy(nombre,argv[3]);//copiamos el nombre
	int simple=atoi(argv[5]);
	
	if (modo!=1) M=leer_red_new(&N,nombre); //necesaario para que el tamaño de los vectores L sea el correcto
	if (modo==2) liberar_nodos(M,N); //necesario para poder vovler a rellenarla luego
	//print_red(M,N);
	
	//return 0;
	
	//delta=get_double_links_prob(M,N);
	//printf("gamma+1-gamma=%lf\n",1.-delta);
	//return 0;
	//Ya tengo N
	L=new double [N+1]; Lp=new double [N+1]; Ln=new double [N+1];
	
	
	//number=get_net_number(nombre);
	//printf("N=%d mu=%f gamma=%f Lmax=%d nombre=%s\n",N,mu,gamma,contmax,nombre);
//Cosas de Loops-------------------
double *LDr_med,*LD2r_med,*sigmar_LDmed; //dirigidos
double *LNDr_med,*LND2r_med,*sigmar_LNDmed; //no dirigidos
double *Loopinesr_med,*Loopines2r_med,*sigmar_Loop; //fraccion
double *Loopines_t;
double *LD_t,*LND_t,*Loop_t;
double *Ldir,*Lsim;//dirigidos y no dirigidos
int cont_ciclos; //longitud maxima de los loops


LDr_med=new double [N+1];LD2r_med=new double[N+1];sigmar_LDmed=new double [N+1];//randomizaciones


LNDr_med=new double [N+1];LND2r_med=new double[N+1];sigmar_LNDmed=new double [N+1];//randomizaciones
Loopinesr_med=new double [N+1];Loopines2r_med=new double [N+1];sigmar_Loop=new double [N+1];
Loopines_t=new double [N+1];
for (int i=0;i<N+1;i++) {LNDr_med[i]=LND2r_med[i]=0.; LDr_med[i]=LD2r_med[i]=0.; Loopinesr_med[i]=0.;}
double *Lreal,*Lrand;
Lreal=new double [N+1]; Lrand=new double [N+1];
LD_t=new double [N+1];	LND_t=new double [N+1]; Loop_t=new double [N+1];
	
//Histogramas de jerarquia-----------------
int *Histo_nivel_exp;
int *Histo_nivel_rand;
int nbins=N;//TODO
Histo_nivel_exp=(int *)calloc(nbins,sizeof(int));
Histo_nivel_rand=(int *)calloc(nbins,sizeof(int));
double a=0.; double b=N*1.; double dx=1.;

//Estadistica -----------------------------
double *Histo_nivel_rand_med,*Histo_nivel_rand2_med,*sigma_rand;
double Nsources_med,Ndrains_med,gamma_med,gammar_med;
double Nsources2_med,Ndrains2_med,gamma2_med,gamma2r_med;
double sigma_Nsources,sigma_Ndrains,sigma_gamma,sigmar_gamma;
double Nsourcesr_med,Ndrainsr_med,Nsources2r_med,Ndrains2r_med,sigmar_Nsources,sigmar_Ndrains;
Histo_nivel_rand_med=(double *)calloc(nbins,sizeof(double));
Histo_nivel_rand2_med=(double *)calloc(nbins,sizeof(double));
sigma_rand=(double *)calloc(nbins,sizeof(double));
Nsources_med=Nsources2_med=Ndrains=Ndrains2_med=gamma_med=gamma2_med=gammar_med=gamma2r_med=Nsourcesr_med=Ndrainsr_med=Nsources2r_med=Ndrainsr_med=0.0;

//Salida de informacion------------------
FILE *vfich;


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////             PROGRAMA PRINSIPAL            ////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	//printf("calculo teorico: Esto solo saca los elementos de la amtriz A\n");
//contmax=5;	
int **A;
//A=(int**)malloc(contmax*sizeof(int*));
//for (int i=0;i<contmax;i++) {
//	A[i]=(int*)calloc(contmax,sizeof(int));
//}
goto aqui;
A=new int *[contmax+1];
for (int i=0;i<contmax+1;i++) A[i]=new int [contmax+1]; 
for (int i=0;i<contmax+1;i++) for (int j=0;j<contmax+1;j++) A[i][j]=0;

int l;

for (int j=3;j<contmax+1;j++){
	printf("i=%d\n",j);
       gsl_permutation * perm = gsl_permutation_alloc (j);
     
       gsl_permutation_init (perm);
     
       do 
        {
        	//gsl_permutation_fprintf (stdout, perm, " %u");
        	//printf ("\n");
        
        	l=0;//ponemos el contador a cero
           	for (int i=0;i<j;i++){//ahora vemos cuantos a favor hay
           		if ( gsl_permutation_get(perm,(i+1)%j) > gsl_permutation_get(perm,(i%j)) ) l ++;
           		 
           	}
           	//printf("l=%d\n",l);
           	A[l][j] ++;
           //printf ("\n");
        }
       while (gsl_permutation_next(perm) == GSL_SUCCESS);
     
       gsl_permutation_free (perm);
          
	//for (int i=0;i<j+1;i++) printf("MAtriz[%d][%d]=%d\n",i,j,A[i][j])     ;
}
     printf("ya tengo los elementos de A\n");
       
        //calculo teorico:
        //lambda=0.2;
//       for (int i=3;i<contmax;i++){
//        //printf("longitud k=%d:\n",i);
//        Loopines_t[i]=0;
//        unsigned long int norma;
//        norma=fact(i,1);
//        //printf("%d!=%d %lf Loopines[%d,%.1lf]: \n",i,norma,1./(norma*1.),i,lambda);
//        	for (int l=1;l<i;l++){
//        		//printf("+ A[%d][%d]/%d!(pow(%.1lf,%d)*pow(%.1lf,%d)+pow(%.1lf,%d)*pow(%.1lf,%d))\n %d/%d [%lf*%lf + %lf*%lf] =",l,i,i,lambda,l,1.-lambda,i-l,lambda,i-l,1.-lambda,l,A[l][i],norma,pow(lambda,l),pow(1.-lambda,i-l),pow(lambda,i-l),pow(1.-lambda,l));
//        		//printf("%lf*[%lf*%lf + %lf*%lf]=",(A[l][i]*1.)/(norma*1.),pow(lambda,l),pow(1.-lambda,i-l),pow(lambda,i-l),pow(1.-lambda,l));
//        		Loopines_t[i] += (pow(lambda,l)*pow(1.-lambda,i-l) + pow(lambda,i-l)*pow(1.-lambda,l))*((A[l][i]*1.)/(norma*1.));
//        		//printf("%lf \n",Loopines_t[i]);
//        	} 
//        	//printf("%lf \n",Loopines_t[i]);       
//        }  
//        
//        for (int i=3;i<contmax;i++) printf("L[%d]=%lf\n",i,Loopines_t[i]);
//          
//        //return 0;
        
        
        
        //calculo prediccion teorica dirigida, solo necesario para las random del modo2
//	printf("calculo teorico dirigido\n");
//	//p=mu/(N*1.);
//	printf("N=%d p=%f\n",N,p);
//	for (int i=3;i<N+1;i++) LD_t[i]=0.;
//	LD_t[2]=pow(p,2)*(exp(lgamma(N+1)-lgamma(N-2+1))); //tenemos que multiplicar por 2, ya que si no no cuenta bien los de 2 nodos
//	printf("L[%d]=%f\n",2,LD_t[2]);
//	for (int i=3;i<contmax+1;i++){
//		LD_t[i]=pow(p,i)*(exp(lgamma(N+1)-lgamma(N-i+1)))/i; //con la función gamma ya no explotan!
//		printf("LD_t[%d]=%f\n",i,LD_t[i]);
//	}
//	
//	printf("calculo teórico no dirigido\n");
//	printf("p=%f\n",p);
//	////printf("calculo teórico no dirigido\n");
//	////printf("p=%f\n",p);
//	if ((2*p-p*p)<1.) p=(2*p-p*p);
//	else p=1.;
//	for (int i=3;i<N+1;i++) LND_t[i]=0.;
//	LND_t[2]=pow(p,2)*(exp(lgamma(N+1)-lgamma(N-2+1)))/2; //aqui no dividimos por 2*k porque no lo estabamos considerando en los dirigidos
//	printf("L[%d]=%f\n",2,LND_t[2]);
//	for (int i=3;i<contmax+1;i++){
//		LND_t[i]=pow(p,i)*(exp(lgamma(N+1)-lgamma(N-i+1)))/(2*i); //con la función gamma ya no explotan!
//		printf("LD_t[%d]=%f\n",i,LND_t[i]);
//	}
//	
//	for (int i=3;i<contmax;i++) {
//		Loop_t[i]=max(0.,LD_t[i]/LND_t[i]);
//		printf("Loopines_t[%d]=%lf\n",i,Loop_t[i]);
//	}
//	
//	 //CALCULO TEORICO
//	 lambda=0.5;
//        for (int i=3;i<contmax;i++){
//        //printf("longitud k=%d:\n",i);
//        Loopines_t[i]=0;
//        unsigned long int norma;
//        norma=fact(i,1);
//        //printf("%d!=%d %lf Loopines[%d,%.1lf]: \n",i,norma,1./(norma*1.),i,lambda);
//        	for (int l=1;l<i;l++){
//        		//printf("+ A[%d][%d]/%d!(pow(%.1lf,%d)*pow(%.1lf,%d)+pow(%.1lf,%d)*pow(%.1lf,%d))\n %d/%d [%lf*%lf + %lf*%lf] =",l,i,i,lambda,l,1.-lambda,i-l,lambda,i-l,1.-lambda,l,A[l][i],norma,pow(lambda,l),pow(1.-lambda,i-l),pow(lambda,i-l),pow(1.-lambda,l));
//        		//printf("%lf*[%lf*%lf + %lf*%lf]=",(A[l][i]*1.)/(norma*1.),pow(lambda,l),pow(1.-lambda,i-l),pow(lambda,i-l),pow(1.-lambda,l));
//        		Loopines_t[i] += (pow(lambda,l)*pow(1.-lambda,i-l) + pow(lambda,i-l)*pow(1.-lambda,l))*((A[l][i]*1.)/(norma*1.));
//        		//printf("%lf \n",Loopines_t[i]);
//        	} 
//        	//printf("%lf \n",Loopines_t[i]);       
//        }  
//        
//        for (int i=3;i<contmax;i++) printf("Loopines_0.5[%d]=%lf\n",i,Loopines_t[i]);
//          
//	//return 0;
//	
//	///////////////////////////////////////////////////////////////////////////////////////
        
        
	//mu=2.*mu; //esto es aqui porque la hago no dirigida apra luego dirigirla!!!! TODO ojo con esto
	
	
if (modo==1){
	
	lambda=0.;
	double dlambda=0.1;
	
printf("comenzamos el bucle\n");
//------------------------------------------------//
FILE *vfich5;
        char arch_name5[60];
        sprintf(arch_name5,"gammaVSlambda.dat");
        vfich5=fopen(arch_name5,"w");
//-----------------------------------------------//
while(lambda <= 0.01){
//lambda =0.5;
gamma_rand=0.; gamma2_rand=0.; sigma_gammarand=0.; 
printf("lambda=%.2lf\n",lambda);
char arch_name1[60];
sprintf(arch_name1,"loops_%.2lf.dat",lambda);
   for (int i=0;i<N+1;i++) {LNDr_med[i]=LND2r_med[i]=0.; LDr_med[i]=LD2r_med[i]=0.;} //no qeuremos q se mezclen q cada vez es con un gamma
      
      int realizaciones=100;
      for (int vez=0;vez<realizaciones;vez++){
      printf("vez %d\n",vez);
	
	//printf("red random con N=%d mu=%lf lambda=%lf \n",N,mu,lambda);
	Mrand=red_random_new(N,mu,r,0,0); //es no pesada y dirigida por defecto	
	//print_red(Mrand,N);
	//return 0;
	//printf("la dirigimos\n");
	Mrand=dirigir_red(Mrand,N,lambda,r);//borra la matriz entrante
	printf("-----------------------------\n");
	//print_red(Mrand,N);

	//Mrand=red_random_new(N,mu,r,1,0);
	
	// SACAMOS LOS LOOPS DIRIGIDOS
	//---------------------------------------------------------------------------------------
	printf("sacando loops dirigidos\n");
				if (contmax > N) contmax=N;	 			
	 			for (int i=0;i<N+1;i++) L[i]=Lp[i]=Ln[i]=Lreal[i]=0;
	 			
				#pragma omp parallel for schedule(dynamic,1) num_threads(num_proc) 
     				for (int i=0;i<N;i++){ 
       				//printf("Desde %d\n",i);
				//int this_thread = omp_get_thread_num(), num_threads = omp_get_num_threads();
				////printf("Hola, soy la hebra %d de un total de %d\n",this_thread,num_threads);
				int cont=0;
				int signo=1;//el signo de un loop particular
				int *rastro;
				int *estado;
				estado=new int[N];
				for(int j=0;j<N;j++) estado[j]=1;
	
				////printf("encontramos loops desde %d \n",i);
				rastro=(int *)malloc(sizeof(int)*1000);	

	  			 loops(i,Mrand,cont,signo,i,rastro,estado,contmax);
	
				delete []estado;
				free(rastro);

    				 }//End quitar nodos --------------------------------------fin cilco principal---- a partir de aqui no hay mas M!! ahora si,pqno la borro
				//return 0;
				//return 0;
				L[2]=L[2]*2; //estos solo aparecen en una direccion por culpa del algoritmo de busqueda, asi q multiplicamos por dos
				for (int i=0;i<contmax+1;i++) {
					LDr_med[i]+=L[i];
					LD2r_med[i]+=L[i]*L[i];
				}
				//for (int i=3;i<contmax+1;i++) printf("LD[%d]=%lf  \n",i,L[i]);	
	//---------------------------------------------------------------------------------------//
	
	Msim=simetrizar(Mrand,N);//TODO volver a poner si queremos los loops
	// SACAMOS LOOPS NO DIRIGIDOS --------------------------------------------------------
	printf("sacamos loops no dirigidos\n");
				for (int i=0;i<N+1;i++) L[i]=Lp[i]=Ln[i]=Lreal[i]=0;
	 			
				#pragma omp parallel for schedule(dynamic,1) num_threads(num_proc) 
     				for (int i=0;i<N;i++){ 
       				//printf("Desde %d\n",i);
				//int this_thread = omp_get_thread_num(), num_threads = omp_get_num_threads();
				////printf("Hola, soy la hebra %d de un total de %d\n",this_thread,num_threads);
				int cont=0;
				int signo=1;//el signo de un loop particular
				int *rastro;
				int *estado;
				estado=new int[N];
				for(int j=0;j<N;j++) estado[j]=1;
	
				////printf("encontramos loops desde %d \n",i);
				rastro=(int *)malloc(sizeof(int)*1000);	

	  			loops(i,Msim,cont,signo,i,rastro,estado,contmax);
	
				delete []estado;
				free(rastro);

    				 }//End quitar nodos --------------------------------------fin cilco principal---- a partir de aqui no hay mas M!! ahora si,pqno la borro
				//return 0;
				//return 0;
				L[2]=L[2]*2; //estos solo aparecen en una direccion por culpa del algoritmo de busqueda, asi q multiplicamos por dos
				for (int i=0;i<contmax+1;i++) L[i]=L[i]/2.;
				for (int i=0;i<contmax+1;i++) {
					LNDr_med[i]+=L[i];
					LND2r_med[i]+= L[i]*L[i];
				}
				//for (int i=3;i<contmax+1;i++) printf("LND[%d]=%lf  \n",i,L[i]);	
				
	//--------------------------------------------------------------------------------------------------------
	
	for (int i=0;i<contmax;i++){
		Loopinesr_med[i]+=max(0.,LDr_med[i]/LNDr_med[i]);
		
	}
	double gamma;
	double *vector;
	gamma=get_Lmax(Mrand,N,0.,vector);
	gamma_rand += gamma ;
	gamma2_rand += gamma*gamma;
	printf("Lmax=%lf\n",gamma);
	liberar_nodos(Mrand,N);
	liberar_nodos(Msim,N);
      }//END estdistica--------------------
      
      for (int i=0;i<contmax;i++){
      	LDr_med[i]=LDr_med[i]/(realizaciones*1.);
        LD2r_med[i]=LD2r_med[i]/(realizaciones*1.);
        //no dirigidos RANDOMIZADOS
        LNDr_med[i]=LNDr_med[i]/(realizaciones*1.);
        LND2r_med[i]=LND2r_med[i]/(realizaciones*1.);
        sigmar_LDmed[i]=sqrt(fabs(LD2r_med[i]-(LDr_med[i]*LDr_med[i])));
        sigmar_LNDmed[i]=sqrt(fabs(LND2r_med[i]-(LNDr_med[i]*LNDr_med[i])));
        Loopinesr_med[i]=Loopinesr_med[i]/(realizaciones*1.);
        //propagacion de errores
        double a,b;
        a=(sigmar_LDmed[i]*sigmar_LDmed[i])*max(0.,1./(LNDr_med[i]*LNDr_med[i]));
        b=(sigmar_LNDmed[i]*sigmar_LNDmed[i])*(LDr_med[i]*LDr_med[i])*max(0.,1./(LNDr_med[i]*LNDr_med[i]*LNDr_med[i]*LNDr_med[i]));        
        sigmar_Loop[i]=sqrt(fabs(a+b));
      }
      gamma_rand=gamma_rand/(realizaciones*1.);
      gamma2_rand=gamma2_rand/(realizaciones*1.);
      sigma_gammarand=sqrt(fabs(gamma2_rand-gamma_rand*gamma_rand));
      printf("con lambda %lf Lmax=%lf +- %lf \n",lambda,gamma_rand,sigma_gammarand);
        
        //for (int i=3;i<contmax;i++) printf("Loopines_r[%d]=%lf\n",i,Loopinesr_med[i]);
        
        //for (int i=3;i<contmax;i++) printf("Loopines_t[%d]=%lf\n",i,Loop_t[i]);
        
        //CALCULO TEORICO
        for (int i=3;i<contmax+1;i++){
        //printf("longitud k=%d:\n",i);
        Loopines_t[i]=0;
        unsigned long int norma;
        norma=fact(i,1);
        //printf("%d!=%d %lf Loopines[%d,%.1lf]: \n",i,norma,1./(norma*1.),i,lambda);
        	for (int l=1;l<i;l++){
        		//printf("+ A[%d][%d]/%d!(pow(%.1lf,%d)*pow(%.1lf,%d)+pow(%.1lf,%d)*pow(%.1lf,%d))\n %d/%d [%lf*%lf + %lf*%lf] =",l,i,i,lambda,l,1.-lambda,i-l,lambda,i-l,1.-lambda,l,A[l][i],norma,pow(lambda,l),pow(1.-lambda,i-l),pow(lambda,i-l),pow(1.-lambda,l));
        		//printf("%lf*[%lf*%lf + %lf*%lf]=",(A[l][i]*1.)/(norma*1.),pow(lambda,l),pow(1.-lambda,i-l),pow(lambda,i-l),pow(1.-lambda,l));
        		//Loopines_t[i] += (pow(lambda,l)*pow(1.-lambda,i-l) + pow(lambda,i-l)*pow(1.-lambda,l))*((A[l][i]*1.)/(norma*1.));
        		//printf("%lf \n",Loopines_t[i]);
        	} 
        	//printf("%lf \n",Loopines_t[i]);       
        }  
        
        //for (int i=3;i<contmax+1;i++) printf("Loopines_0.5[%d]=%lf\n",i,Loopines_t[i]);
          
        //return 0;
        
        //Sacar la realizacion con este valor del parametro: Saco la fraccion de loops de cada longitud para un alpha determinado        
        //vfich=fopen(arch_name1,"w");
        //fprintf(vfich,"LAMBDA 'l=%.2lf'\n",lambda);
        //for (int i=3;i<contmax+1;i++){
        //	fprintf(vfich,"%d %lf %lf %lf %lf \n",i,max(0.,LDr_med[i]/LNDr_med[i]),sigmar_Loop[i],lambda,Loopines_t[i]);
        //}
        //fclose(vfich);
        
        
        
        fprintf(vfich5,"%lf %lf %lf\n",lambda,gamma_rand,sigma_gammarand);
        
        
        lambda += dlambda;    
      
}//END busqueda de parametros WHILE lambda <1 

fclose(vfich5);

}//END modo ==1
aqui:
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////7
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 if (modo ==2){
 
 //printf("entro en el modo2\n");
 
 
  //CALCULO TEORICO
  int K;
  lambda=0.; double dlambda=0.001;
  K=(int)(1./dlambda);
//  printf("dlambda =%lf K=%d y la mitad de k es %d\n",dlambda,K,(int)(K/2.));
//  
//  //genero una matriz con los valores teoricos almacenados
 double **Loopines_tm;
// Loopines_tm=(double **)malloc((contmax+1)*sizeof(double*));
// for (int i=0;i<contmax+1;i++){
// 	Loopines_tm[i]=(double*)calloc(K,sizeof(double));
// }  
// 
//	if (simple==1) {
//		printf("simple\n");
//		delta=0.;
//	}
//  K=0;
//  while (lambda<1.){  
//  
//        for (int i=3;i<contmax+1;i++){
//        //printf("longitud k=%d lambda=%lf:\n",i,K*dlambda);
//        //printf("longitud k=%d:\n",i);
//        Loopines_tm[i][K]=0;
//        unsigned long int norma;
//        norma=fact(i,1);
//        //printf("%d!=%d %lf Loopines[%d,%.1lf]: \n",i,norma,1./(norma*1.),i,lambda);
//        	for (int l=1;l<i;l++){
//        		//printf("+ A[%d][%d]/%d!(pow(%.1lf,%d)*pow(%.1lf,%d)+pow(%.1lf,%d)*pow(%.1lf,%d))\n %d/%d [%lf*%lf + %lf*%lf] =",l,i,i,lambda,l,1.-lambda,i-l,lambda,i-l,1.-lambda,l,A[l][i],norma,pow(lambda,l),pow(1.-lambda,i-l),pow(lambda,i-l),pow(1.-lambda,l));
//        		//printf("%lf*[%lf*%lf + %lf*%lf]=",(A[l][i]*1.)/(norma*1.),pow(lambda,l),pow(1.-lambda,i-l),pow(lambda,i-l),pow(1.-lambda,l));
//        		//Loopines_tm[i][K] += (pow(lambda,l)*pow(1.-lambda,i-l) + pow(lambda,i-l)*pow(1.-lambda,l))*((A[l][i]*1.)/(norma*1.));
//        		//Ahora con el parametro de loops dobles:
//        		Loopines_tm[i][K] += (pow(lambda*(1.-delta)+delta,l)*pow((1.-lambda)*(1.-delta)+delta,i-l) + pow(lambda*(1.-delta)+delta,i-l)*pow((1.-lambda*(1.-delta)+delta),l))*((A[l][i]*1.)/(norma*1.));
//        		//Loopines_t[i] += (pow(lambda,l)*pow(1.-lambda,i-l) + pow(lambda,i-l)*pow(1.-lambda,l))*((A[l][i]*1.)/(norma*1.));
//        		//printf("%lf \n",Loopines_t[i]);
//        		//printf("K=%d Loopines_t[%d][%f]=%f\n",K,i,K*dlambda,Loopines_tm[i][K]);
//        	} 
//        	//printf("%lf \n",Loopines_t[i]);       
//        }  
//  
//  K+=1;
//  lambda += dlambda;
//  
//  }//END calculo teorico lambda parameter , debemos tener una amtriz con los valores teoricos guardados
// 
// for (int i=3;i<contmax+1;i++){
// printf("Loopines_t[%d]=",i);
// 	for (int j=0;j<K;j++){
// 		printf("[%.2lf] %.2lf, ",j*dlambda,Loopines_tm[i][j]);
// 	}
// 	printf("\n");
// }
 
 //return 0;
 
 M=leer_red_new(&N,nombre);
 //tenemos que sacar la loopines para cada tamaño
 
 // Saca la loopines de las redes y las guarda en loopines_nombre.dat
 	FILE *vfich3;
 	char nombre_arch3[70];
 	sprintf(nombre_arch3,"loopines_%s_exp.dat",nombre);
 	//printf("archivo: %s\n",nombre_arch3);
 
 for (int i=0;i<contmax+1;i++) LDr_med[i]=LNDr_med[i]=0.;
 // SACAMOS LOS LOOPS DIRIGIDOS
	//---------------------------------------------------------------------------------------
	//printf("sacando loops dirigidos\n");
				if (contmax > N) contmax=N;	 			
	 			for (int i=0;i<N+1;i++) L[i]=Lp[i]=Ln[i]=0;
	 			
				#pragma omp parallel for schedule(dynamic,1) num_threads(num_proc) 
     				for (int i=0;i<N;i++){ 
       				//printf("Desde %d\n",i);
				//int this_thread = omp_get_thread_num(), num_threads = omp_get_num_threads();
				////printf("Hola, soy la hebra %d de un total de %d\n",this_thread,num_threads);
				int cont=0;
				int signo=1;//el signo de un loop particular
				int *rastro;
				int *estado;
				estado=new int[N];
				for(int j=0;j<N;j++) estado[j]=1;
	
				////printf("encontramos loops desde %d \n",i);
				rastro=(int *)malloc(sizeof(int)*1000);	

	  			//loops(i,M,cont,signo,i,rastro,estado,contmax);
	
				delete []estado;
				free(rastro);

    				 }//End quitar nodos --------------------------------------fin cilco principal---- a partir de aqui no hay mas M!! ahora si,pqno la borro
				//return 0;
				//return 0;
				L[2]=L[2]*2; //estos solo aparecen en una direccion por culpa del algoritmo de busqueda, asi q multiplicamos por dos
				for (int i=0;i<contmax+1;i++) {
					LDr_med[i] += L[i];
					LD2r_med[i] += L[i]*L[i];
				}
				//for (int i=3;i<contmax+1;i++) printf("LD[%d]=%lf  \n",i,L[i]);
				//for (int i=3;i<contmax+1;i++) printf("LDmed[%d]=%lf  \n",i,LDr_med[i]);	
	//---------------------------------------------------------------------------------------//
	//printf("simetrizar\n");
	Msim=simetrizar(M,N);
	//Msim=leer_red_new(&N,nombre);
	//printf("termino de sim\n");
	int todo_cero=1;
	for (int i=3;i<contmax+1;i++) if (L[i]!=0) todo_cero=0; 
	
	if (todo_cero==1) {
		//printf("todo cero, voy al final\n");
		//goto saltar_no_dirigidos;
	}
	
	// SACAMOS LOOPS NO DIRIGIDOS --------------------------------------------------------
	//printf("sacamos loops no dirigidos\n");
				for (int i=0;i<(N+1);i++) L[i]=Lp[i]=Ln[i]=0;
				//#pragma omp parallel for schedule(dynamic,1) num_threads(num_proc) 
     				for (int i=0;i<N;i++){ 
     				
       				//printf("Desde %d\n",i);
				//int this_thread = omp_get_thread_num(), num_threads = omp_get_num_threads();
				////printf("Hola, soy la hebra %d de un total de %d\n",this_thread,num_threads);
				int cont=0;
				int signo=1;//el signo de un loop particular
				int *rastro;
				int *estado;
				estado=new int[N];
				for(int j=0;j<N;j++) estado[j]=1;
	
				////printf("encontramos loops desde %d \n",i);
				rastro=(int *)malloc(sizeof(int)*1000);	

	  			loops(i,Msim,cont,signo,i,rastro,estado,contmax);
	
				delete []estado;
				free(rastro);

    				 }//End quitar nodos --------------------------------------fin cilco principal---- a partir de aqui no hay mas M!! ahora si,pqno la borro
				//return 0;
				//return 0;
				L[2]=L[2]*2; //estos solo aparecen en una direccion por culpa del algoritmo de busqueda, asi q multiplicamos por dos
				for (int i=0;i<contmax+1;i++) L[i]=L[i]/2.;
				for (int i=0;i<contmax+1;i++) {
					LNDr_med[i]+=L[i];
					LND2r_med[i]+= L[i]*L[i];
				}
				//for (int i=3;i<contmax+1;i++) printf("LND[%d]=%lf  \n",i,L[i]);	
				//for (int i=3;i<contmax+1;i++) printf("LD[%d]=%lf  \n",i,LNDr_med[i]);	
				
	//--------------------------------------------------------------------------------------------------------
	
saltar_no_dirigidos:
	
	vfich3=fopen(nombre_arch3,"w");
	fprintf(vfich3,"#Loopines #LD #LND\n");
	for (int i=3;i<contmax+1;i++){
		//printf("Loopines[%d]=%lf/%lf=%lf\n",i,LDr_med[i],LNDr_med[i],Loopinesr_med[i]);
		Loopinesr_med[i]+=max(0.,LDr_med[i]/LNDr_med[i]);
		//printf("Loopines[%d]=%lf/%lf=%lf\n",i,LDr_med[i],LNDr_med[i],Loopinesr_med[i]);
		fprintf(vfich3,"%d %e %e %e\n",i,Loopinesr_med[i],LDr_med[i],LNDr_med[i]);
	}
	
	
	
	
	return 0;
	//Ahora ya tengo el valor de loopines para cada tamaño, tengo q buscar el gamma q mas se acerca
	double gamma_exp;
	double close,closer,dif;  
	

//	 for (int i=3;i<contmax+1;i++){
//	 printf("longitud=%d: Lreal=%lf=\n",i,Loopinesr_med[i]);
//	  	for (int j=0;j<K;j++){
//	 		printf("[%.2lf]=%.2lf, ",j*dlambda,Loopines_tm[i][j]);
//	 	}
//	 	printf("\n");
//	 }
	//return 0;
	//--------------------------------------------------//
	// Saca los gamma fitteados para cada red y los va añadiendo en un archivo
	FILE *vfich2;
	vfich2=fopen("gamma_parameter_basura_boltzman.dat","a");
//	 
	for (int i=3;i<contmax+1;i++){ //para cada longitud
	dif=0.;close=1.; closer=0.;
	printf("Longitud %d Lreal=%lf=\n",i,Loopinesr_med[i]);
		for (int j=0;j<(int)(K/2.)+1;j++){
				dif=fabs(Loopinesr_med[i]-Loopines_tm[i][j]);
				//printf("dif[%lf]=%lf (%lf - %lf)\n",j*dlambda,dif,Loopinesr_med[i],Loopines_tm[i][j]);
				//printf("Loop_t[%d][%f]=%f- L[%d]=%f =%lf \n",i,dlambda*j,Loopines_tm[i][j],Loopinesr_med[i],dif);
				if (dif < close) {
				printf("dif[%lf]=%lf (%lf - %lf)\n",j*dlambda,dif,Loopinesr_med[i],Loopines_tm[i][j]);
//				//printf("Loop_t[%d][%f]=%f- L[%d]=%f =%lf \n",i,dlambda*j,Loopines_tm[i][j],Loopinesr_med[i],dif);
				printf("este es menor j=%d gamma=%lf\n",j,j*dlambda);
					close=dif;
					closer=j*dlambda;
				}
		}
		gamma_exp += closer;
		printf("para longitud %d el mas cercano es gamma=%lf\n",i,closer);
		fprintf(vfich2,"%d %lf %d #%s\n",i,closer,number,nombre);	
	}
	gamma_exp=gamma_exp/(contmax+1);
	printf("el gamma_exp mas cercano es %lf\n",gamma_exp);
	//Tenemos que ver a cual de los que tenemos se parece mas
	gamma_exp=((int)(gamma_exp/dlambda))*(dlambda);
//	printf("el gamma_exp mas cercano es %lf\n",gamma_exp);
	fprintf(vfich2,"\n\n");
// 
 	fclose(vfich2);
 	//--------------------------------------------------//
 	//--------------------------------------------------//
 	
 	vfich3=fopen(nombre_arch3,"w");
 	
 	double *sigma_loop;
 	sigma_loop=(double*)malloc((contmax+1)*sizeof(double));
 	for (int i=3;i<contmax+1;i++) sigma_loop[i]=sqrt(fabs(Loopinesr_med[i]/LNDr_med[i]*1.));
 	
 	for (int i=3;i<contmax+1;i++){
 		fprintf(vfich3,"%d %lf %lf %lf %lf\n",i,Loopinesr_med[i],LDr_med[i],LNDr_med[i],sigma_loop[i]);
 		printf("%d %lf %lf %lf %lf\n",i,Loopinesr_med[i],LDr_med[i],LNDr_med[i],sigma_loop[i]);
 	}
 	fclose(vfich3);
 	//--------------------------------------------------//
 	//--------------------------------------------------//
// 	FILE *vfich4;
// 	char nombre_arch4[60];
// 	sprintf(nombre_arch4,"gamma_%s.dat",nombre);
// 	vfich4=fopen(nombre_arch4,"w");
// 	int index1,index2,index3;
// 	for (int i=0;i<contmax+1;i++) {
// 		
// 	}
 	
 for (int i=0;i<contmax+1;i++) delete []A[i];
 delete []A;	
 
 }//END MODO 2




return 0;
}
