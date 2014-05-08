#include "general.h"
//#include "multicanonical.h"

// Esto lo q hace es remitirme al fichero general donde estan todos los include q necesito, y se comparten tanto por el programa general como por las subrutinas a parte
//ojo, tengo q incluir el general.h en TODOS los archivos .c (funcion) o .cpp (main) para que pueda acceder a las librerias
// En los .h tengo q tener cuidado de poner el #if inicial siempre, para que no pete si los llama mas de una vez

extern int **ciclos;
extern int cont_ciclos;
extern double *L;
extern double *Lp;
extern double *Ln;
extern double *L_estr;
extern int error;

//------------ RED REGULAR ------------------------------------------------------------------------//
void red_regular(nodo *M,int N,int K,gsl_rng *r){
int nodo2,*posibles;
int *cont,nodov;
for (int i=0;i<N;i++) {M[i].kout=0; M[i].kin=K; M[i].label=1; M[i].myself=i;}//pongo todas las conectividades out a cero
for(int i=0;i<N;i++) {
	M[i].in_nodos=new int [K];
	M[i].in_sig=new double [K];
	for(int j=0;j<M[i].kin;j++){
		p_random(M[i].in_sig,M[i].kin,r);//rellenamos M[i].in_sig con kin numeros q sumen 1		
	}
}
printf("K=%d\n",K);

posibles=new int[N];

for (int i=0;i<N;i++){

int cont=0;//contado de los  lleva pillados
for(int j=0;j<N;j++) posibles[j]=j;
//printf("nodo[%d]:",i);

	while (cont < K) {
		do {
			nodo2=gsl_rng_uniform_int(r,N);
			
		} while ((nodo2 == i) or (posibles[nodo2]==-1));
		//printf("nodo2=%d\n",nodo2);
		posibles[nodo2]=-1;
		M[i].in_nodos[cont]=nodo2;
		cont+=1;
		M[nodo2].kout+=1;
	}
}

cont=new int [N]; for (int i=0;i<N;i++) cont[i]=0; //pongo a cero los contadores de k out de los nodos

for (int i=0;i<N;i++) { //preparamos los vectores de nodos out
	M[i].out_nodos=new int [M[i].kout];
	M[i].out_sig=new double [M[i].kout];
	for(int j=0;j<M[i].kout;j++){
		M[i].out_sig[j]=1.;
	}
}

for (int i=0;i<N;i++) {
	for (int j=0;j<M[i].kin;j++){
		nodov=M[i].in_nodos[j]; //este es el nodo victima
		M[nodov].out_nodos[cont[nodov]]=i; //ponemos al que se lo come en su vector de outs
		cont[nodov]+=1;//sumamos uno al contador del nodo victima
	}
}


/*for (int i=0;i<N;i++){*/
/*printf("Nodo %d , se come a :",i);*/
/*	for (int j=0;j<M[i].kin;j++){*/
/*		printf("%d (%f),",M[i].in_nodos[j],M[i].in_sig[j]);*/
/*		*/
/*	}*/
/*	printf("\n");*/
/*printf("es comido por :");*/
/*for (int j=0;j<M[i].kout;j++){*/
/*		printf("%d ",M[i].out_nodos[j]);*/
/*		*/
/*	}*/
/*	printf("\n");*/
/*}*/


delete []posibles;
delete []cont;

return;
}

//-------------------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------------------//
//            Ordena la red de manera q pone primero a las especies q estan dentro del vector my_source
void ordenar_red(nodo *my_node,int N,int *my_source,int Nsources)
{	
	int cont_sources=0;
	int cont_normal=Nsources;
	int source;
	

	for (int i=0;i<N;i++){
		source=0;//llevo 0 sources puestas de Nsources
		for (int j=0;j<Nsources;j++){
			if ( i == my_source[j] ) source=1; //si es una fuente porque está en el vector
		
		}
		if (source == 1) {
			my_node[i].label = cont_sources; //le pongo un label de cont_sources < Nsources :Lo pongo al principio
			cont_sources ++;
		}
		else {
			my_node[i].label = cont_normal;// si no es una source, empiezo a ponerle labels > Nsources
			cont_normal ++;
		}
	}
	
	//for (int i=0;i<N;i++) printf("[%d]='%d'\n",i,my_node[i].label);

	return;
}
//---------------------------------------------------------------------------------------//
double power_law(int N,float gamma,gsl_rng *r)
{
	double x,u;
	double xMAX=sqrt(N);
	double xMIN=1.;
	double norma,c;
	
	u=gsl_rng_uniform(r);
	
	if (gamma == 1.) {
		norma=log(xMAX/xMIN);
		x=xMIN*exp(u*norma);
		
	}
	
	else {
		norma=(gamma -1.)/(pow(xMIN,1.-gamma)-pow(xMAX,1.-gamma));
		
		x=pow(u * ((1.-gamma)/norma) + pow(xMIN,1.-gamma) , 1./(1.-gamma) );
	
	}
	
	
	return x;
}//-------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------
//--------------------------------------------------------------------------//
void loops(int NODO,nodo *M,int cont,int signo,int final,int * rastro,int * estado,int contmax){//saca los loops

int nodo2;
//printf("nivel %d : %d   ",cont,NODO);

//M[NODO].label=0;
estado[NODO]=0;
rastro[cont]=NODO;
cont=cont+1;

for (int i=0;i<M[NODO].kout;i++){
	nodo2=M[NODO].out_nodos[i];
	//printf("%d ",nodo2);
	if (nodo2 == final) {
		//printf("hemos encontrado uno de tamaño %d: \n",cont);
		#pragma omp atomic
		L[cont]++;
		
		if (signo > 0){ //si el loop es positivo lo metemos en el contador positivo
		#pragma omp atomic
		Lp[cont]++;
		}
		
		else if (signo < 0){ //si el loop es negativo, lo metemos en el contador negativo
		#pragma omp atomic 
		Ln[cont]++;
		}
		
		if (cont==contmax){//-----------IMPRIME EL LOOP------------------
			//ciclos[cont_ciclos]=new int [cont];
			for(int o=0; o<cont; o++){
				//ciclos[cont_ciclos][o]=rastro[o];
				printf("%d ",rastro[o]+1);
			}
//			printf("%d",final+1);
/*			//printf("(%d)",signo);*/
			printf("\n");
			//sleep(1);
			//cont_ciclos ++;
		}//------------------------------------------------------------
		
	}
	else {
		if ((estado[nodo2] == 1) and (nodo2 > final) and (cont < contmax)) {//lo visitamos ;(nodo2 > final) ->eswto es como si eliminasemos los nodos de los q sabemos todos los loops
			//M[nodo2].label=0;
			//printf("level %d :",cont);
			signo=signo*M[NODO].out_sig[i];
			loops(nodo2,M,cont,signo,final,rastro,estado,contmax);
			//printf("%d*%d    ",M[NODO].out_sig[i],nodo2);
		}
		else {	

		}
	}

}
//M[NODO].label=1;
estado[NODO]=1;

return;
}

//--------------------------------------------------------------------------//
void loops_record(int NODO,nodo *M,edge *Links,int cont,int signo,int final,int *rastro,int *estado,int contmax){//saca los loops
int nodo2;
//printf("nivel %d : %d   ",cont,NODO);

//M[NODO].label=0;
estado[NODO]=0;
rastro[cont]=NODO;
cont ++ ;

for (int i=0;i<M[NODO].kout;i++){
	nodo2=M[NODO].out_nodos[i];
	//printf("%d ",nodo2);
	if (nodo2 == final) {
		//printf("hemos encontrado uno de tamaño %d: \n",cont);
		#pragma omp atomic
		L[cont]++;		
		
		if (signo > 0){ //si el loop es positivo lo metemos en el contador positivo
		#pragma omp atomic
		Lp[cont]++;
		}		
		
		else if (signo < 0){ //si el loop es negativo, lo metemos en el contador negativo
		#pragma omp atomic 
		Ln[cont]++;
		}
	
		rastro[cont]=final;
		//cont ++;
		//printf("ahora cont=%d y rastro[%d]=%d rastro[%d]=%d\n",cont,cont,final,cont,rastro[cont-1]);
/*		for(int o=0; o<cont+1; o++){*/
/*			printf("%d ",rastro[o]);*/
/*		}*/
/*		printf("\n");*/
		for (int o=0; o < cont; o++){//para cada elemento del loop
			//localizo la x
			int x=0;
			for (int m=0; m < rastro[o]; m++){
				//printf("+ Kout[%d]=%d, ",m,M[m].kout);
				x += M[m].kout;
			}
			int y=0;//esto es lo q tnego q sumar para llear al q es
			for (int m=0; m<M[rastro[o]].kout; m++){
				if (rastro[o+1] == M[rastro[o]].out_nodos[m]) {
					y=m;
					break;
				}
			}
			
			#pragma omp atomic
			Links[x+y].importance ++; 
			
			//printf("es el elemento %d->%d ,x+y=%d: %d\n",rastro[o],rastro[o+1],x+y,Links[x+y].importance);
		}
		//Ahora tenemos q sumar a los Edges correspondientes uno, porque aparecen en el loop
		//printf("(%d)",signo);
		//sleep(1);
		//cont_ciclos ++;
		
	}
	else {
		if ((estado[nodo2] == 1) and (nodo2 > final) and (cont < contmax)) {//lo visitamos
			//M[nodo2].label=0;
			//printf("level %d :",cont);
			signo=signo*M[NODO].out_sig[i];
			loops_record(nodo2,M,Links,cont,signo,final,rastro,estado,contmax);
			//printf("%d*%d    ",M[NODO].out_sig[i],nodo2);
		}
		else {	
			
		}
	}

}

estado[NODO]=1;

return;
}


void quitar_nodo(int NODO,nodo * M,int N){
//printf("entramos\n");
int nodo1,nodo2,indice;

for (int i=0;i<M[NODO].kin;i++){//busco entre los que se come
	nodo2=M[NODO].in_nodos[i];
	for (int j=0;j<M[nodo2].kout;j++){//busco entre los depredadores del nodo2
		if(M[nodo2].out_nodos[j]==NODO) indice=j;
	}
		//printf("IN en %d presa estaba en %d lugar de %d\n",nodo2,indice,M[nodo2].kout);
		//printf("cambio %d por %d\n",M[nodo2].out_nodos[indice],M[nodo2].out_nodos[ M[nodo2].kout -1]);
		M[nodo2].out_nodos[indice]=M[nodo2].out_nodos[M[nodo2].kout -1];//los intercambio
		M[nodo2].out_nodos[M[nodo2].kout]=NODO;
		M[nodo2].kout --;//Le quitamos el depredador NODO al bich	
}
for (int i=0;i<M[NODO].kout;i++){
	nodo1=M[NODO].out_nodos[i];// busco entre los depredadores de NODO
	//printf("%d ",nodo1);
	for (int j=0;j<M[nodo1].kin;j++){
		if (M[nodo1].in_nodos[j]==NODO) indice=j;
	}
		//printf("OUT en %d depredador estaba en %d lugar\n",nodo1,indice);
		//printf("cambio %d por %d\n",M[nodo1].in_nodos[indice],M[nodo1].in_nodos[ M[nodo1].kin -1]);
		M[nodo1].in_nodos[indice]=M[nodo1].in_nodos[M[nodo1].kin -1];//los intercambio
		M[nodo1].in_nodos[M[nodo1].kin]=NODO;
		M[nodo1].kin --;
}

/*for (int i=NODO+1;i<N;i++){*/
/*printf("Nodo %d , se come a de %d:",i,M[i].kin);*/
/*	for (int j=0;j<M[i].kin;j++){*/
/*		printf("%d ",M[i].in_nodos[j]);*/
/*		*/
/*	}*/
/*	printf("\n");*/
/*printf("es comido por %d :",M[i].kout);*/
/*for (int j=0;j<M[i].kout;j++){*/
/*		printf("%d ",M[i].out_nodos[j]);*/
/*		*/
/*	}*/
/*	printf("\n");*/
/*}*/


return;
}//------------este es el dirigido

void quitar_nodo_sim(int NODO,nodo * M,int N){//ya no se usa, pero va quitando el nodo i cada vez q ha apsado por el en la siemtrica
//printf("entramos\n");
int nodo1,nodo2,indice;

/*for (int i=0;i<M[NODO].kin;i++){//busco entre los que se come*/
/*	nodo2=M[NODO].in_nodos[i];*/
/*	for (int j=0;j<M[nodo2].kout;j++){//busco entre los depredadores del nodo2*/
/*		if(M[nodo2].out_nodos[j]==NODO) indice=j;*/
/*	}*/
/*		//printf("IN en %d presa estaba en %d lugar de %d\n",nodo2,indice,M[nodo2].kout);*/
/*		//printf("cambio %d por %d\n",M[nodo2].out_nodos[indice],M[nodo2].out_nodos[ M[nodo2].kout -1]);*/
/*		M[nodo2].out_nodos[indice]=M[nodo2].out_nodos[M[nodo2].kout -1];//los intercambio*/
/*		M[nodo2].out_nodos[M[nodo2].kout]=NODO;*/
/*		M[nodo2].kout --;//Le quitamos el depredador NODO al bich	*/
/*}*/
for (int i=0;i<M[NODO].kout;i++){
	nodo1=M[NODO].out_nodos[i];// busco entre los depredadores de NODO
	//printf("está en %d ",nodo1);
	for (int j=0;j<M[nodo1].kout;j++){
		if (M[nodo1].out_nodos[j]==NODO) indice=j;
	}
		//printf("en %d lugar de %d \n",indice,M[i].kout);
		//printf("cambio %d por %d\n",M[nodo1].out_nodos[indice],M[nodo1].out_nodos[ M[nodo1].kout -1]);
		M[nodo1].out_nodos[indice]=M[nodo1].out_nodos[M[nodo1].kout -1];//los intercambio
		M[nodo1].out_nodos[M[nodo1].kout]=NODO;
		M[nodo1].kout --;
}

//for (int i=NODO+1;i<N;i++){
//printf("Nodo %d , se come a de %d:",i,M[i].kin);
/*	for (int j=0;j<M[i].kin;j++){*/
/*		printf("%d ",M[i].in_nodos[j]);*/
/*		*/
/*	}*/
/*	printf("\n");*/
/*printf("nodo %d es comido por %d :",i,M[i].kout);*/
/*for (int j=0;j<M[i].kout;j++){*/
/*		printf("%d ",M[i].out_nodos[j]);*/
/*		*/
/*	}*/
/*	printf("\n");*/
/*}*/


return;
}
//----------------------------------------------------------------------------------------------------------------//
void simetrizar(nodo *M,nodo *Msim,int N){ //crea una lista simetrica a aprtir de la lista de nodos no simetrica //TODO OLD
vector < vector <int> > aux_in;
vector < vector <int> > aux_out;
aux_in.resize(N);//"2N vectores para almacenar las conexiones q faltan
aux_out.resize(N);

int repe;
printf("entramos a simemtrizar\n");
for (int i=0;i<N;i++) {Msim[i].kout=0; Msim[i].kin=0; Msim[i].label=1; Msim[i].myself=i;}//pongo las banceritas hacia arriba

//printf("simetrizamos\n");

//printf("copiamos los vectores a los auxiliares\n");
for (int i=0;i<N;i++){
	for (int j=0;j<M[i].kin;j++){
		aux_in[i].push_back(M[i].in_nodos[j]);	
		//printf("M[%d].in_modos[%d]=%d aux_in[%d][%d]=%d\n",i,j,M[i].in_nodos[j],i,j,aux_in[i][j]);	
	}
	for (int j=0;j<M[i].kout;j++){
		aux_out[i].push_back(M[i].out_nodos[j]);
		//printf("M[%d].out_modos[%d]=%d aux_out[%d][%d]=%d\n",i,j,M[i].out_nodos[j],i,j,aux_out[i][j]);
	}
}

/*for (int i=0;i<N;i++){*/
/*printf("Nodo %d , se come a :",i);*/
/*	for (int j=0;j<aux_in[i].size();j++){*/
/*		printf("%d ",aux_in[i][j]);*/
/*		*/
/*	}*/
/*	printf("\n");*/
/*printf("es comido por :");*/
/*for (int j=0;j<aux_out[i].size();j++){*/
/*		printf("%d ",aux_out[i][j]);*/
/*		*/
/*	}*/
/*	printf("\n");*/
/*}*/

//printf("ponemos los q faltan\n");
int nodo2;
for (int i=0;i<N;i++){
	for (int j=0;j<M[i].kin;j++){
		//tenemos que ver si j se come a i
		repe=0;
		nodo2=M[i].in_nodos[j];
		for(int k=0;k<aux_in[nodo2].size();k++){
			if(i==aux_in[nodo2][k]) repe=1;
		}
		if (repe==0) aux_in[nodo2].push_back(i);
	}
	for (int j=0;j<M[i].kout;j++){
		//tenemos que ver si por el q es comido es comido por i ya
		repe=0;
		nodo2=M[i].out_nodos[j];
		for (int k=0;k<aux_out[nodo2].size();k++){
			if (i==aux_out[nodo2][k]) repe=1;
		}
		if (repe==0) aux_out[nodo2].push_back(i);
	}
}

/*for (int i=0;i<N;i++){*/
/*printf("Nodo %d , se come a :",i);*/
/*	for (int j=0;j<aux_in[i].size();j++){*/
/*		printf("%d ",aux_in[i][j]);*/
/*		*/
/*	}*/
/*	printf("\n");*/
/*printf("es comido por :");*/
/*for (int j=0;j<aux_out[i].size();j++){*/
/*		printf("%d ",aux_out[i][j]);*/
/*		*/
/*	}*/
/*	printf("\n");*/
/*}*/

//printf("finalmente tenemos que copiarlos de vuelta\n");
for (int i=0;i<N;i++){
/*Msim[i].kin=aux_in[i].size();*/
/*Msim[i].in_nodos=new int [Msim[i].kin];*/
/*	for (int j=0;j<Msim[i].kin;j++){*/
/*		Msim[i].in_nodos[j]=aux_in[i][j];*/
/*	}*/
Msim[i].kout=aux_out[i].size();
Msim[i].out_nodos=new int [Msim[i].kout];
Msim[i].out_sig=new double [Msim[i].kout];
	for (int j=0;j<Msim[i].kout;j++){
		Msim[i].out_nodos[j]=aux_out[i][j];
	}
}

for (int i=0;i<N;i++){
/*	for (int j=0;j<Msim.kin[i];j++){*/ //solo dejamos esta aprte apra q no salga doble
/*		Msim.in_sig[j]=1;*/
/*	}*/
	for (int j=0;j<Msim[i].kout;j++){
		Msim[i].out_sig[j]=1;
	}
}

//ponemos los signos
for (int i=0;i<N;i++){
	delete [] M[i].in_sig;
	M[i].in_sig=new double [N];
	for (int j=0;j<M[i].kin;j++) M[i].in_sig[j]=1.;
	
	delete [] M[i].out_sig;
	M[i].out_sig=new double [N];
	for (int j=0;j<M[i].kout;j++) M[i].out_sig[j]=1;
}


/*printf("matriz simetrica\n");*/
/*for (int i=0;i<N;i++){*/
/*printf("Nodo %d , se come a %d bischos :",i,Msim[i].kin);*/
/*	for (int j=0;j<Msim[i].kin;j++){*/
/*		printf("%d ",Msim[i].in_nodos[j]);*/
/*		*/
/*	}*/
/*	printf("\n");*/
/*printf("%d es comido por %d bichos:",i,Msim[i].kout);*/
/*for (int j=0;j<Msim[i].kout;j++){*/
/*		printf("%d ",Msim[i].out_nodos[j]);*/
/*		*/
/*	}*/
/*	printf("\n");*/
/*}*/
aux_in.clear();
aux_out.clear();


return;
}
// ----------------------- HISTOGRAMS ---------------------------------------------------------------------------------------------
void add_to_histogram(double x, double A, double DELTA_X, int *histogram)
{
	
	histogram[ (int)( (x-A)/DELTA_X ) ]++;
	return;
}//--------------------------------------------------------------------------------
void make_histogram(double *x,int *histo_x,int N,double A,double B,double DX) //coge un vector de numeros y les hace elm hsitpgrama
{

	for (int i=0;i<N;i++) {//para cada elemento del vector x
		if((x[i]>A) and (x[i]<B)) { //si esté en el rango
			histo_x[(int)((x[i]-A)/(DX))] ++;
		}
	}
	
	return;
} //-------------------------------------------------------------------------------
void print_histogram(int nbins, double A, double DELTA_X, int events,int *histogram)
{
	printf("plot '-' w boxes lw 2\n");

	// print the histogram properly normalized
	
	events=0;
	for(int i=0;i<nbins;i++) events += histogram[i];
	
	//for(int i=0; i<nbins; i++) printf("%f %f\n",A+i*DELTA_X+(DELTA_X/2.),histogram[i]*1./events);
	for(int i=0; i<nbins; i++) printf("%f %f\n",A+i*DELTA_X,histogram[i]*1.);
	
	printf("e\n");
	
	//printf("#pause 1\n");
	
	return;
}//----------------------------------------------------------------------------------
void print_histogram_double(int nbins, double A, double DELTA_X, int events, double *histogram)
{
	printf("plot '-' w boxes lw 2\n");

	// print the histogram properly normalized
	
	//printf("set yr[0:1000]\n");
	
	//for(int i=0; i<nbins; i++) printf("%f %lf\n",A+i*DELTA_X,histogram[i]*1./events);
	for(int i=0; i<nbins; i++) printf("%f %lf\n",A+i*DELTA_X,histogram[i]*1.);
	
	printf("e\n");
	
	//printf("#pause 1\n");
	
	return;
}//----------------------------------------------------------------------------------
//---------------------------- HIERARCHY SECTION ----------------------------------------------------------------//
void calc_hierarchy(nodo *my_node,int N,int Nbb)
{
	//printf("#entramos N=%d Nsource=%d\n",N,Nbb);
	//Generamos la lista de correspondencia
	double *label; label =new double[N]; //el etiquetado label de los nodos (Nsources delante)
	size_t *label_index; label_index=new size_t [N]; //ordeno las especies en funcion de su label: al label 1, que especie original corresponde
	for(int i=0;i<N;i++) label[i]=my_node[i].label; //etiquetado nuevo de lso nodos, puesto con "ordenar"
	
	gsl_sort_index(label_index,label,1,N); //odenamos las etiquetas propias segun el nivel trofico del nodo
	//si originalmente los elementos estaban de 1,...N, ,, al reordenar segun los labels, cómo se queda esa ordenacion inicial?:Se guarda en label_index
	
	//for (int i=0;i<N;i++) printf("v[%d]=%d i[%d]=%d\n",i,(int)label[i],i,label_index[i]);//el nodo "inicial(label)" es el nuevo (label_index)
	
	//printf("#creamos la matriz\n");
	//Generamos la amtriz de la red
	int **A;
	A=new int *[N];
	for(int i=0;i<N;i++) A[i]=new int [N]; 
	for (int i=0;i<N;i++) for (int j=0;j<N;j++) A[i][j]=0;
	for (int i=0;i<N;i++) {
		for (int j=0;j<my_node[i].kin;j++) {
			int x=my_node[i].label;
			int nodoy=my_node[i].in_nodos[j];
			int y=my_node[nodoy].label;
			//printf("elemento [%d][%d] = [%d][%d]\n",i,my_node[i].in_nodos[j],x,y);
			A[x][y]=1;
		}
	}
	
	
/*	for (int i=0;i<min(10,N);i++) {*/
/*		for (int j=0;j<min(10,N);j++) {*/
/*			printf("#%d ",A[i][j]);*/
/*		}*/
/*		printf("#\n");*/
/*	}*/
	
	//rellenamos las matrices a resolver
	gsl_matrix * B = gsl_matrix_calloc (N, N);
	gsl_vector * b = gsl_vector_calloc (N);
	gsl_vector * X = gsl_vector_calloc (N);

	for (int i=0;i<Nbb;i++){ // A todas las fuentes les estoy asociando un valor l=1 by the facer //TODO y borro sus conexiones IN ¿tiene esto efecto al calcular?
		for (int j=0;j<N;j++){
			gsl_matrix_set (B, i, j,0. );
		}
		gsl_matrix_set (B, i, i, 1.);
		gsl_vector_set (b,i,1.);
	}
	//printf("copiamoslos animales\n");
	for (int i=Nbb;i<N;i++){
		for (int j=0;j<N;j++){
			gsl_matrix_set (B, i, j, -(A[i][j]*1.)/(my_node[label_index[i]].kin*1.));
		}
		gsl_matrix_set (B, i, i, 1.-(A[i][i]*1.)/(my_node[label_index[i]].kin*1.));
		gsl_vector_set (b,i,1.);
	}
	
	//printf("resolvemos	\n");
	gsl_error_handler_t * handler_viejuno = gsl_set_error_handler (&error_S);
	lu_solve(B,b,X);
	gsl_set_error_handler(handler_viejuno);

	for(int i=0;i<N;i++){
	        my_node[i].s = gsl_vector_get(X,i);
		//printf("%f ,",gsl_vector_get(X,i));
	}
//	printf("\n");

	//sacamos la Smedia
	float Smed=0.;
	for (int i=0;i<N;i++){
		Smed=Smed+my_node[i].s;
	}
	Smed=Smed/(N*1.);
	
	for (int i=0;i<N;i++) my_node[i].label=label_index[i];
	
	for (int i=0;i<N;i++) {
		delete [] A[i];
	}
	delete [] A;
	delete [] label;
	delete [] label_index;
	gsl_vector_free(b);
	gsl_vector_free(X);
	gsl_matrix_free(B);
	
	return;
}


//------------------------------------------//
//RESOLUCION SIST ECUACUIONES
void lu_solve(gsl_matrix *m_a, gsl_vector *v_b, gsl_vector *v_x){
/***********************************************/
/* T.Kouya's GSL sample program collection     */
/*          Solving Linear System of Equations */
/*                            LU decomposition */
/*                   Written by Tomonori Kouya */
/*                                             */
/* Version 0.01: 2007-10-02                    */
/***********************************************/
// INPUT: matriz y vector gsl con A y B
// OUTPUT: almacena el resultado en un vector gsl X

	int signum;
	gsl_permutation *perm;
	
	/* LU decomposition and forward&backward substition */
	perm = gsl_permutation_alloc(m_a->tda);
	gsl_linalg_LU_decomp(m_a, perm, &signum);
	gsl_linalg_LU_solve(m_a, perm, v_b, v_x);
	gsl_permutation_free(perm);

	
	
	return;
}
void error_S(const char* reason,const char* file,int line,int gsl_errno)
{

	error=1;

	return;
}
//-------------- END HIERARCHY--------------------------------------------//
//------------------------------------------------------------------------//

void loops_estrada(int N,int **A){//saca los loops a la estrada
double **Aaux;
double **A2;

Aaux=new double *[N];
A2=new double *[N];

for (int i=0;i<N;i++) {
	Aaux[i]=new double [N];
	A2[i]=new double [N];
	for (int j=0;j<N;j++) {
		//Aaux[i][j]=0;
		A2[i][j]=0.;
	}
}

for (int i=0;i<N;i++) for (int j=0;j<N;j++) Aaux[i][j]=A[i][j];

for(int k=2;k<N+1;k++){ //para cada longitud
        //calculo la matriz A2
	for (int i=0;i<N;i++){
		for (int j=0;j<N;j++){
		A2[i][j]=0;
			for (int m=0;m<N;m++){
				A2[i][j]=A2[i][j]+A[i][m]*Aaux[m][j];
			}
		}
	}
        //---ya tengo A2
        //actualizo la matriz auxiliar con la A2
        for (int i=0;i<N;i++) for (int j=0;j<N;j++) Aaux[i][j]=A2[i][j];
/*	printf("A^%d\n",k);*/
/*	for (int i=0;i<N;i++){*/
/*		for (int j=0;j<N;j++) {*/
/*			printf("%.1f ",A2[i][j]);*/
/*		}*/
/*		printf("\n");*/
/*	}*/
	for (int i=0;i<N;i++) L_estr[k]=L_estr[k]+A2[i][i];
	L_estr[k]=L_estr[k]/(exp(lgamma(k+1)));
	//printf("L_estr[%d]=%f\n",k,L_estr[k]);



}


return;
}

void estrada_index(int N,int **A){ //saca el estrada index de una red
double *evalR,*evalC;
float EE;
evalR=new double [N];
evalC=new double [N];
//guardamos sitio para la matriz gsl
gsl_matrix * m = gsl_matrix_alloc (N, N);

//copiamos la matriz a diagonalizar
for (int i = 0; i < N; i++) for (int j = 0; j < N; j++)
	 gsl_matrix_set (m, i, j, 1.*A[i][j] );
	 
//CALCULO SIN AUTOVECTORES:
gsl_eigen_nonsymm_workspace * w = gsl_eigen_nonsymm_alloc (N);
gsl_vector_complex *eval = gsl_vector_complex_alloc (N);
gsl_eigen_nonsymm(m,eval,w);
gsl_eigen_nonsymm_free (w);  
//gsl_eigen_symm_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC); //ESTO ORDENA LOS AV


for (int i=0;i<N;i++) {
	evalR[i]=GSL_REAL(gsl_vector_complex_get(eval,i));
	evalC[i]=GSL_IMAG(gsl_vector_complex_get(eval,i));
	//printf("av[%d]=%f+%fi\n",i,evalR[i],evalC[i]);
}
printf(" \n");

gsl_vector_complex_free(eval);

EE=0.;
for (int k=0;k<N;k++) EE=EE+exp(evalR[k]);
//printf("EE=%f\n",EE);

return;
}

void randomizar(nodo *M,nodo *Mrand ,int N,gsl_rng *r){
int RELINKS=100*N;
int nodo1,nodo2;
int nodo1p,nodo2p;
int cambios=0;
//copio todo a Mrand
for (int i=0;i<N;i++){
	Mrand[i].kout=M[i].kout;
	Mrand[i].kin=M[i].kin;
	Mrand[i].label=1;
	Mrand[i].myself=i;
	
	Mrand[i].in_nodos= new int [Mrand[i].kin];
	for (int j=0;j<Mrand[i].kin;j++){
		Mrand[i].in_nodos[j]=M[i].in_nodos[j];
	}
	Mrand[i].out_nodos =new int [Mrand[i].kout];
	for (int j=0;j<Mrand[i].kout;j++){
		Mrand[i].out_nodos[j]=M[i].out_nodos[j];
	}
}
//--------------------
while (cambios < RELINKS){
	//elijo nodo1
	nodo1=gsl_rng_uniform_int(r,N);//elijo uno de los N nodos
	while(Mrand[nodo1].kin == 0) nodo1=gsl_rng_uniform_int(r,N);
	//printf("nodo1=%d con k1=%d \n",nodo1,M[nodo1].kin);
	//elijo nodo 2
	nodo2=gsl_rng_uniform_int(r,N); //elijo el segundo nodo
	while ((nodo2==nodo1) or (Mrand[nodo2].kin == 0)) nodo2=gsl_rng_uniform_int(r,N);
	//printf("nodo2=%d con k2=%d \n",nodo2,M[nodo2].kin);
	//elijo el padre del nodo1
	int k1=gsl_rng_uniform_int(r,Mrand[nodo1].kin);
	nodo1p=Mrand[nodo1].in_nodos[k1];
	int k1p;
		for (int j=0;j<Mrand[nodo1p].kout;j++) {
			if (nodo1 == Mrand[nodo1p].out_nodos[j]) k1p=j;
		}
	//printf("k1=%d nodo1p=%d y en out es el %d \n",k1,nodo1p,k1p);
	//elijo el padre del nodo2
	int k2=gsl_rng_uniform_int(r,Mrand[nodo2].kin);
	nodo2p=Mrand[nodo2].in_nodos[k2];
	int k2p;
		for (int j=0;j<Mrand[nodo2p].kout;j++) {
			if (nodo2 == Mrand[nodo2p].out_nodos[j]) k2p=j;
		}
	//printf("k2=%d nodo2p=%d y en out es el %d \n",k2,nodo2p,k2p);
	//tenemos que ver si estan ya cogidos
	int repe=0;
	for (int i=0;i<Mrand[nodo1].kin;i++){
		if (nodo2p == Mrand[nodo1].in_nodos[i]) repe=1;
	}
	for (int i=0;i<Mrand[nodo2].kin;i++){
		if (nodo1p == Mrand[nodo2].in_nodos[i]) repe=1;
	}
	
	if(repe==0){//si se pueden intercambiar,tenemos q cambiar la lista de nodos in de nodo1 y nodo2 y las de out de nodop1 y nodop2
		cambios ++;
	//	printf("cambios %d\n",cambios);
	//	printf("antes %d->%d %d->%d \n ahora %d->%d %d->%d\n",nodo1p,nodo1,nodo2p,nodo2,nodo2p,nodo1,nodo1p,nodo2);
		Mrand[nodo1].in_nodos[k1]=nodo2p;
		Mrand[nodo2].in_nodos[k2]=nodo1p;
		Mrand[nodo1p].out_nodos[k1p]=nodo2;
		Mrand[nodo2p].out_nodos[k2p]=nodo1;
	}
	


}
return;}

void randomize_over(nodo *M,int N,gsl_rng *r,int times) //Es similar solo que sobreescribe la amtriz original
{
	int RELINKS=times*N; //cambio al menos N veces antes de sacar el lambda
	int nodo1,nodo2;
	int nodo1p,nodo2p;
	int cambios=0;
	
	while (cambios < RELINKS){
		//elijo nodo1
		nodo1=gsl_rng_uniform_int(r,N);//elijo uno de los N nodos
		while(M[nodo1].kin == 0) nodo1=gsl_rng_uniform_int(r,N);
		//printf("nodo1=%d con k1=%d \n",nodo1,M[nodo1].kin);
		//elijo nodo 2
		nodo2=gsl_rng_uniform_int(r,N); //elijo el segundo nodo
		while ((nodo2==nodo1) or (M[nodo2].kin == 0)) nodo2=gsl_rng_uniform_int(r,N);
		//printf("nodo2=%d con k2=%d \n",nodo2,M[nodo2].kin);
		//elijo el padre del nodo1
		int k1=gsl_rng_uniform_int(r,M[nodo1].kin);
		nodo1p=M[nodo1].in_nodos[k1];
		int k1p;
			for (int j=0;j<M[nodo1p].kout;j++) {
				if (nodo1 == M[nodo1p].out_nodos[j]) k1p=j;
			}
		//printf("k1=%d nodo1p=%d y en out es el %d \n",k1,nodo1p,k1p);
		//elijo el padre del nodo2
		int k2=gsl_rng_uniform_int(r,M[nodo2].kin);
		nodo2p=M[nodo2].in_nodos[k2];
		int k2p;
			for (int j=0;j<M[nodo2p].kout;j++) {
				if (nodo2 == M[nodo2p].out_nodos[j]) k2p=j;
			}
		//printf("k2=%d nodo2p=%d y en out es el %d \n",k2,nodo2p,k2p);
		//tenemos que ver si estan ya cogidos
		int repe=0;
		for (int i=0;i<M[nodo1].kin;i++){
			if (nodo2p == M[nodo1].in_nodos[i]) repe=1;
		}
		for (int i=0;i<M[nodo2].kin;i++){
			if (nodo1p == M[nodo2].in_nodos[i]) repe=1;
		}
	
		if(repe==0){//si se pueden intercambiar,tenemos q cambiar la lista de nodos in de nodo1 y nodo2 y las de out de nodop1 y nodop2
			cambios ++;
		//	printf("cambios %d\n",cambios);
		//	printf("antes %d->%d %d->%d \n ahora %d->%d %d->%d\n",nodo1p,nodo1,nodo2p,nodo2,nodo2p,nodo1,nodo1p,nodo2);
			M[nodo1].in_nodos[k1]=nodo2p;
			M[nodo2].in_nodos[k2]=nodo1p;
			M[nodo1p].out_nodos[k1p]=nodo2;
			M[nodo2p].out_nodos[k2p]=nodo1;
		}
	


	}


	return;
}//----------------------------------------------------------------------------------------------------------------------

void rewire_directional(nodo *my_node,int N,gsl_rng *r, int *OK_label)
{
	int fail=0;
	int cambios,repe,direction_fail;

	int nodo1p,nodo1h; //nodop.out_nodes[k1p] -> nodoh
	int nodo2p,nodo2h; //nodoh.in_nodes[k1h] -> nodop 
	int node,nodeh,k2p,k2h;

	vector <int> padres;
	vector <int> hijos;
	
/*	for (int i=0;i<N;i++){*/
/*		printf("nivel M[%d]=%.2lf\n",i,my_node[i].s);*/
/*	}	*/

       //--------------------------------------------------------//	
	do {
		nodo1p=gsl_rng_uniform_int(r,N);
		
		
	} while (my_node[nodo1p].kout == 0); //no permitimos coger un sumidero
	int k1p=gsl_rng_uniform_int(r,my_node[nodo1p].kout);
		nodo1h=my_node[nodo1p].out_nodos[k1p];
		int k1h;
			for (int j=0;j<my_node[nodo1h].kin;j++) {
				if (nodo1p == my_node[nodo1h].in_nodos[j]) k1h=j;
			}
	//printf("%d(%.2lf) -> %d(%.2lf)\n",nodo1p,my_node[nodo1p].s,nodo1h,my_node[nodo1h].s);
/*	for (int i=0;i<my_node[nodo1p].kout;i++) printf("%d ",my_node[nodo1p].out_nodos[i]);*/
/*	printf("\n");*/
	//-------------------------------------------------------//
	
	//BUSQUEDA DIRIGIDA
	//Ahora vamos a elegir un nodop proporcionalmente a su parecido con nodo1p, y un nodoh entre sus vecinos out proporcionalmente a su parecido con nodo1h
/*	double total=0.; double pi; double x; double T; //temperatura */
/*	T=4.;*/
/*	for (int i=0;i<N;i++) if ((i != nodo1p) or (my_node[i].kout == 0))*/
/*	{*/
/*		pi= min(1.,1./(T*fabs(my_node[nodo1p].s - my_node[i].s)));*/
/*		total += pi;*/
/*		//printf("pi[%d]=%lf\n",i,pi);*/
/*				*/
/*	}*/
/*		if (total == 0.0) {//si no hay ningun nodo elegible nos piramos*/
/*			fail=1;*/
/*			goto EXIT;*/
/*		}*/
/*	x=gsl_rng_uniform(r)*total; //elijo un numero aleatorio y veo en q nodo cae*/
/*	*/
/*	total=0.; node=-1;*/
/*	while (total < x)  //no podemos coger nial mismo q el primero, ni ninguno q no tenga vecinos hacia afuera!*/
/*	{	*/
/*		node ++;*/
/*		if ( (node != nodo1p) or (my_node[node].kout==0)) total += min(1.,1./fabs(my_node[nodo1p].s - my_node[node].s));*/
/*		*/
/*	}*/
/*	nodo2p=node;*/
/*	printf("nodo2p=%d(%.2f):",node,my_node[node].s);*/
/*	//for (int i=0;i<my_node[nodo2p].kout;i++) printf("%d ",my_node[nodo2p].out_nodos[i]);*/
/*	printf("\n");*/
/*	//Ahroa escogemos uno de sus vecinos out q se aparezca a nodo1h*/
/*	total=0.;	*/
/*	for (int i=0;i<my_node[node].kout;i++){*/
/*		nodeh=my_node[node].out_nodos[i];*/
/*		pi= min(1.,1./(T*fabs(my_node[nodo1h].s - my_node[nodeh].s)));*/
/*		total += pi;*/
/*		//printf("pi[%d]=%lf\n",nodeh,pi);*/
/*	}*/
/*		if (total == 0.0) {//si no hay ningun nodo elegible nos piramos*/
/*			fail=1;*/
/*			goto EXIT;*/
/*		}*/
/*	x=gsl_rng_uniform(r)*total;*/
/*	*/
/*	total=0.; node=-1;*/
/*	while (total < x) */
/*	{	*/
/*		node ++;*/
/*		nodo2h=my_node[nodo2p].out_nodos[node];*/
/*		if ( nodo2h != nodo1h ) total += min(1.,1./fabs(my_node[nodo1h].s - my_node[nodeh].s));*/
/*		*/
/*	}*/
/*	printf("nodo2h=%d(%.2f)\n",nodo2h,my_node[nodo2h].s);*/
/*	k2p=node; //es en este en el orden q esta nodo2h en my_node[nodo2p].out_nodos[node]=nodo2h*/
/*		for (int j=0;j<my_node[nodo2h].kin;j++) {*/
/*			if (nodo2p == my_node[nodo2h].in_nodos[j]) k2h=j;*/
/*		}*/
	//Ahora tenemos que buscarlos entre si
	//-------------------------------------------------------------------------------------------------
	
	
	//BUSQUEDA POR NIVELES
	
	//miro a todos los del mismo nivel y elijo uno aleatoriamente, que no sea el nodo1p, claro . No puede ser uno con kout==0 !
	for (int i=0;i<N;i++) if ((i!= nodo1p) and ((int)my_node[nodo1p].s == (int)my_node[i].s) and (my_node[i].kout!=0) and (i != nodo1h)) padres.push_back(i);
	//printf("intentando elegir uno entre %d padres\n",padres.size());
	fail=0;
	if (padres.size() >0 ) node=gsl_rng_uniform_int(r,padres.size());
	else {
		fail=1;
		//printf("EXIT\n");
		goto EXIT;
	}
	nodo2p=padres[node]; //tengo uno aleatorio con el mismo nivel
	
	//miro todos los q salen de ese y seran elegibles todos los q no sean nodo1h
	for (int i=0;i<my_node[nodo2p].kout;i++) {
		nodeh=my_node[nodo2p].out_nodos[i];
		if ((my_node[nodo2p].out_nodos[i] != nodo1h) and ((int)my_node[nodo1h].s == (int)my_node[nodeh].s)) hijos.push_back(my_node[nodo2p].out_nodos[i]);
	}
	fail=0;
	//printf("intentando elegir uno entre %d hijos\n",hijos.size());
	if (hijos.size() > 0) node=gsl_rng_uniform_int(r,hijos.size());
	else {
		fail=1;
		//printf("EXIT\n");
		goto EXIT;
	}
	nodo2h=hijos[node];
	
		//no se nos puede olvidar ver en q indice estan los nodos para cambiarlos
		for (int i=0;i<my_node[nodo2p].kout;i++){
			if (nodo2h == my_node[nodo2p].out_nodos[i]) k2p=i;
		}
	
		for (int j=0;j<my_node[nodo2h].kin;j++) {
			if (nodo2p == my_node[nodo2h].in_nodos[j]) k2h=j;
		}
	
	direction_fail=0;
	//if ( (my_node[nodo2h].s < my_node[nodo1p].s ) or (my_node[nodo1h] < my_node[nodo2p]) ) direction_fail=1;
	if (my_node[nodo2p].s > my_node[nodo1h].s) direction_fail=1;
	if (my_node[nodo1p].s > my_node[nodo2h].s) direction_fail=1;
	//BUSQUEDA NO DIRIGIDA
	
	
	//CAMBIO EN LA TOPOLOGIA
	//Tenemos que revisar que los nuevos no esten ya linkeados, apra que no haya links dobles!! q eso cambiaria la conectividad
	repe=0;
	for (int i=0;i<my_node[nodo1p].kout;i++) if ( nodo2h == my_node[nodo1p].out_nodos[i]) repe=1; // Ya es su vecino!
	for (int i=0;i<my_node[nodo2p].kout;i++) if ( nodo1h == my_node[nodo2p].out_nodos[i]) repe=1; // Ya es su vecino!
	
	if ((repe ==0) and (direction_fail==0) ){//si se pueden intercambiar,tenemos q cambiar la lista de nodos in de nodo1 y nodo2 y las de out de nodop1 y nodop2
	//if ( repe==0){
		//printf("%d(%.2lf) -> %d(%.2lf) // ",nodo1p,my_node[nodo1p].s,nodo1h,my_node[nodo1h].s);
		//printf("%d(%.2lf) -> %d(%.2lf)\n",nodo2p,my_node[nodo2p].s,nodo2h,my_node[nodo2h].s);
		
			my_node[nodo1p].out_nodos[k1p]=nodo2h;
			my_node[nodo2p].out_nodos[k2p]=nodo1h;
			my_node[nodo1h].in_nodos[k1h]=nodo2p;
			my_node[nodo2h].in_nodos[k2h]=nodo1p;
			*OK_label=1;
		
		//printf("------------------------\n");
	}

	
	EXIT:
	//if (fail==1) printf("hubo un error, volvemos a elegir\n");
	//if (repe==1) printf("repe, try again!\n");
	padres.clear();
	hijos.clear();
	
	return;

}

//-----------------------------------------------------------------------------------------------------------------------
void erase_important_links(nodo *M,int N,edge **my_edge,int *how_many)
{
	//encontramos el link mas popular
	int most_important=0;
	int cont=0;
	for(int i=0;i<N;i++) {//TODO N por Nlinks
		if ((*my_edge)[i].importance > cont) {
			cont=(*my_edge)[i].importance;
			most_important=i;
		} 
	}
	printf("el link mas famoso es %d [%d -> %d]\n",most_important,(*my_edge)[most_important].in,(*my_edge)[most_important].out);
	//localizado el link, lo eliminamos tanto de la lista de edges como de M
	*how_many=(*my_edge)[most_important].importance;
	(*my_edge)[most_important].importance=0;
	int nodo_in,nodo_out;
	nodo_in=(*my_edge)[most_important].in;
	nodo_out=(*my_edge)[most_important].out;
	printf("antes:M[%d]=",nodo_in); for (int i=0;i<M[nodo_in].kout;i++) printf("%d ",M[nodo_in].out_nodos[i]);printf("\n");
	for(int i=0;i<M[nodo_in].kout;i++){
		if(M[nodo_in].out_nodos[i] == nodo_out) {
			int aux;
			M[nodo_in].out_nodos[i]=M[nodo_in].out_nodos[M[nodo_in].kout-1];			
			M[nodo_in].out_nodos[M[nodo_in].kout-1]=-1;
			//break;
		}
	}
	M[nodo_in].kout --;
	printf("ahora:M[%d]=",nodo_in); for (int i=0;i<M[nodo_in].kout;i++) printf("%d ",M[nodo_in].out_nodos[i]);printf("\n");
	
	
	printf("antes:M[%d]=",nodo_out); for (int i=0;i<M[nodo_out].kin;i++) printf("%d ",M[nodo_out].in_nodos[i]);printf("\n");
	for (int i=0;i<M[nodo_out].kin;i++){
		if (M[nodo_out].in_nodos[i]==nodo_in) {
			M[nodo_out].in_nodos[i]=M[nodo_out].in_nodos[M[nodo_out].kin -1];
			M[nodo_out].in_nodos[M[nodo_out].kin -1 ]=-1;
			//break;
		}
	}
	M[nodo_out].kin --;
	printf("ahora:M[%d]=",nodo_out); for (int i=0;i<M[nodo_out].kin;i++) printf("%d ",M[nodo_out].in_nodos[i]);printf("\n");
	
	
	return;
}


void erase_determined_links(nodo *M,int N,edge **my_edge,int who)
{
	//eliminamos el link "who"
	int nodo_in,nodo_out;
	nodo_in=(*my_edge)[who].in;
	nodo_out=(*my_edge)[who].out;
	//printf("antes:Mout[%d]=",nodo_in); for (int i=0;i<M[nodo_in].kout;i++) printf("%d ",M[nodo_in].out_nodos[i]);printf("\n");
	for(int i=0;i<M[nodo_in].kout;i++){
		if(M[nodo_in].out_nodos[i] == nodo_out) {
			M[nodo_in].out_nodos[i]=M[nodo_in].out_nodos[M[nodo_in].kout-1];			
			M[nodo_in].out_nodos[M[nodo_in].kout-1]=-1;
			//break;
		}
	}
	M[nodo_in].kout --;
	//printf("ahora:Mout[%d]=",nodo_in); for (int i=0;i<M[nodo_in].kout;i++) printf("%d ",M[nodo_in].out_nodos[i]);printf("\n");
	
	//printf("antes:Min[%d]=",nodo_out); for (int i=0;i<M[nodo_out].kin;i++) printf("%d ",M[nodo_out].in_nodos[i]);printf("\n");
	//if (nodo_in != nodo_out){	
	for (int i=0;i<M[nodo_out].kin;i++){
		if (M[nodo_out].in_nodos[i]==nodo_in) {
			M[nodo_out].in_nodos[i]=M[nodo_out].in_nodos[M[nodo_out].kin -1];
			M[nodo_out].in_nodos[M[nodo_out].kin -1 ]=-1;
			//break;
		}
	//}
	
	}//end if son diferentes, es q si son iguales quitaba dos veces links
	M[nodo_out].kin --;
	//printf("ahora:Min[%d]=",nodo_out); for (int i=0;i<M[nodo_out].kin;i++) printf("%d ",M[nodo_out].in_nodos[i]);printf("\n");
	
	
	return;
}

void how_many_to_remove(edge *my_edge,long double Nloops_rand, long double Nloops_exp,int contmax,int *who,int *how_many,int Nedges)
{
/*	printf("estoy dentro\n");*/
/*	return;*/
/*	//vamos a ir restando los que se eliminarian hasta q el numero de loops coincida mas o menos*/
/*	double *Links_in_order;*/
/*	int *order;*/
/*	Links_in_order=new double [Nedges];*/
/*	order=new int [Nedges];*/
/*	*/
/*	for(int i=0;i< Nedges; i++) Links_in_order[i]=my_edge[i].importance*1.;*/
/*	*/
/*	gsl_sort(Links_in_order, 1, Nedges);*/
/*	//for (int i=0;i<Nedges;i++) printf("imp[%d]=%f\n",i,Links_in_order[i]);*/
/*	*/
/*	int cont=Nedges-1;*/
/*	int cont_inverso=0;*/
/*	int ocupacion=1;*/
/*	vector <int> to_erase;*/
/*	//printf("cont=%d Nloops_exp=%Lf Nloops_rand=%Lf \n",cont,Nloops_exp,Nloops_rand);*/
/*	*/
/*	while( (Nloops_exp < Nloops_rand) && (cont >= 0 ) && (ocupacion>0) ) {//mientras haya mas en la red randomizada, pero no me salga del vector Links_in_order*/
/*		*/
/*		//printf("probamos con el %d:\n",cont);*/
/*		//printf("Nloops_exp=%Lf Nloops_rand=%Lf - %f/%d = ",Nloops_exp,Nloops_rand,Links_in_order[cont],contmax-1);*/
/*		*/
/*			Nloops_rand -= Links_in_order[cont]/(contmax*1.-1.);*/
/*			//printf("%Lf\n",Nloops_rand);*/
/*			*/
/*			if (Nloops_rand < Nloops_exp) {//significaria q nos hemos colado*/
/*				*/
/*				//printf("nos colariamos quitando, así q probamos con el siguiente menor: \n");*/
/*				//Deshacemos el cambio*/
/*				Nloops_rand += Links_in_order[cont]/(contmax*1.-1.);//primero hago q vuelva a tener el valor de antes*/
/*				*/
/*			}*/
/*		*/
/*			else {	//significaria q esta bien quitarlo*/
/*				//printf(":aceptado\n");*/
/*				to_erase.push_back(cont);*/
/*			}*/
/*		*/
/*		cont --;*/
/*		ocupacion=Links_in_order[cont];	*/
/*	}//end while buscando links a quitar*/
/*	*/
/*	//printf("acabamos con Nloops_rand=%Lf\n",Nloops_rand);*/
/*	*/
/*	printf("tendriamos q quitar:");*/
/*	for (int i=0;i<to_erase.size();i++) printf("%d, ",to_erase[i]);*/
/*	printf("\n");*/
/*	*/
/*	//(*who)=new int [to_erase.size()];*/
/*	int indice=0;*/
/*	for (int i=0;i<to_erase.size();i++){*/
/*		//para cada uno de los q vamos a quitar los tenemos q localizar*/
/*		for (int j=0;j<Nedges;j++){*/
/*			if(Links_in_order[to_erase[i]] == my_edge[j].importance){*/
/*				//lo hemos encontrado*/
/*				(*who)[indice]=j;*/
/*				indice ++;*/
/*			}*/
/*		}*/
/*	}*/
/*	*/
/*	printf("que significa q quitar:");*/
/*	for (int i=0;i<to_erase.size();i++) printf("%d, ",(*who)[i]);*/
/*	printf("\n");*/
/*	(*how_many)=to_erase.size();*/
/*	*/
	//estos generan algun fallo de memoria al liberarse parece
	//to_erase.clear();
	//delete [] Links_in_order;
	//delete [] order;
	
	return;
}
//----------------------------------------------------------------------------------------ES DE MUTICANONICAL
void mutar(nodo *M2,nodo *M1,int N,gsl_rng *r){
//por ahora ,lo que hace es copiarla, por ejemplo
for (int i=0;i<N;i++){
M2[i].kin=M1[i].kin;
M2[i].kout=M1[i].kout;
M2[i].in_nodos=new int [M2[i].kin];
M2[i].out_nodos=new int [M2[i].kout];
}
return;
}

void get_network(int n,nodo *M,double *w)
{

	
	FILE *vfich_error;
	vfich_error=fopen("dsaupderror.net","w");
	
	for (int i=0;i<n;i++){
		for (int j=0;j<M[i].kin;j++){
				fprintf(vfich_error,"%d %d\n",i,M[i].in_nodos[j]);
		}	
		
	}
	
	
	
/*	fprintf(vfich_error,"conectividades\n");*/
/*	for(int i=0;i<n;i++){*/
/*		fprintf(vfich_error,"Kin[%d]=%d , Kout[%d]=%d\n",i,M[i].kin,i,M[i].kout);*/
/*	}*/
	
/*	//veamos la densidad:*/
/*	printf("density:\n");*/
/*	for (int i=0;i<n;i++){*/
/*		fprintf(vfich_error,"%d %f\n",i,w[i]);*/
/*	}*/
	
/*	for (int i=0;i<n;i++){*/
/*		fprintf(vfich_error,"Nodo %d , se come a :",i);*/
/*		for (int j=0;j<M[i].kin;j++){*/
/*			fprintf(vfich_error,"%d (%f),",M[i].in_nodos[j],M[i].in_sig[j]);*/
/*		*/
/*		}*/
/*		fprintf(vfich_error,"\n");*/
/*		fprintf(vfich_error,"es comido por :");*/
/*		for (int j=0;j<M[i].kout;j++){*/
/*			fprintf(vfich_error,"%d ",M[i].out_nodos[j]);*/
/*		*/
/*		}*/
/*		fprintf(vfich_error,"\n");*/
/*	}*/
	
	fclose(vfich_error);
	return;
}


void print_red(nodo *M,int N)
{
	for (int i=0;i<N;i++){
		printf("M[%d]in:",i);
		for (int j=0;j<M[i].kin;j++){
			printf("%d ",M[i].in_nodos[j]);
		}
		printf("\n");
		printf("M[%d]out:",i);
		for (int j=0;j<M[i].kout;j++){
			printf("%d ",M[i].out_nodos[j]);
		}
		printf("\n");
	}
	
return;
}

unsigned long int fact(int n,int stop){
if (n>stop) return n*fact(n-1,stop);
else return 1;

}

void imp_grafo_levels(char *nombre,int N, nodo * my_node,int mode,gsl_rng *r)
{
	double coordx,Smax;
	Smax=0;
	for (int i=0;i<N;i++) if (my_node[i].s > Smax) Smax=my_node[i].s; 

	char color[70];
	FILE *vfich;
	char nombre_archivo[70];
	if (mode==1) sprintf(nombre_archivo,"%s_rand_niveles.net",nombre);
	else if (mode==0) sprintf(nombre_archivo,"%s_niveles.net",nombre);
	vfich=fopen(nombre_archivo,"w");

	fprintf(vfich,"*Vertices %d\n",N);
	for (int i=0;i<N;i++){
		colorea(my_node[i].s,i,color);
		coordx=gsl_rng_uniform(r);
		fprintf(vfich,"%d \"%d\" %.4f   %.4f   %.4f ic %s\n",i+1,i+1,coordx,my_node[i].s/Smax,0.5,color);
	}

	fprintf(vfich,"*Arcs\n");
	for (int i=0;i<N;i++){
		for (int j=0;j<my_node[i].kout;j++){
			//if (A[i][j]==1){
				fprintf(vfich," %2d     %2d   \n",my_node[i].out_nodos[j]+1,i+1);
			//}
		}
	}
	
	fclose(vfich);
	return;
}

void imp_grafo_components(char *nombre,int N, nodo * my_node,int mode,gsl_rng *r,int *sources,int *drains,int Nsources, int Ndrains)
{
	double coordx,Smax;
	Smax=0;
	for (int i=0;i<N;i++) if (my_node[i].s > Smax) Smax=my_node[i].s; 

	char color[70];
	FILE *vfich;
	char nombre_archivo[70];
	if (mode==1) sprintf(nombre_archivo,"%s_rand_compo.net",nombre);
	else if (mode==0) sprintf(nombre_archivo,"%s_compo.net",nombre);
	vfich=fopen(nombre_archivo,"w");

	fprintf(vfich,"*Vertices %d\n",N);
	for (int i=0;i<N;i++){
		define_tipo(sources,Nsources,drains,Ndrains,i,color);
		coordx=gsl_rng_uniform(r);
		fprintf(vfich,"%d \"%d\" %.4f   %.4f   %.4f ic %s\n",i+1,i+1,coordx,my_node[i].s/Smax,0.5,color);
	}

	fprintf(vfich,"*Arcs\n");
	for (int i=0;i<N;i++){
		for (int j=0;j<my_node[i].kout;j++){
			//if (A[i][j]==1){
				fprintf(vfich," %2d     %2d   \n",my_node[i].out_nodos[j]+1,i+1);
			//}
		}
	}
	
	fclose(vfich);
	return;
}

void define_tipo(int *sources,int Nsources,int* drains,int Ndrains,int nodo,char *color)
{
	int tipo=-1;
	for (int i=0;i<Nsources;i++){
	
		if (nodo==sources[i]) {
			tipo=0;
			printf("el nodo %d es base\n",i);
		}
	}
	printf("\n");
	
	for (int i=0;i<Ndrains;i++){
	printf("%d ",drains[i]);
		if(nodo==drains[i]) {
			printf("el nodo %d es apice\n",i);
			if (tipo=-1) tipo=1;
			else tipo=2;
		}
	}
	
	if (tipo==0) sprintf(color,"%s","ForestGreen");
	else if (tipo==1) sprintf(color,"%s","WildStrawberry");
	else if (tipo==2) sprintf(color,"%s","BurntOrange");
	else if (tipo==-1) sprintf(color,"%s","Yellow");
	return;
}

void colorea(float s,int nodo,char *color)
{
	int NIVEL;
	//printf("coloreo?\n");
	NIVEL=(int)s;
	switch(NIVEL){ //para poner el bombre al archivo
    case 1:
    sprintf(color,"%s","ForestGreen");
    break;

    case 2:
    sprintf(color,"%s","LimeGreen");
    break;
   
    case 3:
    sprintf(color,"%s","SpringGreen");
    break;
   
    case 4:
    sprintf(color,"%s","Yellow");
    break;
   
    case 5:
    sprintf(color,"%s","BurntOrange");
    break;
   
    case 6:
    sprintf(color,"%s","WildStrawberry");
    break;
   
    case 7:
    sprintf(color,"%s","RoyalBlue");
    break;
   
    case 8:
    sprintf(color,"%s","RedViolet");
    break;
   
    case 9:
    sprintf(color,"%s","OliveGreen");
    break;
   
    case 10:
    sprintf(color,"%s","Gray15");
    break;   
	}

	if (NIVEL > 10) {
		sprintf(color,"%s","Gray15");
	}

	return;
}






















