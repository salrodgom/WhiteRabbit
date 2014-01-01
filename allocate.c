#include "general.h"
//--------------------LEER RED-----------------------------------------------//
nodo * leer_red_new(int *N,char *nombre){
int nodoi,nodoj,indice,cosa,cont,ord_c;
float qq;
double xvar1;
char signo[5];
int signoint;
char nombre_arch[20];
char basura[4];
vector <vector <int> > aux_in; //auxiliar de entrada para cada nodo
vector <vector <int> > sig_in;
vector <vector <int> > aux_out; //auxiliar de salida para cada nodo
vector <vector <int> > sig_out;
int n;
FILE *fichero;

*N=0;// aun no sabemos cuantos nodos tiene
sprintf(nombre_arch,"%s.red",nombre); //generamos el nombre del archivo a usar
fichero=fopen(nombre_arch,"r");
	while(!feof(fichero)) { //sacamos el numero de nodos
	fscanf(fichero,"%s\n",basura);
	if( (strcmp(basura,"+")) or (strcmp(basura,"-"))) {
	   fscanf(fichero,"%lf\n", &xvar1);
	   nodoi=(int)xvar1;
	   //printf("leido %lf %d\n",xvar1,nodoi);
	   if (nodoi>*N){*N=nodoi;}
	   }
	}
	aux_in.resize(*N+1);//creo espacio para guardar las listas
	sig_in.resize(*N+1);
	aux_out.resize(*N+1);//creo espacio para guardar las listas
	sig_out.resize(*N+1);
	rewind(fichero);
	//printf("#ya tenemos N=%d\n",*N);
	while(!feof(fichero)) { // creamos la red
		fscanf(fichero,"%lf\n", &xvar1);
        	nodoi=(int)xvar1-1;
		//printf("nodoi=%i   ",nodoi);
       		fscanf(fichero,"%lf\n", &xvar1);
      		nodoj=(int)xvar1-1;
      		//printf("nodoj=%i   ",nodoj);
      		fscanf(fichero,"%s\n",signo);
      		//printf("signo %s\n",signo);
      		if (strcmp(signo,"-")== 0) signoint=-1;
      		else if (strcmp(signo,"+")== 0) signoint=1;
      		else if (strcmp(signo,"0")== 0) signoint=-1;//esto tengo q mirarlo!!!
      		else (signoint=8);
      		//printf("signo=%d\n",signoint);
      		//printf("itnento meter en auxin[%d] auxout[%d] de %d:\n",nodoi,nodoj,*N+1);     		
      		aux_in[nodoi].push_back(nodoj);
      		sig_in[nodoi].push_back(signoint);
      		aux_out[nodoj].push_back(nodoi);
      		sig_out[nodoj].push_back(signoint);
      		
      		signoint=0;
      		//printf("el nodo %d se come a ",nodoi);
      		//for (int i=0;i<aux_in[nodoi].size();i++) printf("%d ",aux_in[nodoi][i]);
      		//printf("\n");
      		
   	}
//printf("va bien la asignacion? \n");
//return;
/*for(int i=0;i<*N;i++){*/
/*printf("M[%d]in: ",i);*/
/*	for (int j=0;j<aux_in[i].size();j++){*/
/*		//printf("%d*",sig_in[i][j]);*/
/*		printf("%d ",aux_in[i][j]);*/
/*	}*/
/*	printf("\n");*/
/*printf("M[%d]out: ",i);*/
/*	for (int j=0;j<aux_out[i].size();j++){*/
/*		//printf("%d*",sig_out[i][j]);*/
/*		printf("%d ",aux_out[i][j]);*/
/*	}*/
/*	printf("\n");*/
/*}*/


n=*N;
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//printf("#N=%d n=%d\n",*N,n);

	nodo *my_node;
	my_node = new nodo [*N];
	for (int i=0;i<*N;i++) { my_node[i].kout=0; my_node[i].label=1; my_node[i].myself=i;}//pongo todas las conectividades out a cero

//printf("lo metemos dentro de la final\n");
for (int i=0;i<*N;i++){
	my_node[i].kin=aux_in[i].size();
	my_node[i].in_nodos=new int [my_node[i].kin];
	my_node[i].in_sig=new double [my_node[i].kin];
//printf("Kin[%d]=%d %d \n ",i,aux_in[i].size(),my_node[i].kin);
//printf("M[%d]in=",i);
	for (int j=0;j<my_node[i].kin;j++){
		my_node[i].in_nodos[j]=aux_in[i][j];
		my_node[i].in_sig[j]=sig_in[i][j];
//		printf("%d*(%.0f) ",my_node[i].in_nodos[j],my_node[i].in_sig[j]);
	}
//	printf("\n");
	
	my_node[i].kout=aux_out[i].size();
//printf("kout=%d %d \n",aux_out[i].size(),my_node[i].kout);
//printf("M[%d]out=",i);
	my_node[i].out_nodos=new int [my_node[i].kout];
	my_node[i].out_sig=new double [my_node[i].kout];
	for (int j=0;j<my_node[i].kout;j++){
		my_node[i].out_nodos[j]=aux_out[i][j];
		my_node[i].out_sig[j]=sig_out[i][j];
//		printf("%d*(%.0f) ",my_node[i].out_nodos[j],my_node[i].out_sig[j]);
	}
//	printf("\n");
}
//printf("hemos copiado --------------------------------------------  \n");
/*for (int i=0;i<*N;i++){*/
/*printf("Nodo %d , Kin=%d se come a :",i,my_node[i].kin);*/
/*	for (int j=0;j<my_node[i].kin;j++){*/
/*		printf("(%.0f)*%d ",my_node[i].in_sig[j],my_node[i].in_nodos[j]+1);*/
/*		*/
/*	}*/
/*	printf("\n");*/
/*printf("es comido por %d bichos:",my_node[i].kout);*/
/*for (int j=0;j<my_node[i].kout;j++){*/
/*		printf("(%.0f)*%d ",my_node[i].out_sig[j],my_node[i].out_nodos[j]+1);*/
/*		*/
/*	}*/
/*	printf("\n");*/
/*}*/
/*printf("terminamos?\n");*/

	aux_in.clear();
	aux_out.clear();
	
	//for (int i=0;i<N;i++) my_node[i].label=i;
	
	return my_node;
}//------------------------------------------------------------------------//

nodo * copiar_nodos(nodo * my_original,int N)
{
	nodo *my_copy;
	my_copy= new nodo [N];
	
	for (int i=0;i<N;i++){
		my_copy[i].kin=my_original[i].kin;
		my_copy[i].kout=my_original[i].kout;
		my_copy[i].label=my_original[i].label;
		my_copy[i].myself=my_original[i].myself;
		my_copy[i].s=my_original[i].s;
		my_copy[i].in_nodos=new int [my_copy[i].kin];
		my_copy[i].in_sig=new double [my_copy[i].kin];
		for (int j=0;j<my_copy[i].kin;j++){
			my_copy[i].in_nodos[j]=my_original[i].in_nodos[j];
			my_copy[i].in_sig[j]=my_original[i].in_sig[j];
		}
		my_copy[i].out_nodos=new int [my_copy[i].kout];
		my_copy[i].out_sig=new double [my_copy[i].kout];
		for (int j=0;j<my_copy[i].kout;j++){
			my_copy[i].out_nodos[j]=my_original[i].out_nodos[j];
			my_copy[i].out_sig[j]=my_original[i].out_sig[j];
		}	
	}	
	return my_copy;
}




//-------------------------------------------------------------------------//
nodo * leer_matriz(int *Na,int *Np,char *nombre) //TODO lee "matrices"" no dirigidas!!!
{

	FILE *fichero;
	char nombre_arch[30];
	sprintf(nombre_arch,"%s.red",nombre); //generamos el nombre del archivo a usar
	fichero=fopen(nombre_arch,"r");
	
	nodo *my_node=NULL;
	
	char *bufferlinea;
	int TOPE=1000;//El N maximo de la red
	bufferlinea=(char*)malloc(sizeof(char)*TOPE);
	
	vector <vector <int> > aux_out; //punto de vista de los pollinators, de ellos van a las plantas
	vector <vector <double> > sig_out;
	vector <vector <int> > aux_in;//punto de vista de las plantas, q pollinators vienen a ellas
	vector <vector <double> > sig_in; 
	vector <int> aux;
	
	
	printf("tamaño de entero=%d\n",sizeof(int));
	
	int num_lineas=0;
	int num_columnas=0;
		
	while(fgets(bufferlinea,TOPE,fichero)!=NULL){//mientras no llegue leyendo al final del fichero, q es cuando da NULL
		
		//digo q voy por una linea
		
		int nums_now, bytes_now;
		int bytes_consumed = 0, nums_read = 0;
		
		int aux2;
		
		aux_out.resize(num_lineas+1);//añadimos una fila mas
		sig_out.resize(num_lineas+1);// a este tambien, y dejamos las anteriores
		
		num_columnas=0;
		while ( ( nums_now = sscanf( bufferlinea + bytes_consumed, "%d%n", &aux2, &bytes_now )) > 0 ) {
    			bytes_consumed += bytes_now;
    			nums_read += nums_now;
    			
    			aux.push_back(aux2);
    			num_columnas++;
    			
    			if (aux2 != 0) {
    				aux_out[num_lineas].push_back(nums_read-1);
    				sig_out[num_lineas].push_back(aux2);
    			}

    			
		}// leo una linea------------------------------------------------------------------------------
		
		
		for (int i=0;i<aux.size();i++){
			printf("%d ",aux[i]);
		}
		printf("\n");
		
		num_lineas++;
		aux.resize(0);
		
		
	} // termino de leer el archivo--------------------------------------------------------------------------------------------

printf("--------------------------------------------\n");
	aux_in.resize(num_lineas);
	sig_in.resize(num_lineas);
	for (int i=0;i<num_lineas;i++){
		for (int j=0;j<aux_out[i].size();j++){
			//printf("%d(%.0lf) ",aux_out[i][j],sig_out[i][j]);
			aux_in[aux_out[i][j]].push_back(i);
			sig_in[aux_out[i][j]].push_back(sig_out[i][j]);
		}
		//printf("\n");
	}
printf("--------------------------------------------\n");
	//como estas son no dirigidas me vale con copiar todo a la aprte out:


	for (int i=0;i<num_lineas;i++){
		for (int j=0;j<aux_out[i].size();j++){
			printf("%d(%.0lf) ",aux_out[i][j],sig_out[i][j]);
		}
		printf("\n");
	}

printf("num lineas=%d num columnas=%d\n",num_lineas,num_columnas);

	*Na=num_lineas;//El numero de pollinators es el numero de lineas
	*Np=num_columnas;
	int N = max(*Na,*Np);
	//El numero de plantas es el numero de columnas
	my_node = new nodo [N]; //estos son los num_lineas pollinators del ecosistema
	for (int i=0;i<N;i++) { my_node[i].kout=0; my_node[i].kin=0; my_node[i].label=1; my_node[i].myself=i; }//pongo todas las conectividades out a cero
	
for (int i=0;i<N;i++){
// Esta direccion es la red desde las plantas
	my_node[i].kin=aux_in[i].size();
	my_node[i].in_nodos=new int [my_node[i].kin];
	my_node[i].in_sig=new double [my_node[i].kin];
//		printf("Kin[%d]=%d %d \n ",i,aux_in[i].size(),my_node[i].kin);
//		printf("M[%d]in=",i);
	for (int j=0;j<my_node[i].kin;j++){
		my_node[i].in_nodos[j]=aux_in[i][j];
		my_node[i].in_sig[j]=sig_in[i][j];
//		printf("%d*(%.0f) ",my_node[i].in_nodos[j],my_node[i].in_sig[j]);
	}
//	printf("\n");
//esta relacion es la red desde los animales 
	my_node[i].kout=aux_out[i].size();
//		printf("kout=%d %d \n",aux_out[i].size(),my_node[i].kout);
//		printf("M[%d]out=",i);
	my_node[i].out_nodos=new int [my_node[i].kout];
	my_node[i].out_sig=new double [my_node[i].kout];
	for (int j=0;j<my_node[i].kout;j++){
		my_node[i].out_nodos[j]=aux_out[i][j];
		my_node[i].out_sig[j]=sig_out[i][j];
//		printf("%d*(%.0f) ",my_node[i].out_nodos[j],my_node[i].out_sig[j]);
	}
//	printf("\n");
}



printf("hemos copiado --------------------------------------------  \n");
	
	for (int i=0;i<N;i++){
printf("Nodo %d , Kin=%d se come a :",i,my_node[i].kin);
	for (int j=0;j<my_node[i].kin;j++){
		printf("%d(%.0f) ",my_node[i].in_sig[j],my_node[i].in_nodos[j]+1);
		
	}
	printf("\n");
printf("es comido por %d bichos:",my_node[i].kout);
for (int j=0;j<my_node[i].kout;j++){
		printf("%d(%.0f)",my_node[i].out_sig[j],my_node[i].out_nodos[j]+1);
		
	}
	printf("\n");
}
	
	aux_out.clear();
	aux_in.clear();
	sig_out.clear();
	sig_in.clear();
	aux.clear();
	
	return my_node;
}//-------------------------------------------------------------------------------------------

//------------ RED RANDOM ---------------------------------------------------------------------
nodo * red_random_new(int N,double mu,gsl_rng *r,int directed,int weighted)
{
	//reservo memoria para la red
	nodo *my_node;
	nodo *my_node_sim;
	my_node=new nodo [N];
/*	*my_node =(nodo *) malloc(sizeof(nodo) *N);*/
/*	if (*my_node==NULL) {printf("No se pudo reservar\n"); exit (1);}*/

	//if (directed == 0) mu=mu/2.;
	
	
	int *cont,nodov;
	double p;
	vector <int> elegidos;
	cont=new int [N];
	
	for (int i=0;i<N;i++) cont[i]=0; //pongo a cero los contadores de k out de los nodos

	for (int i=0;i<N;i++) {my_node[i].kout=0; my_node[i].label=1; my_node[i].myself=i;}//pongo todas las conectividades out a cero

        //-------HACEMOS LA RED-------------------------
again:
for (int i=0;i<N;i++){

	for (int j=0;j<N;j++){
		p=gsl_rng_uniform(r);
		if(((N*p)< mu) and (i!=j)) {
			elegidos.push_back(j);//lo metemos a la lista de especies comidas
			//for(int k=0;k<elegidos.size();k++) printf("%d\n ",elegidos[k]);
			my_node[j].kout += 1; //hay alguien que se lo come
		}
		
	}
	my_node[i].kin=elegidos.size();
	my_node[i].in_nodos=new int [my_node[i].kin];
	for (int l=0;l<elegidos.size();l++){
		my_node[i].in_nodos[l]=elegidos[l]; //copiamos el vector
	}
	
	elegidos.clear();
}//------------------------------------------
//printf("terminamos\n");
for (int i=0;i<N;i++) { //preparamos los vectores de nodos out
	my_node[i].out_nodos=new int [my_node[i].kout];
	for (int j=0;j<my_node[i].kout;j++) {
		my_node[i].out_nodos[j]=0;
	}
}
//printf("lista out\n");
for (int i=0;i<N;i++) {
	for (int j=0;j<my_node[i].kin;j++){
		nodov=my_node[i].in_nodos[j]; //este es el nodo victima
		my_node[nodov].out_nodos[cont[nodov]]=i; //ponemos al que se lo come en su vector de outs
		cont[nodov]+=1;//sumamos uno al contador del nodo victima
	}
}
//añado dos filas de 1 que no aportan nada, pero mantienen la esturctura
for (int i=0;i<N;i++){
	my_node[i].out_sig=new double [my_node[i].kout];
	my_node[i].in_sig=new double [my_node[i].kin];
	for(int j=0;j<my_node[i].kin;j++){
		
		if (weighted == 1) p_random(my_node[i].in_sig,my_node[i].kin,r);//rellenamos my_node[i].in_sig con kin numeros q sumen 1
		else my_node[i].in_sig[j]=1;
	}
	for(int j=0;j<my_node[i].kout;j++){
		my_node[i].out_sig[j]=1;
	}
}

double kmed=0.;
for(int i=0;i<N;i++) kmed += my_node[i].kin;

kmed = kmed/(N*1.);
//printf("kmed=%lf\n",kmed);

/*for (int i=0;i<N;i++){*/
/*printf("#Nodo %d , se come a :",i);*/
/*	for (int j=0;j<my_node[i].kin;j++){*/
/*		//printf("%d (%f),",my_node[i].in_nodos[j],my_node[i].in_sig[j]);*/
/*		printf("#%d ",my_node[i].in_nodos[j]);*/
/*	}*/
/*	printf("\n");*/
/*printf("#es comido por :");*/
/*for (int j=0;j<my_node[i].kout;j++){*/
/*		printf("#%d ",my_node[i].out_nodos[j]);*/
/*		*/
/*	}*/
/*	printf("#\n");*/
/*}*/
	elegidos.clear();
	delete [] cont;
	
	if (directed == 0) {
		
		my_node_sim=simetrizar(my_node,N);
		liberar_nodos(my_node,N);
		return my_node_sim;
		
	}	
	//for(int i=0;i<N;i++) my_node[i].label=i;
	else { 
		
		return my_node;
	}
	
}//-----------------------------------------------------------------------------------
//---------
void p_random(double *p,int NSTATES,gsl_rng *r) //genera una lsita de numeros que suman 1
{
   p[0]=0.;

    for(int j=1; j<NSTATES; j++)
        p[j]=gsl_rng_uniform(r);
        
    // ordena la lista
    gsl_sort(p, 1, NSTATES); // el 1 es porque los coge de 1 en uno

    // calcula las diferencias
    for(int j=0; j<(NSTATES-1); j++)
        p[j]=p[j+1]-p[j];

    p[NSTATES-1]=1.-p[NSTATES-1];

    return;
}//----------------------------------------------------------------------------------

//------- RED CONFIGURACIONAL ------------------------------------------------------
//------- RED CONFIGURACIONAL -------------------------------------------------------------------------------
nodo * red_config(int N,double mu,double sigma,gsl_rng *r,int directed,int weighted,char *type)
{
retry:
	int double_links;
	int autoloop;
	int veces=0;
	int SAFE=100000;
	
	//printf("red configuracional %s\n",type);
	nodo *my_node;
	my_node=(nodo *)malloc(N*sizeof(nodo));
/*	//vamos a generar una red segun el modelo configuracional: necesitamos una serie de nodos con sus conectividades, que luego iremos eligiendo:*/
/*	//sacamos la listya de conectividad: N extracciones de una distribucion tipo "Type"*/
	vector <int> nodes;
	
	vector <vector <int> > aux_in;
	vector <vector <int> > aux_out;
	
	aux_in.resize(N);
	aux_out.resize(N);
	
	printf("red %s\n",type);
	int k;
	int ktot=0.;
	for(int i=0;i<N;i++) {
		if (strcmp(type,"GA")==0) { //Gaussian degree dsitributed
			k= gsl_ran_gaussian(r,3.)+mu;
			do{
				k=(int)gsl_ran_gaussian(r,3.)+mu;
						
			}while (k <= 0);
			//printf("k=%d\n",k);
			ktot += k;
		}
		else if (strcmp(type,"PO")==0) { //Poisson degree distributed
			k=gsl_ran_poisson(r,mu);
			do{
				k=gsl_ran_poisson(r,mu);
			
			}while (k <= 0);
			//printf("k=%d\n",k);
			ktot += k;
		}
		else if (strcmp(type,"SC")==0){ //scale free degree distributed, gamma=1.5
			k=power_law(N,3.,r);
			do {
				k=power_law(N,3.,r);
			} while ((k > (int) sqrt(N)) or (k<=0));
			ktot += k;
		}	
		
		//printf("node %d k=%d\n",i,k);
		if (directed == 0) {
			for (int j=0;j<k;j++){
				//printf("meto %d\n",i);
				nodes.push_back(i);
			}
		}
		else {
			for (int j=0;j<2*k;j++){
				//printf("meto %d\n",i);
				nodes.push_back(i);
			}		
		}
	}
	printf("kmed=%lf\n",ktot/(N*1.));
	//printf("#el numero de links es %d\n",nodes.size());
	//for (int i=0;i<nodes.size();i++) printf("%d ",nodes[i]);
	//printf("#\n");
	

	int cont=0;
	//printf("nlinks=%d npair=%d\n",nodes.size(),pair_total);	
	while ((nodes.size()>1) and (veces < SAFE)){
	
	//printf("#eleccion %d\n",cont);
		//vamos sorteando y creando los links
		int nodo1,nodo2; nodo2=-1;
	again:
		nodo1=gsl_rng_uniform_int(r,nodes.size());
		//printf("#nodo1=%d, ",nodes[nodo1]);
		do {
			nodo2=gsl_rng_uniform_int(r,nodes.size());
			
		} while ( (nodo2 == nodo1) );
		veces++;
		if (veces ==SAFE) goto exit;
		//printf("#nodo2=%d ",nodes[nodo2]);
		//printf("es decir %d -> %d size=%d\n ",nodes[nodo2],nodes[nodo1],nodes.size());

		//vigilemos que no esté repetido
		double_links=0;
		for (int i=0;i<aux_in[nodes[nodo1]].size();i++){
			if (nodes[nodo2] == aux_in[nodes[nodo1]][i]) {
				double_links=1;
				//if (nodes.size()>2) {
					//printf("elijo otra vez\n");
					goto again;	
				//}
			}
		}
		autoloop=0;
		if (nodes[nodo2]==nodes[nodo1]) {
			autoloop=1;
				goto again;
		}
		
		
		aux_in[nodes[nodo1]].push_back(nodes[nodo2]);
		aux_out[nodes[nodo2]].push_back(nodes[nodo1]);
	
		
		
		if (directed == 0) {//si es no dirigida , incluimos los reversos tb
		
			aux_in[nodes[nodo2]].push_back(nodes[nodo1]); //el 2 tb se come al 1
			aux_out[nodes[nodo1]].push_back(nodes[nodo2]); //el 1 es comido por el 2
		}
		else {//si es dirigida no lo tenemos preparado
			//printf("en construccion\n");
			//return my_node;
		}
		
		int max=nodo1;
		int min;
		if (nodo2 > nodo1) { max=nodo2; min=nodo1;}
		else {min=nodo2;}
		//printf("tamaño:%d ; eliminamos el elemento %d y luego el %d",nodes.size(),max,min);
		nodes.erase (nodes.begin()+max);//eliminamos los elementos par q no se puedan eleigr otra vez
		nodes.erase (nodes.begin()+min);
		//printf("y el tamaño es %d\n",nodes.size());
		//for (int i=0;i<nodes.size();i++) printf("%d ",nodes[i]);printf("\n");
		cont ++;
	}//mientras queden nodos sin 
	//printf("#parejas elegidas\n");
	printf("autoloop=%d 2links=%d veces=%d\n",autoloop,double_links,veces);
	exit:
	if (veces==SAFE) {//reintentarla entera de nuevo
		aux_in.clear();
		aux_out.clear();
		nodes.clear();
		goto retry;
	}
	//fflush(stdout);
	for(int i=0;i<N;i++) {
	  //printf("nodo %d:",i);
	
		my_node[i].kin=aux_in[i].size();
		//printf("se come a %d bichos y ",my_node[i].kin);
		my_node[i].kout=aux_out[i].size();
		//printf("es comido por %d bichos\n",my_node[i].kout);
		my_node[i].in_nodos=new int [my_node[i].kin];
		my_node[i].in_sig=new double [my_node[i].kin];
		my_node[i].out_nodos=new int [my_node[i].kout];
		my_node[i].out_sig=new double [my_node[i].kout];
		for (int j=0;j<my_node[i].kin;j++){
			my_node[i].in_nodos[j]=aux_in[i][j];
			my_node[i].in_sig[j]=1;
		}
		for (int j=0;j<my_node[i].kout;j++){
			my_node[i].out_nodos[j]=aux_out[i][j];
			my_node[i].out_sig[j]=1;
		}
	}
/*	printf("copiado\n");*/
	
	
	if( nodes.size() > 0) {
		nodes.clear();
		//printf("borro pq es mayor q cero");
	}
	
	for (int i=0;i<N;i++) {
		if (aux_in[i].size() > 0) aux_in[i].clear();
		if (aux_out[i].size() > 0) aux_out[i].clear();
	}
	aux_in.clear();
	aux_out.clear();
	
	for(int i=0;i<N;i++){
		for(int j=0;j<my_node[i].kin;j++){
		int nodoc=my_node[i].in_nodos[j];
			for(int k=j+1;k<my_node[i].kin;k++){
				if (nodoc == my_node[i].in_nodos[k]) printf("nodo %d:%d\n",nodoc,my_node[i].in_nodos[k]);
			}
		}
	}
	
		for (int i=0;i<N;i++){
	printf("Nodo[%d].kin: %d :",i,my_node[i].kin);
		for (int j=0;j<my_node[i].kin;j++){
			printf("%d ",my_node[i].in_nodos[j]);
		
		}
		printf("\n");
	printf("Nodo[%d],kout:%d :",i,my_node[i].kout);
		for (int j=0;j<my_node[i].kout;j++){
			printf("%d ",my_node[i].out_nodos[j]);
		
		}
		printf("\n");
	}
	double km;
	for (int i=0;i<N;i++) {
		km += my_node[i].kout;
	}
	km=km/N;
	printf("<k>=%lf\n",km);
	
/*	//control conexion*/
/*	for (int i=0;i<N;i++){*/
/*		for (int j=0;j<my_node[i].kin;j++){		*/
/*		*/
/*	int k1=j;//elijo uno de los padres del nodo1*/
/*	int nodo1p=my_node[i].in_nodos[k1];//el padre del nodo1*/
/*	//printf("nodo1p=%d\n",nodo1p);*/
/*	int k1p;//el sitio q ocupa el nodo1 en el out_nodos del padre*/
/*	//printf("el padre de %d=%d debe estar entre :",nodo1p,i);*/
/*		for (int k=0;k<my_node[nodo1p].kout;k++) {*/
/*	//	printf("%d ",my_node[nodo1p].out_nodos[k]);*/
/*			if (i == my_node[nodo1p].out_nodos[k]) k1p=k;*/
/*		}*/
/*	//	printf("\n");*/
/*	if((k1p > my_node[nodo1p].kout) or (k1p < 0)) printf("chungo!\n");*/
/*	//printf("#k1=%d nodo1p=%d y en out es el %d \n",k1,nodo1p,k1p);*/
/*	*/
/*		}//end for j*/
/*	}//end forN*/
	

	
	return my_node;
}
//----------------------------------------------------------------------------
//------- RED CONFIG2 -------------------------------------------------------------------------------
nodo * red_config2(int N,double mu,double sigma,gsl_rng *r,int directed)
{
	//#define funcion(mu) gsl_ran_poisson(r,mu)
/*	//#define funcion(mu) gsl_ran_gaussian(r,1.)*mu*/
/*	//#define funcion(mu) power_law(N,mu) //ojo, aqui mu es el exponente*/
/*	*/
	nodo *my_node;
	my_node=(nodo *)malloc(N*sizeof(nodo));
/*	//vamos a generar una red segun el modelo configuracional: necesitamos una serie de nodos con sus conectividades, que luego iremos eligiendo:*/
/*	//sacamos la listya de conectividad: N extracciones de una distribucion tipo "Type"*/
	vector <int> nodes;
	
	vector <vector <int> > aux_in;
	vector <vector <int> > aux_out;
	
	aux_in.resize(N);
	aux_out.resize(N);
	
	int k;
	for(int i=0;i<N;i++) {
		k= gsl_ran_gaussian(r,10.)+mu;
		//k=gsl_ran_poisson(r,mu);
		do{
			k=(int)gsl_ran_gaussian(r,10.)+mu;
			//k=gsl_ran_poisson(r,mu);
			
		}while (k <= 0);
		//printf("k=%d\n",k);
		
		//k=power_law(N,1.5,r);
		//do {
		//	k=power_law(N,1.5,r);
		//} while ((k > (int) sqrt(N)) or (k<=0));
		
		//printf("node %d k=%d\n",i,k);
		for (int j=0;j<k;j++){
		//printf("meto %d\n",i);
			nodes.push_back(i);
		}
	}
	//printf("#el numero de links es %d\n",nodes.size());
	//for (int i=0;i<nodes.size();i++) printf("%d ",nodes[i]);
	//printf("#\n");
	

	int cont=0;
	//printf("nlinks=%d npair=%d\n",nodes.size(),pair_total);	
	while (nodes.size()>1){
	
	//printf("#eleccion %d\n",cont);
		//vamos sorteando y creando los links
		int nodo1,nodo2; nodo2=-1;
		nodo1=gsl_rng_uniform_int(r,nodes.size());
		//printf("#nodo1=%d, ",nodo1);
		do {
			nodo2=gsl_rng_uniform_int(r,nodes.size());
			
		} while ( nodo2 == nodo1 );
		//printf("#nodo2=%d ",nodo2);
		//printf("es decir %d -> %d\n",nodes[nodo2],nodes[nodo1]);

		//vigilemos que no esté repetido

		aux_in[nodes[nodo1]].push_back(nodes[nodo2]);
		aux_out[nodes[nodo2]].push_back(nodes[nodo1]);
		
		if (directed == 0) {//si es no dirigida , incluimos los reversos tb
		
			aux_in[nodes[nodo2]].push_back(nodes[nodo1]); //el 2 tb se come al 1
			aux_out[nodes[nodo1]].push_back(nodes[nodo2]); //el 1 es comido por el 2
		}
		
		int max=nodo1;
		int min;
		if (nodo2 > nodo1) { max=nodo2; min=nodo1;}
		else {min=nodo2;}
		//printf("tamaño:%d ; eliminamos el elemento %d y luego el %d",nodes.size(),max,min);
		nodes.erase (nodes.begin()+max);//eliminamos los elementos par q no se puedan eleigr otra vez
		nodes.erase (nodes.begin()+min);
		//printf("y el tamaño es %d\n",nodes.size());
		//for (int i=0;i<nodes.size();i++) printf("%d ",nodes[i]);printf("\n");
		cont ++;
	}//mientras queden nodos sin 
	//printf("#parejas elegidas\n");
	//fflush(stdout);
	for(int i=0;i<N;i++) {
	  //printf("nodo %d:",i);
	
		my_node[i].kin=aux_in[i].size();
		//printf("se come a %d bichos y ",my_node[i].kin);
		my_node[i].kout=aux_out[i].size();
		//printf("es comido por %d bichos\n",my_node[i].kout);
		my_node[i].in_nodos=new int [my_node[i].kin];
		my_node[i].out_nodos=new int [my_node[i].kout];
		for (int j=0;j<my_node[i].kin;j++){
			my_node[i].in_nodos[j]=aux_in[i][j];
		}
		for (int j=0;j<my_node[i].kout;j++){
			my_node[i].out_nodos[j]=aux_out[i][j];
		}
	}
/*	printf("copiado\n");*/
	
	
	if( nodes.size() > 0) {
		nodes.clear();
		//printf("borro pq es mayor q cero");
	}
	
	for (int i=0;i<N;i++) {
		if (aux_in[i].size() > 0) aux_in[i].clear();
		if (aux_out[i].size() > 0) aux_out[i].clear();
	}
	aux_in.clear();
	aux_out.clear();
	//only in_nodos y out_nodos, pero no in_sig y out_sig
	
	return my_node;
}
//----------------------------------------------------------------------------
nodo * simetrizar(nodo *M,int N)
{	//concretamente pone en una direccion todas las conexiones q habia antes
	vector < vector <int> > aux_in;
	vector < vector <int> > aux_out;
	aux_in.resize(N+1);
	aux_out.resize(N+1);
	
	int repe;
	printf("entramos a simemtrizar N=%d 1\n",N);
	

//printf("simetrizamos\n");

	printf("copiamos los vectores a los auxiliares\n");
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

printf("ponemos los q faltan\n");
	int nodo2;
	for (int i=0;i<N;i++){
		for (int j=0;j<M[i].kin;j++){
			//tenemos que ver si j se come a i
			repe=0;
			nodo2=M[i].in_nodos[j];
			for(int k=0;k<aux_in[nodo2].size();k++){
				if(i==aux_in[nodo2][k]) repe=1;
			}
			if (repe==0) aux_in[nodo2].push_back(i); //El q tiene toda la informacion es aux_in
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
	printf("voy a borrar\n");



	printf("voy a borrar\n");

	nodo *my_node;
	my_node = new nodo [N];
//	return my_node;
	
	
	for (int i=0;i<N;i++){
		my_node[i].kin=aux_in[i].size();
		my_node[i].in_nodos=new int [my_node[i].kin];
		my_node[i].in_sig=new double [my_node[i].kin];
		for (int j=0;j<my_node[i].kin;j++){
			my_node[i].in_nodos[j]=aux_in[i][j];
		}
		my_node[i].kout=aux_out[i].size();
		my_node[i].out_nodos=new int [my_node[i].kout];
		my_node[i].out_sig=new double [my_node[i].kout];
		for (int j=0;j<my_node[i].kout;j++){
			my_node[i].out_nodos[j]=aux_out[i][j];
		}
	}
	
//	printf("matriz simetrica\n");
/*	for (int i=0;i<N;i++){*/
/*	printf("Nodo %d , se come a %d bischos :",i,my_node[i].kin);*/
/*		for (int j=0;j<my_node[i].kin;j++){*/
/*			printf("%d ",my_node[i].in_nodos[j]);*/
/*		*/
/*		}*/
/*		printf("\n");*/
/*	printf("%d es comido por %d bichos:",i,my_node[i].kout);*/
/*		for (int j=0;j<my_node[i].kout;j++){*/
/*			printf("%d ",my_node[i].out_nodos[j]);*/
/*		*/
/*		}*/
/*		printf("\n");*/
/*	}*/
/*	printf("liberamos memoria\n");*/
// liberamos memoria ---------------------------------------------
	for (int i=0;i<N;i++) {
		if (aux_in[i].size() > 0) aux_in[i].clear();
		if (aux_out[i].size() > 0) aux_out[i].clear();
	}
	aux_in.clear();
	aux_out.clear();
	
	//printf("y nos vamos\n");

	return my_node;
}

//---------------------------------------------------------------------------------
nodo * randomizar(nodo *M,int N,gsl_rng *r,int times) //Genera una red randomizada con signo en Mrand
{
	nodo *Mrand;
	//printf("entro a randomizar\n");
	int RELINKS=times*N;
	int nodo1,nodo2;
	int nodo1p,nodo2p;
	int cambios=0;
	//copio todo a Mrand
	
	//Mrand=(nodo *)malloc(N*sizeof(nodo));
	Mrand= new nodo [N];
	
	for (int i=0;i<N;i++){
	//printf("i=%d\n",i);
		Mrand[i].kout=M[i].kout;
		Mrand[i].kin=M[i].kin;
		Mrand[i].label=1;
		Mrand[i].myself=i;
	
		Mrand[i].in_nodos= new int [Mrand[i].kin];
		Mrand[i].in_sig=new double [Mrand[i].kin];
		for (int j=0;j<Mrand[i].kin;j++){
			Mrand[i].in_nodos[j]=M[i].in_nodos[j];
			Mrand[i].in_sig[j]=M[i].in_sig[j];
		}
		Mrand[i].out_nodos =new int [Mrand[i].kout];
		Mrand[i].out_sig =new double [Mrand[i].kout];
		for (int j=0;j<Mrand[i].kout;j++){
			Mrand[i].out_nodos[j]=M[i].out_nodos[j];
			Mrand[i].out_sig[j]=M[i].out_sig[j];
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
   }//END WHILE CAMBIOS < RELINKS
   	//printf("salgo\n");
/*   	printf("matriz randomizada\n");*/
/*	for (int i=0;i<N;i++){*/
/*	printf("Nodo %d , se come a %d bischos :",i,Mrand[i].kin);*/
/*		for (int j=0;j<Mrand[i].kin;j++){*/
/*			printf("%d ",Mrand[i].in_nodos[j]);*/
/*		*/
/*		}*/
/*		printf("\n");*/
/*	printf("%d es comido por %d bichos:",i,Mrand[i].kout);*/
/*		for (int j=0;j<Mrand[i].kout;j++){*/
/*			printf("%d ",Mrand[i].out_nodos[j]);*/
/*		*/
/*		}*/
/*		printf("\n");*/
/*	}*/
	return Mrand;
}

//------------------------------------------------------------------------------------
edge * get_edges(nodo *my_node,int N,int *Nlinks,int contmax)
{
	printf("#entro\n");
	int cont_edges=0;
	
	for (int i=0;i<N;i++) {
		cont_edges += my_node[i].kout;
	}
	printf("#edges=%d\n",cont_edges);
	
	*Nlinks=cont_edges;
	
	edge *my_edge;
//	return my_edge;
	
	//my_edge = new edge [cont_edges];
	printf("#size=%d size=%d\n",sizeof(edge),sizeof(int));
	my_edge = (edge *) malloc((*Nlinks*2)*sizeof(edge));
	if (my_edge == NULL) {printf("No se pudo reservar\n"); }  
	
	printf("#allocated\n");
	
	int edges_number=0;
	
	for (int i=0;i<N;i++){
		for (int j=0;j<my_node[i].kout;j++){
			my_edge[edges_number].in=i;
			my_edge[edges_number].out=my_node[i].out_nodos[j];
			//my_edge[edges_number].loop=new[contmax];
			edges_number += 1;			
		}
	}
	printf("#edges_number=%d\n",edges_number);		
	
	
	for (int i=0;i<*Nlinks;i++) my_edge[i].importance=0;
	
	
	return my_edge;
}

//-------------------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------------------//
int * get_sources(nodo *M,int N,gsl_rng *r,int *Nsources)
{
	int *my_source;
	vector <int> source_aux;
	
	int Nmax,Tmax;//parametros del caminante
	Tmax=1000;//a partir del que empieza a recolecar infomracion
	Nmax=1000;//numero de caminantes que soltamos en las plantas
	int repe,bicho1,bicho2;//bichos por los que va el camino
	int *veces;
	veces=new int[N]; for (int i=0;i<N;i++) veces[i]=0;
	
	//primero vamos a por las evidentes
	for (int i=0;i<N;i++) if ( M[i].kin == 0 ) {
		//source_aux.push_back(i);
		//printf("#M[%d].kin=%d\n",i,M[i].kin);
	}
	
	
	//ahora sacamos el resto con caminantes aleatorios:los debo ir poniendo en toda la red y evolucionar AL CONTRARIO para encontrar las fuentes, podemos intentar encontrar las evidentes
	for (int i=0;i<N;i++){
	//printf("bicho_0=%d\n",i);
		for(int k=0;k<Nmax;k++){//repetimos el lanzamiento Nmax veces
		bicho1=i;
			for (int t=0;t<(int)(2*Tmax);t++){ //Las veces que va a avanzar
				if (M[bicho1].kin == 0) {//si no podemos seguir retrocediendo es una FUENTE. final
					veces[bicho1]=veces[bicho1]+5000;//sumo esto por ejemplo
					bicho1=i;
					break;//empezamos con otro caminante
				}
				else { //si podemos seguri retrocediendo lo hacemos
					int k=gsl_rng_uniform_int(r,M[bicho1].kin);
					bicho2=M[bicho1].in_nodos[k];
					if(t>Tmax){//si ha pasado mas de la mitad del tiempo empezamos a contar
						veces[bicho2] ++;				
					}
					bicho1=bicho2;				
				}//END IF		
			}//END FOR time		
		}//END FOR NMAX walkers
	}//END FOR N nodes

//for (int i=0;i<N;i++) printf("veces[%d]=%d\n",i,veces[i]); //imprime cuantas veces ha apsado el caminante

	//nosquedamos con los que han tenido al caminante mas tiempo, esto quizas haya q generalizarlo
	for (int i=0;i<N;i++) {
		if( (veces[i]!=0) or (M[i].kin == 0)) {
			//printf("spurce!\n");
			source_aux.push_back(i);
		}
	}

	*Nsources=source_aux.size();
	my_source=new int [source_aux.size()];
	for (int i=0;i<source_aux.size();i++) my_source[i]=source_aux[i];
	
	source_aux.clear();
	delete []veces;
	
	return my_source;
}
//-------------------------------------------------------------------------------------------------------------------------//
int * get_drains(nodo *M,int N,gsl_rng *r,int *Ndrains)
{
	int *my_drain;
	vector <int> drain_aux;
	
	int Nmax,Tmax;//parametros del caminante
	Tmax=1000;//a partir del que empieza a recolecar infomracion
	Nmax=1000;//numero de caminantes que soltamos en las plantas
	int repe,bicho1,bicho2;//bichos por los que va el camino
	int *veces;
	veces=new int[N]; for (int i=0;i<N;i++) veces[i]=0;
	
	//ahora sacamos con caminantes aleatorios:los debo ir poniendo en toda la red y evolucionar para encontrar los sumideros, podemos intentar encontrar las evidentes
	for (int i=0;i<N;i++){
	//printf("bicho_0=%d\n",i);
		for(int k=0;k<Nmax;k++){//repetimos el lanzamiento Nmax veces
		bicho1=i;
			for (int t=0;t<(int)(2*Tmax);t++){ //Las veces que va a avanzar
				if (M[bicho1].kout == 0) {//si no podemos seguir retrocediendo es una FUENTE. final
					veces[bicho1]=veces[bicho1]+5000;//sumo esto por ejemplo
					bicho1=i;
					break;//empezamos con otro caminante
				}
				else { //si podemos seguri avanzando lo hacemos
					int k=gsl_rng_uniform_int(r,M[bicho1].kout);
					bicho2=M[bicho1].out_nodos[k];
					if(t>Tmax){//si ha pasado mas de la mitad del tiempo empezamos a contar
						veces[bicho2] ++;				
					}
					bicho1=bicho2;				
				}//END IF		
			}//END FOR time		
		}//END FOR NMAX walkers
	}//END FOR N nodes

for (int i=0;i<N;i++) printf("#veces[%d]=%d\n",i,veces[i]); //imprime cuantas veces ha apsado el caminante

	for (int i=0;i<N;i++) {
		if( ( veces[i]!=0 ) or ( M[i].kout==0 ) ) {
			//printf("drain!\n");
			drain_aux.push_back(i);
		}
	}
	
	for (int i=0;i<N;i++) if (M[i].kout==0) printf("#M[%d]=%d\n",i,M[i].kout);
	
	*Ndrains=drain_aux.size();
	my_drain=new int [drain_aux.size()];
	for (int i=0;i<drain_aux.size();i++) my_drain[i]=drain_aux[i];
	
	delete [] veces;
	drain_aux.clear();
	return my_drain;
}
//-----------------------------------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------------//
int ** crear_matriz_new(int N,nodo * my_node)//meto un vector de nodos y me da una matriz NxN
{
	int ** my_matrix;
	my_matrix=(int **)malloc(N*sizeof(int*));
	for (int i=0;i<N;i++) my_matrix[i]=(int*)calloc(N,sizeof(int));
	
	for (int i=0;i<N;i++){
		for(int j=0;j<my_node[i].kin;j++){
			my_matrix[i][my_node[i].out_nodos[j]]=1;
		}
	}
	
	
	return my_matrix;
}

void liberar_matriz(int N, int** my_matrix)
{
	for (int i=0;i<N;i++) free(my_matrix[i]);
	free(my_matrix);
	
	return;
}//---------------------------------------------------------------------------------------------------//

loop_stat set_loop_stat(int filas,int columnas)
{
	//EL numero de columans da la longitud de los loops q medimos
	//El numero de filas es el binneado en gamma
	loop_stat my_loops;
	
	my_loops.nbins=filas;//(numero de filas, binning de gamma)
	my_loops.contmax=columnas; //(numero de columnas, profunddad del loop)
	
	//genero nbins filas de Lp
	my_loops.redes = (int*)calloc(filas,sizeof(int));
	my_loops.Lp = (double**)malloc(filas*sizeof(double*));
	my_loops.L2p = (double**)malloc(filas*sizeof(double*));
	my_loops.sigma_Lp = (double**)malloc(filas*sizeof(double*));
	for (int i=0;i<filas;i++){
		my_loops.Lp[i]=(double *)calloc(columnas,sizeof(double));
		my_loops.L2p[i]=(double *)calloc(columnas,sizeof(double));
		my_loops.sigma_Lp[i]= (double*)calloc(columnas,sizeof(double*));
	}	
	return my_loops;
}

void liberar_loop_stat(loop_stat my_L)
{
	printf("entro a liberar nbins=%d\n",my_L.nbins);
/*	for (int i=0;i<my_L.nbins;i++){*/
/*		for(int j=0;j<my_L.contmax;j++){*/
/*			printf("L[%d][%d]=%lf\n",i,j,my_L.Lp[i][j]);*/
/*		}*/
/*		printf("\n");*/
/*	}*/
	for (int i=0;i<my_L.nbins;i++){
	//printf("libero %d\n",i);
		free(my_L.Lp[i]);
		free(my_L.L2p[i]);
		free(my_L.sigma_Lp[i]);
	}
	
	free (my_L.Lp);
	free (my_L.L2p);
	free (my_L.sigma_Lp);
	
	return;
}

void liberar_nodos(nodo *my_node,int N)
{
	for (int i=0;i<N;i++){
		delete [] my_node[i].in_nodos;
		delete [] my_node[i].out_nodos;
	//	delete [] my_node[i].in_sig;
	//	delete [] my_node[i].out_sig;
	}
	
	free(my_node);
	
	return;
}

void liberar_edges(edge *my_edge)
{
	free (my_edge);
	return;
}

void liberar_int(int *vector)
{
	delete [] vector;
	return;
}

