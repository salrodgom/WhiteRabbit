#include "general.h"
extern int error;
extern int **ciclos;
extern int cont_ciclos;
extern double *L;
extern double *Lp;
extern double *Ln;
extern double *L_estr;

void get_hierarchy(nodo *my_node,int N,int Nsources,int Ndrains)
{
//si el numero de fuentes permanece continuo, esta parte no tenemos que hacerla!
/*	int *sources;*/
/*	sources=get_sources(M,N,r,&Nsources);//vecotr de fuentes randomizado*/
/*	printf("tenemos %d fuentes:",Nsources);*/
//	for (int i=0;i<Nsources;i++) printf("%d ",sources[i]);
	
	//Veamos como afecta lo q hicimos al orden amiento, o en verdad al no cambiar las fuentes no tenemos q cambair la matriz??
	//ordenar_red(my_node,N,sources,Nsources);//pone las Nsourcs con labels €[0,Nsources], y los demas con labels €[Nsources,N] al alimon

	error=0;
	//printf("calculo jerarquia\n");
	calc_hierarchy(my_node,N,Nsources); //sacamos los niveles
	//printf("done\n");

	//delete []sources;

	return;
}

nodo * crop_network(nodo *my_node,int N,edge *my_edge,int Nlinks,long double Nloops_rand,long double Nloops_exp,int contmax)
{
	nodo *crop_node;
	double *importance,*importance2;
	size_t *importance_index;
	vector <int> aux;
	edge *new_edge;
	int Nlinks_eff;
	
	
	
	for (int i=0;i<Nlinks;i++) my_edge[i].label=-1;
	
/*	for(int i=0;i<Nlinks;i++) {*/
/*		printf("Link[%d] (%d->%d).importance=%d\n",i,my_edge[i].in,my_edge[i].out,my_edge[i].importance);*/
/*	}*/
	
	//printf("defino un nuevo vector de links\n");
	for(int i=0;i<Nlinks;i++) if (my_edge[i].importance != 0) aux.push_back(i);
	new_edge=(edge*)malloc(aux.size()*sizeof(edge));
	for(int i=0;i<aux.size();i++) {
		new_edge[i].in = my_edge[aux[i]].in;
		new_edge[i].out = my_edge[aux[i]].out;
		new_edge[i].importance = my_edge[aux[i]].importance;
		new_edge[i].label= aux[i];
	}
	Nlinks_eff=aux.size();
	
/*	for(int i=0;i<Nlinks_eff;i++) {*/
/*		printf("new_links[%d] (%d->%d).importance=%d\n",i,new_edge[i].in,new_edge[i].out,new_edge[i].importance);*/
/*	}*/
	
	importance=(double*)calloc(Nlinks_eff,sizeof(double));
	importance2=(double*)calloc(Nlinks_eff,sizeof(double));
	importance_index=(size_t*)malloc(Nlinks_eff*sizeof(size_t));
	
	for (int i=0;i<Nlinks_eff;i++) importance[i]=importance2[i]=new_edge[i].importance*1.;
	gsl_sort(importance2,1,Nlinks_eff);	
	gsl_sort_index(importance_index,importance,1,Nlinks_eff); //ya tenemos ordenados los edges por importancia, y su antiguo nombre esta en index
	
	//printf("primero establecemos los labels en la lista larga!\n");
	for (int i=0;i<Nlinks_eff;i++) {
	//printf("para el new_link[%d][%d->%d]:\n",i,new_edge[i].in,new_edge[i].out);
		int index=0;
		for (int j=0;j<new_edge[i].in;j++){
		
			index += my_node[j].kout;
	//		printf("+M[%d].kout(%d) ",j,my_node[j].kout);
		}
	//	printf("\n");
		int nodoi=new_edge[i].in;
		for (int j=0;j<my_node[nodoi].kout;j++){
			if (new_edge[i].out==my_node[nodoi].out_nodos[j]) break;
			index ++;
	//		printf("+ 1(%d)",my_node[nodoi].out_nodos[j]);		
		}
	//	printf("\n");
		//printf("en el largo es Links[%d][%d->%d]\n\n",index,my_edge[index].in,my_edge[index].out);
		my_edge[index].label=i;
		new_edge[i].label2=index;
	}

	
	//for (int i=0;i<Nlinks_eff;i++) printf("importance2[%d]=%.0lf label2[%d]=%d i[%d]=%d %d->%d \n",i,importance2[i],i,new_edge[i].label2,i,importance_index[i],new_edge[i].in,new_edge[i].out);//el nodo "inicial(label)" es el nuevo
	//for (int i=0;i<Nlinks_eff;i++) printf("importance[%d]=%.0lf i[%d]=%d  \n",i,importance[i],i,importance_index[i]);//el nodo "inicial(label)" es el nuevo (label_index)
	int cont=Nlinks_eff; printf("Nlinks_eff=%d\n",Nlinks_eff);
	//printf("entramos al bucle\n");
	
	resta(Nloops_exp,Nloops_rand,Nlinks_eff,new_edge,importance_index,my_node,N,Nlinks_eff,contmax);
	
	//printf("salimos del bucle\n");
	
	aux.clear();
	free(new_edge);
	free(importance);
	free(importance2);
	free(importance_index);
	return crop_node;
}

void resta(long double Nloopsfin,long double Nloops,int cont,edge * my_edge,size_t *imp_index,nodo * my_node,int N,int Nedges,int contmax)
{
	double aux;
	int nodo_in,nodo_out;
	aux=Nloops-my_edge[imp_index[cont-1]].importance; //cont da el tamaño del vector de links eff donde se almacenan los q no son cero
	//restamos usando los nuevos, que tienen menos, pero tenemos que tratar con este vector a traves del extendido con los ceros
	//printf("vamosa restar edge[%ld].importance=%d %d->%d\n",imp_index[cont-1],my_edge[imp_index[cont-1]].importance,my_edge[imp_index[cont-1]].in,my_edge[imp_index[cont-1]].out);
	printf("aux=%lf = %Lf - %d //Nloop_exp=%Lf\n",aux,Nloops,my_edge[imp_index[cont-1]].importance,Nloopsfin);

	if((aux>Nloopsfin) and (cont>=1)) {
		//printf("\n");
		Nloops=aux;
		nodo_in=my_edge[imp_index[cont-1]].in;
		nodo_out=my_edge[imp_index[cont-1]].out;
		//printf(" quitando loops donde aparece el link %d -%d // %lf de %Lf\n",nodo_in,nodo_out,aux,Nloopsfin);
		quitar_link(my_node,nodo_in,nodo_out);
		recalcular_links(my_edge,my_node,N,nodo_in,nodo_out,Nedges,contmax);//elimono los links q ya he recortado indirectamente
		//cont --;siempre intenta restar el mayor despues de reordenar	
		//reordeno la lista de links por importancia!
		if (my_edge[imp_index[cont-1]].importance > 0){
			my_edge[imp_index[cont-1]].importance=0;//quito tb la importance del loop q acabo de romper, tengo q haber quitado todos esos loops
			reordenar_importance(my_edge,imp_index,Nedges);//devuelve el vector imp_index reordenado
			resta(Nloopsfin,Nloops,cont,my_edge,imp_index,my_node,N,Nedges,contmax);
			cont=Nedges-1;
			//printf("     cont=%d ",cont);
		} //si la importancia es cero ya jno lo puedo quitar!!
		

		
	}
	else if ((aux<Nloopsfin) and (cont>=1)){
		//printf("demasié! cont=%d de %d\n",cont,Nedges);
		cont --;//si no puede con el mayor, intenta el siguiente
		resta(Nloopsfin,Nloops,cont,my_edge,imp_index,my_node,N,Nedges,contmax);
	}
	else if ((aux == Nloopsfin) or (cont =1)){
		//printf("yata!\n");
	}	
	
	return;
}
void recalcular_links(edge *my_edge,nodo *my_node,int N,int nodo_in,int nodo_out,int Nedges,int contmax)
{
/*	for (int i=0;i<N+1;i++) L[i]=Lp[i]=Ln[i]=0;*/
/*	int depth=0;*/
/*	int *rastro;*/
/*	int *estado;*/
/*	estado=new int[N];*/
/*	for(int j=0;j<N;j++) estado[j]=1;*/
/*	//printf("encontramos loops desde %d \n",nodo_out);*/
/*	rastro=(int *)malloc(sizeof(int)*1000);	*/

/*	//miro los caminos que pasan por ese link, y vo quitando ese numero de loops de my_edge[i].importance!*/
/*	loops_crop(nodo_out,my_node,my_edge,0,nodo_in,rastro,estado,contmax-1,Nedges);//la longitud debe ser una menos, porque 1 edge YAlo tenemos!!*/
/*	*/
/*	delete []estado;*/
/*	free(rastro);*/
}

void reordenar_importance(edge *my_edge,size_t *imp_index,int Nedges)//devuelve el vector de indices de importancia cambiado
{
	double *importance;
	importance=(double*)calloc(Nedges,sizeof(double));
	for (int i=0;i<Nedges;i++) importance[i]=my_edge[i].importance*1.;
	gsl_sort_index(imp_index,importance,1,Nedges); //ya tenemos ordenados los edges por importancia, y su antiguo nombre esta en index
	free(importance);
	return;
}
void quitar_link(nodo *my_node,int in,int out)
{	
	int aux,index,k; 
	k=my_node[in].kout;
	for (int i=0;i<k;i++){
		if (out==my_node[in].out_nodos[i]) {
			index=i;
			aux=my_node[in].out_nodos[k-1];
			my_node[in].out_nodos[k-1]=out;
			my_node[in].out_nodos[index]=aux;
			break;			
		}
	}
	my_node[in].kout --; //He borrado la conexion de la red en este sentido
	k=my_node[out].kin;
	for (int i=0;i<k;i++){
		if(in==my_node[out].in_nodos[i]){
			index=i;
			aux=my_node[out].in_nodos[k-1];
			my_node[out].in_nodos[k-1]=in;
			my_node[out].in_nodos[index]=aux;
			break;
		}
	}
	my_node[out].kin --; //HE borrado la conexion tb en este sentido
	return;
}


void loops_crop(int NODO,nodo *M,edge *Links,int cont,int final,int *rastro,int *estado,int contmax,int Nedges){//saca los loops
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
		rastro[cont]=final;
/*		for(int o=0; o<cont+1; o++){*/
/*			printf("%d ",rastro[o]);*/
/*		}*/
/*		printf("\n");*/
		//cont ++;
		//printf("ahora cont=%d y rastro[%d]=%d rastro[%d]=%d\n",cont,cont,final,cont,rastro[cont-1]);

		for (int o=0; o < cont; o++){//para cada elemento del loop
		int index;
			for (int k=0;k<Nedges;k++){
				if ((rastro[o]==Links[k].in) and (rastro[o+1]==Links[k].out)) {
					index=k;
					break;
				}
			}
			
		#pragma omp atomic
		Links[index].importance --; 
		//printf("Es el link[%d](%d->%d).importance=%d\n",index,Links[index].in,Links[index].out,Links[index].importance);
		}
		
		
	}
	else {
		if ((estado[nodo2] == 1) and (cont < contmax)) {//lo visitamos
			//M[nodo2].label=0;
			//printf("level %d :",cont);
			//signo=signo*M[NODO].out_sig[i];
			loops_crop(nodo2,M,Links,cont,final,rastro,estado,contmax,Nedges);
			//printf("%d*%d    ",M[NODO].out_sig[i],nodo2);
		}
		else {	
			
		}
	}

}

estado[NODO]=1;

return;
}

nodo * dirigir_red(nodo *my_node,int N,double lambda,gsl_rng *r)
{
	//nodo* directed_node;
	double chance;
	int nodo1, nodo2;
	vector <vector <int> > aux_out; //punto de vista de los pollinators, de ellos van a las plantas
	vector <vector <double> > sig_out;
	vector <vector <int> > aux_in;//punto de vista de las plantas, q pollinators vienen a ellas
	vector <vector <double> > sig_in; 
	//Esto asigna las flechas con probabilidad lambda del nodo menor al mayor
	
	aux_in.resize(N);
	aux_out.resize(N);
	sig_in.resize(N);
	sig_out.resize(N);
	
	for (int i=0;i<N;i++){
		for (int j=0;j<my_node[i].kout;j++){
			if (my_node[i].out_nodos[j] > i){
				//determinar si va al derecho o al reves
				chance=gsl_rng_uniform(r);
				if (chance < lambda) { //va de menor a mayor
			   		if( i != my_node[i].out_nodos[j]){
						nodo1=min(i,my_node[i].out_nodos[j]);
						nodo2=max(i,my_node[i].out_nodos[j]);
			   		}
			 	 //printf("%d->%d\n",nodo1,nodo2);
			   	aux_in[nodo2].push_back(nodo1);
			   	aux_out[nodo1].push_back(nodo2);
				}
				else { //va de mayor a menor
			   		if(i!=my_node[i].out_nodos[j]){
						nodo1=max(i,my_node[i].out_nodos[j]);
						nodo2=min(i,my_node[i].out_nodos[j]);
			   		}
			   	//printf("%d->%d\n",nodo1,nodo2);
			   	aux_in[nodo2].push_back(nodo1);
			   	aux_out[nodo1].push_back(nodo2);			    
				}
			}
		}//for j
	}//for i
	
	//ahora pasamos a nuestra red
	//--------------------------------------
	for (int i=0;i<N;i++){
		delete [] my_node[i].in_nodos;
		delete [] my_node[i].out_nodos;
		delete [] my_node[i].in_sig;
		delete [] my_node[i].out_sig;
	}
	//---------------------------------------
	
	for (int i=0;i<N;i++) {
		my_node[i].kin=0; my_node[i].kout=0;
	}
	//printf("lo metemos dentro de la final\n");
for (int i=0;i<N;i++){
	my_node[i].kin=aux_in[i].size();
	my_node[i].in_nodos=new int [my_node[i].kin];
	my_node[i].in_sig=new double [my_node[i].kin];
//printf("Kin[%d]=%d %d \n ",i,aux_in[i].size(),my_node[i].kin);
//printf("M[%d]in=",i);
	for (int j=0;j<my_node[i].kin;j++){
		my_node[i].in_nodos[j]=aux_in[i][j];
		//my_node[i].in_sig[j]=sig_in[i][j];
		my_node[i].in_sig[j]=1;
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
		//my_node[i].out_sig[j]=sig_out[i][j];
		my_node[i].out_sig[j]=1;
//		printf("%d*(%.0f) ",my_node[i].out_nodos[j],my_node[i].out_sig[j]);
	}
//	printf("\n");
}
	
	
	aux_out.clear();
	aux_in.clear();
	sig_out.clear();
	sig_in.clear();

	
	return my_node;
}

int get_net_number(char *nombre)
{
	int number;
	//Foodwebs
	if (strcmp (nombre,"coachella") == 0) number=1;
	else if (strcmp (nombre,"reef") == 0) number=2;
	else if (strcmp (nombre,"shelf") == 0) number=3;
	else if (strcmp (nombre,"el_verde") == 0) number=4;
	else if (strcmp (nombre,"little_rock") == 0) number=5;
	//social networks
	else if (strcmp (nombre,"twitter_ego") == 0) number=6;
	else if (strcmp (nombre,"polblog") == 0) number=7;
	//else if (strcmp (nombre,"adjacency_spanish") == 0) number=7;
	else if (strcmp (nombre,"citation_dblp") == 0) number=8;
	else if (strcmp (nombre,"advocato") == 0) number=8;
	else if (strcmp (nombre,"slashdot_replies") == 0) number=9;
	else if (strcmp (nombre,"digg_replies") == 0) number=10;
	else if (strcmp (nombre,"wikivote") == 0) number=11;
	else if (strcmp (nombre,"ownership") == 0) number=12;
	else if (strcmp (nombre,"consulting") == 0) number=13;
	else if (strcmp (nombre,"kaitiaki") == 0) number=14;
	//technical networks
	//else if (strcmp (nombre,"aeropuertos") == 0) number=17;
	else if (strcmp (nombre,"p2p-Gnutella08") == 0) number=17;
	//biological networks
	else if (strcmp (nombre,"celegans_neural") == 0) number=19;
	else if (strcmp (nombre,"tuberculosis_trans") == 0) number=20;
	else if (strcmp (nombre,"ecoli-iJO1366") == 0) number=21;
	else if (strcmp (nombre,"ecoli-iAF1260") == 0) number=22;
	else if (strcmp (nombre,"saureus-bigg") == 0) number=23;
	else if (strcmp (nombre,"ecoli_bmc") == 0) number=24;
	else if (strcmp (nombre,"subtilis") == 0) number=25;
	else if (strcmp (nombre,"tuberculosis_bmc") == 0) number=26;
	//else if (strcmp (nombre,"") == 0) number=;

	else number = -1;
	return number;

}

double get_double_links_prob(nodo *my_node,int N)
{
	double delta;
	int double_links,Ktot;
	printf("N=%d ",N);
	Ktot=0;
	for (int i=0;i<N;i++) Ktot += my_node[i].kout;
	printf("numero de edges=%d\n",Ktot);

	//sacamos cuantos links dobles hay:
	double_links=0;
	for (int i=0;i<N;i++){
	//printf("nodo %d [%d]:",i,my_node[i].kout);
		for (int j=0;j<my_node[i].kout;j++){
		//printf("%d ",my_node[i].out_nodos[j]);
			int nodo2=my_node[i].out_nodos[j];
			for (int k=0;k<my_node[nodo2].kout;k++){
				if (my_node[nodo2].out_nodos[k]==i) {
					double_links ++;
					//printf("%d->%d \n",i,my_node[i].out_nodos[j]);
				}
			}
		}
		//printf("\n");
	}
	printf("numero de links dobles=%d\n",double_links);
	//delta=(double_links*1.)/(Ktot*1. - double_links*0.5);
	double_links=double_links/2.;
	delta=(double_links*1.)/(Ktot*1.-double_links);
	
	printf("proporcion de loops dobles=%lf\n",delta);

	return delta;
}

double calc_coherence(nodo *my_node,int N,int Nsources)
{
	double coherence,q,q2,q_i,sum;
	
	
	q=0.;
sum=0.;
int linko=0;
for (int i=0; i<N; i++)
{	
	//printf("nodo %d\n",i);	
	q_i=0.;
	//printf("%f \n",my_node[i].s);
	for (int j=0;j<my_node[i].kin;j++)
	{	
		//printf("(+ %f) =%f\n",fabs(my_node[i].s - my_node[my_node[i].in_nodos[j]].s - 1.),q_i);
		q_i = q_i + fabs(my_node[i].s - my_node[my_node[i].in_nodos[j]].s - 1.);
		//printf("[%d] |%f - %f -1 |=%f \n",my_node[i].in_nodos[j],my_node[i].s,my_node[my_node[i].in_nodos[j]].s,fabs(my_node[i].s - my_node[my_node[i].in_nodos[j]].s - 1.));
		sum+=fabs(my_node[i].s - my_node[my_node[i].in_nodos[j]].s - 1.);
		linko++;			
	}
	
	q_i=max(0.,((q_i)/(my_node[i].kin*1.)));
	
	q2 += q_i;
	//printf(" q[%d]=%f\n",i,q_i);	
}
q2=q2/(N - Nsources)*1.;
q=sum/(1.*linko);
//printf("q=%f q2=%f\n",q,q2);
coherence=q2;	

	return coherence;
}


