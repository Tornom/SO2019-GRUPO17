#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gmp.h> 
#include <pthread.h>

#define num_repet 100000

void f_borwein();

void *fator_1();

void *fator_4();

void *fator_5();

void *fator_6();

mpf_t pow_aux,var_fator_1,var_fator_4,var_fator_5,var_fator_6,soma_aux,var_k8,pi, aux; //variaveis de execucao do algoritmo borwein
mpf_t var_1, var_2, var_4, var_8, var_16; // constantes do algoritmo borwein

int main (){ //chama as tres funcoes paralelizadas de calculo do pi 
	f_borwein();
    return 0;
}


void *fator_1(){
	mpf_add_ui(var_fator_1,var_k8, 1);
	mpf_div(var_fator_1, var_4, var_fator_1);
	pthread_exit(0);
}

void *fator_4(){
	mpf_add_ui(var_fator_4, var_k8, 4);
	mpf_div(var_fator_4, var_2, var_fator_4);
	pthread_exit(0);
}

void *fator_5(){
	mpf_add_ui(var_fator_5, var_k8, 5);
	mpf_div(var_fator_5, var_1, var_fator_5);
	pthread_exit(0);
}

void *fator_6(){
	mpf_add_ui(var_fator_6, var_k8, 6);
	mpf_div(var_fator_6, var_1, var_fator_6);
	pthread_exit(0);
}

void f_borwein(){
	int i;
	
	pthread_t vetor_var_fator_1[num_repet]; // vetores para a geracao das threads divididas entre as 4 operacoes
    pthread_t vetor_var_fator_4[num_repet];
    pthread_t vetor_var_fator_5[num_repet]; 
    pthread_t vetor_var_fator_6[num_repet];  
	
    mpf_set_default_prec(pow(10,5));

	mpf_init(var_k8); //inicia as variaveis de big numbers
	
	mpf_init(var_fator_1);
	
	mpf_init(var_fator_4);
	
	mpf_init(var_fator_5);
	
	mpf_init(var_fator_6);
	
	mpf_init(pi);
	
	mpf_init(aux);
	
	mpf_init(pow_aux);
	
	mpf_init(soma_aux);
    
	
	/*
		constantes para serem utilizadas no algoritmo
	*/
	
	mpf_init_set_ui(var_1,1); //salva valor 1 em big numbers
	
	mpf_init_set_ui(var_2,2); //salva valor 2 em big numbers
	
	mpf_init_set_ui(var_4,4); //salva valor 4 em big numbers
	
	mpf_init_set_ui(var_8,8); //salva valor 8 em big numbers
	
	mpf_init_set_ui(var_16,16); //salva valor 16 em big numbers


	for(i = 0; i < num_repet; i++){
		mpf_mul_ui(var_k8, var_8, i);
		pthread_create (&vetor_var_fator_1[i], NULL, fator_1,&i);
		pthread_create (&vetor_var_fator_4[i], NULL, fator_4,&i);
		pthread_create (&vetor_var_fator_5[i], NULL, fator_5,&i);
		pthread_create (&vetor_var_fator_6[i], NULL, fator_6,&i);
		

		pthread_join (vetor_var_fator_1[i], NULL); //espera todas as threads acabarem para terminar de juntar seus resultados
		pthread_join (vetor_var_fator_4[i], NULL);
		pthread_join (vetor_var_fator_5[i], NULL);
		pthread_join (vetor_var_fator_6[i], NULL);
		
		mpf_pow_ui(pow_aux, var_16, i);
		mpf_div(pow_aux, var_1, pow_aux);

		mpf_sub(soma_aux,var_fator_1,var_fator_4);
		
		mpf_sub(soma_aux,soma_aux,var_fator_5);
		
		mpf_sub(soma_aux,soma_aux,var_fator_6);
		
		mpf_mul(aux,pow_aux,soma_aux);
		
		mpf_add(pi, pi, aux);
	}

    gmp_printf("\n||Resultado => Borwein= %.6Ff      ||\n", pi);
	
	mpf_clear(var_k8);
	mpf_clear(var_fator_1);
	mpf_clear(var_fator_4);
	mpf_clear(var_fator_5);
	mpf_clear(var_fator_6);
	mpf_clear(aux);
	mpf_clear(pow_aux);
	mpf_clear(soma_aux);
}
