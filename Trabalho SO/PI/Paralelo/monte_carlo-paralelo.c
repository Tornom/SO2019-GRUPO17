#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gmp.h> 
#include <pthread.h>

#define num_repet 100000

void f_monte_carlo();

void *calculo_monte_carlo(void *thread_atual);// funcao paralelizada monte carlo

mpf_t num_ptos_cir[4], pi, pontos_dentro_do_circulo;//variaveis de operacao monte carlo

int main (){ //chama as tres funcoes paralelizadas de calculo do pi 
	srand(time(NULL));
	f_monte_carlo();
    return 0;
}


void *calculo_monte_carlo(void *thread_atual){
	int i;
	int qual_thread = *((int *)thread_atual);
	int num_iteracoes_monte_carlo = num_repet/4;
	mpf_set_default_prec(pow(10,5));
	mpf_init(num_ptos_cir[qual_thread]);
	double x, y;
	

	for(i = 0; i < num_iteracoes_monte_carlo; i++){
		x = (double)(rand() / (double)RAND_MAX);
		y = (double)(rand() / (double)RAND_MAX);
		if( ((x*x) + (y*y)) <= 1){
			mpf_add_ui(num_ptos_cir[qual_thread], num_ptos_cir[qual_thread],1);
		}
	}
	pthread_exit(0);
}
	
void f_monte_carlo(){
	int i,j,k;
	
	mpf_set_default_prec(pow(10,5));
	
	pthread_t vetor_threads[4];
    
    mpf_t quantos_pontos_big;
	
    mpf_init(pontos_dentro_do_circulo);
   	mpf_init(pi);
	mpf_init_set_ui(quantos_pontos_big, num_repet);
	
	int* vetor = (int*)malloc(sizeof (int));
	vetor[0] = 0;
	pthread_create (&vetor_threads[0], NULL, calculo_monte_carlo, vetor); // cria quatro threads para rodarem o algoritmo de monte carlo em paralelo
	vetor[0] = 1;
	pthread_create (&vetor_threads[1], NULL, calculo_monte_carlo, vetor);
	vetor[0] = 2;
	pthread_create (&vetor_threads[2], NULL, calculo_monte_carlo, vetor);
	vetor[0] = 3;
	pthread_create (&vetor_threads[3], NULL, calculo_monte_carlo, vetor);
	
	pthread_join (vetor_threads[0], NULL);	//espera todas as threads acabarem para realizar os  calculos
	pthread_join (vetor_threads[1], NULL);
	pthread_join (vetor_threads[2], NULL);
	pthread_join (vetor_threads[3], NULL);

	j = 0;
	while(j < 4){ //soma todos os pontos obtidos na mesma variavel
		mpf_add(pontos_dentro_do_circulo, pontos_dentro_do_circulo, num_ptos_cir[j]);
		j++;
	}
	
	mpf_mul_ui(pi, pontos_dentro_do_circulo, 4);
	mpf_div(pi,pi,quantos_pontos_big);


	free(vetor);
    gmp_printf("\n||Resultado => Monte Carlo= %.6Ff      ||\n\n\n", pi);
    mpf_clear(pi);
    mpf_clear(quantos_pontos_big);
    mpf_clear(pontos_dentro_do_circulo);
}