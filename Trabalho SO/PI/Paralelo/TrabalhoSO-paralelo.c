#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gmp.h> 
#include <pthread.h>

#define num_repet 100000



void f_gauss_legendre(); // funcoes para paralelizar gauss-legendre

void *calcular_elemento_a();

void *calcular_elemento_b();

void f_monte_carlo();

void *calculo_monte_carlo(void *thread_atual);// funcao paralelizada monte carlo

void f_borwein();

void *fator_1();

void *fator_4();

void *fator_5();

void *fator_6();


mpf_t an,bn,an_1,bn_1; //variaveis de iteracao gauss-legendre
mpf_t an_an_1,pn,tn,tn_1; //variaveis de operacao gauss-legendre

mpf_t num_ptos_cir[4], pi, pontos_dentro_do_circulo;//variaveis de operacao monte carlo

mpf_t pow_aux,var_fator_1,var_fator_4,var_fator_5,var_fator_6,soma_aux,var_k8,pi, aux; //variaveis de execucao do algoritmo borwein
mpf_t var_1, var_2, var_4, var_8, var_16; // constantes do algoritmo borwein

int main (){ //chama as tres funcoes paralelizadas de calculo do pi 
	srand(time(NULL));
    f_gauss_legendre();
	f_borwein();
	f_monte_carlo();
    return 0;
}

void *calcular_elemento_a(){
	
    mpf_add(an_1, an, bn);
    mpf_div_ui(an_1,an_1,2); //(an+bn) /2

	
    mpf_sub(an_an_1, an, an_1);
	
    mpf_pow_ui(an_an_1, an_an_1, 2); //salva o valor de (an - (an+1))^2 para usar posteriormente para o calculo de tn+1
    
    pthread_exit(0); // funcao para sair da thread com o status de 0 (terminou a execucao corretamente)
}

void *calcular_elemento_b(){
	
    mpf_mul(bn_1, an, bn);
	
    mpf_sqrt(bn_1, bn_1); // (an*bn)^(1/2)
    
    pthread_exit(0); // funcao para sair da thread com o status de 0 (terminou a execucao corretamente)
}

void f_gauss_legendre(){
	int i;
	
	pthread_t vetor_threads_a[num_repet]; // cria os vetores onde serao armazenados os ids das threads
    pthread_t vetor_threads_b[num_repet];
	
    mpf_set_default_prec(pow(10,5));

    mpf_t pi; //resultado
	mpf_t aux_mult, aux_soma;

    mpf_init(pi); // inicia as variaveis para big numbers
	
    mpf_init(aux_mult);
	
	mpf_init(aux_soma);
	
	mpf_init(an_1);
	
    mpf_init(bn_1);
	
    mpf_init(tn_1);
	
    mpf_init(an_an_1);

    mpf_init_set_d(an,1); // guarda os valores iniciais para cada um dos iteradores
	
	mpf_init_set_d(pn,1);
    
    mpf_init_set_d(bn,1/sqrt(2));
    
    mpf_init_set_d(tn,1/4.0);
	
	//printf("chegou");
    
    
    for(i = 0; i < num_repet; i++){
        pthread_create (&vetor_threads_a[i], NULL, calcular_elemento_a, &i); //cria as threads no vetor de threads
		
        pthread_create (&vetor_threads_b[i], NULL, calcular_elemento_b, &i);
		
		

        pthread_join(vetor_threads_a[i],NULL); // espera as threads acabarem sua execucao e retorna para a parte sequencial
		
        pthread_join(vetor_threads_b[i],NULL);

        mpf_mul(aux_mult, pn, an_an_1); // tn-pn*(an-an+1)^2
		/*
			gmp_printf("||mult %.15Ff||", aux_mult);
			gmp_printf("||pn  %.15Ff||", aux_mult);
			
		*/
        mpf_sub(tn_1, tn, aux_mult); 

        mpf_mul_ui(pn, pn, 2); // atualiza os proximos varoles para os valores atuais
        mpf_set(an,an_1);
        mpf_set(bn,bn_1);
        mpf_set(tn,tn_1);
    }

    mpf_add(aux_soma, an, bn);

    mpf_mul(aux_mult,aux_soma,aux_soma); //aux_soma ^2

    mpf_mul_ui(tn,tn, 4);//tn*4

    mpf_div(pi, aux_mult,tn);
	
	
    gmp_printf("\n||Resultado => Gauss-Legendre= %.6Ff      ||\n", pi);

    mpf_clear(pi);
    
	
    mpf_clear(an_an_1);
    mpf_clear(aux_mult);
	mpf_clear(aux_soma);
	mpf_clear(an_1);
    mpf_clear(bn_1);
    mpf_clear(tn_1);
	mpf_clear(an);
    mpf_clear(bn);
    mpf_clear(tn);
    mpf_clear(pn);
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

/*
	quebrar as quatro partes da operacao de borwein para a paralelizacao
*/

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
