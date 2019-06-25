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

mpf_t an,bn,an_1,bn_1; //variaveis de iteracao gauss-legendre
mpf_t an_an_1,pn,tn,tn_1; //variaveis de operacao gauss-legendre


int main (){ //chama as tres funcoes paralelizadas de calculo do pi 
    f_gauss_legendre();
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