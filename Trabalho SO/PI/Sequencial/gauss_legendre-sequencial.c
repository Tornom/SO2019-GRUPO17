#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gmp.h> 


void f_gauss_legendre(int num_repet);

int main (){
    int num_interacoes = 100000;
    f_gauss_legendre(num_interacoes);
    return 0;
}


void f_gauss_legendre(int num_repet){
	int i;

    mpf_set_default_prec(pow(10,5));

    mpf_t pi; // resultado
	
	mpf_t auxiliar,auxiliar1,aux_somas,aux_mults,t_x4; // variaveis auxiliares
    
	mpf_t a,b,t,p,proximo_a,proximo_b,proximo_t; // variaveis de iteracao
	
	mpf_init(proximo_a);  //inicia as variaveis de big numbers
    
    mpf_init(proximo_b);
    
    mpf_init(proximo_t);
	
    mpf_init(pi);
    
    mpf_init(auxiliar);
    
    mpf_init(auxiliar1);
    
    mpf_init(aux_somas);
    
    mpf_init(aux_mults);
    
    mpf_init(t_x4);
    
    mpf_init_set_d(a,1);
	
	mpf_init_set_d(p,1);
    
    mpf_init_set_d(b,1/sqrt(2));
    
    mpf_init_set_d(t,1/4.0);
	
    
    for(i = 0; i < num_repet; i++){ //comeca a fazer as iteracoes do algoritmo de gauss num_repet vezes seguindo a ordem an+1, bn+1, tn+1 e pn+1
        mpf_add(auxiliar, a, b); 
        mpf_div_ui(proximo_a,auxiliar,2); //an+1

        mpf_mul(aux_mults, a, b);
        mpf_sqrt(proximo_b, aux_mults); //bn+1

        mpf_sub(auxiliar1, a, proximo_a); 
        mpf_mul(aux_mults, auxiliar1, auxiliar1);
        mpf_mul(aux_mults, p, aux_mults);
        mpf_sub(proximo_t, t, aux_mults); //tn+1

        mpf_mul_ui(p, p, 2); //pn+1
        
        
		/*
		printf("||%d||",i);
		gmp_printf("||aux_mults %.10Ff||", aux_mult);
		*/
		
        mpf_set(t,proximo_t);// atualiza os proximos varoles para os valores atuais
		mpf_set(a,proximo_a); 
        mpf_set(b,proximo_b);
    }
    
    mpf_set_ui(aux_mults,0);

    mpf_add(aux_somas, a, b);

    mpf_mul(aux_mults, aux_somas, aux_somas); //(a+b)^2

    mpf_mul_ui(t_x4,t, 4);// 4tn+1

    mpf_div(pi, aux_mults,t_x4); // dividir por 4tn+1
    gmp_printf("\n||Resultado => Gauss-Legendre= %.6Ff      ||\n",pi);
    
    mpf_clear(pi);
	mpf_clear(aux_somas);
    mpf_clear(aux_mults);
    mpf_clear(t_x4);
    mpf_clear(proximo_a);
    mpf_clear(proximo_b);
    mpf_clear(proximo_t);
    mpf_clear(auxiliar);
    mpf_clear(auxiliar1);
}