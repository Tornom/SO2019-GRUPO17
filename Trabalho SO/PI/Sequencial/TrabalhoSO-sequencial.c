#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gmp.h> 

void f_borwein(int num_repet);

void f_gauss_legendre(int num_repet);

void f_monte_carlo(int num_repet);


int main (){
	srand(time(NULL));
    int num_interacoes = 100000;
    f_gauss_legendre(num_interacoes);
    f_borwein(num_interacoes);
   	f_monte_carlo(num_interacoes);
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

void f_borwein(int num_repet){
	int k;
    mpf_set_default_prec(pow(10,5));
    
	mpf_t pi; // resultado
	mpf_t var_k8, var_final, var_x, var_y, var_z,var_16, var_1, var_4, var_5, var_6; //variaveis operacionais
	
	mpf_t valor_1, valor_2, valor_4, valor_8, valor_16; //constantes para operacao
	
	mpf_init(pi); //inicia as variaveis de big numbers
	
	mpf_init(var_k8);
	
	mpf_init(var_final);
	
	mpf_init(var_16);
	
	mpf_init(var_1);
	
	mpf_init(var_4);
	
	mpf_init(var_5);
	
	mpf_init(var_6);
	
	mpf_init(var_x);
	
	mpf_init(var_y);
	
	mpf_init(var_z);

	
	/*
	if(pi == NULL){
		printf("ponteiro deu null");
	}
	*/
	
	
	mpf_init_set_ui(valor_1,1); //salva valor 1 em big numbers
	
	mpf_init_set_ui(valor_2,2); //salva valor 2 em big numbers
	
	mpf_init_set_ui(valor_4,4); //salva valor 4 em big numbers
	
	mpf_init_set_ui(valor_8,8); //salva valor 8 em big numbers
	 
	mpf_init_set_ui(valor_16,16); //salva valor 16 em big numbers

	for(k = 0; k < num_repet; k++){ //repete o somatorio das multiplicacoes das parcelas 16,1,4,5,6

		mpf_pow_ui(var_16, valor_16, k);
		mpf_div(var_16, valor_1, var_16); // parcela 16
		
		mpf_mul_ui(var_k8, valor_8, k);
		mpf_add_ui(var_1,var_k8, 1);
		mpf_div(var_1, valor_4, var_1); // parcela 1

		mpf_add_ui(var_4, var_k8, 4);
		mpf_div(var_4, valor_2, var_4); // parcela 4

		mpf_add_ui(var_5, var_k8, 5);
		mpf_div(var_5, valor_1, var_5); // parcela 5

		mpf_add_ui(var_6, var_k8, 6);
		mpf_div(var_6, valor_1, var_6); // parcela 6

		mpf_sub(var_x,var_1,var_4);
		mpf_sub(var_y,var_x,var_5);
		mpf_sub(var_z,var_y,var_6);
		mpf_mul(var_final,var_16,var_z); // junta todas as parcelas 

		mpf_add(pi, pi, var_final);//somando os resultados obtidos
	}

    gmp_printf("\n||Resultado => Borwein= %.6Ff      ||\n",pi); //mostra o resultado final com 6 casas decimais

    mpf_clear(pi);
	mpf_clear(var_k8);
	mpf_clear(var_final);
	mpf_clear(var_x);
	mpf_clear(var_y);
	mpf_clear(var_z);
	mpf_clear(var_16);
	mpf_clear(var_1);
	mpf_clear(var_4);
	mpf_clear(var_5);
	mpf_clear(var_6);

}

void f_monte_carlo(int num_repet){
	int i;
	double x,y;
	
    mpf_set_default_prec(pow(10,5));
	mpf_t pi; // resultado
	mpf_t pontos_no_circ, pontos_rand;
	
	mpf_init(pontos_rand); //inicia as variaveis de big numbers
	
	mpf_init(pi);
	
	mpf_init(pontos_no_circ);
	
	mpf_init_set_ui(pontos_rand, num_repet);
	
	
	for(i = 0; i < num_repet; i++){
		x = (double)(rand()/(double)RAND_MAX);
		y = (double)(rand()/(double)RAND_MAX);
		
		/*
			printf("valor x = %d, valor y = %d\n",x,y);
		*/
		
		if((x*x )+ (y*y) <= 1){
			mpf_add_ui(pontos_no_circ, pontos_no_circ, 1);
		}
	}
	
	mpf_div(pi, pontos_no_circ, pontos_rand);
	
	mpf_mul_ui(pi,pi,4); // como o valor encontrado e' pi/4 multiplica-se o valor por quatro para encontrar o valor de pi

    gmp_printf("\n||Resultado => Monte Carlo= %.6Ff      ||\n",pi);
    
    
	mpf_clear(pontos_no_circ);
	mpf_clear(pontos_rand);
	mpf_clear(pi);
}
