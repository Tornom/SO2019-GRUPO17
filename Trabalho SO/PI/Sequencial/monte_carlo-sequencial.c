#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gmp.h> 

void f_monte_carlo(int num_repet);


int main (){
	srand(time(NULL));
    int num_interacoes = 100000;
   	f_monte_carlo(num_interacoes);
    return 0;
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
