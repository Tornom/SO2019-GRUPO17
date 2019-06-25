#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gmp.h> 

void f_borwein(int num_repet);


int main (){
    int num_interacoes = 100000;
    f_borwein(num_interacoes);
    return 0;
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