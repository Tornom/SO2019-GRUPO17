#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void f_black_scholes();

int main(){
    srand(time(NULL));
	f_black_scholes();
	return 0;
}

void f_black_scholes(){
	int i, num_int;
	double var_mult_t,sig_t,operacao,aux,valor_t;
	double soma_media, soma_variancia;
	double vari, conf_min, conf_max;
    double raiz_vari, maximo, intervalo;
    double S,E,r,sig,T;

    soma_media = soma_variancia = 0;




    /* 
    FILE* fp;
    fp = fopen("entrada_blackscholes.txt","rw");
    fscanf(fp,"%ld",&S);
    fscanf(fp,"%ld",&E);
    fscanf(fp,"%ld",&r);
    fscanf(fp,"%ld",&sig);
    fscanf(fp,"%d",&num_int);

    printf("%ld \n %ld \n %ld \n %ld \n %ld \n %d \n",S,E,r,sig,T,num_int);
    fclose(fp);
    */

    printf("Coloque os valores de entrada S, E, r, Sigma, T e M na seguinte forma:\nS\nE\nr\nSigma\nT\nM\n\nPor fim digite ok e pressione enter\n");

	scanf("%lf",&S);
    scanf("%lf",&E);
    scanf("%lf",&r);
    scanf("%lf",&sig);
    scanf("%lf",&T);
    scanf("%d ",&num_int);

	double testes[num_int], media_testes;
	
	for(i = 0; i < num_int; i++)// realiza o algoritmo de black scholes iterando na quantia de ?num_int interacoes
    { 
		sig_t = sig * sqrt(T); 
        var_mult_t = T * (r - (1 / 2.0) * pow(sig,2));
		sig_t = sig_t * (double) (rand()/(double)RAND_MAX);
		operacao = var_mult_t + sig_t;

		valor_t = S * pow(M_E, operacao);
		
		if(valor_t - E > 0)
        {
        	maximo = valor_t-E;
		}
        else
        {
        	maximo = 0;
		}

		aux = pow(M_E, - r * T);

		testes[i] = aux * maximo;

	}

    int k = 0;

    while(k < num_int) // Calcula a media dos testes
    {
        soma_media += testes[k];
        k++;
    }

    media_testes = soma_media / (double) num_int;

    k =0;

	while(k < num_int) // Calcula a variancia dos testes
    {
        soma_variancia += pow((testes[k] - media_testes),2);
        k++;
    }
    vari = soma_variancia / (double) num_int;

    raiz_vari = sqrt(vari); // Raiz da variancia

	intervalo = raiz_vari*1.96/(double)sqrt(num_int); // Por fim o intervalo minimo e maximo da confianca

	conf_min = media_testes - intervalo;

	conf_max = media_testes + intervalo;


	printf("|| S . . . . . . . . %0.lf ||\n",S);
	printf("|| E . . . . . . . . %0.lf ||\n",E);
	printf("|| r . . . . . . . . %0.lf ||\n",r);
	printf("|| sig . . . . . . . %0.lf ||\n",sig);
	printf("|| T . . . . . . . . %0.lf ||\n",T);
	printf("|| M . . . . . . . . %d || \n",num_int);
	printf("|| Intervalo de confianca = { Min :%lf, Max :%lf} ||\n\n", conf_min,conf_max);


}
