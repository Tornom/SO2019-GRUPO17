#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>

int f_black_scholes();

void *calculo_black_scholes(void* qual_int);

double S, E, r, sig, T;

int num_int;

double var_mult_t[4],sig_t[4],operacao[4],aux[4],valor_t[4],maximo[4];

double testes[4][100000];

int main()
{
    srand(time(NULL));
    f_black_scholes();
}

int f_black_scholes()
{
	int i = 0;
	double soma_media[4], soma_variancia[4];
    double vari[4], conf_min, conf_max;
    double raiz_vari[4], maximo, intervalo[4];
    double media_testes[4];

    for(i = 0; i < 4; i++)
    soma_media[i] = soma_variancia[i] = 0;

    printf("Coloque os valores de entrada S, E, r, Sigma, T e M na seguinte forma:\nS\nE\nr\nSigma\nT\nM\n\nPor fim digite ok e pressione enter\n");


    scanf("%lf",&S);
    scanf("%lf",&E);
    scanf("%lf",&r);
    scanf("%lf",&sig);
    scanf("%lf",&T);
    scanf("%d ",&num_int);

    pthread_t vetor_threads[4];

    
    int var_1 = 0;
    int var_2 = 1;
    int var_3 = 2;
    int var_4 = 3;

    pthread_create (&vetor_threads[var_1], NULL, calculo_black_scholes, &var_1);
    pthread_create (&vetor_threads[var_2], NULL, calculo_black_scholes, &var_2);
    pthread_create (&vetor_threads[var_3], NULL, calculo_black_scholes, &var_3);
    pthread_create (&vetor_threads[var_4], NULL, calculo_black_scholes, &var_4);

 
	for (i = 0; i< 4 ; i++){
		pthread_join (vetor_threads[i], NULL);
    }	

    
	int k;
    for(i = 0; i < 4; i++)
    {
        k = 0;
        while(k < num_int) // Calcula a media dos testes
        {
            soma_media[i] += testes[i][k];
            k++;
        }

        media_testes[i] = soma_media[i] / (double) num_int;

        k = 0;

        while(k < num_int) // Calcula a variancia dos testes
        {
            soma_variancia[i] += pow((testes[i][k] - media_testes[i]),2);
            k++;
        }
        vari[i] = soma_variancia[i] / (double) num_int;

        raiz_vari[i] = sqrt(vari[i]); // Raiz da variancia

        intervalo[i] = raiz_vari[i]*1.96/(double)sqrt(num_int); // Por fim o intervalo minimo e maximo da confianca
    }
    
    int j = 0;
    double med_in = 0;
    double med_med = 0;
    for(j = 0; j< 4;j++)
    {
        med_in += intervalo[j]/4.0;
        med_med += media_testes[j]/4.0;
    }

	conf_min = med_med - med_in;

	conf_max = med_med + med_in;

	printf("|| S . . . . . . . . %0.lf ||\n",S);
	printf("|| E . . . . . . . . %0.lf ||\n",E);
	printf("|| r . . . . . . . . %0.lf ||\n",r);
	printf("|| sig . . . . . . . %0.lf ||\n",sig);
	printf("|| T . . . . . . . . %0.lf ||\n",T);
	printf("|| M . . . . . . . . %d || \n",num_int);
	printf("|| Intervalo de confianca = { Min :%lf, Max :%lf} ||\n\n", conf_min,conf_max);

}

void *calculo_black_scholes(void* qual_int)
{
    
    int i = 0;
    int numero_thread = (*(int*)qual_int);

    for(i = 0; i < num_int; i++)// realiza o algoritmo de black scholes iterando na quantia de ?num_int interacoes
    { 
		sig_t[numero_thread] = sig * sqrt(T); 
        var_mult_t[numero_thread] = T * (r - (1 / 2.0) * pow(sig,2));
		sig_t[numero_thread] = sig_t[numero_thread] * (double) (rand() / (double) RAND_MAX);
		operacao[numero_thread] = var_mult_t[numero_thread] + sig_t[numero_thread];
		valor_t[numero_thread] = S * pow(M_E, operacao[numero_thread]);
		if(valor_t[numero_thread] - E > 0)
        {
        	maximo[numero_thread] = valor_t[numero_thread]-E;
		}
        else
        {
        	maximo[numero_thread] = 0;
		}
		aux[numero_thread] = pow(M_E, - r * T);
		testes[numero_thread][i] = aux[numero_thread] * maximo[numero_thread];
	}   
}


