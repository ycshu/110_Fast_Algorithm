#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <string.h>


void FFT2(complex *y, complex *x, int N){
	complex *Nx, *Ny;
	int i;
	complex wn, w;
	
	Nx = (complex*) malloc(N*sizeof(complex));
	Ny = (complex*) malloc(N*sizeof(complex));
	
	for(i = 0; i < N/2; i++){
		Nx[i]       = x[2*i];
		Nx[i + N/2] = x[2*i + 1];
	}
	
	
	if(N == 2){
		y[0] = x[0] + x[1];
		y[1] = x[0] - x[1];
	}
	else{
		FFT2(Ny, Nx, N/2);
		FFT2(Ny + N/2, Nx + N/2, N/2);
		
		wn = 1; w = cexp(-I*2*M_PI/N);
		
		for(i = 0; i < N/2; i++){
			y[i]       = Ny[i] + wn*Ny[i + N/2];
			y[i + N/2] = Ny[i] - wn*Ny[i + N/2];
			wn *= w;
		}
	}
	free(Nx); free(Ny);
}

void iFFT2(complex *y, complex *x, int N){
	complex *Nx, *Ny;
	int i;
	complex wn, w;
	
	Nx = (complex*) malloc(N*sizeof(complex));
	Ny = (complex*) malloc(N*sizeof(complex));
	
	for(i = 0; i < N/2; i++){
		Nx[i]       = x[2*i];
		Nx[i + N/2] = x[2*i + 1];
	}
	
	
	if(N == 2){
		y[0] = x[0] + x[1];
		y[1] = x[0] - x[1];
	}
	else{
		iFFT2(Ny, Nx, N/2);
		iFFT2(Ny + N/2, Nx + N/2, N/2);
		
		wn = 1; w = cexp(I*2*M_PI/N);
		
		for(i = 0; i < N/2; i++){
			y[i]       = Ny[i] + wn*Ny[i + N/2];
			y[i + N/2] = Ny[i] - wn*Ny[i + N/2];
			wn *= w;
		}
	}
	
	free(Nx); free(Ny);
}


void FFT_p(complex *y, complex *x, int N){
	int p;
	if(N%2 == 0) p = 2;
	else if(N%3 == 0) p = 3;
	else if(N%5 == 0) p = 5;
	else p = N;
	
	
	complex *Nx, *Ny;
	int i, j, k, bias, idx;	
	
	bias = N/p;		//bias  = N/p  : 偏移量 
	
	Nx = (complex*) malloc(N*sizeof(complex));
	Ny = (complex*) malloc(N*sizeof(complex));
	
	for(j = 0; j < p; j++){
		for(k = 0; k < bias; k++){
			Nx[k + j*bias] = x[p*k + j];
		}
	}
	
	// F[p][p] 
	
	complex F[p][p];
	for(i = 0; i < p; i++){
		F[0][i] = 1;
		F[i][0] = 1;
	}
	for(i = 1; i < p/2+1; i++){
		for(j = 1; j < p; j++){
			F[i][j] = cexp(-2*M_PI*i*j*I/p);
		}
	}
	for(i = p/2+1; i < p; i++){	//有反對稱
		for(j = 1; j < p; j++){
			F[i][j] = F[p-i][p-j];
		}
	}
	
	
	//recursive
	if(N == p){
		for(i = 0; i < p; i++){
			y[i] = 0;
			for(j = 0; j < p; j++){
				y[i] += F[i][j]*x[j];
			} 
		}
	}
	else{
		
		for(j = 0; j < p; j++){
			FFT_p(Ny + j*bias, Nx + j*bias, bias);
		}
		
		for(i = 0; i < bias; i++){
			for(k = 0; k < p; k++){
				Ny[i + k*bias] *= cexp(-2*M_PI*k*i*I/N);
			}
			
			for(j = 0; j < p; j++){
				idx = i + j*bias;
				y[idx] = 0;
				for(k = 0; k < p; k++){
					y[idx] += F[j][k]*Ny[i + k*bias];
				}
				
			}
		}
		
	}
	
	free(Nx); free(Ny);
} 



void iFFT_p(complex *y, complex *x, int N){
	int p;
	if(N%2 == 0) p = 2;
	else if(N%3 == 0) p = 3;
	else if(N%5 == 0) p = 5;
	else p = N;
	
	
	complex *Nx, *Ny;
	int i, j, k, bias, idx;	
	
	bias = N/p;		//bias  = N/p  : 偏移量 
	
	Nx = (complex*) malloc(N*sizeof(complex));
	Ny = (complex*) malloc(N*sizeof(complex));
	
	for(j = 0; j < p; j++){
		for(k = 0; k < bias; k++){
			Nx[k + j*bias] = x[p*k + j];
		}
	}
	
	// F[p][p] 
	
	complex F[p][p];
	for(i = 0; i < p; i++){
		F[0][i] = 1;
		F[i][0] = 1;
	}
	for(i = 1; i < p/2+1; i++){
		for(j = 1; j < p; j++){
			F[i][j] = cexp(2*M_PI*i*j*I/p);
		}
	}
	for(i = p/2+1; i < p; i++){	//有反對稱
		for(j = 1; j < p; j++){
			F[i][j] = F[p-i][p-j];
		}
	}
	
	
	//recursive
	if(N == p){
		for(i = 0; i < p; i++){
			y[i] = 0;
			for(j = 0; j < p; j++){
				y[i] += F[i][j]*x[j];
			} 
		}
	}
	else{
		
		for(j = 0; j < p; j++){
			iFFT_p(Ny + j*bias, Nx + j*bias, bias);
		}
		
		for(i = 0; i < bias; i++){
			for(k = 0; k < p; k++){
				Ny[i + k*bias] *= cexp(2*M_PI*k*i*I/N);
			}
			
			for(j = 0; j < p; j++){
				idx = i + j*bias;
				y[idx] = 0;
				for(k = 0; k < p; k++){
					y[idx] += F[j][k]*Ny[i + k*bias];
				}
				
			}
		}
		
	}
	
	free(Nx); free(Ny);
} 






int main(){
	complex *x, *y, *u, *v, *ans_t;
	int i , N;
	
	N = 3;

		
	x = (complex*) malloc(2*N*sizeof(complex));
	y = (complex*) malloc(2*N*sizeof(complex));
	u = (complex*) malloc(2*N*sizeof(complex));
	v = (complex*) malloc(2*N*sizeof(complex));
	ans_t = (complex*) malloc(2*N*sizeof(complex));
	
	char temp[N], a;
	printf("enter the first number x : \n");
	fgets(temp, N*2, stdin);
	for(i = 0; i < N; i++){
		a = temp[i];
		x[N-1-i] = atoi(&a);
	} 
	for(i = N; i < 2*N; i++) x[i] = 0;
	
	printf("\n");
	printf("enter the second number y : \n");
	fgets(temp, N*2, stdin);
	for(i = 0; i < N; i++){
		a = temp[i];
		y[N-1-i] = atoi(&a);
	} 
	for(i = N; i < 2*N; i++) y[i] = 0;
	
//	for(i = 0; i < 2*N; i++){
//		printf("%f + %f i\n", creal(x[i]), cimag(x[i]));
//	}
//	printf("\n");
//	
//	for(i = 0; i < 2*N; i++){
//		printf("%f + %f i\n", creal(y[i]), cimag(y[i]));
//	}
	
	//u: x after fft
	//v: y after fft
	
	FFT_p(u, x, 2*N);
	FFT_p(v, y, 2*N);
	
//	for(i = 0; i < 2*N; i++){
//		printf("%f + %f i\n", creal(u[i]), cimag(u[i]));
//	}
//	printf("\n");
//	
//	for(i = 0; i < 2*N; i++){
//		printf("%f + %f i\n", creal(v[i]), cimag(v[i]));
//	}
//	printf("\n");
	
	for(i = 0; i < 2*N; i++){
		u[i]*=v[i];
		ans_t[i] = 0;
	} 
	
//	for(i = 0; i < 2*N; i++){
//		printf("%f + %f i\n", creal(u[i]), cimag(u[i]));
//	}
//	printf("\n");
	
	iFFT_p(ans_t, u, 2*N);
	for(i = 0; i < 2*N; i++) ans_t[i] /= (2*N);
	
	for(i = 0; i < 2*N; i++) printf("%f\n", creal(ans_t[i]));
	
	float ans;
	ans = 0;
	for(i = 0; i < 2*N; i++){
		ans += pow(10, i)*creal(ans_t[i]);
	}
	
	
	printf("ans is %f\n", ans);
	
	
	return 0;
} 

