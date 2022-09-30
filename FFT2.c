#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>


// FFT_2(complex *y, complex *x, int N)： FFT2 by bit reverse
// FFT_3(complex *y, complex *x, int N)： FFT3 by bit reverse
// FFT_5(complex *y, complex *x, int N)： FFT5 by bit reverse
// FFT2(complex *y, complex *x, int N)： FFT with base 2 
// FFT3(complex *y, complex *x, int N)： FFT with base 3 
// FFT5(complex *y, complex *x, int N)： FFT with base 5 
// FFTp(complex *y, complex *x, int N)： FFT with 2^q * 3^q * 5^r

// FFT_p(complex *y, complex *x, int N) : FFT with 2^q * 3^q * 5^r , and others are FT
//										  y: put the result after FFT
//										  x: the input of FFT
//										  N: the size of input array

void FFT_2(complex *x, int N){
	int i, j, k, p, M, L, temp;
	complex u, w, t;
	
	M = N/2;
	j = 0;
	for(i = 1; i < N-1; i++){
		L = M;
		while(j >= L){
			j -= L;
			L /= 2;
		}
		j += L;
		if(j > i){
//			printf("%d <-> %d\n",i,j);
			t    = x[i];
			x[i] = x[j];
			x[j] = t;
		}
	}
	j = 1;
	while(j < N){
		w = 1.0;
		u = cexp(-I*M_PI/j);
		for(i = 0; i < j; i++){
			for(k = 0; k < N; k += 2*j){
				x[k+i+j] *= w;
				t = x[k + i];
				x[k + i] = x[k + i] + x[k + i + j];
				x[k+i+j] = t        - x[k + i + j];
			}
			w *= u;
		}
		j = j << 1;
	}
	
}

void FFT_3(complex *x, int N){
	int i, j, k, p, M, L, temp;
	complex u, w, t1, t2, t3, w1, w2;
	
	M = N/3;
	j = 0;	
	for(i = 1; i < N-1; i++){
		L = M;
		while(j >= 2*L){
			j -= 2*L;
			L /= 3;
		}
		j += L;
		if(j > i){
			printf("%d <-> %d\n",i,j);
			t1   = x[i];
			x[i] = x[j];
			x[j] = t1;
		}
	}
	
	j = 1;
	w1 = cexp(-I*M_PI*2/3); // w4 == w1
	w2 = cexp(-I*M_PI*4/3);
	
	while(j < N){
		w = 1.0;
		u = cexp(-I*2*M_PI/N);
		for(i = 0; i < j; i++){
			for(k = 0; k < N; k += 3*j){
				
				t1 = x[k+i];
				t2 = x[k+i+j]*w;
				t3 = x[k+i+2*j]*cpow(w, 2);
				
				x[k+i]     = t1 +    t2 +    t3;
				x[k+i+j]   = t1 + w1*t2 + w2*t3;
				x[k+i+2*j] = t1 + w2*t2 + w1*t3; 
			}
			w *= u;
		}
		j = j*3;
	}
	
}

void FFT_5(complex *x, int N){
	int i, j, k, p, M, L, temp;
	complex u, w, t1, t2, t3, t4, t5, w1, w2, w3, w4;
	
	M = N/5;
	j = 0;	
	for(i = 1; i < N-1; i++){
		L = M;
		while(j >= 4*L){
			j -= 4*L;
			L /= 5;
		}
		j += L;
		if(j > i){
//			printf("%d <-> %d\n",i,j);
			t1   = x[i];
			x[i] = x[j];
			x[j] = t1;
		}
	}
	
	j = 1;
	w1 = cexp(-I*M_PI*2/5); //	w16 == w6 == w1
	w2 = cexp(-I*M_PI*4/5); //  w12 == w2
	w3 = cexp(-I*M_PI*6/5); //	w8 == w3
	w4 = cexp(-I*M_PI*8/5); //	w9 == w4
	
	while(j < N){
		w = 1.0;
		u = cexp(-I*2*M_PI/N);
		for(i = 0; i < j; i++){
			for(k = 0; k < N; k += 5*j){
				
				t1 = x[k+i];
				t2 = x[k+i+j]*w;
				t3 = x[k+i+2*j]*cpow(w, 2);
				t4 = x[k+i+3*j]*cpow(w, 3);
				t5 = x[k+i+4*j]*cpow(w, 4);
				
				x[k+i]     = t1 +    t2 +    t3 +    t4 +    t5;
				x[k+i+j]   = t1 + w1*t2 + w2*t3 + w3*t4 + w4*t5;
				x[k+i+2*j] = t1 + w2*t2 + w4*t3 + w1*t4 + w3*t5; 
				x[k+i+3*j] = t1 + w3*t2 + w1*t3 + w4*t4 + w2*t5;
				x[k+i+4*j] = t1 + w4*t2 + w3*t3 + w2*t4 + w1*t5;
			}
			w *= u;
		}
		j = j*5;
	}
	
}

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

void FFT3(complex *y, complex *x, int N){
	complex *Nx, *Ny, w, wn, w1, w2;
	complex wn1, wn2;
	int i;
	Nx = (complex*) malloc(N*sizeof(complex));
	Ny = (complex*) malloc(N*sizeof(complex));
	
	for(i = 0; i < N/3; i++){
		Nx[i]         = x[3*i];
		Nx[i + N/3]   = x[3*i + 1];
		Nx[i + 2*N/3] = x[3*i + 2];
	}
	
	w1 = cexp(-I*M_PI*2/3); // w4 == w1
	w2 = cexp(-I*M_PI*4/3);
	if(N == 3){
		
		y[0] = x[0] +    x[1] +    x[2];
		y[1] = x[0] + w1*x[1] + w2*x[2];
		y[2] = x[0] + w2*x[1] + w1*x[2];
		
	}
	else{
		FFT3(Ny, Nx, N/3);
		FFT3(Ny + N/3, Nx + N/3, N/3);
		FFT3(Ny + 2*N/3, Nx + 2*N/3, N/3);
		
		wn = 1; w = cexp(-I*2*M_PI/N);
		for(i = 0; i < N/3; i++){
			
			wn1 = wn*Ny[i + N/3];
			wn2 = cpow(wn, 2)*Ny[i + 2*N/3];
			
			y[i]         = Ny[i] +    wn1 +    wn2;
			y[i + N/3]   = Ny[i] + w1*wn1 + w2*wn2;
			y[i + 2*N/3] = Ny[i] + w2*wn1 + w1*wn2;
			wn *= w;
		}
	}
	free(Nx); free(Ny);
	
	
}

void FFT5(complex *y, complex *x, int N){
	complex *Nx, *Ny, w, wn, w1, w2, w3, w4;
	complex wn1, wn2, wn3, wn4;
	int i;
	
	Nx = (complex*) malloc(N*sizeof(complex));
	Ny = (complex*) malloc(N*sizeof(complex));
	
	for(i = 0; i < N/5; i++){
		Nx[i]         = x[5*i];
		Nx[i + N/5]   = x[5*i + 1];
		Nx[i + 2*N/5] = x[5*i + 2];
		Nx[i + 3*N/5] = x[5*i + 3];
		Nx[i + 4*N/5] = x[5*i + 4];
	}
	
	w1 = cexp(-I*M_PI*2/5); //	w16 == w6 == w1
	w2 = cexp(-I*M_PI*4/5); //  w12 == w2
	w3 = cexp(-I*M_PI*6/5); //	w8 == w3
	w4 = cexp(-I*M_PI*8/5); //	w9 == w4
	
	if(N == 5){
		
		y[0] = x[0] +    x[1] +    x[2] +    x[3] +    x[4];
		y[1] = x[0] + w1*x[1] + w2*x[2] + w3*x[3] + w4*x[4];
		y[2] = x[0] + w2*x[1] + w4*x[2] + w1*x[3] + w3*x[4];
		y[3] = x[0] + w3*x[1] + w1*x[2] + w4*x[3] + w2*x[4];
		y[4] = x[0] + w4*x[1] + w3*x[2] + w2*x[3] + w1*x[4];
		
	}
	else{
		FFT5(Ny, Nx, N/5);
		FFT5(Ny + N/5, Nx + N/5, N/5);
		FFT5(Ny + 2*N/5, Nx + 2*N/5, N/5);
		FFT5(Ny + 3*N/5, Nx + 3*N/5, N/5);
		FFT5(Ny + 4*N/5, Nx + 4*N/5, N/5);
		
		wn = 1; w = cexp(-I*2*M_PI/N);
		
		
		for(i = 0; i < N/5; i++){
			wn1 = wn*Ny[i + N/5];
			wn2 = cpow(wn, 2)*Ny[i + 2*N/5];
			wn3 = cpow(wn, 3)*Ny[i + 3*N/5];
			wn4 = cpow(wn, 4)*Ny[i + 4*N/5];
			
			y[i]         = Ny[i] +    wn1 +    wn2 +    wn3 +    wn4;
			y[i + N/5]   = Ny[i] + w1*wn1 + w2*wn2 + w3*wn3 + w4*wn4;
			y[i + 2*N/5] = Ny[i] + w2*wn1 + w4*wn2 + w1*wn3 + w3*wn4;
			y[i + 3*N/5] = Ny[i] + w3*wn1 + w1*wn2 + w4*wn3 + w2*wn4;
			y[i + 4*N/5] = Ny[i] + w4*wn1 + w3*wn2 + w2*wn3 + w1*wn4;
			wn *= w;
		}
	}
	free(Nx); free(Ny);
	
	
}


void FFTp(complex *y, complex *x, int N){
	complex *Nx, *Ny;
	int i;
	complex w, w1, w2, w3, w4, wn, wn1, wn2, wn3, wn4;
	
	Nx = (complex*) malloc(N*sizeof(complex));
	Ny = (complex*) malloc(N*sizeof(complex));
	
	if(N % 2 == 0){
		for(i = 0; i < N/2; i++){
			Nx[i]       = x[2*i];
			Nx[i + N/2] = x[2*i + 1];
		}
		
		
		if(N == 2){
			y[0] = x[0] + x[1];
			y[1] = x[0] - x[1];
		}
		else{
			FFTp(Ny, Nx, N/2);
			FFTp(Ny + N/2, Nx + N/2, N/2);
			
			wn = 1; w = cexp(-I*2*M_PI/N);
			
			for(i = 0; i < N/2; i++){
				y[i]       = Ny[i] + wn*Ny[i + N/2];
				y[i + N/2] = Ny[i] - wn*Ny[i + N/2];
				wn *= w;
			}
		}
	}
	if(N % 3 == 0){
		for(i = 0; i < N/3; i++){
			Nx[i]         = x[3*i];
			Nx[i + N/3]   = x[3*i + 1];
			Nx[i + 2*N/3] = x[3*i + 2];
		}
		
		w1 = cexp(-I*M_PI*2/3); // w4 == w1
		w2 = cexp(-I*M_PI*4/3);
		if(N == 3){
			
			y[0] = x[0] +    x[1] +    x[2];
			y[1] = x[0] + w1*x[1] + w2*x[2];
			y[2] = x[0] + w2*x[1] + w1*x[2];
			
		}
		else{
			FFTp(Ny, Nx, N/3);
			FFTp(Ny + N/3, Nx + N/3, N/3);
			FFTp(Ny + 2*N/3, Nx + 2*N/3, N/3);
			
			wn = 1; w = cexp(-I*2*M_PI/N);
			for(i = 0; i < N/3; i++){
				wn1 = wn*Ny[i + N/3];
				wn2 = cpow(wn, 2)*Ny[i + 2*N/3];
				y[i]         = Ny[i] +    wn1 +    wn2;
				y[i + N/3]   = Ny[i] + w1*wn1 + w2*wn2;
				y[i + 2*N/3] = Ny[i] + w2*wn1 + w1*wn2;
				wn *= w;
			}
		}
	}
	if(N % 5 == 0){
		for(i = 0; i < N/5; i++){
			Nx[i]         = x[5*i];
			Nx[i + N/5]   = x[5*i + 1];
			Nx[i + 2*N/5] = x[5*i + 2];
			Nx[i + 3*N/5] = x[5*i + 3];
			Nx[i + 4*N/5] = x[5*i + 4];
		}
		
		w1 = cexp(-I*M_PI*2/5); //	w16 == w6 == w1
		w2 = cexp(-I*M_PI*4/5); //  w12 == w2
		w3 = cexp(-I*M_PI*6/5); //	w8 == w3
		w4 = cexp(-I*M_PI*8/5); //	w9 == w4
		
		if(N == 5){
			
			y[0] = x[0] +    x[1] +    x[2] +    x[3] +    x[4];
			y[1] = x[0] + w1*x[1] + w2*x[2] + w3*x[3] + w4*x[4];
			y[2] = x[0] + w2*x[1] + w4*x[2] + w1*x[3] + w3*x[4];
			y[3] = x[0] + w3*x[1] + w1*x[2] + w4*x[3] + w2*x[4];
			y[4] = x[0] + w4*x[1] + w3*x[2] + w2*x[3] + w1*x[4];
			
		}
		else{
			FFTp(Ny, Nx, N/5);
			FFTp(Ny + N/5, Nx + N/5, N/5);
			FFTp(Ny + 2*N/5, Nx + 2*N/5, N/5);
			FFTp(Ny + 3*N/5, Nx + 3*N/5, N/5);
			FFTp(Ny + 4*N/5, Nx + 4*N/5, N/5);
			
			wn = 1; w = cexp(-I*2*M_PI/N);
			
			
			for(i = 0; i < N/5; i++){
				wn1 = wn*Ny[i + N/5];
				wn2 = cpow(wn, 2)*Ny[i + 2*N/5];
				wn3 = cpow(wn, 3)*Ny[i + 3*N/5];
				wn4 = cpow(wn, 4)*Ny[i + 4*N/5];
				
				y[i]         = Ny[i] +    wn1 +    wn2 +    wn3 +    wn4;
				y[i + N/5]   = Ny[i] + w1*wn1 + w2*wn2 + w3*wn3 + w4*wn4;
				y[i + 2*N/5] = Ny[i] + w2*wn1 + w4*wn2 + w1*wn3 + w3*wn4;
				y[i + 3*N/5] = Ny[i] + w3*wn1 + w1*wn2 + w4*wn3 + w2*wn4;
				y[i + 4*N/5] = Ny[i] + w4*wn1 + w3*wn2 + w2*wn3 + w1*wn4;
				wn *= w;
			}
		}
	}
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



int main(){
	complex *x, *y;
	int i , N, temp;
	time_t T;
//	N = 30;
	N = pow(2, 26)*pow(3, 0)*pow(5, 0);
	x = (complex*) malloc(N*sizeof(complex));
	y = (complex*) malloc(N*sizeof(complex));
	
	printf("N = %d\n", N);
	for(i = 0; i < N; i++) x[i] = i;
	
	T = clock();
	FFT2(y, x, N);
	T = clock() - T;
	printf("FFT2 time cost ：%d ms\n", T);
	
	
	for(i = 0; i < N; i++) x[i] = i;
	T = clock();
	FFT_2(x, N);
	T = clock() - T;
	printf("FFT_2 time cost ：%d ms\n", T);
	
//	printf("after FFT\n");
//	for(i = 0; i < N; i++) printf("%f + %f i\n", creal(y[i]), cimag(y[i]));
		
	
	return 0;
} 

