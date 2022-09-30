#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <string.h>
int FinSize = 0;

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


int* FastMul(char x[], char y[]){
	complex *Cx, *Cy, *u, *v, *ans_t;
	int i, j, size_x = 0, size_y = 0;
	
	// find the max size of these two numbers
	
	while(x[size_x] != '\0') size_x++;
	while(y[size_y] != '\0') size_y++;
	
	
	// add 0s at the first of the small one
	if(size_x > size_y){
		FinSize = size_x;
		char new_y[FinSize];
		for(i = 0; i < size_x-size_y; i++) new_y[i] = '0';
		for(i = size_x-size_y; i < size_x; i++) new_y[i] = y[i-size_x+size_y];
		y = new_y;
	}
	else if(size_x < size_y){
		FinSize = size_y;
		char new_x[FinSize];
		for(i = 0; i < size_y-size_x; i++) new_x[i] = '0';
		for(i = size_y-size_x; i < size_y; i++) new_x[i] = x[i-size_y+size_x];
		x = new_x;
	}
	else FinSize = size_x;
	
	
	// use FFT and iFFT to find the answer 
	Cx = (complex*) malloc(2*FinSize*sizeof(complex));
	Cy = (complex*) malloc(2*FinSize*sizeof(complex));
	u = (complex*) malloc(2*FinSize*sizeof(complex));
	v = (complex*) malloc(2*FinSize*sizeof(complex));
	ans_t = (complex*) malloc(2*FinSize*sizeof(complex));
	
	
	char a; //用來代替char，以便轉成整數 
	for(i = 0; i < FinSize; i++){
		a = x[i];
		Cx[FinSize-1-i] = atoi(&a); // change the char to the int
	} 
	for(i = FinSize; i < 2*FinSize; i++) Cx[i] = 0;
	

	for(i = 0; i < FinSize; i++){
		a = y[i];
		Cy[FinSize-1-i] = atoi(&a);
	} 
	for(i = FinSize; i < 2*FinSize; i++) Cy[i] = 0;
	
	//u: x after fft
	//v: y after fft
	
	//after this , only use the in fact 2*N 
	FinSize = 2*FinSize;
	FFT_p(u, Cx, FinSize);
	FFT_p(v, Cy, FinSize);
	
	for(i = 0; i < FinSize; i++){
		u[i]*=v[i];
		ans_t[i] = 0;
	} 
	
	
	iFFT_p(ans_t, u, FinSize);
	for(i = 0; i < FinSize; i++) ans_t[i] /= FinSize;

	
	//將結果個別放入array中 
	int *res = (int*)malloc(FinSize*sizeof(int));
	for(i = 0; i < FinSize; i++) res[i] = 0;
	
	char ans[3];
	int temp;
	for(i = 0; i < FinSize; i++){
		temp = (int)(creal(ans_t[i])+0.5); //+0.5 為了可以無條件捨去小數部分 //這部分可加速 
		itoa(temp, ans, 10);
		
		int size_ans;
		//size of this number
		if(ans[1] == '\0' && ans[2] == '\0') size_ans = 1;
		else if(ans[2] == '\0') size_ans = 2;
		else size_ans = 3;
		
		for(j = 0; j < size_ans; j++){
			
			a = ans[j];	
			int idx = FinSize-size_ans-i+j;
			res[idx] += atoi(&a);
			while(res[idx] >= 10){
				res[idx] -= 10;
				res[idx-1] += 1;
			}
			
		}
		
	}
	
	int idx = 0;
	while(res[idx] == 0) idx++;
	// the number is begin from idx
	FinSize = FinSize-idx;
	int *res2 = (int*)malloc((FinSize)*sizeof(int));
	for(i = 0; i < (FinSize); i++){
		res2[i] = res[i+idx];
	}
	
	
	
	return res2;
}


int main(){
	int i ;

	char a[6] = {'1', '2', '3', '4', '5', '6'};
	char b[5] = {'5', '6', '7', '8', '9'};

	int *ans = FastMul(a, b);
	
	for(i = 0; i < FinSize; i++) printf("%d", ans[i]);
	
	
	
	return 0;
} 

