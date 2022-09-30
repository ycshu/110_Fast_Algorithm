#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <string.h>


int Initialize(double *u, int N){
	int i;
	for(i = 0; i < N; ++i) u[i] = 0.0;
	
	return 0;
} 

int Evolution(double *u, double *un, int N){
	int i;
	for(i = 0; i < N; ++i) u[i] += un[i];
	return 0;
}

double Residual(double *r, double *u, double *f, int N){
	//compute r = f-Au
	int i;
	double v = 0.0;
	int m = (int)sqrt(N);
	
	
	for(i = 1; i < N; i++){
		if(i < m){
			r[i] = f[i] - (4*u[i] - u[i-1] - u[i+1] - u[i+m]);
		}
		else if(i >= m*(m-1)){
			r[i] = f[i] - (4*u[i] - u[i+1] - u[i-1] - u[i-m]);
		}
		else{
			r[i] = f[i] - (4*u[i] - u[i-1] - u[i-m] - u[i+1] - u[i+m]);
		}
		if(fabs(r[i]) > v) v = fabs(r[i]);
	}
	
	return v;
}

int Smooth(double *un, double *u, double *f, double omega, int N){
	int i;
	double v = 0.0;
	int m = (int)sqrt(N);
	
	
	for(i = 1; i < N; i++){
		if(i < m){
			un[i] = u[i] + omega*(0.5)*(4*u[i] - u[i-1] - u[i+1] - u[i+m] - f[i]);
		}
		else if(i >= m*(m-1)){
			un[i] = u[i] - omega*(0.5)*(4*u[i] - u[i+1] - u[i-1] - u[i-m] - f[i]);
		}
		else{
			un[i] = u[i] - omega*(0.5)*(4*u[i] - u[i-1] - u[i-m] - u[i+1] - u[i+m] - f[i]);
		}
	}
	
	for(i = 1; i < N; i++){
		u[i] = un[i];
	}
	
	return 0;
}

int Restriction(double *r_f, double *r_c, int Nc){
	// go to coarse grid
	//r_c:coarse
	//r_f:fine
	
	int i, j;
	int m = ((int)sqrt(Nc))/2;
	
	j = 0;
	for(i = m+1; i < Nc; ++i){
		//跳過第一列 & 最後一列 
		if(i%m == m-1){
			//跳過相鄰的一列 
			i = i+m+1;
			continue;
		}
		if(i&m == 0){
			i = i+m;
			continue;
		}
		//i 為中心的8個方向都要算 ,然後除16 
		r_c[j] = 0.0625*(4*r_f[i] + 2*(r_f[i-m] + r_f[i-1] + r_f[i+1] + r_f[i+m]) + r_f[i-m-1] + r_f[i-m+1] + r_f[i+m-1] + r_f[i+m+1]);
		++j;
		++i;
	}
	
	return 0;
}

int Interprolation(double *u_f, double *u_c, int Nf){
	//go to fine grid
	//r_c:coarse
	//r_f:fine
	
	
	int i, j;
	int mf = (int)sqrt(Nf);
	int mc = mf/2;
	
	if(mf%2 == 0){
		//偶數
		int row = 1;
		
		for(i = 0; i < mc; ++i){
			for(j = 0; j < mc; ++j){
				int InFineIdx = mf*row + 2*j + 1;
				int InCoaresIdx = mc*i+j;
				
				u_f[InFineIdx] = u_c[InCoaresIdx];//中心 
				
				if(j != mc-1){
					u_f[InFineIdx+1] += u_c[InCoaresIdx]/2;//右 
					u_f[InFineIdx+1-mf] += u_c[InCoaresIdx]/4;//右上 
					u_f[InFineIdx+1+mf] += u_c[InCoaresIdx]/4;//右下 
				}
				
				if(i != mc-1){
					u_f[InFineIdx+mf] += u_c[InCoaresIdx]/2;//下 
					u_f[InFineIdx+1+mf] += u_c[InCoaresIdx]/4;//左下 
					u_f[InFineIdx-1+mf] += u_c[InCoaresIdx]/4;//右下 
				}
				
				u_f[InFineIdx-mf] += u_c[InCoaresIdx]/2;//上 
				u_f[InFineIdx-1] += u_c[InCoaresIdx]/2;//左
				u_f[InFineIdx-1-mf] += u_c[InCoaresIdx]/4;//左上 
				
				if(i == mc-1 && j == mc-1){
					u_f[InFineIdx-1+mf] -= u_c[InCoaresIdx]/4;//右下被多算一次 
				}
			}
			row += 2;
		}
		
	}
	else{
		//奇數
		int row = 1;
		for(i = 0; i < mc; ++i){
			for(j = 0; j < mc; ++j){
				int InFineIdx = mf*row + 2*j + 1;
				int InCoaresIdx = mc*i+j;
				
				u_f[InFineIdx-mf-1] = u_c[InCoaresIdx]/4;//左上 
				u_f[InFineIdx-mf] = u_c[InCoaresIdx]/2;//上 
				u_f[InFineIdx-mf+1] = u_c[InCoaresIdx]/4;//右上 
				u_f[InFineIdx-1] = u_c[InCoaresIdx]/2;//左
				u_f[InFineIdx] = u_c[InCoaresIdx];//中心 
				u_f[InFineIdx+1] = u_c[InCoaresIdx]/2;//右 
				u_f[InFineIdx+mf-1] = u_c[InCoaresIdx]/4;//左下 
				u_f[InFineIdx+mf] = u_c[InCoaresIdx]/2;//下 
				u_f[InFineIdx+mf+1] = u_c[InCoaresIdx]/4;//右下 
			}
		}
	}
	
	return 0;
}

void Source(double **F, int N){
	// make F
	
	int i, j, k;
	double x, y, h;
	h = 1.0/N;
	
	for(i = 0; i < N-1; ++i){
		for(j = 0; j < N-1; ++j){
			x = (i+1)*h;
			y = (j+1)*h;
			F[i][j] = -(1.0+4.0)*h*h*M_PI*M_PI*sin(M_PI*x)*sin(2*M_PI*y);
		}
	}
}



int main(){
	int i, j, k, s, *Level, N, M;
	double **MEM, **f, **u, **un, **r;
	double res, res0, w, h;
	
	// solving Au = f
	N = 4;
	M = (N-1)*(N-1);
	
	int levelSize = 0;
	int tempN = N;
	while(tempN != 0){
		++levelSize;
		tempN /= 2;
	}
//	printf("size is %d\n", levelSize);
	//open memory
	Level = (int*)malloc((levelSize+1)*sizeof(int));
	MEM = (double**)malloc(5*levelSize*sizeof(double*));
	u = MEM + levelSize;
	un = u + levelSize;
	f = un + levelSize;
	r = f + levelSize;
	
	N = 2;
	for(j = 1; j < levelSize; ++j){
		printf("%d\n", N);
		// open memory once
		MEM[j-1] = (double*)malloc(4*(N+1)*sizeof(double));
		
		u[j] = MEM[j-1];
		un[j] = u[j] + N + 1;
		f[j] = un[j] + N+1;
		r[j] = f[j] + N+1;
		
		h = 1.0/N;
		
		
	}
	
	
	
	
	
	// end of open memory
	
	
	//make F
	
	
	
	
//	for(i = 0; i < levelSize; i++) printf("%d\n", Level[i]);
	
	
	//begin the multigrid
//	res0 = Residual(r, u, f, Level[levelSize]);
//	printf("%f\n", res0);
	
//	for(i = 0; i < M; i++) printf("%f  %f  %f\n",r[i], u[i], b[i]);
	
	
//	w = 2.0/3;
//	
//	for(k = 0; k < 10; ++k){
//		for(j = levelSize; j > 0; --j){
//			Smooth(un, u, f, w, Level[j]);
//			Residual(r, u, f, Level[j]);
//			
//			
//		}
//	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	return 0;
}
