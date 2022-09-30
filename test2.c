#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double Residual(double **r, double **u, double **f, int N);

int main(){
	int i, j, k, s, L = 5, *M, N, N2;
	double ***Mem, ***u, ***un, ***f, ***r, res, res0, h, w; 
	time_t t1, t2;
	
	//					   [ 0, -1,  0] 
	// solving Au = f, A = [-1,  4, -1]
	//                     [ 0, -1,  0]
	 
	// for memory allocation!	
	M   = (int *) malloc((L+1)*sizeof(int));
	Mem = (double ***) malloc(5*L*sizeof(double*));
	u   = Mem + L - 1;	// -1 makes u[0~L-1] to u[1~L]
	un  = u   + L;
	f   = un  + L;
	r   = f   + L;

	N = 2;
	for(k = 1; k <= L; ++k) {
		printf("%d\n", N); // N = 2^(j+1), j=1..L
		N2 = (N+1)*(N+1);
		Mem[k-1] = (double **) malloc(4*(N+1)*sizeof(double*)); // open memory once
		Mem[k-1][0] = (double*) malloc(4*N2*sizeof(double));
		
		u[k]   = Mem[k-1];
		un[k]  = u[k]  + (N+1);
		f[k]   = un[k] + (N+1);
		r[k]   = f[k]  + (N+1);
		
		u[k][0]  = Mem[k-1][0];
		un[k][0] = u[k][0]  + N2;
		f[k][0]  = un[k][0] + N2;
		r[k][0]  = f[k][0]  + N2;
		
		for(i = 1; i <= N; ++i){
			u[k][i]  = u[k][i-1] + N + 1;
			un[k][i] = un[k][i-1]  + N + 1;
			f[k][i]  = f[k][i-1]   + N + 1;
			r[k][i]  = r[k][i-1]   + N + 1;
		}
		
		h = 1.0/N;
		
		for(i = 0; i <= N; ++i){
			for(j = 0; j <= N; ++j){
				u[k][i][j]  = 0.0;
				un[k][i][j] = 0.0;
				r[k][i][j]  = 0.0;
				f[k][i][j]  = -1.0*h*h;
			}
			
		}
		M[k] = N;
		N = N*2;
	}	
	k = L;
	res0 = Residual(r[k],u[k],f[k],M[k]);
	w = 0.9;
	for(k = 0 ; k < 10; ++k){
		for(j = L; j > 1; --j){
			SmootherG(un[j],u[j],f[j],w,M[j]);
			Residual(r[j],u[j],f[j],M[j]);
			Restriction(r[j], f[j-1], M[j-1]);
			Initial(u[j-1],M[j-1]);
		}
		SmootherJ(un[j],u[j],f[j],1.0,M[j]);
		for(j = 2; j <= L; ++j){
			Interpolation(un[j], u[j-1], M[j]);
			Evolution(u[j],un[j],M[j]);
			SmootherG(un[j],u[j],f[j],w,M[j]);			
		}
		j = L;
		res = Residual(r[j],u[j],f[j],M[j]);
		printf("%d Res: %e Res_reduce: %.3f\n",k,res,res/res0);
		res0 = res;
	}
	
	
	for(i=0;i<L;++i) free(Mem[i]);
 	free(Mem);
	free(M);

	return 0;
	
}

int Initial(double **u, int N){
	int i, j;
	for(i = 1; i < N; ++i){
		for(j = 1; j < N; ++j){
			u[i][j] = 0.0;
		} 
	} 
	return 0;
}

int Evolution(double **u, double **un, int N){
	int i, j;
	for(i = 1; i < N; ++i){
		for(j = 1; j < N; ++j){
			u[i][j] += un[i][j];
		}
	} 
	return 0;
}
double Residual(double **r, double **u, double **f, int N){
	// compute r = f - Au
	double maxres = 0.0;
	int i, j;
	for(i = 1; i < N; ++i){
		for(j = 1; j < N; ++j){
			r[i][j] = f[i][j] - (u[i+1][j] - 4.0*u[i][j] + u[i-1][j] + u[i][j+1] + u[i][j-1]);
		if(fabs(r[i][j])>maxres) maxres = fabs(r[i][j]);
		}
		
		
	}
	return maxres;
}

int SmootherJ(double **un, double **u, double **f, double omega, int N){
	int i, j;

	for(i = 1; i < N; ++i){
		for(j = 1; j < N; ++j){
			un[i][j] = u[i][j] + omega*(0.25)*( (u[i+1][j] - 4.0*u[i][j] + u[i-1][j] + u[i][j+1] + u[i][j-1]) - f[i][j] );

		}
	}
	for(i = 1; i < N; ++i){
		for(j = 1; j < N; ++j){
			u[i][j] = un[i][j];
		}
	} 
	
	return 0;
}

int SmootherG(double **un, double **u, double **f, double omega, int N){
	int i, j;
	omega = 1.0;
	
	for(i = 1; i < N; ++i){
		for(j = 1; j < N; ++j){
			u[i][j] = u[i][j] + omega*(0.25)*( (u[i+1][j] - 4.0*u[i][j] + u[i-1][j] + u[i][j+1] + u[i][j-1]) - f[i][j] );

		}
	}
	
	return 0;
}

int Restriction(double **r_f, double **r_c, int Nc) {
	// 
	int i, j;
	for(i = 1; i < Nc; ++i){
		for(j = 1; j < Nc; ++j){
			r_c[i][j] = ( r_f[2*i-1][2*j-1]   + 2.0*r_f[2*i][2*j-1] + r_f[2*i+1][2*j-1]+
						  2.0*r_f[2*i-1][2*j] + 4.0*r_f[2*i][2*j]   + 2.0*r_f[2*i+1][2*j]+
						  r_f[2*i-1][2*j+1]   + 2.0*r_f[2*i][2*j+1] + r_f[2*i+1][2*j+1] )/4.0;
			
		}
		
	}
	return 0;
}

int Interpolation(double **u_f, double **u_c, int Nf) {
	// 
	int i, j;
	
	for(i = 1; i < Nf; ++i){
		for(j = 1; j < Nf; ++j){
			
			if(i%2 == 0 && j%2 == 0){
				u_f[i][j] = u_c[i/2][j/2];
			}
			else if(i%2 == 1 && j%2 == 0){
				u_f[i][j] = 0.5*(u_c[i/2][j/2] + u_c[i/2+1][j/2]);
			}
			else if(i%2 == 0 && j%2 == 1){
				u_f[i][j] = 0.5*(u_c[i/2][j/2] + u_c[i/2][j/2+1]);
			}
			else{
				u_f[i][j] = 0.25*( u_c[i/2][j/2]   + u_c[i/2][j/2+1] +
								   u_c[i/2+1][j/2] + u_c[i/2+1][j/2+1]);
			}
		}
	}
	
	return 0;
}


