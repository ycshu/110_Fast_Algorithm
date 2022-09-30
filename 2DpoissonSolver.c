#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <string.h>


void Exact_Discretization(double **A, int N){
	int i,j,k;
	for(i=0;i<N-1;++i){
		for(j=0;j<N-1;++j){
			
			k = i + j*(N-1);
			
			A[k][k] = -4;
			
			if(j>0) A[k][k-(N-1)] = 1;
			if(i>0) A[k][k-1] = 1;
			if(i<N-2) A[k][k+1] = 1;
			if(j<N-2) A[k][k+(N-1)] = 1;
		}
	}
	
}

void ExactSolution(double **U, int N){
	//the ExactSoltion : matrix U
	
	int i, j, k;
	double x, y, h;
	h = 1.0/N;
	for(i = 0; i < N-1; ++i){
		for(j = 0; j < N-1; j++){
			x = (i+1)*h;
			y = (j+1)*h;
			U[i][j] = sin(M_PI*x)*sin(2*M_PI*y);
		}
	}
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

void BandGaussianElimination(double **A, double *x, double *b, int N){
	//N by N matrix
	//solve Ax = b
	int i, j, k, n;
	double v;
	
	double **M;
	M = (double**)malloc(N*sizeof(double*));
	M[0] = (double*) malloc(N*N*sizeof(double));
	for(i = 1; i < N; ++i){
		M[i] = M[i-1] + N;
	}
	
	for(i = 0; i < N; ++i){
		for(j = 0; j < N; ++j){
			M[i][j] = A[i][j];
		}
	}
	
	for(i = 0; i < N; ++i){
		x[i] = b[i];
	}
	// start to elimination
	// up to down
	for(k = 0; k < N; ++k){
		for(i = k+1; i < N; ++i){
			v = M[i][k]/M[k][k];
			x[i] -= v*x[k];
			
			for(j = k; j < N; ++j){
				M[i][j] -= v*M[k][j];
			}
		}
	}
	
	//backward
	
	for(i = N-1; i >= 0; --i){
		for(j = i+1; j < N; ++j){
			x[i] -= M[i][j]*x[j];
		}
		x[i] /= M[i][i];
	}
		
	// free memory
	free(M[0]);
	free(M);
}

double Residual(double *r, double **A, double *x, double *b, int N){
	//compute r = b-Ax
	int i, j;
	double v = 0.0;
	for(i = 0; i < N; ++i){
		r[i] = b[i];
		
		for(j = 0; j < N; ++j){
			r[i] -= A[i][j]*x[j];
		}
		
		if(fabs(r[i]) > v) v = fabs(r[i]);
	}
	return v;
}

double error(double *x, double *u, int N){
	int i;
	double e, v = 0.0;
	
	for(i = 0; i < N; ++i){
		e = fabs(x[i]-u[i]);
		if(e > v) v = e;
	}
	return v;
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
	
	for(j = 0; j < p; ++j){
		for(k = 0; k < bias; ++k){
			Nx[k + j*bias] = x[p*k + j];
		}
	}
	
	// F[p][p] 
	
//	complex F[p][p];

	complex **F;
	F = (complex**)malloc(p*sizeof(complex*));
	for(i = 0; i < p; ++i){
		F[i] = (complex*)malloc(p*sizeof(complex));
	}
	
	for(i = 0; i < p; ++i){
		F[0][i] = 1;
		F[i][0] = 1;
	}
	for(i = 1; i < p/2+1; ++i){
		for(j = 1; j < p; ++j){
			F[i][j] = cexp(-2*M_PI*i*j*I/p);
		}
	}
	for(i = p/2+1; i < p; ++i){	//有反對稱
		for(j = 1; j < p; ++j){
			F[i][j] = F[p-i][p-j];
		}
	}
	
	
	//recursive
	if(N == p){
		for(i = 0; i < p; ++i){
			y[i] = 0;
			for(j = 0; j < p; ++j){
				y[i] += F[i][j]*x[j];
			} 
		}
	}
	else{
		
		for(j = 0; j < p; ++j){
			FFT_p(Ny + j*bias, Nx + j*bias, bias);
		}
		
		for(i = 0; i < bias; ++i){
			for(k = 0; k < p; ++k){
				Ny[i + k*bias] *= cexp(-2*M_PI*k*i*I/N);
			}
			
			for(j = 0; j < p; ++j){
				idx = i + j*bias;
				y[idx] = 0;
				for(k = 0; k < p; ++k){
					y[idx] += F[j][k]*Ny[i + k*bias];
				}
				
			}
		}
		
	}
	
	free(Nx); free(Ny);
} 

void DST(double *y, double *x, int N){
	//y is the answer
	
	int i;
	N = 2*N + 2;
	complex *Nx, *Ny;
	Nx = (complex*) malloc(N*sizeof(complex));
	Ny = (complex*) malloc(N*sizeof(complex));
	
	for(i = 0; i < N; ++i){
		Nx[i] = 0.0;
	}
	
	for(i = 0; i < N/2-1; ++i){
		Nx[i+1] = x[i];
		Nx[(i+1) + N/2] = -x[N/2-(i+2)];
	}
	
	
	FFT_p(Ny, Nx, N);
	
	
	for(i = 0; i < N/2-1; ++i){
		y[i] = -0.5*cimag(Ny[i+1]);
	}
	
	free(Nx); free(Ny);
	
}

void iDST(double *y, double *x, int N){
	//y is the answer
	int i;
	N = 2*N + 2;
	complex *Nx, *Ny;
	Nx = (complex*) malloc(N*sizeof(complex));
	Ny = (complex*) malloc(N*sizeof(complex));
	
	for(i = 0; i < N; ++i){
		Nx[i] = 0.0;
	}
	
	for(i = 0; i < N/2-1; ++i){
		Nx[i+1] = x[i];
		Nx[(i+1) + N/2] = -x[N/2-(i+2)];
	}
	
	
	FFT_p(Ny, Nx, N);
	
	
	for(i = 0; i < N/2-1; ++i){
		y[i] = (-2.0/N)*cimag(Ny[i+1]);
	}
	
	free(Nx); free(Ny);
}

void Transpose(double **A, int N){
	//transpose a NxN matrix
	int i , j;
	double temp;
	
	for(i = 0; i < N; ++i){
		for(j = i+1; j < N; ++j){
			temp = A[i][j];
			A[i][j] = A[j][i];
			A[j][i] = temp;
		}
	}
}

void DST2D(double **X, int N){
	int i;
	
	for(i = 0; i < N; ++i){
		DST(X[i], X[i], N);
	}
	
	Transpose(X, N);
	
	for(i = 0; i < N; ++i){
		DST(X[i], X[i], N);
	}
	
	Transpose(X, N);
}

void iDST2D(double **X, int N){
	int i;
	
	for(i = 0; i < N; ++i){
		iDST(X[i], X[i], N);
	}
	
	Transpose(X, N);
	
	for(i = 0; i < N; ++i){
		iDST(X[i], X[i], N);
	}
	
	Transpose(X, N);
}

void FastPoissonSolver(double **F, int N){
	int i, j;
	double h;
	
	h = 1.0/(N+1);
	DST2D(F, N);
	
	
	
	for(i = 0; i < N; ++i){
		for(j = 0; j < N; ++j){
			F[i][j]= F[i][j]/((2*(cos(M_PI*(i+1)*h)-1))+(2*(cos(M_PI*(j+1)*h)-1)));
		}
	}
	
	
	iDST2D(F, N);
	
}

int main(){
	//solve AU = F, where r = b - Ax is the residue
	// U is unknown, F is known
	int i, j, k, N, M;
	double **A, *x, *u, **U, *b, **F, *r;
	time_t t;
	
	N = 32;
	
	M = (N-1)*(N-1);
	//create memory, 將2維array開成一列 
	A = (double**) malloc(M*sizeof(double*));
	A[0] = (double*) malloc(M*M*sizeof(double));
	for(i = 1; i < M; ++i){
		A[i] = A[i-1] + M;
	}
	x = (double*) malloc(M*sizeof(double));
	r = (double*) malloc(M*sizeof(double));
	
	
	b = (double*) malloc(M*sizeof(double));
	F = (double**)malloc((N-1)*sizeof(double*));
	F[0] = b;
	for(i = 1; i < N-1; ++i){
		F[i] = F[i-1] + N - 1;
	}
	
	
	u = (double*) malloc(M*sizeof(double));
	U = (double**) malloc((N-1)*sizeof(double*));
	U[0] = u;
	for(i = 1; i < N-1; ++i){
		U[i] = U[i-1] + N - 1;
	}
	//end of open memory
	
	
	//initialize
	for(i = 0; i < M*M; ++i){
		A[0][i] = 0.0;
	}
	//end of initialize
	
	//make the matrix A
	Exact_Discretization(A, N);

	
	//make F
	Source(F, N);
	
	//case1. use Gauss elimination solve
	t = clock();
	BandGaussianElimination(A, x, b, M);

	t = clock() - t;
	printf("use Gauss elimination: %f sec\n", 1.0*t/CLOCKS_PER_SEC);
	printf("the residual : %e\n", Residual(r, A, x, b, M));
	
	//calculate the true answer in order to find the error
	ExactSolution(U, N);
	printf("error is: %e\n", error(x, u, M));
	
	
	//case 2. use DST fast poisson solver
	t = clock();
	FastPoissonSolver(F, N-1);
	t = clock() - t;
	
	printf("use DST fast poisson solver: %f sec\n", 1.0*t/CLOCKS_PER_SEC); 
	printf("error is: %e\n", error(b, u, M));
	
	
	//free memory
	free(A[0]);
	free(A);
	free(x);
	free(r);
	free(b);
	free(u);
	free(U);
	free(F);
	
	
	return 0;
}
