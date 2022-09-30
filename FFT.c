#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
int FFT(double *yr, double *yi, double *xr, double *xi, int N){
	
	if(N == 2){
		yr[0] = xr[0] + xr[1];
		yi[0] = xi[0] + xi[1];
		yr[1] = xr[0] - xr[1];
		yi[1] = xi[0] - xi[1];
		return 0;
	}
	else{
		int n, k, c, s;
		double *yEr, *yEi, *yOr, *yOi, *xEr, *xEi, *xOr, *xOi;
		
		xEr = (double*) malloc((N/2)*sizeof(double));
		xEi = (double*) malloc((N/2)*sizeof(double));
		yEr = (double*) malloc((N/2)*sizeof(double));
		yEi = (double*) malloc((N/2)*sizeof(double));
		xOr = (double*) malloc((N/2)*sizeof(double));
		xOi = (double*) malloc((N/2)*sizeof(double));
		yOr = (double*) malloc((N/2)*sizeof(double));
		yOi = (double*) malloc((N/2)*sizeof(double));
		for(n = 0; n < N/2; n++){
			xEr[n] = xr[2*n];
			xEi[n] = xi[2*n];
			xOr[n] = xr[2*n + 1];
			xOi[n] = xi[2*n + 1];
		}
		FFT(yEr, yEi, xEr, xEi, N/2);
		FFT(yEr, yEi, xEr, xEi, N/2);
		
		for(k = 0; k < N/2; k++){
			c =  cos(2*M_PI*k/N);
			s = -sin(2*M_PI*k/N);
			
			yr[k]       = yEr[k] + (yOr[k]*c - yOi[k]*s);
			yi[k]       = yEr[k] + (yOr[k]*s + yOi[k]*c);
			yr[N/2 + k] = yEr[k] - (yOr[k]*c - yOi[k]*s);
			yi[N/2 + k] = yEi[k] - (yOr[k]*s + yOi[k]*c);
		}
		free(xEr); free(xEi);
		free(yEr); free(yEi);
		free(xOr); free(xOi);
		free(yOr); free(yOi);
	}
}





int main(){
	FILE *fp;
	char ch;
	int L = 0, N = 1;
	int k, n;
	double *xr, *xi, *yr, *yi, *zr, *zi, *ur, *ui; 
	time_t T; 
	
	
	fp = fopen("hw6.csv", "r");
	while(!feof(fp)){
		ch = fgetc(fp);
		if(ch == '\n') L++; 
	}
	fclose(fp);
	
	while(N < L) N = N * 2;
	N = N/2;
	
	xr = (double*) malloc(N*sizeof(double));
	xi = (double*) malloc(N*sizeof(double));
	yr = (double*) malloc(N*sizeof(double));
	yi = (double*) malloc(N*sizeof(double));
	zr = (double*) malloc(N*sizeof(double));
	zi = (double*) malloc(N*sizeof(double));
	ur = (double*) malloc(N*sizeof(double));
	ui = (double*) malloc(N*sizeof(double));
	
	fp = fopen("hw6.csv", "r");
	for(k = 0; k < N; k++){
		fscanf(fp, "%d %d", &n, &L);
		xr[k] = 1.0*L;
		xi[k] = 0.0;
	}
	fclose(fp);
	
	T = clock();
	FFT(zr, zi, xr, xi, N);
	T = clock() - T;
	printf("time = %d ms \n", T);
	
	
	return 0;
} 
