#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int FFT(double *yr, double *yi, double *xr, double *xi, int N);

int main(){
    int k,n,N;
    double *xr, *xi, *yr, *yi, *zr, *zi, *ur, *ui, t, a, b, an, bn, c, s;
    time_t T;

    FILE *fp;
    char ch;
    int L = 0;
    fp = fopen("wave.csv","r");
    while(!feof(fp)){
        ch = fgetc(fp);
        if(ch=='\n'){
            L++;
        }
    }
    fclose(fp);
    printf("Total Lines: %d \n",L);

    N = 1;
    while(N<L){
        N = N * 2;
    }
    N = N / 2;
    printf("Line: %d, Use: %d \n",L,N);

    xr = (double *) malloce(N*sizeof(double));
    xi = (double *) malloce(N*sizeof(double));
    yr = (double *) malloce(N*sizeof(double));
    yi = (double *) malloce(N*sizeof(double));
    zr = (double *) malloce(N*sizeof(double));
    zi = (double *) malloce(N*sizeof(double));
    ur = (double *) malloce(N*sizeof(double));
    ui = (double *) malloce(N*sizeof(double));

    printf("Memory Done! \n");

    fp = fopen("file","r");
    for(k = 0;k<N;k++){
        fscanf(fp,"%d %d", &n, &L);
        xr[k] = 1.0*L;
        xi[k] = 0.0;
    }
    fclose(fp);

    for(k = 0;k<5;++k){
        printf("%d: %f \n",k,xr[k]);
    }

    T = clock();

    a = cos(2 * M_PI / N); b = sin(2 * M_PI / N);
    an = 1; bn = 0;
    for(n = 0; n < N; n++){
        yr[n] = 0.0;
        yi[n] = 0.0;
        c = 1; s = 0;
        for(k = 0; k<N; ++k){
            yr[n] += c*xr[k];
            yi[n] -= s*xr[k];
            t = c * an - bn * b;
            s = c * bn + s * an;
            c = t;
        }
        t = an * a - bn * b;
        bn = an * b + bn * a;
        an = t;
    }

    T = clock() - T;
    printf("%d ms for DFT of %d elements\n",T,N);

    T = clock();
    FFT(zr, zi, xr, xi, N);
    T = clock() - T;
    printf("%d ms for FFT of %d elements\n", T, N);

    t = 0;
    for(k = 0; k<N; ++k){
        t += pow((yr[k] - zr[k]),2) + pow((yi[k] - zi[k]),2);
    }
    printf("Error: %e\n",sqrt(t/N));

    t = 0;
    for(k = 0; k<N; ++k){
        if(sqrt(zr[k]*zr[k]+zi[k]*zi[k])>t){
            n = k;
            t = sqrt(zr[k]*zr[k]+zi[k]*zi[k]);
        }
    }
    printf("mazimum at %d, value = %f\n", n, t/N*2);
    printf("BPM %f\n",60.0*n/(N/240)); //240 = Sample Rate

    free(xr);free(xi);free(yr);free(yi);

    return 0;

}

int FFT(double *yr, double *yi, double *xr, double *xi, int N){
    if(N == 2){
        yr[0] = xr[0]+xr[1];
        yi[0] = xi[0]+xi[1];
        yr[1] = xr[0]-xr[1];
        yr[1] = xi[0]-xi[1];
    } else {
        int k;
        double *yEr, *yEi, *yOr, *yOi;
        double *xEr, *xEi, *xOr, *xOi;
        double c, s;

        xEr = (double *) malloc((N/2)*sizeof(double));
        xEi = (double *) malloc((N/2)*sizeof(double));
        yEr = (double *) malloc((N/2)*sizeof(double));
        yEi = (double *) malloc((N/2)*sizeof(double));
        xOr = (double *) malloc((N/2)*sizeof(double));
        xOi = (double *) malloc((N/2)*sizeof(double));
        yEr = (double *) malloc((N/2)*sizeof(double));
        yOi = (double *) malloc((N/2)*sizeof(double));
        for(k = 0; k<N/2; ++k){
            xEr[k] = xr[2*k];
            xEi[k] = xi[2*k];
            xOr[k] = xr[2*k+1];
            xOi[k] = xi[2*k+1];

        }
        FFT(yEr, yEi, xEr, xEi, N/2);
        FFT(yOr, yOi, xOr, xOi, N/2);
        for(k = 0; k < N/2; ++k){
            c = cos(2*M_PI*k/N);
            s = -sin(2*M_PI*k/N);
            yr[k] = yEr[k] + (yOr[k]*c - yOi[k]*s);
            yi[k] = yEi[k] + (yOr[k]*s + yOi[k]*c);
            yr[N/2+k] = yEr[k] - (yOr[k]*c - yOi[k]*s);
            yi[N/2+k] = yEi[k] - (yOr[k]*s + yOi[k]*c);
        }
        
        free(xEr);free(xEi);
        free(xOr);free(xOi);
        free(yEr);free(yEi);
        free(yOr);free(yOi);
    }
    return 0;
}