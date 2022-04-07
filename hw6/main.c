#include <stdio.h>
#include <stdlib.h>
#include <math.h>


struct __array_t{
    unsigned int size;
    unsigned int capacity;
    double *re;
    double *im;
};
typedef struct __array_t array_t;

#define ARRAY_INIT (array_t){.size = 0, .capacity = 0, .re=NULL, .im=NULL}

void push_data(array_t * ar, int data){
    if(ar == NULL){
        return; 
    }

    if(ar->size + 1 >= ar->capacity){
        ar->re = realloc(ar->re, sizeof(double)*(ar->capacity + 1000));
        ar->im = realloc(ar->im, sizeof(double)*(ar->capacity + 1000));
        ar->capacity += 1000;
    }

    if(ar->re && ar->im){
        ar->re[ar->size] = data;    
        ++ar->size;
    }
}

void prealloc_array(array_t *ar, unsigned int num){
    ar->capacity = num;
    ar->re = malloc(sizeof(double)*num);
    ar->im = malloc(sizeof(double)*num);

    if(!ar->re || !ar->im){
        printf("Malloc failed in prealloc_array function\n");
        exit(-1);
    }
}

void free_array(array_t ar){
    if(ar.re && ar.im){
        free(ar.re);
        free(ar.im);
    }
}

void FFT2(double *xre, double *xim, double *yre, double *yim){
    // y[0] = x[0] + x[1];
    yre[0] = xre[0] + xre[1];
    yim[0] = xim[0] + xim[1];

    // y[1] = x[0] - x[1];
    yre[1] = xre[0] - xre[1];
    yim[1] = xim[0] - xim[1];
}

void FFTN(array_t x, array_t y, int N){

    double cos_value, sin_value;

    if(N & 2){
        FFT2(x.re, x.im, y.re, y.im);
    }else{
        array_t xodd, xeven;
        array_t yodd, yeven;
        prealloc_array(&xodd, N>>1l);
        prealloc_array(&xeven, N>>1l);

        prealloc_array(&yodd, N>>1l);
        prealloc_array(&yeven, N>>1l);
        

        // place the data
        for(int i = 0; i < N/2; ++i){
            xodd.re[i] = x.re[(i<<1) | 1];
            xodd.im[i] = x.im[(i<<1) | 1];

            xeven.re[i] = x.re[i<<1];
            xeven.im[i] = x.im[i<<1];
        }

        FFTN(xeven, yeven, N>>1l);
        FFTN(xodd, yodd, N>>1l);

        // collect the data
        for(int i = 0; i <N/2; ++i){
            cos_value = cos(2*M_PI*i/N);
            sin_value = -sin(2*M_PI*i/N);
            
            // for the value k < N/2
            y.re[i] = yeven.re[i] + (yodd.re[i] *cos_value - yodd.im[i] * sin_value);
            y.im[i] = yeven.im[i] + (yodd.re[i] *sin_value + yodd.im[i] * cos_value);
            // for the value k >= N/2
            y.re[N/2 + i] = yeven.re[i] - (yodd.re[i] *cos_value - yodd.im[i] * sin_value);
            y.im[N/2 + i] = yeven.im[i] + (yodd.re[i] *sin_value + yodd.im[i] * cos_value);
        }

        free_array(xeven);
        free_array(xodd);
        free_array(yeven);
        free_array(yodd);
    }
}


int main(){
    FILE * file = fopen("hw6.csv", "r");
    if(file == NULL){
        perror("Can't open file hw6.csv");
        exit(-1);
    }
    
    array_t array = ARRAY_INIT;
    
    int index, data;
    while(fscanf(file,"%d,%d", &index, &data) != EOF){
        push_data(&array, data);
    }
    fclose(file);

    for(int i = 0; i < array.size; ++i){
        printf("%d,%.2f\n", i, array.re[i]); 
    }
    
    free_array(array);
    

    // truncate the data
    array.size = 1l << (31l - __builtin_clz(array.size));
    array_t array_res;
    prealloc_array(&array_res, array.size);
    
    FFTN(array, array_res, array.size);
    array_res.size = array.size; 
    for(int i = 0; i < array_res.size; ++i){
        printf("%d,%.3f%s%.3f\n", i, array_res.re[i], (array_res.im[i]>0) ? "+": "", array_res.im[i]);
    }
}
