#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "complex.h"

int roundup(int a, int b){
    int c;
    if(a%b){
        a += (b-a%b);
    }
    c = a/b;
    return c;
}

int main(){
	int i, j, k, tmp, N = 2;     //i,j,k for index using

    //prime factorization
    //delcare some variables
    int nextp, nextN, psize;
    int *prime, *ptime;     //record prime and square number
    int Ntmp = N, index = 0, t, ctrl, final;

    //memory allocation for prime and its square number
    psize = roundup(N,2);
    prime = (int *) malloc( psize * sizeof(int));
    ptime = (int *) malloc( psize * sizeof(int));
        
    for(i = 2; i <= N; i++){
        t = 0;
        while(Ntmp%i == 0){
            Ntmp /= i;
            t++;                //sum up square number
            ctrl = 1;           //switch the position of the array
        }

        if(ctrl == 1){
            prime[index] = i;   //prime        
            ptime[index] = t;   //square number of the prime resp.
            index++;
        }
        ctrl = 0;
    }
    prime[index] = N;
    ptime[index] = 1;
    printf("prime done\n");
    if((prime[0] == N && ptime[0] == 1)||(!index)){  //as N = p, do the final FFT
        printf("one time\n");
        final = 1;              
    }else if(ptime[index-1] > 1){
    	printf("two time\n");
        nextp = prime[index-1];
    }else{
    	printf("three time\n");
        nextp = prime[index-2];
    }
}