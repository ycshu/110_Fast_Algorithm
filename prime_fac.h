#include <stdio.h>
#include <stdlib.h>

void p_fac(int N, int &p, int &q, int &r, int &ctrl){
	int O = N;

	while(N%5==0){
		r++;
		N /= 5;
	}
	while(N%3==0){
		q++;
		N /= 3;
	}
	while(N%2==0){
		p++;
		N /= 2;
	}

	if( ( p==0 & q==0 & r==0 ) || N!=1 ){
		ctrl = 1;
		printf("The number %d doesn't have prime factorization in 2,3,5\n",O);
	}
}

int roundup(int a, int b){
    int c;
    if(a%b){
        a += (b-a%b);
    }
    c = a/b;
    return c;
}
