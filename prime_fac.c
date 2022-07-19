#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void p_fac(int N, int &p, int &q, int &r);

int main(){
	int N = 10800000, p = 0, q = 0, r = 0;
	
	p_fac(N, p, q, r);
	
	printf("p = %d, q = %d, r = %d\n", p, q, r);
	
	return 0;
}

void p_fac(int N, int &p, int &q, int &r){
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
}