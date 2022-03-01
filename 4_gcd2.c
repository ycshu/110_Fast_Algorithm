#include <stdio.h>
#include <time.h>

int main() {
	int a, b, r, k, N=1 << 25;
	time_t t1, t2;
	
	printf("Please input a b : ");
	scanf("%d %d", &a, &b);
	
	srand( time(NULL) );
	t1 = clock();
	for(k=0;k<N;++k){
		a = rand();
		b = rand();
		d = GCD(a,b);
		//printf("GCD=%d\n", GCD(a,b));
	}
	t1 = clock()-t1;

	/*	
	while ( r > 0 ) {
		// a = bq + r
		// q = a / b;
		r = a % b;
		// a --> b, b --> r 
		a = b;
		b = r;
		printf("(%d,%d)\n", a,b);
	}
	printf("GCD=%d\n",a);
	*/
	return 100;
}
int GCD(int a, int b){
	int r;
	r = a % b;
	if(r == 0) {
		return b;
	}
	else {
		printf("(%d,%d)\n",b,r);
		return GCD(b, r);
	}
}
