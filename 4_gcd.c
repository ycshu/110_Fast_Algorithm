#include <stdio.h>
#include <time.h>	// for random seed
unsigned int GCD(unsigned int a, unsigned int b);

int main() {
	unsigned int d, k, m, n, a, b, r, t, s, N=1 << 25;
	time_t t1, t2;

	srand( time(NULL) );
	m = (rand()/2)*(rand()/2)+1;
	n = (rand()/2)*(rand()/2)+1;
	printf("%d, %d\n",m,n);
	if (n>m) {
		t = m;
		m = n;
		n = t;
	}
	// Method 1
	t1 = clock();
	for (k=1;k<N;++k) {
		a = m; b = n;
		while (b != 0) {
			r = a % b;
			a = b;
			b = r;
		}
	}
	t1 = clock()-t1;
	printf("Method 1: %d ms\n",t1);
	// Method 2
	a = m; b = n;
	t1 = clock();
	for(k=1;k<N;++k){
		d = GCD(a,b);
	}
	t1 = clock()-t1;
	printf("Method 2: %d ms\n",t1);

	// Method 3
	t2 = clock();
	for (k=1;k<N;++k) {
		a = m; b = n; d = 1;
		while (b != 0) {
			// If a and b are both even, take 2 from them and sent to gcd
			while((a % 2 == 0) && (b % 2 == 0)) {
				d = d * 2;
				a = a / 2; 
				b = b / 2;
			}
			while((a % 2 == 0) && (b % 2 == 1)) {
				a = a / 2; 
			}
			while((a % 2 == 1) && (b % 2 == 0)) {
				b  = b/2; 
			}
			r = a % b;
			a = b;
			b = r;
		}	
		a = a * d;
		//printf("2: GCD of %d and %d is %d\n",m,n,a);
	}
	t2 = clock()-t2;
	printf("Method 3: %d ms\n", t2); 
	return 100;
}
unsigned int GCD(unsigned int a, unsigned int b) {
	if (b == 0) return a;
	unsigned int r = a % b;
	if (r == 0) {
		return b;
	} else {
		return GCD(b,r);
	}
}
