#include <stdio.h>

int main() {
	int a, b, q, r;
	
	printf("Please input a b : ");
	scanf("%d %d", &a, &b); 
	r = 1;
	while(r > 0) {
		q = a / b;
		r = a % b;
		a = b;
		b = r;
		printf("(%d,%d)\n",a,b);
	}
	printf("GCD: %d\n", a);
	return 100;
}
