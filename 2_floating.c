#include <stdio.h>

int main() {
	int k;
	unsigned char *c;
	float  f = 6.0;
	double d = 6.0;
	c = (unsigned char*) &f;
	c = c+sizeof(float);
	//c = c+sizeof(double);
	for(k=0;k<sizeof(float);k++) {
		c--;
		printf("%02X", *c);
	}
	printf("\n");
	return 0;
}
