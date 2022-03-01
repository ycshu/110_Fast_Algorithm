#include <stdio.h>

int main() {
	int k;
	unsigned char *c;
	double d = 6.0;
	// read double by uchar pointer (1 byte)
	c = (unsigned char*) &d;
	for(k=sizeof(d)-1;k>=0;k--) {
		printf("%02X", c[k]);
		if(k % 4 == 0) printf(" ");
	}
	printf("\n");

	unsigned int  *n;
	double v[2] = {6.0, 12.0};
	// read double by uint pointer (4 bytes)
	n = v;
	for(k=sizeof(v)/sizeof(unsigned int)-1;k>=0;k--) {
		printf("%08X ", *(n+k));
	}
	printf("\n");
	
	return 0;
}
