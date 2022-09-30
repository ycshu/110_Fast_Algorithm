#include <stdio.h>
#include <stdlib.h>


void PrintBaseP(int a, int p){ // input 1個整數轉為p進位，並印出此p進位及反序 
	if(a < p) printf("k = %d\nj = %d\n", a, a);
	else{
		int i, cnt, temp, *output, size;
		temp = a;
		cnt = 0;
		while(temp > p-1){
			cnt++;
			temp /= p;
		}
		cnt++;
		size = cnt;
		printf("cnt = %d\n", cnt);
		output = (int*)malloc(cnt*sizeof(int));//0 ~ cnt
		
		
		while(cnt > 0){
			printf("put in : %d\n", a%p);
			output[--cnt] = a%p;
			a /= p;
		}
		printf("size = %d\n",sizeof(output));
		
		printf("k = ");
		for(i = 0; i < size; i++){
			printf("%d", output[i]);
		}
		printf("\n");
		
		printf("j = ");
		for(i = size-1; i >= 0; i--){
			printf("%d", output[i]);
		}
		printf("\n");
		
		
		
	}
	
	
	
	
	
	
	 
}






int main(){
	base_p_1(100, 4);
	
	
	
	
	
	return 0;
} 
