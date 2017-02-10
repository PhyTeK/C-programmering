#include <stdio.h>
#define MAX 1000

int main(){
    int i,sum=0;
    
    
    for(i=1;i<=MAX;i++){
        sum = sum + i;
    }

    printf("Sum f.o.m %d t.o.m %d is: %d\n",1,1000,sum);

    return 0;
}
