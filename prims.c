#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <errno.h>

/*#define MAX 10000000*/

int main(int argc,char *argv[]){
    int i,j,count=0;
    int prime,max;
    char *eptr;
    long long N;
    
    printf("Args: %d %s\n",argc, argv[0]);
    if(argc == 2){
         N=strtoll(argv[1],&eptr,10);
    }else{
        printf("%s\n","Argument should be an integer");
        return 1;
    }

    printf("N=%lld\n",N);
    for(i=1;i<=N;i++){
        prime = 1;
        max = (int)sqrt((double)i);
        /*printf("%d\n",max);*/
        for(j=2;j<=max;j++){
            if(i%j == 0){                
                prime=0;
                break;
            }            
        }
    
        if (prime == 1)
            count++;
           
        
    }
    printf("There is %d primes numbers between %d and %lld\n",count-1,2,N);
    return 0;
}
