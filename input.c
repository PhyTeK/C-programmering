/* Test Scanf and printf */
#include <stdio.h>

int main(){
    int a,b;
    /* &a is the internal adress of a 
     * scanf gets two integers
    */
    scanf("%d%d",&a,&b); 
    printf("I got %d and %d",a,b);
    
    return 0;
}
