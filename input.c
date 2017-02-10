/* Test Scanf and printf */
#include <stdio.h>

int main(){
    int a,b;
    /* &a is the internal adress of a 
     * scanf gets two integers
    */
    scanf("enter two integers\n%d%d",&a,&b);
    printf("address for a is: %d",&a);
    printf("I got %d and %d\n",a,b);
    printf("The sum is equal to %d",a+b);
    
    return 0;
}
