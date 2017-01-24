/* Test of pointers i C */
#include <stdio.h>

int main(){
    
    int a=2;
    int *pta; /* A pointer to a */
    pta = &a; /* set to the address of a */
    printf("Address of a is: %d\n",&a);
    
    /* test that *pta and a are the same value */
    
    printf("%d %d",a,*pta);
    
    return 0;
}
