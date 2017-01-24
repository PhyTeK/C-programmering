/* Pointers and arrays */
#include <stdio.h>

int main(){
    int a[]={5,10,15,20,25};
    
    int *pA;
    
    pA = &a[2];
    
    printf("%d %d\n",*pA, *(pA + 2));
    printf("Second test\n");
    
    /* Another test */
    
    int i;
    double *pB;
    double b[] = {5.1,10.2,15.3,20.4,25.5};
    
    pB = b;
    
    printf("pB value is: %lf\n",*pB);
    
    /* print ut all elements of b[] */
    
    
    for(i=0;i<5;i++){
        printf("%.2lf ",*(pB+i));
    }
    
    return 0;
}
