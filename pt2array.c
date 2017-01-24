/* Pointers and arrays */
#include <stdio.h>

int main(){
    int a[]={5,10,15,20,25};
    
    int *pA;
    
    pA = &a[2];
    
    printf("%d %d\n",*pA, *(pA + 2));
    printf("Second test\n");
    
    /* Another test */
    
    int i,row,col;
    double *pB;
    double b[] = {5.1,10.2,15.3,20.4,25.5};
    int c[2][4] = {
        {10,11,12,13},
        {14,15,16,17}
    }; /* 2D array declaration and initialization*/
    
    pB = b;
    
    printf("pB value is: %lf\n",*pB);
    
    /* print ut all elements of b[] */
    
    
    for(i=0;i<5;i++){
        printf("%.2lf ",*(pB+i));
    }
    
    printf("\n\nPrint a 2D array:\n");
    
    for(row=0;row<2;row++){
        for(col=0;col<4;col++){
            printf("%d ",c[row][col]);
        }
        printf("\n");
    }
    
    printf("\n");
    
    return 0;
}
