/* Pointers and arrays */
#include <stdio.h>

int main(){
    int i,j,row,col;
    int *pA;
    float *pB;
    int *pC;
    
    int a[]={5,10,15,20,25};
    
    pA = &a[0];
    
    printf("%d\n",a[0]);
    printf("%d %d\n",*pA, *(pA+1));
    printf("Second test\n");
    
    /* Another test */
    
    float b[] = {5.1,10.2,15.3,20.4,25.5};
      
    pB = b;
    
    printf("pB value is: %f\n",*pB);
    
    int c[3][4] = {
        {10,11,12,13},
        {14,15,16,17},
        {18,19,20,21}
    }; /* 2D array declaration and initialization*/
  
    pC = &c[0][0];
    
    printf("pC address= %d\n",*pC);
    
    /* print ut all elements of b[] */
    
    for(i=0;i<5;i++){
        printf("%.2lf ",*(pB+i));
    }
    
    printf("\n\nPrint a 2D array:\n");
    
    for(row=0;row<3;row++){
        for(col=0;col<4;col++){
            printf("%d ",c[row][col]);
        }
        printf("\n");
    }
    
    printf("\n");
    
    printf("Using pointers\n");
    int d=0;
    for(i=0;i<3;i++){
        for(j=0;j<4;j++){
            printf("%d ",*(pC+i+j+d));

        }
        d += 3;
        printf("\n");
    }
    return 0;
}
