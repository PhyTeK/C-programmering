/*
 * arraytest.c
 * 
 * print out an array using a pointer
 */


#include <stdio.h>

int main()
{
	int arr[5]={5,10,15,20,25};
    int *parr;
    int i;
    
    parr=&arr[0];
    
    for(i=0;i<5;i++){
        printf("%d ",*(parr+i));
    }
    printf("\n");
	return 0;
}

