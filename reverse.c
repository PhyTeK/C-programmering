#include <stdio.h>
#include <string.h>

int main(int argc,char *argv[]){
    
    char str[20];
    int i,l;
    scanf("%s",&str[0]);
    l=strlen(str);

    for(i=l;i>=0;i--){
        printf("%c",str[i]);
        
    }
    printf("\n");
    
    return 0;
    
}
