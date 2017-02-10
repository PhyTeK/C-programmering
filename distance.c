#include <stdio.h>
#include <math.h>

struct point{
    double x,y;
};


int main(){
    
    struct point a,b;
    double d;
    
    scanf("%lf%lf%lf%lf",&a.x,&a.y,&b.x,&b.y);
    
    d=sqrt(pow(a.x - b.x,2) + pow(a.y - b.y,2));
    
    printf("%lf\n",d);
    
    return 0;
}
