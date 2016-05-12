#include <stdio.h>

int main () {
    int k,n;
    for (n=3,k=0;n<101;n++,k++) {
        fprintf(stdout,"%d,0,%d,%d,0\n",k,n-1,n);
    }
    return 0;
}
