#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
    int l,n,p;
    int count=0;

    for (l=2;l<50;l++) {
        for (n=2;n<50;n++) {
            if (l==n)
                continue;
            for (p=2;p<50;p++) {
                printf("%d,%d,%d,%d\n",count,l,n,p);
                count++;
            }
        }
    }
    return 0;

}
