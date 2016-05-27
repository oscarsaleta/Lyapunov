#include <stdio.h>

int main() {
    int count=0;
    int k,l,m,n;

    for (k=0;k<20;k++) {
        for (l=0;l<20;l++) {
            for (m=0;m<20;m++) {
                for (n=0;n<20;n++) {
                    if (k+l <= m+n)
                        continue;
                    if (k==m && l==n)
                        continue;
                    if (k-l-1==0 || m-n-1==0)
                        continue;
                    fprintf(stdout,"%d,%d,%d,%d,%d\n",count,k,l,m,n);
                    count++;
                }
            }
        }
    }
    return 0;
}
