#include <stdio.h>

int main() {
    int count=0;
    int k,l,m,n;

    for (k=0;k<25;k++) {
        for (l=0;l<25;l++) {
            for (m=0;m<25;m++) {
                for (n=0;n<25;n++) {
                    if (k+l<2 || n+m<2)
                        continue;
                    if (k+l <= m+n)
                        continue;
                    if (k==m && l==n)
                        continue;
                    if (k==n && k==2 && l==m && l==0)
                        continue;
                    if (l==n && l==0)
                        continue;
                    if (k==m && l!=n && (k-l-1)!=0)
                        continue;
                    if (k-l-1==0) {
                        fprintf(stdout,"%d,%d,%d,%d,%d\n",count,m,n,k,l);
                        count++;
                        continue;
                    }
                    fprintf(stdout,"%d,%d,%d,%d,%d\n",count,k,l,m,n);
                    count++;
                }
            }
        }
    }
    return 0;
}
