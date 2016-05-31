#include <stdio.h>

int main() {
    int count=0;
    int k,l,m,n;

    for (k=0;k<20;k++) {
        for (l=0;l<20;l++) {
            for (m=0;m<k;m++) {
                for (n=0;n<l;n++) {
                    if (k+l <= m+n)
                        continue;
                    if (k==n && k==2 && l==m && l==0)
                        continue;
                    if (l==n && l==0)
                        continue;
                    if (k==m && l==n)
                        continue;
                    if ((k-l-1)*(m-n-1)==0)
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
