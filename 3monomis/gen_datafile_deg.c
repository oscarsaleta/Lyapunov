#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[]) {
    int k,l,m,n,p,q;
    int maxdeg=atoi(argv[1]);
    int count=0;

    for (k=0;k<maxdeg;k++) {
        for (l=0;l<maxdeg;l++) {
            if (k+l > maxdeg)
                continue;
            for (m=0;m<maxdeg;m++) {
                for (n=0;n<maxdeg;n++) {
                    if (m+n > maxdeg || (m==l && n==k))
                        continue;
                    for (p=0;p<maxdeg;p++) {
                        for (q=0;q<maxdeg;q++) {
                            if (p+q > maxdeg || (p==m && q==n) || (p==k && q==l))
                                continue;
                            printf("%d,%d,%d,%d,%d,%d,%d\n",count,k,l,m,n,p,q);
                            count++;
                        }
                    }
                }
            }
        }
    }

    return 0;
}
