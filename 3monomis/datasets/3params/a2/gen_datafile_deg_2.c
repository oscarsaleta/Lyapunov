#include <stdio.h>
#include <stdlib.h>

int max(int a, int b) {
    return(a>b?a:b);
}

int main (int argc, char *argv[]) {
    int k,l,m,n,p,q;
    int maxdeg=atoi(argv[1]);
    int count=0;

    for (k=0;k<=maxdeg;k++) {
        for (l=0;l<=maxdeg;l++) {
            if (k+l > maxdeg || (k==0 && l==0) || k+l<2)
                continue;
            for (m=0;m<=maxdeg;m++) {
                for (n=0;n<=maxdeg;n++) {
                    if (m+n > maxdeg || (m==k && n==l) || (m==0 && n==0) || m+n<2)
                        continue;
                    for (p=0;p<=maxdeg;p++) {
                        for (q=0;q<=maxdeg;q++) {
                            if (p+q > maxdeg || (p==m && q==n) || (p==k && q==l)
                                    || (p==0 && q==0) || p+q<2)
                                continue;
                            if (max(max(k+l,m+n),p+q)<maxdeg)
                                continue;
                            if (k==m && m==p && (k-l-1!=0))
                                continue;
                            if (l==n && n==q)
                                continue;
                            if (k-l-1==0)
                                continue;
                            // mirar nous centres monomis de 2 en 2
                            if (p-q-1==0 && k==m) {
                                printf("%d,%d,%d,%d,%d,%d\n",count,count,k,l,p,q);
                                printf("%d,%d,%d,%d,%d,%d\n",count+1,count+1,m,n,p,q);
                            } else
                                continue;
                            count+=2;
                        }
                    }
                }
            }
        }
    }

    return 0;
}
