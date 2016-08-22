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
            if (k+l > maxdeg || k+l<2)
                continue;
            for (m=0;m<=maxdeg;m++) {
                for (n=0;n<=maxdeg;n++) {
                    if (m+n > maxdeg || (m==k && n==l) || m+n<2)
                        continue;
                    for (p=0;p<=maxdeg;p++) {
                        for (q=0;q<=maxdeg;q++) {
                            if (p+q > maxdeg || (p==m && q==n) || (p==k && q==l) || p+q<2)
                                continue;
                            if (max(max(k+l,m+n),p+q)<maxdeg)
                                continue;
                            if (k==q && q==2 && m==n && n==1 && l==p && p==0)
                                continue;
                            if (l==n && n==q && q==0)
                                continue;
                            if (k==m && m==p && l!=n && l!=q && n!=q && ( (k-l-1)!=0 || (m-n-1)!=0 || (p-q-1)!=0) )
                                continue;
                            if (k-l-1!=0)
                                printf("%d,%d,%d,%d,%d,%d,%d,%d\n",count,count,k,l,m,n,p,q);
                            else
                                continue;
                            count++;
                        }
                    }
                }
            }
        }
    }

    return 0;
}
