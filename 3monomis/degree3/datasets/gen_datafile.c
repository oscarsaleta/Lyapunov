#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
    int d=atoi(argv[1]);
    int k,l;
    int m,n;
    int p,q;
    int tasknum=0;

    for(k=0;k<=d;k++) {
        for(l=0;l<=d;l++) {
            for(m=0;m<=d;m++) {
                for(n=0;n<=d;n++) {
                    for(p=0;p<=d;p++) {
                        for(q=0;q<=d;q++) {
                           if (k+l<2 || m+n<2 || p+q<2)
                               continue;
                           if (k+l>d || m+n>d || p+q>d)
                               continue;
                           if (k+l!=d && m+n!=d && p+q!=d)
                               continue;
                           if (m+n<p+q) // per no repetir casos
                               continue;
                           if ( (k==m && l==n) || (k==p && l==q) || (m==p && n==q) )
                               continue; // no serien 3 monomis
                           if (k-l-1==0)
                               continue;
                           tasknum++;
                           printf("%d,%d,%d,%d,%d,%d,%d,%d\n",tasknum,tasknum,k,l,m,n,p,q);
                        }
                    }
                }
            }
        }
    }
    
    return 0;
}
