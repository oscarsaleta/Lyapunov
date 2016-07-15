#include <stdio.h>

int main() {
    int count=0;
    int k,l,m,n;

    for (k=0;k<25;k++) {
        for (l=0;l<25;l++) {
            for (m=0;m<25;m++) {
                for (n=0;n<25;n++) {
                    if (k+l<2 || n+m<2) {
                        fprintf(stderr,"%d,%d,%d,%d - Not valid (R has to be quadratic at least)\n",k,l,m,n);
                        continue;
                    } else if (k+l <= m+n) {
                        fprintf(stderr,"%d,%d,%d,%d - Repeated (k+l <= m+n)\n",k,l,m,n);
                        continue;
                    } else if (k==m && l==n) {
                        fprintf(stderr,"%d,%d,%d,%d - Same monomial twice\n",k,l,m,n);
                        continue;
                    } else if (k==n && k==2 && l==m && l==0) {
                        fprintf(stderr,"%d,%d,%d,%d - Center: (a) quadratic Darboux\n",k,l,m,n);
                        continue;
                    } if (l==n && l==0) {
                        fprintf(stderr,"%d,%d,%d,%d - Center: (b) holomorphic\n",k,l,m,n);
                        continue;
                    } else if (k==m && l!=n && (k-l-1)!=0) {
                        fprintf(stderr,"%d,%d,%d,%d - Center: (d) Hamiltonian or new Darboux\n",k,l,m,n);
                        continue;
                    }
                    if (k-l-1==0) {
                        fprintf(stderr,"%d,%d,%d,%d - Center: (c) reversible\n",k,l,m,n);
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
