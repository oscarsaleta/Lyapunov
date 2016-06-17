#include <stdio.h>
#include <stdlib.h>

int max(int a, int b) {
    return a>b?a:b;
}

int main(int argc, char *argv[]) {

    if (argc<2) {
        fprintf(stderr,"%s grau\n",argv[0]);
        return -1;
    }

    int k,l,m,n,p,q;
    int maxdeg=atoi(argv[1]);
    int count=0;

    for (k=0;k<=maxdeg;k++) {
        for (l=0;l<=maxdeg;l++) {
            if (k+l > maxdeg /* grau major que el desitjat */
                    || (k==0 && l==0) /* graus=0 es una constant */
                    || k+l < 2) /* ha de ser grau>2 perque i*z es grau 1 */
                continue;

            for (m=0;m<=maxdeg;m++) {
                for (n=0;n<=maxdeg;n++) {
                    if (m+n > maxdeg /* grau major que el desitjat */
                            || (m==0 && n==0) /* graus=0 es una constant */
                            || m+n < 2) /* ha de ser grau>2 perque i*z es grau 1 */
                        continue;
                    if (k+l <= m+n)
                        continue; /* monomis repetits */
                    else if (k==m && l==n)
                        continue; /* mateix monomi 2 cops */
                    else if (!(k==n && k==2 && l==m && l==0)/* agafem centres quadratics Darboux */
                            && !(l==n && l==0)/* agafem centres holomorfs */
                            && !(k==m && l!=n && (k-l-1)!=0))/* agafem centres Hamiltonians o new Darboux */
                        continue; 

                    for (p=0;p<=maxdeg;p++) {
                        for (q=0;q<=maxdeg;q++) {
                            if (p+q > maxdeg /* grau major que el desitjat */
                                    || (p==0 && q==0) /* graus=0 es una constant */
                                    || p+q < 2) /* ha de ser grau>2 perque i*z es grau 1 */
                                continue;

                            if (max(max(k+l,m+n),p+q)<maxdeg)
                                continue; /* només volem grau maxdeg */
                            if (k==m && m==p && (k-l-1!=0))
                                continue; /* hamiltonian/new darboux amb 3 monomis */
                            if (l==n && n==q)
                                continue; /* holomorfic */

                            if (k-l-1!=0)
                                fprintf(stdout,"%d,%d,%d,%d,%d,%d,%d,%d\n",count,count,k,l,m,n,p,q);
                            else
                                fprintf(stdout,"%d,%d,%d,%d,%d,%d,%d,%d\n",count,count,m,n,k,l,p,q);
                            count++;
                        }
                    }
                }
            }
        }
    }
    return 0;
}

