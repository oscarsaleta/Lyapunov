#include <stdio.h>
#include <stdlib.h>

int max(int a, int b) {
    return(a>b?a:b);
}

int main (int argc, char *argv[]) {
    int k,l,m,n,p,q;
    int maxdeg=atoi(argv[1]);
    int homoflag=atoi(argv[2]);
    int count=0;

    for (k=0;k<=maxdeg;k++) {
        for (l=0;l<=maxdeg;l++) {
            if (k+l > maxdeg || k+l<2) //eliminar graus majors del que volem o <2
                continue;
            for (m=0;m<=maxdeg;m++) {
                for (n=0;n<=maxdeg;n++) {
                    if (m+n > maxdeg || m+n<2) //eliminar graus majors del que volem o <2
                        continue;
                    for (p=0;p<=maxdeg;p++) {
                        for (q=0;q<=maxdeg;q++) {
                            if (p+q > maxdeg || p+q<2)
                                continue;
                            if (homoflag && max(max(k+l,m+n),p+q)<maxdeg) //forçar homogenis?
                                continue;
                            else
                                printf("%d,%d,%d,%d,%d,%d,%d,%d\n",count,count,k,l,m,n,p,q);
                            count++;
                        }
                    }
                }
            }
        }
    }

    return 0;
}
