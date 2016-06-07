#include <stdio.h>

int main() {
    int k,l;
    int m,n;
    int p,q;
    int tasknum=0;

    for(k=0;k<20;k++) {
        for(l=0;l<20;l++) {
            for(m=0;m<20;m++) {
                for(n=0;n<20;n++) {
                    for(p=0;p<20;p++) {
                        for(q=0;q<20;q++) {
                           if (k+l<=m+n || m+n<=p+q) // per no repetir casos
                               continue;
                           if (l==n && n==q && q==0)
                               continue; //holomorfs
                           if (k==q && q==m+n && m==n && l==p && p==0)
                               continue; //darboux?
                           if (k==m && m==p && l!=n && l!=q && n!=q && k-l-1!=0)
                               continue; //new darboux?
                        }
                    }
                }
            }
        }
    }
    
    return 0;
}
