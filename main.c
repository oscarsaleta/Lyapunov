#include <stdio.h>
#include <pari/pari.h>

#define PARI_STACK 1000000

int main(int argc, char const *argv[]) {
    
    int taskId;
    int m,n,k,l;

    /* taskId m n k l */
    if (argc != 6
        || sscanf(argv[1],"%d",&taskId)!=1
        || sscanf(argv[2],"%d",&m)!=1
        || sscanf(argv[3],"%d",&n)!=1
        || sscanf(argv[4],"%d",&k)!=1
        || sscanf(argv[5],"%d",&l)!=1
       ) {
        fprintf(stderr,"%s:: taskId m n k l\n",argv[0]);
        return -1;
    }

    /* Quanta memoria assignar? */
    pari_init(PARI_STACK,2);
    init_polops();
    fer();
    pari_close();

    return 0;
}
