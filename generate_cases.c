#include <stdio.h>
#include <stdlib.h>

inline int max(int a, int b) { return (a > b ? a : b); }

// taken from http://www.netlib.no/netlib/toms/382
int twiddle(int *x, int *y, int *z, int *p) {
    register int i, j, k;
    j = 1;
    while (p[j] <= 0)
        j++;
    if (p[j - 1] == 0) {
        for (i = j - 1; i != 1; i--)
            p[i] = -1;
        p[j] = 0;
        *x = *z = 0;
        p[1] = 1;
        *y = j - 1;
    } else {
        if (j > 1)
            p[j - 1] = 0;
        do
            j++;
        while (p[j] > 0);
        k = j - 1;
        i = j;
        while (p[i] == 0)
            p[i++] = -1;
        if (p[i] == -1) {
            p[i] = p[k];
            *z = p[k] - 1;
            *x = i - 1;
            *y = k - 1;
            p[k] = -1;
        } else {
            if (i == p[0])
                return (1);
            else {
                p[j] = p[i];
                *z = p[i] - 1;
                p[i] = 0;
                *x = j - 1;
                *y = i - 1;
            }
        }
    }
    return (0);
}

void inittwiddle(int m, int n, int *p) {
    int i;
    p[0] = n + 1;
    for (i = 1; i != n - m + 1; i++)
        p[i] = 0;
    while (i != n + 1) {
        p[i] = i + m - n;
        i++;
    }
    p[n + 1] = -2;
    if (m == 0)
        p[1] = 1;
}

int main(int argc, char *argv[]) {
    int k, l, m, n, p, q;       // exponents
    int taskid = 0;             // case counter
    int maxdeg;                 // max degree
    int count;                  // number of (m,n) cases
    int *exp_z, *exp_w;         // vectors of exponents for z and w
    int x, y, z, *permutations; // variables for twiddle
    int i;                      // for loops
    int *indexes;               // vector [1,2,...,count]
    int selection[2];           // pair of indexes that form a combination

    if (argc != 2 || sscanf(argv[1], "%d", &maxdeg) != 1) {
        fprintf(stderr, "%s: max_degree\n", argv[0]);
        return 1;
    }

    for (k = 0; k <= maxdeg; k++) {
        for (l = 0; l <= maxdeg; l++) {
            // discard incorrect degrees or resonant first monomial
            if (k + l < 2 || k + l > maxdeg || k - l - 1 == 0)
                continue;
            // count how many cases there are for other two monomials
            // if maxdeg=d, there are d^2 combinations of exponents
            // -1: when (m,n)==(k,l)
            // -3: (0,0),(1,0),(0,1) of degree <2
            // last term: counts the cases where m+n>d, it's the sum
            // 1+2+...+(maxdeg-1)+maxdeg
            count = (maxdeg + 1) * (maxdeg + 1) - 1 - 3 -
                    (maxdeg) * (maxdeg + 1) / 2;
            // allocate vector with every possible pair of exponents
            exp_z = malloc(count * sizeof(int));
            exp_w = malloc(count * sizeof(int));
            // fill exponent pairs vector
            i = 0;
            for (m = 0; m <= maxdeg; m++) {
                for (n = 0; n <= maxdeg; n++) {
                    if ((m == k && n == l) || m + n < 2 || m + n > maxdeg)
                        continue;
                    exp_z[i] = m;
                    exp_w[i] = n;
                    i++;
                }
            }
            // define vector of indexes
            indexes = malloc(count * sizeof(int));
            for (i = 0; i < count; i++)
                indexes[i] = i;
            // allocate vector of permutations
            permutations = malloc((count + 2) * sizeof(int));
            inittwiddle(2, count, permutations);
            selection[0] = count - 2;
            selection[1] = count - 1;
            m = exp_z[selection[0]];
            n = exp_w[selection[0]];
            p = exp_z[selection[1]];
            q = exp_w[selection[1]];
            taskid++;
            printf("%d,%d,%d,%d,%d,%d,%d,%d\n", taskid, taskid, k, l, m, n, p,
                   q);
            while (!twiddle(&x, &y, &z, permutations)) {
                selection[z] = indexes[x];
                m = exp_z[selection[0]];
                n = exp_w[selection[0]];
                p = exp_z[selection[1]];
                q = exp_w[selection[1]];
                taskid++;
                printf("%d,%d,%d,%d,%d,%d,%d,%d\n", taskid, taskid, k, l, m, n,
                       p, q);
            }
            free(permutations);
            free(indexes);
            free(exp_z);
            free(exp_w);
        }
    }
    return 0;
}
