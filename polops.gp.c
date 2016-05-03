/*-*- compile-command: "/usr/bin/gcc -c -o polops.gp.o -O3 -Wall -fno-strict-aliasing -fomit-frame-pointer -fPIC -I"/usr/local/include" polops.gp.c && /usr/bin/gcc -o polops.gp.so -shared -O3 -Wall -fno-strict-aliasing -fomit-frame-pointer -fPIC -Wl,-shared polops.gp.o -lc -lm -L/usr/local/lib -lpari"; -*-*/
#include <pari/pari.h>
#include <stdio.h>
#define PARI_STACK 2048000000
/*
GP;install("init_polops","v","init_polops","./polops.gp.so");
GP;install("pol2vec","D0,G,D0,G,D0,G,D0,G,","pol2vec","./polops.gp.so");
GP;install("vpolmult","D0,G,D0,G,","vpolmult","./polops.gp.so");
GP;install("vpoldz","D0,G,","vpoldz","./polops.gp.so");
GP;install("vpoldw","D0,G,","vpoldw","./polops.gp.so");
GP;install("diagmat","D0,G,","diagmat","./polops.gp.so");
GP;install("indcoef","D0,G,D0,G,D0,G,","indcoef","./polops.gp.so");
GP;install("lyapunov","D0,G,D0,G,","lyapunov","./polops.gp.so");
GP;install("firstlyapunov","D0,G,","firstlyapunov","./polops.gp.so");
GP;install("ferP","lD0,G,","ferP","./polops.gp.so");
GP;install("genfield","D0,G,D0,G,D0,G,D0,G,","genfield","./polops.gp.so");
GP;install("fer","v","fer","./polops.gp.so");
*/
void init_polops(void);
GEN pol2vec(GEN P, GEN n, GEN vx, GEN vy);
GEN vpolmult(GEN P, GEN Q);
GEN vpoldz(GEN P);
GEN vpoldw(GEN P);
GEN diagmat(GEN ord);
GEN indcoef(GEN deg, GEN H, GEN R);
GEN lyapunov(GEN N, GEN R);
GEN firstlyapunov(GEN R);
long ferP(GEN N);
GEN genfield(GEN m, GEN n, GEN k, GEN l);
void fer(void);
/*End of prototype*/

void
init_polops(void)	  /* void */
{
  return;
}

GEN
pol2vec(GEN P, GEN n, GEN vx, GEN vy)
{
  GEN aux = gen_0, p1;
  GEN p2;	  /* vec */
  p1 = gaddgs(n, 1);
  {
    long l3;
    p2 = cgetg(gtos(p1)+1, t_VEC);
    for (l3 = 1; gcmpsg(l3, p1) <= 0; ++l3)
      gel(p2, l3) = gen_0;
  }
  aux = p2;
  {
    GEN i;
    for (i = gen_1; gcmp(i, n) <= 0; i = gaddgs(i, 1))
    {
      gel(aux, gtos(i)) = polcoeff0(P, gtos(gaddgs(gsub(n, i), 1)), gvar(vx));
      if (gcmpgs(i, 1) > 0)
        gel(aux, gtos(i)) = polcoeff0(gel(aux, gtos(i)), gtos(i), gvar(vy));
    }
  }
  return aux;
}

/* multiplies 2 homogeneous polynomials as vectors */
GEN
vpolmult(GEN P, GEN Q)
{
  GEN len = gen_0, res = gen_0, aux = gen_0;
  GEN p1;	  /* vec */
  long l2;
  len = stoi((glength(P) + glength(Q)) - 1);
  {
    long l3;
    p1 = cgetg(gtos(len)+1, t_VEC);
    for (l3 = 1; gcmpsg(l3, len) <= 0; ++l3)
      gel(p1, l3) = gen_0;
  }
  res = p1;
  l2 = glength(P);
  {
    long j;
    for (j = 1; j <= l2; ++j)
    {
      GEN p4;	  /* vec */
      long l5;
      {
        long l6;
        p4 = cgetg(gtos(len)+1, t_VEC);
        for (l6 = 1; gcmpsg(l6, len) <= 0; ++l6)
          gel(p4, l6) = gen_0;
      }
      aux = p4;
      l5 = glength(Q);
      {
        long i;
        for (i = 1; i <= l5; ++i)
          gel(aux, (i + j) - 1) = gmul(gel(P, j), gel(Q, i));
      }
      res = gadd(res, aux);
    }
  }
  return res;
}

/* differentiates a homogeneous polynomial with resp. to z */
GEN
vpoldz(GEN P)
{
  GEN deg = gen_0, res = gen_0;
  GEN p1;	  /* vec */
  deg = stoi(glength(P) - 1);
  {
    long l2;
    p1 = cgetg(gtos(deg)+1, t_VEC);
    for (l2 = 1; gcmpsg(l2, deg) <= 0; ++l2)
      gel(p1, l2) = gen_0;
  }
  res = p1;
  {
    GEN i;
    for (i = gen_1; gcmp(i, deg) <= 0; i = gaddgs(i, 1))
      gel(res, gtos(i)) = gmul(gaddgs(gsub(deg, i), 1), gel(P, gtos(i)));
  }
  return res;
}

/* differentiates a homogeneous polynomial with resp. to w */
GEN
vpoldw(GEN P)
{
  GEN deg = gen_0, res = gen_0;
  GEN p1;	  /* vec */
  deg = stoi(glength(P) - 1);
  {
    long l2;
    p1 = cgetg(gtos(deg)+1, t_VEC);
    for (l2 = 1; gcmpsg(l2, deg) <= 0; ++l2)
      gel(p1, l2) = gen_0;
  }
  res = p1;
  {
    GEN i;
    for (i = gen_1; gcmp(i, deg) <= 0; i = gaddgs(i, 1))
      gel(res, gtos(gaddgs(gsub(deg, i), 1))) = gmul(gaddgs(gsub(deg, i), 1), gel(P, gtos(gaddgs(gsubsg(glength(P), i), 1))));
  }
  return res;
}

/* Create diagonal matrix of coefficients */
GEN
diagmat(GEN ord)
{
  GEN aux = gen_0, p1;
  GEN p2;	  /* vec */
  GEN p3;
  p1 = gaddgs(ord, 1);
  {
    long l4;
    p2 = cgetg(gtos(p1)+1, t_VEC);
    for (l4 = 1; gcmpsg(l4, p1) <= 0; ++l4)
      gel(p2, l4) = gen_0;
  }
  aux = p2;
  p3 = gaddgs(ord, 1);
  {
    GEN i;
    for (i = gen_1; gcmp(i, p3) <= 0; i = gaddgs(i, 1))
      gel(aux, gtos(i)) = gsub(ord, gmulsg(2, gsubgs(i, 1)));
  }
  return aux;
}

/* H and R are lists with the vectors needed to compute up to desired
* order ord (assume #H = #R)
*/
GEN
indcoef(GEN deg, GEN H, GEN R)
{
  GEN res = gen_0, aux = gen_0, p1;
  GEN p2;	  /* vec */
  long l3;
  p1 = gaddgs(deg, 1);
  {
    long l4;
    p2 = cgetg(gtos(p1)+1, t_VEC);
    for (l4 = 1; gcmpsg(l4, p1) <= 0; ++l4)
      gel(p2, l4) = gen_0;
  }
  res = p2;
  l3 = glength(H);
  {
    long i;
    for (i = 1; i <= l3; ++i)
    {
      long l5;
      l5 = glength(R);
      {
        long j;
        for (j = 1; j <= l5; ++j)
        {
          if (gequalsg((glength(gel(H, i)) + glength(gel(R, j))) - 3, deg))
          {
            long l6;
            GEN p7;	  /* vec */
            long l8;
            l6 = glength(gel(R, j));
            {
              long l9;
              p7 = cgetg(l6+1, t_VEC);
              for (l9 = 1; l9 <= l6; ++l9)
                gel(p7, l9) = gen_0;
            }
            /* Inverse order of conjugate vector */
            aux = p7;
            l8 = glength(aux);
            {
              long k;
              for (k = 1; k <= l8; ++k)
                gel(aux, k) = gconj(gel(gel(R, j), (glength(gel(R, j)) - k) + 1));
            }
            res = gadd(res, gadd(vpolmult(vpoldz(gel(H, i)), gel(R, j)), vpolmult(vpoldw(gel(H, i)), aux)));
          }
        }
      }
    }
  }
  return gmul(gen_I(), res);
}

/* Calcula les N primeres constants de Lyapunov del sistema donat */
GEN
lyapunov(GEN N, GEN R)
{
  GEN lastdg = gen_0, H = gen_0, L = gen_0;
  GEN p1, p2;	  /* vec */
  GEN p3, g = pol_x(fetch_user_var("g")), d = pol_x(fetch_user_var("d")), h = pol_x(fetch_user_var("h"));
  lastdg = gmulsg(2, gaddgs(N, 1));
  p1 = cgetg(2, t_VEC);
  p2 = cgetg(4, t_VEC);
  gel(p2, 1) = gen_0;
  gel(p2, 2) = gen_1;
  gel(p2, 3) = gen_0;
  gel(p1, 1) = p2;
  H = gtolist(p1);
  L = listcreate();
  p3 = gsubgs(lastdg, 1);
  {
    GEN i;
    for (i = stoi(3); gcmp(i, p3) <= 0; i = gaddgs(i, 2))
    {
      GEN p4;
      GEN p5;	  /* vec */
      GEN p6, p7;
      GEN p8;	  /* vec */
      GEN p9;
      pari_printf("%Ps\n", i);
      /* Part senar */
      g = indcoef(i, H, R);
      d = diagmat(i);
      p4 = gaddgs(i, 1);
      {
        long l10;
        p5 = cgetg(gtos(p4)+1, t_VEC);
        for (l10 = 1; gcmpsg(l10, p4) <= 0; ++l10)
          gel(p5, l10) = gen_0;
      }
      h = p5;
      p6 = gaddgs(i, 1);
      {
        GEN j;
        for (j = gen_1; gcmp(j, p6) <= 0; j = gaddgs(j, 1))
          gel(h, gtos(j)) = gdiv(gel(g, gtos(j)), gel(d, gtos(j)));
      }
      listput(H, h, 0);
      /* Part parella */
      g = indcoef(gaddgs(i, 1), H, R);
      d = diagmat(gaddgs(i, 1));
      p7 = gaddgs(i, 2);
      {
        long l11;
        p8 = cgetg(gtos(p7)+1, t_VEC);
        for (l11 = 1; gcmpsg(l11, p7) <= 0; ++l11)
          gel(p8, l11) = gen_0;
      }
      h = p8;
      listput(L, gdiv(gel(g, gtos(gaddgs(gdivgs(gaddgs(i, 1), 2), 1))), gen_I()), 0);
      p9 = gaddgs(i, 2);
      {
        GEN j;
        for (j = gen_1; gcmp(j, p9) <= 0; j = gaddgs(j, 1))
        {
          if (!gequalgs(gel(d, gtos(j)), 0))
            gel(h, gtos(j)) = gdiv(gel(g, gtos(j)), gel(d, gtos(j)));
        }
      }
      listput(H, h, 0);
    }
  }
  return L;
}

/* Calcula la primera constant de Lyapunov no nul·la i la retorna,
nomes busca fins la constant N
*/
GEN
firstlyapunov(GEN R)
{
  GEN lastdg = gen_0, H = gen_0, N = gen_0, g = gen_0, d = gen_0, h = gen_0, k = gen_0;
  long l1;
  GEN maxL = pol_x(fetch_user_var("maxL"));
  GEN p2, p3;	  /* vec */
  GEN p4;
  N = gen_0;
  k = gen_0;
  l1 = glength(R);
  {
    long i;
    for (i = 1; i <= l1; ++i)
      N = gmaxgs(N, glength(gel(R, i)));
  }
  maxL = gsubgs(gadd(gmul(N, N), gmulsg(3, N)), 7);
  lastdg = gmulsg(2, gaddgs(maxL, 1));
  p2 = cgetg(2, t_VEC);
  p3 = cgetg(4, t_VEC);
  gel(p3, 1) = gen_0;
  gel(p3, 2) = gen_1;
  gel(p3, 3) = gen_0;
  gel(p2, 1) = p3;
  H = gtolist(p2);
  p4 = gsubgs(lastdg, 1);
  {
    GEN i;
    for (i = stoi(3); gcmp(i, p4) <= 0; i = gaddgs(i, 2))
    {
      GEN p5;
      GEN p6;	  /* vec */
      GEN p7, p8;
      GEN p9;	  /* vec */
      GEN p10;
      /* Part senar */
      g = indcoef(i, H, R);
      d = diagmat(i);
      p5 = gaddgs(i, 1);
      {
        long l11;
        p6 = cgetg(gtos(p5)+1, t_VEC);
        for (l11 = 1; gcmpsg(l11, p5) <= 0; ++l11)
          gel(p6, l11) = gen_0;
      }
      h = p6;
      p7 = gaddgs(i, 1);
      {
        GEN j;
        for (j = gen_1; gcmp(j, p7) <= 0; j = gaddgs(j, 1))
          gel(h, gtos(j)) = gdiv(gel(g, gtos(j)), gel(d, gtos(j)));
      }
      listput(H, h, 0);
      /* Part parella */
      k = gaddgs(k, 1);
      g = indcoef(gaddgs(i, 1), H, R);
      if (!gequalgs(gel(g, gtos(gaddgs(gdivgs(gaddgs(i, 1), 2), 1))), 0))
      {
        GEN p12;	  /* vec */
        p12 = cgetg(3, t_VEC);
        gel(p12, 1) = gcopy(k);
        gel(p12, 2) = gdiv(gel(g, gtos(gaddgs(gdivgs(gaddgs(i, 1), 2), 1))), gen_I());
        return p12;
      }
      d = diagmat(gaddgs(i, 1));
      p8 = gaddgs(i, 2);
      {
        long l13;
        p9 = cgetg(gtos(p8)+1, t_VEC);
        for (l13 = 1; gcmpsg(l13, p8) <= 0; ++l13)
          gel(p9, l13) = gen_0;
      }
      h = p9;
      p10 = gaddgs(i, 2);
      {
        GEN j;
        for (j = gen_1; gcmp(j, p10) <= 0; j = gaddgs(j, 1))
        {
          if (!gequalgs(gel(d, gtos(j)), 0))
            gel(h, gtos(j)) = gdiv(gel(g, gtos(j)), gel(d, gtos(j)));
        }
      }
      listput(H, h, 0);
    }
  }
  return strtoGENstr("Centre");
}

long
ferP(GEN N)
{
  GEN r1 = gen_0, r2 = gen_0, R = gen_0;
  GEN p1;	  /* vec */
  GEN p2;
  GEN p3, p4;	  /* vec */
  {
    long l5;
    p1 = cgetg(gtos(N)+1, t_VEC);
    for (l5 = 1; gcmpsg(l5, N) <= 0; ++l5)
      gel(p1, l5) = gen_0;
  }
  /*,NN*/r1 = p1;
  gel(r1, gtos(N)) = gen_1;
  p2 = gaddgs(N, 1);
  {
    long l6;
    p3 = cgetg(gtos(p2)+1, t_VEC);
    for (l6 = 1; gcmpsg(l6, p2) <= 0; ++l6)
      gel(p3, l6) = gen_0;
  }
  r2 = p3;
  gel(r2, 1) = gen_1;
  p4 = cgetg(3, t_VEC);
  gel(p4, 1) = gcopy(r1);
  gel(p4, 2) = gcopy(r2);
  R = gtolist(p4);
  /*NN=(N-1)*(N-1);*/
  gettime();
  pari_printf("%Ps\n", firstlyapunov(R));
  return gettime();
}

/* Generar pol z^m*w^n+z^k*w^l en notacio vectorial */
GEN
genfield(GEN m, GEN n, GEN k, GEN l)	  /* list */
{
  GEN v1 = gen_0, v2 = gen_0, p1;
  GEN p2;	  /* vec */
  GEN a1 = pol_x(fetch_user_var("a1")), b1 = pol_x(fetch_user_var("b1")), p3;
  GEN p4;	  /* vec */
  GEN a2 = pol_x(fetch_user_var("a2")), b2 = pol_x(fetch_user_var("b2"));
  GEN p5;	  /* vec */
  p1 = gaddgs(gadd(m, n), 1);
  {
    long l6;
    p2 = cgetg(gtos(p1)+1, t_VEC);
    for (l6 = 1; gcmpsg(l6, p1) <= 0; ++l6)
      gel(p2, l6) = gen_0;
  }
  v1 = p2;
  gel(v1, gtos(gaddgs(n, 1))) = gadd(a1, gmul(gen_I(), b1));
  p3 = gaddgs(gadd(k, l), 1);
  {
    long l7;
    p4 = cgetg(gtos(p3)+1, t_VEC);
    for (l7 = 1; gcmpsg(l7, p3) <= 0; ++l7)
      gel(p4, l7) = gen_0;
  }
  v2 = p4;
  gel(v2, gtos(gaddgs(l, 1))) = gadd(a2, gmul(gen_I(), b2));
  p5 = cgetg(3, t_VEC);
  gel(p5, 1) = gcopy(v1);
  gel(p5, 2) = gcopy(v2);
  return gtolist(p5);
}

void
fer(void)	  /* void */
{
  {
    long m;
    for (m = 1; m <= 3; ++m)
    {
      long n;
      for (n = 1; n <= 3; ++n)
      {
        long k;
        for (k = 1; k <= 3; ++k)
        {
          long l;
          for (l = 1; l <= 3; ++l)
          {
            pari_printf("%ld,%ld,%ld,%ld\n", m, n, k, l);
            pari_printf("%Ps\n", firstlyapunov(genfield(stoi(m), stoi(n), stoi(k), stoi(l))));
          }
        }
      }
    }
  }
  return;
}

int main(int argc, char const *argv[]) {
    
    int taskId;
    GEN m,n,k,l;

    /* taskId m n k l */
    if (argc != 6 || sscanf(argv[1],"%d",&taskId)!=1
       ) {
        fprintf(stderr,"%s:: taskId m n k l\n",argv[0]);
        return -1;
    }
    m = gp_read_str(argv[2]);
    n = gp_read_str(argv[3]);
    k = gp_read_str(argv[4]);
    l = gp_read_str(argv[5]);

    /* Quanta memoria assignar? */
    pari_init(PARI_STACK,2);
    init_polops();
    firstlyapunov(genfield(m,n,k,l));
    pari_close();

    return 0;
}