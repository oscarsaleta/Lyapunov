
pol2vec(P,n,vx,vy)=
{
    local(i,aux);
    aux = vector(n+1);
    for(i=1,n,
        aux[i]=polcoeff(P,n-i+1,vx);
        if(i>1,
            aux[i]=polcoeff(aux[i],i,vy);
        );
    );
    return(res);
}


/* multiplies 2 homogeneous polynomials as vectors */
vpolmult(P,Q)=
{
    local(i,j,len,res,aux);
    len = #P + #Q - 1;
    res = vector(len);

    for(j=1,#P,
        aux = vector(len);
        for (i=1,#Q,
            aux[i+j-1] = P[j]*Q[i];
            );
        res += aux;
       );
    res
};

/* differentiates a homogeneous polynomial with resp. to z */
vpoldz(P)=
{
    local(i,deg,res);
    deg = #P-1;
    res = vector(deg);

    for(i=1,deg,
        res[i] = (deg-i+1)*P[i];
       );
    res
};

/* differentiates a homogeneous polynomial with resp. to w */
vpoldw(P)=
{
    local(i,deg,res);
    deg = #P-1;
    res = vector(deg);

    for(i=1,deg,
        res[deg-i+1] = (deg-i+1)*P[#P-i+1];
       );
    res
};

/* Create diagonal matrix of coefficients */
diagmat(ord)=
{
    local(i,aux);
    aux = vector(ord+1);
    for(i=1,ord+1,aux[i]=ord-2*(i-1););
    aux
};

/* H and R are lists with the vectors needed to compute up to desired
 * order ord (assume #H = #R)
 */
indcoef(deg,H,R)=
{
    local(i,j,k,res);
    res = vector(deg+1);

    for(i=1,#H,
        for(j=1,#R,
            if(#H[i]+#R[j]-3==deg,
                /* Inverse order of conjugate vector */
                aux = vector(#R[j]);
                for(k=1,#aux,
                    aux[k] = conj(R[j][#R[j]-k+1]);
                );
                res += vpolmult(vpoldz(H[i]),R[j])+vpolmult(vpoldw(H[i]),aux);
            );
        );
    );
    I*res
};


/* Calcula les N primeres constants de Lyapunov del sistema donat */
lyapunov(N,R)=
{
    local(lastdg,H,L,i,j,k,kmax);
    lastdg = 2*(N+1);
    H=List([[0,1,0]]);
    L=List();
    forstep(i=3,lastdg-1,2,
            print(i);
        /* Part senar */
        g=indcoef(i,H,R);
        d=diagmat(i);
        h=vector(i+1);
        for(j=1,i+1,
            h[j]=g[j]/d[j];
        );
        listput(H,h);
        /* Part parella */
        g=indcoef(i+1,H,R);
        d=diagmat(i+1);
        h=vector(i+2);
        listput(L,g[((i+1)/2)+1]/I);
        for(j=1,i+2,
            if(d[j]!=0,
                h[j]=g[j]/d[j];
            );
        );
        listput(H,h);
    );
    return(L);
};

/* Calcula la primera constant de Lyapunov no nul·la i la retorna,
   nomes busca fins la constant N
 */
firstlyapunov(R)=
{
    local(lastdg,H,L,i,j,k,kmax,N,g,d,h);
    N=0;
    k=0;
    for(i=1,#R,
        N = max(N,#R[i]);
    );
    maxL = N*N+3*N-7;
    lastdg = 2*(maxL+1);
    H=List([[0,1,0]]);
    forstep(i=3,lastdg-1,2,
        /* Part senar */
        g=indcoef(i,H,R);
        d=diagmat(i);
        h=vector(i+1);
        for (j=1,i+1,
            h[j]=g[j]/d[j];
        );
        listput(H,h);
        /* Part parella */
        k++;
        g=indcoef(i+1,H,R);
        if(g[((i+1)/2)+1]!=0,
            return([k,g[((i+1)/2)+1]/I]);
        );
        d=diagmat(i+1);
        h=vector(i+2);
        for(j=1,i+2,
            if(d[j]!=0,
                h[j]=g[j]/d[j];
            );
        );
        listput(H,h);
    );
    return("Centre");
};

ferP(N)=
{
    local(r1,r2,R,NN);
    r1=vector(N);r1[N]=1;
    r2=vector(N+1);r2[1]=1;
    R=List([r1,r2]);
    NN=(N-1)*(N-1);
    gettime();
    print(firstlyapunov(R));
    gettime()
};

