
/*pol2vec(P,n,vx,vy)=
{
    local(aux);
    aux = vector(n+1);
    for(i=1,n,
        aux[i]=polcoeff(P,n-i+1,vx);
        if(i>1,
            aux[i]=polcoeff(aux[i],i,vy);
        );
    );
    return(aux);
}*/


/* multiplies 2 homogeneous polynomials as vectors */
vpolmult(P:vec,Q:vec)=
{
    local(len:int,res:vec,aux:vec);
    len = #P + #Q - 1;
    res = vector(len);

    for(j=1,#P,
        aux = vector(len);
        for (i=1,#Q,
            aux[i+j-1] = P[j]*Q[i];
            );
        res += aux;
       );
    return(res)
};

/* differentiates a homogeneous polynomial with resp. to z */
vpoldz(P)=
{
    local(deg:int,res:vec);
    deg = #P-1;
    res = vector(deg);

    for(i=1,deg,
        res[i] = (deg-i+1)*P[i];
       );
    return(res)
};

/* differentiates a homogeneous polynomial with resp. to w */
vpoldw(P:vec)=
{
    local(deg:int,res:vec);
    deg = #P-1;
    res = vector(deg);

    for(i=1,deg,
        res[deg-i+1] = (deg-i+1)*P[#P-i+1];
       );
    return(res)
};

/* Create diagonal matrix of coefficients */
diagmat(ord:int)=
{
    local(aux:vec);
    aux = vector(ord+1);
    for(i=1,ord+1,aux[i]=ord-2*(i-1););
    return(aux)
};

/* H and R are lists with the vectors needed to compute up to desired
 * order ord (assume #H = #R)
 */
indcoef(deg:int,H:vec,R:vec)=
{
    local(res:vec,aux:vec);
    res = vector(deg+1);

    for(i=1,#H,
        if(H[i]==0,next);
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
    return(I*res)
};


/* Calcula les N primeres constants de Lyapunov del sistema donat */
lyapunov(N:int,R:vec)=
{
    local(lastdg:int,H:vec,L:vec,g:vec,d:vec,h:vec);
    lastdg = 2*(N+1);
    H=vector(lastdg);
    L=vector(N);
    H[1]=[0,1,0];
    forstep(i=3,lastdg-1,2,
        /*print(i);*/
        /* Part senar */
        g=indcoef(i,H,R);
        d=diagmat(i);
        h=vector(i+1);
        for(j=1,i+1,
            h[j]=g[j]/d[j];
        );
        H[i-1]=h;
        /* Part parella */
        g=indcoef(i+1,H,R);
        d=diagmat(i+1);
        h=vector(i+2);
        L[(i-1)/2]=g[((i+1)/2)+1]/I;
        for(j=1,i+2,
            if(d[j]!=0,
                h[j]=g[j]/d[j];
            );
        );
        H[i]=h;
    );
    return(L);
};

/* Calcula la primera constant de Lyapunov no nul·la i la retorna,
   nomes busca fins la constant N
 */
firstlyapunov(R:vec)=
{
    local(lastdg:int,H:vec,N:int,g:vec,d:vec,h:vec,k:int,maxL:int);
    N=0;
    k=0;
    for(i=1,#R,
        N = max(N,#R[i]);
    );
    maxL = N*N+3*N-7;
    lastdg = 2*(maxL+1);
    H=vector(lastdg);
    H[1]=[0,1,0];
    forstep(i=3,lastdg-1,2,
        /* Part senar */
        g=indcoef(i,H,R);
        d=diagmat(i);
        h=vector(i+1);
        for (j=1,i+1,
            h[j]=g[j]/d[j];
        );
        H[i-1]=h;
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
        H[i]=h;
    );
    return("Centre");
};

ferP(N)=
{
    local(r1,r2,R/*,NN*/);
    r1=vector(N);r1[N]=1;
    r2=vector(N+1);r2[1]=1;
    R=List([r1,r2]);
    /*NN=(N-1)*(N-1);*/
    gettime();
    print(firstlyapunov(R));
    gettime()
};

/* Generar pol z^m*w^n+z^k*w^l en notacio vectorial */
genfield(m:int,n:int,k:int,l:int)=
{
    local(v1:vec,v2:vec);
    v1=vector(m+n+1);
    v1[n+1]=a1+I*b1;
    v2=vector(k+l+1);
    v2[l+1]=a2+I*b2;
    return([v1,v2]);
}

fer()=
{
    for(m=1,3,
        for(n=1,3,
            for(k=1,3,
                for(l=1,3,
                    print(m,",",n,",",k,",",l);
                    print(firstlyapunov(genfield(m,n,k,l)));
                );
            );
        );
    );
}    

