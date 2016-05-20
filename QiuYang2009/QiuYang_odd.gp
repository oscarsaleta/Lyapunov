\r ../polops.gp

n=taskArgs[1];

oddfield(n)=
{
    local(v);
    v=vector(n+1);
    v[1]=n/(n-2);
    v[n]=1;
    v[n+1]=1+I*a;
    return(List([v]));
}

L=firstlyapunovN(2,evenfield(n));
print(L[1][1]);
print(L[1][2]);
/*s=solve

\q
