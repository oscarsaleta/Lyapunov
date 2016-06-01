\r polops.gp

m=taskArgs[1];
n=taskArgs[2];
k=taskArgs[3];
l=taskArgs[4];

genfield(k,l,m,n)=
{
    local(v1,v2);
    v1=vector(k+l+1);
    v1[l+1]=1;
    v2=vector(m+n+1);
    v2[n+1]=a+b*I;
    return(List([v1,v2]));
}

print(firstlyapunov(genfield(k,l,m,n)));

\q
