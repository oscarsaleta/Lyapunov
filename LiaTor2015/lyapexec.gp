\r polops.gp

m=taskArgs[1];
n=taskArgs[2];
k=taskArgs[3];
l=taskArgs[4];

gencleanfield(m,n,k,l)=
{
    local(v1,v2);
    v1=vector(m+n+1);
    v1[n+1]=1;
    v2=vector(k+l+1);
    v2[l+1]=1;
    return(List([v1,v2]));
}

print(firstlyapunov(gencleanfield(m,n,k,l)));

\q
