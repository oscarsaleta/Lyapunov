\r polops_complex.gp

k=taskArgs[1];
l=taskArgs[2];
m=taskArgs[3];
n=taskArgs[4];
p=taskArgs[5];
q=taskArgs[6];

genR_3mono(k,l,m,n,p,q)=
{
    local(v1,v2,v3);
    v1=vector(k+l+1);
    v1[l+1]=1;
    v2=vector(m+n+1);
    v2[n+1]=r1;
    v3=vector(p+q+1);
    v3[q+1]=r2;
    return(List([v1,v2,v3]));
}

genCR_3mono(k,l,m,n,p,q)=
{
    local(v1,v2,v3);
    v1=vector(k+l+1);
    v1[l+1]=-1;
    v2=vector(m+n+1);
    v2[n+1]=s1;
    v3=vector(p+q+1);
    v3[q+1]=s2;
    return(List([v1,v2,v3]));
}

R=genR_3mono(k,l,m,n,p,q);
CR=genCR_3mono(l,k,n,m,q,p);
