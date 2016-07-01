\r polops.gp

k=taskArgs[1];
l=taskArgs[2];
m=taskArgs[3];
n=taskArgs[4];
p=taskArgs[5];
q=taskArgs[6];

genfield1(k,l,p,q)=
{
    local(v1,v3);
    v1=vector(k+l+1);
    v1[l+1]=1;
    v3=vector(p+q+1);
    v3[q+1]=b2*I;
    return(List([v1,v3]));
}
genfield2(k,l,p,q)=
{
    local(v2,v3);
    v2=vector(m+n+1);
    v2[n+1]=a1+b1*I;
    v3=vector(p+q+1);
    v3[q+1]=b2*I;
    return(List([v2,v3]));
}

R1=genfield1(k,l,p,q);
R2=genfield2(m,n,p,q);
