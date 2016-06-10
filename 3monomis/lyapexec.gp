\r polops.gp

k=taskArgs[1];
l=taskArgs[2];
m=taskArgs[3];
n=taskArgs[4];
p=taskArgs[5];
q=taskArgs[6];

gentrifield(k,l,m,n,p,q)=
{
    local(v1,v2,v3);
    v1=vector(k+l+1);
    v1[l+1]=1;
    v2=vector(m+n+1);
    v2[n+1]=a1+b1*I;
    v3=vector(p+q+1);
    v3[q+1]=a2+b2*I;
    return(List([v1,v2,v3]));
}

lyap=firstlyapunovN(3,gentrifield(k,l,m,n,p,q));
{
if (lyap=="Centre?",
        print(k,",",l,",",m,",",n,",",p,",",q,",",max(max(k+l,m+n),p+q),",",lyap);
        ,
        print(k,",",l,",",m,",",n,",",p,",",q,",",max(max(k+l,m+n),p+q),",",lyap[1],lyap[2],lyap[3]);
   );
}
