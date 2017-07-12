\r polops.gp

ntask=taskArgs[1];
k=taskArgs[2];
l=taskArgs[3];
m=taskArgs[4];
n=taskArgs[5];
p=taskArgs[6];
q=taskArgs[7];

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

ordres=11;
lyap=firstlyapunovN(ordres,gentrifield(k,l,m,n,p,q));

{
if (lyap=="Centre?",
        print(ntask,",",k,",",l,",",m,",",n,",",p,",",q,",",max(max(k+l,m+n),p+q),",",lyap);
        ,
        if(( (k+l)==(m+n) && (m+n)==(p+q) ),
            print(ntask,",",k,",",l,",",m,",",n,",",p,",",q,",",max(max(k+l,n+m),p+q),",",
                lyap[1][1],",",lyap[2][1],",",lyap[3][1],",",lyap[4][1],",",lyap[5][1],",",lyap[6][1],
                ",","'Homogeni'");
            ,
            print(ntask,",",k,",",l,",",m,",",n,",",p,",",q,",",max(max(k+l,n+m),p+q),",",
                lyap[1][1],",",lyap[2][1],",",lyap[3][1],",",lyap[4][1],",",lyap[5][1],",",lyap[6][1],
                ",","'No homogeni'");
          );

        for (i=1,ordres,
            print("L",lyap[i][1],":=",lyap[i][2],":");
        );
   );
}
