\r polops.gp

n=taskArgs[1];

evenfield(n)=
{
    local(v);
    v=vector(n+1);
    v[1]=-n/(n-2);
    v[n]=1;
    v[n+1]=I*a;
    return(List([v]));
}

{
    L=firstlyapunov(evenfield(n));

    print("QiuYang n = ",n," (expected first nonzero constant: ",(n*n+n-2),")");
    print(" - First nonzero Lyapunov const: \n\tL[",L[1],"] = ",L[2]);

    facts=factor(L[2]);
    j=2;
    tau=-polcoeff(facts[j,1],0,a)/polcoeff(facts[j,1],2,a);

    print("\nWith a² = ",tau);
    for(i=2,5,
        Li=firstlyapunovN(i,evenfield(n))[i];
        if(substpol(Li[2],a^2,tau)!=0,
            break;
        ,
            Li=0;
        );
    );
      
    if(Li!=0,
        print(" - First nonzero Lyapunov const:\n\t L[",Li[1],"] = ",substpol(Li[2],a^2,tau));
        print("\tNumerical approximation = ",substpol(Li[2],a,sqrt(tau)));
    );
    if(Li==0,
        print(" - Suspected center, further analysis required\n");
    );

}
\q
