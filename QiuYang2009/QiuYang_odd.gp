\r polops.gp

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

{
    L=firstlyapunov(oddfield(n));
    print("QiuYang n = ",n," (expected first nonzero constant: ",(n*n+n-2)/2,")");
    print(" - First nonzero Lyapunov const:\n\tL[",L[1][1],"] = ",L[1][2]);
}

{
    fctr=factor(L[1][2]);
    print("\tFactorisation: ",fctr);
    /*for(k=1,#fctr~,
        mypol=fctr[k,1];
        if(poldegree(mypol,a)==2,
            tau=-polcoeff(mypol,0,a)/polcoeff(mypol,2,a);
            print("\nWith a² = ",tau);
            for(i=L[1],n*n+3*n-7,
                Li=lyapunov(i,oddfield(n))[i];
                if(substpol(Li,a^2,tau)!=0,
                    ind=i;
                    break;
                ,
                    Li=0;
                );
            );
            if(Li!=0,
                print(" - First nonzero Lyapunov const:\n\tL[",ind,"] = ",substpol(Li,a^2,tau));
                print("\tNumerical approximation = ",substpol(Li,a,sqrt(tau)));
            );
            if(Li==0,
                print(" - Suspected center, further analysis required\n");
            );
        );
    );*/
}

\q
