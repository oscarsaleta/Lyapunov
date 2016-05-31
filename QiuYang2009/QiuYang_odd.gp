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

L=firstlyapunov(oddfield(n));
print("Degree = ",n);
print("First nonzero Lyapunov const: ",L[1],"th = ",L[2]);

tau=-polcoeff(factor(L[1][2])[2,1],0,a)/polcoeff(factor(L[1][2])[2,1],2,a);

print("\nWith a²=",tau);
print(L[1][1],"th Lyapunov const = ",substpol(L[1][2],a^2,tau));
print(L[2][1],"th Lyapunov const = ",substpol(L[2][2],a^2,tau));


/*\q*/
