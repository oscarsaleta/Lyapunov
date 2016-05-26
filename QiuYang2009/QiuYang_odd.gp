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

L=firstlyapunovN(2,oddfield(n));
print(L[1][1],",",L[1][2]);
print(L[2][1],",",L[2][2]);

write("results_odd/to_solve.txt",n,",",L[1][2])

\q
