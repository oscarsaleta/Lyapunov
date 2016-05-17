/*\r /home/osr/pari/Lyapunov/polops.gp*/
\r polops.gp

n=taskArgs[1];

evenfield(n)=
{
    local(v);
    v=vector(n+1);
    v[1]=n/(n-2);
    v[n]=1;
    v[n+1]=I*a;
    return(List([v]));
}

print(firstlyapunovN(2,evenfield(n)));

\q
