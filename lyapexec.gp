\r /home/osr/pari/Lyapunov/liblyap.gp

m=taskArgs[1];
n=taskArgs[2];
k=taskArgs[3];
l=taskArgs[4];

print(firstlyapunov(gencleanfield(m,n,k,l)));

\q
