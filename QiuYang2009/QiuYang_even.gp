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

L=firstlyapunov(evenfield(n));

print("Degree = ",n);
print("\tFirst nonzero Lyapunov const: ",L[1],"th = ",L[2]);

facts=factor(L[2])

for(j=1,#(facts~),tau=-polcoeff(facts)[j,1],0,a)/polcoeff(facts)[j,1],2,a);

print("\nWith a² = ",tau);
for(i=1,5,Li=firstlyapunovN(i,evenfield(n))[i];if(substpol(Li[2],a^2,tau)!=0,break;,Li="Center";););
  
if(Li!="Center",print("\tFirst nonzero Lyapunov const: ",Li[1],"th = ",Li[2]);print("\tNumerical approximation = ",substpol(Li[2],a,sqrt(tau))););
if(Li=="Center",print("\tUs)


/*\q*/
