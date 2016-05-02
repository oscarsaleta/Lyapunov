read("polops.gp");

/* Generar pol z^m*w^n+z^k*w^l en notacio vectorial */
genfield(m,n,k,l)=
{
    local(v1,v2);
    v1=vector(m+n+1);
    v1[n+1]=a1+I*b1;
    v2=vector(k+l+1);
    v2[l+1]=a2+I*b2;
    return(List([v1,v2]));
}

