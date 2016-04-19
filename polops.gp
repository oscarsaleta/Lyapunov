/* multiplies 2 homogeneous polynomials as vectors */
vpolmult(P,Q)=
{
    local(i,j,len,res,aux);
    len = #P + #Q - 1;
    res = vector(len);
    aux = vector(len);

    for(j=1,#P,
        for (i=1,j-1,
            aux[i] = 0;
            );
        for (i=1,#Q,
            aux[i+j-1] = P[j]*Q[i];
            );
        for (i=j+#Q,len,
            aux[i] = 0;
            );
        res += aux;
       );
    res
}

/* differentiates a homogeneous polynomial with resp. to z */
vpoldz(P)=
{
    local(i,deg,res);
    deg = #P-1;
    res = vector(deg);

    for(i=1,deg,
        res[i] = (deg-i+1)*P[i];
       );
    res
}

/* differentiates a homogeneous polynomial with resp. to w */
vpoldw(P)=
{
    local(i,deg,res);
    deg = #P-1;
    res = vector(deg);

    for(i=1,deg,
        res[deg-i+1] = (deg-i+1)*P[#P-i+1];
       );
    res
}

/* Create diagonal matrix of coefficients */
diagmat(ord)=
{
    local(i,aux);
    aux = vector(ord+1);
    for(i=1,ord+1,aux[i]=1;);
    matdiagonal(vpolmult(vpoldz(aux),[1,0])+vpolmult(vpoldw(aux),[0,-1]))
}

/* H and R are lists with the vectors needed to compute up to desired
 * order ord (assume #H = #R)
 */
indcoef(H,R)=
{
    local(i,j,res);
    res = vector((#H[1]-1)+#R[#R]-1);

    for(i=1,#H,
        ii = #R-i+1;
        /* Inverse order of conjugate vector */
        aux = vector(#R[ii]);
        for (j=1,#R[ii],
            aux[j] = conj(R[ii][#R[ii]-j+1]);
            );
        res += vpolmult(vpoldz(H[i]),R[#R-i+1])+vpolmult(vpoldw(H[i]),aux);
       );
    res/I
}

/* Test per ordre 3 */
print(diagmat(3))
H = List([[0,1,0]]);
R = List([[a20+I*b20,a11+I*b11,a02+I*b02]]);
print(indcoef(H,R)~)

/* Test per ordre 4 */
print(diagmat(4))
listput(H,[h30,h21,h12,h03]);
listput(R,[a30+I*b30,a21+I*b21,a12+I*b12,a03+I*b03]);
print(indcoef(H,R))

/* Test per ordre 5 */
print(diagmat(5))
listput(H,[h40,h31,h22,h13,h04]);
listput(R,[a40+I*b40,a31+I*b31,a22+I*b22,a13+I*b13,a04+I*b04]);
print(indcoef(H,R))

/* Test per ordre 6 */
print(diagmat(6))
listput(H,[h50,h41,h32,h23,h14,h05]);
listput(R,[a50+I*b50,a41+I*b41,a32+I*b32,a23+I*b23,a14+I*b14,a05+I*b05]);
print(indcoef(H,R))
