/* multiplies 2 homogeneous polynomials as vectors */
vpolmult(P,Q)=
{
    local(i,j,len,res,aux);
    len = #P + #Q - 1;
    res = vector(len);

    for(j=1,#P,
        aux = vector(len);
        for (i=1,#Q,
            aux[i+j-1] = P[j]*Q[i];
            );
        res += aux;
       );
    res
};

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
};

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
};

/* Create diagonal matrix of coefficients */
diagmat(ord)=
{
    local(i,aux);
    aux = vector(ord+1);
    for(i=1,ord+1,aux[i]=ord-2*(i-1););
    aux
};

/* H and R are lists with the vectors needed to compute up to desired
 * order ord (assume #H = #R)
 */
indcoef(deg,H,R)=
{
    local(i,j,k,res);
    res = vector(deg+1);

    for(i=1,#H,
        for(j=1,#R,
            if(#H[i]+#R[j]-3==deg,
                /* Inverse order of conjugate vector */
                aux = vector(#R[j]);
                for(k=1,#aux,
                    aux[k] = conj(R[j][#R[j]-k+1]);
                );
                res += vpolmult(vpoldz(H[i]),R[j])+vpolmult(vpoldw(H[i]),aux);
            );
        );
    );
    I*res
};


/* Test per ordre 3 */
/*print(diagmat(3))
H = List([[0,1,0]]);
R = List([[a20+I*b20,a11+I*b11,a02+I*b02]]);
print(indcoef(3,H,R)~)*/

/* Test per ordre 4 */
/*print(diagmat(4))
listput(H,[h30,h21,h12,h03]);
listput(R,[a30+I*b30,a21+I*b21,a12+I*b12,a03+I*b03]);
print(indcoef(4,H,R))*/

/* Test per ordre 5 */
/*print(diagmat(5))
listput(H,[h40,h31,h22,h13,h04]);
listput(R,[a40+I*b40,a31+I*b31,a22+I*b22,a13+I*b13,a04+I*b04]);
print(indcoef(H,R))*/

/* Test per ordre 6 */
/*print(diagmat(6))
listput(H,[h50,h41,h32,h23,h14,h05]);
listput(R,[a50+I*b50,a41+I*b41,a32+I*b32,a23+I*b23,a14+I*b14,a05+I*b05]);
print(indcoef(H,R))*/


/* Calcula les N primeres constants de Lyapunov del sistema donat */
lyapunov(N,Rtot)=
{
    local(lastdg,H,L,i,j,k,kmax);
    lastdg = 2*(N+1);
    H=List([[0,1,0]]);
    L=List();
    for(i=3,lastdg,
        g=indcoef(i,H,Rtot);
        d=diagmat(i);
        h=vector(i+1);
        if(i%2!=0,
            for(j=1,i+1,
                h[j]=g[j]/d[j];
            );
            listput(H,h);
        ,
            listput(L,g[(i/2)+1]/I);
            print(i,i/2,i/2+1);
            for(j=1,i+1,
                if(d[j]!=0,
                    h[j]=g[j]/d[j];
                );
            );
            listput(H,h);
        );
    );
    L
};

/* Calcula la primera constant de Lyapunov no nul·la i la retorna,
   nomes busca fins la constant N
 */
firstlyapunov(Rtot)=
{
    local(lastdg,H,L,i,j,k,kmax);
    N=0;
    for(i=1,#Rtot,
        N = max(N,#Rtot[i]);
    );
    maxL = N*N+N-2;
    lastdg = 2*(maxL+1);
    H=List([[0,1,0]]);
    for(i=3,lastdg,
        g=indcoef(i,H,Rtot);
        if(i%2==0 && g[(i/2)+1]!=0,
            print(g[(i/2)+1]/I);
            return;
        );
        d=diagmat(i);
        h=vector(i+1);
        if(i%2!=0,
            for(j=1,i+1,
                h[j]=g[j]/d[j];
            );
            listput(H,h);
        ,
            for(j=1,i+1,
                if(d[j]!=0,
                    h[j]=g[j]/d[j];
                );
            );
            listput(H,h);
        );
    );
};

ferP(N)=
{
    local(r1,r2,R,NN);
    r1=vector(N);r1[N]=1;
    r2=vector(N+1);r2[1]=1;
    R=List([r1,r2]);
    NN=(N-1)*(N-1);
    gettime();
    print(firstlyapunov(R));
    gettime()
};

/*
? d3=diagmat(3)
    %9 = [3, 1, -1, -3]
    ? g3=indcoef([[0,1,0]],[[a20+I*b20,a11+I*b11,a02+I*b02]])
    %10 = [-I*a02 - b02, -I*a20 + (b20 + (-I*a11 - b11)), -I*a20 + (-b20 + (-I*a11 + b11)), -I*a02 + b02]
    ? h3=vector(4)
    %11 = [0, 0, 0, 0]
    ? for(i=1,#g3,h3[i]=g3[i]/d3[i];);
    ? h3
    %13 = [-1/3*I*a02 - 1/3*b02, -I*a20 + (b20 + (-I*a11 - b11)), I*a20 + (b20 + (I*a11 - b11)), 1/3*I*a02 - 1/3*b02]
    ? H=List([[0,1,0],h3])
    %14 = List([[0, 1, 0], [(-1/3*I*a02 - 1/3*b02), -I*a20 + (b20 + (-I*a11 - b11)), I*a20 + (b20 + (I*a11 - b11)), (1/3*I*a02 - 1/3*b02)]])
    ? R=List([[a20+I*b20,a11+I*b11,a02+I*b02],[a30+I*b30,a21+I*b21,a12+I*b12,a03+I*b03]])
    %15 = List([[a20 + I*b20, a11 + I*b11, a02 + I*b02], [a30 + I*b30, a21 + I*b21, a12 + I*b12, a03 + I*b03]])
    ? indcoef(H,R)
    %16 = [(-2*a02 + 2*I*b02)*a20 + ((-2*I*a02 - 2*b02)*b20 + ((-a02 + I*b02)*a11 + ((I*a02 + b02)*b11 + (-I*a03 - b03)))), -2*a20^2 + (-4*I*b20 + (-3*a11 + (3*I*b11 + (2*a02 - 2*I*b02))))*a20 + (2*b20^2 + (-3*I*a11 + (-3*b11 + (-2*I*a02 - 2*b02)))*b20 + (-a11^2 + (2*I*b11 + (a02 - I*b02))*a11 + (b11^2 + (I*a02 + b02)*b11 + (-I*a30 + (b30 + (-I*a12 - b12)))))), -2*I*b11*a20 + (-2*I*a11*b20 - 2*I*a21), 2*a20^2 + (-4*I*b20 + (3*a11 + (3*I*b11 + (-2*a02 - 2*I*b02))))*a20 + (-2*b20^2 + (-3*I*a11 + (3*b11 + (-2*I*a02 + 2*b02)))*b20 + (a11^2 + (2*I*b11 + (-a02 - I*b02))*a11 + (-b11^2 + (I*a02 - b02)*b11 + (-I*a30 + (-b30 + (-I*a12 + b12)))))), (2*a02 + 2*I*b02)*a20 + ((-2*I*a02 + 2*b02)*b20 + ((a02 + I*b02)*a11 + ((I*a02 - b02)*b11 + (-I*a03 + b03))))]
    ? g4=indcoef(H,R)
    %17 = [(-2*a02 + 2*I*b02)*a20 + ((-2*I*a02 - 2*b02)*b20 + ((-a02 + I*b02)*a11 + ((I*a02 + b02)*b11 + (-I*a03 - b03)))), -2*a20^2 + (-4*I*b20 + (-3*a11 + (3*I*b11 + (2*a02 - 2*I*b02))))*a20 + (2*b20^2 + (-3*I*a11 + (-3*b11 + (-2*I*a02 - 2*b02)))*b20 + (-a11^2 + (2*I*b11 + (a02 - I*b02))*a11 + (b11^2 + (I*a02 + b02)*b11 + (-I*a30 + (b30 + (-I*a12 - b12)))))), -2*I*b11*a20 + (-2*I*a11*b20 - 2*I*a21), 2*a20^2 + (-4*I*b20 + (3*a11 + (3*I*b11 + (-2*a02 - 2*I*b02))))*a20 + (-2*b20^2 + (-3*I*a11 + (3*b11 + (-2*I*a02 + 2*b02)))*b20 + (a11^2 + (2*I*b11 + (-a02 - I*b02))*a11 + (-b11^2 + (I*a02 - b02)*b11 + (-I*a30 + (-b30 + (-I*a12 + b12)))))), (2*a02 + 2*I*b02)*a20 + ((-2*I*a02 + 2*b02)*b20 + ((a02 + I*b02)*a11 + ((I*a02 - b02)*b11 + (-I*a03 + b03))))]
    ? d4=diagmat(4)
    %18 = [4, 2, 0, -2, -4]
    ? h4=vector(5)
    %19 = [0, 0, 0, 0, 0]
    ? for(i=1,5,if(d4[i]!=0,h4[i]=g4[i]/d4[i];););
    ? h4
    %21 = [(-1/2*a02 + 1/2*I*b02)*a20 + ((-1/2*I*a02 - 1/2*b02)*b20 + ((-1/4*a02 + 1/4*I*b02)*a11 + ((1/4*I*a02 + 1/4*b02)*b11 + (-1/4*I*a03 - 1/4*b03)))), -a20^2 + (-2*I*b20 + (-3/2*a11 + (3/2*I*b11 + (a02 - I*b02))))*a20 + (b20^2 + (-3/2*I*a11 + (-3/2*b11 + (-I*a02 - b02)))*b20 + (-1/2*a11^2 + (I*b11 + (1/2*a02 - 1/2*I*b02))*a11 + (1/2*b11^2 + (1/2*I*a02 + 1/2*b02)*b11 + (-1/2*I*a30 + (1/2*b30 + (-1/2*I*a12 - 1/2*b12)))))), 0, -a20^2 + (2*I*b20 + (-3/2*a11 + (-3/2*I*b11 + (a02 + I*b02))))*a20 + (b20^2 + (3/2*I*a11 + (-3/2*b11 + (I*a02 - b02)))*b20 + (-1/2*a11^2 + (-I*b11 + (1/2*a02 + 1/2*I*b02))*a11 + (1/2*b11^2 + (-1/2*I*a02 + 1/2*b02)*b11 + (1/2*I*a30 + (1/2*b30 + (1/2*I*a12 - 1/2*b12)))))), (-1/2*a02 - 1/2*I*b02)*a20 + ((1/2*I*a02 - 1/2*b02)*b20 + ((-1/4*a02 - 1/4*I*b02)*a11 + ((-1/4*I*a02 + 1/4*b02)*b11 + (1/4*I*a03 - 1/4*b03))))]
    ? L=List([[h4[3]]])
    %22 = List([[0]])
    ? L=List([[h4[3]]])
    %23 = List([[0]])
    ? L=List([[g4[3]]])
    %24 = List([[-2*I*b11*a20 + (-2*I*a11*b20 - 2*I*a21)]])
    ? L[1]
    %25 = [-2*I*b11*a20 + (-2*I*a11*b20 - 2*I*a21)]

*/
