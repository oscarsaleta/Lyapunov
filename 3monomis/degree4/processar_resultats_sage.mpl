restart;
Digits:=50;

# Define if we are in rationals (primer=0) or finite field
primer:=0;

# Read PBala arguments
deg:=taskArgs[1];
taskId:=taskArgs[2];

# Pick correct filenames depending on working field
if primer <> 0 then
    fname:=cat("results/task",taskId,"_stdout.txt");
    result_fname:=cat("results_maple/",taskId,"output.txt");
else
    fname:=cat("results_no_fin/task",taskId,"_stdout.txt");
    result_fname:=cat("results_maple_no_fin/",taskId,"output.txt");
end if;

# Open output file
fd:=fopen(result_fname,WRITE):

# Define A and conj(A) instead of A=a+b*i
a1:=(a+ca)/2:
b1:=(a-ca)/2/I:
a2:=(b+cb)/2:
b2:=(b-cb)/2/I:

# Max number of Lyapunov constants computed
N:=deg^2+3*deg-7;

# Initialise all Lyapunov constants to 0
for i from 1 to N do
    L||i:=0:
end do:

# Read data from Sage PBala execution
read fname:
fprintf(fd,"# Differential system:\n");
fprintf(fd,"R:=I*z+z^%d*w^%d+a*z^%d*w^%d+b*z^%d*w^%d;\n",k,l,m,n,p,q):

# Simplify and print Lyapunov constants
fprintf(fd,"\n# Lyapunov constants:\n");
for i from 1 to N do
    l||i:=simplify(factor(L||i));
    fprintf(fd,"L%d:=%a;\n",i,l||i);
end do:

# Compute Lyapunov center conditions
eqs:={seq(l||i,i=1..N)} minus {0};

# Compute reversible center conditions
if primer <> 0 then
    conds:={c0 mod primer,c1 mod primer,c2 mod primer};
else
    conds:={c0,c1,c2};
end if;
#fprintf(fd,"\nconds:=%a\n",conds);

# Solve previous equations using Groebner basis
with(Groebner):

# Define and solve Lyapunov equations
eqs:={seq(l||i,i=1..N)} minus {0};
lsols:=Solve(eqs,{a,b});

# Simplify and write Lyapunov solution sets
fprintf(fd,"\n# Lyapunov center conditions:\n");
n_lsols:=0:
for i1 in lsols do
    eq:={op(i1[1])};
    if nops(i1[3])>0 then
        eq:=eq union {op(i1[3])<>0};
    end if;
    ii1:=solve(eq,{a,b});
    if whattype(ii1)=exprseq then
        for v in ii1 do
            vv:=allvalues(v);
            if whattype(vv)=exprseq then
                for vvv in vv do
                    lsols||n_lsols:=vvv;
                    fprintf(fd,"%a;\n",vvv);
                    n_lsols:=n_lsols+1;
                end do;
            else
                lsols||n_lsols:=v;
                fprintf(fd,"%a;\n",v);
                n_lsols:=n_lsols+1;
            end if;
        end do;
    else
        v:=allvalues(ii1);
        if whattype(v)=exprseq then
            for vv in v do
                lsols||n_lsols:=vv;
                fprintf(fd,"%a;\n",vv);
                n_lsols:=n_lsols+1;
            end do;
        else
            lsols||n_lsols:=ii1;
            fprintf(fd,"%a;\n",ii1);
            n_lsols:=n_lsols+1;
        end if;
    end if;
end do;

# Define and solve center equations
conds:={c0,c1,c2}:
csols:=Solve(conds,{a,b,x}):

# Simplify and write reversible center solution sets
fprintf(fd,"\n# Reversible center conditions:\n");
with(ListTools):
n_csols:=0;
for i2 in csols do
    eq:={op(i2[1])};
    if nops(i2[3])>0 then
        eq:=eq union {op(i2[3])<>0};
    end if;
    if primer <> 0 then
        ii2:=solve(eq,{a,b,x}) mod primer;
    else
        ii2:=solve(eq,{a,b,x});
    end if;
    ii2:=ii2 minus {SelectLast([op(ii2)])};
    if whattype(ii2)=exprseq then
        for v in ii2 do
            vv:=allvalues(ii2);
            if whattype(vv)=exprseq then
                for vvv in vv do
                    csols||n_csols:=vvv;
                    fprintf(fd,"%a;\n",csols||n_csols);
                    n_csols:=n_csols+1;
                end do;
            else
                csols_n_csols:=v;
                fprintf(fd,"%a;\n",csols||n_csols);
                n_csols:=n_csols+1;
            end if;
        end do;
    else
        v:=allvalues(ii2);
        if whattype(v)=exprseq then
            for vv in v do
                csols||n_csols:=vv;
                fprintf(fd,"%a;\n",csols||n_csols);
                n_csols:=n_csols+1;
            end do;
        else
            csols||n_csols:=ii2;
            fprintf(fd,"%a;\n",csols||n_csols);
            n_csols:=n_csols+1;
        end if;
    end if;
end do;

# Create sets of all conditions for Lyapunov and reversible center
LSOLS:={simplify(expand(simplify(lsols0)))};
for i from 1 to n_lsols-1 do
    LSOLS:=LSOLS union {simplify(expand(simplify(lsols||i)))};
end do;
CSOLS:={simplify(expand(simplify(csols0)))};
for i from 1 to n_csols-1 do
    CSOLS:=CSOLS union {simplify(expand(simplify(csols||i)))};
end do;
#fprintf(fd,"\nLSOLS:=%a\nCSOLS:=%a\n",LSOLS,CSOLS);

# Compare conditions and remove Lyapunov ones that are reversible
for c in CSOLS do
    for l in LSOLS do
        if l subset c and c subset l then
            LSOLS:=LSOLS minus {l};
            break;
        elif abs(coeff(rhs(evalf(l)[1])-rhs(evalf(c)[1]),ca))<1e-49 and abs(coeff(rhs(evalf(l)[2])-rhs(evalf(c)[2]),cb))<1e-49 then
            LSOLS:=LSOLS minus {l};
            break;
        end if;
    end do;
end do;

# Print remaining Lyapunov center conditions
if nops(LSOLS)>0 then
    fprintf(fd,"\n# Non-reversible center conditions:\n");
    n_nrconds:=0;
    for l in LSOLS do
        c||n_nrconds:=l;
        fprintf(fd,"c%d:=%a;\n",n_nrconds,l);
        n_nrconds:=n_nrconds+1;
    end do;
else
    n_nrconds:=-1;
    fprintf(fd,"\n# All center conditions are reversible\n");
end if;

# Test conditions for Hamiltonian/easy Darboux integrability
a1:='a1';a2:='a2';b1:='b1';b2:='b2';
defs:={a=a1+a2*I,ca=a1-a2*I,b=b1+b2*I,cb=b1-b2*I};
for i from 0 to n_nrconds-1 do
    fprintf(fd,"\n# Non-reversible condition c%d\n", i);
    RR:=subs(solve(c||i,{a,ca,b,cb}),R);
    expand(subs(z=x+y*I,w=x-y*I,defs,RR));
    P,Q:=coeff(%,I,0),coeff(%,I,1);
    fprintf(fd,"P:=%a;\nQ:=%a;\n",P,Q);
    if diff(P,x)+diff(Q,y)=0 then
        fprintf(fd,"# The center is Hamiltonian\n");
    else
        #fprintf(fd,"# The center is not Hamiltonian\n");
        Nn:=1;
        K:=b00+b10*x+b01*y+b20*x^2+b11*x*y+b02*y^2;
        #fprintf(fd,"K:=%a;\n",K);
        F:=1+sum(sum(a[j1,j2-j1]*x^j1*y^(j2-j1),j1=0..j2),j2=1..Nn);
        #fprintf(fd,"F:=%a;\n",F);
        var:=indets(F*K) minus ({x,y});
        #fprintf(fd,"var:=%a;\n",var);
        expand(diff(F,x)*P+diff(F,y)*Q-F*K):
        FKsols:={solve({coeffs(%,[x,y]),b20<>0},var)};
        #fprintf(fd,"FKsols:=%a;\n",FKsols);
        if nops(FKsols)>0 then
            fprintf(fd,"# The system could be Darboux integrable\n");
            fprintf(fd,"# TODO: find cofactor\n");
        else
            fprintf(fd,"# Could not find algebraic curves of degree %d;\n",Nn);
        end if;
    end if;
end do;

fclose(fd);
