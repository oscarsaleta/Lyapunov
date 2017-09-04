restart;

# Define if we are in rationals (primer=0) or finite field
primer:=0:

# Read PBala arguments
deg:=taskArgs[1]:
taskId:=taskArgs[2]:

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
N:=deg^2+3*deg-7:
fprintf(result_fname,"\t\"Number of Lyapunov constants\": %d,\n",N):

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
    l||i:=simplify(factor(L||i)):
    fprintf(fd,"L%d:=%a;\n",i,l||i):
end do:
fprintf(result_fname,"\t],\n"):

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
eqs:={seq(l||i,i=1..N)} minus {0}:
lsols:=Solve(eqs,{a,b}):

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
fprintf(result_fname,"\t],\n"):

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
        ii2:=solve(eq,{a,b,x}) mod primer:
    else
        ii2:=solve(eq,{a,b,x}):
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
LSOLS:={simplify(expand(lsols0))}:
for i from 1 to n_lsols-1 do
    LSOLS:=LSOLS union {simplify(expand(lsols||i))};
end do;
CSOLS:={simplify(expand(csols0))}:
for i from 1 to n_csols-1 do
    CSOLS:=CSOLS union {simplify(expand(csols||i))};
end do;
#fprintf(fd,"\nLSOLS:=%a\nCSOLS:=%a\n",LSOLS,CSOLS);

# Compare conditions and remove Lyapunov ones that are reversible
for c in CSOLS do
    for l in LSOLS do
        if l subset c and c subset l then
            LSOLS:=LSOLS minus {l};
            break;
        end if;
    end do;
end do;

# Print remaining Lyapunov center conditions
if nops(LSOLS)>0 then
    fprintf(fd,"\n# Non-reversible center conditions:\n");
    for l in LSOLS do
        fprintf(fd,"%a\n",l);
    end do;
else
    fprintf(fd,"\n# All center conditions are reversible\n");
end if;

fclose(fd);

