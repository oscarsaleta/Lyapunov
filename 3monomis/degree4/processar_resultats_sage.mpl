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

# Compute Lyapunov center conditions
eqs:={seq(l||i,i=1..N)} minus {0};

# Compute reversible center conditions
if primer <> 0 then
    conds:={c0 mod primer,c1 mod primer,c2 mod primer};
else
    conds:={c0,c1,c2};
end if;

# Solve previos equations using Groebner basis
with(Groebner):
lsols:=Solve(eqs,{a,b}):
csols:=Solve(conds,{a,b,x}):

# Simplify and write Lyapunov solution sets
fprintf(fd,"\n# Lyapunov center conditions:\n");
n_lsols:=0:
for i1 in lsols do
    ii1:=solve(i1[1]);
    fprintf(fd,"%a\n",ii1);
    for v in allvalues(ii1) do
        if whattype(v)=set then
            lsols||n_lsols:=v;
            fprintf(fd,"%a\n",v):
        else
            lsols||n_lsols:=ii1;
        end if;
        n_lsols:=n_lsols+1:
    end do;
end do;

# Simplify and write reversible center solution sets
fprintf(fd,"\n# Reversible center conditions:\n");
with(ListTools):
n_csols:=0;
for i2 in csols do
    if primer <> 0 then
        ii2:=solve(i2[1]) mod primer:
    else
        ii2:=solve(i2[1]):
    end if;
    ii2:=ii2 minus {SelectLast([op(ii2)])};
    for v in allvalues(ii2) do
        if whattype(v)=set then
            csols||n_csols:=v;
            fprintf(fd,"%a\n",v):
        else
            csols||n_csols:=ii2;
        end if;
        n_csols:=n_csols+1:
    end do;
end do;

# Compare Lyapunov and reversible soltion sets
LSOLS:={lsols0}:
for i from 1 to n_lsols-1 do
    LSOLS:=LSOLS union {lsols||i};
CSOLS:={csols0}:
for i from 1 to n_csols-1 do
    CSOLS:=CSOLS union {csols||i};
end do;

for c in CSOLS do
    for l in LSOLS do
        if l subset c and c subset l then
            LSOLS:=LSOLS minus l;
            break;
        end if;
    end do;
end do;

# Print remaining Lyapunov center conditions
if nops(LSOLS)>0 then
    fprintf(fd,"\n# Conflicting Lyapunov conditions:\n");
    for l in LSOLS do
        fprintf(fd,"%a\n",l);
    end do;
else
    fprintf(fd,"\n# All Lyapunov center conditions are reversible\n");
end if;

fclose(fd);

