restart;

primer:=0:

deg:=taskArgs[1]:
taskId:=taskArgs[2]:
if primer <> 0 then
    fname:=cat("results/task",taskId,"_stdout.txt");
    result_fname:=cat("results_maple/",taskId,"output.txt");
else
    fname:=cat("results_no_fin/task",taskId,"_stdout.txt");
    result_fname:=cat("results_maple_no_fin/",taskId,"output.txt");
end if;

fd:=fopen(result_fname,WRITE):

a1:=(a+ca)/2:
b1:=(a-ca)/2/I:
a2:=(b+cb)/2:
b2:=(b-cb)/2/I:

N:=deg^2+3*deg-7:

for i from 1 to N do
    L||i:=0:
end do:

read fname:
fprintf(fd,"# Differential system:\n");
fprintf(fd,"R:=I*z+z^%d*w^%d+a*z^%d*w^%d+b*z^%d*w^%d;\n",k,l,m,n,p,q):

fprintf(fd,"\n# Lyapunov constants:\n");
for i from 1 to N do
    l||i:=simplify(factor(L||i)):
    fprintf(fd,"L%d:=%a;\n",i,l||i):
end do:

eqs:={seq(l||i,i=1..N)} minus {0};
if primer <> 0 then
    conds:={c0 mod primer,c1 mod primer,c2 mod primer};
else
    conds:={c0,c1,c2};
end if;

with(Groebner):
lsols:=Solve(eqs,{a,b}):
csols:=Solve(conds,{a,b,x}):

fprintf(fd,"\n# Lyapunov center conditions:\n");
for i1 in lsols do
    fprintf(fd,"%a\n",solve(i1[1])):
end do;

fprintf(fd,"\n# Reversible center conditions:\n");
for i2 in csols do
    if primer <> 0 then
        fprintf(fd,"%a\n",solve(i2[1]) mod primer):
    else
        fprintf(fd,"%a\n",solve(i2[1])):
    end if;
end do;

fclose(fd);

