restart;

deg:=taskArgs[1]:
taskId:=taskArgs[2]:
fname:=cat("results/task",taskId,"_stdout.txt");

result_fname:=cat("results_maple/",taskId,"output.txt");
writeto(result_fname);

a1:=(a+ca)/2:
b1:=(a-ca)/2/I:
a2:=(b+cb)/2:
b2:=(b-cb)/2/I:

N:=deg^2+3*deg-7:

for i from 1 to N do
    L||i:=0:
end do:

read fname:
print(k,l,m,n,p,q):

for i from 1 to N do
    l||i:=simplify(factor(L||i)):
    print(cat("L",i,":=",l||i)):
end do:

eqs:={seq(l||i,i=1..N)} minus {0}:
conds:={c0,c1,c2}:

with(Groebner):
lsols:=Solve(eqs,{a,b}):
csols:=Solve(conds,{a,b,x}):

for i1 in lsols do
    print(solve(i1[1])):
end do;

for i2 in csols do
    print(solve(i2[1])):
end do;


