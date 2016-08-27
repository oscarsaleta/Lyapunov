restart;

deg:=taskArgs[1]:
taskId:=taskArgs[2]:
fname:=cat("resultats_graus/g",deg,"/task",taskId,"_stdout.txt");

for i from 1 to 11 do
    L||i:=0:
end do:

a1:=(a+ca)/2:
b1:=(a-ca)/2/I:
a2:=(b+cb)/2:
b2:=(b-cb)/2/I:

read fname:
print(k,l,m,n,p,q);

for i from 1 to 11 do
    l||i:=simplify(factor(L||i)):
    print(cat("L",i,":=",l||i));
end do:

eqs:={l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11} minus {0}:
conds:={c0,c1,c2}:

with(Groebner):
lsols:=Solve(eqs,{a,b}):
csols:=Solve(conds,{a,b,x}):

for i1 in lsols do
    print(solve(i1[1]));
end do;

for i2 in csols do
    print(solve(i2[1]));
end do;


