restart;

# Get degree and tasknumber for reading data
deg:=taskArgs[1]:
taskId:=taskArgs[2]:
fname:=cat("results/task",taskId,"_stdout.txt");

# Set output file
result_fname:=cat("results_maple/",taskId,"output.txt");

# Start JSON output
fprintf(result_fname,"{\n"):

# Print taskId and degree
fprintf(result_fname,"\t\"Id\": %d,\n",taskId):
fprintf(result_fname,"\t\"Degree\": %d,\n",deg):

# Variable conversion for faster computations (sometimes)
a1:=(a+ca)/2:
b1:=(a-ca)/2/I:
a2:=(b+cb)/2:
b2:=(b-cb)/2/I:

# Compute max number of Lyapunov constants
N:=deg^2+3*deg-7:
fprintf(result_fname,"\t\"Number of Lyapunov constants\": %d,\n",N):

# Define variables for Lyapunov constants and initialise to 0
for i from 1 to N do
    L||i:=0:
end do:

# Read data (exponents, Lyapunov constants and reversibility equations)
read fname:

# Print exponents
fprintf(result_fname,"\t\"Exponents\": [%d,%d,%d,%d,%d,%d],\n",k,l,m,n,p,q):

# Print Lyapunov constants
fprintf(result_fname,"\t\"Lyapunov constants\": [\n"):
for i from 1 to N do
    l||i:=simplify(factor(L||i)):
    fprintf(result_fname,"\t\t{ \"L%d\": \"%a\" },\n",i,l||i):
    if (l||i <> 0) then
        maxlyap := i:
    fi:
end do:
fprintf(result_fname,"\t],\n"):

# Print order of max nonzero Lyapunov constant (focus degree)
fprintf(result_fname,"\t\"Maximum focus degree\": %d,\n", maxlyap):

# Include Groebner functions
with(Groebner):

# Define and solve Lyapunov equations
eqs:={seq(l||i,i=1..N)} minus {0}:
lsols:=Solve(eqs,{a,b}):

# Print Lyapunov center conditions
fprintf(result_fname,"\t\"Lyapunov center conditions\": [\n"):
for i1 in lsols do
    s:=solve(i1[1]):
    fprintf(result_fname,"\t\t{ \"%a\" },\n",s):
end do;
fprintf(result_fname,"\t],\n"):

# Define and solve center equations
conds:={c0,c1,c2}:
csols:=Solve(conds,{a,b,x}):

# Print reversible center conditions
fprintf(result_fname,"\t\"Reversible center conditions\": [\n"):
for i2 in csols do
    s := solve(i2[1]):
    fprintf(result_fname,"\t\t{ \"%a\" },\n",s):
end do;
fprintf(result_fname,"\t]\n"):

# End JSON output
fprintf(result_fname,"}\n"):
