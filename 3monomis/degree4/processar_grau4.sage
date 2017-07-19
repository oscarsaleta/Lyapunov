# -*- coding: utf-8 -*-
set_verbose(-1)

# Llegir dades
taskId = taskArgs[0]
k = taskArgs[1]
l = taskArgs[2]
m = taskArgs[3]
n = taskArgs[4]
p = taskArgs[5]
q = taskArgs[6]

status = ""+str(taskId)+","+str(k)+","+str(l)+"," + \
    str(m)+","+str(n)+","+str(p)+","+str(q)

print("#R = I*z + z^"+str(k)+"*w^"+str(l)+" + (a1+b1*I)*z^" +
      str(m)+"*w^"+str(n)+" + (a2+b2*I)*z^"+str(p)+"*w^"+str(q))
grau = max(max(k+l, m+n), p+q)
print("#Field degree = "+str(grau))

# Check if it is a trivial center before computing
# Holomorfic
if l == 0 and n == 0 and q == 0:
    print(status+", HOLOMORPHIC CENTER")
    sys.exit(0)

# New Darboux?
if k >= 1 and k == m and p >= k+1 and q == p-k:
    print(status+", DARBOUX CENTER with integrating factor 1/(z*w)")
    sys.exit(0)
elif k >= 1 and k == p and m >= k+1 and n == m-k:
    print(status+", DARBOUX CENTER with integrating factor 1/(z*w)")
    sys.exit(0)

# Hamiltonian
import sympy as sp
a1,b1,a2,b2 = sp.symbols('a1 b1 a2 b2')
i, z, w, x, y = sp.symbols('i z w x y')

f = i*z + z**k*w**l + (a1+b1*i)*z**m*w**n + (a2+b2*i)*z**p*w**q
g = f.subs({z:x+i*y,w:x-i*y}).expand().subs(i, sp.I).subs(sp.I, i)
P = g.coeff(i, 0)
Q = g.coeff(i, 1)
if diff(P, x) - diff(Q, y) == 0:
    print(status+", HAMILTONIAN CENTER")
    sys.exit(0)


# Center not included in previous filters
print("\n#Computing Lyapunov constants...")

# Load data and functions in PARI
gp("taskArgs=["+str(k)+","+str(l)+"," +
   str(m)+","+str(n)+","+str(p)+","+str(q)+"]")
gp("read(\"lyaps.gp\")")

# Compute first nonzero constant
gp("l=nextlyapunov(R);")
if sage_eval(gp.eval("l==-1")) == 1:
    print("\n#"+status+", Center")
    sys.exit()

# Reduce Lyapunov constants in Groebner bases
primer = 32003
R = singular.ring(primer, '(a1,b1,a2,b2)', 'dp')
if (primer != 0):
    print("#Using ring on a finite field modulo 32003")

lyaps = []
lyaps.append(gp.eval("l[1][1][2]"))
ordres = []
ordres.append(gp.eval("l[1][1][1]"))

print("\nL"+ordres[0]+":="+lyaps[0]+":\n")
ordre = ordres[0]

i = 1
reduct = 0
while (int(ordres[0]) <= grau*grau+3*grau-7):
    # Compute next nonzero constant
    gp("l=nextlyapunov(R,l[2],l[1]);")
    f = gp.eval("l[1]["+str(i+1)+"][2]")
    o = gp.eval("l[1]["+str(i+1)+"][1]")
    # Generate ideal with previous constants
    I = singular.ideal(lyaps)
    # Reduce new constant with respect to the previous
    B = I.groebner()
    # Stopping condition: 2 reduce in a row or we computed more than n*n+3*n-7
    g = singular(f).sage().reduce(B.sage())
    if g == 0:
        reduct += 1
        print("L"+o+":="+str(g)+": #reduced\n")
        if int(o) > grau*(grau+2)-1 or reduct > 3/2*grau:
            break
    else:
        # If it does not reduce, store the order
        lyaps.append(str(g))
        print("L"+o+":="+str(g)+":\n")
        ordres.append(o)
        ordre = o
        reduct = 0
    i += 1


# Reversible center conditions
print("# Computing reversible center conditions\n")
if primer != 0:
    S.<i> = GF(primer)[]
    SS.<I> = S.quotient(i^2+1)
else:
    SS.<I> = QQ[I]
K.<a1,a2,b1,b2,x> = PolynomialRing(SS)
c0 = numerator(1+x^(k-l-1))
c1 = numerator((a1+b1*I)+(a1-b1*I)*x^(m-n-1))
c2 = numerator((a2+b2*I)+(a2-b2*I)*x^(p-q-1))
#condicions = singular.facstd(singular.ideal(c0,c1,c2))
print("c0:="+str(c0)+"\n")
print("c1:="+str(c1)+"\n")
print("c2:="+str(c2)+"\n")
