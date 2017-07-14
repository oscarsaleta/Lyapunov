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
print("\n")

# Comprovar si es centre d'algun tipus conegut abans de fer cap calcul
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
from sympy import *

a1 = Symbol('a1')
b1 = Symbol('b1')
a2 = Symbol('a2')
b2 = Symbol('b2')
i = Symbol('i')
z = Symbol('z')
w = Symbol('w')
x = Symbol('x')
y = Symbol('y')

f = i*z + z**k*w**l + (a1+b1*i)*z**m*w**n + (a2+b2*i)*z**p*w**q
f = f.subs(z,x+i*y).subs(w,x-i*y).expand().subs(i,I).subs(I,i)
P = f.coeff(i,0)
Q = f.coeff(i,1)
if diff(P,x) - diff(Q,y) == 0:
    print(status+", HAMILTONIAN CENTER")
    sys.exit(0)
else:
    print(status+", h? "+diff(P,x) - diff(Q,y))


print(status+", NEED TO COMPUTE FURTHER")
sys.exit(0)
