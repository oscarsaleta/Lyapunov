set_verbose(-1)

# Llegir dades
taskId=taskArgs[0]
k=taskArgs[1]
l=taskArgs[2]
m=taskArgs[3]
n=taskArgs[4]
p=taskArgs[5]
q=taskArgs[6]

print("\n"+str(taskId)+","+str(k)+","+str(l)+","+str(m)+","+str(n)+","+str(p)+","+str(q)+"\n")

print("R1 = i*z + z^"+str(k)+"*w^"+str(l)+" + (b2*i)*z^"+str(p)+"*w^"+str(q))
grau1 = max(k+l,p+q)
print("deg(R1) = "+str(grau1))
print("R2 = i*z"+" + (a1+b1*i)*z^"+str(m)+"*w^"+str(n)+" + (b2*i)*z^"+str(p)+"*w^"+str(q))
grau2 = max(m+n,p+q)
print("deg(R2) = "+str(grau2))

# Carregar dades i funcions en Pari
gp("taskArgs=["+str(k)+","+str(l)+","+str(m)+","+str(n)+","+str(p)+","+str(q)+"]")
gp("read(\"lyap_a2_2.gp\")")

# Calcular primera constant no nul·la
gp("l=nextlyapunov(R1);")
if sage_eval(gp.eval("l==-1"))==1:
    print("\nR1 is a center")
else:
    print("\nR1 is NOT a center")

gp("l=nextlyapunov(R2);")
if sage_eval(gp.eval("l==-1"))==1:
    print("\nR2 is a center")
else:
    print("\nR2 is NOT a center")
