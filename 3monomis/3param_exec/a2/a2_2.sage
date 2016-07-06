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

print("R1 = i*z + 1/2*z^"+str(k)+"*w^"+str(l)+" + (a1+b1*i)/2*z^"+str(m)+"*w^"+str(n))
grau1 = max(k+l,m+n)
print("deg(R1) = "+str(grau1))
print("R2 = i*z"+" + 1/2*z^"+str(k)+"*w^"+str(l)+" + (b2*i)/2*z^"+str(p)+"*w^"+str(q))
grau2 = max(k+l,p+q)
print("deg(R2) = "+str(grau2))
print("R3 = i*z"+" + (a1+b1*i)/2*z^"+str(m)+"*w^"+str(n)+" + (b2*i)/2*z^"+str(p)+"*w^"+str(q))
grau3 = max(m+n,p+q)
print("deg(R3) = "+str(grau3))

# Carregar dades i funcions en Pari
gp("taskArgs=["+str(k)+","+str(l)+","+str(m)+","+str(n)+","+str(p)+","+str(q)+"]")
gp("read(\"lyap_a2_2.gp\")")

txt="\n"+str(taskId)+","+str(k)+","+str(l)+","+str(m)+","+str(n)+","+str(p)+","+str(q)+","

# Calcular primera constant no nul·la
gp("l1=nextlyapunov(R1);")
if sage_eval(gp.eval("l1==-1"))==1:
    txt+="SI,"
else:
    txt+="NO,"

gp("l2=nextlyapunov(R2);")
if sage_eval(gp.eval("l2==-1"))==1:
    txt+="SI,"
else:
    txt+="NO,"

gp("l3=nextlyapunov(R3);")
if sage_eval(gp.eval("l3==-1"))==1:
    txt+="SI"
else:
    txt+="NO"

print(txt)
