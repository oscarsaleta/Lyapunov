set_verbose(-1)

# Llegir dades
taskId=taskArgs[0]
k=taskArgs[1]
l=taskArgs[2]
m=taskArgs[3]
n=taskArgs[4]
p=taskArgs[5]
q=taskArgs[6]

print("R = i*z + z^"+str(k)+"*w^"+str(l)+" + (a1)*z^"+str(m)+"*w^"+str(n)+" + (a2+b2*i)*z^"+str(p)+"*w^"+str(q))
grau = max(max(k+l,m+n),p+q)
print("Field degree = "+str(grau))

# Carregar dades i funcions en Pari
gp("taskArgs=["+str(k)+","+str(l)+","+str(m)+","+str(n)+","+str(p)+","+str(q)+"]")
gp("read(\"lyap_b1.gp\")")

# Calcular primera constant no nul·la
gp("l=nextlyapunov(R);")
if sage_eval(gp.eval("l==-1"))==1:
    print("\nR is a center")
    print("\n"+str(taskId)+","+str(k)+","+str(l)+","+str(m)+","+str(n)+","+str(p)+","+str(q)+",Center")
    sys.exit()

R = singular.ring(32003,'(a1,a2,b2)','dp')

lyaps = []
lyaps.append(gp.eval("l[1][1][2]"))
ordres = []
ordres.append(gp.eval("l[1][1][1]"))

print("\nFirst constant: L"+ordres[0]+" = "+lyaps[0])
ordre = ordres[0]

i = 1
reduct = 0
while (int(ordres[0])<=grau*grau+3*grau-7):
    # Calcular constant no nul·la
    gp("l=nextlyapunov(R,l[2],l[1]);")
    f=gp.eval("l[1]["+str(i+1)+"][2]")
    lyaps.insert(0,f)
    ordres.insert(0,gp.eval("l[1]["+str(i+1)+"][1]"))
    # Generar ideal amb constants anteriors
    I = singular.ideal(lyaps[1:])
    # Reduir nova constant resp les anteriors
    B = I.groebner()
    # Si redueix, pararem si en portem 2 seguides o passem de n*n+n-2
    if singular(lyaps[0]).sage().reduce(B.sage())==0:
        print("L"+ordres[0]+" reduces with previous: YES")
        reduct += 1
        if int(ordres[0])>grau*(grau+1)-2 or reduct>grau-1:
            break
    else:
        # Si no redueix, guardem l'ultim ordre
        print("L"+ordres[0]+" reduces with previous: NO")
        ordre = ordres[0]
        reduct = 0
    i += 1

print("\nR is a weak focus of order "+ordre)
print("\n"+str(taskId)+","+str(k)+","+str(l)+","+str(m)+","+str(n)+","+str(p)+","+str(q)+","+ordre)
