set_verbose(-1)

# Llegir dades
taskId=taskArgs[0]
k=taskArgs[1]
l=taskArgs[2]
m=taskArgs[3]
n=taskArgs[4]
p=taskArgs[5]
q=taskArgs[6]

print("R = i*z + z^"+str(k)+"*w^"+str(l)+" + (b1*i)*z^"+str(m)+"*w^"+str(n)+" + (a2+b2*i)*z^"+str(p)+"*w^"+str(q))
grau = max(max(k+l,m+n),p+q)
print("Field degree = "+str(grau))

# Carregar dades i funcions en Pari
gp("taskArgs=["+str(k)+","+str(l)+","+str(m)+","+str(n)+","+str(p)+","+str(q)+"]")
gp("read(\"lyap_a1.gp\")")

# Calcular primera constant no nul·la
gp("l=nextlyapunov(R);")
if sage_eval(gp.eval("l==-1"))==1:
    print("\nR is a center")
    print("\n"+str(taskId)+","+str(k)+","+str(l)+","+str(m)+","+str(n)+","+str(p)+","+str(q)+",Center")
    sys.exit()

R = singular.ring(32003,'(a1,b1,a2,b2)','dp')
print("Using ring on a finite field modulo 32003")

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
    o=gp.eval("l[1]["+str(i+1)+"][1]")
    # Generar ideal amb constants anteriors
    I = singular.ideal(lyaps)
    # Reduir nova constant resp les anteriors
    B = I.groebner()
    # Si redueix, pararem si en portem 2 seguides o passem de n*n+3*n-7
    g = singular(f).sage().reduce(B.sage())
    if g==0 or singular('('+f+')^2').sage().reduce(B.sage())==0:
        print("reduce(L"+o+", "+str(["L"+str(x) for x in ordres])+") = 0")
        reduct += 1
        if int(o)>grau*(grau+3)-7 or reduct>2*grau-1:
            break
    else:
        # Si no redueix, guardem l'ultim ordre
        lyaps.append(str(g))
        print("reduce(L"+o+", "+str(["L"+str(x) for x in ordres])+") != 0")
        ordres.append(o)
        ordre = o
        reduct = 0
    i += 1

print("\nR is a weak focus of order "+ordre)
print("\n"+str(taskId)+","+str(k)+","+str(l)+","+str(m)+","+str(n)+","+str(p)+","+str(q)+","+ordre)
