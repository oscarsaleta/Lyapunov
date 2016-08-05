set_verbose(-1)

# Llegir dades
k=taskArgs[0]
l=taskArgs[1]
m=taskArgs[2]
n=taskArgs[3]
p=taskArgs[4]
q=taskArgs[5]

status=""+taskId+","+k+","+l+","+m+","+n+","+p+","+q


print("#R = i*z + z^"+str(k)+"*w^"+str(l)+" + (a1+b1*i)*z^"+str(m)+"*w^"+str(n)+" + (a2+b2*i)*z^"+str(p)+"*w^"+str(q))
grau = max(max(k+l,m+n),p+q)
print("#Field degree = "+str(grau))
print("\n")

# Comprovar si es centre d'algun tipus abans de fer cap calcul
# (a) Quadratic Darboux (Zoladek)
if k==2 and q==2 and m==1 and n==1 and l==0 and p==0:
    print(status+"CENTRE A")
    sys.exit(0)
# (b) Holomorphic centers
if l==0 and n==0 and q==0:
    print(status+"CENTRE B")
    sys.exit(0)
# (d) Hamiltonian/new Darboux
if k==m and m==p and l!=n and n!=q and k-l-1!=0:
    print(status+"CENTRE D")
    sys.exit(0)

# (c) Si hem arribat fins aquí, tenim un possible centre reversible
print("\n#Possible centre reversible, calculant constants de Liapunov...")

# Carregar dades i funcions en Pari
gp("taskArgs=["+str(k)+","+str(l)+","+str(m)+","+str(n)+","+str(p)+","+str(q)+"]")
gp("read(\"lyap.gp\")")

# Calcular primera constant no nul·la
gp("l=nextlyapunov(R);")
if sage_eval(gp.eval("l==-1"))==1:
    print("\nR is a center")
    print("\n"+str(taskId)+","+str(k)+","+str(l)+","+str(m)+","+str(n)+","+str(p)+","+str(q)+",Center")
    sys.exit()

R = singular.ring(32003,'(a1,b1,a2,b2)','dp')
print("#Using ring on a finite field modulo 32003")

lyaps = []
lyaps.append(gp.eval("l[1][1][2]"))
ordres = []
ordres.append(gp.eval("l[1][1][1]"))

print("\n#First constant: \nL"+ordres[0]+":="+lyaps[0]+":")
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
    if g==0:
        reduct += 1
        if int(o)>grau*(grau+3)-7 or reduct>grau-1:
            break
    else:
        # Si no redueix, guardem l'ultim ordre
        lyaps.append(str(g))
        print("L"+o+":="+g+":\n")
        ordres.append(o)
        ordre = o
        reduct = 0
    i += 1

test = singular.facstd(I)
print(test)

print("\nR is a weak focus of order "+ordre)
print("\n"+str(taskId)+","+str(k)+","+str(l)+","+str(m)+","+str(n)+","+str(p)+","+str(q)+","+ordre)
