set_verbose(-1)

# Llegir dades
k=taskArgs[0]
l=taskArgs[1]
m=taskArgs[2]
n=taskArgs[3]
p=taskArgs[4]
q=taskArgs[5]

print("R = i*z + z^"+str(k)+"*w^"+str(l)+" + (a1+b1*i)*z^"+str(m)+"*w^"+str(n)+" + (a2+b2*i)*z^"+str(p)+"*w^"+str(q))
grau = max(max(k+l,m+n),p+q)
print("Field degree = "+str(grau))

# Carregar dades i funcions en Pari
gp("taskArgs=["+str(k)+","+str(l)+","+str(m)+","+str(n)+","+str(p)+","+str(q)+"]")
gp("read(\"lyap_3mono.gp\")")

# Calcular primera constant no nul·la
gp("l=nextlyapunov(R);")
if sage_eval(gp.eval("l==-1"))==1:
    print("\nR is a center")
    sys.exit()

#K = GF(32003)
#R = PolynomialRing(K,'a1,a2,b1,b2')
R = PolynomialRing(QQbar,'a1,a2,b1,b2')
a1,a2,b1,b2 = R.gens()
R.inject_variables(verbose=False)
#str_cmd="R.<a1,b1,a2,b2>=PolynomialRing(GF(32003))"
str_cmd="R.<a1,b1,a2,b2>=QQbar[]"

lyaps = []
lyaps.append(sage_eval(gp.eval("l[1][1]"),cmds=str_cmd)[1])
ordres = []
ordres.append(gp.eval("l[1][1][1]"))

print("\nFirst constant: L"+ordres[0]+" = "+str(lyaps[0]))
ordre = ordres[0]

i = 1
reduct = 0
while (int(ordres[0])<=grau*grau+3*grau-7):
    # Calcular constant no nul·la
    gp("l=nextlyapunov(R,l[2],l[1]);")
    f=sage_eval(gp.eval("l[1]["+str(i+1)+"]"),cmds=str_cmd)[1]
    lyaps.insert(0,f)
    ordres.insert(0,gp.eval("l[1]["+str(i+1)+"][1]"))
    # Generar ideal amb constants anteriors
    I = ideal(*lyaps[1:])
    # Reduir nova constant resp les anteriors
    B = I.groebner_basis()
    # Si redueix, pararem si en portem 2 seguides o passem de n*n+n-2
    if lyaps[0].reduce(B)==0:
        print("reduce(L"+ordres[0]+", "+str(["L"+str(x) for x in ordres[i:0:-1]])+") = 0")
        reduct += 1
        if int(ordres[0])>grau*(grau+1)-2 or reduct>2:
            break
    else:
        # Si no redueix, guardem l'ultim ordre
        print("reduce(L"+ordres[0]+", "+str(["L"+str(x) for x in ordres[i:0:-1]])+") != 0")
        ordre = ordres[0]
        reduct = 0
    i += 1

print("\nR is a weak focus of order "+ordre)
print("\n"+str(taskId)+","+str(k)+","+str(l)+","+str(m)+","+str(n)+","+str(p)+","+str(q)+","+ordre)
