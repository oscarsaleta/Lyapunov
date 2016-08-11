# -*- coding: utf-8 -*-
set_verbose(-1)

# Llegir dades
taskId=1
k=0
l=4
m=1
n=1
p=1
q=4

status=""+str(taskId)+","+str(k)+","+str(l)+","+str(m)+","+str(n)+","+str(p)+","+str(q)


print("#R = i*z + z^"+str(k)+"*w^"+str(l)+" + (a1+b1*i)*z^"+str(m)+"*w^"+str(n)+" + (a2+b2*i)*z^"+str(p)+"*w^"+str(q))
grau = max(max(k+l,m+n),p+q)
print("#Field degree = "+str(grau))
print("\n")

# Comprovar si es centre d'algun tipus abans de fer cap calcul
# (a) Quadratic Darboux (Zoladek)
if k==2 and q==2 and m==1 and n==1 and l==0 and p==0:
    print(status+",CENTRE A")
    sys.exit(0)
# (b) Holomorphic centers
if l==0 and n==0 and q==0:
    print(status+",CENTRE B")
    sys.exit(0)
# (d) Hamiltonian/new Darboux
if k==m and m==p and l!=n and n!=q and k-l-1!=0:
    print(status+",CENTRE D")
    sys.exit(0)

# (c) Si hem arribat fins aqui, tenim un possible centre reversible
print("\n#Possible reversible centre, computing Liapunov constants...")

# Carregar dades i funcions en Pari
gp("taskArgs=["+str(k)+","+str(l)+","+str(m)+","+str(n)+","+str(p)+","+str(q)+"]")
gp("read(\"lyap.gp\")")

# Calcular primera constant no nul·la
gp("l=nextlyapunov(R);")
if sage_eval(gp.eval("l==-1"))==1:
    print("\n#"+status+"Center")
    sys.exit()

#primer = 32003
primer = 0
R = singular.ring(primer,'(a1,b1,a2,b2)','dp')
if (primer!=0):
    print("#Using ring on a finite field modulo 32003")

lyaps = []
lyaps.append(gp.eval("l[1][1][2]"))
ordres = []
ordres.append(gp.eval("l[1][1][1]"))

print("\nL"+ordres[0]+":="+lyaps[0]+":\n")
ordre = ordres[0]

i = 1
reduct = 0
while (int(ordres[0])<=grau*grau+3*grau-7 and len(ordres)<4):
    # Calcular constant no nul·la
    gp("l=nextlyapunov(R,l[2],l[1]);")
    f=gp.eval("l[1]["+str(i+1)+"][2]")
    o=gp.eval("l[1]["+str(i+1)+"][1]")
    # Generar ideal amb constants anteriors
    Ii = singular.ideal(lyaps)
    # Reduir nova constant resp les anteriors
    B = Ii.groebner()
    # Si redueix, pararem si en portem 2 seguides o passem de n*n+3*n-7
    g = singular(f).sage().reduce(B.sage())
    if g==0:
        reduct += 1
        print("L"+o+":="+str(g)+": #reduced\n")
        if int(o)>grau*(grau+3)-7 or reduct>grau:
            break
    else:
        # Si no redueix, guardem l'ultim ordre
        lyaps.append(str(g))
        print("L"+o+":="+str(g)+":\n")
        ordres.append(o)
        ordre = o
        reduct = 0
    i += 1

test = singular.facstd(Ii).sage()
print(test)

y = var('y')
valors_x = solve(1+y^(k-l-1)==0,y)
print(valors_x)


# Busquem les condicions de centre
if primer!=0:
    S.<i>=GF(primer)
    RR.<I>=S.quotient(i^2+1)
else:
    RR.<I>=QQ[I]
K.<a1,a2,b1,b2,x>=PolynomialRing(RR)
#singular.eval("ring R = (0,I),(a1,a2,b1,b2,x),dp;")
#singular.eval("minpoly = I^2+1;")
c0=numerator(1+x^(k-l-1))
#singular.eval("poly c0 = "+str(c0)+";")
c1=numerator((a1+b1*I)+(a1-b1*I)*x^(m-n-1))
#singular.eval("poly c1 = "+str(c1)+";")
c2=numerator((a2+b2*I)+(a2-b2*I)*x^(p-q-1))
#singular.eval("poly c2 = "+str(c2)+";")
condicions = singular.facstd(singular.ideal(c0,c1,c2))
print(condicions)


#print("\nR is a weak focus of order "+ordre)
#print("\n"+str(taskId)+","+str(k)+","+str(l)+","+str(m)+","+str(n)+","+str(p)+","+str(q)+","+ordre)
