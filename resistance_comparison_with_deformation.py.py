import numpy as np

import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix, csc_matrix
from scipy.sparse.linalg import spsolve
import math

n, m, L, l = 50, 50,  0.0196, 0.00655
x = L / (n-1)
y = l / (m-1)
C = 0.000018
P1 = 10.5
P2 = 10
Y = y**2
X = x**2
def debit(u, n, m, l):
    D = []
    for i in range(n):
        Um = 0
        for j in range(m):
            Um += u[i + j * n]
        D.append(Um * y)
    return D
def directe(i, j, n):
    return (j) * n + i

def inverse(k, n):
    return (k % n, k // n)
def defo(Xmin, Xmax, Ymax):
        L = []
        x = Xmin
        y = 0
        while Xmin <= x <= Xmax and y <= Ymax:
            L.append((x, y))
            if x == Xmax:
                x = Xmin
                y += 1
                continue
            x += 1
        return L

def remplir(xmin, xmax, ymax):
    # Dimensions et paramètres


    # Créer la matrice et le vecteur b
    matrix = lil_matrix((3 * m * n, 3 * m * n))
    b = np.zeros((3 * m * n, 1))

    # Définir les indices de déformation


    Q = defo(xmin, xmax, ymax)

    # Fonctions d'indexation

    # Coins de la matrice
    I = [0, n - 1, m * n - n, m * n - 1]
    LX, LY, LZ = [], [], []

    # Remplissage des conditions
    for k in range(m * n):
        if k not in I and inverse(k, n) not in Q:
            i, j = inverse(k, n)

            if j == 0:  # UX=UY=0
                matrix[directe(i, j, n), directe(i, j, n)] = 1
                LX.append(directe(i, j, n))
                matrix[directe(i, j, n) + m * n, directe(i, j, n) + m * n] = 1
                LY.append(directe(i, j, n) + m * n)

            if j == m - 1:  # UX=UY=0
                matrix[directe(i, j, n), directe(i, j, n)] = 1
                LX.append(directe(i, j, n))
                matrix[directe(i, j, n) + m * n, directe(i, j, n) + m * n] = 1
                LY.append(directe(i, j, n) + m * n)

            if i == 0:  # P1, périodicité UY=0
                matrix[k,directe(i + 1, j, n)] = 1 / X
                matrix[k,directe(n - 2, j, n)] = 1 / X
                matrix[k,directe(i, j + 1, n)] = 1 / Y
                matrix[k,directe(i, j - 1, n)] = 1 / Y
                matrix[k,k] = -2 * ((1 / Y) + (1 / X))
                matrix[k,directe(i + 1, j, n) + 2 * n * m] = -1 / (C * x)
                matrix[k,directe(i, j, n) + 2 * n * m] = 1 / (C * x)
                LX.append(directe(i, j, n))
                matrix[directe(i, j, n) + 2 * m * n, directe(i, j, n) + 2 * m * n] = 1
                b[directe(i, j, n) + 2 * m * n, 0] = P1
                LZ.append(directe(i, j, n) + 2 * m * n)
                matrix[directe(i, j, n) + m * n, directe(i, j, n) + m * n] = 1
                LY.append(directe(i, j, n) + m * n)

            if i == n - 1:  # P2, périodicité UY=0
                matrix[k,directe(1, j, n)] = 1 / X
                matrix[k,directe(i - 1, j, n)] = 1 / X
                matrix[k,directe(i, j + 1, n)] = 1 / Y
                matrix[k,directe(i, j - 1, n)] = 1 / Y
                matrix[k,k] = -2 * ((1 / Y) + (1 / X))
                matrix[k,directe(i, j, n) + 2 * n * m] = -1 / (C * x)
                matrix[k,directe(i - 1, j, n) + 2 * n * m] = 1 / (C * x)
                LX.append(directe(i, j, n))
                matrix[directe(i, j, n) + m * n, directe(i, j, n) + m * n] = 1
                LY.append(directe(i, j, n) + m * n)
                matrix[directe(i, j, n) + 2 * m * n, directe(i, j, n) + 2 * m * n] = 1
                b[directe(i, j, n) + 2 * m * n, 0] = P2
                LZ.append(directe(i, j, n) + 2 * m * n)
            if(j==0 and i!=0 and i!=n-1):
                matrix[k+2*m*n,directe(i+1, j, n)] = 1/(2*x)
                matrix[k+2*m*n,directe(i-1, j, n)] = -1/(2*x)
                matrix[k+2*m*n,directe(i, j+1, n)+m*n] = 1/y
                matrix[k+2*m*n,directe(i, j, n)+m*n] = -1/y
                LZ.append(directe(i,j,n)+2*m*n)
            if(j==m-1 and i!=0 and i!=n-1):
                matrix[k+2*m*n,directe(i+1, j, n)] = 1/(2*x)
                matrix[k+2*m*n,directe(i-1, j, n)] = -1/(2*x)
                matrix[k+2*m*n,directe(i, j, n)+m*n] = 1/y
                matrix[k+2*m*n,directe(i, j-1, n)+m*n] = -1/y
                LZ.append(directe(i,j,n)+2*m*n)
    # Remplissage de l'intérieur
    for k in range(m * n):
        if inverse(k, n) not in Q:
            if k not in LX and k not in I:
                i, j = inverse(k, n)
                matrix[k,directe(i + 1, j, n)] = 1 / X
                matrix[k,directe(i - 1, j, n)] = 1 / X
                matrix[k,directe(i, j + 1, n)] = 1 / Y
                matrix[k,directe(i, j - 1, n)] = 1 / Y
                matrix[k,k] = -2 * ((1 / Y) + (1 / X))
                matrix[k,directe(i + 1, j, n) + 2 * n * m] = -1 / (2 * C * x)
                matrix[k,directe(i - 1, j, n) + 2 * n * m] = 1 / (2 * C * x)

            w = k + m * n
            if w not in LY and k not in I:
                i, j = inverse(k, n)
                matrix[w,directe(i + 1, j, n) + m * n] = 1 / X
                matrix[w,directe(i - 1, j, n) + m * n] = 1 / X
                matrix[w,directe(i, j + 1, n) + m * n] = 1 / Y
                matrix[w,directe(i, j - 1, n) + m * n] = 1 / Y
                matrix[w,w] = -2 * ((1 / Y) + (1 / X))
                matrix[w,directe(i, j + 1, n) + 2 * n * m] = -1 / (2 * C * y)
                matrix[w,directe(i, j - 1, n) + 2 * n * m] = 1 / (2 * C * y)

            v = k + 2 * m * n
            if v not in LZ and k not in I:
                i, j = inverse(k, n)
                matrix[v, directe(i + 1, j, n)] = 1 / (2 * x)
                matrix[v, directe(i - 1, j, n)] = -1 / (2 * x)
                matrix[v, directe(i, j + 1, n) + m * n] = 1 / (2 * y)
                matrix[v, directe(i, j - 1, n) + m * n] = -1 / (2 * y)
    for i in I:

      matrix[i,i]=1
      matrix[i+m*n,i+m*n]=1
      matrix[i+2*m*n,i+2*m*n]=1
      b[0+2*m*n,0]=P1
      b[m*n-n+2*m*n,0]=P1
      b[n-1+2*m*n,0]=P2
      b[m*n-1+2*m*n,0]=P2

    # Remplissage des indices dans Q
    for k in range(m * n):
        i, j = inverse(k, n)
        if inverse(k, n) in Q:
            matrix[k,k] = 1
            matrix[k + m * n,k + m * n] = 1
        if j < ymax and xmin < i < xmax:
            matrix[directe(i, j, n) + 2 * m * n,directe(i, j, n) + 2 * m * n] = 1
        if j == ymax and xmin < i < xmax:
            matrix[k + 2 * m * n,directe(i, j + 1, n) + m * n] = 1
            matrix[k + 2 * m * n,directe(i, j, n) + m * n] = -1
        if i==xmin and 0<j<ymax:
            matrix[k+2*m*n,directe(i,j,n)] =1
            matrix[k+2*m*n,directe(i-1,j,n)] = -1
        if i==xmax and 0<j<ymax:
            matrix[k+2*m*n,directe(i,j,n)] =-1
            matrix[k+2*m*n,directe(i+1,j,n)] = 1
        if (i==xmin and j==ymax ):
      # (Ux(k-1)-Ux(k+1))/dx+(Uy(k-N)-Uy(k+N))/dy=0 le divergent
            matrix[k+2*m*n,directe(i-1,j,n)] = 1/(2*x)
            matrix[k+2*m*n,directe(i+1,j,n)] = - 1/(2*x)
            matrix[k+2*m*n,directe(i,j-1,n)+m*n] = 1/(2*y)
            matrix[k+2*m*n,directe(i,j+1,n)+m*n] = -1/(2*y)
        if ( i==xmax and j==ymax):
            matrix[k+2*m*n,directe(i-1,j,n)] = 1/(2*x)
            matrix[k+2*m*n,directe(i+1,j,n)] = - 1/(2*x)
            matrix[k+2*m*n,directe(i,j-1,n)+m*n] = 1/(2*y)
            matrix[k+2*m*n,directe(i,j+1,n)+m*n] = -1/(2*y)
        if (i==xmin and j==0 ):
            matrix[k+2*m*n,directe(i-1,j,n)+2*m*n] = 1/X
            matrix[k+2*m*n,directe(i,j+1,n)+2*m*n] =  1/Y
            matrix[k+2*m*n,directe(i,j,n)+2*m*n] = -1*((1/Y)+(1/X))
            #matrix[k+2*m*n][directe(i,j,n)+2*m*n] = 1
            #matrix[k+2*m*n][directe(i-1,j,n)+2*m*n] = -1
        if (i==xmax and j==0 ):
            matrix[k+2*m*n,directe(i+1,j,n)+2*m*n] = 1/X
            matrix[k+2*m*n,directe(i,j+1,n)+2*m*n] =  1/Y
            matrix[k+2*m*n,directe(i,j,n)+2*m*n] = -1*((1/Y)+(1/X))
            #matrix[k+2*m*n][directe(i,j,n)+2*m*n] = 1
            #matrix[k+2*m*n][directe(i+1,j,n)+2*m*n] = -1

    return matrix, b

debitm_list = []
ymax_values = range(m-5)  # ymax variant de 0 à 28
pxmax=[]

pxmin=[]
r1=r2=r3=0
R=0
R_eq=[]
Rnumeric=0
R_numeric=[]

for ymax in ymax_values:
    pmax=[] # pression max pour ymax
    pmin=[]
    xmin, xmax = 17, 33
    matrix, b = remplir(xmin, xmax, ymax)
    b = csc_matrix(b)
    matrix=csc_matrix(matrix)
    A=matrix
    solution = spsolve(A, b)
# Extraction de Ux, Uy et P de la solution
    Ux = solution[:m * n].reshape((m, n))
    Uy = solution[m * n:2 * m * n].reshape((m, n))
    P = solution[2 * m * n:].reshape((m, n))
    ux = solution[:m * n]
    uy = solution[m * n:2 * m * n]
    p = solution[2 * m * n:]
    for j in range (m):
      pmax.append(p[directe(xmax,j,n)])
      pmin.append(p[directe(xmin,j,n)])

    U = []

    D = debit(ux, n, m, l)
    debitm = sum(D) / len(D)  # Moyenne du débit pour chaque ymax
    debitm_list.append(debitm)
    pmaxm = sum(pmax) / len(pmax)#pression moyenne pour chaque ymax danx x=xmax
    pminm = sum(pmin) / len(pmin)
    pxmax.append(pmaxm)
    pxmin.append(pminm)
    #calcul de resistance theorique
    r2=12*C*(xmax-xmin)*x/(((m-1-ymax)*y)**3)
    r1=12*C*(n-1-xmax)*x/(l**3)
    r3=12*C*xmin*x/(l**3)
    R=(r1+r2+r3)
    R_eq.append(R)
    RPnumeric=(P1-P2)/debitm
    R_numeric.append(RPnumeric)
    #calcul de resistance numerique

pxmax=[item for item in pxmax]
pxmin=[item for item in pxmin]
#RP=12*c*L/((n-1-ymax)*x)
R_poiseuille=12*C*L/l**3
print(pxmax)
print(pxmin)
print(R_eq)
print(R_numeric)
print(R_poiseuille)
resultat = [a / b for a, b in zip(R_eq, R_numeric) ]
Rap_list = [R_poiseuille / r for r in R_numeric]
log10_resultat = np.log10(resultat)
log10_Rap_list = np.log10(Rap_list)

c=0
abs=[]

for e in range(m-5):
  c=e/n
  c=c*100
  abs.append(c)

# Plot R_list / RP against ymax_values
plt.plot(abs, log10_Rap_list, marker='o', linestyle='-', label='log(R-poiseille / R_numeric)')
plt.plot(abs, log10_resultat, marker='o', linestyle='-', label='log(R-eq / R_numeric)')
plt.xlabel('la valeur de ymax en pourecentage de la largeur totale')
plt.ylabel('Rapport des resistances')
plt.title('rapport du resistance theorique et numeric')
plt.legend()
plt.grid()
plt.show()