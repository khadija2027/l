import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.sparse import lil_matrix, csc_matrix
from scipy.sparse.linalg import spsolve

# Dimensions et paramètres

n, m, L, l = 50, 50,  0.0196, 0.00655
x = L / (n-1)
y = l / (m-1)
C = 0.000018
P1 = 10.5
P2 = 10
Y = y**2
X = x**2
# Créer la matrice sparse et le vecteur b
matrix = lil_matrix((3 * m * n, 3 * m * n))  # Création d'une matrice creuse
b = np.zeros((3 * m * n, 1))
# Listes pour suivre les indices de remplissage
LX = []
LY = []
LZ = []
def defo(Xmin,Xmax,Ymax):
    L=[]
    x=Xmin
    y=0
    while Xmin <= x <= Xmax and y <= Ymax :

        L.append((x,y))
        if x==Xmax :
            x=Xmin
            y+=1
            continue
        x+=1
    return L
xmin, xmax, ymax = 17, 33, 23
Q = defo(xmin, xmax, ymax)

# Fonctions auxiliaires
def directe(i, j, n):
    return j * n + i

def inverse(k, n):
    return k % n, k // n
I = [0, n-1, m*n-n, m*n-1]
for k in range(m*n)  :
    if k not in I and inverse(k, n) not in Q :
        i=inverse(k, n)[0]
        j=inverse(k, n)[1]
        if(j==0):#UX=UY=0
            matrix[directe(i,j,n),directe(i,j,n)]=1
            LX.append(directe(i,j,n))
            matrix[directe(i,j,n)+m*n,directe(i,j,n)+m*n]=1
            LY.append(directe(i,j,n)+m*n)
        if(j==m-1):#UX=UY=0
            matrix[directe(i,j,n),directe(i,j,n)]=1
            LX.append(directe(i,j,n))
            matrix[directe(i,j,n)+m*n,directe(i,j,n)+m*n]=1
            LY.append(directe(i,j,n)+m*n)
        if(i==0):#P1 , PERIDICITE UY=0
            matrix[k,directe(i+1, j, n)] = 1/X
            matrix[k,directe(n-2, j, n)] = 1/X
            matrix[k,directe(i, j+1, n)] = 1/Y
            matrix[k,directe(i, j-1, n)] = 1/Y
            matrix[k,k] = -2*((1/Y)+(1/X))
            matrix[k,directe(i+1, j, n)+2*n*m] = -1/(C*x)
            matrix[k,directe(i, j, n)+2*n*m] = 1/(C*x)
            LX.append(directe(i,j,n))
            matrix[directe(i,j,n)+2*m*n,directe(i,j,n)+2*m*n]=1
            b[directe(i,j,n)+2*m*n,0]=P1
            LZ.append(directe(i,j,n)+2*m*n)
            matrix[directe(i,j,n)+m*n,directe(i,j,n)+m*n]=1
            LY.append(directe(i,j,n)+m*n)

        if(i==n-1 ):#P2 DERIVE DE PERIDICITE UY=0
            matrix[k,directe(1, j, n)] = 1/X
            matrix[k,directe(i-1, j, n)] = 1/X
            matrix[k,directe(i, j+1, n)] = 1/Y
            matrix[k,directe(i, j-1, n)] = 1/Y
            matrix[k,k] = -2*((1/Y)+(1/X))
            matrix[k,directe(i, j, n)+2*n*m] = -1/(C*x)
            matrix[k,directe(i-1, j, n)+2*n*m] = 1/(C*x)
            LX.append(directe(i,j,n))
            matrix[directe(i,j,n)+m*n,directe(i,j,n)+m*n]=1
            LY.append(directe(i,j,n)+m*n)
            matrix[directe(i,j,n)+2*m*n,directe(i,j,n)+2*m*n]=1
            b[directe(i,j,n)+2*m*n,0]=P2
            LZ.append(directe(i,j,n)+2*m*n)
            # peridicite + laplacien de u en haut et bas
        if(j==m-1 and i!=0 and i!=n-1):
            matrix[k+2*m*n,directe(i+1, j, n)] = 1/(2*x)
            matrix[k+2*m*n,directe(i-1, j, n)] = -1/(2*x)
            matrix[k+2*m*n,directe(i, j, n)+m*n] = 1/y
            matrix[k+2*m*n,directe(i, j-1, n)+m*n] = -1/y
            LZ.append(directe(i,j,n)+2*m*n)
        if(j==0 and i!= n-1 and i!=0 ):
            matrix[k+2*m*n,directe(i-1, j, n)] = -1/(2*x)
            matrix[k+2*m*n,directe(i+1, j, n)] = 1/(2*x)
            matrix[k+2*m*n,directe(i, j+1, n)+m*n] = 1/y
            matrix[k+2*m*n,directe(i, j, n)+m*n] = -1/y
            LZ.append(directe(i,j,n)+2*m*n)

# remplissage de l interieur
for k in range(m * n):
  if inverse(k, n) not in Q :
    if k not in LX  and k not in I :
        i=inverse(k, n)[0]
        j=inverse(k, n)[1]
        matrix[k,directe(i+1, j, n)] = 1/X
        matrix[k,directe(i-1, j, n)] = 1/X
        matrix[k,directe(i, j+1, n)] = 1/Y
        matrix[k,directe(i, j-1, n)] = 1/Y
        matrix[k,k] = -2*((1/Y)+(1/X))
        matrix[k,directe(i+1, j, n)+2*n*m] = -1/(2*C*x)
        matrix[k,directe(i-1, j, n)+2*n*m] = 1/(2*C*x)
    w=k+m*n
    if w not in LY  and k not in I :
        i=inverse(k, n)[0]
        j=inverse(k, n)[1]
        matrix[k+m*n,directe(i+1, j, n)+m*n] = 1/X
        matrix[k+m*n,directe(i-1, j, n)+m*n] = 1/X
        matrix[k+m*n,directe(i, j+1, n)+m*n] = 1/Y
        matrix[k+m*n,directe(i, j-1, n)+m*n] = 1/Y
        matrix[k+m*n,k+m*n] = -2*((1/Y)+(1/X))
        matrix[k+m*n,directe(i, j+1, n)+2*n*m] = -1/(2*C*y)
        matrix[k+m*n,directe(i, j-1, n)+2*n*m] = 1/(2*C*y)
    v=k+2*m*n
    if v not in LZ  and k not in I  :
        i=inverse(k, n)[0]
        j=inverse(k, n)[1]
        matrix[k+2*m*n, directe(i+1, j, n)] = 1 / (2*x)
        matrix[k+2*m*n, directe(i-1, j, n)] = -1 / (2*x)
        matrix[k+2*m*n, directe(i, j+1, n)+m*n] = 1 / (2*y)#fix
        matrix[k+2*m*n, directe(i, j-1, n)+m*n] = -1 / (2*y)#fix
#remplissage des coins
for i in I:
    print(i)
    matrix[i,i]=1
    matrix[i+m*n,i+m*n]=1
    matrix[i+2*m*n,i+2*m*n]=1
    b[0+2*m*n,0]=P1
    b[m*n-n+2*m*n,0]=P1
    b[n-1+2*m*n,0]=P2
    b[m*n-1+2*m*n,0]=P2

# la deformation
for k in range(m*n):
    i=inverse(k, n)[0]
    j=inverse(k, n)[1]
    if inverse(k, n) in Q:
            matrix[k,k] = 1
            matrix[k+m*n,k+m*n] = 1
        # l interieur ou le tous est null
    if j < ymax and xmin<i<xmax:
            matrix[directe(i,j,n)+2*m*n,directe(i,j,n)+2*m*n] = 1
    if j == ymax and  xmin<i<xmax:
            matrix[k+2*m*n,directe(i,j+1,n)+m*n] = 1
            matrix[k+2*m*n,directe(i,j,n)+m*n] = -1
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

matrix = csc_matrix(matrix)

# Résolution du système
solution = spsolve(matrix, b)

# Extraction de Ux, Uy et P de la solution
Ux = solution[:m * n].reshape((m, n))
Uy = solution[m * n:2 * m * n].reshape((m, n))
P = solution[2 * m * n:].reshape((m, n))
ux = solution[:m * n]
uy = solution[m * n:2 * m * n]
p = solution[2 * m * n:]
U = []


# Définition de la fonction de débit
def debit(u, n, m, l):
    D = []
    for i in range(n):
        Um = 0
        for j in range(m):
            Um += u[i + j * n]
        D.append(Um * y)
    return D

# Calcul du débit
D = debit(ux, n, m, l)
ecart_type = np.std(D)
print(ecart_type)
debitm = sum(D) / len(D)
print(debitm)
# Création des points sur l'axe X
x = np.linspace(0, L, n)

# Tracé du débit (figure séparée)
plt.figure(figsize=(10, 6))
plt.plot(x, D, marker='o', label='Débit')
plt.title("Débit en fonction de la position x")
plt.xlabel("Position (x)")
plt.ylabel("Débit (D)")
plt.text(0.0075, debitm, f"Écart-type : {ecart_type:.8f}", fontsize=12, color='blue')
plt.grid(True)
plt.legend()
plt.show()
# Création d'une figure avec un agencement de grille
fig, axs = plt.subplot_mosaic([['Ux', 'Uy', 'P']], figsize=(18, 6))

# Affichage de Ux
img1 = axs['Ux'].imshow(Ux.T, cmap='viridis', aspect='auto', extent=[0,L,0,l])
fig.colorbar(img1, ax=axs['Ux'])
axs['Ux'].set_title('Variation de Ux')
axs['Ux'].set_xlabel('Index Y')
axs['Ux'].set_ylabel('Index X')
# Affichage de Uy
img2 = axs['Uy'].imshow(Uy.T, cmap='viridis', aspect='auto', extent=[0,L,0,l])
fig.colorbar(img2, ax=axs['Uy'])
axs['Uy'].set_title('Variation de Uy')
axs['Uy'].set_xlabel('Index Y')
axs['Uy'].set_ylabel('Index X')
cmap = plt.cm.get_cmap('viridis')
cmap.set_under(color='red')
P_min = 10.5
P_max = 10
img3 = axs['P'].imshow(P.T, cmap=cmap, aspect='auto', vmin=P_min, vmax=P_max , extent=[0,L,0,l])
fig.colorbar(img3, ax=axs['P'])
axs['P'].set_title('Pression P')
axs['P'].set_xlabel('Index Y')
axs['P'].set_ylabel('Index X')
np.set_printoptions(linewidth=200, threshold=np.inf, suppress=True)
plt.show()
