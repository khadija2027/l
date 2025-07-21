import numpy as np
import matplotlib.pyplot as plt

# Dimensions et paramètres
n, m, L, l = 50, 50,  0.0196, 0.00655
x = L / n
y = l / m
C = 0.000018
Y = y**2
X = x**2
P1 = 10.5
P2 = 10

# Créer la matrice et le vecteur b
matrix = np.zeros((3 * m * n, 3 * m * n))
b = np.zeros((3 * m * n, 1))

# Listes pour suivre les indices de remplissage
LX = []
LY = []
LZ = []

def directe(i, j, n):
    return (j)*n + i

def inverse(k, n):
    return (k % n, k // n)
I = [0, n-1, m*n-n, m*n-1]
for k in range(m*n)  :
    if k not in I:
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
            matrix[k][directe(i+1, j, n)] = 1/X
            matrix[k][directe(n-2, j, n)] = 1/X
            matrix[k][directe(i, j+1, n)] = 1/Y
            matrix[k][directe(i, j-1, n)] = 1/Y
            matrix[k][k] = -2*((1/Y)+(1/X))
            matrix[k][directe(i+1, j, n)+2*n*m] = -1/(C*x)
            matrix[k][directe(i, j, n)+2*n*m] = 1/(C*x)
            LX.append(directe(i,j,n))
            matrix[directe(i,j,n)+2*m*n,directe(i,j,n)+2*m*n]=1
            b[directe(i,j,n)+2*m*n,0]=P1
            LZ.append(directe(i,j,n)+2*m*n)
            matrix[directe(i,j,n)+m*n,directe(i,j,n)+m*n]=1
            LY.append(directe(i,j,n)+m*n)

        if(i==n-1 ):#P2 DERIVE DE PERIDICITE UY=0
            matrix[k][directe(1, j, n)] = 1/X
            matrix[k][directe(i-1, j, n)] = 1/X
            matrix[k][directe(i, j+1, n)] = 1/Y
            matrix[k][directe(i, j-1, n)] = 1/Y
            matrix[k][k] = -2*((1/Y)+(1/X))
            matrix[k][directe(i, j, n)+2*n*m] = -1/(C*x)
            matrix[k][directe(i-1, j, n)+2*n*m] = 1/(C*x)
            LX.append(directe(i,j,n))
            matrix[directe(i,j,n)+m*n,directe(i,j,n)+m*n]=1
            LY.append(directe(i,j,n)+m*n)
            matrix[directe(i,j,n)+2*m*n,directe(i,j,n)+2*m*n]=1
            b[directe(i,j,n)+2*m*n,0]=P2
            LZ.append(directe(i,j,n)+2*m*n)
            # peridicite + laplacien de u en haut et bas
        if(j==0 and i!=0 and i!=n-1):
            matrix[k+2*m*n][directe(i+1, j, n)] = 1/(2*x)
            matrix[k+2*m*n][directe(i-1, j, n)] = -1/(2*x)
            matrix[k+2*m*n][directe(i, j+1, n)+m*n] = 1/y
            matrix[k+2*m*n][directe(i, j, n)+m*n] = -1/y
            LZ.append(directe(i,j,n)+2*m*n)
        if(j==m-1 and i!=0 and i!=n-1):
            matrix[k+2*m*n][directe(i+1, j, n)] = 1/(2*x)
            matrix[k+2*m*n][directe(i-1, j, n)] = -1/(2*x)
            matrix[k+2*m*n][directe(i, j, n)+m*n] = 1/y
            matrix[k+2*m*n][directe(i, j-1, n)+m*n] = -1/y
            LZ.append(directe(i,j,n)+2*m*n)
# remplissage de l interieur
for k in range(m * n):
    if k not in LX  and k not in I :
        i=inverse(k, n)[0]
        j=inverse(k, n)[1]
        matrix[k][directe(i+1, j, n)] = 1/X
        matrix[k][directe(i-1, j, n)] = 1/X
        matrix[k][directe(i, j+1, n)] = 1/Y
        matrix[k][directe(i, j-1, n)] = 1/Y
        matrix[k][k] = -2*((1/Y)+(1/X))
        matrix[k][directe(i+1, j, n)+2*n*m] = -1/(2*C*x)
        matrix[k][directe(i-1, j, n)+2*n*m] = 1/(2*C*x)
    w=k+m*n
    if w not in LY  and k not in I :
        i=inverse(k, n)[0]
        j=inverse(k, n)[1]
        matrix[k+m*n][directe(i+1, j, n)+m*n] = 1/X
        matrix[k+m*n][directe(i-1, j, n)+m*n] = 1/X
        matrix[k+m*n][directe(i, j+1, n)+m*n] = 1/Y
        matrix[k+m*n][directe(i, j-1, n)+m*n] = 1/Y
        matrix[k+m*n][k+m*n] = -2*((1/Y)+(1/X))
        matrix[k+m*n][directe(i, j+1, n)+2*n*m] = -1/(2*C*y)
        matrix[k+m*n][directe(i, j-1, n)+2*n*m] = 1/(2*C*y)
    v=k+2*m*n
    if v not in LZ  and k not in I  :
        i=inverse(k, n)[0]
        j=inverse(k, n)[1]
        matrix[k+2*m*n, directe(i+1, j, n)] = 1 / (2*x)
        matrix[k+2*m*n, directe(i-1, j, n)] = -1 / (2*x)
        matrix[k+2*m*n, directe(i, j+1, n)+m*n] = 1 / y#fix
        matrix[k+2*m*n, directe(i, j, n)+m*n] = -1 / y#fix
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


# Calcul de la solution
solution, _, _, _ = np.linalg.lstsq(matrix, b, rcond=None)

# Extraction de Ux, Uy et P de la solution
Ux = solution[:m * n].reshape((m, n))
Uy = solution[m * n:2 * m * n].reshape((m, n))
P = solution[2 * m * n:].reshape((m, n))

# Définir les dimensions de la grille
x = np.linspace(0, L, n)  # 30 points de 0 à 0.07 (largeur)
y = np.linspace(0, l, m)  # 15 points de 0 à 0.02 (hauteur)
X, Y = np.meshgrid(x, y)

# Calculer la norme de la vitesse pour le color mapping
speed = np.sqrt(Ux**2 + Uy**2)

# Création du graphique
plt.figure(figsize=(8, 6))
plt.quiver(X, Y, Ux, Uy, speed, cmap='viridis')
plt.colorbar(label="Norme de la vitesse")

# Titre et légendes
plt.title("Les vecteurs vitesse pour une pression d'entrée 10.5 et une pression de sortie 10")
plt.xlabel("Position x(m)")
plt.ylabel("Position y(m)")

# Afficher le graphique
plt.show()