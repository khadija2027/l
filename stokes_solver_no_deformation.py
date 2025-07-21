import numpy as np
import matplotlib.pyplot as plt

# Dimensions et paramètres
n, m, L, l = 50, 50,  0.0196, 0.00655 #L et l en m
x = L / n
y = l / m
C = 0.000018  # C est la viscosité dynamique en Pa.s
P1 = 10.5
P2 = 10    #P1 et P2 en Pa
Y = y**2
X = x**2

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


# Création d'une figure avec un agencement de grille
fig, axs = plt.subplot_mosaic([['Ux', 'Uy', 'P']], figsize=(18, 6))

# Affichage de Ux
img1 = axs['Ux'].imshow(Ux, cmap='viridis', aspect='auto', extent=[0, L, 0, l])

fig.colorbar(img1, ax=axs['Ux'])
axs['Ux'].set_title('Variation de Ux(m/s)')
axs['Ux'].set_xlabel('Index X(m)')
axs['Ux'].set_ylabel('Index Y(m)')

# Affichage de Uy
img2 = axs['Uy'].imshow(Uy, cmap='viridis', aspect='auto', extent=[0, L, 0, l])

fig.colorbar(img2, ax=axs['Uy'])
axs['Uy'].set_title('Variation de Uy(m/s)')
axs['Uy'].set_xlabel('Index X(m)')
axs['Uy'].set_ylabel('Index Y(m)')

# Affichage de la pression P
img3 = axs['P'].imshow(P, cmap='viridis', aspect='auto', extent=[0, L, 0, l])

fig.colorbar(img3, ax=axs['P'])
axs['P'].set_title('Pression P(Pa)')
axs['P'].set_xlabel('Index X(m)')
axs['P'].set_ylabel('Index Y(m)')



np.set_printoptions(linewidth=200, threshold=np.inf, suppress=True)
plt.show()
