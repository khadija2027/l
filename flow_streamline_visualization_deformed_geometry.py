# Remplissage des conditions pour Ux
import numpy as np
import matplotlib.pyplot as plt

# Dimensions et paramètres
n, m, L, l = 50, 50 , 0.0196 ,0.00655
x = L / (n-1)
y = l / (m-1)
C = 0.000018
P1 = 10.5
P2 = 10
Y = y**2
X = x**2


# Listes pour suivre les indices de remplissage
LX = []
LY = []
LZ = []

def rotate_matrix_90_antihour(matrix):
    # Effectue une rotation de 90° antihoraire
    return matrix[::-1 ]

fig, axes = plt.subplots(1, 3, figsize=(15, 4))
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

xmin, xmax = 17,33
J = [15,30,45]

for idx, ymax in enumerate(J):
    Q = defo(xmin, xmax, ymax)
    matrix = np.zeros((3 * m * n, 3 * m * n))
    b = np.zeros((3 * m * n,1 ))
    def directe(i, j, n):
        return (j)*n + i

    def inverse(k, n):
        return (k % n, k // n)
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
                matrix[k][directe(i+1, j, n)] = 1/X
                matrix[k][directe(n-2, j, n)] = 1/X
                matrix[k][directe(i, j+1, n)] = 1/Y
                matrix[k][directe(i, j-1, n)] = 1/Y
                matrix[k][k] = -2*((1/Y)+(1/X))
                matrix[k][directe(i+1, j, n)+2*n*m] = -1/(C*x)
                matrix[k][directe(i, j, n)+2*n*m] = 1/(C*x)
                LX.append(directe(i,j,n))
                matrix[directe(i,j,n)+2*m*n,directe(i,j,n)+2*m*n]=1
                b[directe(i,j,n)+2*m*n]=P1
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
     if inverse(k, n) not in Q :
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
                matrix[k][k] = 1
                matrix[k+m*n][k+m*n] = 1
            # l interieur ou le tous est null
        if j < ymax and xmin<i<xmax:
                matrix[directe(i,j,n)+2*m*n][directe(i,j,n)+2*m*n] = 1
        if j == ymax and  xmin<i<xmax:
                matrix[k+2*m*n][directe(i,j+1,n)+m*n] = 1
                matrix[k+2*m*n][directe(i,j,n)+m*n] = -1
        if i==xmin and 0<j<ymax:
                matrix[k+2*m*n][directe(i,j,n)] =1
                matrix[k+2*m*n][directe(i-1,j,n)] = -1
        if i==xmax and 0<j<ymax:
                matrix[k+2*m*n][directe(i,j,n)] =-1
                matrix[k+2*m*n][directe(i+1,j,n)] = 1
        if (i==xmin and j==ymax ):
        # (Ux(k-1)-Ux(k+1))/dx+(Uy(k-N)-Uy(k+N))/dy=0 le divergent
                matrix[k+2*m*n][directe(i-1,j,n)] = 1/(2*x)
                matrix[k+2*m*n][directe(i+1,j,n)] = - 1/(2*x)
                matrix[k+2*m*n][directe(i,j-1,n)+m*n] = 1/(2*y)
                matrix[k+2*m*n][directe(i,j+1,n)+m*n] = -1/(2*y)
        if ( i==xmax and j==ymax):
                matrix[k+2*m*n][directe(i-1,j,n)] = 1/(2*x)
                matrix[k+2*m*n][directe(i+1,j,n)] = - 1/(2*x)
                matrix[k+2*m*n][directe(i,j-1,n)+m*n] = 1/(2*y)
                matrix[k+2*m*n][directe(i,j+1,n)+m*n] = -1/(2*y)
        if (i==xmin and j==0 ):
                matrix[k+2*m*n][directe(i-1,j,n)+2*m*n] = 1/X
                matrix[k+2*m*n][directe(i,j+1,n)+2*m*n] =  1/Y
                matrix[k+2*m*n][directe(i,j,n)+2*m*n] = -1*((1/Y)+(1/X))
                #matrix[k+2*m*n][directe(i,j,n)+2*m*n] = 1
                #matrix[k+2*m*n][directe(i-1,j,n)+2*m*n] = -1
        if (i==xmax and j==0 ):
                matrix[k+2*m*n][directe(i+1,j,n)+2*m*n] = 1/X
                matrix[k+2*m*n][directe(i,j+1,n)+2*m*n] =  1/Y
                matrix[k+2*m*n][directe(i,j,n)+2*m*n] = -1*((1/Y)+(1/X))
                #matrix[k+2*m*n][directe(i,j,n)+2*m*n] = 1
                #matrix[k+2*m*n][directe(i+1,j,n)+2*m*n] = -1

    A=matrix
    # Calcul de la solution
    solution, _, _, _ = np.linalg.lstsq(A, b, rcond=None)

    # Extraction de Ux, Uy et P de la solution
    Ux = solution[:m * n].reshape((m, n))
    Uy = solution[m * n:2 * m * n].reshape((m, n))
    P = solution[2 * m * n:].reshape((m, n))

    # Calculer la norme de la vitesse pour le color mapping
    speed = np.sqrt(Ux**2 + Uy**2)
    ax = axes[idx % 3] 
    X_grid, Y_grid = np.meshgrid(np.linspace(0, l, m), np.linspace(0, L, n))
    stream = ax.streamplot(X_grid, Y_grid, Ux, Uy, color='black', linewidth=0.5)

    # Affichage de la vitesse avec imshow
    st = ax.imshow(rotate_matrix_90_antihour(speed), cmap='viridis', aspect='auto', extent=[0, l, 0, L])

    # Configurer le sous-graphe
    ax.set_title(f"Lignes de courant (ymax={ymax})")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    # Ajouter la barre de couleur
    cbar_image = fig.colorbar(st, ax=ax)
    cbar_image.set_label("Vitesse")  # L'étiquette de la barre de couleur
print("Taille de A :", A.shape)
print("Taille de b :", b.shape)
plt.tight_layout()
plt.show()