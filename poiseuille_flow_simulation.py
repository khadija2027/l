import matplotlib.pyplot as plt
import numpy as np

# Parameters

N, M, L, l = 50, 50,  0.0196, 0.00655
C = 0.000018
P1 = 10.5
P2 = 10

# Generate spatial grids
x = np.linspace(0, L, N)  # Corrected range for x
y = np.linspace(0, l, M)  # Corrected range for y

# Meshgrid for calculations
X, Y = np.meshgrid(x, y)

# Prepare figure
fig, ax = plt.subplots(1, 3, figsize=(14, 6))

# Ux plot
Ux = ((P2 - P1) / (2*C * L)) * (Y ** 2 - l * Y)
img1 = ax[0].imshow(Ux.T, cmap='viridis', aspect='auto', extent=(0,l, L, 0))
fig.colorbar(img1, ax=ax[0])
ax[0].set_title('Variation de Ux en m/s')
ax[0].set_xlabel('Position Y en m')
ax[0].set_ylabel('Position X en m')

# Uy plot (zero field)
Uy = np.zeros((N, M))  # Matrix of zeros with shape (N, M)
img2 = ax[1].imshow(Uy.T, cmap='viridis', aspect='auto', extent=(0,l, L, 0))
fig.colorbar(img2, ax=ax[1])
ax[1].set_title('Variation de Uy en m/s')
ax[1].set_xlabel('Position Y en m')
ax[1].set_ylabel('Position X en m')

# Pressure (P) plot
P = ((P2 - P1) / L) * X + P1
img3 = ax[2].imshow(P.T, cmap='viridis', aspect='auto', extent=(0,l, L, 0))
fig.colorbar(img3, ax=ax[2])
ax[2].set_title('Variation de la pression en Pa ')
ax[2].set_xlabel( 'Position Y en m')
ax[2].set_ylabel('Position X en m')

plt.tight_layout()
plt.show()