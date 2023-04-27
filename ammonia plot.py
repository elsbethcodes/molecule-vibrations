""" This code plots the vibrational modes of ammonia in 3D. The eigenvector you wish to plot must be uncommented first, one at a time. """

import matplotlib.pyplot as plt
import numpy as np

def normalise(vec):
    """
    Normalises a vector using its norm then halves its magnitude in the interest of the normal modes being
    plotted shorter than the ammonia bonds which are length sqrt(3).
    
    Args:
    vec: a list representing a column vector.
    
    Returns:
    normalised_vec/2: a numpy array representing vec normalised and with magnitude 1/2.
    """
    vec = np.array(vec)
    normalised_vec = vec/np.linalg.norm(vec)
    return normalised_vec/2

# Sets up the figure with an axis.
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Defines the coordinates of the nitrogen and hydrogen atoms.
nitrogen = [0, 0, np.sqrt(3)]
hydrogen1 = [1, 0, 0]
hydrogen2 = [-1/2, np.sqrt(3)/2, 0]
hydrogen3 = [-1/2, -np.sqrt(3)/2, 0]

# Sets the spring constants (a for N-H and b for H-H).
a = 0.05
b = 1

# Defines the normalisation constants in the eigenvectors.
N1 = np.sqrt(144*a**2 - 264*a*b + 169*b**2)
N2 = np.sqrt(144*a**2 - 216*a*b + 121*b**2)

# Defines the eigenvectors for all six vibrational modes of ammonia.
# Uncomment the mode you'd like to plot.
# c1
Vector_n = [np.sqrt(3)*(N1 + 12*a - 11*b)/(18*b), 0, 1/3]
Vector_h1 = [-np.sqrt(3)*(N1 + 12*a - 11*b)/(36*b), (N1 + 12*a - 11*b)/(12*b), 1/3]
Vector_h2 = [-np.sqrt(3)*(N1 + 12*a - 11*b)/(36*b), -(N1 + 12*a - 11*b)/(12*b), 1/3]
Vector_h3 = [0, 0, -1]

# c2
"""Vector_n = [np.sqrt(3)*(-N1 + 12*a - 11*b)/(18*b), 0, 1/3]
Vector_h1 = [-np.sqrt(3)*(-N1 + 12*a - 11*b)/(36*b), (-N1 + 12*a - 11*b)/(12*b), 1/3]
Vector_h2 = [-np.sqrt(3)*(-N1 + 12*a - 11*b)/(36*b), -(-N1 + 12*a - 11*b)/(12*b), 1/3]
Vector_h3 = [0, 0, -1]"""

# c5
"""Vector_n = [-(N2 + 12*a - 7*b)/(6*b), 0, -2*np.sqrt(3)/3]
Vector_h1 = [(N2 + 12*a - 13*b)/(12*b), np.sqrt(3)*(N2 + 12*a -9*b)/(12*b), np.sqrt(3)/3]
Vector_h2 = [(N2 + 12*a -13*b)/(12*b), -np.sqrt(3)*(N2 + 12*a - 9*b)/(12*b), np.sqrt(3)/3]
Vector_h3 = [1, 0, 0]"""

# c9
"""Vector_n = [0, (N2 + 12*a -11*b)/(6*b), 0]
Vector_h1 = [np.sqrt(3)*(N2 + 12*a -9*b)/(12*b), -(N2 + 12*a -5*b)/(12*b), -1]
Vector_h2 = [-np.sqrt(3)*(N2 + 12*a -9*b)/(12*b), -(N2 + 12*a -5*b)/(12*b), 1]
Vector_h3 = [0,1,0]"""

# c6
"""Vector_n = [-(-N2 + 12*a - 7*b)/(6*b), 0, -2*np.sqrt(3)/3]
Vector_h1 = [(-N2 + 12*a - 13*b)/(12*b), np.sqrt(3)*(-N2 + 12*a -9*b)/(12*b), np.sqrt(3)/3]
Vector_h2 = [(-N2 + 12*a -13*b)/(12*b), -np.sqrt(3)*(-N2 + 12*a - 9*b)/(12*b), np.sqrt(3)/3]
Vector_h3 = [1, 0, 0]"""

# c10
"""Vector_n = [0, (-N2 + 12*a -11*b)/(6*b), 0]
Vector_h1 = [np.sqrt(3)*(-N2 + 12*a -9*b)/(12*b), -(-N2 + 12*a -5*b)/(12*b), -1]
Vector_h2 = [-np.sqrt(3)*(-N2 + 12*a -9*b)/(12*b), -(-N2 + 12*a -5*b)/(12*b), 1]
Vector_h3 = [0,1,0]"""

# Normalises (and halves the magnitude of) the eigenvectors. Note that the normalised (and shortened) eigenvectors have lowercase v.
vector_n = normalise(Vector_n)
vector_h1 = normalise(Vector_h1)
vector_h2 = normalise(Vector_h2)
vector_h3 = normalise(Vector_h3)

# Plots the atoms where Nitrogen is in blue 4 times the size of the Hydrogens in white. Note that these are toy sizes, not realistic.
ax.scatter(*nitrogen, color='b', s=200)
ax.scatter(*hydrogen1, color='w', s=50)
ax.scatter(*hydrogen2, color='w', s=50)
ax.scatter(*hydrogen3, color='w', s=50)

# Plots the eigenvectors as green arrows.
ax.quiver(nitrogen[0], nitrogen[1], nitrogen[2], vector_n[0], vector_n[1], vector_n[2], color='g', arrow_length_ratio=0.2)
ax.quiver(hydrogen1[0], hydrogen1[1], hydrogen1[2], vector_h1[0], vector_h1[1], vector_h1[2], color='g', arrow_length_ratio=0.2)
ax.quiver(hydrogen2[0], hydrogen2[1], hydrogen2[2], vector_h2[0], vector_h2[1], vector_h2[2], color='g', arrow_length_ratio=0.2)
ax.quiver(hydrogen3[0], hydrogen3[1], hydrogen3[2], vector_h3[0], vector_h3[1], vector_h3[2], color='g', arrow_length_ratio=0.2)

# Draws the bonds as black lines between atoms. 
ax.plot([nitrogen[0], hydrogen1[0]], [nitrogen[1], hydrogen1[1]], [nitrogen[2], hydrogen1[2]], color='k')
ax.plot([nitrogen[0], hydrogen2[0]], [nitrogen[1], hydrogen2[1]], [nitrogen[2], hydrogen2[2]], color='k')
ax.plot([nitrogen[0], hydrogen3[0]], [nitrogen[1], hydrogen3[1]], [nitrogen[2], hydrogen3[2]], color='k')

# Adjusts the axes.
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_zlim(0, np.sqrt(3) + 0.5)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

# Removes the axis numbering and grid.
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
ax.grid(False)

# Shows the plot.
plt.show()
