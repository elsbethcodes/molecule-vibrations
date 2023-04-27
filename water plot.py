""" This code plots the vibrational modes of the water molecule in 3D. The eigenvector you wish to plot must be uncommented first, one at a time. """
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

# # Sets up the figure with an axis.
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Defines the coordinates of the oxygen and hydrogen atoms.
oxygen = [0, 0, 0]
hydrogen1 = [0, -np.sqrt(3)/2, -0.5]
hydrogen2 = [0, np.sqrt(3)/2, -0.5]

# Defines the eigenvectors for all three vibrational modes of water.
# Uncomment the mode you'd like to plot.
# c1
Vector_o = [0, 0, 2]
Vector_h1 = [0, 1.851366, -1]
Vector_h2 = [0, -1.851366, -1]
# c2
"""Vector_o = [0, 0, 2]
Vector_h1 = [0, -1.620425, -1]
Vector_h2 = [0, 1.620425, -1]"""
# c3
"""Vector_o = [0, -0.774597, 0]
Vector_h1 = [0, 0.774597/2, 0.223607]
Vector_h2 = [0, 0.774597/2, -0.223607]"""

# Normalises (and halves the magnitude of) the eigenvectors. Note that the normalised (and shortened) eigenvectors have lowercase v.
vector_o = normalise(Vector_o)
vector_h1 = normalise(Vector_h1)
vector_h2 = normalise(Vector_h2)

# Plots the atoms where Oxygen is in red 2 times the size of the Hydrogens in white. Note that these are toy sizes, not realistic.
ax.scatter(*oxygen, color='r', s=100)
ax.scatter(*hydrogen1, color='w', s=50)
ax.scatter(*hydrogen2, color='w', s=50)

# Plots the eigenvectors as green arrows.
ax.quiver(oxygen[0], oxygen[1], oxygen[2], vector_o[0], vector_o[1], vector_o[2], color='g', arrow_length_ratio=0.2)
ax.quiver(hydrogen1[0], hydrogen1[1], hydrogen1[2], vector_h1[0], vector_h1[1], vector_h1[2], color='g', arrow_length_ratio=0.2)
ax.quiver(hydrogen2[0], hydrogen2[1], hydrogen2[2], vector_h2[0], vector_h2[1], vector_h2[2], color='g', arrow_length_ratio=0.2)

# Draws the bonds as black lines between atoms.
ax.plot([oxygen[0], hydrogen1[0]], [oxygen[1], hydrogen1[1]], [oxygen[2], hydrogen1[2]], color='k')
ax.plot([oxygen[0], hydrogen2[0]], [oxygen[1], hydrogen2[1]], [oxygen[2], hydrogen2[2]], color='k')

# Adjusts the axes.
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
ax.set_zlim(-1, 1)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Removes the axis numbering and grid (feel free to undo).
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
ax.grid(False)

# Reveals the plot.
plt.show()
