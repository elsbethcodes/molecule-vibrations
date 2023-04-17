# Force matrix calculator

# SymPy enables symbolic manipulation
import sympy as sp
from sympy import *
import pprint as pp


# R = equilibrium coordinates for n atoms.
"""H2O: central atom O is atom 1, alpha is the strength of O-H and beta is the strength of H-H."""
mol = 'H2O'
R = (0,0,0,0,-sp.sqrt(3)/2,-1/2,0,sp.sqrt(3)/2,-1/2)
centre=1
"""NH3: central atom N is atom 4, alpha is the strength of H-H and beta is the strength of N-H."""
#mol = 'NH3'
#R = [1,0,0,-1/2,sp.sqrt(3)/2,0,-1/2,-sp.sqrt(3)/2,0,0,0,sp.sqrt(3)]
#centre=4

#n is the number of the atoms and is found by dividing the number of equilibrium coordinates into 3.
n = int(len(R)/3)

# d is a dictionary of equilibrium coordinates for each atom.
d = {}
# R1 is the first 3 values of R, R2 is the second 3 values, and so on.
for x in range(0,n):
    d["R{}".format(x+1)] = Matrix(R[3*x:3*x+3])

# |Ri - Rj|
def length(i,j):
    diff = d["R{}".format(i)] - d["R{}".format(j)]
    sumofsquares = 0
    for x in diff:
        sumofsquares += x**2
    return sp.sqrt(sumofsquares)

# Spring constants function.
A = sp.Symbol('a')
B = sp.Symbol('b')
def spring(i,j):
    if i == centre or j == centre:
        if mol == 'H2O':
            return A
        if mol == 'NH3':
            return B
    else:
        if mol == 'H2O':
            return B
        if mol == 'NH3':
            return A    

# Creates an empty 3n dimensional force matrix.
F = sp.zeros(3*n,3*n)

for i in range(1,n+1):
    for j in range(1, n+1):
        if i <= j:
            for a in range(1,4):
                for b in range (1,4):
                    # Row and column numbers of F based on atoms i,j and dimensions a,b.
                    r = (3*i -3 + a) -1
                    c = (3*j -3 + b) -1

                    # Indexing for Ria.
                    # The -1 is because python indexing starts at 0 but we start labelling at 1.
                    ia = (3*i -3 + a) -1
                    ja = (3*j -3 + a) -1
                    ib = (3*i -3 + b) -1
                    jb = (3*j -3 + b) -1

                    if r<=c:
                        if i == j:
                            ans = 0
                            for p in range(1,n+1):
                                if i != p:
                                    ans += spring(i,p)*(R[ia]-R[3*p -3 + a -1])*(R[ib]-R[3*p -3 + b -1])/length(i,p)**2
                            F[c,r] = ans
                            F[r,c] = ans
                        else:
                            ans = -spring(i,j)*(R[ia]-R[ja])*(R[ib]-R[jb])/length(i,j)**2
                            F[c,r] = ans
                            F[r,c] = ans

print(F)

