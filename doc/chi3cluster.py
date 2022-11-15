#
# \chi^{(3)} of a Hubbard cluster calculated using pycommute
#

from itertools import product
import numpy as np

from pycommute.expression import c, c_dag, n
from pycommute.models import tight_binding, dispersion, hubbard_int
from pycommute.loperator import HilbertSpace, LOperatorR, make_matrix

from networkx.generators.lattice import grid_2d_graph
from networkx.linalg.graphmatrix import adjacency_matrix

# 3x1 Hubbard plaquette
Nx, Ny = 3, 1

# Parameters
beta = 5.0
t = 1.0
U = 4.0
mu = 0.6 * U

# Lattice and hopping matrix
lat = grid_2d_graph(Nx, Ny, periodic=True)
hopping_matrix = -t * adjacency_matrix(lat).todense()

# Lists of indices for electronic spin-up and spin-down operators.
indices_up = [(x, y, "up") for x, y in lat.nodes()]
indices_dn = [(x, y, "dn") for x, y in lat.nodes()]
indices = indices_up + indices_dn

# Hamiltonian of the model

H = tight_binding(hopping_matrix, indices=indices_up) \
  + tight_binding(hopping_matrix, indices=indices_dn)
H += dispersion(-mu * np.ones(Nx*Ny), indices=indices_up) \
  + dispersion(-mu * np.ones(Nx*Ny), indices=indices_dn)
H += hubbard_int(U * np.ones(Nx*Ny),
                 indices_up=indices_up, indices_dn=indices_dn)

# Hilbert space
hs = HilbertSpace(H)

# Operator form of the Hamiltonian
H_op = LOperatorR(H, hs)

# Matrix form of the Hamiltonian
H_mat = np.matrix(make_matrix(H_op, hs))
np.testing.assert_allclose(H_mat.H, H_mat)

def to_eigenbasis(O):
    return U_mat.H @ O @ U_mat

# Diagonalization
E, U_mat = np.linalg.eigh(H_mat)
print("E =", E)
np.testing.assert_allclose(U_mat.H @ U_mat, np.eye(hs.dim), atol=1e-14)
np.testing.assert_allclose(to_eigenbasis(H_mat), np.diag(E), atol=1e-13)

w = np.exp(-beta * E) # Statistical weights of levels
Z = np.sum(w) # Partition function
w /= Z

# Construct matrix representation of creation and annihilation operators
C_mats = {ind: to_eigenbasis(make_matrix(LOperatorR(c(*ind), hs), hs)) for ind in indices}
Cdag_mats = {ind: to_eigenbasis(make_matrix(LOperatorR(c_dag(*ind), hs), hs)) for ind in indices}
for ind in indices:
    np.testing.assert_allclose(C_mats[ind].H, Cdag_mats[ind], atol=1e-13)

# Function f_{234} from Dominik Kiese's notes
def f(i, j, k, w1, w2):
    res = (w[j] + w[i]) / (E[i] - E[j] - 1j*w1)
    if np.isclose(E[i], E[k], atol=1e-13):
        res += beta * w[i] * np.isclose(w1, -w2, atol=1e-13)
    else:
        res += (w[k] - w[i]) / (E[i] - E[k] - 1j*w1 - 1j*w2)
    res /= -(E[j] - E[k] - 1j*w2)
    return res

# \chi^{(3)}_{pp}
def chi3_pp(x1p, x1, x2p, x2, w1, w2):
    Delta = C_mats[x1] @ C_mats[x2]
    res = 0
    for i, j, k in product(range(hs.dim), repeat=3):
        res += f(i, j, k, w1, w2) * Cdag_mats[x1p][i, j] * Cdag_mats[x2p][j, k] * Delta[k, i]
        res += -f(i, j, k, w2, w1) * Cdag_mats[x2p][i, j] * Cdag_mats[x1p][j, k] * Delta[k, i]
    return res

# \chi^{(3)}_{ph}
def chi3_ph(x1p, x1, x2p, x2, w1, w2):
    N = Cdag_mats[x2p] @ C_mats[x2]
    res = 0
    for i, j, k in product(range(hs.dim), repeat=3):
        res += -f(i, j, k, w1, -w2) * Cdag_mats[x1p][i, j] * C_mats[x1][j, k] * N[k, i]
        res += f(i, j, k, -w2, w1) * C_mats[x1][i, j] * Cdag_mats[x1p][j, k] * N[k, i]
    return res

# \chi^{(3)}_{xph}
def chi3_xph(x1p, x1, x2p, x2, w1, w2):
    barN = -Cdag_mats[x2p] @ C_mats[x1]
    res = 0
    for i, j, k in product(range(hs.dim), repeat=3):
        res += -f(i, j, k, w1, -w2) * Cdag_mats[x1p][i, j] * C_mats[x2][j, k] * barN[k, i]
        res += f(i, j, k, -w2, w1) * C_mats[x2][i, j] * Cdag_mats[x1p][j, k] * barN[k, i]
    return res

# Matsubara frequency
def nu(n):
    return np.pi*(2*n+1)/beta

print("Particle-particle")
x1p = (0, 0, "up")
x1 = (2, 0, "up")
x2p = (0, 0, "dn")
x2 = (2, 0, "dn")

for n1, n2 in product(range(-1, 2), repeat=2):
    print(n1, n2, chi3_pp(x1p, x1, x2p, x2, nu(n1), nu(n2)))

print("Particle-hole")
x1p = (0, 0, "up")
x1 = (0, 0, "up")
x2p = (2, 0, "dn")
x2 = (2, 0, "dn")

for n1, n2 in product(range(-1, 2), repeat=2):
    print(n1, n2, chi3_ph(x1p, x1, x2p, x2, nu(n1), nu(n2)))

print("Crossed particle-hole")
x1p = (0, 0, "up")
x1 = (2, 0, "up")
x2p = (2, 0, "dn")
x2 = (0, 0, "dn")

for n1, n2 in product(range(-1, 2), repeat=2):
    print(n1, n2, chi3_xph(x1p, x1, x2p, x2, nu(n1), nu(n2)))
