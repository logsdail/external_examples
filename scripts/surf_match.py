from susmost.cell import enum_cells
from susmost.symmetry import get_2d_sym_ops
from ase.build import fcc111, make_supercell
from ase.io import write, read
from ase.visualize import view
from itertools import product
import numpy as np, sys

def transform_cell(a, T):
  assert T.shape==(3,3)
  b = a.copy()
  b.set_cell((T @ a.cell.T).T, scale_atoms=True)
  return b

def enum_2d_cells(a, max_n):
  symops = get_2d_sym_ops((a.cell, a.get_scaled_positions(), a.numbers))
  res = []
  for i in range(1, max_n+1):
    res = res + [make_supercell(a, c) for c in enum_cells(i, a.cell, symops)]
  return res

a = fcc111('Au', (1,1, 5), vacuum=10)
b = fcc111('Ni', (1,1, 5), vacuum=10)

assert (a.pbc==[True, True, False]).all()
assert (b.pbc==[True, True, False]).all()
assert (a.cell[2, :2]==0).all()
assert (b.cell[2, :2]==0).all()

print (a.cell)
print (b.cell)

a.cell[2,2] = max(a.cell[2,2], b.cell[2,2])
b.cell[2,2] = a.cell[2,2]
#sys.exit()
a_supercells = enum_2d_cells(a, 10)
b_supercells = enum_2d_cells(b, 10)
print ("Supercells #:", len(a_supercells), len(b_supercells))
best_fit_strain = np.inf
for ai,bi in product(a_supercells, b_supercells):
  T = ai.cell.T @ np.linalg.inv(bi.cell.T) # a = T @ b => a @ b^-1 = T
  assert np.linalg.norm(ai.cell.T - T @ bi.cell.T) < 1e-6
  U,S,V = np.linalg.svd(T) # a = U @ S @ V @ b => U^-1 @ a = S @ V @ b => S^-1/2 @ U^-1 @ a = S^1/2 @ V @ b
  max_strain = np.max([S, 1/S])
  if (max_strain < best_fit_strain):
    best_fit_strain = max_strain
    best_pair = (ai, bi, U,S,V,T)

ai, bi, U,S,V,T = best_pair
Uinv = np.linalg.inv(U)
Shalf = S**0.5
Shalf_inv = Shalf**-1.0

ai_1 = transform_cell(ai, Uinv)
ai_2 = transform_cell(ai_1, np.diag(Shalf_inv))

bi_1 = transform_cell(bi, V)
bi_2 = transform_cell(bi_1, np.diag(Shalf))

view([ai, ai_1, ai_2])
view([bi, bi_1, bi_2])


