from ase.build import fcc111, add_adsorbate # ASE's utilities to build the surface
from clusterx.parent_lattice import ParentLattice
from clusterx.structures_set import StructuresSet
from clusterx.super_cell import SuperCell
from random import randint
import numpy as np
import os
np.random.seed(10)
from ase.build import surface
from ase import Atoms

scell_x, scell_y = 3, 3
layer_height = 7

a = 2.874
c = 3.661

Alloy= Atoms('AuCu', scaled_positions=[(0, 0, 0),
                                (0.5, 0.5, 0.5),
                                ],
              cell=[a, a, c],
              pbc=True)

pri = surface(Alloy, (1, 1, 0), layer_height)
pri.center(vacuum=20.0, axis=2)

symbols = pri.get_chemical_symbols() # Get Chemical symbols of slab
z_coords = pri.get_positions()[:,2] # Get z-coordinate of atomic positions in slab

print("{0:<19s}|{1:<19s}|{2:<19s}".format("Atom index","Chemical symbol","z coordinate")) # Print headers
for i, (symbol, z_coord) in enumerate(zip(symbols,z_coords)):
    print("{0:<19d}|{1:<19s}|{2:<19.3f}".format(i,symbol,z_coord)) # Print atom indexes, symbols and z_coordinat

list_of_elements = [["Au", "Cu"]] * len(pri)
platt = ParentLattice(pri, symbols= list_of_elements)
scell = SuperCell(platt,[scell_x, scell_y])
scell.get_sublattice_types(pretty_print=True)
scell.serialize("scell.json")
sset = StructuresSet(platt)
z_coords_2 = scell.get_positions()[:]

nstruc = 60 # we are generating 60 random structures
for i in range(nstruc):
    concentration = {0:[randint(1,len(scell))]} # Pick a random concentration of "Zn" substitutions starting from 1 to 4*4*3
    sset.add_structure(scell.gen_random(concentration)) # Generate and add a random structure to the StructuresSet

# Write out the file structure and list of file locations
sset.write_files(prefix="sset")
sset.serialize("sset.json") # Write JSON db file for visualization with ASE's GUI.

structure_locations = sset.get_folders()
sset_write_locations = open("folder_paths.txt", "w")
for i in structure_locations:
        sset_write_locations.write(i+" ")
sset_write_locations.close()

# Generate reference structures and the file structure
refs = StructuresSet(platt)
refs.add_structure(scell.gen_random({0:[0]})) # Pristine
refs.add_structure(scell.gen_random({0:[0],0:[len(scell)]})) # Full substitution
refs.write_files(prefix="refs")
refs.serialize("refs.json")

structure_locations = refs.get_folders()
refs_write_locations = open("folder_paths.txt", "a")
for i in structure_locations:
    if i is not structure_locations[-1]:
        refs_write_locations.write(i+" ")
    else:
        refs_write_locations.write(i) # avoid ending file with the denominator

refs_write_locations.close()

from clusterx.clusters.clusters_pool import ClustersPool
import os
if not os.path.exists("cpool.json"):
    cpool = ClustersPool(platt, npoints=[0,1,2], radii=[0, -1,-1], super_cell=scell)
    print(len(cpool)," clusters were generated.")
    cpool.serialize(db_name="cpool.json")

"""
# this can be applied once energies are calculated and retrieved
plot_property_vs_concentration(sset, site_type=0, property_name="total_energy_emt",refs=ref_en,scale=0.6)
"""
