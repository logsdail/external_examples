from ase.build import bulk
from clusterx.parent_lattice import ParentLattice
from clusterx.structures_set import StructuresSet
from clusterx.super_cell import SuperCell
from clusterx.clusters.clusters_pool import ClustersPool
from random import randint
import numpy as np
import os

# Generate Parent lattice and supercell
pri = bulk("Pd", a=3.89)
list_of_elements = [["Pd", "Cu"]] * len(pri)
platt = ParentLattice(pri, symbols=list_of_elements)
scell = SuperCell(platt, [2,2,2])
sset = sset = StructuresSet(platt)

# Generate random structures
nstruc = 10
for i in range(nstruc):
    concentration = {0:[randint(1,len(scell))]}
    sset.add_structure(scell.gen_random(concentration))

# Generate pool of clusters
cpool = ClustersPool(platt, npoints=[0,1,2], radii=[0, -1,-1], super_cell=scell)

# Calculate EMT energy
from ase.calculators.emt import EMT
sset.set_calculator(EMT())
sset.calculate_property()

# Build the training model
from clusterx.model import ModelBuilder, Model
mb = ModelBuilder(selector_type="linreg",
                  selector_opts={'clusters_sets':'size'},
                  estimator_type="skl_LinearRegression",
                  estimator_opts={"fit_intercept":False})

cemodel = mb.build(sset, cpool, "energy") # This returns a Model class object
cemodel.serialize("cemodel.json") # Model object serialization
################ DEBUG ##################################
# TODO: The MonteCarlo does not work when CE_model is read from .json
# When the line 44 is commented, the script works as normal
# (i.e. Model generated in the same script as MC or predictions, inefficient)
cemodel = Model(json_db_filepath="cemodel.json") # Model object reinitialisation from .json
#########################################################

# Energy predictions do not work also
print("Energy predictions:", sset.get_predictions(cemodel))

# Perform MonteCarlo
from clusterx.thermodynamics.monte_carlo import MonteCarlo
from ase.io.trajectory import Trajectory
from ase.io import write, read
kb = float(8.6173303*10**(-5)) # Boltzmann constant in eV/K
temp = 300 # Temperature in K

# Analyse for 4 substitutions
nsubs = {0:[4]}
nsites = len(scell.get_substitutional_atoms())
# Initialization of a MonteCarlo object
mc = MonteCarlo(cemodel,
                scell,
                ensemble = "canonical",
                nsubs=nsubs,
                predict_swap = True,
                filename = 'MC'+str(nsubs)+'.json')

# Execution of a Metropolis Monte-Carlo sampling
traj = mc.metropolis(no_of_sampling_steps = 100,
                    temperature = 500, # Kelvin
                    boltzmann_constant = kb,
                    scale_factor = [1/(1.0*nsites)])
lowest_energy = traj.get_lowest_energy_structure()
sset.add_structure(lowest_energy, write_db=True)
