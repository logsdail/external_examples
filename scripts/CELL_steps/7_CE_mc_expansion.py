from clusterx.structures_set import StructuresSet
from clusterx.clusters.clusters_pool import ClustersPool
from clusterx.super_cell import SuperCell
from clusterx.parent_lattice import ParentLattice
from clusterx.model import ModelBuilder

platt = ParentLattice(json_db_filepath="platt.json")
# Define a new, larger super cell
scell = SuperCell(platt, [3, 3, 7])

# Build a Cluster Expansion model based on structures from previous calculations
# and an optimised pool of clusters
sset = StructuresSet(db_fname="sset_mc.json")
cpool = ClustersPool(json_db_filepath="cpool.json") # can't use cpool_opt?

mb = ModelBuilder(selector_type="linreg",
                  selector_opts={'clusters_sets':'size'},
                  estimator_type="skl_LinearRegression",
                  estimator_opts={"fit_intercept":False})

cemodel = mb.build(sset, cpool, "energy") # Build CE model using the training data set
cpool_opt_mc = mb.get_opt_cpool()
cemodel.report_errors(sset)
cpool_opt_mc.display_info(ecis=cemodel.get_ecis())
cpool_opt_mc.write_clusters_db(db_name="cpool_opt_mc.json")


nsites = len(scell.get_substitutional_atoms())

from clusterx.thermodynamics.monte_carlo import MonteCarlo
from clusterx.model import Model
from ase.io.trajectory import Trajectory
from ase.io import write, read
kb = float(8.6173303*10**(-5)) # Boltzmann constant in eV/K
temp = 300 # Temperature in K

list_of_mc_structures = []

for n in range(1, nsites, int(nsites/10)):
    print("Now analysing systems with", n, "substitutions.")
    nsubs = {0:[n]}

    # Initialization of a MonteCarlo object

    mc = MonteCarlo(cemodel,
                    scell,
                    ensemble = "canonical",
                    nsubs=nsubs,
                    predict_swap = True,
                    filename = 'MC_CE'+str(nsubs)+'.json')

    # Execution of a Metropolis Monte-Carlo sampling
    traj = mc.metropolis(no_of_sampling_steps = 100, # increase if necessary
                        temperature = 500, # Kelvin
                        boltzmann_constant = kb,
                        scale_factor = [1/(1.0*nsites)])
    # do we need to run DFT on MC structures?
    lowest_energy = traj.get_lowest_energy_structure() # get lowest energy non duplicated
    list_of_mc_structures += [lowest_energy]
    sset.add_structure(lowest_energy, write_db=True)      # adding the structures to the sset? supercell has no add structure

# Write the database once all MonteCarlo simulations have finished
sset.serialize("sset_ce.json")

# Write out the new folder structure and generate a file cointaining the list of folders
sset.write_files(prefix="sset_ce")
# Plot the data
from ase.visualize import view
view(list_of_mc_structures)
refs = StructuresSet(db_fname="refs.json")
from clusterx.visualization import plot_property_vs_concentration
property_name = "energy"
plot_property_vs_concentration(sset, site_type=0, property_name=property_name,refs=ref_en,scale=0.6)
