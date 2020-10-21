#building the CE model using linear regression method (other methods are available too)
from clusterx.model import ModelBuilder
from clusterx.clusters.clusters_pool import ClustersPool
from clusterx.structures_set import StructuresSet
from clusterx.super_cell import SuperCell

# very inconsistent keyword naming scheme for retrieving data from .json files
sset = StructuresSet(db_fname="sset.json")
cpool = ClustersPool(json_db_filepath="cpool.json")

mb = ModelBuilder(selector_type="linreg",
                  selector_opts={'clusters_sets':'size'},
                  estimator_type="skl_LinearRegression",
                  estimator_opts={"fit_intercept":False})
cemodel = mb.build(sset, cpool, "energy") # Build CE model using the training data set
cpool_opt1 = mb.get_opt_cpool()

cemodel.report_errors(sset)
cpool_opt1.display_info(ecis=cemodel.get_ecis())
cpool_opt1.write_clusters_db(db_name="cpool_opt.json")

# Show sublattice properties
scell = SuperCell(json_db_filepath="scell.json")
scell.get_sublattice_types(pretty_print=True)
sites_dict = scell.get_nsites_per_type()
for key in sites_dict.keys():
    print("Number of atoms in sublattice "+str(key)+":", sites_dict[key])

nsites = len(scell.get_substitutional_atoms())

from clusterx.thermodynamics.monte_carlo import MonteCarlo
from clusterx.model import Model
from ase.io.trajectory import Trajectory
from ase.io import write, read
kb = float(8.6173303*10**(-5)) # Boltzmann constant in eV/K
temp = 300 # Temperature in K

for n in range(1, len(scell)):
    nsubs = {0:[n]}

    # Initialization of a MonteCarlo object

    mc = MonteCarlo(cemodel,
                    scell,
                    ensemble = "canonical",
                    nsubs=nsubs,
                    predict_swap = True,
                    filename = 'MC'+str(nsubs)+'.json')

    # Execution of a Metropolis Monte-Carlo sampling
    traj = mc.metropolis(no_of_sampling_steps = 1000,
                        temperature = 500, # Kelvin
                        boltzmann_constant = kb,
                        scale_factor = [1/(1.0*nsites)])
    # do we need to run DFT on MC structures?
    lowest_energy = traj.get_lowest_energy_structure() # get lowest energy non duplicated
    sset.add_structure(lowest_energy, write_db=True)      #adding the structures to the sset? supercell has no add structure

# Write out the new folder structure and generate a file cointaining the list of folders
sset.write_files(prefix="sset_mc")


# Write the database once all MonteCarlo simulations have finished
sset.serialize("sset_mc.json")

structure_locations = sset.get_folders()
sset_write_locations = open("folder_paths.txt", "w")
for i in structure_locations:
    if i is not structure_locations[-1]:
        sset_write_locations.write(i+" ")
    else:
        sset_write_locations.write(i) # avoid ending file with the denominator

sset_write_locations.close()

# Go back to step 3 to recalculate the energies
# And then to step 4 to generate the updated plots
