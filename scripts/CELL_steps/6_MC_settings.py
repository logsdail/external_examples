# Starting MC simulations
# Show sublattice properties
from clusterx.model import ModelBuilder, Model
from clusterx.clusters.clusters_pool import ClustersPool
from clusterx.structures_set import StructuresSet
from clusterx.super_cell import SuperCell


# Read the previously generated structures_set from .json database file
sset = StructuresSet(db_fname="sset.json")

################ TODO: Raise with Santiago #################################
cpool = ClustersPool(json_db_filepath="cpool.json")
mb = ModelBuilder(selector_type="linreg",
                  selector_opts={'clusters_sets':'size'},
                  estimator_type="skl_LinearRegression",
                  estimator_opts={"fit_intercept":False})
cemodel = mb.build(sset, cpool, "energy") # Build CE model using the training data set
cpool_opt1 = mb.get_opt_cpool()
############################################################################
#cemodel = Model(json_db_filepath="CE_model.json") # TODO: Calling a serialized model does not work

scell = SuperCell(json_db_filepath="scell.json")
scell.get_sublattice_types(pretty_print=True)
sites_dict = scell.get_nsites_per_type()
for key in sites_dict.keys():
    print("Number of atoms in sublattice "+str(key)+":", sites_dict[key])

nsites = len(scell.get_substitutional_atoms())
print("len scell", len(scell))

from clusterx.thermodynamics.monte_carlo import MonteCarlo
from ase.io.trajectory import Trajectory
from ase.io import write, read
kb = float(8.6173303*10**(-5)) # Boltzmann constant in eV/K
temp = 300 # Temperature in K

for n in range(1, len(scell)):
    print("Now analysing systems with", n, "substitutions.")
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
# Include previously calculated values in energy.dat file
for i in range(0, len(sset)):
        property_name = "energy"
        e_model = sset.get_property_values(property_name)[i]
        f = open(sset.get_folders()[i]+"/"+property_name+".dat", "w")
        f.write(str(e_model))
        f.close()



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
# Then update CE model in step 5
