# Starting MC simulations

# Show sublattice properties
from clusterx.model import ModelBuilder, Model
from clusterx.clusters.clusters_pool import ClustersPool
from clusterx.structures_set import StructuresSet
from clusterx.super_cell import SuperCell

def read_property(i, folder, structure=None, **kwargs):
    import os
    f = open(os.path.join(folder, "energy.dat"), 'r')
    erg = float(f.readlines()[0])
    f.close()
    print(folder, erg)
    return erg

property_name = "energy"
# ref
#refs = StructuresSet(db_fname="refs.json")
#refs.read_property_values(property_name, write_to_file=False, read_property=read_property)
#refs.serialize("refs.json")

# Read the previously generated structures_set from .json database file
sset = StructuresSet(db_fname="sset_mc.json")

def read_property(i, folder, structure=None, **kwargs):
    import os
    f = open(os.path.join(folder, "energy.dat"), 'r')
    erg = float(f.readlines()[0])
    f.close()
    print(folder, erg)
    return erg

property_name = "energy"
# ref
refs = StructuresSet(db_fname="refs.json")
refs.read_property_values(property_name, write_to_file=False, read_property=read_property)
refs.serialize("refs.json")
ref_en = refs.get_property_values(property_name)
     
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


from clusterx.thermodynamics.monte_carlo import MonteCarlo
from ase.io.trajectory import Trajectory
from ase.io import write, read


scell = SuperCell(json_db_filepath="scell.json")
for n in range(1, len(scell)):
    print("Now analysing systems with", n, "substitutions.")
    nsubs = {0:[n]}
    mc = MonteCarlo(cemodel,
                    scell,
                    ensemble = "canonical",
                    nsubs=nsubs,
                    predict_swap = True,
                    filename = 'MC1'+str(nsubs)+'.json')

    cp_list = []
    temp_list =[]
# Simulated annealing
    for i,temp in enumerate(range(1400,0,-100)):
       if i == 0:
           init_structure = scell.gen_random(nsubs=nsubs)
       traj = mc.metropolis(no_of_sampling_steps = 1000, \
                         temperature = temp, \
                         boltzmann_constant = kb, \
                         scale_factor = [1/(1.0*nsites)], \
                         initial_decoration = init_structure.get_atomic_numbers())
       cp = traj.calculate_average_property(prop_name = "C_p", \
                         no_of_equilibration_steps = 200)

       print("Sampling at temperature ", temp, "K gives a specific heat of ", cp)
       temp_list.append(temp)
       cp_list.append(cp)
       init_structure = traj.get_structure(-1)

plot_property(temp_list, cp_list, \
             prop_name = "Specific heat", \
             xaxis_label = "Temperature [K]", yaxis_label = r'C$_p$ [k$_B$/#sites]')

