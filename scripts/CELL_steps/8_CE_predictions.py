from clusterx.structures_set import StructuresSet
from clusterx.clusters.clusters_pool import ClustersPool
from clusterx.super_cell import SuperCell
from clusterx.parent_lattice import ParentLattice
from clusterx.model import ModelBuilder, Model
from random import randint

platt = ParentLattice(json_db_filepath="platt.json")
# Define a new, larger super cell
scell = SuperCell(platt, [3, 3, 7])


##################################################################
# Build a Cluster Expansion model based on structures from previous calculations
sset = StructuresSet(db_fname="sset.json")
cpool = ClustersPool(json_db_filepath="cpool.json")

mb = ModelBuilder(selector_type="linreg",
                  selector_opts={'clusters_sets':'size'},
                  estimator_type="skl_LinearRegression",
                  estimator_opts={"fit_intercept":False})

cemodel = mb.build(sset, cpool, "energy") # Build CE model using the training data set
cpool_opt_mc = mb.get_opt_cpool()
cemodel.report_errors(sset)
cpool_opt_mc.display_info(ecis=cemodel.get_ecis())
cpool_opt_mc.write_clusters_db(db_name="cpool_opt_mc.json")
###################################################################


#cemodel = Model(json_db_filepath="CE_model.json")

sset = StructuresSet(platt)
property_name= "energy"
z_coords = scell.get_positions()[:]

nstruc = 1 # we are generating 60 random structures
for i in range(nstruc):
    concentration = {0:[randint(1,len(scell))]} # Pick a random concentration of "Zn"
    sset.add_structure(scell.gen_random(concentration)) # Generate and add a random structure to the StructuresSet
print("Structure generation done now getting predictions")
# Get predictions
predictions = sset.get_predictions(cemodel)
#refs = StructuresSet(db_fname="refs.json")
#ref_en = refs.get_property_values(property_name)
#from clusterx.visualization import plot_property_vs_concentration
#plot_property_vs_concentration(sset, site_type=0, property_name=property_name,refs=ref_en,scale=0.6)

from ase.visualize import view
print(predictions)
view(sset.get_structures())
sset.serialize("sset_predict.json")
sset.write_input_files()
