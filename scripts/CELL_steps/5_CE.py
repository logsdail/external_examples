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

from clusterx.visualization import plot_optimization_vs_number_of_clusters
from clusterx.visualization import plot_predictions_vs_target
plot_optimization_vs_number_of_clusters(mb.get_selector(),scale=0.5)
plot_predictions_vs_target(sset,cemodel,"total_energy_emt",scale=0.5)
plot_property_vs_concentration(sset, site_type=0, property_name="total_energy_emt",cemodel=cemodel,refs=ref_en,scale=0.5)


