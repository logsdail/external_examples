# Read previously generated ParentLattice and SuperCell
from clusterx.super_cell import SuperCell
from clusterx.parent_lattice import ParentLattice
platt = ParentLattice(json_db_filepath="platt.json")
scell = SuperCell(json_db_filepath="scell.json")

from clusterx.clusters.clusters_pool import ClustersPool
import os
if not os.path.exists("cpool.json"):
    cpool = ClustersPool(platt, npoints=[0,1,2,3,4], radii=[0, -1,-1, -1, -1], super_cell=scell)
    print(len(cpool)," clusters were generated.")
    cpool.serialize(db_name="cpool.json")
