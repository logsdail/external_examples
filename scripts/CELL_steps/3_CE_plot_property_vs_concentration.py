import os
from clusterx.structures_set import StructuresSet
from clusterx.visualization import plot_property_vs_concentration


# need a way of retrieving custom values, current template is not working as intended
def read_property(i, folder, structure=None, **kwargs):
    import os
    f = open(os.path.join(folder, "energy.dat"), 'r')
    erg = float(f.readlines()[0])
    f.close()
    print(folder, erg)
    return erg


property_name = "energy"
sset = StructuresSet(db_fname="sset.json")
sset.read_property_values(property_name, write_to_file=False, read_property=read_property)
sset.serialize("sset.json")


refs = StructuresSet(db_fname="refs.json")
refs.read_property_values(property_name, write_to_file=False, read_property=read_property)
refs.serialize("refs.json")
ref_en = refs.get_property_values(property_name)

#plot_property_vs_concentration(sset, site_type=0, property_name=property_name,refs=ref_en,scale=0.6)

# See the most stable structure
from ase.visualize import view
import numpy as np
array = sset.get_property_values(property_name)
lowest_energy_structure_index = np.where(array == np.amin(array))[0][0]
# GET STRUCTURES/ GET STRUCTURE IS NOT SHOWING THE CORRECT STRUCTURE
lowest_energy_structure = sset.get_structures()[lowest_energy_structure_index]
view(lowest_energy_structure)
