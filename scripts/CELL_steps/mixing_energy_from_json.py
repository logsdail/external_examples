import json
from collections import Counter
from clusterx.structures_set import StructuresSet
from clusterx.clusters.clusters_pool import ClustersPool
from clusterx.structures_set import StructuresSet
from clusterx.super_cell import SuperCell
from clusterx.visualization import plot_property

from clusterx.structures_set import StructuresSet
from clusterx.super_cell import SuperCell
sset = StructuresSet(db_fname="sset_mc_new_added_refs_1.json")

list_of_energies = []
list_of_energies_cell = []

with open('sset_mc_new_added_refs_1.json') as f:
    data = json.load(f)
    #print(data)
    list1=range(1,54)
    list2=(list(map(str,list1)))
    #print(list1)
    for n in list2:
        print("number of structure is :", n)
        #print(list2)
        data_numbers = data[n]["numbers"]
        print("the ordering of atooms is :", data_numbers)
        data_numbers_count = Counter(data_numbers)
        print("number_count :", data_numbers_count)
        Pd = "Pd"
        Zn= "Zn"
       
        Pd_numbers=(data_numbers_count[46])
        Zn_numbers=(data_numbers_count[30])
        #print(Pd_numbers)
        y=Pd_numbers + Zn_numbers
        print("number of atoms is :", y)
        Pd_concentration = (Pd_numbers)/y
        formula = Pd_numbers*str(Pd) + Zn_numbers*str(Zn)
        print("formula is:", formula )
        Zn_concentration=(Zn_numbers)/y   #16 is the total number
        print("Zn_numbers_concentration is:", Zn_concentration)
        energy=data[n]['data']['properties']['energy']
        print("total_energy :", energy)
        refs0 = -2224065.94000551    
        refs1 = -786860.0251030402
        energy_mixing= energy-((Pd_concentration*refs0)+(Zn_concentration*refs1))
        energy_mixing_atom = (energy-((Pd_concentration*refs0)+(Zn_concentration*refs1)))/y
        print("energy_mixing :", energy_mixing)
        print("energy_mixing_atom :", energy_mixing_atom)
        
        list_of_energies += [energy_mixing_atom]
        list_of_energies_cell += [energy_mixing]

#appending the energies to json file
sset.set_property_values(property_name="energy_mixing_atom", property_vals=list_of_energies)
sset.set_property_values(property_name="energy_mixing", property_vals=list_of_energies_cell)


      
