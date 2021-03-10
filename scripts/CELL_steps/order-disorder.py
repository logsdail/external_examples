from ase.io import read, write
from ase.build import molecule
from ase import Atoms
from ase.constraints import FixAtoms
from ase.utils import hsv
from collections import Counter

atoms=read('concentration-05875-new.traj')  #structure to check the disorder
chemical= atoms.get_chemical_symbols()
numbers=Counter(chemical)
print("number_count_atoms :", numbers)

atoms_ordered=read('concentration-0585-222.traj')  # the reference structure
chemical_ordered=atoms_ordered.get_chemical_symbols()
numbers=Counter(chemical)

#definition of the disorder factor with Pdzn is the Pd replacing Zn in ordered structure and ZnPd is the Zn replacing Pd in prderd structure
#disorder= (Pdzn+ZnPd)/(PdPd+ZnZn) as described in https://doi.org/10.1080/01614940.2013.796192


n=0 #index of chemical and chemical_ordered
m=0 #number of disoredred
l=0 #number of ordered
#print(len(chemical))

while n < len(chemical):
    if chemical[n] == chemical_ordered[n]:
        l= l+1
    else:
       m= m+1
    n = n+1

print('number of ordered site is :',l)
print('number of disordered site is :', m)
print('disorder factor is :', m/len(chemical))



