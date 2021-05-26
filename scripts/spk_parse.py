##################################################################################### READ ME ####################################################################################################

### This script parses the data (including energies and forces) for the trajectory files of a given directory and converts them into a form which can be read by Schnetpack.
## Here, the script uses the principles highlighted on the online documentation (https://schnetpack.readthedocs.io/en/stable/tutorials/tutorial_01_preparing_data.html#Preparing-your-own-data,
# as of May 2021). However, in this script, the energies and forces for each step of the trajectory file are extracted and appended to a list of dictionaries (property_list) before being
# passed into a database.
## This script was originally created in a PyCharm virtual environment with the following dependencies:
# Dependencies: Python (3.8), Schnetpack (version 0.3), PyYaml (5.4.1), numpy (1.20.3), h5py (3.2.1), tqdm (4.60.0), ASE (3.21.1), PyTorch (1.8.1), tensorboardX (2.2).

## NOTE THAT THIS SCRIPT HAS NOT BEEN TESTED FOR TRAINING MODELS HOWEVER THE DATABASE FORMAT IS IN-LINE WITH THE ONLINE DOCUMENTATION.

##################################################################################################################################################################################################

### INITIALISES RELEVANT MODULES

import schnetpack
from ase.io.trajectory import Trajectory
from ase import Atoms
from ase.io import read, write
from schnetpack.data import AtomsData
import numpy as np
import os


### CREATES A "EMPTY" PYTORCH DATASET TO PARSE THE DATA TO

## TO BE DEFINED BY USER - defines file location for trajectory files and the database name
dirpath = 'INSERT_DIRECTORY_PATH_HERE'
dbname = "test.db"

## If a database of the same name already exists, it is removed as otherwise it will cause this script to fail
if os.path.isfile(os.path.join(dirpath, dbname)):
    os.remove(os.path.join(dirpath, dbname))

## Creates the schnetpack database (spk_db) in the given directory and defines the properties we are interested in
spk_db = AtomsData(os.path.join(dirpath, dbname), available_properties=['energy', 'forces'])  # note that {name}.db must not previously exist in said directory


### APPENDING PROPERTIES TO THE DATABASE

## Parses the energy and forces for every image of every trajectory file in the given directory to the previously defined database
for root, dirs, files in os.walk(dirpath):
    for name in files:
        ## Defines the trajectory and requests all images
        trajectory = read(os.path.join(dirpath, name + "@:"))

        ## Extracts the energies and forces for said trajectory file and puts them in a list of dictionaries
        property_list = [{"energy": np.array([atoms.get_potential_energy()], dtype=np.float32)
                         ,"forces":np.array([atoms.get_forces()], dtype=np.float32)
                        } for atoms in trajectory]

        ## Adds the system geometry (trajectory) and the properties (property_list) of every trajectory image to the previously defined database
        spk_db.add_systems(trajectory, property_list)



### VISUALIZING PROPERTIES OF THE DATABASE

## Shows what properties the database has
#for p in spk_db.available_properties:
#    print('-', p)

## Gives an output showing the properties of the database
print('Number of reference calculations:', len(spk_db))
print('Available properties:')
for p in spk_db.available_properties:
    print('-', p)
print()

## Displays the available properties and their sizes as a torch array
example = spk_db[0]
print('Properties of the 0th database entry:')
for k, v in example.items():
    print('-', k, ':', v.shape)
print()

## Displays the first entry in the database
print('Below the 0th database entry is shown:\n', spk_db.get_properties(0))




