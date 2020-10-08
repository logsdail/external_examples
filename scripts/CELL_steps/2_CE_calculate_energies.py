# Looping over folders and calculating energies
def calculation_function(structure_locations, property_name, calculator, socket=False):
    '''TODO:Description needed'''
    import os
    from ase.io import read

    parent_directory = os.path.abspath(".")
    parent_directory = str(parent_directory) # path to current folder saved for later

    if socket is True:
        for paths in structure_locations:
            os.chdir(paths)
            if not os.path.exists(str(property_name)+".dat"): # do not perform unnecessary calculations
                model = read("geometry.json")
                with calculator as calc:
                    model.set_calculator(calc) #  socket FHI-aims calculator goes here
                    e_model = model.get_potential_energy()

                    f = open(str(property_name)+".dat", "w")
                    f.write(str(e_model))
                    f.close()
            else:
                pass

            os.chdir(parent_directory) # go back to the parent_directory to finish the loop

    elif socket is False:
        for paths in structure_locations:
            os.chdir(paths)
            if not os.path.exists(str(property_name)+".dat"): # do not perform unnecessary calculations
                model = read("geometry.json")
                model.set_calculator(calculator) # non-socket FHI-aims calculator goes here
                e_model = model.get_potential_energy()

                f = open(str(property_name)+".dat", "w") # TODO: make sure calculations are not repeated if already present
                f.write(str(e_model))
                f.close()

            os.chdir(parent_directory) # go back to the parent_directory to finish the loop

def aims_socket(k_grid):
    # New method that gives a default calculator
    from carmm.run.aims_calculator import get_aims_and_sockets_calculator
    # On machines where ASE and FHI-aims are run separately (e.g. ASE on login node, FHI-aims on compute nodes)
    # we need to specifically state what the name of the login node is so the two packages can communicate
    import socket
    sockets_calc, fhi_calc = get_aims_and_sockets_calculator(dimensions=2, k_grid=k_grid)
    # remove previous xc argument to ensure libxc warning override is first
    fhi_calc.parameters.pop("xc")
    fhi_calc.set(override_warning_libxc="true",
                 xc_pre=['pbe', '10'],
                 xc='libxc MGGA_X_MBEEF+GGA_C_PBE_SOL',
                 spin='none',
                 relativistic=('atomic_zora','scalar'),
                 compute_forces="true",
                 #output=['mulliken'],
                 use_dipole_correction='True',
                 final_forces_cleaned='true',
                 #sc_accuracy_etot=1e-6,
                 #sc_accuracy_eev=1e-3,
                 #sc_accuracy_rho=1e-6,
                 #sc_accuracy_forces=1e-4
                 )

    return sockets_calc

#############################################
from ase.calculators.emt import EMT
with open("folder_paths.txt", "r") as f:
    structure_locations = f.read().split(' ')

calculation_function(structure_locations, property_name="energy", calculator=EMT())

#from carmm.run.aims_path import set_aims_command
#set_aims_command(hpc='thomas', basis_set='tight')
#calculation_function(structure_locations, property_name="energy", calculator=aims_socket((8,8,1)), socket=True)
