def mpb_solv(results, solvation_energy, rhomin, rhomax): 

    ## Initialises modules
    import os
    import numpy as np
    from decimal import Decimal
    import pandas as pd

    ## Values to be set be user
    # Parametrisation Variables (adjust these to the values in the original python script used to run the mpb).
    rhomin = rhomin
    rhomax = rhomax
    # Set 'solute_energy' to be the energy of the gas-phase optimised solute (without the solvation model)
    sol_energy =  solvation_energy # optimized energy of water using pbe0 and a tight basis  

    ## Script
    data = {'rhomin': [], 
            'rhomax': [],
            'Free Solvation Energy / eV': [],
            }
    dataframe = pd.DataFrame(data)

    for i in rhomin: 
        for j in rhomax:
            i = Decimal(format(i, '.5f'))  # converts formats i to 5 decimal places
            j = Decimal(format(j, '.4f'))  # formats j to 4 decimal places
            folder = str(i)+'_'+str(j) # creates a name for folder in the form {rhomin}_{rhomax}
            pwd = os.path.dirname(__file__) # gets the directory of this python file
            try:
                os.getcwd() # gets current working directory
                os.chdir(pwd+'\\..\\Data\\'+results) # moves the current directory to where the results are located
                os.chdir(os.getcwd()+'\\'+folder) # moves the directory to the defined folder for this loop
            except:
                pass

            ## Appends all lines containing the phase '| Free Energy in Electrolyte' to a list
            is_converged = False
            with open('aims.out','r') as file:
                for n, line in enumerate(file):
                    if 'Have a nice day.' in line:
                        is_converged = True
                        break
            if is_converged:
                line_list = []
                with open('aims.out','r') as file:
                    for n, line in enumerate(file):
                        if '| Free Energy in Electrolyte' in line:
                            line_list.append(line)
							
                    FE = line_list[-1] # Extracts the final value (converged free energy) in the file
                    FE = FE.split()[8] # Extracts the energy in eV (9th string in line)
                    FSE = float(FE) - float(sol_energy) # Calculates the free solvation energy for the given rhomin, rhomax values
                    dataframe = dataframe.append({'rhomin':i, 'rhomax':j, 'Free Solvation Energy / eV':FSE}, ignore_index=True)
            os.chdir(pwd) #goes back one directory to 'pwd'

    ## Returns dataframe
    return dataframe

def contour(x, y, z):

    ## Initialises modules
    import pandas as pd
    from matplotlib import rcParams
    import matplotlib.pyplot as plt
    import numpy as np

    ## Initalises plot objects
    rcParams['figure.figsize'] = 6, 5 # sets plot size
    fig = plt.figure()
    left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
    ax = fig.add_axes([left, bottom, width, height])


    # Make plot and customize axes
    cp = plt.contour(x, y, z, cmap='viridis')
    ax.clabel(cp, inline=True, fontsize=10) 
    plt.xticks(np.arange(0.0010, 0.0045, step=0.0005))
    plt.yticks(np.arange(0.0001, 0.00045, step=0.00005))
    ax.set_xlabel('rhomax')
    ax.set_ylabel('rhomin')
