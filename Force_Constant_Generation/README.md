## Generation of Force constants using Finite difference method

    
###  Usage
   
    mpiexec -np <no. of cores> get_force_constant -i <lammps input file>   
        
    [ Optional arguments: -o Output file, 
                          -d Displacement of each atom 
                          -s Turn on symmetrization 
                          -t number of symmetry iterations]  
    
###  Theory  

We displace every atom by a finite distance and compute the derivative of the Forces on the atoms due to the displacement. 
$$\phi_{\alpha,\beta}^{j,j'} = -\frac{\partial F_{\beta}^{j'}}{\partial r_{\alpha}^{j}}$$


## Known issue with Drip Potential

For potentials that work between different molecular types, one has to define the criteria for defining the molecular types in a replicated layout. This is an artefact of how the Drip and other potentials that define interactions between different molecular types work as mentioned [here](https://github.com/lammps/lammps/issues/3047).  
For using the script with these potentials and replicated structures, uncomment lines 209-214.  
Note that if the user wants to work with multilayered graphene, cutoffs for defining each layer has to be defined seperately and appropriate groups have to be created in LAMMPS as per requirements.  
