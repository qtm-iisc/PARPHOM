## Generation of Force constants using Finite difference method

    
###  Usage
   
    mpiexec -np <no. of cores> get_force_constant -i <lammps input file>   
        
    [ Optional arguments: -o Output file, 
                          -d Displacement of each atom 
                          -s Turn on symmetrization 
                          -t number of symmetry iterations]  
    
###  Theory  

We displace every atom by a finite distance and compute the derivative of the Forces on the atoms due to the displacement. $$\phi_{\alpha,\beta}^{j,j'} = -\frac{\partial F_{\beta}^{j'}}{\partial r_{\alpha}^{j}}$$

