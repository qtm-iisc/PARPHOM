## Generation of Force constants using Finite difference method

    
###  Usage
   
    mpiexec -np <no. of cores> get_force_constant -i <lammps input file>   
        
    [ Optional arguments: -o Output file, -d Displacement of each atom]  
    
###  Theory  

We displace every atom by a finite distance and compute the derivative of the Forces on the atoms due to the displacement.  
```
        $\phi_{a,b}^{i,j} = -\frac{dF_{b,j}}{dr_{a,i}}$
```
