# PARPHOM

Package to calculate the phonon frequencies in moiré systems. 

## Prerequisites

The prerequisites to run this code are:

>
>   MPI  
>   Parallel HDF5 (checked for v1.12.0.0 and v1.12.1)   
>   Scalapack (checked with Intel MKL and AMD Scalapack)  
>

## Instructions for usage

### For the Force constant generation

#### Step 1: 
Prepare a LAMMPS input file
#### Step 2:
Make sure you have the serial LAMMPS python wrapper (the MPI python wrapper is not supported for now). For instructions on how to install the LAMMPS serial python wrapper please refer to [the LAMMPS documentation](https://docs.lammps.org/Python_install.html)
#### Step 3:
Generate the force constant using 
```
mpiexec -np <no. of cores> Force_Constant_Generation/get_force_constant.py -i <LAMMPS input File> -o <Name of the output File> -d <displacement of each atom>
```

### For the Phonon frequency calculation

#### Step 1: 
Build the code using:
```
make -f Makefile.<dist> (supported dists = Intel MKL, AOCL)
```
While compilation, the user can turn on the `-D__DEBUG` flag to get debugging information. Please note that this will generate a very large log file. Use it with appropriate caution.


#### Step 2:  
Prepare a q-point file using the examples in the Utilities folder. For more information refer to README in the Utilities folder

#### Step 3:  
Prepare an input file for the phonon code. A sample input file would look like this:  
```
# ---------------------------------------------------
# Input file for Phonon Calculations in Moire Systems
# ---------------------------------------------------

# LAMMPS Structure
# ----------------
lammps file location         : 
lammps file name             :
natom                        :
atom types                   :
atom style                   :          # LAMMPS atom style used in the data file
                                        # Supported styles are: atomic and molecular 

# Q Points  
# --------
q file location              : 
q file name                  : 
nqpt                         :          # Number of q points

# Force constants  
# ---------------
force constant file location : 
force constant file name     : 
force constant dataset       :          # Name of the dataset containing the force constants
                                        # Dimension of this dataset should be (N,N,3,3)


# Calculation Parameters
# ----------------------
mb                           :                          # (mb,nb) are the block sizes for the block cyclic
nb                           :                          # distribution of the matrix
abstol                       :                          # Tolerance for eigenvalue convergence
orfac                        :                          # Tolerance for eigenvector orthogonaliztion

# Quantities to be calculated
# ---------------------------
compute eigvecs              : < yes/no, default= no >
range                        : < A/I/V, default = A >   # A = compute all evals
                                                        # I = compute eigenvalues from min index to max index
                                                        # V = compute eigenvalues from a min-max range
min eigval                   :                          # Minimum value of the eigenvalue to be found, ignored if compute eigvecs = A/I
max eigval                   :                          # Maximum value of the eigenvalue to be found, ignored if compute eigvecs = A/I
min index                    :                          # lowest index of the eigenvalue to be found, ignored if range = A/V
max index                    :                          # highest index of the eigenvalue to be found,ignored if range = A/V
group velocity               : < yes/no, default = no >
velocity method              : < A/CD, default = A >    # A = Analytic method of computing dynamical derivatives
                                                        # CD = Dynamical matrix derivatives generated from central differences


# Output File
# -----------
output file name             : 
output_file_location         : 
```

#### Step 4:  
Run the code as 
```
mpiexec -np <no. of cores> parphom.x <input file>
```


**For post processing refer to the README and examples in the Utilities folder**


## Known Issues

q-pools are not supported yet. Will be supported in v2.0
