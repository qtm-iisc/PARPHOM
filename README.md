# Moire Phonons

Code to calculate the phonon frequencies in Moire systems. 

## Prerequisites

The prerequisites to run this code are:

>
>   Parallel HDF5 (checked for v1.12.0.0 and v1.12.1)   
>   Scalapack (checked with Intel MKL and AMD Scalapack)  
>

## Instructions for usage

### For the Force constant generation

### For the Phonon frequency calculation

Build the code using:
```
make -f Makefile.<dist> (supported dists = Intel MKL, AOCL)
```

While compilation, the user can turn on the `-D__DEBUG` flag to get debugging information. Please note that this will generate a ver large log file. Use it with appropriate caution.

A sample input file would look like this:
```
# ---------------------------------------------------
# Input file for Phonon Calculations in Moire Systems
# ---------------------------------------------------


lammps file location         : 
lammps file name             :
natom                        :
atom types                   :
atom style                   :          # LAMMPS atom style used in the data file
                                        # Supported styles are: atomic and molecular 


q file location              : 
q file name                  : 
nqpt                         :          # Number of q points


force constant file location : 
force constant file name     : 
force constant dataset       :          # Name of the dataset containing the force constants
                                        # Dimension of this dataset should be (N,N,3,3)

num qpools                   :          # Number of q-pools --- Not implemented in v1.0


compute eigvecs              : < yes/no >
range                        : < A/I/V >   # A = compute all evals
                                           # I = compute eigenvalues from min index to max index
                                           # V = compute eigenvalues from a min-max range

min eigval                   :          # Minimum value of the eigenvalue to be found
                                        # ignored if compute eigvecs = A/I
max eigval                   :          # Maximum value of the eigenvalue to be found
                                        # ignored if compute eigvecs = A/I
min index                    :          # lowest index of the eigenvalue to be found
                                        # ignored if range = A/V
max index                    :          # highest index of the eigenvalue to be found
                                        # ignored if range = A/V

mb                           :          # (mb,nb) are the block sizes for the block cyclic
nb                           :          # distribution of the matrix

abstol                       :          # Tolerance for eigenvalue convergence
orfac                        :          # Tolerance for eigenvector orthogonaliztion

output file name             : 
output_file_location         : 
```


## Known Issues

Q-pools are not supported yet. Will be supported in v2.0

