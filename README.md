# PARPHOM

> **A package to calculate phonon frequencies in moiré systems**

---

## 🚦 Prerequisites

- [**MPI**](https://www.mpich.org/static/downloads/4.0.2/mpich-4.0.2-installguide.pdf)
- [**Parallel HDF5**](https://support.hdfgroup.org/documentation/hdf5/latest/_cookbook.html)
- [**ScaLAPACK**](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html)

---

## ⚡ Usage Instructions

### 🛠️ Force Constant Generation

1. **Prepare a LAMMPS input file.**
2. **Install the serial LAMMPS Python wrapper.**  
   _(MPI Python wrapper is not supported.)_  
   See the [LAMMPS documentation](https://docs.lammps.org/Python_install.html) for installation instructions.
3. **Generate the force constants:**

   ```bash
   mpiexec -np <no. of cores> Force_Constant_Generation/get_force_constant.py -i <LAMMPS input file> -o <output file> -d <displacement>
   ```

---

### 🎵 Phonon Frequency Calculation

1. **Build the code:**
   ```bash
   make -f Makefile.<dist>   # supported: Intel MKL, AOCL
   ```
   _Tip: Use the `-D__DEBUG` flag for debugging (generates large log files)._  
   
   💡 **Note:**
   - In the `Makefile.intel` (if building with INTEL compilers), you can change the Fortran compiler to the appropriate compiler wrapped with HDF5.
   - If you are using native Fortran compilers, make sure to link the HDF5 libraries separately.
   - Change the `Makefile.aocl` appropriately for AOCC compilers.
   - You can also build the code with any other compiler of your choice. To do that, copy one of the provided Makefiles and replace the Fortran compiler (`FC`) with the compiler of your choice. Put the paths to the BLAS, LAPACK, and ScaLAPACK libs in `BLAS_LIBS`, `LAPACK_LIBS`, and `SCALAPACK_LIBS` respectively. Also update your Fortran compiler flags (`FCFLAGS`) and the linked library flags (`FLFLAGS`) appropriately. Then build with the Makefile.

2. **Prepare a q-point file** using the Utilities folder examples.  
   _See the Utilities/README for details._  
   You can use `Utilities/Examples/gen_band_q_pt.py` and `Utilities/Examples/gen_dos_q_pt.py` as examples for generating q-points along high-symmetry lines and across the BZ, respectively.

3. **Create an input file** for the phonon code. Use the following template:

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

4. **Run the code:**
   ```bash
   mpiexec -np <no. of cores> parphom.x <input file>
   ```

---

## 📝 Notes
- Ensure all required input and potential files are present before running calculations.
- The workflow is compatible with LAMMPS and the phonon calculation tools used in this project.

---

## 🐞 Known Issues
- q-pools are not supported yet. Will be supported in v2.0.
