# Force Constant Generation Example: 21.8° Twisted Bilayer Graphene

Generate the force constants for 21.8° twisted bilayer graphene using the provided script and input files.

## How to Run

```
mpiexec -np <number_of_cores> python ../../get_force_constant.py -i in.phonon -o FORCE_CONSTANTS -d 0.0001
```
- Replace `<number_of_cores>` with the number of CPU cores you want to use.
- Adjust the displacement value (`-d 0.0001`) as needed for your calculation accuracy.

## Files in This Folder
- **in.phonon**: Input file for the force constant calculation
- **FORCE_CONSTANTS**: Output file containing the generated force constants
- **C.drip, CH.rebo, BNC.tersoff**: Potential files 

## Tips
- Make sure all required input and potential files are present before running the script.
- For different systems or displacement values, simply modify the input arguments as needed.

---
For more details, refer to the main project documentation or contact the authors.
