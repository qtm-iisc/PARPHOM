program compute_electronic_bands

   implicit none

   call initialize_mpi()
   call read_input()
   call blacs_grid_initialization() 
   call read_lammps_data()
   call read_q_file()
   call allocate_distributed_arrays()
   call read_force_constants()
   call diagonalize_and_write()
   call close_mpi()

end program
