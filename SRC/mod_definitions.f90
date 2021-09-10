MODULE global_variables
      
      INTEGER :: allocate_status

      ! MPI Variables    

      TYPE mpi_group
          INTEGER :: comm
          INTEGER :: color
          INTEGER :: key
          INTEGER :: size_
          INTEGER :: rank
      END TYPE mpi_group
      
      TYPE(mpi_group) :: mpi_local, mpi_global

      INTEGER :: no_of_groups 
      INTEGER :: mpierr

      ! BLACS variables

      TYPE blacs_info
          INTEGER, ALLOCATABLE :: context(:)
          INTEGER :: nprow
          INTEGER :: npcol
          INTEGER :: myprow
          INTEGER :: mypcol
      END TYPE blacs_info

      TYPE(blacs_info) :: blacs_grid

      INTEGER, PARAMETER :: BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,     &
                            CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,    &
                            RSRC_ = 7, CSRC_ = 8, LLD_ = 9

      ! Input File read
      ! ---------------
      
      
      TYPE lammps
          CHARACTER(500) :: location, name, suffix, atom_style
          INTEGER        :: at_types
      END TYPE lammps
      
      TYPE(lammps) :: lammps_file
      
      
      
      TYPE q_points
          CHARACTER(500) :: location, name
          INTEGER        :: nqpt
          DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: points
      END TYPE
      
      TYPE(q_points) :: q_file
      
      
      
      TYPE force_constant
          CHARACTER(500) :: location, file_name, dataset
          INTEGER :: size_, lld, ncol
          DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: a
          INTEGER, DIMENSION(DLEN_) :: desca
      END TYPE force_constant
      
      TYPE(force_constant) :: fc
      


      TYPE scalapack_variables
          LOGICAL      :: comp_evec
          INTEGER      :: mb
          INTEGER      :: nb
          CHARACTER(1) :: range_
          INTEGER      :: il_
          INTEGER      :: iu_
          INTEGER      :: vl_
          INTEGER      :: vu_
          INTEGER      :: num_eval_comp
          INTEGER      :: num_evec_comp
      END TYPE scalapack_variables
      
      TYPE(scalapack_variables) :: pzheevx_vars 



      TYPE system
          INTEGER :: natom
          DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: mass
          INTEGER,          ALLOCATABLE, DIMENSION(:)   :: at_types_i
          DOUBLE PRECISION, DIMENSION(3,3)              :: lat, rec_lat
          DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: real_pos, crys
      END TYPE system

      TYPE(system) :: moire


    
      TYPE complex_arrays
          DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:) :: a
          INTEGER, DIMENSION(DLEN_) :: desca
          INTEGER :: size_, lld, ncol
      END TYPE complex_arrays
      
      TYPE(complex_arrays) :: dynmat, evec

      ! Q file read

      INTEGER :: st_q, en_q, q_index 
      INTEGER :: info, ierr(1)
      LOGICAL :: isok

      ! Eigenvectors and eigenvalues

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: W
      
END MODULE
