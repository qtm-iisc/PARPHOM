    PROGRAM PHONON_FREQUENCIES_FINAL

    ! -----------------------------------------------------------------------
    ! Program to read the force constants, distribute it among the processes,
    ! create the dynamical matrix and diagonalize it
    ! -----------------------------------------------------------------------
    USE HDF5 
    IMPLICIT NONE

    ! Defining the variables
    ! ----------------------

        ! **  BLACS Definitions **
 
    INTEGER :: iam, nprocs, nprow, npcol, icontext, myprow, mypcol
    LOGICAL :: isroot

        ! ** Reading Data **

    INTEGER        :: natom                 
    INTEGER        :: at_types              
    CHARACTER(1)   :: atoms_style           
    INTEGER        :: nqpt                  
    INTEGER        :: mb, nb
    CHARACTER(500) :: location, fil, suffix, q_file 
    CHARACTER(LEN=500)  :: FC_file       
    CHARACTER(LEN=500)  :: FC_dataset_name 
        
        ! ** Reading LAMMPS data and q-points **

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)    :: mass
    DOUBLE PRECISION,              DIMENSION(3,3)  :: lat, rec_lat
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)  :: real_pos, crys
    INTEGER,          ALLOCATABLE, DIMENSION(:)    :: at_types_i
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)  :: q_array
    DOUBLE PRECISION,              DIMENSION(3)    :: q_pt


        ! ** Distributing Matrix **

    INTEGER  :: csrc, rsrc, lld, Locp, Locq, size_need
    INTEGER  :: numroc
    EXTERNAL :: numroc
    
    INTEGER  :: AllocateStatus

    DOUBLE PRECISION, ALLOCATABLE :: force_cnst_dist(:)
    INTEGER, PARAMETER :: BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,     &
                          CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,    &
                          RSRC_ = 7, CSRC_ = 8, LLD_ = 9

    INTEGER  :: descforce(DLEN_)
    INTEGER  :: info, ierr(1)
    LOGICAL  :: isok
    
    INTEGER  :: ia_first, ja_first, iastart, jastart, iaend, jaend
    INTEGER  :: ia, ja
    INTEGER  :: lroffset, lcoffset
    
    
        ! ** Iterate over q points **
    
    INTEGER  :: q_index
    DOUBLE COMPLEX, ALLOCATABLE :: dyn_matrix(:)
    DOUBLE COMPLEX :: d_ij
    INTEGER  :: descdyn(DLEN_)
    INTEGER  :: lrindx, lcindx, ipos
    DOUBLE PRECISION, DIMENSION(3) :: rc_i, rc_j
    
        ! ** Diagonalization **
    
    DOUBLE COMPLEX, ALLOCATABLE  :: Z(:), work(:)
    INTEGER  :: descEigvec(DLEN_)
    DOUBLE PRECISION, ALLOCATABLE :: W(:)
    DOUBLE PRECISION, ALLOCATABLE :: rwork(:), gap(:)
    INTEGER, ALLOCATABLE  :: iwork(:), ifail(:), iclustr(:)
    INTEGER  :: lwork, liwork, lrwork, lgap, lifail, liclustr, numval, numvec
    INTEGER, PARAMETER  :: abstol = -1
    INTEGER, PARAMETER  :: orfac = -1
    
        ! ** HDF5 Reader variables **

    INTEGER(HID_T)  :: file_id
    INTEGER(HID_T)  :: dset_id
    INTEGER :: error
    INTEGER(HID_T)  :: space, memspace
    INTEGER(SIZE_T), PARAMETER :: num_ele = 2
    INTEGER(HSIZE_T), DIMENSION(4, num_ele) :: coord
    INTEGER(HSIZE_T), DIMENSION(1) :: count_
    INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
    DOUBLE PRECISION, DIMENSION(num_ele) :: phi_ij
        
       ! ** Print Array **

    INTEGER :: irprnt, icprnt
    INTEGER, PARAMETER :: nout = 16
!    DOUBLE PRECISION, DIMENSION(20000) :: pwork
!    DOUBLE PRECISION, DIMENSION(20000) :: pzwork


       ! ** Read input file **

    CHARACTER(500) :: input_file



    ! ------------ 
    ! Main Program
    ! ------------

    
    
    ! ** Initiate BLACS Environment **
    ! -------------------------------- 
    

    CALL blacs_pinfo(iam, nprocs)
    DO nprow = int(sqrt(real(nprocs)))+1, 1, -1
      npcol = nprocs/nprow
      IF (nprow*npcol.eq.nprocs) GOTO 11
    END DO
11  CONTINUE

    CALL blacs_get(-1,0,icontext)
    CALL blacs_gridinit(icontext, 'Col-major', nprow, npcol)
    CALL blacs_gridinfo(icontext, nprow, npcol, myprow, mypcol)

    isroot = (myprow.eq.0).and.(mypcol.eq.0)
  


    ! ** Read the crystal coordinates, the lattice parameters, list of q points **
    ! ----------------------------------------------------------------------------
    !
    ! Variables -
    !
    ! natom   :: Number of atoms in the system
    ! crys    :: array of atomic positions in crystal coordinates
    ! lat     :: lattice vectors
    ! q_array    :: array of q points 
    ! mb, nb  :: Block size of distributed matrix
    !

    CALL GET_COMMAND_ARGUMENT(1,input_file)
    IF (LEN_TRIM(input_file) == 0) STOP "Enter proper input file"
 
    CALL read_input(input_file, location, fil, suffix, natom, at_types, atoms_style, &
                    mb, nb, nqpt, q_file, FC_file, FC_dataset_name)
    
    ALLOCATE(mass(at_types))
    ALLOCATE(real_pos(natom,3))
    ALLOCATE(crys(natom,3))
    ALLOCATE(at_types_i(natom))
    
    CALL read_param(natom, at_types, atoms_style, location, fil, suffix, mass, lat, rec_lat, &
                      real_pos, crys, at_types_i)


    ALLOCATE(q_array(nqpt,3))    

    CALL read_q_file(q_file, location, nqpt, q_array)



    ! ** Read and distribute force constant matrix
    
    rsrc = 0
    csrc = 0

    Locq = numroc(3*natom,mb,mypcol,csrc,npcol)
    Locq = max(1,Locq)
    Locp = numroc(3*natom,nb,myprow,rsrc,nprow)
    lld = max(1,Locp)

    size_need = 1+lld*Locq

    IF (size_need.le.0) STOP "Memory overflow. Increase number of processors"
    ALLOCATE(force_cnst_dist(size_need), STAT = AllocateStatus)
    IF (AllocateStatus.ne.0) STOP "Insufficient Memory for Force Constants"

    CALL descinit(descforce, 3*natom, 3*natom, mb, nb, rsrc, csrc, icontext, lld, info)

    ierr(1) = 0
    CALL igsum2d(icontext, 'All', ' ', 1, 1, ierr, 1, -1, -1)
    isok = (info.eq.0)
    IF (.not.isok) THEN
        IF (isroot) THEN
            WRITE(*,*) 'descinit for force constant distributed matrix returned info = ', info
        ENDIF
        GOTO 99
    ENDIF

    
    CALL blacs_barrier(icontext, 'All')

    
    !   
    ! ** Read the force constant matrix in block cyclic distributed fashion ** !
    !   
       
    ! -------- Compute the first array index on local processor

    IF (myprow.ge.descforce(RSRC_)) THEN
        ia_first = (myprow - descforce(RSRC_))*descforce(MB_) + 1
    ELSE
        ia_first = (myprow + (nprow-descforce(RSRC_)))*descforce(MB_) + 1
    ENDIF


    IF (mypcol.ge.descforce(CSRC_)) THEN
        ja_first = (mypcol - descforce(CSRC_))*descforce(NB_) + 1
    ELSE
        ja_first = (mypcol + (npcol-descforce(CSRC_)))*descforce(NB_) + 1
    ENDIF

    
    ! ---------- Load the HDF5 file and the dataset in each processor ----- !

    ! open hdf5 interface
    CALL h5open_f(error) 
    
    ! open hdf5 file                         
    CALL h5fopen_f(TRIM(ADJUSTL(location))//TRIM(ADJUSTL(FC_file)), H5F_ACC_RDONLY_F, file_id, error)    
    
    ! open dataset
    CALL h5dopen_f(file_id, TRIM(ADJUSTL(FC_dataset_name)), dset_id, error)

    ! retrieve dataspace
    CALL h5dget_space_f(dset_id, space, error)   
    
    count_(1) = num_ele
    data_dims(1) = num_ele
   
    
    !           skip by npcol*nb and nprow*mb 

    DO jastart = ja_first, descforce(N_), npcol*descforce(NB_)
        DO iastart = ia_first, descforce(M_), nprow*descforce(MB_)
            iaend = min(descforce(M_), iastart+descforce(MB_)-1)
            jaend = min(descforce(N_), jastart+descforce(NB_)-1)

            ! ---------------------------------------------
            ! block (iastart:iaend, jastart:jaend) is 
            ! within the same block on the local processor
            !  
            ! Need to compute local array index for 1st entry
            ! in the local block
            ! ---------------------------------------------
        
            ia = iastart
            ja = jastart
            
            CALL infog2l(ia, ja, descforce, nprow, npcol, myprow, mypcol, lroffset, lcoffset, rsrc, csrc)
            
            DO ja=jastart, jaend
              DO ia=iastart, iaend
                
                lrindx = lroffset + (ia-iastart)
                lcindx = lcoffset + (ja-jastart)
                ipos = lrindx + (lcindx-1)*descforce(LLD_)
                
                CALL twod_to_fourd(ia,ja,coord)
                
                CALL h5sselect_elements_f(space, H5S_SELECT_SET_F, 4, num_ele, coord, error)
                CALL h5screate_simple_f(1,count_,memspace,error)
                CALL h5dread_f(dset_id, H5T_IEEE_F64LE, phi_ij, data_dims, error, memspace, &
                               file_space_id=space)

                force_cnst_dist(ipos) = (phi_ij(1)+phi_ij(2))*0.5
                CALL h5sclose_f(memspace, error)

              END DO
            END DO 
        END DO
    END DO

    CALL h5sclose_f(space, error)
    CALL h5dclose_f(dset_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

    ! ------ Force constant files read complete ----------------
    
    IF (isroot) THEN
      WRITE(*,*) "Force Constant Read"
    END IF

    ! ------------------------
    ! Create dynamical matrix
    ! ------------------------


    DO q_index = 1, nqpt
        
        q_pt(1) = q_array(q_index,1)
        q_pt(2) = q_array(q_index,2)
        q_pt(3) = q_array(q_index,3)
        

        ! ** Allocate dynamical matrix
        
        rsrc = 0
        csrc = 0

        Locq = numroc(3*natom,mb,mypcol,csrc,npcol)
        Locq = max(1,Locq)
        Locp = numroc(3*natom,nb,myprow,rsrc,nprow)
        lld = max(1,Locp)

        size_need = 1+lld*Locq
        ALLOCATE(dyn_matrix(size_need), STAT = AllocateStatus)
        IF (AllocateStatus.ne.0) STOP "Insufficient Memory for Dynamical Matrix"
        CALL descinit(descdyn, 3*natom, 3*natom, mb, nb, rsrc, csrc, icontext, lld, info)
        ierr(1) = 0
        CALL igsum2d(icontext, 'All', ' ', 1, 1, ierr, 1, -1, -1)
        isok = (info.eq.0)
        IF (.not.isok) THEN
            IF (isroot) THEN
                WRITE(*,*) 'descinit for dynamical matrix returned info = ', info
            ENDIF
            GOTO 99
        ENDIF

        CALL blacs_barrier(icontext, 'All')

        ! ** locate the correct global indices and generate the element from the force constant matrix

        IF (myprow.ge.descdyn(RSRC_)) THEN
            ia_first = (myprow - descdyn(RSRC_))*descdyn(MB_) + 1
        ELSE
            ia_first = (myprow + (nprow-descdyn(RSRC_)))*descdyn(MB_) + 1
        ENDIF


        IF (mypcol.ge.descdyn(CSRC_)) THEN
            ja_first = (mypcol - descdyn(CSRC_))*descdyn(NB_) + 1
        ELSE
            ja_first = (mypcol + (npcol-descdyn(CSRC_)))*descdyn(NB_) + 1
        ENDIF

        DO jastart = ja_first, descdyn(N_), npcol*descdyn(NB_)
            DO iastart = ia_first, descdyn(M_), nprow*descdyn(MB_)
                iaend = min(descdyn(M_), iastart+descdyn(MB_)-1)
                jaend = min(descdyn(N_), jastart+descdyn(NB_)-1)
                ia = iastart
                ja = jastart
                CALL infog2l(ia,ja,descdyn,nprow,npcol,myprow,mypcol,lroffset,lcoffset,rsrc,csrc)
                DO ja = jastart, jaend

                    rc_j(1) = crys((ja-1)/3 + 1, 1)
                    rc_j(2) = crys((ja-1)/3 + 1, 2)
                    rc_j(3) = crys((ja-1)/3 + 1, 3)
                        
                    DO ia = iastart, iaend

                        lrindx = lroffset + (ia-iastart)
                        lcindx = lcoffset + (ja-jastart)
                        ipos = lrindx+(lcindx-1)*descdyn(LLD_)
                        
                        phi_ij(1) = force_cnst_dist(ipos)

                        rc_i(1) = crys((ia-1)/3 + 1, 1)
                        rc_i(2) = crys((ia-1)/3 + 1, 2)
                        rc_i(3) = crys((ia-1)/3 + 1, 3)
                      
                        
                        IF (ia.eq.ja) THEN 
                          dyn_matrix(ipos) = phi_ij(1)/mass(at_types_i( (ia-1)/3+1 ))
                        ELSE
                          CALL create_dynamical_ij(phi_ij, q_pt, rc_i, rc_j, lat, d_ij)
                          dyn_matrix(ipos) = d_ij/ SQRT( mass( at_types_i( (ia-1)/3+1 ))* &
                                                         mass( at_types_i( (ja-1)/3+1 )))
                        END IF

                    END DO
                END DO
            END DO
        END DO

        CALL blacs_barrier(icontext,'All')

        ! *** Diagonalize the Dynamical Matrix ***
        
        rsrc = 0
        csrc = 0

        Locq = numroc(3*natom,mb,mypcol,csrc,npcol)
        Locq = max(1,Locq)
        Locp = numroc(3*natom,nb,myprow,rsrc,nprow)
        lld = max(1,Locp)

        size_need = 1+lld*Locq

        ALLOCATE(Z(size_need), STAT=AllocateStatus)
        
        CALL descinit(descEigvec, 3*natom, 3*natom, mb, nb, rsrc, csrc, icontext, lld, info)
        
        ierr(1) = info
        CALL igsum2d(icontext, 'All', ' ', 1,1, ierr, 1,-1,-1)
        isok = (info.eq.0)
        IF (.not.isok) THEN
          IF (isroot) THEN
             WRITE(*,*) 'Descinit for Eigenvectors returns info = ', info
          END IF
          GOTO 99
        END IF

        ALLOCATE(W(3*natom))

        ! Allocate optimum arrays for PZHEEVX

        ia = 1
        ja = 1
        lgap = nprow*npcol
        lifail = 3*natom
        liclustr = 2*nprow*npcol
        ALLOCATE(gap(lgap))
        ALLOCATE(ifail(lifail))
        ALLOCATE(iclustr(liclustr))
        lwork = -1
        lrwork = -1
        liwork = -1
        
        ALLOCATE(work(1))
        ALLOCATE(rwork(1))
        ALLOCATE(iwork(1))

        CALL PZHEEVX('N', 'A', 'U', 3*natom, dyn_matrix, ia, ja, descdyn, 0.0D0, 0.0D0, 1, 3*natom,  &
                     abstol, numval, numvec, W, orfac, Z, 1, 1, descEigvec, work, lwork, rwork, &
                     lrwork, iwork, liwork, ifail, iclustr, gap, info)

        lwork = int(work(1)+1)
        lrwork = int(rwork(1)+1)
        liwork = int(iwork(1)+1)

        DEALLOCATE(work)
        DEALLOCATE(iwork)
        DEALLOCATE(rwork)

        ALLOCATE(work(lwork))
        ALLOCATE(rwork(lrwork))
        ALLOCATE(iwork(liwork))

        CALL PZHEEVX('N', 'A', 'U', 3*natom, dyn_matrix, ia, ja, descdyn, 0.0D0, 0.0D0, 1, 3*natom,  &
                     abstol, numval, numvec, W, orfac, Z, 1, 1, descEigvec, work, lwork, rwork, &
                     lrwork, iwork, liwork, ifail, iclustr, gap, info)

        DEALLOCATE(work)
        DEALLOCATE(iwork)
        DEALLOCATE(rwork)
        DEALLOCATE(gap)
        DEALLOCATE(ifail)
        DEALLOCATE(iclustr)
        DEALLOCATE(Z)
        DEALLOCATE(dyn_matrix)

        IF (isroot) THEN
           WRITE(*,*) q_pt
           DO ia=1,3*natom
              WRITE(*,'(A,I5,A,F16.6)') 'w(',ia,') = ', SQRT(ABS(W(ia)))*15.633302*33.35641* &
                                                        SIGN(1.00,W(ia))
           END DO
        END IF
        DEALLOCATE(W)

    END DO


    DEALLOCATE(force_cnst_dist)

    DEALLOCATE(mass)
    DEALLOCATE(real_pos)
    DEALLOCATE(crys)
    DEALLOCATE(at_types_i)
    DEALLOCATE(q_array)


99  CONTINUE

    ! ** Close BLACS Grid

    CALL blacs_barrier(icontext, 'All')
    CALL blacs_gridexit(icontext)
    CALL blacs_exit(0)

    END PROGRAM
