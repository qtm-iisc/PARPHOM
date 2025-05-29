! Package: PARPHOM
! Authors: Shinjan Mandal, Indrajit Maity, H R Krishnamurthy, Manish Jain
! License: GPL-3.0
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!> \file allocate_arrays.F90
!> \brief Allocates distributed arrays for force constants, dynamical matrices, eigenvalues, eigenvectors, and group velocities.
!> \details
!> This subroutine sets up and allocates all distributed arrays required for phonon calculations using ScaLAPACK/BLACS grid.
!> It initializes matrix descriptors, computes local array sizes, and allocates memory for each process. Optionally, it also allocates arrays for group velocities if requested.
!> Debug output is provided for each allocation if compiled with __DEBUG.
!>
!> - Allocates: force_const, dyn_mat, evec, eval, vel (if comp_vel)
!> - Sets up ScaLAPACK descriptors for each distributed matrix
!> - Handles both eigenvector and non-eigenvector cases
!> - Uses MPI barriers for synchronized debug output
!>
!> \author Shinjan Mandal, Indrajit Maity, H R Krishnamurthy, Manish Jain
!> \ingroup phonon_allocation
!>
!> \subroutine allocate_distributed_arrays
!>     Allocates and initializes all distributed arrays for phonon calculations.
!>     Must be called after MPI and BLACS grid are initialized and before any matrix operations.
!>
!> \par Usage
!>     call allocate_distributed_arrays()
!>
!> \par Example
!>     call allocate_distributed_arrays()
!>
!> \note
!>     This routine assumes all global variables and grid structures are initialized.
!>     Debug output is controlled by the __DEBUG preprocessor flag.
!>     Memory allocation is collective across all MPI ranks.
!>     For eigenvector allocation, checks pzheevx_vars%comp_evec == 'V'.
!>     Group velocity arrays are only allocated if comp_vel is true.
!>     All arrays are deallocated elsewhere.
!> 
!> \see main.F90, mod_variables.F90
!> 
!> \warning
!>     Ensure that the number of atoms and BLACS grid are set correctly before calling this routine.
!>     Incorrect grid or atom count may lead to allocation errors or incorrect results.
!> 
!> \remarks
!>     This file is part of the PARPHOM package for phonon calculations.
!> 
!> \copyright GPL-3.0 Shinjan Mandal, Indrajit Maity, H R Krishnamurthy, Manish Jain
!> 
!> -----------------------------------------------------------------------------
subroutine allocate_distributed_arrays()
    !> \brief Main allocation routine for all distributed arrays used in phonon calculations.
    !> \details
    !> Allocates and initializes distributed arrays for force constants, dynamical matrices, eigenvalues, eigenvectors, and group velocities.
    !> Sets up ScaLAPACK descriptors and computes local sizes for each process. Handles both eigenvector and non-eigenvector cases, and optionally allocates group velocity arrays if requested.
    !> Debug output is provided for each allocation if compiled with __DEBUG.
    !>
    !> \note
    !>     Must be called after MPI and BLACS grid are initialized and before any matrix operations.
    !> \see main.F90, mod_variables.F90
    !> \warning Ensure that the number of atoms and BLACS grid are set correctly before calling this routine.
    !> \par Example
    !>     call allocate_distributed_arrays()
    !> 

    use global_variables
    use mpi
    implicit none

    integer :: rsrc, csrc, info
    integer, external :: numroc
#ifdef __DEBUG
    integer :: i
#endif

    rsrc = 0
    csrc = 0

    force_const%locq = numroc(3*moire%natom,pzheevx_vars%mb,grid%mypcol,csrc,grid%npcol)
    force_const%locq = max(1,force_const%locq)
    force_const%lld  = numroc(3*moire%natom,pzheevx_vars%nb,grid%myprow,rsrc,grid%nprow)
    force_const%lld  = max(1,force_const%lld)

    force_const%size_ = force_const%lld*force_const%locq
    allocate(force_const%mat(force_const%size_))

    call descinit(force_const%desca, 3*moire%natom, 3*moire%natom, pzheevx_vars%mb, &
                  pzheevx_vars%nb, rsrc, csrc, grid%context, force_const%lld, info)
#ifdef __DEBUG
    write(debug_str,'(A)') "Force Constant Matrix Allocated"
    call debug_output(0)
    write(debug_str,'(A)') "\r\nForce Constant Parameters: "
    call debug_output(0)
    call mpi_barrier(mpi_global%comm, mpierr)
    do i=0,mpi_global%size_-1
        if (i==mpi_global%rank) then
            write(*,'(A,I8,3(A,I0))') "Rank: ", mpi_global%rank, &
                                   "    |  Size: ", force_const%size_, &
                                   " Local Rows: ", force_const%lld, &
                                   " Local Cols: ", force_const%locq
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do
    call mpi_barrier(mpi_global%comm, mpierr)

    write(debug_str, '(A)') "\r\nForce Constant Matrix descriptor:"
    call debug_output(0)
    do i=0,mpi_global%size_-1
        if (i==mpi_global%rank) then
            write(*,'(A)') " ---------"
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  DTYPE_: ", force_const%desca(DTYPE_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  CTXT_: ", force_const%desca(CTXT_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  M_: ", force_const%desca(M_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  N_: ", force_const%desca(N_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  MB_: ", force_const%desca(MB_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  NB_: ", force_const%desca(NB_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  RSRC_: ", force_const%desca(RSRC_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  CSRC_: ", force_const%desca(CSRC_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  LLD_: ", force_const%desca(LLD_)
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do


#endif

    dyn_mat%locq = numroc(3*moire%natom, pzheevx_vars%mb, grid%mypcol, csrc, grid%npcol)
    dyn_mat%locq = max(dyn_mat%locq,1)
    dyn_mat%lld = numroc(3*moire%natom, pzheevx_vars%nb, grid%myprow, rsrc, grid%nprow)
    dyn_mat%lld = max(dyn_mat%lld,1)

    dyn_mat%size_ = dyn_mat%locq*dyn_mat%lld

    allocate(dyn_mat%mat(dyn_mat%size_))

    call descinit(dyn_mat%desca, 3*moire%natom, 3*moire%natom, pzheevx_vars%mb, &
                  pzheevx_vars%nb, rsrc, csrc, grid%context, dyn_mat%lld, info)

#ifdef __DEBUG
    write(debug_str,'(A)') "Dynamical Matrix Allocated"
    call debug_output(0)
    write(debug_str,'(A)') "\r\nDynamical Matrix Parameters: "
    call debug_output(0)
    call mpi_barrier(mpi_global%comm, mpierr)
    do i=0,mpi_global%size_-1
        if (i==mpi_global%rank) then
            write(*,'(A,I8,3(A,I0))') "Rank: ", mpi_global%rank, &
                                   "    |  Size: ", dyn_mat%size_, &
                                   " Local Rows: ", dyn_mat%lld, &
                                   " Local Cols: ", dyn_mat%locq
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do
    call mpi_barrier(mpi_global%comm, mpierr)

    write(debug_str, '(A)') "\r\nDynamical Matrix descriptor:"
    call debug_output(0)
    do i=0,mpi_global%size_-1
        if (i==mpi_global%rank) then
            write(*,'(A)') " ---------"
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  DTYPE_: ", dyn_mat%desca(DTYPE_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  CTXT_: ", dyn_mat%desca(CTXT_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  M_: ", dyn_mat%desca(M_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  N_: ", dyn_mat%desca(N_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  MB_: ", dyn_mat%desca(MB_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  NB_: ", dyn_mat%desca(NB_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  RSRC_: ", dyn_mat%desca(RSRC_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  CSRC_: ", dyn_mat%desca(CSRC_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  LLD_: ", dyn_mat%desca(LLD_)
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do
#endif
    
    evec%locq = numroc(3*moire%natom, pzheevx_vars%mb, grid%mypcol, csrc, grid%npcol)
    evec%locq = max(evec%locq,1)
    evec%lld = numroc(3*moire%natom, pzheevx_vars%nb, grid%myprow, rsrc, grid%nprow)
    evec%lld = max(evec%lld,1)

    evec%size_ = evec%locq*evec%lld
    call descinit(evec%desca, 3*moire%natom, 3*moire%natom, pzheevx_vars%mb, &
                  pzheevx_vars%nb, rsrc, csrc, grid%context, evec%lld, info)
    
    if (pzheevx_vars%comp_evec=='V') then
        allocate(evec%mat(evec%size_))

#ifdef __DEBUG
    write(debug_str,'(A)') "Eigenvectors Allocated"
    call debug_output(0)
    write(debug_str,'(A)') "\r\nEigenvector Parameters: "
    call debug_output(0)
    call mpi_barrier(mpi_global%comm, mpierr)
    do i=0,mpi_global%size_-1
        if (i==mpi_global%rank) then
            write(*,'(A,I8,3(A,I0))') "Rank: ", mpi_global%rank, &
                                   "    |  Size: ", evec%size_, &
                                   " Local Rows: ", evec%lld, &
                                   " Local Cols: ", evec%locq
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do
    call mpi_barrier(mpi_global%comm, mpierr)

    write(debug_str, '(A)') "\r\nEigenvector descriptor:"
    call debug_output(0)
    do i=0,mpi_global%size_-1
        if (i==mpi_global%rank) then
            write(*,'(A)') " ---------"
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  DTYPE_: ", evec%desca(DTYPE_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  CTXT_: ", evec%desca(CTXT_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  M_: ", evec%desca(M_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  N_: ", evec%desca(N_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  MB_: ", evec%desca(MB_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  NB_: ", evec%desca(NB_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  RSRC_: ", evec%desca(RSRC_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  CSRC_: ", evec%desca(CSRC_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  LLD_: ", evec%desca(LLD_)
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do
#endif
    else
        allocate(evec%mat(1))
    end if

    
    allocate(eval(3*moire%natom))
#ifdef __DEBUG
    write(debug_str,'(A)') "Eigenvalue Allocated"
    call debug_output(0)
#endif

    if (comp_vel) then
        vel%locq = numroc(3*moire%natom, pzheevx_vars%mb, grid%mypcol, csrc, grid%npcol)
        vel%locq = max(vel%locq,1)
        vel%lld = numroc(3*moire%natom, pzheevx_vars%nb, grid%myprow, rsrc, grid%nprow)
        vel%lld = max(vel%lld,1)

        vel%size_ = vel%locq*vel%lld
        call descinit(vel%desca, 3*moire%natom, 3*moire%natom, pzheevx_vars%mb, &
                      pzheevx_vars%nb, rsrc, csrc, grid%context, vel%lld, info)
        allocate(vel%mat(vel%size_))
#ifdef __DEBUG
        write(debug_str,'(A)') "Group Velocity Matrix Allocated"
        call debug_output(0)
        write(debug_str,'(A)') "\r\nGroup Velocity Matrix Parameters: "
        call debug_output(0)
        call mpi_barrier(mpi_global%comm, mpierr)
        do i=0,mpi_global%size_-1
            if (i==mpi_global%rank) then
                write(*,'(A,I8,3(A,I0))') "Rank: ", mpi_global%rank, &
                                       "    |  Size: ", vel%size_, &
                                       " Local Rows: ", vel%lld, &
                                       " Local Cols: ", vel%locq
            end if
            call mpi_barrier(mpi_global%comm, mpierr)
        end do
        call mpi_barrier(mpi_global%comm, mpierr)
    
        write(debug_str, '(A)') "\r\nVelocity Matrix descriptor:"
        call debug_output(0)
        do i=0,mpi_global%size_-1
            if (i==mpi_global%rank) then
                write(*,'(A)') " ---------"
                write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                       "    |  DTYPE_: ", vel%desca(DTYPE_)
                write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                       "    |  CTXT_: ", vel%desca(CTXT_)
                write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                       "    |  M_: ", vel%desca(M_)
                write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                       "    |  N_: ", vel%desca(N_)
                write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                       "    |  MB_: ", vel%desca(MB_)
                write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                       "    |  NB_: ", vel%desca(NB_)
                write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                       "    |  RSRC_: ", vel%desca(RSRC_)
                write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                       "    |  CSRC_: ", vel%desca(CSRC_)
                write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                       "    |  LLD_: ", vel%desca(LLD_)
            end if
            call mpi_barrier(mpi_global%comm, mpierr)
        end do
#endif

    end if

end subroutine
