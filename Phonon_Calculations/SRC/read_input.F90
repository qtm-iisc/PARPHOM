! Package: PARPHOM - PARallel PHOnon calculator for Moiré systems
! Authors: Shinjan Mandal, Indrajit Maity, H. R. Krishnamurthy, Manish Jain
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

!> @file read_input.F90
!> @brief Input file parser for PARPHOM phonon calculations
!> @details This file contains routines for reading and validating input parameters,
!>          setting default values, and sanitizing input data for the phonon calculations.

!> @brief Read and process the input file for phonon calculations
!> @details This subroutine reads an input file that contains configuration parameters
!>          for the phonon calculations. It processes the following types of parameters:
!>          - File paths and names for LAMMPS and force constant data
!>          - Calculation parameters (eigenvalues, eigenvectors, velocities)
!>          - Parallel computation settings (block sizes, tolerances)
!>          - Output configuration
!>
!> The input file format uses ':' as a delimiter between parameter names and values.
!> Lines starting with '#' or '!' are treated as comments.
!>
!> @note This routine requires MPI initialization before being called
!> @warning Will exit the program if critical input parameters are missing or invalid
subroutine read_input()
    
    ! Subroutine to read the input file
    ! ----------------------------------

    use global_variables
    implicit none
    
    character(len=char_len) :: input_file
    integer :: error
    character(len=char_len) :: data_ 
    integer :: pos, i
    character(char_len), dimension(2) :: args_
    integer, parameter :: inp_unit = 20
    character(1), external :: capital

    call default_variables()

    call get_command_argument(1,input_file)
    if (len_trim(input_file)==0) then
        err_msg = "\r\n#!#!#!#!#!#!#!#!\r\nError! No input file"
        call error_message()
        call close_mpi()
        call exit
    end if

    call mpi_barrier(mpi_global%comm, mpierr)

    open(unit=inp_unit, file=trim(adjustl(input_file)), iostat=error)

    if (error.ne.0) then
        write(err_msg,'(2A)') 'Error reading input file ', trim(adjustl(input_file))
        call error_message()
        call exit
    end if

    debug_str = ' '
    call debug_output(0)
    write(debug_str,'(2A)') "Reading input file : ", trim(adjustl(input_file))
    call debug_output(0)


    ! -------------------------------------------------------------------------------- !
    !                                                                                  !
    ! Read all lines in the input file.                                                !
    ! If a "#" or "!" is encountered, the line is treated as a comment                 ! 
    ! Ignores empty lines                                                              !
    !                                                                                  !
    ! Splits the lines taking the delimiter as ':'                                     !
    ! The left hand side of the delimiter is used to identify the variable being read  !
    ! The right hand side is stored directly to the necessary variables                !
    !                                                                                  !
    ! -------------------------------------------------------------------------------- !


    do 
100     read(inp_unit,'(A)',iostat = error) data_
        data_ = trim(adjustl(data_))
        if (error.eq.0) then
            if (data_(1:1) == '#'.or.data_(1:1) == '!'.or.data_(1:1) == ' ' ) goto 100 
            
            pos = index(data_, '!')
            if (pos .ne. 0) then
                data_ = data_(1:pos-1)
            end if
            
            pos = index(data_, '#')
            if (pos .ne. 0) then
                data_ = data_(1:pos-1)
            end if

            pos = index(data_,':')
            args_(1) = data_(1:pos-1)
            args_(2) = data_(pos+1:)

            do i=1,len_trim(args_(1))
                args_(1)(i:i) = capital(args_(1)(i:i))
            end do            
            
            select case(trim(args_(1)))
                case ("LAMMPS FILE LOCATION", "LAMMPS_FILE_LOCATION")
                    lammps_file%location = trim(adjustl(args_(2)))
                case ("LAMMPS_FILE_NAME", "LAMMPS FILE NAME")
                    lammps_file%name_ = trim(adjustl(args_(2)))
                case ("NATOM")
                    read(args_(2), *) moire%natom
                case ("ATOM_STYLE", "ATOM STYLE")
                    lammps_file%atom_style = trim(adjustl(args_(2)))
                case ("ATOM_TYPES", "ATOM TYPES")
                    read(args_(2),*) lammps_file%at_types
                case ("FORCE CONSTANT FILE LOCATION", "FORCE_CONSTANT_FILE_LOCATION")
                    force_const%location = trim(adjustl(args_(2)))
                case ("FORCE CONSTANT FILE NAME", "FORCE_CONSTANT_FILE_NAME")
                    force_const%name_ = trim(adjustl(args_(2)))
                case ("FORCE CONSTANT DATASET", "FORCE_CONSTANT_DATASET")
                    force_const%dset_name = trim(adjustl(args_(2)))
                case ("Q_FILE_LOCATION", "Q FILE LOCATION")
                    q_file%location = trim(adjustl(args_(2)))
                case ("Q_FILE_NAME","Q FILE NAME")
                    q_file%name_ = trim(adjustl(args_(2)))
                case ("NQPT")
                    read(args_(2),*) q_file%npt
#ifdef __QPOOL
                case ("NUM_QPOOLS","NUM QPOOLS")
                    read(args_(2), *) num_pools
#endif  
                case ("COMPUTE_EIGVECS","COMPUTE EIGVECS")
                    args_(2) = trim(adjustl(args_(2)))
                    if (args_(2) == 'yes' .or. args_(2) == 'true' .or. &
                        args_(2) == 'Yes' .or. args_(2) == 'True' .or. &
                        args_(2) == 'YES' .or. args_(2) == 'TRUE') then
                        pzheevx_vars%comp_evec = 'V'
                        evec_comp = .true.
                    elseif (args_(2) == 'no' .or. args_(2) == 'false' .or. &
                            args_(2) == 'No' .or. args_(2) == 'False' .or. &
                            args_(2) == 'NO' .or. args_(2) == 'FALSE') then
                        pzheevx_vars%comp_evec = 'N'
                        evec_comp = .false.
                    else
                        write(err_msg, '(3A)') "Could not interpret command ", &
                                               args_(2), &
                                               "\r\nEigenvectors will not be computed."
                        call error_message()
                        pzheevx_vars%comp_evec = 'N'
                        evec_comp = .false.
                    end if
                case ("GROUP VELOCITY", "GROUP_VELOCITY")
                    args_(2) = trim(adjustl(args_(2)))
                    if (args_(2) == 'yes' .or. args_(2) == 'true' .or. &
                        args_(2) == 'Yes' .or. args_(2) == 'True' .or. &
                        args_(2) == 'YES' .or. args_(2) == 'TRUE') then
                        comp_vel = .true.
                    elseif (args_(2) == 'no' .or. args_(2) == 'false' .or. &
                            args_(2) == 'No' .or. args_(2) == 'False' .or. &
                            args_(2) == 'NO' .or. args_(2) == 'FALSE') then
                        comp_vel = .false.
                    else
                        write(err_msg, '(3A)') "Could not interpret command ", &
                                               args_(2), &
                                               "\r\nVelocity will not be computed."
                        call error_message()
                        comp_vel = .false.
                    end if
                case ("VELOCITY METHOD", "VELOCITY_METHOD")
                    args_(2) = trim(adjustl(args_(2)))
                    if (args_(2) == 'cd' .or. args_(2) == 'CD' .or. &
                        args_(2) == 'central_difference' .or.       &
                        args_(2) == 'CENTRAL_DIFFERENCE') then
                        vel_method = 'C'
                    elseif (args_(2) == 'A' .or. args_(2) == 'a' .or. &
                            args_(2) == 'Analytic' .or.               &
                            args_(2) == 'analytic') then
                        vel_method = 'A'
                    else
                        write(err_msg, '(3A)') "Could not interpret command ", &
                                               args_(2), &
                                    "\r\nVelocity will be computed using analytic derivative."
                        call error_message()
                        vel_method = 'A'
                    end if
                case ("RANGE")
                    if (trim(adjustl(args_(2))) == 'A') then
                        pzheevx_vars%range_ = 'A'
                    elseif (trim(adjustl(args_(2))) == 'V') then 
                        pzheevx_vars%range_ = 'V'
                    elseif (trim(adjustl(args_(2))) == 'I') then   
                        pzheevx_vars%range_ = 'I'
                    else 
                        write(err_msg, '(3A)') 'Could not interpret command', &
                                                trim(adjustl(args_(2))), &
                                                '\r\n Setting range to A'
                    end if
                case ("MAX_EIGVAL", "MAX EIGVAL")
                    read(args_(2), *) pzheevx_vars%vu
                case ("MIN_EIGVAL", "MIN EIGVAL")
                    read(args_(2), *) pzheevx_vars%vl
                case ("MAX_INDEX", "MAX INDEX")
                    read(args_(2), *) pzheevx_vars%iu
                case ("MIN_INDEX", "MIN INDEX")
                    read(args_(2), *) pzheevx_vars%il
                case ("ABSTOL")
                    read(args_(2), *) pzheevx_vars%abstol
                case ("ORFAC")
                    read(args_(2), *) pzheevx_vars%orfac 
                case ("MB")
                    read(args_(2),*) pzheevx_vars%mb
                case ("NB")
                    read(args_(2),*) pzheevx_vars%nb
                case ("OUTPUT FILE NAME", "OUTPUT_FILE_NAME")
                    write(output_file_name,'(A)') trim(adjustl(args_(2)))
                case ("OUTPUT FILE LOCATION", "OUTPUT_FILE_LOCATION")
                    write(output_file_location,'(A)') trim(adjustl(args_(2)))
                case ("PRINT PROGRESS","PRINT_PROGRESS")
                    args_(2) = trim(adjustl(args_(2)))
                    if (args_(2) == 'yes' .or. args_(2) == 'true' .or. &
                        args_(2) == 'Yes' .or. args_(2) == 'True' .or. &
                        args_(2) == 'YES' .or. args_(2) == 'TRUE') then
                        print_progress = .true.
                    elseif (args_(2) == 'no' .or. args_(2) == 'false' .or. &
                            args_(2) == 'No' .or. args_(2) == 'False' .or. &
                            args_(2) == 'NO' .or. args_(2) == 'FALSE') then
                        print_progress = .false.
                    else
                        write(err_msg, '(3A)') "Could not interpret command ", &
                                               args_(2), &
                                               "\r\nProgress will be printed."
                        call error_message()
                        print_progress = .true.
                    end if

                case default
                    write(err_msg, '(3A)') "\r\n%%%%%\r\nCommand ", trim(args_(1)), " not recognised"
                    call error_message()
                    err_msg = "Ignoring this command\r\n"
                    call error_message()  
            end select

        elseif (error.gt.0) then
            err_msg = 'Error encountered while reading file'
            call error_message()
            call exit
        else
#ifdef __DEBUG
            debug_str = "\r\nFile read successfully"
            call debug_output(0)
#endif
            exit
        end if
    end do

    close(unit=inp_unit)

    call sanitize_input()

    debug_str = '\r\n========================================================='
    call debug_output(0)
    debug_str = '\r\nCalculations started with the following parameters: '
    call debug_output(0)
    write(debug_str, '( 2A )') "\r\nLAMMPS file location:  ", trim(lammps_file%location)
    call debug_output(0)
    write(debug_str, '(2A)') "LAMMPS file name:  ", trim(lammps_file%name_)
    call debug_output(0)
    write(debug_str, '(A, I0)') "Number of atoms:  ", moire%natom
    call debug_output(0)
    write(debug_str, '(A,I0)') "Atom types:  ", lammps_file%at_types
    call debug_output(0)
    write(debug_str, '(2A)') "LAMMPS atom style:  ", trim(lammps_file%atom_style)
    call debug_output(0)
    write(debug_str, '(2A)') "Q file location:  ", trim(q_file%location)
    call debug_output(0)
    write(debug_str, '(2A)') "Q file name:  ", trim(q_file%name_)
    call debug_output(0)
    write(debug_str, '(A,I0)') "Number of Q-points:  ", q_file%npt
    call debug_output(0)
    write(debug_str, '(2A)') "Force Constant file location:  ", trim(force_const%location)
    call debug_output(0)
    write(debug_str, '(2A)') "Force Constant file name:  ", trim(force_const%name_)
    call debug_output(0)
    write(debug_str, '(2A)') "Force Constant dataset name:  ", trim(force_const%dset_name)
    call debug_output(0)
#ifdef __QPOOL
    write(debug_str, '(A,I0)') "Number of simultaneous Q-point diagonalizations:  ", num_pools
    call debug_output(0)
#endif
    write(debug_str, '(A,L)') "Compute group  velocities: ", comp_vel
    call debug_output(0)
    write (debug_str,'(2A)') "Velocity computation method: ", vel_method
    call debug_output(0)
    write(debug_str, '(A,L)') "Compute Eigenvectors: ", evec_comp
    call debug_output(0)    
    write(debug_str, '(2A)') "Range of the calculation: ", trim(pzheevx_vars%range_)
    call debug_output(0)
    write(debug_str, '(A,F0.6)') "Minimum eigenvalue : ", pzheevx_vars%vl
    call debug_output(0)
    write(debug_str, '(A,F0.6)') "Maximum eigenvalue : ", pzheevx_vars%vu
    call debug_output(0)
    write(debug_str, '(A,I0)') "Minimum eigenvalue index : ", pzheevx_vars%il
    call debug_output(0)
    write(debug_str, '(A,I0)') "Maximum eigenvalue index : ", pzheevx_vars%iu
    call debug_output(0)
    write(debug_str, '(2(A,I0),A)') "Block size for scalapack diagonalization : (", &
                                     pzheevx_vars%mb,',',pzheevx_vars%nb,')'
    call debug_output(0)
    write(debug_str, '(A,E12.4)') "Abstol : ", pzheevx_vars%abstol
    call debug_output(0)
    write(debug_str, '(A,E12.4)') "Orfac : ", pzheevx_vars%orfac
    call debug_output(0)
    write(debug_str, '(2A)') "Output file name : ", trim(output_file_name)
    call debug_output(0)
    write(debug_str, '(2A)') "Output file location : ", trim(output_file_location)
    call debug_output(0)
    write(debug_str, '(A,L)') "Print Progress of calculations: ", print_progress
    call debug_output(0)
    debug_str = '\r\n========================================================='
    call debug_output(0)

    return

end subroutine






!> @brief Convert lowercase character to uppercase
!> @param[in] in_char Input character to be converted
!> @return The uppercase version of the input character if it was lowercase,
!>         otherwise returns the original character unchanged
function capital (in_char)
    implicit none
    character(len=1), intent(in) :: in_char
    character(len=1) :: capital
    character(len=26), parameter :: lower='abcdefghijklmnopqrstuvwxyz', &
                                    upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    integer :: i

    do i=1,26
        if (in_char == lower(i:i)) then
            capital = upper(i:i)
            return
        end if
    end do
    capital = in_char
    return
end function capital



!> @brief Initialize default values for all input parameters
!> @details Sets conservative default values for all configurable parameters
!>          including file paths, calculation options, and numerical parameters.
!>          These values will be used if not explicitly specified in the input file.
subroutine default_variables()
    use global_variables
    implicit none

    lammps_file%location = './'
    lammps_file%name_ = 'lammps.dat'
    lammps_file%atom_style = 'A'
    lammps_file%at_types = 1
    moire%natom = 1
    q_file%location = './'
    q_file%name_ = 'q_file.dat'
    q_file%npt = 1
    q_file%start = 1
    q_file%finish = 1
    evec_comp = .false.
    comp_vel = .false.
    print_progress = .true.
    vel_method = 'A'
    pzheevx_vars%comp_evec = 'N'
    pzheevx_vars%range_ = 'A'
    pzheevx_vars%il = 1
    pzheevx_vars%iu = 1
    pzheevx_vars%vl = 0.0
    pzheevx_vars%vu = 0.0
    pzheevx_vars%mb = -1
    pzheevx_vars%nb = -1
    pzheevx_vars%abstol= 1E-15
    pzheevx_vars%orfac = 1E-15
    force_const%location = './'
    force_const%name_ = 'FC_hdf5'
    force_const%dset_name = 'force_constants'
    output_file_name = 'results'
    output_file_location = './'

    return

end subroutine


!> @brief Validate and adjust input parameters
!> @details Performs the following checks and adjustments:
!> - Ensures file paths end with '/'
!> - Sets appropriate block sizes for parallel calculations
!> - Validates eigenvalue ranges and indices
!> - Ensures compatibility between velocity and eigenvector calculations
!> - Handles output file naming and conflicts
!> @warning Will exit the program if critical parameters are invalid
subroutine sanitize_input()

    use global_variables
    implicit none
    integer :: len_
    character(char_len) :: file_name, full_file_details
    logical :: file_exists
    integer :: counter

    len_ = len_trim(lammps_file%location)

    if (lammps_file%location(len_:len_).ne.'/') then
        write(lammps_file%location,'(2A)') trim(adjustl(lammps_file%location)),'/'
    end if

    len_ = len_trim(q_file%location)

    if (q_file%location(len_:len_).ne.'/') then
        write(q_file%location,'(2A)') trim(adjustl(q_file%location)),'/'
    end if
   
    len_ = len_trim(output_file_location)

    if (output_file_location(len_:len_).ne.'/') then
        write(output_file_location,'(2A)') trim(adjustl(output_file_location)),'/'
    end if

    len_ = len_trim(force_const%location)

    if (force_const%location(len_:len_).ne.'/') then
        write(force_const%location,'(2A)') trim(adjustl(force_const%location)),'/'
    end if
    

    if (pzheevx_vars%mb == -1) then
        pzheevx_vars%mb = min(32,moire%natom)/mpi_global%size_ + 1
    end if

    if (pzheevx_vars%nb == -1) then
        pzheevx_vars%nb = min(32,moire%natom)/mpi_global%size_ + 1
    end if

    select case (pzheevx_vars%range_)
        case("V")
            pzheevx_vars%il = 1
            pzheevx_vars%iu = moire%natom
        case("I")
            pzheevx_vars%vl = 0.0
            pzheevx_vars%vu = 0.0
        case("A")
            pzheevx_vars%il = 1
            pzheevx_vars%iu = moire%natom
            pzheevx_vars%vl = 0.0
            pzheevx_vars%vu = 0.0
    end select

    if ((pzheevx_vars%range_ == 'I').and.(pzheevx_vars%il.gt.pzheevx_vars%iu)) then
        err_msg = "\r\n Invalid input for eigenvalue indices to be computed"
        call error_message()
        call exit
    end if

    if ((pzheevx_vars%range_ == 'V').and.(pzheevx_vars%vl .gt. pzheevx_vars%vu)) then
        err_msg = "\r\n Invalid input for the interval in which eigenvalue is to be computed"
        call error_message()
        call exit
    end if
  
    if (comp_vel) then
        pzheevx_vars%comp_evec='V'
    end if

    write(full_file_details, '(3A)') trim(adjustl(output_file_location)), &
                                     trim(adjustl(output_file_name)),'.hdf5'

    inquire(file=full_file_details, exist=file_exists)

    if (file_exists) then
        write(err_msg,'(5A)') "\r\n File ", trim(adjustl(output_file_name)),'.hdf5', &
                              " already exists at ", trim(adjustl(output_file_location))
        call error_message()

        counter = 0

        do while (file_exists) 
            counter = counter + 1
            write(file_name,'(2A,I0)') trim(adjustl(output_file_name)),'_v',counter
            write(full_file_details, '(3A)') trim(adjustl(output_file_location)), &
                                             trim(adjustl(file_name)),'.hdf5'
            inquire(file=full_file_details, exist=file_exists)
        end do

        write(output_file_name,'(2A)') trim(adjustl(file_name)),'.hdf5'
        write(err_msg,'(2A)') "\r\n Output file being renamed as ", &
                              trim(adjustl(output_file_name))
        call error_message()
    else
        write(output_file_name,'(2A)') trim(adjustl(output_file_name)),'.hdf5'
    end if

    return

end subroutine