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
! Program: read_input (Parse input file and initialize parameters)
!> \file   read_input.F90
!> \brief  Reads and parses the input configuration file for phonon calculations.
!> \details
!>   - Calls `default_variables` to set initial defaults.
!>   - Retrieves and validates the input filename.
!>   - Opens and reads the input file line by line, ignoring comments and blank lines.
!>   - Splits each line at `:` to identify parameter names and values.
!>   - Maps parameter names to global variables via a `select case` statement.
!>   - Handles special parsing for boolean flags and numeric values.
!>   - Performs final sanitation and reports all parameter values via debug output.
subroutine read_input()
    !> Initialize default parameter values
    call default_variables()

    !> Retrieve the input filename from the first command-line argument
    call get_command_argument(1, input_file)
    if (len_trim(input_file) == 0) then
        err_msg = "\r\n#!#!#!#!#!#!#!#!\r\nError! No input file"
        call error_message()
        call close_mpi()
        call exit
    end if

    !> Synchronize all MPI processes before reading
    call mpi_barrier(mpi_global%comm, mpierr)

    !> Open the input file for reading
    open(unit = inp_unit, file = trim(adjustl(input_file)), iostat = error)
    if (error .ne. 0) then
        write(err_msg, '(2A)') 'Error reading input file ', trim(adjustl(input_file))
        call error_message()
        call exit
    end if

    !> Debug: announce start of file reading
    debug_str = ' '
    call debug_output(0)
    write(debug_str,'(2A)') "Reading input file : ", trim(adjustl(input_file))
    call debug_output(0)

    !> Read all lines in the input file, skipping comments and blank lines
    do
100     read(inp_unit, '(A)', iostat = error) data_
        data_ = trim(adjustl(data_))
        if (error == 0) then
            !> Skip comments (#, !, or blank)
            if (data_(1:1) == '#' .or. data_(1:1) == '!' .or. data_(1:1) == ' ') goto 100

            !> Strip inline comments marked by '!' or '#'
            pos = index(data_, '!')
            if (pos > 0) data_ = data_(1:pos-1)
            pos = index(data_, '#')
            if (pos > 0) data_ = data_(1:pos-1)

            !> Split at ':' into name and value
            pos = index(data_, ':')
            args_(1) = data_(1:pos-1)
            args_(2) = data_(pos+1:)

            !> Convert parameter name to uppercase
            do i = 1, len_trim(args_(1))
                args_(1)(i:i) = capital(args_(1)(i:i))
            end do

            !> Assign parameters based on name
            select case (trim(args_(1)))
                case ("LAMMPS FILE LOCATION", "LAMMPS_FILE_LOCATION")
                    lammps_file%location      = trim(adjustl(args_(2)))
                case ("LAMMPS FILE NAME", "LAMMPS_FILE_NAME")
                    lammps_file%name_         = trim(adjustl(args_(2)))
                case ("NATOM")
                    read(args_(2), *) moire%natom
                case ("ATOM STYLE", "ATOM_STYLE")
                    lammps_file%atom_style    = trim(adjustl(args_(2)))
                case ("ATOM TYPES", "ATOM_TYPES")
                    read(args_(2), *) lammps_file%at_types
                case ("FORCE CONSTANT FILE LOCATION", "FORCE_CONSTANT_FILE_LOCATION")
                    force_const%location      = trim(adjustl(args_(2)))
                case ("FORCE CONSTANT FILE NAME", "FORCE_CONSTANT_FILE_NAME")
                    force_const%name_         = trim(adjustl(args_(2)))
                case ("FORCE CONSTANT DATASET", "FORCE_CONSTANT_DATASET")
                    force_const%dset_name     = trim(adjustl(args_(2)))
                case ("Q FILE LOCATION", "Q_FILE_LOCATION")
                    q_file%location           = trim(adjustl(args_(2)))
                case ("Q FILE NAME", "Q_FILE_NAME")
                    q_file%name_              = trim(adjustl(args_(2)))
                case ("NQPT")
                    read(args_(2), *) q_file%npt
#ifdef __QPOOL
                case ("NUM QPOOLS", "NUM_QPOOLS")
                    read(args_(2), *) num_pools
#endif
                case ("COMPUTE EIGVECS", "COMPUTE_EIGVECS")
                    args_(2) = trim(adjustl(args_(2)))
                    if (any([args_(2) == 'yes', args_(2) == 'true', &
                             args_(2) == 'Yes', args_(2) == 'True', &
                             args_(2) == 'YES', args_(2) == 'TRUE'])) then
                        pzheevx_vars%comp_evec = 'V'
                        evec_comp             = .true.
                    elseif (any([args_(2) == 'no', args_(2) == 'false', &
                                 args_(2) == 'No', args_(2) == 'False', &
                                 args_(2) == 'NO', args_(2) == 'FALSE'])) then
                        pzheevx_vars%comp_evec = 'N'
                        evec_comp             = .false.
                    else
                        write(err_msg, '(3A)') "Could not interpret command ", &
                                               args_(2), &
                                               "\r\nEigenvectors will not be computed."
                        call error_message()
                        pzheevx_vars%comp_evec = 'N'
                        evec_comp             = .false.
                    end if
                case ("GROUP VELOCITY", "GROUP_VELOCITY")
                    args_(2) = trim(adjustl(args_(2)))
                    if (any([args_(2) == 'yes', args_(2) == 'true', &
                             args_(2) == 'Yes', args_(2) == 'True', &
                             args_(2) == 'YES', args_(2) == 'TRUE'])) then
                        comp_vel              = .true.
                    elseif (any([args_(2) == 'no', args_(2) == 'false', &
                                 args_(2) == 'No', args_(2) == 'False', &
                                 args_(2) == 'NO', args_(2) == 'FALSE'])) then
                        comp_vel              = .false.
                    else
                        write(err_msg, '(3A)') "Could not interpret command ", &
                                               args_(2), &
                                               "\r\nVelocity will not be computed."
                        call error_message()
                        comp_vel              = .false.
                    end if
                case ("VELOCITY METHOD", "VELOCITY_METHOD")
                    args_(2) = trim(adjustl(args_(2)))
                    if (any([args_(2) == 'cd', args_(2) == 'CD', &
                             args_(2) == 'central_difference', &
                             args_(2) == 'CENTRAL_DIFFERENCE'])) then
                        vel_method            = 'C'
                    elseif (any([args_(2) == 'a', args_(2) == 'A', &
                                 args_(2) == 'analytic', &
                                 args_(2) == 'ANALYTIC'])) then
                        vel_method            = 'A'
                    else
                        write(err_msg, '(3A)') "Could not interpret command ", &
                                               args_(2), &
                                               "\r\nUsing analytic derivative."
                        call error_message()
                        vel_method            = 'A'
                    end if
                case ("RANGE")
                    select case (trim(adjustl(args_(2)))
                        case ('A') pzheevx_vars%range_ = 'A'
                        case ('V') pzheevx_vars%range_ = 'V'
                        case ('I') pzheevx_vars%range_ = 'I'
                        case default
                            write(err_msg, '(3A)') 'Invalid range ', &
                                                   trim(adjustl(args_(2))), &
                                                   '\r\nSetting range to A'
                    end select
                case ("MAX EIGVAL", "MAX_EIGVAL")
                    read(args_(2), *) pzheevx_vars%vu
                case ("MIN EIGVAL", "MIN_EIGVAL")
                    read(args_(2), *) pzheevx_vars%vl
                case ("MAX INDEX", "MAX_INDEX")
                    read(args_(2), *) pzheevx_vars%iu
                case ("MIN INDEX", "MIN_INDEX")
                    read(args_(2), *) pzheevx_vars%il
                case ("ABSTOL")
                    read(args_(2), *) pzheevx_vars%abstol
                case ("ORFAC")
                    read(args_(2), *) pzheevx_vars%orfac
                case ("MB")
                    read(args_(2), *) pzheevx_vars%mb
                case ("NB")
                    read(args_(2), *) pzheevx_vars%nb
                case ("OUTPUT FILE NAME", "OUTPUT_FILE_NAME")
                    write(output_file_name, '(A)') trim(adjustl(args_(2)))
                case ("OUTPUT FILE LOCATION", "OUTPUT_FILE_LOCATION")
                    write(output_file_location, '(A)') trim(adjustl(args_(2)))
                case ("PRINT PROGRESS", "PRINT_PROGRESS")
                    args_(2) = trim(adjustl(args_(2)))
                    if (any([args_(2) == 'yes', args_(2) == 'true', &
                             args_(2) == 'Yes', args_(2) == 'True', &
                             args_(2) == 'YES', args_(2) == 'TRUE'])) then
                        print_progress        = .true.
                    elseif (any([args_(2) == 'no', args_(2) == 'false', &
                                 args_(2) == 'No', args_(2) == 'False', &
                                 args_(2) == 'NO', args_(2) == 'FALSE'])) then
                        print_progress        = .false.
                    else
                        write(err_msg, '(3A)') "Could not interpret command ", &
                                               args_(2), &
                                               "\r\nDefaulting to print progress."
                        call error_message()
                        print_progress        = .true.
                    end if
                case default
                    write(err_msg, '(3A)') "Command not recognised: ", trim(args_(1)), "\r\nIgnoring."
                    call error_message()
            end select

        elseif (error > 0) then
            !> End of file or read error
            err_msg = 'Error encountered while reading file'
            call error_message()
            call exit
        else
#ifdef __DEBUG
            !> Debug: file read completed
            debug_str = "\r\nFile read successfully"
            call debug_output(0)
#endif
            exit
        end if
    end do

    !> Close input and finalize parameters
    close(unit = inp_unit)
    call sanitize_input()

    !> Debug: output final parameter summary
    debug_str = '\r\n========================================================='
    call debug_output(0)
    debug_str = '\r\nCalculations started with the following parameters: '
    call debug_output(0)
    write(debug_str,'(2A)') "LAMMPS file location:  ", trim(lammps_file%location)
    call debug_output(0)
    write(debug_str,'(2A)') "LAMMPS file name:      ", trim(lammps_file%name_)
    call debug_output(0)
    write(debug_str,'(A,I0)') "Number of atoms:       ", moire%natom
    call debug_output(0)
    write(debug_str,'(A,I0)') "Atom types:           ", lammps_file%at_types
    call debug_output(0)
    write(debug_str,'(2A)') "LAMMPS atom style:     ", trim(lammps_file%atom_style)
    call debug_output(0)
    write(debug_str,'(2A)') "Q file location:       ", trim(q_file%location)
    call debug_output(0)
    write(debug_str,'(2A)') "Q file name:           ", trim(q_file%name_)
    call debug_output(0)
    write(debug_str,'(A,I0)') "Number of Q-points:    ", q_file%npt
    call debug_output(0)
#ifdef __QPOOL
    write(debug_str,'(A,I0)') "Number of QPools:      ", num_pools
    call debug_output(0)
#endif
    write(debug_str,'(A,L)') "Compute group velocities: ", comp_vel
    call debug_output(0)
    write(debug_str,'(2A)') "Velocity method:      ", vel_method
    call debug_output(0)
    write(debug_str,'(A,L)') "Compute eigenvectors:  ", evec_comp
    call debug_output(0)
    write(debug_str,'(2A)') "Range:                 ", trim(pzheevx_vars%range_)
    call debug_output(0)
    write(debug_str,'(A,F0.6)') "Min eigenvalue:        ", pzheevx_vars%vl
    call debug_output(0)
    write(debug_str,'(A,F0.6)') "Max eigenvalue:        ", pzheevx_vars%vu
    call debug_output(0)
    write(debug_str,'(A,I0)') "Min eigen index:        ", pzheevx_vars%il
    call debug_output(0)
    write(debug_str,'(A,I0)') "Max eigen index:        ", pzheevx_vars%iu
    call debug_output(0)
    write(debug_str,'(2(A,I0),A)') "Block size (mb, nb):   (", pzheevx_vars%mb, ',', pzheevx_vars%nb, ")"
    call debug_output(0)
    write(debug_str,'(A,E12.4)') "Abstol:                ", pzheevx_vars%abstol
    call debug_output(0)
    write(debug_str,'(A,E12.4)') "Orfac:                 ", pzheevx_vars%orfac
    call debug_output(0)
    write(debug_str,'(2A)') "Output file:           ", trim(output_file_name)
    call debug_output(0)
    write(debug_str,'(2A)') "Output location:       ", trim(output_file_location)
    call debug_output(0)
    write(debug_str,'(A,L)') "Print progress:        ", print_progress
    call debug_output(0)

    return
end subroutine read_input

!> \brief Convert a single character to uppercase.
!> \param[in] in_char  Single character to convert.
!> \return  Uppercase equivalent if alphabetic, otherwise returns input.
function capital(in_char)
    implicit none
    character(len=1), intent(in) :: in_char
    character(len=1) :: capital
    character(len=26), parameter :: lower = 'abcdefghijklmnopqrstuvwxyz', upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    integer :: i
    do i = 1, 26
        if (in_char == lower(i:i)) then
            capital = upper(i:i)
            return
        end if
    end do
    capital = in_char
    return
end function capital

!> \brief Set default input parameter values.
subroutine default_variables()
    use global_variables
    implicit none

    ! Default file locations and names
    lammps_file%location        = './'
    lammps_file%name_           = 'lammps.dat'
    lammps_file%atom_style      = 'A'
    lammps_file%at_types        = 1
    moire%natom                = 1
    q_file%location             = './'
    q_file%name_                = 'q_file.dat'
    q_file%npt                  = 1
    q_file%start                = 1
    q_file%finish               = 1

    ! Computation flags
    evec_comp                   = .false.
    comp_vel                    = .false.
    print_progress              = .true.
    vel_method                  = 'A'
    pzheevx_vars%comp_evec      = 'N'
    pzheevx_vars%range_         = 'A'

    ! Eigenvalue index and range defaults
    pzheevx_vars%il             = 1
    pzheevx_vars%iu             = 1
    pzheevx_vars%vl             = 0.0
    pzheevx_vars%vu             = 0.0

    ! Scalapack block size defaults
    pzheevx_vars%mb             = -1
    pzheevx_vars%nb             = -1
    pzheevx_vars%abstol         = 1E-15
    pzheevx_vars%orfac          = 1E-15

    ! Force constant file settings
    force_const%location        = './'
    force_const%name_           = 'FC_hdf5'
    force_const%dset_name       = 'force_constants'

    ! Output file settings
    output_file_name            = 'results'
    output_file_location        = './'

    return
end subroutine default_variables
