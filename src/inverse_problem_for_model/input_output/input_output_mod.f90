!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!###################################################################################################################################
!                                                                                                                                  !
!                                                                                                                                  !
!    MODULE input_output,  PURPOSE : READ ACQUISITON FOR FWI                                                                       !
!                                                                                                                                  !
!                                                                                                                                  !
!                   -- read all sources   (tractions, moment, forces, ...)                                                         !
!                   -- read all stations  (list of stations and position for each event )                                          !
!                   -- read all data      (components to be inverted)                                                              !
!                                                                                                                                  !
!    all information is stored in type : acqui                                                                                     !
!                                                                                                                                  !
!                                                                                                                                  !
!                                                                                                                                  !
! NOTE : MPI is based on WORLD define in specfem : my_local_mpi_comm_world and my_local_mpi_comm_for_bcast                         !
!                                                                                                                                  !
!  August 2016                                                                                                                     !
!  Vadim Monteiller, CNRS Marseille, France                                                                                        !
!                                                                                                                                  !
!##################################################################################################################################!

module input_output

  !! IMPORT VARIABLES FROM SPECFEM ---------------------------------------------------------
  use constants, only: mygroup,NDIM,IMAIN,SOURCES_CAN_BE_BURIED

  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS, BROADCAST_SAME_MESH_AND_MODEL, ANISOTROPY, &
                               NSOURCES, NSOURCES_STF, NSTEP, NSTEP_STF

  use specfem_par, only: CUSTOM_REAL, HUGEVAL, NGNOD, NUM_ITER, NPROC, MAX_STRING_LEN, &
                         NGLLX, NGLLY, NGLLZ, NDIM, NSPEC_AB, NGLOB_AB, MIDX, MIDY, MIDZ, &
                         LOCAL_PATH, xigll, yigll, zigll, DT, &
                         ibool, xstore, ystore, zstore, &
                         myrank, USE_SOURCES_RECEIVERS_Z, &
                         USE_FORCE_POINT_SOURCE,USE_EXTERNAL_SOURCE_FILE, USE_RICKER_TIME_FUNCTION

  use specfem_par, only: nrec,islice_selected_rec,ispec_selected_rec, &
                         xi_receiver,eta_receiver,gamma_receiver,nu_rec, &
                         station_name,network_name, &
                         stlat,stlon,stele,stbur, &
                         x_target_station,y_target_station,z_target_station

  use specfem_par, only: hdur,hdur_Gaussian,force_stf,tshift_src,min_tshift_src_original,t0, &
                         islice_selected_source,ispec_selected_source

  use specfem_par_elastic, only: ispec_is_elastic
  use specfem_par_acoustic, only: ispec_is_acoustic
  !!----------------------------------------------------------------------------------------

  use inverse_problem_par
  use mesh_tools
  use IO_model
  use Teleseismic_IO_mod
  use rotations_mod

  implicit none

  PUBLIC  :: init_input_output_mod, SetUpInversion, get_mode_running, dump_adjoint_sources, write_bin_sismo_on_disk, &
             WriteOutputs

  PRIVATE :: read_acqui_file, read_inver_file, store_default_acqui_values, is_blank_line,& !bcast_all_acqui, &
             get_stations, read_data_gather, create_name_database_inversion, read_and_distribute_events_for_simultaneous_runs

  !! DEFINITION OF PRIVATE VARIABLES
  integer,                PRIVATE                                   :: NEVENT, NIFRQ, ncomp_read, ncomp_inv
  real(kind=CUSTOM_REAL), PRIVATE, dimension(:), allocatable        :: fl, fh
  real(kind=CUSTOM_REAL), PRIVATE                                   :: nb_traces_tot
  logical,                PRIVATE                                   :: use_band_pass_filter
  logical, dimension(8) :: is_component_read, is_component_inv
  integer, dimension(8) :: id_component_read, id_component_inv
contains

!################################## ROUTINES THAT CAN BE CALLED FROM MAIN  #########################################################

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------------------------------------------------------------------------------------------------------
!> \brief read input data, sources, receiver and store all in acqui_simu and send to the concerned MPI slices
!-------------------------------------------------------------------------------------------------------------------
  subroutine init_input_output_mod(inversion_param, acqui_simu, myrank)

    integer,                                 intent(in)    ::  myrank
    type(acqui),  dimension(:), allocatable, intent(inout) ::  acqui_simu
    type(inver),                             intent(inout) ::  inversion_param
    ! locals
    integer                                   :: ievent
    character(len=MAX_LEN_STRING)             :: name_file
    character(len=MAX_LEN_STRING)             :: acqui_file, inver_file
    real(kind=CUSTOM_REAL)                    :: elemsize_min_glob,elemsize_max_glob
    real(kind=CUSTOM_REAL)                    :: distance_min_glob,distance_max_glob

    if (myrank == 0) then
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*) '*********************************************'
      write(INVERSE_LOG_FILE,*) '***       READING INPUT PARAMETERS        ***'
      write(INVERSE_LOG_FILE,*) '*********************************************'
      write(INVERSE_LOG_FILE,*)
      call flush_iunit(INVERSE_LOG_FILE)
    endif

    acqui_file = inversion_param%input_acqui_file
    inver_file = inversion_param%input_inver_file

    !! select input file that we want to read ----------------------------------------
    select case(trim(adjustl(type_input)))
    case('exploration')

      call read_acqui_file(acqui_file, acqui_simu, myrank)
      call read_inver_file(inver_file, inversion_param, myrank)
      call get_stations(acqui_simu)     !! stations are based in specfem format
      call get_point_source(acqui_simu) !! sources can be read in specfem format

    case('teleseismic')

      !! this is purely teleseismic case : we don not have sources in
      !! mesh so we don't need to read and locate sources
      call read_acqui_teleseismic_file(acqui_file, acqui_simu, myrank)
      call read_inver_file(inver_file, inversion_param, myrank)
      call setup_teleseismic_stations(acqui_simu, myrank)

    case default

      write(*,*) ' Abort : Invalid input type :', trim(adjustl(type_input))
      stop

    end select
    !! -------------------------------------------------------------------------------

    if (myrank == 0) call flush_iunit(INVERSE_LOG_FILE)

    !! no need to read data for only forward simulation
    if (.not. inversion_param%only_forward) then

      select case (trim(adjustl(type_input)))
      case('teleseismic')
        call read_pif_data_gather(acqui_simu, inversion_param, myrank)
        inversion_param%nb_traces_tot = nb_traces_tot
      case default
        call read_data_gather(acqui_simu, inversion_param, myrank)
        inversion_param%nb_traces_tot = nb_traces_tot
      end select

    endif

    !! create name for outputs
    do ievent = 1,acqui_simu(1)%nevent_tot
      call create_name_database_inversion(acqui_simu(ievent)%prname_inversion, myrank, ievent, LOCAL_PATH)
    enddo

    if (myrank == 0) call flush_iunit(INVERSE_LOG_FILE)

    !!! We can read input model from external files instead of get model stored in specfem databases -----------------------------
    !!! otherwise the model read from mesher cubit or meshfem3D is kept.
    if (inversion_param%input_SEM_model .and. inversion_param%input_FD_model) then

      if (myrank == 0) then
        write(*,*) ' ERROR : You want to read input model both in SEM and FD grid '
        write(*,*) ' Please choose one of the solution  '
      endif
      stop

    else

      if (inversion_param%input_FD_model) then

        if (myrank == 0) then
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*) '*********************************************'
          write(INVERSE_LOG_FILE,*) '***      READING EXTERNAL FD MODEL        ***'
          write(INVERSE_LOG_FILE,*) '*********************************************'
          write(INVERSE_LOG_FILE,*)
          call flush_iunit(INVERSE_LOG_FILE)
        endif

        call ReadInputFDmodel(inversion_param)

      endif

      if (inversion_param%input_SEM_model) then

        if (myrank == 0) then
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*) '*********************************************'
          write(INVERSE_LOG_FILE,*) '***      READING EXTERNAL SEM MODEL       ***'
          write(INVERSE_LOG_FILE,*) '*********************************************'
          write(INVERSE_LOG_FILE,*)
          call flush_iunit(INVERSE_LOG_FILE)
        endif

        call ReadInputSEMmodel(inversion_param)

      endif

      if (inversion_param%input_SEM_prior) then

        if (myrank == 0) then
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*) '*********************************************'
          write(INVERSE_LOG_FILE,*) '***   READING EXTERNAL SEM PRIOR MODEL    ***'
          write(INVERSE_LOG_FILE,*) '*********************************************'
          write(INVERSE_LOG_FILE,*)
          call flush_iunit(INVERSE_LOG_FILE)
        endif

        call ReadInputSEMpriormodel(inversion_param)

      endif

    endif

    !! get dimension of mesh
    call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB, &
                              ibool,xstore,ystore,zstore, &
                              inversion_param%xmin,inversion_param%xmax, &
                              inversion_param%ymin,inversion_param%ymax, &
                              inversion_param%zmin,inversion_param%zmax, &
                              elemsize_min_glob,elemsize_max_glob, &
                              distance_min_glob,distance_max_glob)

    !!-------------------------------------------------------------------------------------------------------------------------

    !!! write what we read to debug : check what is stored in type acqui_simu
    if (DEBUG_MODE) then
      ! check what I read
      write(IIDD,*)
      write(IIDD,*) '######################################################'
      write(IIDD,*)
      write(IIDD,*) 'MYRANK :' ,myrank, ' READ : '
      do ievent = 1,acqui_simu(1)%nevent_tot
        write(IIDD,*)
        write(IIDD,*) 'event ',ievent
        write(IIDD,*) trim(acqui_simu(ievent)%event_name)
        write(IIDD,*) trim(acqui_simu(ievent)%source_file)
        write(IIDD,*) trim(acqui_simu(ievent)%station_file)
        write(IIDD,*) trim(acqui_simu(ievent)%data_file_gather)
        write(IIDD,*) 'total station     :', acqui_simu(ievent)%nsta_tot
        write(IIDD,*) 'stations in slice :',acqui_simu(ievent)%nsta_slice
        if (acqui_simu(ievent)%nsta_slice > 0 .and. .not. inversion_param%only_forward) then
          write(IIDD,*) 'Check data stored : ', acqui_simu(ievent)%nsta_slice, acqui_simu(ievent)%Nt_data
          write(name_file,'(a21,i8.8,a5,i8.8)') 'Check_read_data_event',ievent,'_rank',myrank
          write(IIDD,*) 'Data that was read are dumped in file for checking: ',trim(name_file)
          open(IINN,file=trim(name_file),access='direct', &
              recl=CUSTOM_REAL*acqui_simu(ievent)%nsta_slice*acqui_simu(ievent)%Nt_data)
          write(IINN,rec=1)  acqui_simu(ievent)%data_traces(:,:,1)
          write(IINN,rec=2)  acqui_simu(ievent)%data_traces(:,:,2)
          write(IINN,rec=3)  acqui_simu(ievent)%data_traces(:,:,3)
          close(IINN)
        endif
      enddo
      write(IIDD,*) '######################################################'
      write(IIDD,*)
      call flush_iunit(IIDD)
    endif

  end subroutine init_input_output_mod
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------------------------------------------------
!> read input inversion parameter file
!--------------------------------------------------------------------------------------------------------------------
  subroutine SetUpInversion(inversion_param, myrank)
    type(inver),                         intent(inout)    :: inversion_param
    integer,                             intent(in)       :: myrank

    ! local
    character(len=MAX_LEN_STRING)                         :: acqui_file_ref
    character(len=MAX_LEN_STRING)                         :: name_file
    integer                                               :: sizeprocs
    integer                                               :: ier

    ! MPI process output
    call world_size(sizeprocs)
    call synchronize_all()

    ! acquisition file
    if (NUMBER_OF_SIMULTANEOUS_RUNS > 1) then
      write(prefix_to_path,"('run',i4.4,'/')") mygroup + 1
      inversion_param%input_acqui_file = trim(prefix_to_path)//'/DATA/inverse_problem/acquisition.dat'
      acqui_file_ref = './DATA/inverse_problem/acquisition.dat'
      if (myrank == 0) then
        write(*,*) ' DISTRIBUTION OF EVENTS '
        call flush_iunit(6)
        !! only one process must do I/O (myrank=0, mygroup=0)
        if (mygroup == 0) call read_and_distribute_events_for_simultaneous_runs(NUMBER_OF_SIMULTANEOUS_RUNS, acqui_file_ref)
      endif
    else
      prefix_to_path = './'
      inversion_param%input_acqui_file = './DATA/inverse_problem/acquisition.dat'
    endif

    ! inversion parameter file
    inversion_param%input_inver_file = './DATA/inverse_problem/inversion_fwi.dat'

    ! user output
    if (myrank == 0) then
      ! Log output file
      if (INVERSE_LOG_FILE /= ISTANDARD_OUTPUT) then
        open(INVERSE_LOG_FILE,file=trim(prefix_to_path)//'output_inverse.txt',status='unknown',action='write',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening file output_inverse.txt for writing output info')
      endif

      ! header
      call header_output(INVERSE_LOG_FILE)

      if (.not. inversion_param%only_forward) then
        ! Iteration output
        open(OUTPUT_ITERATION_FILE,file=trim(prefix_to_path)//'output_iterations.txt')
        ! header
        write(OUTPUT_ITERATION_FILE,'(a7,"|",a15,"|",a15,"|",a12,"|",a12,"|",a15)') &
            "# iter "," cost function "," gradient norm ", " relat cost "," relat grad ", " nb line search"
        call flush_iunit(OUTPUT_ITERATION_FILE)

        ! FWI Log output file
        open(OUTPUT_FWI_LOG,file=trim(prefix_to_path)//'output_fwi_log.txt')
      endif
    endif

    ! debug file for each process
    if (DEBUG_MODE) then
      write(name_file,'("debug_file_input_output_module_rank",i6.6)') myrank
      open(IIDD,file=trim(prefix_to_path)//trim(name_file))
      write(IIDD,*) '///////////////////////////////////////////////////////////////////////'
      write(IIDD,*) 'Debug output'
      write(IIDD,*) '///////////////////////////////////////////////////////////////////////'
      write(IIDD,*)
      ! debug output rank/group ids
      write(IIDD,*) 'SETUP INVERSION : rank : ', myrank, ' group : ', mygroup
      write(IIDD,*)
      write(IIDD,*) 'number of simultaneous runs : ', NUMBER_OF_SIMULTANEOUS_RUNS
      write(IIDD,*)
      call flush_iunit(IIDD)
    endif

    !! to be sure that all processes in the MPI_COMM_WORLD wait for the end of deining new inputs files
    call synchronize_all_world()

  end subroutine SetUpInversion
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------------------------------------------------
! read acqui_file, count events then distribute over groups and rewrite suitable acqui_file for each group
!--------------------------------------------------------------------------------------------------------------------
  subroutine read_and_distribute_events_for_simultaneous_runs(NUMBER_OF_SIMULTANEOUS_RUNS, acqui_file_ref)

    character(len=MAX_LEN_STRING),       intent(in)       :: acqui_file_ref
    integer,                             intent(in)       :: NUMBER_OF_SIMULTANEOUS_RUNS
    ! local
    integer                                               :: number_of_events_in_acqui_file_ref
    integer                                               :: nevent_per_group, nevent_remained
    integer                                               :: igroup, ievent_in_group, ievent, ievent_global, ier
    integer,  dimension(:),     allocatable               :: nevent_in_group
    character(len=MAX_LEN_STRING)                         :: line, prefix_to_path_tmp


    write(*,*)
    write(*,*)  ' NUMBER OF SIMULTANEOUS RUNS > 0 '
    write(*,*)
    call flush_iunit(6)
    number_of_events_in_acqui_file_ref=0
    open(666,file=trim(acqui_file_ref))
    do
       read(666,'(a)',end=99) line
       !! no significant line
       if (is_blank_line(line)) cycle
       !! new event
       if (INDEX(line,'event_name') > 0) number_of_events_in_acqui_file_ref=number_of_events_in_acqui_file_ref+1
    enddo
99  close(666)


    nevent_per_group = number_of_events_in_acqui_file_ref / NUMBER_OF_SIMULTANEOUS_RUNS
    nevent_remained =  mod( number_of_events_in_acqui_file_ref, NUMBER_OF_SIMULTANEOUS_RUNS)

    allocate(nevent_in_group(NUMBER_OF_SIMULTANEOUS_RUNS),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 392')
    do ievent=1,NUMBER_OF_SIMULTANEOUS_RUNS
       if (ievent <= nevent_remained) then
          nevent_in_group(ievent)= nevent_per_group+1
       else
          nevent_in_group(ievent)= nevent_per_group
       endif
    enddo

    ievent_global = 0
    ievent_in_group=0
    igroup=1
    open(666,file=trim(acqui_file_ref))
    write(prefix_to_path_tmp,"('run',i4.4,'/')") igroup
    open(777, file=trim(prefix_to_path_tmp)//'DATA/inverse_problem/acquisition.dat')
    do
       read(666,'(a)',end=999) line
       !! no significant line
       if (is_blank_line(line)) cycle

       !! new event
       if (INDEX(line,'event_name') > 0) then
           ievent_global = ievent_global + 1
           write(*,*)
           write(*,*) '   next event ', ievent_global
           ievent_in_group=ievent_in_group+1

       endif

       !! write lines related to the current event
       if (ievent_in_group > nevent_in_group(igroup)) then
          igroup=igroup+1
          write(*,*) ' group ', igroup
          write(prefix_to_path_tmp,"('run',i4.4,'/')") igroup
          close(777)
          open(777, file=trim(prefix_to_path_tmp)//'DATA/inverse_problem/acquisition.dat')
          ievent_in_group=1
       endif
       write(777, '(a)') trim(line)
       write(*,*) trim(line)
       call flush_iunit(6)
    enddo
999  close(666)
    close(777)
    deallocate(nevent_in_group)

  end subroutine read_and_distribute_events_for_simultaneous_runs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------------------------------------------------
!  ! create the name of the database for inverse problem
!--------------------------------------------------------------------------------------------------------------------
  subroutine create_name_database_inversion(prname,iproc,ievent,LOCAL_PATH)

    implicit none

    integer,                       intent(in)   :: iproc, ievent
    ! name of the database file
    character(len=MAX_LEN_STRING), intent(inout) :: prname
    character(len=MAX_STRING_LEN), intent(in   ) :: LOCAL_PATH
    ! local
    character(len=MAX_STRING_LEN)                :: procname

    ! create the name for the database of the current slide and region
    write(procname,"('/proc',i6.6,'_',i6.6,'_')") iproc,ievent
    ! create full name with path
    prname = LOCAL_PATH(1:len_trim(LOCAL_PATH)) // procname

  end subroutine create_name_database_inversion

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------------------------------------------------
!  ! write final outputs model in disk
!--------------------------------------------------------------------------------------------------------------------
  subroutine WriteOutputs(inversion_param)

    type(inver),                                                intent(in)    :: inversion_param

    !! FOR NOT ONLY ONE OUPTPUT BUT FOR FURTHER DEV WE WILL ADD OTHER

    !! writing the final solution in SEM mesh
    call WriteOutputSEMmodel(inversion_param)

  end subroutine WriteOutputs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------------------------------------------------
!  ! get which option we want to run FWI or just direct problem
!--------------------------------------------------------------------------------------------------------------------
  subroutine get_mode_running(mode_running, inversion_param)

    use my_mpi
    use specfem_par, only: myrank

    type(inver),                                    intent(inout) :: inversion_param
    character(len=MAX_LEN_STRING),                  intent(inout) :: mode_running
    integer                                                       :: ier
    character(len=MAX_LEN_STRING)                                 :: arg

    ! gets program arguments
    if (myrank == 0) then
      call get_command_argument(1, mode_running)
      if (trim(mode_running) == '') then
          print *,'Usage: ./xinverse_problem_for_model mode( == forward or l-bfgs) [type]'
          print *
          print *,'  where'
          print *,'      mode - run mode for inverse problem, must be either "forward" or "l-bfgs"'
          print *,'      type - (optional) sets type of inverse problem to "exploration" or "teleseismic" [default == exploration]'
          print *
          stop 'Reenter command line options'
        endif
    endif

    call MPI_BCAST(mode_running,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)

    select case(trim(adjustl(mode_running)))
    case('forward')
      inversion_param%only_forward = .true.

    case('l-bfgs')
      inversion_param%only_forward = .false.

    case default
      print *,'Error: called with invalid run mode argument : ', mode_running
      print *
      print *,'       usage: ./xinverse_problem_for_model mode( == forward or l-bfgs) [type == exploration or teleseismic]'
      print *
      stop 'Invalid program argument'
    end select

    !! check if there is something else
    if (myrank == 0) then
      call get_command_argument(2, arg)
    endif
    call MPI_BCAST(arg,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)

    if (len_trim(arg) > 0) then
      type_input = trim(adjustl(arg))
      select case(type_input)
      case('exploration')
        ! ok
      case('teleseismic')
        ! ok
      case default
        print *,'Error: called with invalid type argument : ', type_input
        print *
        print *,'       usage: ./xinverse_problem_for_model mode( == forward or l-bfgs) [type == exploration or teleseismic]'
        print *
        stop 'Invalid program argument'
      end select
    endif

    ! user output
    if (myrank == 0) then
      write(*,*)
      write(*,*) ' INVERSION       : mode = ', trim(mode_running)
      write(*,*) '                   type = ', trim(type_input)
      write(*,*)
      call flush_iunit(6)
    endif

  end subroutine get_mode_running

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------
! main write adjoint sources gather to check
!----------------------------------------------------------------

  subroutine dump_adjoint_sources(ievent, acqui_simu, myrank)

    integer,                    intent(in)    :: ievent, myrank
    type(acqui),  dimension(:), intent(inout) :: acqui_simu
    ! local
    character(len=MAX_LEN_STRING)             :: name_file_tmp, ch_to_add
    integer                                   :: global_iter

    ! combined index for current frequency stage and inversion iteration
    global_iter = iter_frq * 10000  + iter_inverse

    write(ch_to_add,'(a5,i5.5,a4)') '_iter',global_iter,'_adj'
    name_file_tmp = trim(acqui_simu(ievent)%data_file_gather)//trim(adjustl(ch_to_add))

    if (myrank == 0) then
      write(INVERSE_LOG_FILE,*) ' ... Writing adjoint sources gather for event : ', ievent
    endif

    call  write_bin_sismo_on_disk(ievent, acqui_simu, acqui_simu(ievent)%adjoint_sources,  name_file_tmp, myrank)

  end subroutine dump_adjoint_sources
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------
! main write synthetics gather to check
!----------------------------------------------------------------
  subroutine dump_seismograms(ievent, array_to_write,  acqui_simu, myrank)

    integer,                                   intent(in)    :: ievent, myrank
    real(kind=CUSTOM_REAL),  dimension(:,:,:), intent(in)    :: array_to_write
    type(acqui),  dimension(:),                intent(inout) :: acqui_simu
    ! local
    character(len=MAX_LEN_STRING)                            :: name_file_tmp, ch_to_add
    integer                                                  :: global_iter

    ! combined index for current frequency stage and inversion iteration
    global_iter = iter_frq * 10000  + iter_inverse

    if (iter_frq == 1 .and. iter_inverse == 0 .and. iter_wolfe == 0) then
      ! initial model seismograms
      write(ch_to_add,'(a12)') '_initial_syn'
    else
      ! updated model seismograms
      ! single file per sub-iteration
      !write(ch_to_add,'(a5,i5.5,a6,i4.4,a4)') '_iter',global_iter,'_wolfe',iter_wolfe,'_syn'
      ! overrides sub-iterations
      write(ch_to_add,'(a5,i5.5,a4)') '_iter',global_iter,'_syn'
    endif
    name_file_tmp = trim(acqui_simu(ievent)%data_file_gather)//trim(adjustl(ch_to_add))

    if (myrank == 0) then
      write(INVERSE_LOG_FILE,*) ' ... Writing synthetic data gather for event  : ', ievent
    endif

    call  write_bin_sismo_on_disk(ievent, acqui_simu, array_to_write, name_file_tmp, myrank)

  end subroutine dump_seismograms
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------
! main write filtered data gather to check
!----------------------------------------------------------------
  subroutine dump_filtered_data(ievent, array_to_write,  acqui_simu, myrank)

    integer,                                   intent(in)    :: ievent,myrank
    real(kind=CUSTOM_REAL),  dimension(:,:,:), intent(in)    :: array_to_write
    type(acqui),  dimension(:),                intent(inout) :: acqui_simu
    ! local
    character(len=MAX_LEN_STRING)                            :: name_file_tmp, ch_to_add
    integer                                                  :: global_iter

    ! combined index for current frequency stage and inversion iteration
    global_iter = iter_frq * 10000  + iter_inverse

    write(ch_to_add,'(a5,i5.5,a4)') '_iter',global_iter,'_fil'
    name_file_tmp = trim(acqui_simu(ievent)%data_file_gather)//trim(adjustl(ch_to_add))

    if (myrank == 0) then
      write(INVERSE_LOG_FILE,*) ' ... Writing filtered data gather for event   : ', ievent
    endif

    call  write_bin_sismo_on_disk(ievent, acqui_simu, array_to_write, name_file_tmp, myrank)

  end subroutine dump_filtered_data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------
! main write waveform synthetic data gather (need to choose between write_bin_sismo_on_disk or write_gather_on_disk)
!----------------------------------------------------------------
  subroutine write_bin_sismo_on_disk(ievent, acqui_simu, array_to_write, name_file_to_write, myrank)

    use my_mpi             !! module from specfem
    include "precision.h"  !! from specfem

    integer,                                   intent(in)    :: myrank, ievent
    character(len=MAX_LEN_STRING),             intent(in)    :: name_file_to_write
    real(kind=CUSTOM_REAL),  dimension(:,:,:), intent(in)    :: array_to_write
    type(acqui),             dimension(:),     intent(inout) :: acqui_simu
    ! local
    integer                                                  :: idim, NSTA, NSTA_LOC, Nt, irec, irec_local
    integer                                                  :: tag, ier, nsta_irank, irank, icomp
    real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable    :: Gather, Gather_loc
    integer,                dimension(:),     allocatable    :: irec_global
    integer                                                  :: status(MPI_STATUS_SIZE)

    if (myrank == 0) then
      NSTA = acqui_simu(ievent)%nsta_tot
      Nt = acqui_simu(ievent)%Nt_data
      allocate(Gather(NSTA,Nt,NDIM),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 393')
    endif
    ! not sure if need this sync
    call synchronize_all()

    do irank = 1, NPROC-1

      if (myrank == 0) then
        ! count the receiver in slice irank
        nsta_irank = 0
        do irec = 1,  NSTA
          if (acqui_simu(ievent)%islice_selected_rec(irec) == irank) nsta_irank = nsta_irank + 1
        enddo
        if (nsta_irank > 0) then
          !! data to receive
          allocate(Gather_loc(nsta_irank,Nt,NDIM),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 394')
          allocate(irec_global(nsta_irank),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 395')

          irec_local = 0
          tag   = MPI_ANY_TAG
          call MPI_RECV(Gather_loc, Nt*nsta_irank*NDIM, CUSTOM_MPI_TYPE, irank, tag, my_local_mpi_comm_world, status,  ier)
          call MPI_RECV(irec_global, nsta_irank, MPI_INTEGER, irank, tag, my_local_mpi_comm_world, status,  ier)

          do icomp=1,NDIM
            do irec_local = 1, nsta_irank
              Gather(irec_global(irec_local), :, icomp) = Gather_loc(irec_local, :, icomp)
            enddo
          enddo

          deallocate(Gather_loc)
          deallocate(irec_global)
        endif
      else
        if (myrank == irank .and. acqui_simu(ievent)%nsta_slice > 0) then
          NSTA_LOC = acqui_simu(ievent)%nsta_slice
          Nt = acqui_simu(ievent)%Nt_data
          allocate(Gather_loc(NSTA_LOC,Nt,NDIM),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 396')
          allocate(irec_global(NSTA_LOC),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 397')

          do irec_local = 1, NSTA_LOC
            irec_global(irec_local) = acqui_simu(ievent)%number_receiver_global(irec_local)
            !! choose the rigth seismograms_*
            do icomp = 1,NDIM
              select case (trim(acqui_simu(ievent)%component(icomp)))
              case ('PR')
                Gather_loc(irec_local,:,icomp) = array_to_write(1,irec_local,:)
              case ('UX')
                Gather_loc(irec_local,:,icomp) = array_to_write(1,irec_local,:)
              case ('UY')
                Gather_loc(irec_local,:,icomp) = array_to_write(2,irec_local,:)
              case ('UZ')
                Gather_loc(irec_local,:,icomp) = array_to_write(3,irec_local,:)
              end select
            enddo
          enddo

          tag    = 2001
          call MPI_SEND(Gather_loc,  Nt*NSTA_LOC*NDIM, CUSTOM_MPI_TYPE, 0, tag, my_local_mpi_comm_world, ier)
          call MPI_SEND(irec_global, NSTA_LOC, CUSTOM_MPI_TYPE, 0, tag, my_local_mpi_comm_world, ier)

          deallocate(Gather_loc)
          deallocate(irec_global)
        endif

      endif

      ! not sure if need this sync
      call synchronize_all()

    enddo

    !!  write gather file
    if (myrank == 0) then

      do icomp = 1,NDIM
        do irec_local = 1, acqui_simu(ievent)%nsta_slice
          !! choose the rigth seismograms_*
          select case (trim(acqui_simu(ievent)%component(icomp)))
          case ('PR')
            Gather(acqui_simu(ievent)%number_receiver_global(irec_local),:,icomp) = array_to_write(1,irec_local,:)
          case ('UX')
            Gather(acqui_simu(ievent)%number_receiver_global(irec_local),:,icomp) = array_to_write(1,irec_local,:)
          case ('UY')
            Gather(acqui_simu(ievent)%number_receiver_global(irec_local),:,icomp) = array_to_write(2,irec_local,:)
          case ('UZ')
            Gather(acqui_simu(ievent)%number_receiver_global(irec_local),:,icomp) = array_to_write(3,irec_local,:)
          end select
        enddo
      enddo

      open(IINN,file=trim(adjustl(name_file_to_write)), access='direct', recl=CUSTOM_REAL*Nt*NSTA, status='replace')

      !! write only the asked component or pressure
      irec = 0
      do idim = 1,NDIM
        select case (trim(acqui_simu(ievent)%component(idim)))
        case('UX', 'UY', 'UZ', 'PR')
          irec = irec+1
          write(IINN,rec=irec) Gather(:,:,idim)
        end select
      enddo

      close(IINN)
      deallocate(Gather)
    endif

  end subroutine write_bin_sismo_on_disk

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------
! main read waveform data gather and bcast to MPI slice concerned
!----------------------------------------------------------------
  subroutine read_data_gather(acqui_simu, inversion_param, myrank)

    use my_mpi             !! module from specfem
    include "precision.h"  !! from specfem

    integer,                                     intent(in)    :: myrank
    type(acqui), dimension(:),                   intent(inout) :: acqui_simu
    type(inver),                                 intent(inout) :: inversion_param
    ! local
    integer                                                    :: ievent, idim, NSTA, NSTA_LOC, Nt, irec, irec_local
    integer                                                    :: tag, ier, nsta_irank, irank
    real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable      :: Gather, Gather_loc
    integer                                                    :: status(MPI_STATUS_SIZE)
    real(kind=CUSTOM_REAL)                                     :: gather_maximum, dummy_real

    nb_traces_tot = 0._CUSTOM_REAL

    if (myrank == 0) then
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*) ' READING data'
    endif

    do ievent = 1, acqui_simu(1)%nevent_tot

      if (myrank == 0) then

        NSTA = acqui_simu(ievent)%nsta_tot
        Nt = acqui_simu(ievent)%Nt_data

        allocate(Gather(NSTA,Nt,NDIM),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 398')
        Gather(:,:,:) = 0._CUSTOM_REAL
        ! read gather file
        open(IINN,file=trim(adjustl(acqui_simu(ievent)%data_file_gather)), access='direct', &
             recl=CUSTOM_REAL*Nt*NSTA,status='old')
        !! read only the asked component or pressure
        irec = 0
        do idim = 1,NDIM
          select case (trim(acqui_simu(ievent)%component(idim)))
          case('UX', 'UY', 'UZ', 'PR')
            nb_traces_tot = nb_traces_tot + NSTA
            irec = irec + 1
            read(IINN,rec=irec) Gather(:,:,idim)
          end select
        enddo
        close(IINN)

        ! takes maximum for weigthing
        gather_maximum = maxval(abs(Gather))
        if (gather_maximum == 0._CUSTOM_REAL) gather_maximum = 1._CUSTOM_REAL

        !! store data gather in my slice if needed
        NSTA_LOC = acqui_simu(ievent)%nsta_slice

        allocate(acqui_simu(ievent)%data_traces(NSTA_LOC,Nt,NDIM),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 399')
        acqui_simu(ievent)%data_traces(:,:,:) = 0._CUSTOM_REAL

        allocate(acqui_simu(ievent)%adjoint_sources(NDIM, NSTA_LOC, Nt),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 400')
        acqui_simu(ievent)%adjoint_sources(:,:,:) = 0._CUSTOM_REAL

        !! here weight_trace is allocated with nt = 1 and for local stations only
        ! (will be set in adjoint source creation)
        allocate(acqui_simu(ievent)%weight_trace(NDIM, NSTA_LOC, 1),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 401')
        acqui_simu(ievent)%weight_trace(:,:,:) = 1._CUSTOM_REAL

        ! array for storing filtered data output
        if (VERBOSE_MODE .or. DEBUG_MODE) then
          allocate(acqui_simu(ievent)%synt_traces(NDIM, NSTA_LOC, Nt),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 402')
          acqui_simu(ievent)%synt_traces(:,:,:) = 0._CUSTOM_REAL
        endif

        irec_local = 0
        do irec = 1, NSTA
          if (acqui_simu(ievent)%islice_selected_rec(irec) == myrank) then
            irec_local = irec_local + 1
            acqui_simu(ievent)%data_traces(irec_local,:,:) = Gather(irec, :, :)
          endif
        enddo
      endif

      ! send gather to other MPI slices
      do irank = 1, NPROC-1

        if (myrank == 0) then !! then send

          ! count the receiver in slice irank
          nsta_irank = 0
          do irec = 1,  NSTA
            if (acqui_simu(ievent)%islice_selected_rec(irec) == irank) nsta_irank = nsta_irank + 1
          enddo

          ! if there is receiver in slice irank then MPI send data
          if (nsta_irank > 0) then
            !! data to send
            allocate(Gather_loc(nsta_irank,Nt,NDIM),stat=ier)
            if (ier /= 0) call exit_MPI_without_rank('error allocating array 403')
            irec_local = 0
            do irec = 1, NSTA
              if (acqui_simu(ievent)%islice_selected_rec(irec) == irank) then
                irec_local = irec_local + 1
                Gather_loc(irec_local, :, :) = Gather(irec, :, :) !! store data to send
              endif
            enddo
            if (DEBUG_MODE) then
              write(IIDD,*) 'myrank ', myrank , 'send to ', irank, ' :' , nsta_irank, Nt
              call flush_iunit(IIDD)
            endif

            ! MPI send rank gather
            tag    = 2001
            call MPI_SEND(Gather_loc, Nt*nsta_irank*NDIM, CUSTOM_MPI_TYPE, irank, tag, my_local_mpi_comm_world, ier)

            deallocate(Gather_loc)
          endif

        else !! then receive gather

          if (myrank == irank .and. acqui_simu(ievent)%nsta_slice > 0) then
            NSTA_LOC = acqui_simu(ievent)%nsta_slice
            Nt = acqui_simu(ievent)%Nt_data

            allocate(Gather_loc(NSTA_LOC,Nt,NDIM),stat=ier)
            if (ier /= 0) call exit_MPI_without_rank('error allocating array 404')
            Gather_loc(:,:,:) = 0._CUSTOM_REAL

            allocate(acqui_simu(ievent)%data_traces(NSTA_LOC,Nt,NDIM),stat=ier)
            if (ier /= 0) call exit_MPI_without_rank('error allocating array 405')
            acqui_simu(ievent)%data_traces(:,:,:) = 0._CUSTOM_REAL

            allocate(acqui_simu(ievent)%adjoint_sources(NDIM, NSTA_LOC, Nt),stat=ier)
            if (ier /= 0) call exit_MPI_without_rank('error allocating array 406')
            acqui_simu(ievent)%adjoint_sources(:,:,:) = 0._CUSTOM_REAL

            ! weights for local stations only
            allocate(acqui_simu(ievent)%weight_trace(NDIM, NSTA_LOC,1),stat=ier)
            if (ier /= 0) call exit_MPI_without_rank('error allocating array 407')
            acqui_simu(ievent)%weight_trace(:,:,:) = 1._CUSTOM_REAL

            if (VERBOSE_MODE .or. DEBUG_MODE) then
              allocate(acqui_simu(ievent)%synt_traces(NDIM, NSTA_LOC, Nt),stat=ier)
              if (ier /= 0) call exit_MPI_without_rank('error allocating array 408')
              acqui_simu(ievent)%synt_traces(:,:,:) = 0._CUSTOM_REAL
            endif

            if (DEBUG_MODE) then
              write(IIDD,*) 'myrank ',myrank,' wait for 0 :', NSTA_LOC,Nt
              call flush_iunit(IIDD)
            endif

            ! MPI receive rank gather
            tag   = MPI_ANY_TAG
            call MPI_RECV(Gather_loc,Nt*NSTA_LOC*NDIM,CUSTOM_MPI_TYPE, 0, tag, my_local_mpi_comm_world, status,  ier)

            !! store in acqui_simu
            acqui_simu(ievent)%data_traces(:,:,:) = Gather_loc(:,:,:)
            deallocate(Gather_loc)
          endif

        endif

      enddo

      if (myrank == 0) deallocate(Gather)
      call synchronize_all()

      ! gradient weighting
      if (inversion_param%is_src_weight_gradient) then
        ! we take the maximum value of gather
        ! broadcast gather maximum (determined on myrank==0) to all slices
        call MPI_BCAST(gather_maximum,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
        if (acqui_simu(ievent)%nsta_slice > 0) then
          acqui_simu(ievent)%weight_trace(:,:,:) = acqui_simu(ievent)%weight_trace(:,:,:) / gather_maximum
        endif
      endif

      !! set other parameters (in future work need to read any additional files)

      !! get band pass filter values if needed
      if (use_band_pass_filter) then
        acqui_simu(ievent)%Nfrq = NIFRQ
        acqui_simu(ievent)%band_pass_filter = use_band_pass_filter
        allocate(acqui_simu(ievent)%fl_event(acqui_simu(ievent)%Nfrq),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 409')
        allocate(acqui_simu(ievent)%fh_event(acqui_simu(ievent)%Nfrq),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 410')
        acqui_simu(ievent)%fl_event(:) = fl(:)
        acqui_simu(ievent)%fh_event(:) = fh(:)
        !! WARNING WARNING
        !! this is for telesismic case for now only one
        !! frequency is allowed (todo fix it)
        acqui_simu(ievent)%freqcy_to_invert(:,1,:) = fl(1)
        acqui_simu(ievent)%freqcy_to_invert(:,2,:) = fh(1)
      endif

    enddo

    call MPI_BCAST(nb_traces_tot,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    dummy_real = nb_traces_tot
    call sum_all_all_cr_for_simulatenous_runs(dummy_real,nb_traces_tot,1)

    if (myrank == 0) write(INVERSE_LOG_FILE,*) ' READING data passed '

  end subroutine read_data_gather

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------
! main read waveform data gather and bcast to MPI slice concerned
!----------------------------------------------------------------
  subroutine read_pif_data_gather(acqui_simu, inversion_param, myrank)

    use my_mpi             !! module from specfem
    include "precision.h"  !! from specfem

    integer,                                     intent(in)    :: myrank
    type(acqui),  dimension(:),                  intent(inout) :: acqui_simu
    type(inver),                                 intent(inout) :: inversion_param

    integer                                                    :: ievent, idim, NSTA, NSTA_LOC, Nt, irec, irec_local
    integer                                                    :: tag, ier, nsta_irank, irank
    real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable      :: Gather, Gather_loc
    real(kind=CUSTOM_REAL), dimension(:),     allocatable      :: weight_loc
    integer                                                    :: status(MPI_STATUS_SIZE)
    real(kind=CUSTOM_REAL)                                     :: dummy_real, W
    character(len=MAX_LEN_STRING)                              :: filename
    logical                                                    :: data_comp_inv, data_comp_read
    character(len=1)                                           :: data_type_inv, data_type_read
    character(len=3)                                           :: data_sys_inv, data_sys_read
    integer                                                    :: it1, it2, it3, it4
    double precision                                           :: lat0, lon0, azi0

    nb_traces_tot = 0
    W = 1.

    if (myrank == 0) write(INVERSE_LOG_FILE,'(/a17)') '... reading data '

    do ievent = 1, acqui_simu(1)%nevent_tot

      lat0 = acqui_simu(ievent)%Origin_chunk_lat
      lon0 = acqui_simu(ievent)%Origin_chunk_lon
      azi0 = acqui_simu(ievent)%Origin_chunk_azi

      if (myrank == 0) then

        NSTA = acqui_simu(ievent)%nsta_tot
        Nt = acqui_simu(ievent)%Nt_data

        allocate(Gather(NSTA,Nt,NDIM),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 411')
        Gather(:,:,:) = 0._CUSTOM_REAL

        !! Read pif gather file component by conponent
        ncomp_read = 0
        data_type_read = acqui_simu(ievent)%read_data_type

        do idim = 1,NDIM
          ! Check type of data
          data_sys_read  = acqui_simu(ievent)%read_data_sys(idim:idim)
          data_comp_read = acqui_simu(ievent)%read_data_comp(idim)

          write(filename,*) trim(acqui_simu(ievent)%event_rep),'/fsismo_', &
                            trim(adjustl(data_type_read)),trim(adjustl(data_sys_read)),'.bin'

          ! Read data
          if (data_comp_read) then
            call read_binary_data(filename, nsta, nt, gather(:,:,idim))
            nb_traces_tot = nb_traces_tot + nsta
            ncomp_read = ncomp_read + 1
          endif
        enddo

        !! store data gather in my slice if needed
        NSTA_LOC = acqui_simu(ievent)%nsta_slice
        allocate(acqui_simu(ievent)%data_traces(ndim,nsta_loc,nt),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 412')
        acqui_simu(ievent)%data_traces(:,:,:) = 0._CUSTOM_REAL

        allocate(acqui_simu(ievent)%adjoint_sources(ndim,nsta_loc,nt),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 413')
        acqui_simu(ievent)%adjoint_sources(:,:,:) = 0._CUSTOM_REAL

        ! weights for total number of stations
        allocate(acqui_simu(ievent)%weight_trace(ndim,NSTA,nt),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 414')
        acqui_simu(ievent)%weight_trace(:,:,:) = 1._CUSTOM_REAL

        !! manage data taper here
        ! first windowing
        if (acqui_simu(ievent)%is_time_pick) then
          allocate(weight_loc(nt),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 415')
          ! loops over all stations
          ! (time picks provided for all station records in teleseismic acquisition file)
          do irec = 1, NSTA
            weight_loc(:) = 1._CUSTOM_REAL
            it1 = int(floor(acqui_simu(ievent)%time_pick(irec) / acqui_simu(ievent)%dt_data ))
            it2 = int(floor((acqui_simu(ievent)%time_pick(irec) &
                            - acqui_simu(ievent)%time_before_pick) &
                            / acqui_simu(ievent)%dt_data))
            it3 = int(ceiling((acqui_simu(ievent)%time_pick(irec) &
                            + acqui_simu(ievent)%time_after_pick) &
                            / acqui_simu(ievent)%dt_data))
            it4 = int(ceiling((acqui_simu(ievent)%time_pick(irec) &
                            + acqui_simu(ievent)%time_after_pick &
                            + acqui_simu(ievent)%time_before_pick) /  acqui_simu(ievent)%dt_data))
            call taper_window_W(weight_loc,it1,it2,it3,it4,nt,W)
            do idim = 1,ndim
              acqui_simu(ievent)%weight_trace(idim,irec,:) = weight_loc(:)
            enddo
          enddo
        endif

        ! then gradient weighting
        if (inversion_param%is_src_weight_gradient) then ! we take the maximum value of gather
          acqui_simu(ievent)%weight_trace(:,:,:) = acqui_simu(ievent)%weight_trace(:,:,:) / maxval(abs(gather))
        endif

        ! manage inverted data (use weight_trace to select which parameter)
        ncomp_inv = 0
        data_type_inv = inversion_param%inverted_data_type
        if (data_type_inv /= data_type_read) then
          write(*,*)'ERROR: requested type of inverted data is different from observed data'
          write(*,*)'       integration of differentiation of observed not implemented yet'
          write(*,*)'NOW STOP'
          stop
        endif
        if (data_type_inv == 'd') inversion_param%get_synthetic_displacement = .true.
        if (data_type_inv == 'v') inversion_param%get_synthetic_velocity     = .true.
        if (data_type_inv == 'a') inversion_param%get_synthetic_acceleration = .true.
        if (data_type_inv == 'p') inversion_param%get_synthetic_pressure     = .true.

        ! mute non_inverted data
        do idim = 1,ndim
          ! Check type of data
          data_sys_inv  = inversion_param%inverted_data_sys(idim:idim)
          data_comp_inv = inversion_param%inverted_data_comp(idim)

          if (.not. data_comp_inv) then
            acqui_simu(ievent)%weight_trace(idim,:,:) = 0._CUSTOM_REAL
          endif
        enddo

        ! pass from data system to local mesh
        select case (acqui_simu(ievent)%read_data_sys)
        case ('xyz')
          ! nothing to do
        case('enz')
          ! Data are in standard coordinate system
          ! Data rotation required to pass in mesh system (zen -> xyz)
          call define_mesh_rotation_matrix(lat0,lon0,azi0)
          call rotate_comp_glob2mesh(gather(:,:,3), gather(:,:,2), gather(:,:,1), &
                                     acqui_simu(ievent)%read_station_position(1,:), &
                                     acqui_simu(ievent)%read_station_position(2,:), &
                                     nt, nsta, gather(:,:,1), gather(:,:,2), gather(:,:,3))
        case('rtz')
          ! Data are in the souce receiver coordinate system
          ! Data rotation required (baz-azi) (rtz -> zne)
          call rotate_ZRT_to_ZNE(gather(:,:,3), gather(:,:,1), gather(:,:,2), &
                                 gather(:,:,3), gather(:,:,2), gather(:,:,1), &
                                 nsta,nt,acqui_simu(ievent)%baz)
          ! Data rotation required to pass in mesh system (zen -> xyz)
          call define_mesh_rotation_matrix(lat0, lon0, azi0)
          call rotate_comp_glob2mesh(gather(:,:,3), gather(:,:,2), gather(:,:,1), &
                                     acqui_simu(ievent)%read_station_position(1,:), &
                                     acqui_simu(ievent)%read_station_position(2,:), &
                                     nt, nsta, gather(:,:,1), gather(:,:,2), gather(:,:,3))
        case('qtl')
          !! Data are in the ray coordinate system
          !! Data rotation required (baz-azi and incidence angle) (rtz -> zen)
          !call rotate_LQT_to_ZNE(vl,vq,vt,vz,vn,ve,nrec,nt,bazi,inci)
          !
          !! Data rotation required to pass in mesh system (zen -> xyz)
          !call define_mesh_rotation_matrix(lat0,lon0,azi0)
          !call rotate_comp_glob2mesh(vz2, vn, ve, stalat, stalon, nt, nsta, vx, vy, vz)
          write(*,*)'ERROR: qtl is not implemented yet'
          write(*,*)'NOW STOP'
          stop
        end select

        if (VERBOSE_MODE .or. DEBUG_MODE) then
          allocate(acqui_simu(ievent)%synt_traces(NDIM, NSTA_LOC, Nt),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 416')
          acqui_simu(ievent)%synt_traces(:,:,:) = 0._CUSTOM_REAL
        endif

        ! store (rotated) data traces
        irec_local = 0
        do irec = 1, NSTA
          if (acqui_simu(ievent)%islice_selected_rec(irec) == myrank) then
            irec_local = irec_local + 1
            acqui_simu(ievent)%data_traces(irec_local,:,:) = Gather(irec, :, :)
          endif
        enddo
      endif

      ! send data gather to other MPI slices
      do irank = 1, NPROC-1

        if (myrank == 0) then !! then send

          ! count the receiver in slice irank
          nsta_irank = 0
          do irec = 1,  NSTA
            if (acqui_simu(ievent)%islice_selected_rec(irec) == irank) nsta_irank = nsta_irank + 1
          enddo

          ! if there is receiver in slice irank then MPI send data
          if (nsta_irank > 0) then
            !! data to send
            allocate(Gather_loc(nsta_irank,Nt,NDIM),stat=ier)
            if (ier /= 0) call exit_MPI_without_rank('error allocating array 417')
            irec_local = 0
            do irec = 1, NSTA
              if (acqui_simu(ievent)%islice_selected_rec(irec) == irank) then
                irec_local = irec_local + 1
                Gather_loc(irec_local, :, :) = Gather(irec, :, :) !! store data to send
              endif
            enddo
            if (DEBUG_MODE) then
              write(IIDD,*) 'myrank ', myrank , 'send to ', irank, ' :' , nsta_irank, Nt
              call flush_iunit(IIDD)
            endif

            ! MPI send rank gather
            tag    = 2001
            call MPI_SEND(Gather_loc, Nt*nsta_irank*NDIM, CUSTOM_MPI_TYPE, irank, tag, my_local_mpi_comm_world, ier)

            deallocate(Gather_loc)
          endif

        else !! then receive gather

          if (myrank == irank .and. acqui_simu(ievent)%nsta_slice > 0) then
            NSTA = acqui_simu(ievent)%nsta_tot
            NSTA_LOC = acqui_simu(ievent)%nsta_slice
            Nt = acqui_simu(ievent)%Nt_data

            allocate(Gather_loc(NSTA_LOC,Nt,NDIM),stat=ier)
            if (ier /= 0) call exit_MPI_without_rank('error allocating array 418')
            Gather_loc(:,:,:) = 0._CUSTOM_REAL

            allocate(acqui_simu(ievent)%data_traces(NSTA_LOC,Nt,NDIM),stat=ier)
            if (ier /= 0) call exit_MPI_without_rank('error allocating array 419')
            acqui_simu(ievent)%data_traces(:,:,:) = 0._CUSTOM_REAL

            allocate(acqui_simu(ievent)%adjoint_sources(NDIM,NSTA_LOC,Nt),stat=ier)
            if (ier /= 0) call exit_MPI_without_rank('error allocating array 420')
            acqui_simu(ievent)%adjoint_sources(:,:,:) = 0._CUSTOM_REAL

            ! weights for total number of stations
            allocate(acqui_simu(ievent)%weight_trace(NDIM,nsta,nt),stat=ier)
            if (ier /= 0) call exit_MPI_without_rank('error allocating array 421')
            acqui_simu(ievent)%weight_trace(:,:,:) = 1._CUSTOM_REAL

            if (VERBOSE_MODE .or. DEBUG_MODE) then
              allocate(acqui_simu(ievent)%synt_traces(NDIM, NSTA_LOC, Nt),stat=ier)
              if (ier /= 0) call exit_MPI_without_rank('error allocating array 422')
              acqui_simu(ievent)%synt_traces(:,:,:) = 0._CUSTOM_REAL
            endif

            if (DEBUG_MODE) then
              write(IIDD,*) 'myrank ',myrank,' wait for 0 :', NSTA_LOC,Nt
              call flush_iunit(IIDD)
            endif

            ! MPI receive rank gather
            tag   = MPI_ANY_TAG
            call MPI_RECV(Gather_loc,Nt*NSTA_LOC*NDIM,CUSTOM_MPI_TYPE, 0, tag, my_local_mpi_comm_world, status,  ier)

            !! store in acqui_simu
            acqui_simu(ievent)%data_traces(:,:,:) = Gather_loc(:,:,:)

            deallocate(Gather_loc)
          endif

        endif

      enddo

      if (myrank == 0) deallocate(Gather)
      call synchronize_all()

      ! broadcast data weights (determined on myrank==0) to all slices
      NSTA = acqui_simu(ievent)%nsta_tot
      Nt = acqui_simu(ievent)%Nt_data
      call MPI_BCAST(acqui_simu(ievent)%weight_trace,NDIM * NSTA * Nt,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)

      !! set other parameters (in future work need to read any additional files)

      !! get band pass filter values if needed
      if (use_band_pass_filter) then
        acqui_simu(ievent)%Nfrq = NIFRQ
        acqui_simu(ievent)%band_pass_filter = use_band_pass_filter
        allocate(acqui_simu(ievent)%fl_event(acqui_simu(ievent)%Nfrq),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 423')
        allocate(acqui_simu(ievent)%fh_event(acqui_simu(ievent)%Nfrq),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 424')
        acqui_simu(ievent)%fl_event(:) = fl(:)
        acqui_simu(ievent)%fh_event(:) = fh(:)
        !! WARNING WARNING
        !! this is for telesismic case for now only one
        !! frequency is allowed (todo fix it)
        acqui_simu(ievent)%freqcy_to_invert(:,1,:) = fl(1)
        acqui_simu(ievent)%freqcy_to_invert(:,2,:) = fh(1)
      endif

    enddo

    call MPI_BCAST(nb_traces_tot,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    dummy_real = nb_traces_tot
    call sum_all_all_cr_for_simulatenous_runs(dummy_real,nb_traces_tot,1)

    if (myrank == 0) write(INVERSE_LOG_FILE,'(a25//)') '... reading data : passed'

  end subroutine read_pif_data_gather

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! WARNING TO BE CONSISTENT I USED NAME FILE TO WRITE TO CHOOSE BETWEEN DATA AND ADJOINT SOURCES
  subroutine write_pif_data_gather(ievent, acqui_simu, inversion_param, array_to_write, name_file_to_write, myrank)

    use my_mpi             !! module from specfem
    include "precision.h"  !! from specfem

    integer,                                   intent(in)    :: myrank, ievent
    real(kind=CUSTOM_REAL),  dimension(:,:,:), intent(in)    :: array_to_write
    type(acqui),             dimension(:),     intent(inout) :: acqui_simu
    character(len=MAX_LEN_STRING),             intent(in)    :: name_file_to_write
    type(inver),                               intent(in)    :: inversion_param
    ! local
    character(len=MAX_LEN_STRING)                            :: filename
    integer                                                  :: idim, NSTA, NSTA_LOC, Nt, irec, irec_local
    integer                                                  :: tag, ier, nsta_irank, irank, icomp
    real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable    :: Gather, Gather_loc
    integer,                dimension(:),     allocatable    :: irec_global
    integer                                                  :: status(MPI_STATUS_SIZE)
    double precision                                         :: lat0, lon0, azi0
    character(len=1)                                         :: data_type
    logical                                                  :: data_comp
    character(len=3)                                         :: data_sys

    lat0 = acqui_simu(ievent)%Origin_chunk_lat
    lon0 = acqui_simu(ievent)%Origin_chunk_lon
    azi0 = acqui_simu(ievent)%Origin_chunk_azi

    if (myrank == 0) then
      NSTA = acqui_simu(ievent)%nsta_tot
      Nt = acqui_simu(ievent)%Nt_data
      allocate(Gather(NSTA,Nt,NDIM),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 425')
    endif
    ! not sure if need this sync
    call synchronize_all()

    ! collects array data on main process
    do irank = 1, NPROC-1
      ! main process gathers all data from other processes
      if (myrank == 0) then
        ! count the receiver in slice irank
        nsta_irank = 0
        do irec = 1,  NSTA
          if (acqui_simu(ievent)%islice_selected_rec(irec) == irank) nsta_irank = nsta_irank + 1
        enddo

        if (nsta_irank > 0) then
          !! data to receive
          allocate(Gather_loc(nsta_irank,Nt,NDIM),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 426')
          allocate(irec_global(nsta_irank),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 427')
          irec_local = 0

          ! receives gather data
          tag   = MPI_ANY_TAG
          call MPI_RECV(Gather_loc, Nt*nsta_irank*NDIM, CUSTOM_MPI_TYPE, irank, tag, my_local_mpi_comm_world, status,  ier)
          call MPI_RECV(irec_global, nsta_irank, MPI_INTEGER, irank, tag, my_local_mpi_comm_world, status,  ier)

          do icomp=1,NDIM
            do irec_local = 1, nsta_irank
              Gather(irec_global(irec_local), :, icomp) = Gather_loc(irec_local, :, icomp)
            enddo
          enddo

          deallocate(Gather_loc)
          deallocate(irec_global)
        endif
      else
        ! other slice sends data to main process
        if (myrank == irank .and. acqui_simu(ievent)%nsta_slice > 0) then
          NSTA_LOC = acqui_simu(ievent)%nsta_slice
          Nt = acqui_simu(ievent)%Nt_data
          allocate(Gather_loc(NSTA_LOC,Nt,NDIM),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 428')
          allocate(irec_global(NSTA_LOC),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 429')

          do irec_local = 1, NSTA_LOC
            irec_global(irec_local) = acqui_simu(ievent)%number_receiver_global(irec_local)
            !! choose the rigth seismograms_*
            do icomp=1,NDIM
              Gather_loc(irec_local,:,icomp) = array_to_write(icomp,irec_local,:)
            enddo
          enddo
          ! MPI send dat to main process
          tag    = 2001
          call MPI_SEND(Gather_loc,  Nt*NSTA_LOC*NDIM, CUSTOM_MPI_TYPE, 0, tag, my_local_mpi_comm_world, ier)
          call MPI_SEND(irec_global, NSTA_LOC, CUSTOM_MPI_TYPE, 0, tag, my_local_mpi_comm_world, ier)

          deallocate(Gather_loc)
          deallocate(irec_global)
        endif
      endif

      ! not sure if need this sync
      call synchronize_all()
    enddo

    !!  write gather file
    if (myrank == 0) then
      ! main process still needs to gather its own array data
      do icomp=1,NDIM
        do irec_local = 1, acqui_simu(ievent)%nsta_slice
          !! choose the right seismograms_*
          Gather(acqui_simu(ievent)%number_receiver_global(irec_local),:,icomp) = array_to_write(icomp,irec_local,:)
        enddo
      enddo

      !! write only the asked component or pressure
      select case (trim(adjustl(name_file_to_write)))
      case('adjoint_source')
        ! we write adjoint sources in the inverted data system (write all))
        irec = 0
        data_type = inversion_param%inverted_data_type
        do idim = 1,ndim
          data_sys  = inversion_param%inverted_data_sys(idim:idim)
          data_comp = inversion_param%inverted_data_comp(idim)
          write(filename,*)trim(acqui_simu(ievent)%data_file_gather),'/adjoint_src_fsismo_', &
                trim(adjustl(data_type)),trim(adjustl(data_sys)),'.bin'
          call write_binary_data(filename,nsta,nt,gather(:,:,idim))
        enddo

      case('data')
        ! we write data the read data system (write all)
        irec = 0
        data_type = acqui_simu(ievent)%read_data_type

        ! perform rotations
        select case (acqui_simu(ievent)%read_data_sys)
        case ('xyz')
          ! nothing to do
        case('enz')
          ! Data are in standard coordinate system
          ! Data rotation required to pass in mesh system (zen -> xyz)
          call define_mesh_rotation_matrix(lat0,lon0,azi0)
          call rotate_comp_mesh2glob(gather(:,:,1), gather(:,:,2), gather(:,:,3), &
                                     acqui_simu(ievent)%read_station_position(1,:), &
                                     acqui_simu(ievent)%read_station_position(2,:), &
                                     nt, nsta, gather(:,:,3), gather(:,:,2), gather(:,:,1))
        case('rtz')
          ! Data are in the souce receiver coordinate system
          ! Data rotation required (baz-azi) (rtz -> zne)
          ! Data rotation required to pass in mesh system (zen -> xyz)
          call define_mesh_rotation_matrix(lat0,lon0,azi0)
          call rotate_comp_mesh2glob(gather(:,:,1), gather(:,:,2), gather(:,:,3), &
                                     acqui_simu(ievent)%read_station_position(1,:), &
                                     acqui_simu(ievent)%read_station_position(2,:), &
                                     nt, nsta, gather(:,:,3), gather(:,:,2), gather(:,:,1))
          call rotate_ZNE_to_ZRT(gather(:,:,3), gather(:,:,2), gather(:,:,1), &
                                 gather(:,:,3), gather(:,:,1), gather(:,:,2), &
                                 nsta,nt,acqui_simu(ievent)%baz)
        case('qtl')
          !! Data are in the ray coordinate system
          !! Data rotation required (baz-azi and incidence angle) (rtz -> zen)
          !call rotate_LQT_to_ZNE(vl,vq,vt,vz,vn,ve,nrec,nt,bazi,inci)
          !
          !! Data rotation required to pass in mesh system (zen -> xyz)
          !call define_mesh_rotation_matrix(lat0,lon0,azi0)
          !call rotate_comp_glob2mesh(vz2, vn, ve, stalat, stalon, nt, nsta, vx, vy, vz)
          write(*,*)'ERROR: qtl is not implemented yet'
          write(*,*)'NOW STOP'
          stop
        end select

        ! output data to file
        do idim = 1,ndim
          data_sys  = acqui_simu(ievent)%read_data_sys(idim:idim)
          data_comp = acqui_simu(ievent)%read_data_comp(idim)
          write(filename,*)trim(acqui_simu(ievent)%data_file_gather),'/output_fsismo_', &
                trim(adjustl(data_type)),trim(adjustl(data_sys)),'.bin'
          call write_binary_data(filename,nsta,nt,gather(:,:,idim))
        enddo

      case default
        ! write as it stands
        ! we write adjoint sources in the inverted data system (write all))
        irec = 0
        data_type = inversion_param%inverted_data_type
        do idim = 1,ndim
          data_sys  = inversion_param%inverted_data_sys(idim:idim)
          data_comp = inversion_param%inverted_data_comp(idim)
          if (len_trim(adjustl(name_file_to_write)) > 0) then
            filename = trim(name_file_to_write) // '_'
          else
            filename = 'unknown_fsismo_'
          endif
          write(filename,*) filename,idim,'.bin'
          filename = trim(acqui_simu(ievent)%event_rep) // '/' // trim(filename)
          call write_binary_data(filename,nsta,nt,gather(:,:,idim))
        enddo
      end select

      ! free temporary array
      deallocate(Gather)
    endif

    ! to be sure... there is a bug somewhere...
    call synchronize_all()

  end subroutine write_pif_data_gather

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------
! read acquisiton input file and store acqui_simu type
!----------------------------------------------------------------
  subroutine read_acqui_file(acqui_file, acqui_simu, myrank)

    use my_mpi             !! module from specfem
    include "precision.h"  !! from specfem

    character(len=MAX_LEN_STRING),           intent(in)    ::  acqui_file
    integer,                                 intent(in)    ::  myrank
    type(acqui),  dimension(:), allocatable, intent(inout) ::  acqui_simu
    ! locals
    character(len=MAX_LEN_STRING)                          :: line, keyw, line_to_read
    integer                                                :: ipos0, ipos1, ievent
    integer                                                :: ier

    ! debug output
    if (DEBUG_MODE) then
      if (myrank == 0) write(IIDD,*) 'read_acqui_file:'
      call synchronize_all()
      write(IIDD,*) '       MYGROUP  ', mygroup, '    MYRANK ', myrank
      if (myrank == 0) then
        write(IIDD,*)
        call flush_iunit(IIDD)
      endif
    endif

    ! initializes
    NEVENT = 0

    ! log output
    if (myrank == 0) then
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*) ' READING acquisition'
      write(INVERSE_LOG_FILE,*)
      call flush_iunit(INVERSE_LOG_FILE)
    endif

    ! only main reads acqui file
    if (myrank == 0) then

      !! 1/ read to count the number of events
      open(666,file=trim(acqui_file),iostat=ier)
      if (ier /= 0) then
        write(*,*) 'Error opening  ', trim(acqui_file), ' mygroup ', mygroup
        stop
      endif
      ! debug output
      if (DEBUG_MODE) then
        write(IIDD,*) ' opening  ', trim(acqui_file), ' mygroup ', mygroup
        call flush_iunit(IIDD)
      endif

      do

        read(666,'(a)',end=99) line
        if (DEBUG_MODE) write(IIDD,'(a)') trim(line)
        if (is_blank_line(line)) cycle                 !! no significant line
        if (INDEX(line,'event_name') > 0) NEVENT = NEVENT+1  !! new event

      enddo
99    close(666)

      write(INVERSE_LOG_FILE,*) ' ALLOCATE  acquisition structure for ', NEVENT, ' events '

      !! 2/ allocate and store type(acqui) acqui_simu
      if (NEVENT > 0) then
        allocate(acqui_simu(NEVENT),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 430')
      else
        write(*,*) 'ERROR NO EVENTS FOUND IN ACQUISITION FILE ',myrank, mygroup, trim(acqui_file)
        stop
      endif

      ! open event file
      open(666,file=trim(acqui_file))
      ievent = 0
      do    !! loop on all lines

        !! READ AND STORE ALL ITEM RELATED TO EVENT -------------------------------
        do
          read(666,'(a)',end=999) line
          if (is_blank_line(line)) cycle

          !! INDICES TO READ line -----------------------------------------------
          ipos0 = index(line,':')+1
          ipos1 = index(line,'#')-1
          if (ipos1 < 0 ) ipos1 = len_trim(line)

          !! STORE KEYWORD ITEM -------------------------------------------------
          keyw = trim(adjustl(line(1:ipos0-2)))
          if (len_trim(keyw) == 0) cycle

          !! DIFFERENT ITEM TO READ ---------------------------------------------
          select case (trim(keyw))

          case('event_name')
            ievent = ievent + 1
            acqui_simu(ievent)%event_name = trim(adjustl(line(ipos0:ipos1)))
            acqui_simu(ievent)%nevent_tot = NEVENT
            call store_default_acqui_values(acqui_simu, ievent)

          case('moment', 'force')
            acqui_simu(ievent)%source_file = trim(adjustl(line(ipos0:ipos1)))
            acqui_simu(ievent)%source_type = trim(adjustl(keyw))
            acqui_simu(ievent)%adjoint_source_type = 'L2_OIL_INDUSTRY'

          case('shot')
            acqui_simu(ievent)%source_type = trim(adjustl(keyw))
            read(line(ipos0:ipos1),*) acqui_simu(ievent)%xshot, &
                                      acqui_simu(ievent)%yshot, &
                                      acqui_simu(ievent)%zshot, &
                                      acqui_simu(ievent)%shot_ampl
            acqui_simu(ievent)%adjoint_source_type = 'L2_OIL_INDUSTRY'

          case('axisem', 'dsm', 'plane', 'fk')
            acqui_simu(ievent)%source_file = trim(adjustl(line(ipos0:ipos1)))
            acqui_simu(ievent)%source_type = trim(adjustl(keyw))
            acqui_simu(ievent)%adjoint_source_type = 'L2_FWI_TELESEISMIC'

          case('traction_dir')
            acqui_simu(ievent)%traction_dir = trim(adjustl(line(ipos0:ipos1)))

          case('station_file')
            acqui_simu(ievent)%station_file = trim(adjustl(line(ipos0:ipos1)))

          case('data_file')
            acqui_simu(ievent)%data_file_gather = trim(adjustl(line(ipos0:ipos1)))

          case ('component')
            !! initiallize components
            acqui_simu(ievent)%component(1) = '  '
            acqui_simu(ievent)%component(2) = '  '
            acqui_simu(ievent)%component(3) = '  '
            !need to add '00' in case of missing components (because Fortran cannot let the missing component to '  ')
            line_to_read = line(ipos0:ipos1)//' 00 00 00'

            read(line_to_read,*)  acqui_simu(ievent)%component(1), &
                                  acqui_simu(ievent)%component(2), &
                                  acqui_simu(ievent)%component(3)


          case('source_wavelet')
            acqui_simu(ievent)%source_wavelet_file = trim(adjustl(line(ipos0:ipos1)))
            acqui_simu(ievent)%external_source_wavelet = .true.

          case ('NSTEP')
            read(line(ipos0:ipos1),*) acqui_simu(ievent)%Nt_data

          case ('DT')
            read(line(ipos0:ipos1),*) acqui_simu(ievent)%dt_data

          case default
            write(*,*) 'ERROR INVALID KEY WORD : ', trim(keyw), ' in file ', trim(acqui_file)
            stop

          end select

        enddo

      enddo

999   close(666)

    endif ! myrank == 0

    if (myrank == 0) then
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*) ' READING acquisition passed '
      call flush_iunit(INVERSE_LOG_FILE)
    endif

    ! main broadcasts read values
    call MPI_BCAST(NEVENT,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
    if (myrank > 0) then
      allocate(acqui_simu(NEVENT),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 432')
    endif
    do ievent = 1, NEVENT
       acqui_simu(ievent)%nevent_tot = NEVENT
       call MPI_BCAST(acqui_simu(ievent)%traction_dir,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(ievent)%data_file_gather,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(ievent)%source_type,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(ievent)%source_file,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(ievent)%source_wavelet_file, MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(ievent)%adjoint_source_type,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(ievent)%event_name,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(ievent)%component,6,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(ievent)%Nt_data,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(ievent)%dt_data,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(ievent)%xshot,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(ievent)%yshot,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(ievent)%zshot,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
       call MPI_BCAST(acqui_simu(ievent)%shot_ampl,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    enddo

    ! log output
    if (myrank == 0) then
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*) ' Nb tot events ', acqui_simu(1)%nevent_tot
      write(INVERSE_LOG_FILE,*)
      call flush_iunit(INVERSE_LOG_FILE)
    endif

  end subroutine read_acqui_file
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------
! read inversion input file and store acqui_simu type
!----------------------------------------------------------------
  subroutine read_inver_file(inver_file, inversion_param, myrank)

    use my_mpi             !! module from specfem
    include "precision.h"  !! from specfem

    character(len=MAX_LEN_STRING), intent(in)    ::  inver_file
    integer,                       intent(in)    ::  myrank
    type(inver),                   intent(inout) ::  inversion_param
    ! locals
    character(len=MAX_LEN_STRING)                :: line, keyw
    integer                                      :: ipos0,ipos1,ier

    ! only main reads inver_file
    if (myrank == 0) then

      open(666, file=trim(inver_file))
      do
        read(666,'(a)',end=99) line
        if (is_blank_line(line)) cycle

        !! INDICES TO READ line -----------------------------------------------
        ipos0 = index(line,':')+1
        ipos1 = index(line,'#')-1
        if (ipos1 < 0 ) ipos1 = len_trim(line)

        !! STORE KEYWORD ITEM -------------------------------------------------
        keyw = trim(adjustl(line(1:ipos0-2)))

        !! DIFFERENT ITEM TO READ ---------------------------------------------
        select case (trim(keyw))

        case ('Niter')
          read(line(ipos0:ipos1),*)  inversion_param%Niter

        case('Niter_wolfe')
          read(line(ipos0:ipos1),*)  inversion_param%Niter_wolfe

        case('m1_wolfe')
          read(line(ipos0:ipos1),*)  inversion_param%m1

        case('m2_wolfe')
          read(line(ipos0:ipos1),*)  inversion_param%m2

        case('max_history_bfgs')
          read(line(ipos0:ipos1),*)  inversion_param%max_history_bfgs

        case('max_relative_pert')
          read(line(ipos0:ipos1),*)  inversion_param%max_relative_pert

        case('param_family')
          read(line(ipos0:ipos1),*) inversion_param%parameter_family_name

        case('nb_inver')
          read(line(ipos0:ipos1),*) inversion_param%NinvPar

        case('param_to_inv')
          read(line(ipos0:ipos1),*) inversion_param%param_inv_name(1:inversion_param%NinvPar)

        case('use_frequency_band_pass_filter')
          inversion_param%use_band_pass_filter = .true.
          read(line(ipos0:ipos1),*) inversion_param%Nifrq
          allocate(fl(inversion_param%Nifrq),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 433')
          allocate(fh(inversion_param%Nifrq),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 434')

        case('fl')
          read(line(ipos0:ipos1),*) fl(1:inversion_param%Nifrq)

        case('fh')
          read(line(ipos0:ipos1),*) fh(1:inversion_param%Nifrq)

        case('input_sem_model')
          read(line(ipos0:ipos1),*)  inversion_param%input_sem_model

        case('input_sem_prior')
          read(line(ipos0:ipos1),*)  inversion_param%input_sem_prior

        case('output_model')
          read(line(ipos0:ipos1),*)  inversion_param%output_model

        case('input_fd_model')
          read(line(ipos0:ipos1),*)  inversion_param%input_fd_model

        case ('taper')
          inversion_param%use_taper = .true.
          read(line(ipos0:ipos1),*) inversion_param%xmin_taper, inversion_param%xmax_taper, &
                                    inversion_param%ymin_taper, inversion_param%ymax_taper, &
                                    inversion_param%zmin_taper, inversion_param%zmax_taper

        case('shin_precond')
          read(line(ipos0:ipos1),*) inversion_param%shin_precond

        case('energy_precond')
          read(line(ipos0:ipos1),*) inversion_param%energy_precond

        case('z2_precond')
          read(line(ipos0:ipos1),*) inversion_param%z2_precond

        case('z_precond')
          inversion_param%z_precond = .true.
          read(line(ipos0:ipos1),*) inversion_param%aPrc, &
                                    inversion_param%zPrc1, &
                                    inversion_param%zPrc2

        case('relat_grad')
          read(line(ipos0:ipos1),*) inversion_param%relat_grad

        case('relat_cost')
          read(line(ipos0:ipos1),*) inversion_param%relat_cost

        case('use_tk_regularization_weight')
          read(line(ipos0:ipos1),*) inversion_param%weight_Tikhonov

        case('use_tk_sem_regularization')
          inversion_param%use_regularization_SEM_Tikhonov = .true.
          allocate(inversion_param%smooth_weight(inversion_param%NinvPar),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 435')
          read(line(ipos0:ipos1),*) inversion_param%smooth_weight(1:inversion_param%NinvPar)

        case('use_tk_sem_damping')
          inversion_param%use_damping_SEM_Tikhonov = .true.
          allocate(inversion_param%damp_weight(inversion_param%NinvPar),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 436')
          read(line(ipos0:ipos1),*) inversion_param%damp_weight(1:inversion_param%NinvPar)
          !! we have read standard deviation for model
          !! then need to change
          inversion_param%damp_weight(:) = 1._CUSTOM_REAL / inversion_param%damp_weight(:)**2

        case('use_tk_sem_variable_damping_sources')
          inversion_param%use_variable_SEM_damping_sources = .true.
          read(line(ipos0:ipos1),*) inversion_param%min_damp_src, inversion_param%max_damp_src, &
                                    inversion_param%distance_from_source
          ! makes sure SEM damping turned on
          inversion_param%use_damping_SEM_Tikhonov = .true.
          if (.not. allocated(inversion_param%damp_weight)) then
            allocate(inversion_param%damp_weight(inversion_param%NinvPar),stat=ier)
            if (ier /= 0) call exit_MPI_without_rank('error allocating array 436')
            inversion_param%damp_weight(:) = 1._CUSTOM_REAL
          endif

        case('use_tk_sem_variable_damping_stations')
          inversion_param%use_variable_SEM_damping_stations = .true.
          read(line(ipos0:ipos1),*) inversion_param%min_damp_sta, inversion_param%max_damp_sta, &
                                    inversion_param%distance_from_station
          ! makes sure SEM damping turned on
          inversion_param%use_damping_SEM_Tikhonov = .true.
          if (.not. allocated(inversion_param%damp_weight)) then
            allocate(inversion_param%damp_weight(inversion_param%NinvPar),stat=ier)
            if (ier /= 0) call exit_MPI_without_rank('error allocating array 436')
            inversion_param%damp_weight(:) = 1._CUSTOM_REAL
          endif

        case('use_tk_fd_regularization')
          read(line(ipos0:ipos1),*) inversion_param%use_regularization_FD_Tikhonov

        case('prior_data_std')
          read(line(ipos0:ipos1),*) inversion_param%prior_data_std

        case('data_to_invert_type')
          read(line(ipos0:ipos1),*) inversion_param%inverted_data_type

        case('data_to_invert_system')
          read(line(ipos0:ipos1),*) inversion_param%inverted_data_sys

        case('data_to_invert_is_component')
          read(line(ipos0:ipos1),*) inversion_param%inverted_data_comp

        case('apply_src_weighting_to_gradient')
          read(line(ipos0:ipos1),*) inversion_param%is_src_weight_gradient

        case('convolve_synth_with_wavelet')
          read(line(ipos0:ipos1),*) inversion_param%convolution_by_wavelet

        case('dump_model_at_each_iteration')
          read(line(ipos0:ipos1),*) inversion_param%dump_model_at_each_iteration

        case('dump_gradient_at_each_iteration')
          read(line(ipos0:ipos1),*) inversion_param%dump_gradient_at_each_iteration

        case('dump_descent_direction_at_each_iteration')
          read(line(ipos0:ipos1),*) inversion_param%dump_descent_direction_at_each_iteration

        case('dump_image_at_each_iteration')
          read(line(ipos0:ipos1),*) inversion_param%dump_image_at_each_iteration

        case default
          write(*,*) 'ERROR KEY WORD NOT MATCH : ', trim(keyw), ' in file ', trim(inver_file)
          exit

        end select

      enddo

99    close(666)

      ! log output
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*) ' READ  ', trim(inver_file)
      write(INVERSE_LOG_FILE,*)
      call flush_iunit(INVERSE_LOG_FILE)
    endif

    ! main broadcasts read values
    ! stopping criteria
    call MPI_BCAST(inversion_param%Niter,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%Niter_wolfe,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%m1,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%m2,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%relat_grad,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%relat_cost,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    ! L-BFGS
    call MPI_BCAST(inversion_param%max_history_bfgs,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%max_relative_pert,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    ! inversion parameters
    call MPI_BCAST(inversion_param%NinvPar,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%parameter_family_name,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%param_inv_name, 50*MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
    ! band pass filter
    call MPI_BCAST(inversion_param%use_band_pass_filter,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    if (inversion_param%use_band_pass_filter) then
      call MPI_BCAST(inversion_param%Nifrq, 1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
      if (myrank > 0) then
        allocate(fl(inversion_param%Nifrq), fh(inversion_param%Nifrq),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 437')
      endif
      call MPI_BCAST(fl,inversion_param%Nifrq,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
      call MPI_BCAST(fh,inversion_param%Nifrq,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    endif
    ! external models
    call MPI_BCAST(inversion_param%output_model,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%input_fd_model,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%input_sem_model,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%input_sem_prior,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    ! preconditioning
    call MPI_BCAST(inversion_param%use_taper,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%xmin_taper,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%xmax_taper,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%ymin_taper,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%ymax_taper,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%zmin_taper,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%zmax_taper,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%shin_precond,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%energy_precond,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%z2_precond,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%z_precond, 1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    ! regularization
    call MPI_BCAST(inversion_param%use_regularization_FD_Tikhonov,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%use_regularization_SEM_Tikhonov,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%aPrc,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%zPrc1,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%zPrc2,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%relat_cost,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%weight_Tikhonov,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%use_damping_SEM_Tikhonov,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%use_variable_SEM_damping_sources,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%min_damp_src,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%max_damp_src,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%distance_from_source,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%use_variable_SEM_damping_stations,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%min_damp_sta,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%max_damp_sta,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%distance_from_station,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    if (myrank > 0 .and. inversion_param%use_regularization_SEM_Tikhonov ) then
      allocate(inversion_param%smooth_weight(inversion_param%NinvPar),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 438')
    endif
    if (inversion_param%use_regularization_SEM_Tikhonov) &
      call MPI_BCAST(inversion_param%smooth_weight(1),inversion_param%NinvPar,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    if (myrank > 0 .and. inversion_param%use_damping_SEM_Tikhonov ) then
      allocate(inversion_param%damp_weight(inversion_param%NinvPar),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 439')
    endif
    if (inversion_param%use_damping_SEM_Tikhonov) &
      call MPI_BCAST(inversion_param%damp_weight(1),inversion_param%NinvPar,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    ! data processing
    call MPI_BCAST(inversion_param%prior_data_std,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
    !call MPI_BCAST(inversion_param%component,6,mpi_character,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%inverted_data_comp,3,mpi_logical,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%inverted_data_sys,3,mpi_character,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%inverted_data_type,1,mpi_character,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%is_src_weight_gradient,1,mpi_logical,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%convolution_by_wavelet,1,mpi_logical,0,my_local_mpi_comm_world,ier)
    ! visualization
    call MPI_BCAST(inversion_param%dump_model_at_each_iteration, 1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%dump_gradient_at_each_iteration, 1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%dump_descent_direction_at_each_iteration, 1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
    call MPI_BCAST(inversion_param%dump_image_at_each_iteration, 1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)

    ! debugging overrides dump flags
    if (DEBUG_MODE) then
      inversion_param%dump_model_at_each_iteration = .true.
      inversion_param%dump_gradient_at_each_iteration = .true.
      inversion_param%dump_descent_direction_at_each_iteration = .true.
      inversion_param%dump_image_at_each_iteration = .true.
    endif

    ! safety check
    ! need at least 1 frequency band
    if (inversion_param%Nifrq == 0) then
      write(*,*) 'Error: Invalid number of frequency bands, must be at least 1'
      stop
    endif

    ! Wolfe condition constants should be 0 < m1 < m2 < 1
    ! (see Nocedal & Wright, Numerical Optimization, 2006. chapter 3, page 34)
    if (inversion_param%m1 <= 0.0 .or. inversion_param%m1 >= 1.0) then
      write(*,*) 'Error: Invalid Wolfe constant m1, must be 0 < m1 < 1'
      stop
    endif
    if (inversion_param%m2 <= 0.0 .or. inversion_param%m2 >= 1.0) then
      write(*,*) 'Error: Invalid Wolfe constant m2, must be 0 < m2 < 1'
      stop
    endif
    if (inversion_param%m1 >= inversion_param%m2) then
      write(*,*) 'Error: Invalid Wolfe constants m1 & m2, must be 0 < m1 < m2 < 1'
      stop
    endif

    if (inversion_param%use_regularization_FD_Tikhonov) then
      write(*,*) " ABORT :: BROKEN OPTION : FD Tikhonov - not fully implemented and tested yet"
      stop
    endif

    !! set private values
    use_band_pass_filter = inversion_param%use_band_pass_filter
    NIFRQ = inversion_param%Nifrq

  end subroutine read_inver_file
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------
! define default values
!----------------------------------------------------------------
  subroutine store_default_acqui_values(acqui_simu, ievent)

    type(acqui), dimension(:), intent(inout)  :: acqui_simu
    integer,                   intent(in)     :: ievent

    acqui_simu(ievent)%source_file         = 'none'
    acqui_simu(ievent)%traction_dir        = 'none'
    acqui_simu(ievent)%data_file_gather    = 'none'
    acqui_simu(ievent)%source_type         = 'none'
    acqui_simu(ievent)%station_file        = 'none'
    acqui_simu(ievent)%adjoint_source_type = 'none'

  end subroutine store_default_acqui_values

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------
! true if blank line or commented line (ie that begins with #)
!----------------------------------------------------------------
  logical function is_blank_line(line)

    character(len=MAX_LEN_STRING), intent(in) :: line
    is_blank_line = .false.
    if (len(trim(adjustl(line))) == 0) is_blank_line = .true.
    if (INDEX(trim(adjustl(line)),'#') == 1) is_blank_line = .true.

  end function is_blank_line
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------
! true if blank line or commented line (ie that begins with #)
!----------------------------------------------------------------
  subroutine remove_blank_in_line(line, line_without_blank)

    character(len=MAX_LEN_STRING), intent(in)    :: line
    character(len=MAX_LEN_STRING), intent(inout) :: line_without_blank
    integer                                      :: n, i, k

    n = len_trim(adjustl(line))
    k = 0
    do i = 1, n
       if (line(i:i+1) == ' ') then
          ! it's blanc character then do nothing
       else
          k = k+1
          line_without_blank(k:k+1)=line(i:i+1)
       endif
    enddo

  end subroutine remove_blank_in_line
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------------------------------------------------
! read station file in specfem format
!------------------------------------------------------------
  subroutine get_stations(acqui_simu)

    type(acqui), dimension(:), intent(inout)               :: acqui_simu

    integer                                                :: ievent, irec, nsta, nrec_loc, ier
    character(len=MAX_LEN_STRING)                          :: rec_filename,filtered_rec_filename

    ! log output
    if (myrank == 0) then
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*) ' READING and LOCATING stations '
      call flush_iunit(INVERSE_LOG_FILE)
    endif

    ! loop on all events
    do ievent = 1, NEVENT
      ! only slice 0 will read the STATIONS files
      if (myrank == 0) then
        rec_filename          = trim(adjustl(acqui_simu(ievent)%station_file))
        filtered_rec_filename = rec_filename(1:len_trim(rec_filename))//'_FILTERED'
      else
        rec_filename          = 'dummy_string'
        filtered_rec_filename = 'dummy_string'
      endif

      ! output_solver.txt
      if (myrank == 0) then
        write(IMAIN,*) '**********************************'
        write(IMAIN,*) 'Inverse problem for model - event: ',ievent
        write(IMAIN,*) '**********************************'
        write(IMAIN,*)
        write(IMAIN,*) 'getting stations from: ',trim(rec_filename)
        write(IMAIN,*)
        call flush_IMAIN()
      endif

      call station_filter(rec_filename,filtered_rec_filename,nsta)
      acqui_simu(ievent)%nsta_tot = nsta

      ! log output
      if (myrank == 0) then
        write(INVERSE_LOG_FILE,*) '   event : ',ievent,' has stations : ',nsta
        call flush_iunit(INVERSE_LOG_FILE)
      endif

      allocate(acqui_simu(ievent)%station_name(nsta),acqui_simu(ievent)%network_name(nsta),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 440')
      allocate(acqui_simu(ievent)%position_station(3,nsta),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 441')
      allocate(acqui_simu(ievent)%xi_rec(nsta),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 442')
      allocate(acqui_simu(ievent)%eta_rec(nsta),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 443')
      allocate(acqui_simu(ievent)%gamma_rec(nsta),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 444')
      allocate(acqui_simu(ievent)%islice_selected_rec(nsta),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 445')
      allocate(acqui_simu(ievent)%ispec_selected_rec(nsta),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 446')
      allocate(acqui_simu(ievent)%number_receiver_global(nsta),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 447')
      allocate(acqui_simu(ievent)%nu_rec(NDIM,NDIM,nsta),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 448')

      nrec = nsta

      ! re-allocate memory for receiver arrays, i.e. stations given in STATIONS file
      if (allocated(islice_selected_rec)) then
        deallocate(islice_selected_rec,ispec_selected_rec, &
                   xi_receiver,eta_receiver,gamma_receiver,nu_rec)
        deallocate(station_name,network_name, &
                   stlat,stlon,stele,stbur)
      endif
      ! allocation with correct nrec sizes
      allocate(islice_selected_rec(nrec), &
              ispec_selected_rec(nrec), &
              xi_receiver(nrec), &
              eta_receiver(nrec), &
              gamma_receiver(nrec), &
              nu_rec(NDIM,NDIM,nrec), stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2067')
      if (ier /= 0) stop 'error allocating arrays for receivers'
      ! initializes arrays
      islice_selected_rec(:) = -1
      ispec_selected_rec(:) = 0
      xi_receiver(:) = 0.d0; eta_receiver(:) = 0.d0; gamma_receiver(:) = 0.d0
      nu_rec(:,:,:) = 0.0d0

      allocate(station_name(nrec), &
               network_name(nrec), &
               stlat(nrec), &
               stlon(nrec), &
               stele(nrec), &
               stbur(nrec),stat=ier)
       if (ier /= 0 ) call exit_MPI(myrank,'Error allocating receiver arrays')
       stlat(:) = 0.d0; stlon(:) = 0.d0; stele(:) = 0.d0; stbur(:) = 0.d0
       station_name(:) = ""; network_name(:) = ""

      ! reads STATIONS_FILTERED file, locates receivers in the mesh and compute Lagrange interpolators
      call locate_receivers(filtered_rec_filename,1.0d0,1.0d0)

      acqui_simu(ievent)%islice_selected_rec(:) = islice_selected_rec(:)
      acqui_simu(ievent)%ispec_selected_rec(:) = ispec_selected_rec(:)
      acqui_simu(ievent)%xi_rec(:) = xi_receiver(:)
      acqui_simu(ievent)%eta_rec(:) = eta_receiver(:)
      acqui_simu(ievent)%gamma_rec(:) = gamma_receiver(:)
      acqui_simu(ievent)%station_name(:) = station_name(:)
      acqui_simu(ievent)%network_name(:) = network_name(:)
      acqui_simu(ievent)%nu_rec(:,:,:) = nu_rec(:,:,:)
      ! stores positions
      acqui_simu(ievent)%position_station(1,:) = x_target_station(:)
      acqui_simu(ievent)%position_station(2,:) = y_target_station(:)
      acqui_simu(ievent)%position_station(3,:) = z_target_station(:)

      ! local to global receiver numbering
      nrec_loc = 0
      do irec = 1, nsta
        if (myrank == acqui_simu(ievent)%islice_selected_rec(irec)) then
          nrec_loc = nrec_loc+1
          acqui_simu(ievent)%number_receiver_global(nrec_loc) = irec
        endif
      enddo

      acqui_simu(ievent)%nsta_slice = nrec_loc
      if (acqui_simu(ievent)%nsta_slice > 0) then
        allocate(acqui_simu(ievent)%hxi    (NGLLX,nrec_loc),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 449')
        allocate(acqui_simu(ievent)%heta   (NGLLY,nrec_loc),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 450')
        allocate(acqui_simu(ievent)%hgamma (NGLLZ,nrec_loc),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 451')
        allocate(acqui_simu(ievent)%hpxi   (NGLLX,nrec_loc),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 452')
        allocate(acqui_simu(ievent)%hpeta  (NGLLY,nrec_loc),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 453')
        allocate(acqui_simu(ievent)%hpgamma(NGLLZ,nrec_loc),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 454')
        allocate(acqui_simu(ievent)%freqcy_to_invert(NDIM,2,nrec_loc),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 455')
      else
        allocate(acqui_simu(ievent)%hxi    (1,1),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 456')
        allocate(acqui_simu(ievent)%heta   (1,1),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 457')
        allocate(acqui_simu(ievent)%hgamma (1,1),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 458')
        allocate(acqui_simu(ievent)%hpxi   (1,1),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 459')
        allocate(acqui_simu(ievent)%hpeta  (1,1),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 460')
        allocate(acqui_simu(ievent)%hpgamma(1,1),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 461')
        allocate(acqui_simu(ievent)%freqcy_to_invert(NDIM,2,1),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 462')
      endif

      nrec_loc = 0
      do irec = 1, acqui_simu(ievent)%nsta_tot
        if (myrank == acqui_simu(ievent)%islice_selected_rec(irec)) then
          nrec_loc = nrec_loc + 1
          ! compute Lagrange polynomials at the receiver location
          call lagrange_any(acqui_simu(ievent)%xi_rec(irec),NGLLX,xigll, &
                            acqui_simu(ievent)%hxi(1,nrec_loc),acqui_simu(ievent)%hpxi(1,nrec_loc))
          call lagrange_any(acqui_simu(ievent)%eta_rec(irec),NGLLY,yigll, &
                            acqui_simu(ievent)%heta(1,nrec_loc),acqui_simu(ievent)%hpeta(1,nrec_loc))
          call lagrange_any(acqui_simu(ievent)%gamma_rec(irec),NGLLZ,zigll, &
                            acqui_simu(ievent)%hgamma(1,nrec_loc),acqui_simu(ievent)%hpgamma(1,nrec_loc))

        endif
      enddo
    enddo

    if (myrank == 0) then
      write(INVERSE_LOG_FILE,*) ' READING and LOCATING stations passed '
      write(INVERSE_LOG_FILE,*)
      call flush_iunit(INVERSE_LOG_FILE)
    endif

  end subroutine get_stations

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------------------------------------------------------------------------
! read source parameter file and compute (xi,eta,gamma) local position for each source
!-------------------------------------------------------------------------------------
  subroutine get_point_source(acqui_simu)

    use my_mpi             !! module from specfem
    include "precision.h"  !! from specfem

    type(acqui), dimension(:),              intent(inout)     :: acqui_simu
    ! locals
    character(len=MAX_STRING_LEN)                             :: filename
    integer                                                   :: ievent,isrc,ier,nsrc_loc
    double precision                                          :: min_tshift, elevation
    double precision, dimension(:), allocatable               :: x_target_source,y_target_source,z_target_source
    double precision, dimension(:), allocatable               :: Mxx,Myy,Mzz,Mxy,Mxz,Myz,lat,long,depth
    double precision, dimension(:,:), allocatable             :: moment_tensor
    double precision, dimension(:), allocatable               :: factor_force_source,Fx,Fy,Fz
    double precision, dimension(:), allocatable               :: xi_source,eta_source,gamma_source
    double precision, dimension(:,:,:), allocatable           :: nu_source
    real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: interparray
    double precision                                          :: final_distance_source,final_distance_squared
    integer                                                   :: idomain

    ! location search
    real(kind=CUSTOM_REAL) :: distance_min_glob,distance_max_glob
    real(kind=CUSTOM_REAL) :: elemsize_min_glob,elemsize_max_glob
    real(kind=CUSTOM_REAL) :: x_min_glob,x_max_glob
    real(kind=CUSTOM_REAL) :: y_min_glob,y_max_glob
    real(kind=CUSTOM_REAL) :: z_min_glob,z_max_glob

    real(kind=CUSTOM_REAL) :: dt_dummy
    integer :: it

    ! date
    integer :: yr,jda,ho,mi,mo,da
    double precision :: sec

    if (myrank == 0) then
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*) ' READING and LOCATING sources  '
      call flush_iunit(INVERSE_LOG_FILE)
    endif

    ! compute typical size of elements
    ! gets mesh dimensions
    call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                              x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
                              elemsize_min_glob,elemsize_max_glob, &
                              distance_min_glob,distance_max_glob)

    do ievent = 1,acqui_simu(1)%nevent_tot

      ! 1/ Get the number of sources
      select case (trim(acqui_simu(ievent)%source_type))

        case('moment','force')

          ! only slice 0 will read the sources files
          if (myrank == 0) then
            filename = trim(adjustl(acqui_simu(ievent)%source_file))
            call count_number_of_sources(NSOURCES,filename)
          endif
          ! NSOURCES has been updated in count_number_of_sources for slice 0, thus
          ! we broadcast it to other slices
          call MPI_BCAST(NSOURCES,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)

         case('shot')
            !! in case of shot we assume that we have only one point source
            NSOURCES = 1
            !! and we manatoru use external stf
            USE_EXTERNAL_SOURCE_FILE = .true.
            USE_RICKER_TIME_FUNCTION = .false.
         case default
          !! define here the way the number of sources is obtained
          write(*,*) 'Error: source type ',trim(acqui_simu(ievent)%source_type),' not implemented yet'
          stop
      end select

      ! output_solver.txt
      if (myrank == 0) then
        write(IMAIN,*) '**********************************'
        write(IMAIN,*) 'Inverse problem for model - event: ',ievent
        write(IMAIN,*) '**********************************'
        write(IMAIN,*)
        write(IMAIN,*) 'source type               : ',trim(acqui_simu(ievent)%source_type)
        if (trim(acqui_simu(ievent)%source_type) == 'moment' .or. trim(acqui_simu(ievent)%source_type) == 'force') &
          write(IMAIN,*) 'getting point sources from: ',trim(filename)
        write(IMAIN,*) 'number of point sources   : ',NSOURCES
        write(IMAIN,*)
        call flush_IMAIN()
      endif

      ! log output
      if (myrank == 0) then
        write(INVERSE_LOG_FILE,*) '   event : ',ievent,' has sources : ',NSOURCES
        call flush_iunit(INVERSE_LOG_FILE)
      endif

      ! 2/ Allocate the acqui_simu structure
      acqui_simu(ievent)%nsources_tot = NSOURCES

      allocate(acqui_simu(ievent)%islice_selected_source(NSOURCES),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 463')
      allocate(acqui_simu(ievent)%ispec_selected_source(NSOURCES),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 464')
      allocate(acqui_simu(ievent)%tshift(NSOURCES),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 465')
      allocate(acqui_simu(ievent)%hdur(NSOURCES),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 466')
      allocate(acqui_simu(ievent)%hdur_Gaussian(NSOURCES),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 467')
      allocate(acqui_simu(ievent)%sourcearrays(NSOURCES,NDIM,NGLLX,NGLLY,NGLLZ),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 468')
      allocate(acqui_simu(ievent)%Xs(NSOURCES),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 469')
      allocate(acqui_simu(ievent)%Ys(NSOURCES),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 470')
      allocate(acqui_simu(ievent)%Zs(NSOURCES), stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 471')
      if (ier /= 0) stop 'error allocating arrays for sources'

      ! re-allocates specfem source arrays
      if (allocated(islice_selected_source)) &
        deallocate(islice_selected_source,ispec_selected_source,tshift_src,hdur,hdur_Gaussian)
      ! allocates specfem source arrays with new size
      allocate(islice_selected_source(NSOURCES), &
               ispec_selected_source(NSOURCES), &
               tshift_src(NSOURCES), &
               hdur(NSOURCES), &
               hdur_Gaussian(NSOURCES),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2048')
      if (ier /= 0) stop 'error allocating arrays for sources'
      islice_selected_source(:) = -1
      ispec_selected_source(:) = 0
      tshift_src(:) = 0.d0; hdur(:) = 0.d0; hdur_Gaussian(:) = 0.d0

      ! allocates source arrays with new size
      allocate(Mxx(NSOURCES), &
               Myy(NSOURCES), &
               Mzz(NSOURCES), &
               Mxy(NSOURCES), &
               Mxz(NSOURCES), &
               Myz(NSOURCES), &
               xi_source(NSOURCES), &
               eta_source(NSOURCES), &
               gamma_source(NSOURCES), &
               nu_source(NDIM,NDIM,NSOURCES),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2049')
      if (ier /= 0) stop 'error allocating arrays for sources'
      ! initializes arrays
      Mxx(:) = 0.d0; Myy(:) = 0.d0; Mzz(:) = 0.d0
      Mxy(:) = 0.d0; Mxz(:) = 0.d0; Myz(:) = 0.d0
      xi_source(:) = 0.d0; eta_source(:) = 0.d0; gamma_source(:) = 0.d0
      nu_source(:,:,:) = 0.d0

      ! positions
      allocate(x_target_source(NSOURCES), &
               y_target_source(NSOURCES), &
               z_target_source(NSOURCES),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 480')
      if (ier /= 0) stop 'error allocating utm source arrays'
      x_target_source(:) = 0.d0; y_target_source(:) = 0.d0; z_target_source(:) = 0.d0

      ! point source arrays
      allocate(factor_force_source(NSOURCES),Fx(NSOURCES),Fy(NSOURCES),Fz(NSOURCES),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 485')
      if (ier /= 0) stop 'error allocating arrays for force point sources'
      factor_force_source(:) = 0.d0; Fx(:) = 0.d0; Fy(:) = 0.d0; Fz(:) = 0.d0

      ! re-allocates specfem array
      if (allocated(force_stf)) deallocate(force_stf)
      allocate(force_stf(NSOURCES),stat=ier)
      force_stf(:) = 0

      !! VM VM set the size of user_source_time_function
      if (USE_EXTERNAL_SOURCE_FILE) then
        NSTEP_STF = NSTEP
        NSOURCES_STF = NSOURCES
      else !! We don't need the array user_source_time_function : use a small dummy array
        NSTEP_STF = 1
        NSOURCES_STF = 1
      endif
      !! allocate the array contains the user defined source time function
      allocate(acqui_simu(ievent)%user_source_time_function(NSTEP_STF, NSOURCES_STF),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 487')
      if (ier /= 0) stop 'error allocating arrays for user sources time function'

      ! 3/ Read the file describing the sources
      select case (trim(acqui_simu(ievent)%source_type))

        case('moment','force')

          allocate(lat(NSOURCES),long(NSOURCES),depth(NSOURCES),moment_tensor(6,NSOURCES),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 488')
          ! read all the sources
          if (USE_FORCE_POINT_SOURCE) then
            ! point forces
            if (myrank == 0) then
              ! only main process reads in FORCESOLUTION file
              call get_force(filename,acqui_simu(ievent)%tshift,acqui_simu(ievent)%hdur, &
                             lat,long,depth,NSOURCES,min_tshift,force_stf,factor_force_source, &
                             Fx,Fy,Fz,acqui_simu(ievent)%user_source_time_function)
            endif
            ! broadcasts specific point force infos
            call bcast_all_i(force_stf,NSOURCES)
            call bcast_all_dp(factor_force_source,NSOURCES)
            call bcast_all_dp(Fx,NSOURCES)
            call bcast_all_dp(Fy,NSOURCES)
            call bcast_all_dp(Fz,NSOURCES)
          else
            ! CMT moment tensors
            if (myrank == 0) then
              ! only main process reads in CMTSOLUTION file
              call get_cmt(filename,yr,jda,mo,da,ho,mi,sec, &
                           acqui_simu(ievent)%tshift,acqui_simu(ievent)%hdur,lat,long,depth,moment_tensor, &
                           DT,NSOURCES,min_tshift,acqui_simu(ievent)%user_source_time_function)
            endif
            ! broadcasts specific moment tensor infos
            call bcast_all_dp(moment_tensor,6*NSOURCES)
          endif

          ! broadcasts general source information read on the main to the nodes
          call bcast_all_dp(acqui_simu(ievent)%tshift,NSOURCES)
          call bcast_all_dp(acqui_simu(ievent)%hdur,NSOURCES)
          call bcast_all_dp(lat,NSOURCES)
          call bcast_all_dp(long,NSOURCES)
          call bcast_all_dp(depth,NSOURCES)
          call bcast_all_singledp(min_tshift)
          call bcast_all_cr(acqui_simu(ievent)%user_source_time_function,NSOURCES_STF*NSTEP_STF)

          ! get the moment tensor
          Mzz(:) = + moment_tensor(1,:)
          Mxx(:) = + moment_tensor(3,:)
          Myy(:) = + moment_tensor(2,:)
          Mxz(:) = + moment_tensor(5,:)
          Myz(:) = - moment_tensor(4,:)
          Mxy(:) = - moment_tensor(6,:)

          ! loop on all the sources
          do isrc = 1,NSOURCES

            ! get z target coordinate, depending on the topography
            if (.not. USE_SOURCES_RECEIVERS_Z) depth(isrc) = depth(isrc)*1000.0d0
            call get_elevation_and_z_coordinate(long(isrc),lat(isrc),x_target_source(isrc),y_target_source(isrc), &
                                                       z_target_source(isrc),elevation,depth(isrc))

          enddo

          deallocate(lat,long,depth,moment_tensor)

        case ('shot')

           !! set point source
           acqui_simu(ievent)%Xs(:) = acqui_simu(ievent)%xshot
           acqui_simu(ievent)%Ys(:) = acqui_simu(ievent)%yshot
           acqui_simu(ievent)%Zs(:) = acqui_simu(ievent)%zshot
           x_target_source(:) = acqui_simu(ievent)%xshot
           y_target_source(:) = acqui_simu(ievent)%yshot
           z_target_source(:) = acqui_simu(ievent)%zshot

           !! define shot mechanism
           Mzz(:) = acqui_simu(ievent)%shot_ampl * 1.d-7
           Mxx(:) = acqui_simu(ievent)%shot_ampl * 1.d-7
           Myy(:) = acqui_simu(ievent)%shot_ampl * 1.d-7
           Mxz(:) = 0.
           Myz(:) = 0.
           Mxy(:) = 0.

           !! read source time function
           if (myrank == 0) then
              open(IINN, file=trim(acqui_simu(ievent)%source_wavelet_file))
              do it=1,acqui_simu(ievent)%Nt_data
                 read(IINN, *) dt_dummy, acqui_simu(ievent)%user_source_time_function(it,1)
              enddo
              close(IINN)
           endif

           call MPI_BCAST(acqui_simu(ievent)%user_source_time_function,acqui_simu(ievent)%Nt_data, &
                CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
           USE_FORCE_POINT_SOURCE=.false.

        case default
          !! define here reading of source file, and define accordingly the arrays :
          !  - factor_force_source,Fx,Fy,Fz  if USE_FORCE_POINT_SOURCE==.true.
          !  - Mxx,Myy,Mzz,Mxz,Myz,Mxy       if USE_FORCE_POINT_SOURCE==.false.
          !  - x_target_source,y_target_source,z_target_source,min_tshift,user_source_time_function,hdur
           write(*,*) " Your source is not implemented "
           stop

      end select

      !4/ Locate the sources in the mesh

      ! loop on all the sources
      do isrc = 1,NSOURCES

        call locate_point_in_mesh(x_target_source(isrc), y_target_source(isrc), z_target_source(isrc), &
                SOURCES_CAN_BE_BURIED, elemsize_max_glob, &
                acqui_simu(ievent)%ispec_selected_source(isrc), xi_source(isrc), eta_source(isrc), gamma_source(isrc), &
                acqui_simu(ievent)%Xs(isrc), acqui_simu(ievent)%Ys(isrc), acqui_simu(ievent)%Zs(isrc), &
                idomain,nu_source(:,:,isrc),final_distance_squared)

        ! synchronize all the processes to make sure all the estimates are available
        call synchronize_all()

        call locate_MPI_slice_and_bcast_to_all_single(x_target_source(isrc), y_target_source(isrc), z_target_source(isrc), &
                   acqui_simu(ievent)%Xs(isrc), acqui_simu(ievent)%Ys(isrc),acqui_simu(ievent)%Zs(isrc), &
                   xi_source(isrc), eta_source(isrc), gamma_source(isrc), &
                   acqui_simu(ievent)%ispec_selected_source(isrc),acqui_simu(ievent)%islice_selected_source(isrc), &
                   final_distance_source, idomain,nu_source(:,:,isrc))

      enddo ! end of loop on all the sources

      tshift_src(:) = acqui_simu(ievent)%tshift(:)
      hdur(:) = acqui_simu(ievent)%hdur(:)
      min_tshift_src_original = min_tshift
      islice_selected_source(:) = acqui_simu(ievent)%islice_selected_source(:)
      ispec_selected_source(:) = acqui_simu(ievent)%ispec_selected_source(:)

      call setup_stf_constants()

      acqui_simu(ievent)%hdur(:) = hdur(:)
      acqui_simu(ievent)%hdur_Gaussian = hdur_Gaussian(:)
      acqui_simu(ievent)%tshift(:) = tshift_src(:)
      min_tshift = min_tshift_src_original
      acqui_simu(ievent)%t0 = t0

      nsrc_loc = 0
      do isrc=1, NSOURCES
        if (myrank == acqui_simu(ievent)%islice_selected_source(isrc)) then
          nsrc_loc = nsrc_loc + 1
          !! Warning in this subroutine you must add your case for source source
          call compute_source_coeff(xi_source(isrc), eta_source(isrc), gamma_source(isrc), &
                                    acqui_simu(ievent)%ispec_selected_source(isrc), &
                                    interparray,Mxx(isrc),Myy(isrc),Mzz(isrc),Mxy(isrc), &
                                    Mxz(isrc),Myz(isrc),factor_force_source(isrc),Fx(isrc),Fy(isrc),Fz(isrc), &
                                    acqui_simu(ievent)%source_type,nu_source(:,:,isrc))
          acqui_simu(ievent)%sourcearrays(isrc,:,:,:,:) = interparray(:,:,:,:)
        endif
      enddo

      acqui_simu(ievent)%nsources_local = nsrc_loc
      acqui_simu(ievent)%nsources_tot = NSOURCES

      deallocate(x_target_source,y_target_source,z_target_source)
      deallocate(islice_selected_source,ispec_selected_source)
      deallocate(Mxx,Myy,Mzz,Mxy,Mxz,Myz,xi_source,eta_source,gamma_source,nu_source)
      deallocate(tshift_src,hdur,hdur_Gaussian)
      deallocate(factor_force_source,Fx,Fy,Fz)
      deallocate(force_stf)

    enddo ! loop over events

    if (myrank == 0) then
      write(INVERSE_LOG_FILE,*) ' READING and LOCATING sources passed '
      call flush_iunit(INVERSE_LOG_FILE)
    endif

  end subroutine get_point_source

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------
! MPI BCAST of acqui_simu
!---------------------------------------------------------------

!  subroutine bcast_all_acqui(acqui_simu, inversion_param, myrank)
!    use my_mpi             !! module from specfem
!    include "precision.h"  !! from specfem
!    type(acqui), dimension(:), intent(inout)  :: acqui_simu
!    type(inver),               intent(inout)  :: inversion_param
!    integer,                   intent(in)     :: myrank
!    integer                                   :: i,ier

!    if (myrank == 0) NEVENT=acqui_simu(1)%nevent_tot
!    call MPI_BCAST(NEVENT,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)

!    do i=1,NEVENT

       ! events
!       acqui_simu(i)%nevent_tot=NEVENT
       ! bcast file paths
!       call MPI_BCAST(acqui_simu(i)%source_file,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(acqui_simu(i)%traction_dir,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(acqui_simu(i)%data_file_gather,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(acqui_simu(i)%source_type,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(acqui_simu(i)%adjoint_source_type,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(acqui_simu(i)%station_file,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(acqui_simu(i)%event_name,MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(acqui_simu(i)%component, 6,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)

!       call MPI_BCAST(acqui_simu(i)%Nt_data,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(acqui_simu(i)%dt_data,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)

       ! bcast source parameters
!       call MPI_BCAST(acqui_simu(i)%Xs,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(acqui_simu(i)%Ys,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(acqui_simu(i)%Zs,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)

!       call MPI_BCAST(acqui_simu(i)%Mxx,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(acqui_simu(i)%Myy,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(acqui_simu(i)%Mzz,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(acqui_simu(i)%Mxy,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(acqui_simu(i)%Mxz,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(acqui_simu(i)%Myz,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)

!       call MPI_BCAST(acqui_simu(i)%Fx,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(acqui_simu(i)%Fy,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(acqui_simu(i)%Fz,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)

!       call MPI_BCAST(acqui_simu(i)%t_shift,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(acqui_simu(i)%hdur,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(acqui_simu(i)%source_wavelet_file, MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)

!       call MPI_BCAST(acqui_simu(i)%external_source_wavelet, 1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
!       if (acqui_simu(i)%external_source_wavelet) then
!         if (myrank > 0) then
!           allocate(acqui_simu(i)%source_wavelet(acqui_simu(i)%Nt_data,1))
!         endif
!         call MPI_BCAST(acqui_simu(i)%source_wavelet,acqui_simu(i)%Nt_data, CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
!       endif

       ! stations

!       if (myrank == 0) nsta_tot=acqui_simu(i)%nsta_tot
!       call  MPI_BCAST(nsta_tot,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
!       acqui_simu(i)%nsta_tot=nsta_tot
!       if (myrank > 0) then
!          allocate(acqui_simu(i)%station_name(nsta_tot),acqui_simu(i)%network_name(nsta_tot))
!          allocate(acqui_simu(i)%position_station(3,nsta_tot))
!       endif

!       call MPI_BCAST(acqui_simu(i)%station_name, nsta_tot* MAX_LENGTH_STATION_NAME,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(acqui_simu(i)%network_name, nsta_tot* MAX_LENGTH_NETWORK_NAME,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(acqui_simu(i)%position_station, 3*nsta_tot,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)


       ! inverse params
!       call MPI_BCAST(inversion_param%Niter,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(inversion_param%Niter_wolfe,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(inversion_param%max_history_bfgs,1,MPI_INTEGER,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(fl,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(fh,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(inversion_param%max_relative_pert,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(inversion_param%relat_grad,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(inversion_param%relat_cost,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(inversion_param%output_model,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(inversion_param%input_fd_model,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(inversion_param%input_sem_model,1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(inversion_param%use_taper, 1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(inversion_param%shin_precond, 1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(inversion_param%energy_precond, 1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(inversion_param%z2_precond, 1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(inversion_param%dump_model_at_each_iteration, 1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(inversion_param%dump_gradient_at_each_iteration, 1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(inversion_param%dump_descent_direction_at_each_iteration, 1,MPI_LOGICAL,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(inversion_param%param_family, MAX_LEN_STRING,MPI_CHARACTER,0,my_local_mpi_comm_world,ier)

       ! user taper on gradient (MASK)
!       call MPI_BCAST(inversion_param%xmin_taper,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(inversion_param%xmax_taper,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(inversion_param%ymin_taper,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(inversion_param%ymax_taper,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(inversion_param%zmin_taper,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)
!       call MPI_BCAST(inversion_param%zmax_taper,1,CUSTOM_MPI_TYPE,0,my_local_mpi_comm_world,ier)

!    enddo

!  end subroutine bcast_all_acqui

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine header_output(iunit)

  implicit none
  integer, intent(in) :: iunit

  ! header ANSI picture
  write(iunit,*) '//////////////////////////////////////////////////////////////////////////////'
  write(iunit,*)
  write(iunit,*) '  cONXXXXXX  x0x:   ckx:   ck0d;  dXNX'
  write(iunit,*) '  cKMNkdddo  OWNd  :kWNx;  xNWk;  xWWK'
  write(iunit,*) '  cKMX       dNWk  xNMMNd  OMNo   xWWO'
  write(iunit,*) '  cKMWNXXXK  lKM0ldNW00WXolKMKc   xWMW'
  write(iunit,*) '  cKMNkdddo   kWN0XWOcc0WX0NWk,   xWWK'
  write(iunit,*) '  cKMK        dNMMWKc  lKWMMNd    xWWk'
  write(iunit,*) '  cKMN        cKMMXo    lXMMKc    xWWk'
  write(iunit,*) '  ox00         xNXd      oKNk     oOOo               Inverse Problem for Model'
  write(iunit,*)
  write(iunit,*) '//////////////////////////////////////////////////////////////////////////////'
  write(iunit,*)
  call flush_iunit(iunit)

  end subroutine header_output

end module input_output
