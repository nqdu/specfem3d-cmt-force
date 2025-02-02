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


  subroutine compute_seismograms()

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,ZERO

  use specfem_par, only: SIMULATION_TYPE,NGLOB_AB,NSPEC_AB,ibool,NGLOB_ADJOINT, &
    seismo_current, &
    ispec_selected_source,ispec_selected_rec, &
    number_receiver_global,nrec_local, &
    nu_source,nu_rec, &
    hxir_store,hetar_store,hgammar_store, &
    USE_TRICK_FOR_BETTER_PRESSURE, &
    SAVE_SEISMOGRAMS_DISPLACEMENT,SAVE_SEISMOGRAMS_VELOCITY,SAVE_SEISMOGRAMS_ACCELERATION,SAVE_SEISMOGRAMS_PRESSURE

  ! seismograms
  use specfem_par, only: seismograms_d,seismograms_v,seismograms_a,seismograms_p

  ! wavefields
  use specfem_par_acoustic, only: ispec_is_acoustic,potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
    b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic
  use specfem_par_elastic, only: ispec_is_elastic,displ,veloc,accel, &
    b_displ,b_veloc,b_accel
  use specfem_par_poroelastic, only: ispec_is_poroelastic,displs_poroelastic,velocs_poroelastic,accels_poroelastic, &
    b_displs_poroelastic,b_velocs_poroelastic,b_accels_poroelastic

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ):: displ_element,veloc_element,accel_element

  ! interpolated wavefield values
  double precision :: dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd,pd
  double precision,dimension(NDIM,NDIM) :: rotation_seismo

  ! receiver Lagrange interpolators
  double precision,dimension(NGLLX) :: hxir
  double precision,dimension(NGLLY) :: hetar
  double precision,dimension(NGLLZ) :: hgammar

  integer :: irec_local,irec
  integer :: ispec

  ! flag to indicate that traces for kernel runs are taken from adjoint wavefields instead of backward/reconstructed wavefields;
  ! useful for debugging.
  ! default (.false.) is to output backward/reconstructed wavefield
  logical, parameter :: OUTPUT_ADJOINT_WAVEFIELD = .false.

  ! loops over local receivers
  do irec_local = 1,nrec_local

    ! initializes wavefield values
    dxd = ZERO
    dyd = ZERO
    dzd = ZERO

    vxd = ZERO
    vyd = ZERO
    vzd = ZERO

    axd = ZERO
    ayd = ZERO
    azd = ZERO

    pd  = ZERO

    ! gets local receiver interpolators
    ! (1-D Lagrange interpolators)
    hxir(:) = hxir_store(:,irec_local)
    hetar(:) = hetar_store(:,irec_local)
    hgammar(:) = hgammar_store(:,irec_local)

    ! gets global number of that receiver
    irec = number_receiver_global(irec_local)

    ! spectral element in which the receiver is located
    if (SIMULATION_TYPE == 2) then
      ! adjoint "receivers" are located at CMT source positions
      ! note: we take here xi_source,.. when FASTER_RECEIVERS_POINTS_ONLY is set
      ispec = ispec_selected_source(irec)
    else
      ! receiver located at station positions
      ispec = ispec_selected_rec(irec)
    endif

    ! calculates interpolated wavefield values at receiver positions
    select case (SIMULATION_TYPE)
    case (1,2)
      ! forward simulations & pure adjoint simulations
      ! wavefields stored in displ,veloc,accel

      ! elastic wave field
      if (ispec_is_elastic(ispec)) then
        ! interpolates displ/veloc/accel at receiver locations
        call compute_interpolated_dva_viscoelast(displ,veloc,accel,NGLOB_AB, &
                                                 ispec,NSPEC_AB,ibool, &
                                                 hxir,hetar,hgammar, &
                                                 dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd,pd)
      endif ! elastic

      ! acoustic wave field
      if (ispec_is_acoustic(ispec)) then
        ! displacement vector
        call compute_gradient_in_acoustic(ispec,potential_acoustic,displ_element)

        ! velocity vector
        call compute_gradient_in_acoustic(ispec,potential_dot_acoustic,veloc_element)

        ! acceleration vector
        call compute_gradient_in_acoustic(ispec,potential_dot_dot_acoustic,accel_element)

        ! interpolates displ/veloc/accel/pressure at receiver locations
        call compute_interpolated_dva_acoust(displ_element,veloc_element,accel_element, &
                                             potential_dot_dot_acoustic,potential_acoustic,NGLOB_AB, &
                                             ispec,NSPEC_AB,ibool, &
                                             hxir,hetar,hgammar, &
                                             dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd,pd,USE_TRICK_FOR_BETTER_PRESSURE)
      endif ! acoustic

      ! poroelastic wave field
      if (ispec_is_poroelastic(ispec)) then
        ! interpolates displ/veloc/accel at receiver locations
        call compute_interpolated_dva_viscoelast(displs_poroelastic,velocs_poroelastic,accels_poroelastic,NGLOB_AB, &
                                                 ispec,NSPEC_AB,ibool, &
                                                 hxir,hetar,hgammar, &
                                                 dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd,pd)
      endif ! poroelastic

    case (3)
      ! adjoint/kernel simulations
      ! reconstructed forward wavefield stored in b_displ, b_veloc, b_accel

      ! elastic wave field
      if (ispec_is_elastic(ispec)) then
        if (OUTPUT_ADJOINT_WAVEFIELD) then
          ! adjoint field: interpolates displ/veloc/accel at receiver locations
          call compute_interpolated_dva_viscoelast(displ,veloc,accel,NGLOB_ADJOINT, &
                                                   ispec,NSPEC_AB,ibool, &
                                                   hxir,hetar,hgammar, &
                                                   dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd,pd)
        else
          ! backward field: interpolates displ/veloc/accel at receiver locations
          call compute_interpolated_dva_viscoelast(b_displ,b_veloc,b_accel,NGLOB_ADJOINT, &
                                                   ispec,NSPEC_AB,ibool, &
                                                   hxir,hetar,hgammar, &
                                                   dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd,pd)
        endif
      endif ! elastic

      ! acoustic wave field
      if (ispec_is_acoustic(ispec)) then
        if (OUTPUT_ADJOINT_WAVEFIELD) then
          ! adjoint field: displacement vector
          call compute_gradient_in_acoustic(ispec,potential_acoustic,displ_element)

          ! adjoint field: velocity vector
          call compute_gradient_in_acoustic(ispec,potential_dot_acoustic,veloc_element)

          ! adjoint field: acceleration vector
          call compute_gradient_in_acoustic(ispec,potential_dot_dot_acoustic,accel_element)

          ! adjoint field: interpolates displ/veloc/accel/pressure at receiver locations
          call compute_interpolated_dva_acoust(displ_element,veloc_element,accel_element, &
                                               potential_dot_dot_acoustic,potential_acoustic,NGLOB_ADJOINT, &
                                               ispec,NSPEC_AB,ibool, &
                                               hxir,hetar,hgammar, &
                                               dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd,pd,USE_TRICK_FOR_BETTER_PRESSURE)
        else
          ! backward field: displacement vector
          call compute_gradient_in_acoustic(ispec,b_potential_acoustic,displ_element)

          ! backward field: velocity vector
          call compute_gradient_in_acoustic(ispec,b_potential_dot_acoustic,veloc_element)

          ! backward field: acceleration vector
          call compute_gradient_in_acoustic(ispec,b_potential_dot_dot_acoustic,accel_element)

          ! backward field: interpolates displ/veloc/accel/pressure at receiver locations
          call compute_interpolated_dva_acoust(displ_element,veloc_element,accel_element, &
                                               b_potential_dot_dot_acoustic,b_potential_acoustic,NGLOB_ADJOINT, &
                                               ispec,NSPEC_AB,ibool, &
                                               hxir,hetar,hgammar, &
                                               dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd,pd,USE_TRICK_FOR_BETTER_PRESSURE)
        endif
      endif ! acoustic

      ! poroelastic wavefield
      if (ispec_is_poroelastic(ispec)) then
        ! outputs wavefield from solid phase: displs/velocs/accels
        if (OUTPUT_ADJOINT_WAVEFIELD) then
          ! adjoint field: interpolates displ/veloc/accel at receiver locations
          call compute_interpolated_dva_viscoelast(displs_poroelastic,velocs_poroelastic,accels_poroelastic,NGLOB_ADJOINT, &
                                                   ispec,NSPEC_AB,ibool, &
                                                   hxir,hetar,hgammar, &
                                                   dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd,pd)
        else
          ! backward field: interpolates displ/veloc/accel at receiver locations
          call compute_interpolated_dva_viscoelast(b_displs_poroelastic,b_velocs_poroelastic,b_accels_poroelastic,NGLOB_ADJOINT, &
                                                   ispec,NSPEC_AB,ibool, &
                                                   hxir,hetar,hgammar, &
                                                   dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd,pd)
        endif
      endif

    end select ! SIMULATION_TYPE

    if (SIMULATION_TYPE == 2) then
      ! adjoint simulations
      ! adjoint "receiver" N/E/Z orientations given by nu_source array
      rotation_seismo(:,:) = nu_source(:,:,irec)
    else
      rotation_seismo(:,:) = nu_rec(:,:,irec)
    endif

    ! we only store if needed
    ! note: current index is seismo_current, this allows to store arrays only up to nlength_seismogram
    !       which could be used to limit the allocation size of these arrays for a large number of receivers
    if (SAVE_SEISMOGRAMS_DISPLACEMENT) &
      seismograms_d(:,irec_local,seismo_current) = real(rotation_seismo(:,1)*dxd + rotation_seismo(:,2)*dyd &
                                                      + rotation_seismo(:,3)*dzd,kind=CUSTOM_REAL)

    if (SAVE_SEISMOGRAMS_VELOCITY) &
      seismograms_v(:,irec_local,seismo_current) = real(rotation_seismo(:,1)*vxd + rotation_seismo(:,2)*vyd &
                                                      + rotation_seismo(:,3)*vzd,kind=CUSTOM_REAL)

    if (SAVE_SEISMOGRAMS_ACCELERATION) &
      seismograms_a(:,irec_local,seismo_current) = real(rotation_seismo(:,1)*axd + rotation_seismo(:,2)*ayd &
                                                      + rotation_seismo(:,3)*azd,kind=CUSTOM_REAL)

    ! only one scalar in the case of pressure
    if (SAVE_SEISMOGRAMS_PRESSURE) &
      seismograms_p(1,irec_local,seismo_current) = real(pd,kind=CUSTOM_REAL)

  enddo ! nrec_local

  end subroutine compute_seismograms

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_seismograms_strain()

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM

  use specfem_par, only: SIMULATION_TYPE,ACOUSTIC_SIMULATION,ELASTIC_SIMULATION, &
    NGLOB_AB,ibool, &
    seismo_current, &
    ispec_selected_source,ispec_selected_rec, &
    number_receiver_global,nrec_local, &
    hxir_store,hetar_store,hgammar_store, &
    SAVE_SEISMOGRAMS_STRAIN

  ! seismograms
  use specfem_par, only: seismograms_eps

  !use specfem_par, only: nu_source,nu_rec

  ! wavefields
  use specfem_par_acoustic, only: ispec_is_acoustic,potential_acoustic,b_potential_acoustic
  use specfem_par_elastic, only: ispec_is_elastic,displ,b_displ
  use specfem_par_poroelastic, only: ispec_is_poroelastic,displs_poroelastic,b_displs_poroelastic

  ! GPU simulations
  use specfem_par, only: Mesh_pointer,GPU_MODE

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ):: displ_element
  ! receiver Lagrange interpolators
  double precision,dimension(NGLLX) :: hxir
  double precision,dimension(NGLLY) :: hetar
  double precision,dimension(NGLLZ) :: hgammar

  ! strain
  real(kind=CUSTOM_REAL),dimension(NDIM,NDIM):: eps_s

  integer :: irec_local,irec
  integer :: ispec,i,j,k,iglob

  ! flag to indicate that traces for kernel runs are taken from adjoint wavefields instead of backward/reconstructed wavefields;
  ! useful for debugging.
  ! default (.false.) is to output backward/reconstructed wavefield
  logical, parameter :: OUTPUT_ADJOINT_WAVEFIELD = .false.

  ! safety check
  if (.not. SAVE_SEISMOGRAMS_STRAIN) return

  ! strain is only implemented on CPU so far, thus we need to transfer wavefields before computing strain
  if (GPU_MODE) then
    ! transfers displacement to the CPU
    if (ELASTIC_SIMULATION) then
      call transfer_displ_from_device(NDIM*NGLOB_AB, displ, Mesh_pointer)
      if (SIMULATION_TYPE == 3 .and. .not. OUTPUT_ADJOINT_WAVEFIELD) &
        call transfer_b_displ_from_device(NDIM*NGLOB_AB,b_displ,Mesh_pointer)
    endif
    if (ACOUSTIC_SIMULATION) then
      call transfer_potential_ac_from_device(NGLOB_AB,potential_acoustic,Mesh_pointer)
      if (SIMULATION_TYPE == 3 .and. .not. OUTPUT_ADJOINT_WAVEFIELD) &
        call transfer_b_potential_ac_from_device(NGLOB_AB,b_potential_acoustic,Mesh_pointer)
    endif
    ! poroelastic simulations: only supported on CPU, so no transfer needed yet...
  endif

  ! loops over local receivers
  do irec_local = 1,nrec_local

    ! initializes strain
    eps_s(:,:) = 0._CUSTOM_REAL

    ! gets local receiver interpolators
    ! (1-D Lagrange interpolators)
    hxir(:) = hxir_store(:,irec_local)
    hetar(:) = hetar_store(:,irec_local)
    hgammar(:) = hgammar_store(:,irec_local)

    ! gets global number of that receiver
    irec = number_receiver_global(irec_local)

    ! spectral element in which the receiver is located
    if (SIMULATION_TYPE == 2) then
      ! adjoint "receivers" are located at the source positions
      ispec = ispec_selected_source(irec)
    else
      ! receiver located at station positions
      ispec = ispec_selected_rec(irec)
    endif

    ! calculates interpolated wavefield values at receiver positions
    select case (SIMULATION_TYPE)
    case (1,2)
      ! forward simulations & pure adjoint simulations
      ! wavefields stored in displ,veloc,accel

      ! strain
      ! elastic wave field
      if (ispec_is_elastic(ispec)) then
        ! gets elements displacement field
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              iglob = ibool(i,j,k,ispec)
              displ_element(:,i,j,k) = displ(:,iglob)
            enddo
          enddo
        enddo
        ! computes strain
        call compute_interpolated_strain(ispec,displ_element,hxir,hetar,hgammar,eps_s)
      endif ! elastic

      ! acoustic wave field
      if (ispec_is_acoustic(ispec)) then
        ! displacement vector
        call compute_gradient_in_acoustic(ispec,potential_acoustic,displ_element)

        ! computes strain
        call compute_interpolated_strain(ispec,displ_element,hxir,hetar,hgammar,eps_s)
      endif ! acoustic

      ! poroelastic wave field
      if (ispec_is_poroelastic(ispec)) then
        ! gets elements displacement field
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              iglob = ibool(i,j,k,ispec)
              displ_element(:,i,j,k) = displs_poroelastic(:,iglob)
            enddo
          enddo
        enddo
        ! computes strain
        call compute_interpolated_strain(ispec,displ_element,hxir,hetar,hgammar,eps_s)
      endif ! poroelastic

    case (3)
      ! adjoint/kernel simulations
      ! reconstructed forward wavefield stored in b_displ, b_veloc, b_accel

      ! strain
      ! elastic wave field
      if (ispec_is_elastic(ispec)) then
        if (OUTPUT_ADJOINT_WAVEFIELD) then
          ! adjoint field: interpolates displ/veloc/accel at receiver locations
          ! gets elements displacement field
          do k = 1,NGLLZ
            do j = 1,NGLLY
              do i = 1,NGLLX
                iglob = ibool(i,j,k,ispec)
                displ_element(:,i,j,k) = displ(:,iglob)
              enddo
            enddo
          enddo
        else
          ! backward field: interpolates displ/veloc/accel at receiver locations
          ! gets elements displacement field
          do k = 1,NGLLZ
            do j = 1,NGLLY
              do i = 1,NGLLX
                iglob = ibool(i,j,k,ispec)
                displ_element(:,i,j,k) = b_displ(:,iglob)
              enddo
            enddo
          enddo
        endif
        ! computes strain
        call compute_interpolated_strain(ispec,displ_element,hxir,hetar,hgammar,eps_s)
      endif ! elastic

      ! acoustic wave field
      if (ispec_is_acoustic(ispec)) then
        if (OUTPUT_ADJOINT_WAVEFIELD) then
          ! adjoint field: displacement vector
          call compute_gradient_in_acoustic(ispec,potential_acoustic,displ_element)
        else
          ! backward field: displacement vector
          call compute_gradient_in_acoustic(ispec,b_potential_acoustic,displ_element)
        endif
        ! computes strain
        call compute_interpolated_strain(ispec,displ_element,hxir,hetar,hgammar,eps_s)
      endif ! acoustic

      ! poroelastic wavefield
      if (ispec_is_poroelastic(ispec)) then
        ! outputs wavefield from solid phase: displs/velocs/accels
        if (OUTPUT_ADJOINT_WAVEFIELD) then
          ! adjoint field: interpolates displ/veloc/accel at receiver locations
          ! gets elements displacement field
          do k = 1,NGLLZ
            do j = 1,NGLLY
              do i = 1,NGLLX
                iglob = ibool(i,j,k,ispec)
                displ_element(:,i,j,k) = displs_poroelastic(:,iglob)
              enddo
            enddo
          enddo
        else
          ! backward field: interpolates displ/veloc/accel at receiver locations
          ! gets elements displacement field
          do k = 1,NGLLZ
            do j = 1,NGLLY
              do i = 1,NGLLX
                iglob = ibool(i,j,k,ispec)
                displ_element(:,i,j,k) = b_displs_poroelastic(:,iglob)
              enddo
            enddo
          enddo
        endif
        ! computes strain
        call compute_interpolated_strain(ispec,displ_element,hxir,hetar,hgammar,eps_s)
      endif

    end select ! SIMULATION_TYPE

    ! todo: check if rotation is required
    !
    ! --- from global code: ----
    ! un-rotated
    !eps_loc_new(:,:) = eps_loc(:,:)
    !
    ! rotates from global x-y-z to the local coordinates (n-e-z):  eps_new = P*eps*P'
    ! nu is the rotation matrix from ECEF to local N-E-UP as defined.
    ! thus, if the nu is the rotation matrix that transforms coordinates from the global system (x,y,z) to the local
    ! coordinate system (N,E,V), e.g., a tensor is transformed as
    ! T_L = \nu * T_g * \nu^T
    !
    ! global -> local (n-e-up)
    ! eps_xx -> eps_nn
    ! eps_yy -> eps_ee
    ! eps_zz -> eps_zz (z in radial direction up)
    ! eps_xy -> eps_ne
    ! eps_xz -> eps_nz
    ! eps_yz -> eps_ez
    !eps_loc_new(:,:) = matmul(matmul(nu_rec(:,:,irec),eps_loc(:,:)), transpose(nu_rec(:,:,irec)))
    !
    ! --- for Cartesian code: ----
    ! this would become:
    !if (SIMULATION_TYPE == 2) then
    !  ! adjoint simulations
    !  ! adjoint "receiver" N/E/Z orientations given by nu_source array
    !  rotation_seismo(:,:) = nu_source(:,:,irec)
    !else
    !  rotation_seismo(:,:) = nu_rec(:,:,irec)
    !endif
    !eps_s(:,:) = matmul(matmul(rotation_seismo(:,:),eps_s(:,:)), transpose(rotation_seismo(:,:)))

    ! stores strain
    seismograms_eps(:,:,irec_local,seismo_current) = eps_s(:,:)

  enddo ! nrec_local

  end subroutine compute_seismograms_strain

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_seismograms_moment_adjoint()

  use constants, only: myrank,CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM

  use specfem_par, only: SIMULATION_TYPE,NGLOB_AB,NSPEC_AB,ibool, &
    deltat,DT,t0,NSTEP,it, &
    seismo_current,seismo_offset,NTSTEP_BETWEEN_OUTPUT_SAMPLE, &
    ispec_selected_source, &
    number_receiver_global,nrec_local, &
    Mxx,Myy,Mzz,Mxy,Mxz,Myz,tshift_src, &
    hprime_xx,hprime_yy,hprime_zz, &
    hxir_store,hetar_store,hgammar_store, &
    hpxir_store,hpetar_store,hpgammar_store, &
    ELASTIC_SIMULATION

  use specfem_par, only: GPU_MODE, Mesh_pointer

  ! strain "seismogram"
  !use specfem_par, only: seismograms_eps

  ! for source derivatives
  use specfem_par, only: Mxx_der,Myy_der,Mzz_der,Mxy_der,Mxz_der,Myz_der,sloc_der

  ! wavefield
  use specfem_par_elastic, only: ispec_is_elastic,displ

  implicit none

  ! local parameters
  integer :: irec_local,irec,idx
  integer :: i,j,k,iglob,ispec
  ! adjoint locals
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ):: displ_element
  real(kind=CUSTOM_REAL),dimension(NDIM,NDIM):: eps_s
  real(kind=CUSTOM_REAL),dimension(NDIM):: eps_m_s
  real(kind=CUSTOM_REAL):: stf_deltat
  double precision :: stf,time_source_dble
  ! receiver Lagrange interpolators
  double precision,dimension(NGLLX) :: hxir
  double precision,dimension(NGLLY) :: hetar
  double precision,dimension(NGLLZ) :: hgammar
  double precision :: hpxir(NGLLX),hpetar(NGLLY),hpgammar(NGLLZ)

  double precision, external :: get_stf_viscoelastic

  ! checks if anything to do
  if (SIMULATION_TYPE /= 2) return
  if (.not. ELASTIC_SIMULATION) return

  ! strain and moment derivatives are computed here on CPU routines,
  ! thus transfers displacement to the CPU
  if (GPU_MODE) call transfer_displ_from_device(NDIM*NGLOB_AB, displ, Mesh_pointer)

  ! checks index bounds
  idx = seismo_offset + seismo_current
  if (idx < 1 .or. idx > NSTEP/NTSTEP_BETWEEN_OUTPUT_SAMPLE) &
    call exit_mpi(myrank,'Error: seismograms_eps has wrong current index')

  ! loops over local receivers
  do irec_local = 1,nrec_local
    ! gets global number of that receiver
    irec = number_receiver_global(irec_local)

    ! spectral element in which the receiver is located
    ! adjoint "receivers" are located at CMT source positions
    ! note: we take here xi_source,.. when FASTER_RECEIVERS_POINTS_ONLY is set
    ispec = ispec_selected_source(irec)

    ! additional calculations for pure adjoint simulations
    ! computes derivatives of source parameters

    ! elastic wave field
    if (ispec_is_elastic(ispec)) then
      ! gets local receiver interpolators
      ! (1-D Lagrange interpolators)
      hxir(:) = hxir_store(:,irec_local)
      hetar(:) = hetar_store(:,irec_local)
      hgammar(:) = hgammar_store(:,irec_local)

      ! gets derivatives of local receiver interpolators
      hpxir(:) = hpxir_store(:,irec_local)
      hpetar(:) = hpetar_store(:,irec_local)
      hpgammar(:) = hpgammar_store(:,irec_local)

      ! strain
      ! gets elements displacement field
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob = ibool(i,j,k,ispec)
            displ_element(:,i,j,k) = displ(:,iglob)
          enddo
        enddo
      enddo
      ! computes strain
      call compute_interpolated_strain(ispec,displ_element,hxir,hetar,hgammar,eps_s)

      ! stores strain value
      ! not needed here, as we already compute strain seismograms in routine compute_seismograms();
      ! left here for checking, just in case.
      !seismograms_eps(:,:,irec_local,seismo_current) = eps_s(:,:)

      ! computes the integrated derivatives of source parameters (M_jk and X_s)
      call compute_adj_source_frechet(displ,NGLOB_AB, &
                                      ispec,NSPEC_AB,ibool, &
                                      Mxx(irec),Myy(irec),Mzz(irec), &
                                      Mxy(irec),Mxz(irec),Myz(irec), &
                                      eps_m_s, &
                                      hxir,hetar,hgammar, &
                                      hpxir,hpetar,hpgammar, &
                                      hprime_xx,hprime_yy,hprime_zz)

      ! source time function value
      time_source_dble = dble(NSTEP-it) * DT - t0 - tshift_src(irec)

      ! determines source time function value
      stf = get_stf_viscoelastic(time_source_dble,irec,NSTEP-it+1)

      stf_deltat = real(stf * deltat * NTSTEP_BETWEEN_OUTPUT_SAMPLE,kind=CUSTOM_REAL)

      ! integrated moment tensor derivatives
      Mxx_der(irec_local) = Mxx_der(irec_local) + eps_s(1,1) * stf_deltat
      Myy_der(irec_local) = Myy_der(irec_local) + eps_s(2,2) * stf_deltat
      Mzz_der(irec_local) = Mzz_der(irec_local) + eps_s(3,3) * stf_deltat
      Mxy_der(irec_local) = Mxy_der(irec_local) + 2 * eps_s(1,2) * stf_deltat
      Mxz_der(irec_local) = Mxz_der(irec_local) + 2 * eps_s(1,3) * stf_deltat
      Myz_der(irec_local) = Myz_der(irec_local) + 2 * eps_s(2,3) * stf_deltat

      ! source location derivative
      sloc_der(:,irec_local) = sloc_der(:,irec_local) + eps_m_s(:) * stf_deltat

    endif ! elastic

  enddo ! nrec_local

  end subroutine compute_seismograms_moment_adjoint









