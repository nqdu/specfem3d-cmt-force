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


  subroutine finalize_simulation()

  use manager_adios
  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic
  use pml_par
  use gravity_perturbation, only: gravity_output, GRAVITY_SIMULATION

  ! hdf5 i/o server
  use io_server_hdf5, only: finalize_io_server

  implicit none

  ! hdf5 i/o server
  ! checks if anything to do
  if (.not. IO_compute_task) then
    ! finalizes MPI subgroup for intercommunication
    if (HDF5_IO_NODES > 0) call finalize_io_server()
    ! all done
    return
  endif

  ! synchronize all processes, waits until all processes have written their seismograms
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'finalizing simulation'
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! write gravity perturbations
  if (GRAVITY_SIMULATION) call gravity_output()

  ! save last frame
  if (SAVE_FORWARD) call save_forward_arrays()

  ! adjoint simulations kernels
  if (SIMULATION_TYPE == 3) call save_adjoint_kernels()

  ! seismograms and source parameter gradients for (pure type=2) adjoint simulation runs
  if (SIMULATION_TYPE == 2) then
    if (nrec_local > 0) then
      ! source gradients  (for sources in elastic domains)
      call save_kernels_source_derivatives()
    endif
  endif

  ! stacey absorbing fields will be reconstructed for adjoint simulations
  ! using snapshot files of wavefields
  if (STACEY_ABSORBING_CONDITIONS) then
    ! closes absorbing wavefield saved/to-be-saved by forward simulations
    if (num_abs_boundary_faces > 0 .and. SAVE_STACEY) then
      ! closes files
      if (ELASTIC_SIMULATION) call close_file_abs(IOABS)
      if (ACOUSTIC_SIMULATION) call close_file_abs(IOABS_AC)
    endif
  endif

  ! coupling
  if (COUPLE_WITH_INJECTION_TECHNIQUE) call couple_with_injection_cleanup()

  ! ADIOS file i/o
  if (ADIOS_ENABLED) then
    call finalize_adios()
  endif

  ! asdf finalizes
  if ((SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) .and. READ_ADJSRC_ASDF) then
    call asdf_cleanup()
  endif

  ! synchronize all
  call synchronize_all()

  ! frees dynamically allocated memory
  call finalize_simulation_cleanup()

  ! close the main output file
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'End of the simulation'
    write(IMAIN,*)
    call flush_IMAIN()
    close(IMAIN)
  endif

  ! synchronize all the processes to make sure everybody has finished
  call synchronize_all()

  ! hdf5 i/o server
  ! finalizes MPI subgroup for intercommunication
  if (HDF5_IO_NODES > 0) call finalize_io_server()

  end subroutine finalize_simulation

!
!-------------------------------------------------------------------------------------------------
!

  subroutine finalize_simulation_cleanup()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_lts

  use fault_solver_dynamic, only: SIMULATION_TYPE_DYN,BC_DYNFLT_free
  use fault_solver_kinematic, only: SIMULATION_TYPE_KIN,BC_KINFLT_free

  !! solving wavefield discontinuity problem with non-split-node scheme
  use wavefield_discontinuity_solver, only: finalize_wavefield_discontinuity

  implicit none

  ! from here on, no gpu data is needed anymore
  ! frees allocated memory on GPU
  if (GPU_MODE) call prepare_cleanup_device(Mesh_pointer,ACOUSTIC_SIMULATION,ELASTIC_SIMULATION,NOISE_TOMOGRAPHY)

  ! C-PML absorbing boundary conditions deallocates C_PML arrays
  if (PML_CONDITIONS) call pml_cleanup()

  ! cleanup for wavefield discontinuity
  if (IS_WAVEFIELD_DISCONTINUITY) call finalize_wavefield_discontinuity()

  ! free arrays
  ! mass matrices
  if (ELASTIC_SIMULATION) then
    deallocate(rmassx)
    deallocate(rmassy)
    deallocate(rmassz)
  endif
  if (ACOUSTIC_SIMULATION) then
    deallocate(rmass_acoustic)
  endif

  ! boundary surfaces
  deallocate(ibelm_xmin)
  deallocate(ibelm_xmax)
  deallocate(ibelm_ymin)
  deallocate(ibelm_ymax)
  deallocate(ibelm_bottom)
  deallocate(ibelm_top)

  ! fixed and roller boundaries
  if(allocated(fixed_bdry_ijk)) then 
    deallocate(fixed_bdry_ijk,fixed_bdry_ispec)
  endif
  if(allocated(roller_bdry_ijk)) then 
     deallocate(roller_bdry_ijk,roller_bdry_ispec,roller_bdry_normal)
  endif

  ! sources
  if (NSOURCES > 0) then
    deallocate(islice_selected_source,ispec_selected_source)
    deallocate(Mxx,Myy,Mzz,Mxy,Mxz,Myz)
    deallocate(xi_source,eta_source,gamma_source)
    deallocate(tshift_src,hdur,hdur_Gaussian)
    deallocate(utm_x_source,utm_y_source)
    deallocate(nu_source)
    deallocate(user_source_time_function)
    if (allocated(is_POINTFORCE)) deallocate(is_POINTFORCE)
    if (allocated(isource_glob2loc)) deallocate(isource_glob2loc)
  endif

  ! receivers
  if (nrec > 0) then
    deallocate(islice_selected_rec,ispec_selected_rec)
    deallocate(xi_receiver,eta_receiver,gamma_receiver)
    deallocate(station_name,network_name)
    deallocate(nu_rec)
  endif

  if (allocated(pm1_source_encoding)) deallocate(pm1_source_encoding)
  if (allocated(sourcearrays)) deallocate(sourcearrays)
  if (allocated(source_adjoint)) deallocate(source_adjoint)

  ! receiver arrays
  if (allocated(number_receiver_global)) deallocate(number_receiver_global)
  if (allocated(hxir_store)) deallocate(hxir_store,hetar_store,hgammar_store)
  if (allocated(hpxir_store)) deallocate(hpxir_store,hpetar_store,hpgammar_store)

  ! adjoint sources
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    if (nadj_rec_local > 0) then
      if (SIMULATION_TYPE == 2) then
        deallocate(number_adjsources_global)
        deallocate(hxir_adjstore,hetar_adjstore,hgammar_adjstore)
      else
        nullify(number_adjsources_global)
        nullify(hxir_adjstore,hetar_adjstore,hgammar_adjstore)
      endif
    endif
  endif

  ! seismograms
  if (allocated(seismograms_d)) deallocate(seismograms_d,seismograms_v,seismograms_a,seismograms_p)
  if (allocated(seismograms_eps)) deallocate(seismograms_eps)

  ! moment tensor derivatives
  if (allocated(Mxx_der)) deallocate(Mxx_der,Myy_der,Mzz_der,Mxy_der,Mxz_der,Myz_der,sloc_der)

  ! mesh
  deallocate(ibool)
  deallocate(irregular_element_number)
  deallocate(xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore,jacobianstore)
  deallocate(deriv_mapping)
  deallocate(xstore,ystore,zstore)
  deallocate(kappastore,mustore,rhostore)
  deallocate(ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic)

  if (ELASTIC_SIMULATION) then
    if (allocated(rho_vp)) deallocate(rho_vp)
    if (allocated(rho_vs)) deallocate(rho_vs)
  endif

  ! kernels
  if (SIMULATION_TYPE == 3) then
    ! acoustic
    if (allocated(rho_ac_kl)) deallocate(rho_ac_kl)
    if (allocated(kappa_ac_kl)) deallocate(kappa_ac_kl)
    if (allocated(rhop_ac_kl)) deallocate(rhop_ac_kl)
    if (allocated(alpha_ac_kl)) deallocate(alpha_ac_kl)
    ! elastic
    if (allocated(rho_kl)) deallocate(rho_kl)
    if (allocated(mu_kl)) deallocate(mu_kl)
    if (allocated(kappa_kl)) deallocate(kappa_kl)
    if (allocated(cijkl_kl)) deallocate(cijkl_kl)
    ! poroelastic
    if (allocated(rhot_kl)) deallocate(rhot_kl)
    if (allocated(rhof_kl)) deallocate(rhof_kl)
    if (allocated(sm_kl)) deallocate(sm_kl)
    if (allocated(eta_kl)) deallocate(eta_kl)
    if (allocated(mufr_kl)) deallocate(mufr_kl)
    if (allocated(B_kl)) deallocate(B_kl)
    if (allocated(C_kl)) deallocate(C_kl)
    if (allocated(M_kl)) deallocate(M_kl)
    if (allocated(rhob_kl)) deallocate(rhob_kl)
    if (allocated(rhofb_kl)) deallocate(rhofb_kl)
    if (allocated(Bb_kl)) deallocate(Bb_kl)
    if (allocated(Cb_kl)) deallocate(Cb_kl)
    if (allocated(Mb_kl)) deallocate(Mb_kl)
    if (allocated(mufrb_kl)) deallocate(mufrb_kl)
    if (allocated(phi_kl)) deallocate(phi_kl)
    if (allocated(rhobb_kl)) deallocate(rhobb_kl)
    if (allocated(rhofbb_kl)) deallocate(rhofbb_kl)
    if (allocated(phib_kl)) deallocate(phib_kl)
    if (allocated(ratio_kl)) deallocate(ratio_kl)
    if (allocated(cs_kl)) deallocate(cs_kl)
    if (allocated(cpI_kl)) deallocate(cpI_kl)
    if (allocated(cpII_kl)) deallocate(cpII_kl)
    ! Hessians
    if (allocated(hess_ac_kl)) deallocate(hess_ac_kl)
    if (allocated(hess_rho_ac_kl)) deallocate(hess_rho_ac_kl)
    if (allocated(hess_kappa_ac_kl)) deallocate(hess_kappa_ac_kl)
    if (allocated(hess_kl)) deallocate(hess_kl)
    if (allocated(hess_rho_kl)) deallocate(hess_rho_kl)
    if (allocated(hess_kappa_kl)) deallocate(hess_kappa_kl)
    if (allocated(hess_mu_kl)) deallocate(hess_mu_kl)
  endif

  ! faults
  if (SIMULATION_TYPE_DYN) call BC_DYNFLT_free()
  if (SIMULATION_TYPE_KIN) call BC_KINFLT_free()

  ! LTS
  if (LTS_MODE) then
    deallocate(p_level,p_level_loops,p_level_steps)
    deallocate(iglob_p_refine,ispec_p_refine)
    deallocate(p_elem,boundary_elem)
  endif

  end subroutine finalize_simulation_cleanup
