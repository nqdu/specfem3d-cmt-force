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

  subroutine save_databases(nspec,nglob, &
                            iMPIcut_xi,iMPIcut_eta, &
                            nodes_coords,ispec_material_id, &
                            nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                            ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top)

  use constants, only: MAX_STRING_LEN,IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC,IDOMAIN_POROELASTIC, &
    NDIM,IMAIN,IIN_DB,myrank

  use constants_meshfem, only: NGLLX_M,NGLLY_M,NGLLZ_M

  use shared_parameters, only: NGNOD,NGNOD2D

  use meshfem_par, only: ibool, &
    addressing,NPROC_XI,NPROC_ETA,iproc_xi_current,iproc_eta_current, &
    prname, &
    NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
    NMATERIALS,material_properties,material_properties_undef, &
    nspec_CPML,is_CPML,CPML_to_spec,CPML_regions

  implicit none

  ! number of spectral elements in each block
  integer, intent(in) :: nspec

  ! number of vertices in each block
  integer, intent(in) :: nglob

  logical, intent(in) :: iMPIcut_xi(2,nspec),iMPIcut_eta(2,nspec)

  ! arrays with the mesh
  double precision, intent(in) :: nodes_coords(nglob,NDIM)
  integer, intent(in) :: ispec_material_id(nspec)

  ! boundary parameters locator
  integer, intent(in) :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax
  integer, intent(in) :: ibelm_xmin(NSPEC2DMAX_XMIN_XMAX),ibelm_xmax(NSPEC2DMAX_XMIN_XMAX)
  integer, intent(in) :: ibelm_ymin(NSPEC2DMAX_YMIN_YMAX),ibelm_ymax(NSPEC2DMAX_YMIN_YMAX)
  integer, intent(in) :: ibelm_bottom(NSPEC2D_BOTTOM)
  integer, intent(in) :: ibelm_top(NSPEC2D_TOP)

  ! local parameters
  ! MPI Cartesian topology
  ! W for West (= XI_MIN), E for East (= XI_MAX), S for South (= ETA_MIN), N for North (= ETA_MAX)
  integer, parameter :: W = 1,E = 2,S = 3,N = 4,NW = 5,NE = 6,SE = 7,SW = 8

  ! CPML
  integer :: nspec_CPML_total,ispec_CPML
  double precision , dimension(17) :: matpropl
  integer :: i,ispec,iglob,ier

  ! for MPI interfaces
  integer ::  nb_interfaces,nspec_interfaces_max
  logical, dimension(8) ::  interfaces
  integer, dimension(8) ::  nspec_interface

  integer, parameter :: IIN_database = IIN_DB

  integer :: ndef,nundef
  integer :: mat_id,domain_id
  ! there was a potential bug here if nspec is big
  !integer,dimension(2,nspec) :: material_index
  integer,dimension(:,:),allocatable :: material_index
  character(len=MAX_STRING_LEN), dimension(6,1) :: undef_mat_prop

  ! temporary array for local nodes (either NGNOD2D or NGNOD)
  integer, dimension(NGNOD) :: loc_node
  integer, dimension(NGNOD) :: anchor_iax,anchor_iay,anchor_iaz
  integer :: ia,inode

  ! assignes material index
  ! format: (1,ispec) = #material_id , (2,ispec) = #material_definition
  allocate(material_index(2,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1346')
  if (ier /= 0) stop 'Error allocating array material_index'
  material_index (:,:) = 0
  do ispec = 1, nspec
    ! material id
    material_index(1,ispec) = ispec_material_id(ispec)
    ! material definition: 1 = interface type / 2 = tomography type
    if (ispec_material_id(ispec) > 0) then
      ! dummy value, not used any further
      material_index(2,ispec) = 1
    else
      ! negative material ids
      ! by default, assumes tomography model material
      ! (note: interface type not implemented yet...)
      material_index(2,ispec) = 2
    endif
  enddo

  ! Materials properties
  ! counts defined/undefined materials
  ndef = 0
  nundef = 0
  do i = 1,NMATERIALS
    ! material_properties(:,:) array:
    !   first dimension  : material_id
    !   second dimension : #rho  #vp  #vs  #Q_Kappa  #Q_mu  #anisotropy_flag  #domain_id  #material_id
    !
    ! material properties format: #rho  #vp  #vs  #Q_Kappa  #Q_mu  #anisotropy_flag  #domain_id  #material_id
    mat_id = int(material_properties(i,8))
    if (mat_id > 0) ndef = ndef + 1
    if (mat_id < 0) nundef = nundef + 1
  enddo
  !debug
  !print *,'materials def/undef: ',ndef,nundef

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'mesh files:'
    write(IMAIN,*) '  saving files: proc***_Database'
    call flush_IMAIN()
  endif

  ! opens database file
  open(unit=IIN_database,file=prname(1:len_trim(prname))//'Database', &
       status='unknown',action='write',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening Database file: ',prname(1:len_trim(prname))//'Database'
    stop 'error opening Database file'
  endif

  ! global nodes
  write(IIN_database) nglob
  do iglob = 1,nglob
    ! format: #id #x #y #z
    write(IIN_database) iglob,nodes_coords(iglob,1),nodes_coords(iglob,2),nodes_coords(iglob,3)
  enddo

  ! materials
  ! format: #number of defined materials #number of undefined materials
  write(IIN_database) ndef, nundef

  ! writes out defined materials
  do i = 1,NMATERIALS
    ! material properties format: #rho  #vp  #vs  #Q_Kappa  #Q_mu  #anisotropy_flag  #domain_id  #material_id
    domain_id = int(material_properties(i,7))
    mat_id = int(material_properties(i,8))
    if (mat_id > 0) then
      ! pad dummy zeros to fill up 17 entries
      matpropl(:) = 0.d0
      select case(domain_id)
      case (IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC)
        ! material properties format:
        !#(1)rho  #(2)vp #(3)vs #(4)Q_Kappa #(5)Q_mu #(6)anisotropy_flag #(7)domain_id #(8)mat_id
        !
        ! output format for xgenerate_database (same as for cubit/trelis inputs):
        !   rho,vp,vs,Q_Kappa,Q_mu,anisotropy_flag,material_domain_id
        !
        ! skipping mat_id, not needed
        matpropl(1:7) = material_properties(i,1:7)
      case (IDOMAIN_POROELASTIC)
        ! material properties format:
        !#(1)rho_s #(2)rho_f #(3)phi #(4)tort #(5)eta #(6)0 #(7)domain_id #(8)mat_id
        !            .. #(9)kxx #(10)kxy #(11)kxz #(12)kyy #(13)kyz #(14)kzz #(15)kappa_s #(16)kappa_f #(17)kappa_fr #(18)mu_fr
        !
        ! output format for xgenerate_database (same as for cubit/trelis inputs):
        !   rhos,rhof,phi,tort,eta,0,material_domain_id,kxx,kxy,kxz,kyy,kyz,kzz,kappas,kappaf,kappafr,mufr
        matpropl(1:7) = material_properties(i,1:7)
        ! skipping mat_id, not needed
        matpropl(8:17) = material_properties(i,9:18)
      end select
      ! writes to database
      write(IIN_database) matpropl(:)
    endif
  enddo

  ! writes out undefined materials
  do i = 1,NMATERIALS
    domain_id = int(material_properties(i,7))
    mat_id = int(material_properties(i,8))
    if (mat_id < 0) then
      ! format:
      ! #material_id #type-keyword #domain-name #tomo-filename #tomo_id #domain_id
      ! format example tomography: -1 tomography elastic tomography_model.xyz 0 2
      undef_mat_prop(:,:) = ''
      ! material id
      write(undef_mat_prop(1,1),*) mat_id
      ! keyword/domain/filename
      undef_mat_prop(2,1) = material_properties_undef(i,1)  ! type-keyword : tomography/interface
      undef_mat_prop(3,1) = material_properties_undef(i,2)  ! domain-name  : acoustic/elastic/poroelastic
      undef_mat_prop(4,1) = material_properties_undef(i,3)  ! tomo-filename: tomography_model**.xyz
      ! checks consistency between domain-name and domain_id
      select case (domain_id)
      case (IDOMAIN_ACOUSTIC)
        if (trim(undef_mat_prop(3,1)) /= 'acoustic') stop 'Error in undef_mat_prop acoustic domain'
      case (IDOMAIN_ELASTIC)
        if (trim(undef_mat_prop(3,1)) /= 'elastic')  stop 'Error in undef_mat_prop elastic domain'
      case (IDOMAIN_POROELASTIC)
        if (trim(undef_mat_prop(3,1)) /= 'poroelastic')  stop 'Error in undef_mat_prop poroelastic domain'
      end select
      ! default name if none given
      if (trim(undef_mat_prop(4,1)) == "") undef_mat_prop(4,1) = 'tomography_model.xyz'
      ! default tomo-id (unused)
      write(undef_mat_prop(5,1),*) 0
      ! domain-id
      write(undef_mat_prop(6,1),*) domain_id
      ! debug
      !print *,'undef mat: ',undef_mat_prop
      ! writes out properties
      write(IIN_database) undef_mat_prop
    endif
  enddo

  ! spectral elements
  !
  ! note: check with routine write_partition_database() to produce identical output
  write(IIN_database) nspec

  ! sets up node addressing
  call hex_nodes_anchor_ijk_NGLL(NGNOD,anchor_iax,anchor_iay,anchor_iaz,NGLLX_M,NGLLY_M,NGLLZ_M)

  do ispec = 1,nspec
    ! format: #ispec #material_id #dummy/tomo_id #iglob1 #iglob2 ..

    ! assumes NGLLX_M == NGLLY_M == NGLLZ_M == 2
    !write(IIN_database) ispec,material_index(1,ispec),material_index(2,ispec), &
    !       ibool(1,1,1,ispec),ibool(2,1,1,ispec),ibool(2,2,1,ispec),ibool(1,2,1,ispec),ibool(1,1,2,ispec), &
    !       ibool(2,1,2,ispec),ibool(2,2,2,ispec),ibool(1,2,2,ispec)

    ! gets anchor nodes
    do ia = 1,NGNOD
      iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
      loc_node(ia) = iglob
    enddo
    ! output
    write(IIN_database) ispec,material_index(1,ispec),material_index(2,ispec),(loc_node(ia),ia = 1,NGNOD)
  enddo


  ! Boundaries
  !
  ! note: check with routine write_boundaries_database() to produce identical output
  write(IIN_database) 1,nspec2D_xmin
  write(IIN_database) 2,nspec2D_xmax
  write(IIN_database) 3,nspec2D_ymin
  write(IIN_database) 4,nspec2D_ymax
  write(IIN_database) 5,NSPEC2D_BOTTOM
  write(IIN_database) 6,NSPEC2D_TOP

  do i = 1,nspec2D_xmin
    ispec = ibelm_xmin(i)
    ! gets anchor nodes on xmin
    inode = 0
    do ia = 1,NGNOD
      if (anchor_iax(ia) == 1) then
        iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
        inode = inode + 1
        if (inode > NGNOD2D) stop 'inode index exceeds NGNOD2D for xmin'
        loc_node(inode) = iglob
      endif
    enddo
    if (inode /= NGNOD2D) stop 'Invalid number of inodes found for xmin'
    write(IIN_database) ispec,(loc_node(inode), inode = 1,NGNOD2D)

    ! assumes NGLLX_M == NGLLY_M == NGLLZ_M == 2
    !write(IIN_database) ispec,ibool(1,1,1,ispec),ibool(1,NGLLY_M,1,ispec), &
    !                          ibool(1,1,NGLLZ_M,ispec),ibool(1,NGLLY_M,NGLLZ_M,ispec)
  enddo

  do i = 1,nspec2D_xmax
    ispec = ibelm_xmax(i)
    ! gets anchor nodes on xmax
    inode = 0
    do ia = 1,NGNOD
      if (anchor_iax(ia) == NGLLX_M) then
        iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
        inode = inode + 1
        if (inode > NGNOD2D) stop 'inode index exceeds NGNOD2D for xmax'
        loc_node(inode) = iglob
      endif
    enddo
    if (inode /= NGNOD2D) stop 'Invalid number of inodes found for xmax'
    write(IIN_database) ispec,(loc_node(inode), inode = 1,NGNOD2D)

    ! write(IIN_database) ibelm_xmax(i),ibool(NGLLX_M,1,1,ibelm_xmax(i)),ibool(NGLLX_M,NGLLY_M,1,ibelm_xmax(i)), &
    !      ibool(NGLLX_M,1,NGLLZ_M,ibelm_xmax(i)),ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_xmax(i))
  enddo

  do i = 1,nspec2D_ymin
    ispec = ibelm_ymin(i)
    ! gets anchor nodes on ymin
    inode = 0
    do ia = 1,NGNOD
      if (anchor_iay(ia) == 1) then
        iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
        inode = inode + 1
        if (inode > NGNOD2D) stop 'inode index exceeds NGNOD2D for ymin'
        loc_node(inode) = iglob
      endif
    enddo
    if (inode /= NGNOD2D) stop 'Invalid number of inodes found for ymin'
    write(IIN_database) ispec,(loc_node(inode), inode = 1,NGNOD2D)

    ! write(IIN_database) ibelm_ymin(i),ibool(1,1,1,ibelm_ymin(i)),ibool(NGLLX_M,1,1,ibelm_ymin(i)), &
    !      ibool(1,1,NGLLZ_M,ibelm_ymin(i)),ibool(NGLLX_M,1,NGLLZ_M,ibelm_ymin(i))
  enddo

  do i = 1,nspec2D_ymax
    ispec = ibelm_ymax(i)
    ! gets anchor nodes on ymax
    inode = 0
    do ia = 1,NGNOD
      if (anchor_iay(ia) == NGLLY_M) then
        iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
        inode = inode + 1
        if (inode > NGNOD2D) stop 'inode index exceeds NGNOD2D for ymax'
        loc_node(inode) = iglob
      endif
    enddo
    if (inode /= NGNOD2D) stop 'Invalid number of inodes found for ymax'
    write(IIN_database) ispec,(loc_node(inode), inode = 1,NGNOD2D)

    ! write(IIN_database) ibelm_ymax(i),ibool(NGLLX_M,NGLLY_M,1,ibelm_ymax(i)),ibool(1,NGLLY_M,1,ibelm_ymax(i)), &
    !      ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_ymax(i)),ibool(1,NGLLY_M,NGLLZ_M,ibelm_ymax(i))
  enddo

  do i = 1,NSPEC2D_BOTTOM
    ispec = ibelm_bottom(i)
    ! gets anchor nodes on bottom
    inode = 0
    do ia = 1,NGNOD
      if (anchor_iaz(ia) == 1) then
        iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
        inode = inode + 1
        if (inode > NGNOD2D) stop 'inode index exceeds NGNOD2D for bottom'
        loc_node(inode) = iglob
      endif
    enddo
    if (inode /= NGNOD2D) stop 'Invalid number of inodes found for bottom'
    write(IIN_database) ispec,(loc_node(inode), inode = 1,NGNOD2D)

    ! write(IIN_database) ibelm_bottom(i),ibool(1,1,1,ibelm_bottom(i)),ibool(NGLLX_M,1,1,ibelm_bottom(i)), &
    !      ibool(NGLLX_M,NGLLY_M,1,ibelm_bottom(i)),ibool(1,NGLLY_M,1,ibelm_bottom(i))
  enddo

  do i = 1,NSPEC2D_TOP
    ispec = ibelm_top(i)
    ! gets anchor nodes on top
    inode = 0
    do ia = 1,NGNOD
      if (anchor_iaz(ia) == NGLLZ_M) then
        iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
        inode = inode + 1
        if (inode > NGNOD2D) stop 'inode index exceeds NGNOD2D for top'
        loc_node(inode) = iglob
      endif
    enddo
    if (inode /= NGNOD2D) stop 'Invalid number of inodes found for top'
    write(IIN_database) ispec,(loc_node(inode), inode = 1,NGNOD2D)

    ! write(IIN_database) ibelm_top(i),ibool(1,1,NGLLZ_M,ibelm_top(i)),ibool(NGLLX_M,1,NGLLZ_M,ibelm_top(i)), &
    !      ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_top(i)),ibool(1,NGLLY_M,NGLLZ_M,ibelm_top(i))
  enddo

  ! compatibility for fixed and roller boundaries (note: not used in xgenerate_database, this option will be enabled in SAVE_CUBIT_MESH = .true.
  write(IIN_database) 8,0
  write(IIN_database) 8,0

  ! CPML
  !
  ! note: check with routine write_cpml_database() to produce identical output
  call sum_all_i(nspec_CPML,nspec_CPML_total)
  call synchronize_all()
  call bcast_all_singlei(nspec_CPML_total)
  call synchronize_all()

  write(IIN_database) nspec_CPML_total

  if (nspec_CPML_total > 0) then
     write(IIN_database) nspec_CPML

     do ispec_CPML = 1,nspec_CPML
        write(IIN_database) CPML_to_spec(ispec_CPML), CPML_regions(ispec_CPML)
     enddo
     do ispec = 1,nspec
        write(IIN_database) is_CPML(ispec)
     enddo
  endif

  ! MPI Interfaces
  !
  ! note: check with routine write_interfaces_database() to produce identical output
  if (NPROC_XI > 1 .or. NPROC_ETA > 1) then
    ! determines number of MPI interfaces for each slice
    nb_interfaces = 4
    interfaces(W:N) = .true.
    interfaces(NW:SW) = .false.

    ! slices at model boundaries
    if (iproc_xi_current == 0) then
      nb_interfaces =  nb_interfaces -1
      interfaces(W) = .false.   ! XI_min
    endif
    if (iproc_xi_current == NPROC_XI-1) then
      nb_interfaces =  nb_interfaces -1
      interfaces(E) = .false.   ! XI_max
    endif
    if (iproc_eta_current == 0) then
      nb_interfaces =  nb_interfaces -1
      interfaces(S) = .false.   ! ETA_min
    endif
    if (iproc_eta_current == NPROC_ETA-1) then
      nb_interfaces =  nb_interfaces -1
      interfaces(N) = .false.   ! ETA_max
    endif

    ! slices in middle of model
    if ((interfaces(N) .eqv. .true.) .and. (interfaces(W) .eqv. .true.)) then
      interfaces(NW) = .true.
      nb_interfaces =  nb_interfaces +1
    endif
    if ((interfaces(N) .eqv. .true.) .and. (interfaces(E) .eqv. .true.)) then
      interfaces(NE) = .true.
      nb_interfaces =  nb_interfaces +1
    endif
    if ((interfaces(S) .eqv. .true.) .and. (interfaces(E) .eqv. .true.)) then
      interfaces(SE) = .true.
      nb_interfaces =  nb_interfaces +1
    endif
    if ((interfaces(S) .eqv. .true.) .and. (interfaces(W) .eqv. .true.)) then
      interfaces(SW) = .true.
      nb_interfaces =  nb_interfaces +1
    endif

    nspec_interface(:) = 0
    if (interfaces(W))  nspec_interface(W) = count(iMPIcut_xi(1,:) .eqv. .true.)    ! XI_min
    if (interfaces(E))  nspec_interface(E) = count(iMPIcut_xi(2,:) .eqv. .true.)    ! XI_max
    if (interfaces(S))  nspec_interface(S) = count(iMPIcut_eta(1,:) .eqv. .true.)   ! ETA_min
    if (interfaces(N))  nspec_interface(N) = count(iMPIcut_eta(2,:) .eqv. .true.)   ! ETA_max
    if (interfaces(NW))  nspec_interface(NW) = count((iMPIcut_xi(1,:) .eqv. .true.) .and. (iMPIcut_eta(2,:) .eqv. .true.))
    if (interfaces(NE))  nspec_interface(NE) = count((iMPIcut_xi(2,:) .eqv. .true.) .and. (iMPIcut_eta(2,:) .eqv. .true.))
    if (interfaces(SE))  nspec_interface(SE) = count((iMPIcut_xi(2,:) .eqv. .true.) .and. (iMPIcut_eta(1,:) .eqv. .true.))
    if (interfaces(SW))  nspec_interface(SW) = count((iMPIcut_xi(1,:) .eqv. .true.) .and. (iMPIcut_eta(1,:) .eqv. .true.))

    nspec_interfaces_max = maxval(nspec_interface)

    ! format: #number_of_MPI_interfaces  #maximum_number_of_elements_on_each_interface
    write(IIN_database) nb_interfaces,nspec_interfaces_max

    ! face elements
    if (interfaces(W)) then
      ! format: #process_interface_id  #number_of_elements_on_interface
      write(IIN_database) addressing(iproc_xi_current-1,iproc_eta_current),nspec_interface(W)
      do ispec = 1,nspec
        if (iMPIcut_xi(1,ispec)) then
          ! format: #(1)spectral_element_id  #(2)interface_type  #(3)node_id1  #(4)node_id2 #(5).. #(6)..
          ! note: face outputs 4 corner points
          write(IIN_database) ispec,4,ibool(1,1,1,ispec),ibool(1,NGLLY_M,1,ispec), &
                                      ibool(1,1,NGLLZ_M,ispec),ibool(1,NGLLY_M,NGLLZ_M,ispec)
        endif
      enddo
    endif

    if (interfaces(E)) then
      ! format: #process_interface_id  #number_of_elements_on_interface
      write(IIN_database) addressing(iproc_xi_current+1,iproc_eta_current),nspec_interface(E)
      do ispec = 1,nspec
        if (iMPIcut_xi(2,ispec)) then
          ! format: #(1)spectral_element_id  #(2)interface_type  #(3)node_id1  #(4)node_id2 #(5).. #(6)..
          write(IIN_database) ispec,4,ibool(NGLLX_M,1,1,ispec),ibool(NGLLX_M,NGLLY_M,1,ispec), &
                                      ibool(NGLLX_M,1,NGLLZ_M,ispec),ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
        endif
      enddo
    endif

    if (interfaces(S)) then
      ! format: #process_interface_id  #number_of_elements_on_interface
      write(IIN_database) addressing(iproc_xi_current,iproc_eta_current-1),nspec_interface(S)
      do ispec = 1,nspec
        if (iMPIcut_eta(1,ispec)) then
          ! format: #(1)spectral_element_id  #(2)interface_type  #(3)node_id1  #(4)node_id2 #(5).. #(6)..
          write(IIN_database) ispec,4,ibool(1,1,1,ispec),ibool(NGLLX_M,1,1,ispec), &
                                      ibool(1,1,NGLLZ_M,ispec),ibool(NGLLX_M,1,NGLLZ_M,ispec)
        endif
      enddo
    endif

    if (interfaces(N)) then
      ! format: #process_interface_id  #number_of_elements_on_interface
      write(IIN_database) addressing(iproc_xi_current,iproc_eta_current+1),nspec_interface(N)
      do ispec = 1,nspec
        if (iMPIcut_eta(2,ispec)) then
          ! format: #(1)spectral_element_id  #(2)interface_type  #(3)node_id1  #(4)node_id2 #(5).. #(6)..
          write(IIN_database) ispec,4,ibool(NGLLX_M,NGLLY_M,1,ispec),ibool(1,NGLLY_M,1,ispec), &
                                      ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ispec),ibool(1,NGLLY_M,NGLLZ_M,ispec)
        endif
      enddo
    endif

    ! edge elements
    if (interfaces(NW)) then
      ! format: #process_interface_id  #number_of_elements_on_interface
      write(IIN_database) addressing(iproc_xi_current-1,iproc_eta_current+1),nspec_interface(NW)
      do ispec = 1,nspec
        if ((iMPIcut_xi(1,ispec) .eqv. .true.) .and. (iMPIcut_eta(2,ispec) .eqv. .true.)) then
          ! note: edge elements output 2 corners and 2 dummy values
          ! format: #(1)spectral_element_id  #(2)interface_type  #(3)node_id1  #(4)node_id2 #(5).. #(6)..
          write(IIN_database) ispec,2,ibool(1,NGLLY_M,1,ispec),ibool(1,NGLLY_M,NGLLZ_M,ispec),-1,-1
        endif
      enddo
    endif

    if (interfaces(NE)) then
      ! format: #process_interface_id  #number_of_elements_on_interface
      write(IIN_database) addressing(iproc_xi_current+1,iproc_eta_current+1),nspec_interface(NE)
      do ispec = 1,nspec
        if ((iMPIcut_xi(2,ispec) .eqv. .true.) .and. (iMPIcut_eta(2,ispec) .eqv. .true.)) then
          ! format: #(1)spectral_element_id  #(2)interface_type  #(3)node_id1  #(4)node_id2 #(5).. #(6)..
          write(IIN_database) ispec,2,ibool(NGLLX_M,NGLLY_M,1,ispec),ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ispec),-1,-1
        endif
      enddo
    endif

    if (interfaces(SE)) then
      ! format: #process_interface_id  #number_of_elements_on_interface
      write(IIN_database) addressing(iproc_xi_current+1,iproc_eta_current-1),nspec_interface(SE)
      do ispec = 1,nspec
        if ((iMPIcut_xi(2,ispec) .eqv. .true.) .and. (iMPIcut_eta(1,ispec) .eqv. .true.)) then
          ! format: #(1)spectral_element_id  #(2)interface_type  #(3)node_id1  #(4)node_id2 #(5).. #(6)..
          write(IIN_database) ispec,2,ibool(NGLLX_M,1,1,ispec),ibool(NGLLX_M,1,NGLLZ_M,ispec),-1,-1
        endif
      enddo
    endif

    if (interfaces(SW)) then
      ! format: #process_interface_id  #number_of_elements_on_interface
      write(IIN_database) addressing(iproc_xi_current-1,iproc_eta_current-1),nspec_interface(SW)
      do ispec = 1,nspec
        if ((iMPIcut_xi(1,ispec) .eqv. .true.) .and. (iMPIcut_eta(1,ispec) .eqv. .true.)) then
          ! format: #(1)spectral_element_id  #(2)interface_type  #(3)node_id1  #(4)node_id2 #(5).. #(6)..
          write(IIN_database) ispec,2,ibool(1,1,1,ispec),ibool(1,1,NGLLZ_M,ispec),-1,-1
        endif
      enddo
    endif
  else
    ! single process execution, no MPI boundaries
    nb_interfaces = 0
    nspec_interfaces_max = 0
    ! format: #number_of_MPI_interfaces  #maximum_number_of_elements_on_each_interface
    write(IIN_database) nb_interfaces,nspec_interfaces_max
  endif

  close(IIN_database)

  deallocate(material_index)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  done mesh files'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine save_databases

