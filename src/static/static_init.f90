submodule (static_module) static_init 

contains 

  subroutine initialize_static_module
    use specfem_par, only: NDIM,NGLLX,NGLLY,NGLLZ 
    use specfem_par, only: nglob => NGLOB_AB, nspec => NSPEC_AB
    implicit none

    ! allocate arrays
    allocate(inv_mult_(nglob))
    allocate(ssol%force_ext(NDIM,nglob))
    allocate(ssol%stress(NGLLX,NGLLY,NGLLZ,nspec,6))
    allocate(ssol%strain(NGLLX,NGLLY,NGLLZ,nspec,6))

    ! read parameters for static solver
    call read_params_static()

    ! create PETSc backend if enabled
    call ssol%init_petsc()

    ! get external force from file and interpolate to global points, store in force_ext
    call read_ext_force()

    ! set inv_mat, used for parallel inner product
    call prepare_inv_mult()

  end subroutine initialize_static_module

  subroutine read_params_static()
    use constants, only: myrank,MAX_STRING_LEN
    implicit none
    integer :: ier

    if(myrank == 0) then 
      call open_parameter_file(ier)
      if (ier /= 0) then
        print*, "Error opening parameter file for static module initialization"
        stop 
      endif

      ! read static solver parameters
      ssol%is_nonlinear = .false. ! default value
      !call read_value_logical(ssol%SAVE_STATIC_FIELD, "SAVE_STATIC_FIELD",ier)
      ssol%SAVE_STATIC_FIELD = .true. ! default value

      ! read USE_PETSC_AS_BACKEND parameter
      call read_value_logical(ssol%USE_PETSC_AS_BACKEND, "USE_PETSC_AS_BACKEND",ier)
      if (ier /= 0) then
        !print*, "no USE_PETSC_AS_BACKEND specified, default to true"
        ssol%USE_PETSC_AS_BACKEND = .false.
      endif
      print*, "USE_PETSC_AS_BACKEND = ", ssol%USE_PETSC_AS_BACKEND


      ! close
      call close_parameter_file()
    endif

    ! broadcast parameters to all ranks
    call bcast_all_singlel(ssol%SAVE_STATIC_FIELD)
    call bcast_all_singlel(ssol%is_nonlinear)
    call bcast_all_singlel(ssol%USE_PETSC_AS_BACKEND)

  end subroutine read_params_static

  subroutine read_ext_force()
    use constants, only: myrank,MAX_STRING_LEN
    use specfem_par, only: nspec => NSPEC_AB
    use specfem_par, only: NDIM,NGLLX,NGLLY,NGLLZ,wxgll,wygll,wzgll 
    use specfem_par, only: ibool,xstore,ystore,zstore,jacobianstore,ROTATION
    use specfem_par, only: NPROC,NGLOB_AB,rhostore,&
                          buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                          num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                          nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                          my_neighbors_ext_mesh,&
                          request_send_vector_ext_mesh,request_recv_vector_ext_mesh
    use specfem_par, only: fixed_bdry_ijk, fixed_bdry_ispec,num_fixed_bdry_faces
    use specfem_par, only: num_roller_bdry_faces,roller_bdry_ijk, &
                            roller_bdry_ispec, roller_bdry_normal
    use petsc_interfaces, only: fill_vec_petsc
    implicit none
    
    integer :: i,j,k,iglob,ier,ispec,iface,igll2 
    integer :: nx,ny,nz
    real(kind=dp) :: xmin,xmax,ymin,ymax,zmin,zmax 
    character(len=MAX_STRING_LEN) :: force_filename 
    integer,parameter :: IO_UNIT = 10
    real(kind=dp),allocatable :: tomo_x(:),tomo_y(:),tomo_z(:)
    real(kind=dp),allocatable :: tomo_force(:,:,:,:)
    real(kind=dp) :: fp(NDIM), temp, normal_vec(NDIM)
    real(kind=CUSTOM_REAL) :: omega(NDIM), rot_org(NDIM),tempx,tempy,tempz 
    real(kind=CUSTOM_REAL) :: xp, yp,zp

    ! element wise force terms
    real(kind=dp), allocatable :: elem_force(:,:,:,:,:) ! shape(NDIM,NGLLX,NGLLY,NGLLZ,nspec)

    ! read and interpolate force_ext from file
    if(myrank == 0) then 
      ! open Par_file
      call open_parameter_file(ier)
      if (ier /= 0) then
        print*, "Error opening parameter file for static module initialization"
        stop 
      endif

      ! read whether the model has a rotation system
      call read_value_string(force_filename, "STATIC_FORCE_FILE",ier)
      if (ier /= 0) then
        print*, "no STATIC_FORCE_FILE specified, assuming zero external force"
        nx = 2; ny = 2; nz = 2
        xmin = 0.0_dp; xmax = 1.0_dp
        ymin = 0.0_dp; ymax = 1.0_dp
        zmin = 0.0_dp; zmax = 1.0_dp

        allocate(tomo_x(nx), tomo_y(ny), tomo_z(nz))
        allocate(tomo_force(nx,ny,nz,NDIM))
        tomo_force(:,:,:,:) = 0.0_dp
      else
        print*, "Reading static force from file: ", trim(force_filename)
        open(unit=IO_UNIT, file=trim(force_filename), status='old', action='read', iostat=ier)
        if (ier /= 0) then
          print*, "Error opening static force file: ", trim(force_filename)
          stop
        endif

        ! read dimensions
        read(IO_UNIT,*) nx,  ny, nz
        read(IO_UNIT,*) xmin, xmax, ymin, ymax, zmin, zmax

        ! allocate arrays
        allocate(tomo_x(nx), tomo_y(ny), tomo_z(nz))
        allocate(tomo_force(nx,ny,nz,NDIM))

        ! read data
        do k = 1,nz; do j = 1,ny; do i = 1,nx
          read(IO_UNIT,*) tomo_force(i,j,k,:)
        enddo; enddo; enddo

        close(IO_UNIT)
      endif 

      ! set tomo_x/y/z based on bounds and dimensions
      do i = 1, nx
        tomo_x(i) = xmin + (i-1)*(xmax-xmin)/(nx-1)
      enddo
      do j = 1, ny
        tomo_y(j) = ymin + (j-1)*(ymax-ymin)/(ny-1)
      enddo
      do k = 1, nz
        tomo_z(k) = zmin + (k-1)*(zmax-zmin)/(nz-1)
      enddo

      ! close Par_file
      call close_parameter_file()
    endif

    ! broadcast dimensions and bounds to all ranks
    call bcast_all_singlei(nx)
    call bcast_all_singlei(ny)
    call bcast_all_singlei(nz)
    call bcast_all_singledp(xmin)
    call bcast_all_singledp(xmax)
    call bcast_all_singledp(ymin)
    call bcast_all_singledp(ymax)
    call bcast_all_singledp(zmin)
    call bcast_all_singledp(zmax)
    
    ! broadcast tomo_x/y/z and tomo_force to all ranks
    if(myrank /= 0) then
      allocate(tomo_x(nx), tomo_y(ny), tomo_z(nz))
      allocate(tomo_force(nx,ny,nz,NDIM))
    endif
    call synchronize_all()

    ! bcast tomo_x/y/z and tomo_force to all ranks
    call bcast_all_dp(tomo_x,size(tomo_x))
    call bcast_all_dp(tomo_y,size(tomo_y))
    call bcast_all_dp(tomo_z,size(tomo_z))
    call bcast_all_dp(tomo_force,size(tomo_force))

    ! init rotation parameters if rotation is enabled
    omega(:) = 0.0_CUSTOM_REAL
    rot_org(:) = 0.0_CUSTOM_REAL

    if(ROTATION) then
      if(myrank == 0) then
        call open_parameter_file(ier)
        call read_value_string(force_filename, "ROTATION_OMEGA",ier)
         if (ier /= 0) then
          print*, "Error reading ROTATION_OMEGA from parameter file, no rotation will be applied"
          stop 
         endif
        read(force_filename,*) omega

        call read_value_string(force_filename, "ROTATION_ORIGIN",ier)
        if (ier /= 0) then
          print*, "Error reading ROTATION_ORIGIN from parameter file, no rotation will be applied"
          stop 
        endif
        read(force_filename,*) rot_org

        ! print rotation parameters
        print*, "Applying rotation with omega = ", omega
        print*, "Rotation origin = ", rot_org

        call close_parameter_file()
      endif

      ! broadcast rotation parameters if needed
      call bcast_all_cr(omega,size(omega))
      call bcast_all_cr(rot_org,size(rot_org))
    endif

    ! set value to parameters for later use
    ssol%omega(:) = omega(:)
    ssol%rot_org(:) = rot_org(:)

    ! allocate elem_force array
    allocate(elem_force(NDIM,NGLLX,NGLLY,NGLLZ,nspec))
    elem_force(:,:,:,:,:) = 0.0_dp

    ! loop each global point and compute interpolated force, store in force_ext
    do ispec = 1, nspec
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool(i,j,k,ispec)
            temp = jacobianstore(i,j,k,ispec) * wxgll(i)  * wygll(j) * wzgll(k)
            call trilinear_interp(nx,ny,nz,tomo_x,tomo_y,tomo_z,tomo_force,&
                                  dble(xstore(iglob)),dble(ystore(iglob)),dble(zstore(iglob)),&
                                  fp)
            elem_force(:,i,j,k,ispec) = fp * temp
          enddo
        enddo
      enddo
    enddo

    ! add rotation force contribution
    do ispec = 1, nspec
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool(i,j,k,ispec)
            temp = jacobianstore(i,j,k,ispec) * wxgll(i)  * wygll(j) * wzgll(k)
            temp = rhostore(i,j,k,ispec) * temp 

            ! shift to rotation origin
            xp = xstore(iglob) - rot_org(1)
            yp = ystore(iglob) - rot_org(2)
            zp = zstore(iglob) - rot_org(3)

            ! omega x r 
            call cross_product(omega(1), omega(2), omega(3), xp, yp, zp, tempx, tempy, tempz)
            xp = tempx; yp = tempy; zp = tempz
            
            ! omega x (omega x r)
            call cross_product(omega(1), omega(2), omega(3), xp, yp, zp, tempx, tempy, tempz)

            ! note centrifugal force is - omega x (omega x r), and Coriolis force is - 2 omega x v, but we only consider centrifugal force here since it's a static problem
            elem_force(1,i,j,k,ispec) = elem_force(1,i,j,k,ispec) - real(tempx * temp, kind=CUSTOM_REAL)
            elem_force(2,i,j,k,ispec) = elem_force(2,i,j,k,ispec) - real(tempy * temp, kind=CUSTOM_REAL)
            elem_force(3,i,j,k,ispec) = elem_force(3,i,j,k,ispec) - real(tempz * temp, kind=CUSTOM_REAL)

          enddo
        enddo
      enddo
    enddo

    ! apply fixed boundary conditions 
    do iface = 1, num_fixed_bdry_faces
    !do iface = 1,0 ! do nothing --- IGNORE ---
      ispec = fixed_bdry_ispec(iface)

      do igll2 = 1, NGLLX*NGLLZ 
        i = fixed_bdry_ijk(1,igll2,iface)
        j = fixed_bdry_ijk(2,igll2,iface)
        k = fixed_bdry_ijk(3,igll2,iface)
        
        elem_force(:,i,j,k,ispec) = 0.0_dp
      enddo
    enddo

    ! apply roller boundary conditions
    do iface = 1, num_roller_bdry_faces
      ispec = roller_bdry_ispec(iface)

      do igll2 = 1, NGLLX*NGLLZ 
        i = roller_bdry_ijk(1,igll2,iface)
        j = roller_bdry_ijk(2,igll2,iface)
        k = roller_bdry_ijk(3,igll2,iface)
        
        normal_vec(:) = roller_bdry_normal(:,igll2,iface)
        fp(:) = elem_force(:,i,j,k,ispec)

        fp = fp -dot_product(fp, normal_vec) * normal_vec ! remove normal component of the force    
        elem_force(:,i,j,k,ispec) = fp
      enddo
    enddo

    ! check if PETSc assembly is enabled, if so, we need to assemble the global force vector using MPI communication, otherwise we can directly store the interpolated force in force_ext and let the Fortran solver handle the assembly
    if(ssol%USE_PETSC_AS_BACKEND) then 
      if(myrank == 0) then
        print*,'-------------------------------'
        print*, "Assembling global force vector for PETSc backend..."
        print*,'-------------------------------'
      end if
      call fill_vec_petsc(ssol%petsc_ptr, elem_force(:,:,:,:,:))
    else 
      ! directly store the interpolated force in force_ext
      ssol%force_ext(:,:) = 0.0_CUSTOM_REAL
      do ispec = 1, nspec
        do k = 1, NGLLZ
          do j = 1, NGLLY
            do i = 1, NGLLX
              iglob = ibool(i,j,k,ispec)
              ssol%force_ext(:,iglob) = ssol%force_ext(:,iglob) + real(elem_force(:,i,j,k,ispec), kind=CUSTOM_REAL)
            enddo
          enddo
        enddo
      enddo

      ! mpi sync 
      call assemble_MPI_vector_async_send(NPROC,NGLOB_AB,ssol%force_ext, &
                                          buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                                          num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                          nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                          my_neighbors_ext_mesh, &
                                          request_send_vector_ext_mesh,request_recv_vector_ext_mesh)
      call assemble_MPI_vector_async_recv(NPROC,NGLOB_AB,ssol%force_ext, &
                                          buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
                                          max_nibool_interfaces_ext_mesh, &
                                          nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                          request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                                          my_neighbors_ext_mesh)
    endif

    ! free tomo arrays
    if (allocated(tomo_x)) then 
      deallocate(tomo_x, tomo_y, tomo_z, tomo_force)
      deallocate(elem_force)
    endif

  end subroutine read_ext_force

  subroutine prepare_inv_mult()
    use specfem_par
    implicit none

    ! set inv_mat, used for parallel inner product
    inv_mult_(:) = 1.0_CUSTOM_REAL
    call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,inv_mult_, &
                                    num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                    nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                    my_neighbors_ext_mesh)
  
    ! inverse
    where(inv_mult_ == 0.0_CUSTOM_REAL)
      inv_mult_ = 1.0_CUSTOM_REAL
    end where
    inv_mult_ = 1.0_CUSTOM_REAL / inv_mult_

  end subroutine prepare_inv_mult

  subroutine trilinear_interp(nx,ny,nz,tomo_x,tomo_y,tomo_z,&
                            tomo_force,xp,yp,zp,fp)
    use constants,only: NDIM
    implicit none
    
    integer, intent(in) :: nx, ny, nz

    real(kind=dp), intent(in) :: tomo_x(nx), tomo_y(ny), tomo_z(nz)
    real(kind=dp), intent(in) :: tomo_force(nx,ny,nz,NDIM)
    real(kind=dp), intent(in) :: xp, yp, zp
    real(kind=dp), intent(out) :: fp(3)

  
    integer :: i,j,k
    real(kind=dp) :: x1,x2,y1,y2,z1,z2
    real(kind=dp) :: xd,yd,zd

    ! find indices for interpolation
    i = max(1, min(nx-1, int((xp - tomo_x(1)) / (tomo_x(2)-tomo_x(1))) + 1))
    j = max(1, min(ny-1, int((yp - tomo_y(1)) / (tomo_y(2)-tomo_y(1))) + 1))
    k = max(1, min(nz-1, int((zp - tomo_z(1)) / (tomo_z(2)-tomo_z(1))) + 1))

    ! get corner points
    x1 = tomo_x(i); x2 = tomo_x(i+1)
    y1 = tomo_y(j); y2 = tomo_y(j+1)
    z1 = tomo_z(k); z2 = tomo_z(k+1)

    ! compute interpolation weights
    xd = (xp - x1) / (x2 - x1)
    yd = (yp - y1) / (y2 - y1)
    zd = (zp - z1) / (z2 - z1)

    ! trilinear interpolation
    fp(:) = 0.0_CUSTOM_REAL
    fp(:) = fp(:) + tomo_force(i,j,k,:) * (1-xd)*(1-yd)*(1-zd)
    fp(:) = fp(:) + tomo_force(i+1,j,k,:) * xd*(1-yd)*(1-zd)
    fp(:) = fp(:) + tomo_force(i,j+1,k,:) * (1-xd)*yd*(1-zd)
    fp(:) = fp(:) + tomo_force(i,j,k+1,:) * (1-xd)*(1-yd)*zd
    fp(:) = fp(:) + tomo_force(i+1,j+1,k,:) * xd*yd*(1-zd)
    fp(:) = fp(:) + tomo_force(i+1,j,k+1,:) * xd*(1-yd)*zd
    fp(:) = fp(:) + tomo_force(i,j+1,k+1,:) * (1-xd)*yd*zd
    fp(:) = fp(:) + tomo_force(i+1,j+1,k+1,:) * xd*yd*zd

  end subroutine trilinear_interp

  subroutine create_petsc_backend(this)
    use constants, only: NGLLX,NGLLY,NGLLZ,NDIM
    use specfem_par, only: nglob => NGLOB_AB, nspec => NSPEC_AB, myrank, ibool
    use specfem_par, only: num_interfaces_ext_mesh, max_nibool_interfaces_ext_mesh, &
                          nibool_interfaces_ext_mesh, ibool_interfaces_ext_mesh, &
                          my_neighbors_ext_mesh,xstore,ystore,zstore
    use petsc_interfaces

    implicit none

    integer, parameter :: NGLL3 = NGLLX*NGLLY*NGLLZ
    class(static_solver_class), intent(inout) :: this

    if(.not. this%USE_PETSC_AS_BACKEND) return

    ! allocate owner ranks
    allocate(this%owner_rank(nglob))

    ! set initial value
    this%owner_rank(:) = myrank

    ! build ownership metadata and create the PETSc backend context
    call setup_petsc(this%petsc_ptr, nglob, myrank, &
                     nspec, NGLL3, NDIM, &
                     num_interfaces_ext_mesh, &
                     nibool_interfaces_ext_mesh, &
                     ibool_interfaces_ext_mesh, &
                     max_nibool_interfaces_ext_mesh, &
                     my_neighbors_ext_mesh, &
                     xstore, ystore, zstore, &
                     ibool, this%owner_rank)

  end subroutine create_petsc_backend

  subroutine destroy_petsc_backend(this)
    use petsc_interfaces
    implicit none

    class(static_solver_class), intent(inout) :: this

    if(.not. this%USE_PETSC_AS_BACKEND) return

    ! destroy PETSc solver context
    call cleanup_petsc(this%petsc_ptr)

    ! deallocate owner ranks
    if (allocated(this%owner_rank)) deallocate(this%owner_rank)

  end subroutine destroy_petsc_backend

  subroutine finalize_static_module
    implicit none

    ! deallocate arrays
    if (allocated(inv_mult_)) deallocate(inv_mult_)
    if (allocated(ssol%force_ext)) deallocate(ssol%force_ext)
    if (allocated(ssol%stress)) deallocate(ssol%stress)
    if (allocated(ssol%strain)) deallocate(ssol%strain)

    ! free petsc objects
    call ssol%free_petsc()

  end subroutine finalize_static_module


end submodule static_init