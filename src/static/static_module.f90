module static_module
  use iso_fortran_env, only: dp => real64
  use constants,only: NDIM
  implicit none
  private

  public:: initialize_static_module, finalize_static_module,solve_static_problem

  integer,parameter :: CUSTOM_REAL = 4

  real(kind=CUSTOM_REAL), dimension(:), allocatable :: inv_mult_ ! shape(NGLOB_AB)

  type static_solver_class 
    ! parameters 
    logical :: is_nonlinear = .false.
    logical :: SAVE_STATIC_FIELD = .false.
    real(kind=CUSTOM_REAL) :: omega(NDIM) = 0.0_CUSTOM_REAL ! rotation rate for centrifugal force calculation, only used when ROTATION is enabled
    real(kind=CUSTOM_REAL) :: rot_org(NDIM) = 0.0_CUSTOM_REAL ! rotation origin for centrifugal force calculation, only used when ROTATION is enabled

    real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: force_ext ! shape(NDIM,NGLOB_AB) external body force
    ! strain and stress tensors
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: stress,strain ! shape(NGLLX,NGLLY,NGLLZ,NSPEC,6)
    end type static_solver_class

  ! GLOBAL variable to hold static solver class
  type(static_solver_class) :: ssol 

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
      end if

      ! read static solver parameters
      ssol%is_nonlinear = .false. ! default value
      !call read_value_logical(ssol%SAVE_STATIC_FIELD, "SAVE_STATIC_FIELD",ier)
      ssol%SAVE_STATIC_FIELD = .true. ! default value

      ! close
      call close_parameter_file()
    end if

    ! broadcast parameters to all ranks
    call bcast_all_singlel(ssol%SAVE_STATIC_FIELD)
    call bcast_all_singlel(ssol%is_nonlinear)

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
    implicit none
    
    integer :: i,j,k,iglob,ier,ispec
    integer :: nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax 
    character(len=MAX_STRING_LEN) :: force_filename 
    integer,parameter :: IO_UNIT = 10
    real(kind=dp),allocatable :: tomo_x(:),tomo_y(:),tomo_z(:)
    real(kind=dp),allocatable :: tomo_force(:,:,:,:)
    real(kind=dp) :: fp(NDIM), temp
    real(kind=CUSTOM_REAL) :: omega(NDIM), rot_org(NDIM),tempx,tempy,tempz 
    real(kind=CUSTOM_REAL) :: xp, yp, zp

    ! read and interpolate force_ext from file
    if(myrank == 0) then 
      ! open Par_file
      call open_parameter_file(ier)
      if (ier /= 0) then
        print*, "Error opening parameter file for static module initialization"
        stop 
      end if

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
        end if

        ! read dimensions
        read(IO_UNIT,*) nx, ny, nz
        read(IO_UNIT,*) xmin, xmax, ymin, ymax, zmin, zmax

        ! allocate arrays
        allocate(tomo_x(nx), tomo_y(ny), tomo_z(nz))
        allocate(tomo_force(nx,ny,nz,NDIM))

        ! read data
        do k = 1,nz; do j = 1,ny; do i = 1,nx
          read(IO_UNIT,*) tomo_force(i,j,k,:)
        end do; end do; end do

        close(IO_UNIT)
      endif 

      ! set tomo_x/y/z based on bounds and dimensions
      do i = 1, nx
        tomo_x(i) = xmin + (i-1)*(xmax-xmin)/(nx-1)
      end do
      do j = 1, ny
        tomo_y(j) = ymin + (j-1)*(ymax-ymin)/(ny-1)
      end do
      do k = 1, nz
        tomo_z(k) = zmin + (k-1)*(zmax-zmin)/(nz-1)
      end do

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
    end if
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
         end if
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
      end if

      ! broadcast rotation parameters if needed
      call bcast_all_cr(omega,size(omega))
      call bcast_all_cr(rot_org,size(rot_org))
    end if

    ! set value to parameters for later use
    ssol%omega(:) = omega(:)
    ssol%rot_org(:) = rot_org(:)

    ! loop each global point and compute interpolated force, store in force_ext
    ssol%force_ext(:,:) = 0.0_CUSTOM_REAL
    do ispec = 1, nspec
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool(i,j,k,ispec)
            temp = jacobianstore(i,j,k,ispec) * wxgll(i)  * wygll(j) * wzgll(k)
            call trilinear_interp(nx,ny,nz,tomo_x,tomo_y,tomo_z,tomo_force,&
                                  dble(xstore(iglob)),dble(ystore(iglob)),dble(zstore(iglob)),&
                                  fp)
            ssol%force_ext(:,iglob) = ssol%force_ext(:,iglob) +  real(fp * temp,kind=CUSTOM_REAL)

          end do
        end do
      end do
    end do

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
            ssol%force_ext(1,iglob) = ssol%force_ext(1,iglob) - real(tempx * temp, kind=CUSTOM_REAL)
            ssol%force_ext(2,iglob) = ssol%force_ext(2,iglob) - real(tempy * temp, kind=CUSTOM_REAL)
            ssol%force_ext(3,iglob) = ssol%force_ext(3,iglob) - real(tempz * temp, kind=CUSTOM_REAL)

          end do
        end do
      end do
    end do

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

    ! free tomo arrays
    if (allocated(tomo_x)) then 
      deallocate(tomo_x, tomo_y, tomo_z, tomo_force)
    end if

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

  subroutine solve_static_problem()
    implicit none
    
    ! call the actual implementation
    call static_problem_impl()

    if(ssol%SAVE_STATIC_FIELD) then 
      call save_static_field_bin()
    end if

    ! save results on recievers 
    call save_field_at_receivers()

  end subroutine solve_static_problem

  subroutine save_static_field_bin()
    use constants, only: CUSTOM_REAL,MAX_STRING_LEN
    use specfem_par, only: NSPEC_AB,NGLLX,NGLLY,NGLLZ,ibool,myrank,LOCAL_PATH
    use specfem_par_elastic,only : displ
    implicit none

    !  local 
    integer :: i,j,k,ispec,ierr
    real(kind=CUSTOM_REAL), allocatable :: buffer(:,:,:,:)
    character(len=MAX_STRING_LEN) :: filename, comps(6)

    ! allocate buffer
    allocate(buffer(NGLLX,NGLLY,NGLLZ,NSPEC_AB)) 


    if(myrank == 0) then
      write(*,*) 
      write(*,*) '----------------------------------------'
      write(*,*) 'Writing static fields to disk ...'
      write(*,*) '----------------------------------------'
    endif

    ! displ_x/r 
    write(filename,'(A,I6.6,A)') trim(LOCAL_PATH) //'/proc',myrank,'_static_displ_x.bin'
    open(unit=20,file=trim(filename),status='replace',form='unformatted',iostat=ierr)
    if (ierr /= 0)  stop 'Error opening output file for static results'
    do ispec = 1,NSPEC_AB 
      do k = 1,NGLLZ; do j = 1,NGLLY; do i = 1,NGLLX
        buffer(i,j,k,ispec) = displ(1,ibool(i,j,k,ispec))
      enddo; enddo; enddo
    enddo
    write(20) buffer
    close(20)

    ! displ_y/phi
    write(filename,'(A,I6.6,A)') trim(LOCAL_PATH) //'/proc',myrank,'_static_displ_y.bin'
    open(unit=20,file=trim(filename),status='replace',form='unformatted',iostat=ierr)
    if (ierr /= 0)  stop 'Error opening output file for static results'
    do ispec = 1,NSPEC_AB
      do k = 1,NGLLZ; do j = 1,NGLLY; do i = 1,NGLLX
        buffer(i,j,k,ispec) =  displ(2,ibool(i,j,k,ispec))
      enddo; enddo; enddo
    enddo
    write(20) buffer
    close(20)

    ! displ_z
    write(filename,'(A,I6.6,A)') trim(LOCAL_PATH) //'/proc',myrank,'_static_displ_z.bin'
    open(unit=20,file=trim(filename),status='replace',form='unformatted',iostat=ierr)
    if (ierr /= 0)  stop 'Error opening output file for static results'
    do ispec = 1,NSPEC_AB
      do k = 1,NGLLZ; do j = 1,NGLLY; do i = 1,NGLLX
        buffer(i,j,k,ispec) =  displ(3,ibool(i,j,k,ispec))
      enddo; enddo; enddo
    enddo
    write(20) buffer
    close(20)

    ! stress tensor 
    comps = [character(len=2) :: 'xx','yy','zz','yz','xz','xy']
    do i = 1,6
      write(filename,'(A,I6.6,A)') trim(LOCAL_PATH) //'/proc',myrank,'_static_stress_'//trim(comps(i)) // '.bin'
      open(unit=20+i,file=trim(filename),status='replace',form='unformatted',iostat=ierr)
      if (ierr /= 0)  stop 'Error opening output file for static results'
      buffer(:,:,:,:) = ssol%stress(:,:,:,:,i)
      write(20+i) buffer
      close(20+i)

      ! strain tensor
      write(filename,'(A,I6.6,A)') trim(LOCAL_PATH) //'/proc',myrank,'_static_strain_'//trim(comps(i)) // '.bin'
      open(unit=30+i,file=trim(filename),status='replace',form='unformatted',iostat=ierr)
      if (ierr /= 0)  stop 'Error opening output file for static results'
      buffer(:,:,:,:) = ssol%strain(:,:,:,:,i)
      write(30+i) buffer
      close(30+i)
    enddo

  end subroutine save_static_field_bin

  subroutine save_field_at_receivers()
    ! for static problem, we can save the field at receivers as well, for easier comparison with observations
    ! this is optional and can be controlled by a parameter
    use specfem_par
    use specfem_par_elastic,only : displ
    implicit none

    ! local 
    integer :: ir,irloc,ispec,i,j,k,iglob,icomp 
    character(len=MAX_STRING_LEN) :: filename,sisname
    ! receiver Lagrange interpolators
    double precision,dimension(NGLLX) :: hxir
    double precision,dimension(NGLLY) :: hetar
    double precision,dimension(NGLLZ) :: hgammar
    double precision,dimension(NGLLX,NGLLY,NGLLZ) :: temp
    real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: seismo_static ! shape(15,nrec_local) in rank 0, where 15 = 3 displ + 6 stress + 6 strain
    
    allocate(seismo_static(15,nrec_local)) 

    do irloc = 1,nrec_local
      ! global receiver index
      ir = number_receiver_global(irloc)
      ispec = ispec_selected_rec(ir)

      ! compute Lagrange interpolators for receiver location
      hxir(:) = hxir_store(:,irloc)
      hetar(:) = hetar_store(:,irloc)
      hgammar(:) = hgammar_store(:,irloc)

      ! interpolate displacements
      do icomp = 1,15 
        do k = 1,NGLLZ; do j = 1,NGLLY; do i = 1,NGLLX
          if(icomp <= 3) then
            iglob = ibool(i,j,k,ispec)
            temp(i,j,k) = displ(icomp,iglob)
          else if (icomp <= 9) then 
            temp(i,j,k) = ssol%stress(i,j,k,ispec,icomp-3)
          else 
            temp(i,j,k) = ssol%strain(i,j,k,ispec,icomp-9)
          end if

          ! add interpolation weights
          temp(i,j,k) = temp(i,j,k) * hxir(i) * hetar(j) * hgammar(k)
        enddo; enddo; enddo

        seismo_static(icomp,irloc) = real(sum(temp), kind=CUSTOM_REAL)
      end do
    end do

    ! write to file on each rank 
    sisname = 'static.sem'
    filename = trim(OUTPUT_FILES) // '/' // trim(sisname)
    do i = 0, NPROC -1
      if(myrank == i) then
        if(myrank == 0) then 
          open(unit=40,file=trim(filename),status='replace',form='formatted')
        else 
          open(unit=40,file=trim(filename),status='old',form='formatted',position='append')
        end if

        do irloc = 1,nrec_local
          ir = number_receiver_global(irloc)
          write(40,12313) station_name(ir),network_name(ir), &
                           seismo_static(:,irloc)
        end do
      endif 
      call synchronize_all()
    end do
    call synchronize_all()

12313 format(2(a,1x),15(g0,1x))

    ! deallocate
    deallocate(seismo_static)
    
  end subroutine save_field_at_receivers

  !> main subroutine for static solution, CG metho is used
  subroutine static_problem_impl()
    use constants, only: MAX_STRING_LEN,NDIM
    use specfem_par, only: NGLOB_AB,NGLLX,NGLLY,NGLLZ,myrank 
    use specfem_par_elastic, only : displ
    implicit none

    ! local
    integer :: i, maxiter
    real(kind=dp) :: alpha, beta, rsnew, temp_sum, rsinit, max_u, max_st
    real(kind=dp) :: rho_old, rho_new ! Variables for the PCG dot product (r . z)
    
    ! Added 'z' for the preconditioned residual
    real(kind=dp), dimension(:, :), allocatable :: r, p, Ap, z, inv_pred, u 
    real(kind=CUSTOM_REAL),dimension(:,:), allocatable :: p_cr, Ap_cr, z_cr,kdotu
    
    allocate(kdotu(NDIM,NGLOB_AB), r(NDIM,NGLOB_AB), p(NDIM,NGLOB_AB), &
             Ap(NDIM,NGLOB_AB), z(NDIM,NGLOB_AB),inv_pred(NDIM,NGLOB_AB))
    allocate(p_cr(NDIM,NGLOB_AB), Ap_cr(NDIM,NGLOB_AB), z_cr(NDIM,NGLOB_AB),&
             u(NDIM,NGLOB_AB))

    ! preconditioner are set to 1
    inv_pred(:,:) = 1.0_dp

    ! -----------------------------------------------------
    ! 1. Initialization
    ! -----------------------------------------------------
    ! initial guess: u = 0
    u = 0.0 
    displ = 0.0
    call compute_forces_static(displ, .false., kdotu, ssol%stress, ssol%strain)
    r = ssol%force_ext - kdotu
    call enforce_fixed_bc(r)  
    call enforce_roller_bc(r)
    
    ! Calculate Initial Norm for convergence check (r . r)
    rsinit = parallel_inner_product(r, r)
    if (myrank == 0) then
      write(*,*)
      write(*,*) '----------------------------------------'
      write (*,*) 'Starting PCG Solver for Static Problem'
      write(*,*) 'Initial residual norm: ', sqrt(rsinit)
      write(*,*) '----------------------------------------'
      write(*,*)
    end if
    
    ! Safety for division by zero if exact solution is 0
    if (rsinit < 1.0e-20) rsinit = 1.0

    ! -----------------------------------------------------
    ! 2. Apply Preconditioner (Initial Step)
    !    z = M^-1 * r
    ! -----------------------------------------------------
    ! We use safe division to handle potential zeros in rmass (e.g. ghosts)
    z = r * inv_pred

    ! p = z (Start direction is the preconditioned residual)
    p = z
    
    ! Calculate rho_old = r . z (This replaces rsold in standard CG)
    rho_old = parallel_inner_product(r, z)

    ! -----------------------------------------------------
    ! 3. PCG Loop
    ! -----------------------------------------------------
    maxiter = NGLOB_AB
    do i = 1, maxiter 
      
      ! Apply BCs to p to be safe
      call enforce_fixed_bc(p)
      call enforce_roller_bc(p)

      ! Matrix-Vector Product: Ap = K * p
      p_cr = real(p, kind=CUSTOM_REAL)
      call compute_forces_static(p_cr, .false., Ap_cr, ssol%stress, ssol%strain)
      Ap = real(Ap_cr, kind=dp)

      ! Calculate Alpha (Step Size)
      ! alpha = (r_old . z_old) / (p . Ap)
      temp_sum = parallel_inner_product(p, Ap)
      
      ! Safety check for zero/negative curvature (indefinite matrix)
      if (temp_sum < 1.0e-30_CUSTOM_REAL) then
          if (myrank == 0) write(*,*) "Warning: Curvature too small/negative, stopping.", i,temp_sum
          exit 
      endif

      alpha = rho_old / temp_sum

      ! Update Solution: x = x + alpha * p
      u = u + alpha * p
      displ = real(u, kind=CUSTOM_REAL)

      ! Update Residual: r = r - alpha * Ap
      if (1 == 0) then 
        ! Restart to kill floating point drift (Increased to 50 for efficiency)
        call compute_forces_static(displ, .false., kdotu, ssol%stress, ssol%strain)
        r = ssol%force_ext - kdotu
      else
        r = r - alpha * Ap
      endif
      call enforce_fixed_bc(r)
      call enforce_roller_bc(r)

      ! Check Convergence (using true residual norm r . r)
      rsnew = parallel_inner_product(r, r)
      
      ! Print status
      if ((mod(i,100) == 0 .or. i == 5 .or. i == maxiter)) then
        temp_sum = maxval(abs(u))
        call max_all_all_dp(temp_sum,max_u)
        if (myrank == 0) &
          write(*,*) 'Iter:', i, ' Rel Resid:', sqrt(rsnew / rsinit), &
                    ' norm of Displ:', max_u
      endif
      
      if (sqrt(rsnew / rsinit) < 1.0e-6_CUSTOM_REAL) exit 

      ! -------------------------------------------------
      ! 4. Apply Preconditioner for Next Step
      !    z = M^-1 * r
      ! -------------------------------------------------
      z = r * inv_pred

      ! Calculate rho_new = r . z
      rho_new = parallel_inner_product(r, z)
      
      ! Calculate Beta (Polak-Ribiere/Fletcher-Reeves equivalent for PCG)
      beta = rho_new / rho_old
      
      ! Update Search Direction: p = z + beta * p
      p = z + beta * p
      
      rho_old = rho_new
    enddo

    ! apply boundary conditions to u
    call enforce_fixed_bc(u)
    call enforce_roller_bc(u)

    ! call again to compute new stress
    displ = real(u, kind=CUSTOM_REAL)
    call compute_forces_static(displ, .false., kdotu, ssol%stress, ssol%strain)

    ! check max strain
    temp_sum = maxval(abs(ssol%strain))
    call max_all_all_dp(temp_sum,max_st)
    if (myrank == 0) write(*,*) 'Max Strain:', max_st

    ! deallocate
    deallocate(kdotu, r, p, Ap, z, inv_pred)
    deallocate(u, p_cr, Ap_cr, z_cr)

  end subroutine static_problem_impl

  subroutine compute_forces_static(displ,is_nonlinear,kdotu,stress_tensor,strain_tensor)
    !! compute forces kdotu = K * u, and compute stress tensor
    !! 
    use specfem_par
    implicit none

    real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB), intent(in) :: displ
    logical, intent(in) :: is_nonlinear
    real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB), intent(out) :: kdotu
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), intent(out) :: stress_tensor,strain_tensor ! shape(NGLLX,NGLLY,NGLLZ,NSPEC,6)

    ! initialize kdotu
    kdotu(:,:) = 0.0_CUSTOM_REAL

    ! compute forces for outer phase
    call compute_forces_phase(displ,1,is_nonlinear,kdotu,stress_tensor,strain_tensor)

    call assemble_MPI_vector_async_send(NPROC,NGLOB_AB,kdotu, &
                                        buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                        my_neighbors_ext_mesh, &
                                        request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

    ! compute forces for inner phase
    call compute_forces_phase(displ,2,is_nonlinear,kdotu,stress_tensor,strain_tensor)


    ! receives MPI buffers
    call assemble_MPI_vector_async_recv(NPROC,NGLOB_AB,kdotu, &
                                        buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
                                        max_nibool_interfaces_ext_mesh, &
                                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                        request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                                        my_neighbors_ext_mesh)


  end subroutine compute_forces_static


  pure subroutine second2firstPK_stress( &
    sigma_xx,sigma_xy,sigma_xz, &
    sigma_yx,sigma_yy,sigma_yz, &
    sigma_zx,sigma_zy,sigma_zz, &
    duxdxl,duxdyl,duxdzl, &
    duydxl,duydyl,duydzl, &
    duzdxl,duzdyl,duzdzl)
    
    use constants, only: CUSTOM_REAL
    implicit none

    real(kind=CUSTOM_REAL), intent(inout) :: sigma_xx,sigma_xy,sigma_xz
    real(kind=CUSTOM_REAL), intent(inout) :: sigma_yx,sigma_yy,sigma_yz
    real(kind=CUSTOM_REAL), intent(inout) :: sigma_zx,sigma_zy,sigma_zz
    real(kind=CUSTOM_REAL), intent(in) :: duxdxl,duxdyl,duxdzl
    real(kind=CUSTOM_REAL), intent(in) :: duydxl,duydyl,duydzl
    real(kind=CUSTOM_REAL), intent(in) :: duzdxl,duzdyl,duzdzl

    ! backup original stresses
    real(kind=CUSTOM_REAL) :: s_xx,s_xy,s_xz
    real(kind=CUSTOM_REAL) :: s_yx,s_yy,s_yz
    real(kind=CUSTOM_REAL) :: s_zx,s_zy,s_zz
    s_xx = sigma_xx; s_xy = sigma_xy; s_xz = sigma_xz;
    s_yx = sigma_yx; s_yy = sigma_yy; s_yz = sigma_yz;
    s_zx = sigma_zx; s_zy = sigma_zy; s_zz = sigma_zz;

    ! update to first PK stress
    ! P_{ij} = S_{ij} + S_{ik} * du_{j}/dx_{k}
    sigma_xx = s_xx + s_xx * duxdxl + s_xy * duxdyl + s_xz * duxdzl
    sigma_xy = s_xy + s_xx * duydxl + s_xy * duydyl + s_xz * duydzl
    sigma_xz = s_xz + s_xx * duzdxl + s_xy * duzdyl + s_xz * duzdzl
    sigma_yx = s_yx + s_yx * duxdxl + s_yy * duxdyl + s_yz * duxdzl
    sigma_yy = s_yy + s_yx * duydxl + s_yy * duydyl + s_yz * duydzl
    sigma_yz = s_yz + s_yx * duzdxl + s_yy * duzdyl + s_yz * duzdzl
    sigma_zx = s_zx + s_zx * duxdxl + s_zy * duxdyl + s_zz * duxdzl
    sigma_zy = s_zy + s_zx * duydxl + s_zy * duydyl + s_zz * duydzl
    sigma_zz = s_zz + s_zx * duzdxl + s_zy * duzdyl + s_zz * duzdzl

  end subroutine second2firstPK_stress


  subroutine compute_forces_phase(displ,iphase,is_nonlinear,kdotu,stress_tensor,strain_tensor)
    use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,&
                        N_SLS,ONE_THIRD,FOUR_THIRDS,m1,m2
    use specfem_par, only: xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                          gammaxstore,gammaystore,gammazstore,jacobianstore, &
                          NGLOB_AB, rhostore,&
                          hprime_xx,hprime_xxT, &
                          hprimewgll_xx,hprimewgll_xxT, &
                          kappastore,mustore,ibool,ANISOTROPY,ROTATION
    use specfem_par, only: wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D,wxgll,wygll,wzgll
    use specfem_par_elastic, only: c11store,c12store,c13store,c14store,c15store,c16store, &
                                  c22store,c23store,c24store,c25store,c26store,c33store, &
                                  c34store,c35store,c36store,c44store,c45store,c46store, &
                                  c55store,c56store,c66store,&
                                  nspec_inner_elastic,nspec_outer_elastic,phase_ispec_inner_elastic

    implicit none
    
    real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB), intent(in) :: displ
    real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB), intent(inout) :: kdotu
    integer, intent(in) :: iphase
    logical, intent(in) :: is_nonlinear
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), intent(out) :: stress_tensor,strain_tensor ! shape(NGLLX,NGLLY,NGLLZ,NSPEC,6)

    ! local 
    integer :: num_elements
    integer :: ispec, i,j,k,iglob, ispec_p
    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc
    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1,tempy1,tempz1
    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx2,tempy2,tempz2
    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx3,tempy3,tempz3
    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: duxdxl,duxdyl,duxdzl
    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: duydxl,duydyl,duydzl
    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: duzdxl,duzdyl,duzdzl
    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: newtempx1,newtempy1,newtempz1
    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: newtempx2,newtempy2,newtempz2
    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: newtempx3,newtempy3,newtempz3
    real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
    real(kind=CUSTOM_REAL) :: duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl
    real(kind=CUSTOM_REAL) :: sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy

    ! local material parameters
    real(kind=CUSTOM_REAL) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                              c33,c34,c35,c36,c44,c45,c46,c55,c56,c66
    real(kind=CUSTOM_REAL) :: lambdal,mul,lambdalplus2mul
    real(kind=CUSTOM_REAL) :: kappal,fac1,fac2,fac3
    real(kind=CUSTOM_REAL) :: strain_xx,strain_yy,strain_zz,strain_xy,strain_xz,strain_yz

    ! choses inner/outer elements
    if (iphase == 1) then
      num_elements = nspec_outer_elastic
    else
      num_elements = nspec_inner_elastic
    endif

    do ispec_p = 1,num_elements
      ispec = phase_ispec_inner_elastic(ispec_p,iphase)
      
      ! gather local displacements
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob = ibool(i,j,k,ispec)
            dummyx_loc(i,j,k) = displ(1,iglob)
            dummyy_loc(i,j,k) = displ(2,iglob)
            dummyz_loc(i,j,k) = displ(3,iglob)
          enddo
        enddo
      enddo

      ! compute gradient
      call mxm5_3comp_singleA(hprime_xx,m1,dummyx_loc,dummyy_loc,dummyz_loc,tempx1,tempy1,tempz1,m2)
      call mxm5_3comp_3dmat_single(dummyx_loc,dummyy_loc,dummyz_loc,m1,hprime_xxT,m1,tempx2,tempy2,tempz2,m1)
      call mxm5_3comp_singleB(dummyx_loc,dummyy_loc,dummyz_loc,m2,hprime_xxT,tempx3,tempy3,tempz3,m1)
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            xixl = xixstore(i,j,k,ispec)
            xiyl = xiystore(i,j,k,ispec)
            xizl = xizstore(i,j,k,ispec)
            etaxl = etaxstore(i,j,k,ispec)
            etayl = etaystore(i,j,k,ispec)
            etazl = etazstore(i,j,k,ispec)
            gammaxl = gammaxstore(i,j,k,ispec)
            gammayl = gammaystore(i,j,k,ispec)
            gammazl = gammazstore(i,j,k,ispec)
            
            ! displ gradients
            duxdxl(i,j,k) = xixl*tempx1(i,j,k) + etaxl*tempx2(i,j,k) + gammaxl*tempx3(i,j,k)
            duxdyl(i,j,k) = xiyl*tempx1(i,j,k) + etayl*tempx2(i,j,k) + gammayl*tempx3(i,j,k)
            duxdzl(i,j,k) = xizl*tempx1(i,j,k) + etazl*tempx2(i,j,k) + gammazl*tempx3(i,j,k)

            duydxl(i,j,k) = xixl*tempy1(i,j,k) + etaxl*tempy2(i,j,k) + gammaxl*tempy3(i,j,k)
            duydyl(i,j,k) = xiyl*tempy1(i,j,k) + etayl*tempy2(i,j,k) + gammayl*tempy3(i,j,k)
            duydzl(i,j,k) = xizl*tempy1(i,j,k) + etazl*tempy2(i,j,k) + gammazl*tempy3(i,j,k)

            duzdxl(i,j,k) = xixl*tempz1(i,j,k) + etaxl*tempz2(i,j,k) + gammaxl*tempz3(i,j,k)
            duzdyl(i,j,k) = xiyl*tempz1(i,j,k) + etayl*tempz2(i,j,k) + gammayl*tempz3(i,j,k)
            duzdzl(i,j,k) = xizl*tempz1(i,j,k) + etazl*tempz2(i,j,k) + gammazl*tempz3(i,j,k)
          enddo
        enddo
      enddo

      ! compute stresses and assemble forces
      do k = 1,NGLLZ;do j = 1,NGLLY; do i = 1,NGLLX
        ! compute strain components
        strain_xx = duxdxl(i,j,k)
        strain_yy = duydyl(i,j,k)
        strain_zz = duzdzl(i,j,k)
        strain_xy = 0.5_CUSTOM_REAL * (duxdyl(i,j,k) + duydxl(i,j,k))
        strain_xz = 0.5_CUSTOM_REAL * (duzdxl(i,j,k) + duxdzl(i,j,k))
        strain_yz = 0.5_CUSTOM_REAL * (duzdyl(i,j,k) + duydzl(i,j,k))

        if (is_nonlinear) then
          ! add finite strain terms here if needed
          strain_xx = strain_xx + 0.5_CUSTOM_REAL * (duxdxl(i,j,k)**2 + &
                                                     duydxl(i,j,k)**2 + &
                                                     duzdxl(i,j,k)**2)                                 
          strain_yy = strain_yy + 0.5_CUSTOM_REAL * (duxdyl(i,j,k)**2 + &
                                                     duydyl(i,j,k)**2 + &
                                                     duzdyl(i,j,k)**2)
                                                     
          strain_zz = strain_zz + 0.5_CUSTOM_REAL * (duxdzl(i,j,k)**2 + &
                                                     duydzl(i,j,k)**2 + &
                                                     duzdzl(i,j,k)**2)
          strain_xy = strain_xy + 0.5_CUSTOM_REAL * (duxdxl(i,j,k) * duxdyl(i,j,k) + &
                                                     duydxl(i,j,k) * duydyl(i,j,k) + &
                                                     duzdxl(i,j,k) * duzdyl(i,j,k))
          strain_xz = strain_xz + 0.5_CUSTOM_REAL * (duxdxl(i,j,k) * duxdzl(i,j,k) + & 
                                                     duydxl(i,j,k) * duydzl(i,j,k) + &
                                                     duzdxl(i,j,k) * duzdzl(i,j,k))
          strain_yz = strain_yz + 0.5_CUSTOM_REAL * (duxdyl(i,j,k) * duxdzl(i,j,k) + &
                                                    duydyl(i,j,k) * duydzl(i,j,k) + &
                                                    duzdyl(i,j,k) * duzdzl(i,j,k))
        endif

        ! save strain/Green strain tensor, voigt notation
        strain_tensor(i,j,k,ispec,1) = strain_xx
        strain_tensor(i,j,k,ispec,2) = strain_yy
        strain_tensor(i,j,k,ispec,3) = strain_zz
        strain_tensor(i,j,k,ispec,4) = strain_yz
        strain_tensor(i,j,k,ispec,5) = strain_xz
        strain_tensor(i,j,k,ispec,6) = strain_xy

        ! precompute some sums to save CPU time
        duxdyl_plus_duydxl = 2.0_CUSTOM_REAL * strain_xy
        duzdxl_plus_duxdzl = 2.0_CUSTOM_REAL * strain_xz
        duzdyl_plus_duydzl = 2.0_CUSTOM_REAL * strain_yz

        ! computes either isotropic or anisotropic element stresses
        if (ANISOTROPY) then
          ! full anisotropic case, stress calculations
          c11 = c11store(i,j,k,ispec)
          c12 = c12store(i,j,k,ispec)
          c13 = c13store(i,j,k,ispec)
          c14 = c14store(i,j,k,ispec)
          c15 = c15store(i,j,k,ispec)
          c16 = c16store(i,j,k,ispec)
          c22 = c22store(i,j,k,ispec)
          c23 = c23store(i,j,k,ispec)
          c24 = c24store(i,j,k,ispec)
          c25 = c25store(i,j,k,ispec)
          c26 = c26store(i,j,k,ispec)
          c33 = c33store(i,j,k,ispec)
          c34 = c34store(i,j,k,ispec)
          c35 = c35store(i,j,k,ispec)
          c36 = c36store(i,j,k,ispec)
          c44 = c44store(i,j,k,ispec)
          c45 = c45store(i,j,k,ispec)
          c46 = c46store(i,j,k,ispec)
          c55 = c55store(i,j,k,ispec)
          c56 = c56store(i,j,k,ispec)
          c66 = c66store(i,j,k,ispec)

          sigma_xx = c11 * strain_xx + c16 * duxdyl_plus_duydxl + c12 * strain_yy + &
                    c15 * duzdxl_plus_duxdzl + c14 * duzdyl_plus_duydzl + c13 * strain_zz
          sigma_yy = c12 * strain_xx + c26 * duxdyl_plus_duydxl + c22 * strain_yy + &
                    c25 * duzdxl_plus_duxdzl + c24 * duzdyl_plus_duydzl + c23 * strain_zz
          sigma_zz = c13 * strain_xx + c36 * duxdyl_plus_duydxl + c23 * strain_yy + &
                    c35 * duzdxl_plus_duxdzl + c34 * duzdyl_plus_duydzl + c33 * strain_zz
          sigma_xy = c16 * strain_xx + c66 * duxdyl_plus_duydxl + c26 * strain_yy + &
                    c56 * duzdxl_plus_duxdzl + c46 * duzdyl_plus_duydzl + c36 * strain_zz
          sigma_xz = c15 * strain_xx + c56 * duxdyl_plus_duydxl + c25 * strain_yy + &
                    c55 * duzdxl_plus_duxdzl + c45 * duzdyl_plus_duydzl + c35 * strain_zz
          sigma_yz = c14 * strain_xx + c46 * duxdyl_plus_duydxl + c24 * strain_yy + &
                    c45 * duzdxl_plus_duxdzl + c44 * duzdyl_plus_duydzl + c34 * strain_zz

        else
          ! isotropic case
          kappal = kappastore(i,j,k,ispec)
          mul = mustore(i,j,k,ispec)

          lambdalplus2mul = kappal + FOUR_THIRDS * mul
          lambdal = lambdalplus2mul - 2._CUSTOM_REAL * mul

          ! compute stress sigma
          sigma_xx = lambdalplus2mul * strain_xx + lambdal * (strain_yy + strain_zz)
          sigma_yy = lambdalplus2mul * strain_yy + lambdal * (strain_xx + strain_zz)
          sigma_zz = lambdalplus2mul * strain_zz + lambdal * (strain_xx + strain_yy)

          sigma_xy = mul * duxdyl_plus_duydxl
          sigma_xz = mul * duzdxl_plus_duxdzl
          sigma_yz = mul * duzdyl_plus_duydzl
        endif ! ANISOTROPY

        ! note stress is the second Piola-Kirchhoff stress if is_nonlinear is true
        ! save stress tensor,voigt notation
        stress_tensor(i,j,k,ispec,1) = sigma_xx
        stress_tensor(i,j,k,ispec,2) = sigma_yy
        stress_tensor(i,j,k,ispec,3) = sigma_zz
        stress_tensor(i,j,k,ispec,4) = sigma_yz
        stress_tensor(i,j,k,ispec,5) = sigma_xz
        stress_tensor(i,j,k,ispec,6) = sigma_xy

        ! symmetric stresses
        sigma_yx = sigma_xy
        sigma_zx = sigma_xz
        sigma_zy = sigma_yz

        if(is_nonlinear) then 
          ! we should use first Piola-Kirchhoff stress for force computation
          call second2firstPK_stress(&
               sigma_xx,sigma_xy,sigma_xz, &
                sigma_yx,sigma_yy,sigma_yz, &
                sigma_zx,sigma_zy,sigma_zz, &
               duxdxl(i,j,k),duxdyl(i,j,k),duxdzl(i,j,k), &
               duydxl(i,j,k),duydyl(i,j,k),duydzl(i,j,k), &
               duzdxl(i,j,k),duzdyl(i,j,k),duzdzl(i,j,k) &
          )
        endif 

        xixl = xixstore(i,j,k,ispec)
        xiyl = xiystore(i,j,k,ispec)
        xizl = xizstore(i,j,k,ispec)
        etaxl = etaxstore(i,j,k,ispec)
        etayl = etaystore(i,j,k,ispec)
        etazl = etazstore(i,j,k,ispec)
        gammaxl = gammaxstore(i,j,k,ispec)
        gammayl = gammaystore(i,j,k,ispec)
        gammazl = gammazstore(i,j,k,ispec)
        jacobianl = jacobianstore(i,j,k,ispec)

        ! form dot product with test vector, non-symmetric form (which is useful in the case of PML)
        tempx1(i,j,k) = jacobianl * (sigma_xx * xixl + sigma_yx * xiyl + sigma_zx * xizl) ! this goes to accel_x
        tempy1(i,j,k) = jacobianl * (sigma_xy * xixl + sigma_yy * xiyl + sigma_zy * xizl) ! this goes to accel_y
        tempz1(i,j,k) = jacobianl * (sigma_xz * xixl + sigma_yz * xiyl + sigma_zz * xizl) ! this goes to accel_z

        tempx2(i,j,k) = jacobianl * (sigma_xx * etaxl + sigma_yx * etayl + sigma_zx * etazl) ! this goes to accel_x
        tempy2(i,j,k) = jacobianl * (sigma_xy * etaxl + sigma_yy * etayl + sigma_zy * etazl) ! this goes to accel_y
        tempz2(i,j,k) = jacobianl * (sigma_xz * etaxl + sigma_yz * etayl + sigma_zz * etazl) ! this goes to accel_z

        tempx3(i,j,k) = jacobianl * (sigma_xx * gammaxl + sigma_yx * gammayl + sigma_zx * gammazl) ! this goes to accel_x
        tempy3(i,j,k) = jacobianl * (sigma_xy * gammaxl + sigma_yy * gammayl + sigma_zy * gammazl) ! this goes to accel_y
        tempz3(i,j,k) = jacobianl * (sigma_xz * gammaxl + sigma_yz * gammayl + sigma_zz * gammazl) ! this goes to accel_z
      enddo; enddo; enddo

      ! dot product with test functions
      call mxm5_3comp_singleA(hprimewgll_xxT,m1,tempx1,tempy1,tempz1,newtempx1,newtempy1,newtempz1,m2)
      call mxm5_3comp_3dmat_single(tempx2,tempy2,tempz2,m1,hprimewgll_xx,m1,newtempx2,newtempy2,newtempz2,m1)
      call mxm5_3comp_singleB(tempx3,tempy3,tempz3,m2,hprimewgll_xx,newtempx3,newtempy3,newtempz3,m1)

      ! scatter local forces to global force vector
      do k = 1,NGLLZ;do j = 1,NGLLY; do i = 1,NGLLX
        iglob = ibool(i,j,k,ispec)

        ! compute rotational forces if needed
        if(ROTATION) then  
          ! add centrifugal forces 
          ! omega x u 
          call cross_product(ssol%omega(1),ssol%omega(2),ssol%omega(3), &
                             dummyx_loc(i,j,k), dummyy_loc(i,j,k),dummyz_loc(i,j,k), &
                             fac1,fac2,fac3)
          
          call cross_product(ssol%omega(1),ssol%omega(2),ssol%omega(3), &
                             fac1,fac2,fac3, &
                             c11,c22,c33)
          ! times rho and weights
          jacobianl = jacobianstore(i,j,k,ispec) * rhostore(i,j,k,ispec) * real(wxgll(i)*wygll(j)*wzgll(k),kind=CUSTOM_REAL)
          c11 = c11 * jacobianl
          c22 = c22 * jacobianl
          c33 = c33 * jacobianl
        else 
          c11 = 0.0_CUSTOM_REAL
          c22 = 0.0_CUSTOM_REAL
          c33 = 0.0_CUSTOM_REAL    
        end if

        fac1 = wgllwgll_yz_3D(i,j,k) !or wgllwgll_yz(j,k)
        fac2 = wgllwgll_xz_3D(i,j,k) !or wgllwgll_xz(i,k)
        fac3 = wgllwgll_xy_3D(i,j,k) !or wgllwgll_xy(i,j)

        ! accumulate forces
        c11 = c11 + fac1 * newtempx1(i,j,k) + fac2 * newtempx2(i,j,k) + fac3 * newtempx3(i,j,k)
        c22 = c22 + fac1 * newtempy1(i,j,k) + fac2 * newtempy2(i,j,k) + fac3 * newtempy3(i,j,k)
        c33 = c33 + fac1 * newtempz1(i,j,k) + fac2 * newtempz2(i,j,k) + fac3 * newtempz3(i,j,k)
        ! c11 = -c11
        ! c22 = -c22
        ! c33 = -c33

        kdotu(1,iglob) = kdotu(1,iglob) + c11
        kdotu(2,iglob) = kdotu(2,iglob) + c22
        kdotu(3,iglob) = kdotu(3,iglob) + c33
      enddo; enddo; enddo

    enddo ! ispec_p
  end subroutine compute_forces_phase

  subroutine finalize_static_module
    implicit none

    ! deallocate arrays
    if (allocated(inv_mult_)) deallocate(inv_mult_)
    if (allocated(ssol%force_ext)) deallocate(ssol%force_ext)
    if (allocated(ssol%stress)) deallocate(ssol%stress)
    if (allocated(ssol%strain)) deallocate(ssol%strain)

  end subroutine finalize_static_module

  function parallel_inner_product(vec1,vec2) result(psum)
    use specfem_par, only: NGLOB_AB,NDIM
    implicit none

    real(kind=dp), dimension(NDIM, NGLOB_AB), intent(in) :: vec1, vec2
    real(kind=dp) :: psum 

    ! local
    real(kind=dp) :: local_sum,temp,a, b 

    ! to avoid underflow/overflow, we scale the vectors by their max values before doing the 
    ! inner product, and then scale back the sum by a*b at the end. 
    ! This does not change the result but helps with numerical stability
    a = maxval(abs(vec1))
    b = maxval(abs(vec2))
    call max_all_all_dp(a,temp)
    a = temp 
    call max_all_all_dp(b,temp)
    b = temp

    ! avoid very small scaling factors which can cause overflow in the inner product
    if (a < 1.0e-20_dp) a = 1.0_dp
    if (b < 1.0e-20_dp) b = 1.0_dp
    
    ! local sum
    local_sum = sum(vec1 / a * vec2 / b * spread(inv_mult_,1,NDIM))
    
    ! mpi reduce
    call sum_all_all_dp(local_sum,temp)
    psum = temp * a * b ! scale back the sum by the max values of the vectors to get the correct inner product

  end function parallel_inner_product

  ! utilities for boundary conditions
  pure subroutine enforce_roller_bc(vec)
    use constants,only: NGLLSQUARE
    use specfem_par,only: NGLOB_AB, NDIM,ibool 
    use specfem_par,only : num_roller_bdry_faces, roller_bdry_ispec, roller_bdry_ijk, roller_bdry_normal

    implicit none

    real(kind=dp),intent(inout) :: vec(NDIM,NGLOB_AB)

    ! local
    integer :: iglob,ispec,i,j,k,iface,igll2
    logical(kind=1) :: mask_nodes(NGLOB_AB)
    real(kind=dp) :: normal_vec(NDIM)

    mask_nodes(:) = .false.
    
    do iface = 1, num_roller_bdry_faces
      ispec = roller_bdry_ispec(iface)
      do igll2 = 1, NGLLSQUARE
        i = roller_bdry_ijk(1,igll2,iface)
        j = roller_bdry_ijk(2,igll2,iface)
        k = roller_bdry_ijk(3,igll2,iface)
        iglob = ibool(i,j,k,ispec)

        if(mask_nodes(iglob)) cycle ! skip if node has already been processed by another face

        normal_vec(:) = roller_bdry_normal(:,igll2,iface)
        vec(:,iglob) = vec(:,iglob) - dot_product(vec(:,iglob),normal_vec)*normal_vec
        mask_nodes(iglob) = .true.
      enddo
    end do 

  end subroutine enforce_roller_bc

  pure subroutine enforce_fixed_bc(vec)
    use constants,only: NGLLSQUARE
    use specfem_par,only: NGLOB_AB, NDIM,ibool
    use specfem_par,only : num_fixed_bdry_faces, fixed_bdry_ijk,fixed_bdry_ispec

    implicit none

    real(kind=dp),intent(inout) :: vec(NDIM,NGLOB_AB)

    !local
    integer :: iglob,ispec,i,j,k,iface,igll2

    do iface = 1, num_fixed_bdry_faces
      ispec = fixed_bdry_ispec(iface)
      do igll2 = 1, NGLLSQUARE
        i = fixed_bdry_ijk(1,igll2,iface)
        j = fixed_bdry_ijk(2,igll2,iface)
        k = fixed_bdry_ijk(3,igll2,iface)
        iglob = ibool(i,j,k,ispec)

        vec(:,iglob) = 0.0_dp
      enddo
    end do

  end subroutine enforce_fixed_bc

  pure subroutine cross_product(x1,x2,x3,y1,y2,y3,z1,z2,z3)
    implicit none

    real(kind=CUSTOM_REAL), intent(in) :: x1,x2,x3
    real(kind=CUSTOM_REAL), intent(in) :: y1,y2,y3
    real(kind=CUSTOM_REAL), intent(out) :: z1,z2,z3

    z1 = x2 * y3 - x3 * y2
    z2 = x3 * y1 - x1 * y3
    z3 = x1 * y2 - x2 * y1

  end subroutine cross_product

end module static_module