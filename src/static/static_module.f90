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

    ! for PETSc parallel matrix assembly
    logical :: USE_PETSC_AS_BACKEND = .false. ! whether to use PETSc as the linear solver backend, if false, a simple CG solver implemented in Fortran will be used, this is mainly for testing and debugging purposes
    integer(kind=8) :: petsc_ptr ! dummy variable to hold PETSc pointers as integers, will be cast to proper types in C
    integer, dimension(:), allocatable :: owner_rank ! shape(NGLOB_AB), stores the owning rank for each global DOF, used for parallel assembly with PETSc

    ! for GPU_MODE: holds the CGSolver* from static_gpu.cu (0 until lazily created by
    ! create_cg_gpu_backend(), cast to/from a C pointer the same way petsc_ptr is)
    integer(kind=8) :: cg_gpu_ptr = 0

    contains
      procedure :: init_petsc => create_petsc_backend
      procedure :: free_petsc => destroy_petsc_backend
      procedure :: execute => solve_static_problem_petsc
  end type static_solver_class


  interface
    module subroutine initialize_static_module()
    end subroutine initialize_static_module

    module subroutine finalize_static_module()
    end subroutine finalize_static_module

    module subroutine solve_static_problem()
    end subroutine solve_static_problem

    module subroutine create_petsc_backend(this)
      class(static_solver_class), intent(inout) :: this
    end subroutine create_petsc_backend

    module subroutine destroy_petsc_backend(this)
      class(static_solver_class), intent(inout) :: this
    end subroutine destroy_petsc_backend

    ! NOTE: everything below is declared here (rather than left as an ordinary contains-only
    ! procedure) purely so it keeps proper external linkage when called from a *different* file
    ! (static_init.f90 / static_impl.f90). Without a "module subroutine"/"module function"
    ! declaration + submodule implementation, gfortran is free to fully inline/eliminate a
    ! private module procedure that's only ever called locally *as far as this translation unit
    ! can see* -- it has no way of knowing another file will call it via host association, so the
    ! symbol silently disappears from the compiled .o and the final link fails with "undefined
    ! reference". See static_impl.f90 for the implementations.

    module subroutine static_problem_impl()
    end subroutine static_problem_impl

    module subroutine create_cg_gpu_backend()
    end subroutine create_cg_gpu_backend

    module subroutine save_static_field_bin()
    end subroutine save_static_field_bin

    module subroutine save_field_at_receivers()
    end subroutine save_field_at_receivers

    module subroutine compute_forces_static(displ,is_nonlinear,kdotu,stress_tensor,strain_tensor)
      use specfem_par, only: NGLOB_AB
      real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB), intent(in) :: displ
      logical, intent(in) :: is_nonlinear
      real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB), intent(out) :: kdotu
      real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), intent(out) :: stress_tensor,strain_tensor
    end subroutine compute_forces_static

    module subroutine compute_forces_phase(displ,iphase,is_nonlinear,kdotu,stress_tensor,strain_tensor)
      use specfem_par, only: NGLOB_AB
      real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB), intent(in) :: displ
      real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB), intent(inout) :: kdotu
      integer, intent(in) :: iphase
      logical, intent(in) :: is_nonlinear
      real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), intent(out) :: stress_tensor,strain_tensor
    end subroutine compute_forces_phase

    module subroutine compute_elemwise_Kxu(ispec,dummyx_loc,dummyy_loc,dummyz_loc,&
                                          force_x,force_y,force_z,&
                                          is_nonlinear,compute_stress_and_strain,&
                                          stress_loc,strain_loc)
      use constants, only: NGLLX,NGLLY,NGLLZ
      integer, intent(in) :: ispec
      real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ), intent(in) :: dummyx_loc,dummyy_loc,dummyz_loc
      real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ), intent(out) :: force_x,force_y,force_z
      real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,6), intent(out) :: stress_loc,strain_loc
      logical, intent(in) :: is_nonlinear
      logical, intent(in) :: compute_stress_and_strain
    end subroutine compute_elemwise_Kxu

    pure module subroutine second2firstPK_stress( &
      sigma_xx,sigma_xy,sigma_xz, &
      sigma_yx,sigma_yy,sigma_yz, &
      sigma_zx,sigma_zy,sigma_zz, &
      duxdxl,duxdyl,duxdzl, &
      duydxl,duydyl,duydzl, &
      duzdxl,duzdyl,duzdzl)
      use constants, only: CUSTOM_REAL
      real(kind=CUSTOM_REAL), intent(inout) :: sigma_xx,sigma_xy,sigma_xz
      real(kind=CUSTOM_REAL), intent(inout) :: sigma_yx,sigma_yy,sigma_yz
      real(kind=CUSTOM_REAL), intent(inout) :: sigma_zx,sigma_zy,sigma_zz
      real(kind=CUSTOM_REAL), intent(in) :: duxdxl,duxdyl,duxdzl
      real(kind=CUSTOM_REAL), intent(in) :: duydxl,duydyl,duydzl
      real(kind=CUSTOM_REAL), intent(in) :: duzdxl,duzdyl,duzdzl
    end subroutine second2firstPK_stress

    pure module subroutine cross_product(x1,x2,x3,y1,y2,y3,z1,z2,z3)
      real(kind=CUSTOM_REAL), intent(in) :: x1,x2,x3
      real(kind=CUSTOM_REAL), intent(in) :: y1,y2,y3
      real(kind=CUSTOM_REAL), intent(out) :: z1,z2,z3
    end subroutine cross_product

    pure module subroutine enforce_roller_bc(vec)
      use specfem_par,only: NGLOB_AB
      real(kind=dp),intent(inout) :: vec(NDIM,NGLOB_AB)
    end subroutine enforce_roller_bc

    pure module subroutine enforce_fixed_bc(vec)
      use specfem_par,only: NGLOB_AB
      real(kind=dp),intent(inout) :: vec(NDIM,NGLOB_AB)
    end subroutine enforce_fixed_bc

    module function parallel_inner_product(vec1,vec2) result(psum)
      use specfem_par, only: NGLOB_AB
      real(kind=dp), dimension(NDIM, NGLOB_AB), intent(in) :: vec1, vec2
      real(kind=dp) :: psum
    end function parallel_inner_product
  end interface


  ! GLOBAL variable to hold static solver class
  type(static_solver_class) :: ssol

  contains

  !> main subroutine for static solution using PETSc, this will be called when ssol%execute() is called, which is set to solve_static_problem_petsc in the type definition of static_solver_class, the actual
  subroutine solve_static_problem_petsc(this)
    use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM, NGLLSQUARE
    use specfem_par, only: NSPEC_AB,ibool,NGLOB_AB,myrank
    use specfem_par_elastic, only : ispec_is_elastic,displ

    use specfem_par, only: fixed_bdry_ijk, fixed_bdry_ispec,num_fixed_bdry_faces
    use specfem_par, only: num_roller_bdry_faces,roller_bdry_ijk, &
                            roller_bdry_ispec, roller_bdry_normal

    use petsc_interfaces, only: fill_mat_petsc
    use petsc_interfaces,only: solve_petsc,extract_petsc

    implicit none

    integer,parameter :: NGLL3 = NGLLX*NGLLY*NGLLZ
    class(static_solver_class), intent(inout) :: this


    real(kind=dp), pointer :: Kloc(:,:,:,:)
    integer :: ibool0_loc(NGLLX,NGLLY,NGLLZ), ispec
    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: ux,uy,uz
    integer :: i, j, k, p, q, r, node_in, node_out
    integer :: ib,dim_in, dim_out, iface,igll2
    real(kind=CUSTOM_REAL) :: force_x(NGLLX, NGLLY, NGLLZ), stress_loc(NGLLX, NGLLY, NGLLZ, 6)
    real(kind=CUSTOM_REAL) :: force_y(NGLLX, NGLLY, NGLLZ), strain_loc(NGLLX, NGLLY, NGLLZ, 6)
    real(kind=CUSTOM_REAL) :: force_z(NGLLX, NGLLY, NGLLZ)
    real(kind=dp) :: Pmat(NDIM,NDIM), norm_vec(NDIM), Ktemp(NDIM,NDIM)
    real(kind=dp), allocatable :: displ_dp(:,:)
    real(kind=CUSTOM_REAL), allocatable :: kdotu(:,:)
    real(kind=dp),target, allocatable :: coo_v(:,:,:,:,:) ! shape(NDIM,NGLL3,NDIM,NGLL3,NSPEC_AB),coo values for stiffness matrix, used for PETSc assembly

    ! boundary arrays
    integer,allocatable :: fixed_bdry_faces(:,:), roller_bdry_faces(:,:)

    ! allocate space
    allocate(coo_v(NDIM,NGLL3,NDIM,NGLL3,NSPEC_AB), &
            displ_dp(NDIM,NGLOB_AB), kdotu(NDIM,NGLOB_AB),&
            fixed_bdry_faces(6,NSPEC_AB), &
            roller_bdry_faces(6,NSPEC_AB))

    ! init boundary arrays
    fixed_bdry_faces(:,:) = 0
    roller_bdry_faces(:,:) = 0
    do i = 1, num_fixed_bdry_faces
      ispec = fixed_bdry_ispec(i)
      do j=1,6
        if(fixed_bdry_faces(j, ispec) == 0) then
          fixed_bdry_faces(j,ispec) = i ! mark the element with fixed boundary condition
          exit
        endif
      enddo
    enddo
    do i = 1, num_roller_bdry_faces
      ispec = roller_bdry_ispec(i)
      do j=1,6
        if(roller_bdry_faces(j, ispec) == 0) then
          roller_bdry_faces(j,ispec) = i ! mark the element with roller boundary condition
          exit
        endif
      enddo
    enddo

    if(myrank == 0) then
      write(*,*)
      write(*,*) '--------------------------------------------------'
      write(*,*) 'element-wise computation of stiffness matrix K ...'
      write(*,*) '--------------------------------------------------'
    endif

    do ispec = 1, NSPEC_AB
      if(.not. ispec_is_elastic(ispec)) cycle
      ibool0_loc(:,:,:) = ibool(:,:,:,ispec) - 1 ! convert to 0-based indexing for C

      ! point to local stiffness matrix for this element
      kloc => coo_v(:,:,:,:,ispec)

      ! compute kloc
      kloc(:,:,:,:) = 0.0_dp
      ux(:,:,:) = 0.0_CUSTOM_REAL; uy(:,:,:) = 0.0_CUSTOM_REAL; uz(:,:,:) = 0.0_CUSTOM_REAL
      do r=1,NGLLZ; do q=1,NGLLY; do p=1,NGLLX
        node_in = (r-1)*NGLLY*NGLLX + (q-1)*NGLLX + p

        ! loop over each direction
        do dim_in = 1,NDIM
          ! apply unit displacement in dim_in direction at node_in
          if (dim_in == 1) then
            ux(p,q,r) = 1.0_CUSTOM_REAL
          else if (dim_in == 2) then
            uy(p,q,r) = 1.0_CUSTOM_REAL
          else
            uz(p,q,r) = 1.0_CUSTOM_REAL
          endif

          ! compute force response at all nodes, store in force_x/y/z
          call compute_elemwise_Kxu(ispec, ux, uy, uz, &
                                    force_x, force_y, force_z,&
                                    .false., .false.,&
                                    stress_loc, strain_loc)
          ! fill kloc based on force response
          do k=1,NGLLZ; do j=1,NGLLY; do i=1,NGLLX
            node_out = (k-1)*NGLLY*NGLLX + (j-1)*NGLLX + i

            kloc(1,node_out,dim_in,node_in) = force_x(i,j,k)
            kloc(2,node_out,dim_in,node_in) = force_y(i,j,k)
            kloc(3,node_out,dim_in,node_in) = force_z(i,j,k)
          enddo; enddo; enddo

          ! reset ux/uy/uz for next iteration
          ux(p,q,r) = 0.0_CUSTOM_REAL
          uy(p,q,r) = 0.0_CUSTOM_REAL
          uz(p,q,r) = 0.0_CUSTOM_REAL
        enddo ! dim_in

      enddo; enddo; enddo

      ! apply boundary conditions to kloc, for fixed boundary condition, we can simply zero out the corresponding rows and columns in kloc,
      ! and set the diagonal entry to 1, for roller boundary condition, we need to remove the normal component of the force response,
      ! which is equivalent to zeroing out the entries in kloc that correspond to the normal direction
      !DEBUG
      do ib = 1,6
      !do ib = 1,0 ! do nothing --- IGNORE ---
        iface = fixed_bdry_faces(ib,ispec)
        if(iface == 0) cycle ! no more fixed boundary condition for this element
        do igll2 = 1, NGLLSQUARE
          i = fixed_bdry_ijk(1,igll2,iface)
          j = fixed_bdry_ijk(2,igll2,iface)
          k = fixed_bdry_ijk(3,igll2,iface)
          node_out = (k-1)*NGLLY*NGLLX + (j-1)*NGLLX + i

          ! zero out rows and columns in kloc
          kloc(:,node_out,:,:)= 0.0_dp
          kloc(:,:,:,node_out) = 0.0_dp

          ! set diagonal entry to 1
          kloc(1,node_out,1,node_out) = 1.0_dp
          kloc(2,node_out,2,node_out) = 1.0_dp
          kloc(3,node_out,3,node_out) = 1.0_dp
        enddo
      enddo

      ! ROLLER BOUNDARY CONDITION
      do ib = 1,6
        iface = roller_bdry_faces(ib,ispec)
        if(iface == 0) cycle ! no more roller boundary condition for this element
        do igll2 = 1, NGLLSQUARE
          i = roller_bdry_ijk(1,igll2,iface)
          j = roller_bdry_ijk(2,igll2,iface)
          k = roller_bdry_ijk(3,igll2,iface)
          node_out = (k-1)*NGLLY*NGLLX + (j-1)*NGLLX + i

          ! copy normal vector for this node
          norm_vec(:) = roller_bdry_normal(:,igll2,iface)
          Pmat = 0.0_dp
          do dim_in = 1,NDIM
            Pmat(dim_in,dim_in) = 1.0_dp
            do dim_out = 1,NDIM
              Pmat(dim_out,dim_in) = Pmat(dim_out,dim_in) - norm_vec(dim_out) * norm_vec(dim_in)
            enddo
          enddo

          ! project K_row to Pmat * K
          do node_in = 1,NGLL3
            Ktemp(:,:) = Kloc(:,node_out,:,node_in)
            Kloc(:,node_out,:,node_in) = matmul(Pmat, Ktemp)
          enddo

          ! K_col
          do node_in = 1,NGLL3
            Ktemp(:,:) = Kloc(:,node_in,:,node_out)
            Kloc(:,node_in,:,node_out) = matmul(Ktemp, Pmat)
          enddo

          ! handle null space
          do dim_in = 1,NDIM
            do dim_out = 1,NDIM
              Kloc(dim_out,node_out,dim_in,node_out) = Kloc(dim_out,node_out,dim_in,node_out) + norm_vec(dim_out) * norm_vec(dim_in)
            enddo
          enddo

        enddo ! igll2
      enddo
    enddo

    if(myrank == 0) then
      write(*,*)
      write(*,*) '----------------------------------------'
      write(*,*) 'global assembly ...'
      write(*,*) '----------------------------------------'
    endif

    ! assemble global matrix in PETSc
    call fill_mat_petsc(this%petsc_ptr, coo_v)

    call synchronize_all();

    if(myrank == 0) then
      write(*,*)
      write(*,*) '----------------------------------------'
      write(*,*) 'run linear solver ...'
      write(*,*) '----------------------------------------'
    endif

    ! solve the linear system using PETSc
    call solve_petsc(this%petsc_ptr)

    ! copy results to displ_dp for later use, note that PETSc uses 0-based indexing while Fortran uses 1-based indexing, so we need to add 1 to the indices when copying
    call extract_petsc(this%petsc_ptr, displ_dp)

    displ(:,:) = real(displ_dp(:,:), kind=CUSTOM_REAL)

    ! compute stress and strain based on the solution, and store in ssol%stress and ssol%strain for later use
    call compute_forces_static(displ,.false., kdotu, ssol%stress, ssol%strain)

    ! free allocated arrays
    deallocate(coo_v, displ_dp, kdotu, fixed_bdry_faces, roller_bdry_faces)

  end subroutine solve_static_problem_petsc

end module static_module
