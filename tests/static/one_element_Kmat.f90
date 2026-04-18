module one_element_static_support
  use iso_fortran_env, only: real64
  implicit none

  integer, parameter :: rk = real64
  integer, parameter :: NDIM = 3
  integer, parameter :: NGLLX = 5
  integer, parameter :: NGLLY = 5
  integer, parameter :: NGLLZ = 5
  integer, parameter :: NSPEC = 1
  integer, parameter :: NGLLCUBE = NGLLX * NGLLY * NGLLZ
  integer, parameter :: m1 = NGLLX
  integer, parameter :: m2 = NGLLX * NGLLY
  real(rk), parameter :: ONE_THIRD = 1.0_rk / 3.0_rk
  real(rk), parameter :: FOUR_THIRDS = 4.0_rk / 3.0_rk
  real(rk), parameter :: TOL_ZERO = 1.0e-10_rk
  real(rk), parameter :: TOL_MATCH = 5.0e-11_rk

    integer, parameter :: VOIGT_LOOKUP(3, 3) = reshape([ &
      0, 5, 4, &
      5, 1, 3, &
      4, 3, 2 ], shape=[3, 3])

  logical :: ANISOTROPY = .false.
  logical :: ROTATION = .false.

  real(rk), dimension(NGLLX) :: xigll, wxgll
  real(rk), dimension(NGLLY) :: yigll, wygll
  real(rk), dimension(NGLLZ) :: zigll, wzgll
  real(rk), dimension(NGLLX, NGLLX) :: hprime_xx, hprime_xxT
  real(rk), dimension(NGLLX, NGLLX) :: hprimewgll_xx, hprimewgll_xxT
  real(rk), dimension(NGLLX, NGLLY, NGLLZ) :: wgllwgll_xy_3D
  real(rk), dimension(NGLLX, NGLLY, NGLLZ) :: wgllwgll_xz_3D
  real(rk), dimension(NGLLX, NGLLY, NGLLZ) :: wgllwgll_yz_3D

  real(rk), dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: xixstore, xiystore, xizstore
  real(rk), dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: etaxstore, etaystore, etazstore
  real(rk), dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: gammaxstore, gammaystore, gammazstore
  real(rk), dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: jacobianstore
  real(rk), dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: kappastore, mustore, rhostore
  real(rk), dimension(NGLLX, NGLLY, NGLLZ) :: xcoord, ycoord, zcoord

  abstract interface
    subroutine displacement_filler(ux_out, uy_out, uz_out)
      import :: rk, NGLLX, NGLLY, NGLLZ
      real(rk), intent(out) :: ux_out(NGLLX, NGLLY, NGLLZ)
      real(rk), intent(out) :: uy_out(NGLLX, NGLLY, NGLLZ)
      real(rk), intent(out) :: uz_out(NGLLX, NGLLY, NGLLZ)
    end subroutine displacement_filler
  end interface

contains

  subroutine initialize_state()
    integer :: i, j, k
    real(rk) :: lx, ly, lz
    real(rk) :: jacobian

    call initialize_gll_data()

    lx = 2.0_rk
    ly = 1.5_rk
    lz = 0.75_rk
    jacobian = lx * ly * lz / 8.0_rk

    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          xcoord(i, j, k) = 0.5_rk * lx * (xigll(i) + 1.0_rk)
          ycoord(i, j, k) = 0.5_rk * ly * (yigll(j) + 1.0_rk)
          zcoord(i, j, k) = 0.5_rk * lz * (zigll(k) + 1.0_rk)

          xixstore(i, j, k, 1) = 2.0_rk / lx
          xiystore(i, j, k, 1) = 0.0_rk
          xizstore(i, j, k, 1) = 0.0_rk
          etaxstore(i, j, k, 1) = 0.0_rk
          etaystore(i, j, k, 1) = 2.0_rk / ly
          etazstore(i, j, k, 1) = 0.0_rk
          gammaxstore(i, j, k, 1) = 0.0_rk
          gammaystore(i, j, k, 1) = 0.0_rk
          gammazstore(i, j, k, 1) = 2.0_rk / lz
          jacobianstore(i, j, k, 1) = jacobian

          kappastore(i, j, k, 1) = 17.5_rk
          mustore(i, j, k, 1) = 6.25_rk
          rhostore(i, j, k, 1) = 1.0_rk
        enddo
      enddo
    enddo
  end subroutine initialize_state

  subroutine initialize_gll_data()
    integer :: i, j, k

    xigll = [ -1.0_rk, -sqrt(3.0_rk / 7.0_rk), 0.0_rk, sqrt(3.0_rk / 7.0_rk), 1.0_rk ]
    yigll = xigll
    zigll = xigll

    wxgll = [ 0.1_rk, 49.0_rk / 90.0_rk, 32.0_rk / 45.0_rk, 49.0_rk / 90.0_rk, 0.1_rk ]
    wygll = wxgll
    wzgll = wxgll

    call build_lagrange_derivative_matrix(xigll, hprime_xx)
    hprime_xxT = transpose(hprime_xx)

    do i = 1, NGLLX
      hprimewgll_xx(i, :) = hprime_xx(i, :) * wxgll(i)
    enddo
    hprimewgll_xxT = transpose(hprimewgll_xx)

    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          wgllwgll_xy_3D(i, j, k) = wxgll(i) * wygll(j)
          wgllwgll_xz_3D(i, j, k) = wxgll(i) * wzgll(k)
          wgllwgll_yz_3D(i, j, k) = wygll(j) * wzgll(k)
        enddo
      enddo
    enddo
  end subroutine initialize_gll_data

  subroutine build_lagrange_derivative_matrix(nodes, deriv)
    real(rk), intent(in) :: nodes(:)
    real(rk), intent(out) :: deriv(size(nodes), size(nodes))

    integer :: i, j, n
    real(rk) :: bary(size(nodes))

    n = size(nodes)
    bary = 1.0_rk

    do i = 1, n
      do j = 1, n
        if (j /= i) bary(i) = bary(i) / (nodes(i) - nodes(j))
      enddo
    enddo

    do i = 1, n
      do j = 1, n
        if (i /= j) then
          deriv(i, j) = bary(j) / (bary(i) * (nodes(i) - nodes(j)))
        else
          deriv(i, j) = 0.0_rk
        endif
      enddo
      deriv(i, i) = -sum(deriv(i, :))
    enddo
  end subroutine build_lagrange_derivative_matrix

  pure function voigt4(i, j, p, q) result(idx)
    integer, intent(in) :: i, j, p, q
    integer :: idx
    integer :: m0, n0, m, n
    
    m0 = VOIGT_LOOKUP(i, j)
    n0 = VOIGT_LOOKUP(p, q)
    
    m = m0 ! to convert to 0-based indexing for the formula, will convert back to 1-based at the end
    n = n0
    if (m0 > n0) then
      m = n0
      n = m0
    endif
    
    idx = m * 6 + n - (m * (m + 1)) / 2 + 1 ! the +1 at the end is to convert back to 1-based indexing
  end function voigt4

subroutine get_kloc_analytical(kloc,grad_phi)
    real(rk), intent(out) :: kloc(NDIM, NGLLCUBE, NDIM, NGLLCUBE)

    integer :: i, j, k, p, q, r, s, x, y, z
    integer :: node_in, node_out, xyz, v_idx
    real(rk) :: mu0, lambda0, lam2mu, w_jac, sum_val
    real(rk) :: dphi(NDIM), Jinv(NDIM, NDIM)
    
    ! ---------------------------------------------------------
    ! PRECOMPUTED ARRAYS (Fits easily in L1/L2 Cache)
    ! ---------------------------------------------------------
    real(rk) :: c21_int(21, NGLLCUBE)
    real(rk),intent(inout) :: grad_phi(NDIM, NGLLCUBE, NGLLCUBE) 
    integer  :: v4(NDIM, NDIM, NDIM, NDIM)
    
    integer, parameter :: ispec = 1

    ! 1. Precompute all possible Voigt index mappings (81 calls instead of millions)
    do s = 1, NDIM; do r = 1, NDIM; do q = 1, NDIM; do p = 1, NDIM
      v4(p,q,r,s) = voigt4(p,q,r,s)
    enddo; enddo; enddo; enddo

    ! 2. Precompute Material Properties & Basis Gradients at all Integration Points
    do z = 1, NGLLZ; do y = 1, NGLLY; do x = 1, NGLLX
      xyz = x + NGLLX * (y - 1) + NGLLX*NGLLY * (z - 1)
      
      ! Material properties combined with integration weights
      mu0 = mustore(x, y, z, ispec)
      lambda0 = kappastore(x, y, z, ispec) - (2.0_rk/3.0_rk) * mu0
      lam2mu = lambda0 + 2.0_rk * mu0
      w_jac = wxgll(x) * wygll(y) * wzgll(z) * jacobianstore(x,y,z,ispec)

      ! Explicitly map non-zero components of the isotropic c21 tensor
      c21_int(:, xyz) = 0.0_rk
      c21_int(1, xyz) = lam2mu * w_jac
      c21_int(2, xyz) = lambda0 * w_jac
      c21_int(3, xyz) = lambda0 * w_jac
      c21_int(7, xyz) = lam2mu * w_jac
      c21_int(8, xyz) = lambda0 * w_jac
      c21_int(12, xyz) = lam2mu * w_jac
      c21_int(16, xyz) = mu0 * w_jac
      c21_int(19, xyz) = mu0 * w_jac
      c21_int(21, xyz) = mu0 * w_jac
      
      ! Inverse Jacobian at this point
      Jinv(1,1) = xixstore(x,y,z,ispec); Jinv(1,2) = etaxstore(x,y,z,ispec); Jinv(1,3) = gammaxstore(x,y,z,ispec)
      Jinv(2,1) = xiystore(x,y,z,ispec); Jinv(2,2) = etaystore(x,y,z,ispec); Jinv(2,3) = gammaystore(x,y,z,ispec)
      Jinv(3,1) = xizstore(x,y,z,ispec); Jinv(3,2) = etazstore(x,y,z,ispec); Jinv(3,3) = gammazstore(x,y,z,ispec)

      ! Gradients of all shape functions at this integration point
      do k = 1, NGLLZ; do j = 1, NGLLY; do i = 1, NGLLX
        node_out = i + NGLLX * (j - 1) + NGLLX*NGLLY * (k - 1)
        
        dphi(1) = 0.0_rk; dphi(2) = 0.0_rk; dphi(3) = 0.0_rk
        if (y == j .and. z == k) dphi(1) = hprime_xx(x, i)
        if (x == i .and. z == k) dphi(2) = hprime_xx(y, j)
        if (x == i .and. y == j) dphi(3) = hprime_xx(z, k)

        grad_phi(1, node_out, xyz) = Jinv(1,1)*dphi(1) + Jinv(1,2)*dphi(2) + Jinv(1,3)*dphi(3)
        grad_phi(2, node_out, xyz) = Jinv(2,1)*dphi(1) + Jinv(2,2)*dphi(2) + Jinv(2,3)*dphi(3)
        grad_phi(3, node_out, xyz) = Jinv(3,1)*dphi(1) + Jinv(3,2)*dphi(2) + Jinv(3,3)*dphi(3)
      enddo; enddo; enddo
    enddo; enddo; enddo

    ! ---------------------------------------------------------
    ! 3. ASSEMBLE MATRIX (Pure Tensor Contraction)
    ! ---------------------------------------------------------
    kloc = 0.0_rk
    
    ! Loop order matches Fortran column-major memory mapping exactly
    do node_in = 1, NGLLCUBE
      do r = 1, NDIM
        do node_out = 1, NGLLCUBE
          do p = 1, NDIM
            sum_val = 0.0_rk
            
            do xyz = 1, NGLLCUBE
              do s = 1, NDIM
                if (grad_phi(s, node_in, xyz) == 0.0_rk) cycle 
                
                do q = 1, NDIM
                  if (grad_phi(q, node_out, xyz) == 0.0_rk) cycle 
                  
                  v_idx = v4(p,q,r,s)
                  if (c21_int(v_idx, xyz) /= 0.0_rk) then
                    sum_val = sum_val + c21_int(v_idx, xyz) * &
                                        grad_phi(q, node_out, xyz) * &
                                        grad_phi(s, node_in, xyz)
                  endif
                enddo
              enddo
            enddo
            
            kloc(p, node_out, r, node_in) = sum_val
            
          enddo
        enddo
      enddo
    enddo

  end subroutine get_kloc_analytical

  subroutine get_kloc_analytical1(kloc)
    real(rk), intent(out) :: kloc(NDIM, NGLLCUBE, NDIM, NGLLCUBE)

    integer :: a,b,c, i,j,k, p,r,q,s 
    integer:: x,y,z 
    real(rk) :: tempq(3), temps(3)
    real(rk) ::  der_q(NDIM,NDIM),der_s(NDIM,NDIM), c21(21),mu0,lambda0 
    integer:: node_out,node_in 
    integer,parameter :: ispec = 1

    ! outer loops 
    do r=1,NDIM; do c = 1,NGLLZ; do b = 1,NGLLY; do a = 1,NGLLX;
        node_in = a + NGLLX * (b - 1) + NGLLX*NGLLY * (c - 1)
        ! inner loops
        do p=1,NDIM; do k=1,NGLLZ; do j=1,NGLLY; do i=1,NGLLX;
            node_out = i + NGLLX * (j - 1) + NGLLX*NGLLY * (k - 1)
            kloc(p, node_out, r, node_in) = 0.0_rk 

            ! sum loop 
            do z = 1,NGLLZ; do y = 1,NGLLY; do x = 1,NGLLX;
              mu0 = mustore(x, y, z,ispec)
              lambda0 = kappastore(x, y, z,ispec) - 2.0_rk/3.0_rk * mu0
              c21(:) = [ lambda0+2.d0*mu0, lambda0, lambda0, 0.d0, 0.d0, 0.d0, lambda0+2.d0*mu0, lambda0, 0.d0, 0.d0, 0.d0, lambda0+2.d0*mu0, 0.d0, 0.d0, 0.d0, mu0, 0.d0, 0.d0, mu0, 0.d0, mu0 ]
              do q=1,NDIM
                if(q == 1) then ! -> x
                  der_q(q,1) = xixstore(x,y,z,ispec); der_q(q,2) = etaxstore(x,y,z,ispec); der_q(1,3) = gammaxstore(x,y,z,ispec)
                else if(q == 2) then ! -> y
                  der_q(q,1) = xiystore(x,y,z,ispec); der_q(q,2) = etaystore(x,y,z,ispec); der_q(2,3) = gammaystore(x,y,z,ispec)
                else ! -> z
                  der_q(q,1) = xizstore(x,y,z,ispec); der_q(q,2) = etazstore(x,y,z,ispec); der_q(3,3) = gammazstore(x,y,z,ispec)
                endif 
                tempq(:) = 0.0_rk
                if(y == j .and. z == k) tempq(1) = hprime_xx(x,i) * der_q(q,1)
                if(x == i .and. z == k) tempq(2) = hprime_xx(y,j) * der_q(q,2)
                if(x == i .and. y == j) tempq(3) = hprime_xx(z,k) * der_q(q,3)
                do s=1,NDIM;
                  if(s == 1) then ! -> x
                    der_s(s,1) = xixstore(x,y,z,ispec); der_s(s,2) = etaxstore(x,y,z,ispec); der_s(1,3) = gammaxstore(x,y,z,ispec)
                  else if(s == 2) then ! -> y
                    der_s(s,1) = xiystore(x,y,z,ispec); der_s(s,2) = etaystore(x,y,z,ispec); der_s(2,3) = gammaystore(x,y,z,ispec)
                  else ! -> z
                    der_s(s,1) = xizstore(x,y,z,ispec); der_s(s,2) = etazstore(x,y,z,ispec); der_s(3,3) = gammazstore(x,y,z,ispec)
                  endif
                  temps(:) = 0.0_rk
                  if(y == b .and. z == c) temps(1) = hprime_xx(x,a) * der_s(s,1)
                  if(x == a .and. z == c) temps(2) = hprime_xx(y,b) * der_s(s,2)
                  if(x == a .and. y == b) temps(3) = hprime_xx(z,c) * der_s(s,3)

                  kloc(p, node_out, r, node_in) = kloc(p, node_out, r, node_in) + &
                    c21(voigt4(p,q,r,s)) * wxgll(x) * wygll(y) * wzgll(z) * &
                    jacobianstore(x,y,z,ispec) * &
                    sum(tempq(:)) * sum(temps(:))
                enddo 
              enddo 
            enddo; enddo; enddo;
        enddo; enddo; enddo; enddo;
    enddo; enddo; enddo; enddo;
  end subroutine get_kloc_analytical1

  subroutine assemble_kloc(kloc)
    real(rk), intent(out) :: kloc(NDIM, NGLLCUBE, NDIM, NGLLCUBE)

    integer :: dim_in, i, j, k, node_in, node_out, p, q, r
    real(rk) :: ux(NGLLX, NGLLY, NGLLZ)
    real(rk) :: uy(NGLLX, NGLLY, NGLLZ)
    real(rk) :: uz(NGLLX, NGLLY, NGLLZ)
    real(rk) :: force_x(NGLLX, NGLLY, NGLLZ)
    real(rk) :: force_y(NGLLX, NGLLY, NGLLZ)
    real(rk) :: force_z(NGLLX, NGLLY, NGLLZ)
    real(rk) :: stress_loc(NGLLX, NGLLY, NGLLZ, 6)
    real(rk) :: strain_loc(NGLLX, NGLLY, NGLLZ, 6)

    kloc = 0.0_rk
    ux = 0.0_rk
    uy = 0.0_rk
    uz = 0.0_rk

    do r = 1, NGLLZ
      do q = 1, NGLLY
        do p = 1, NGLLX
          node_in = flatten_node(p, q, r)
          do dim_in = 1, NDIM
            if (dim_in == 1) then
              ux(p, q, r) = 1.0_rk
            else if (dim_in == 2) then
              uy(p, q, r) = 1.0_rk
            else
              uz(p, q, r) = 1.0_rk
            endif

            call compute_elemwise_Kxu(1, ux, uy, uz, force_x, force_y, force_z, .false., .false., stress_loc, strain_loc)

            do k = 1, NGLLZ
              do j = 1, NGLLY
                do i = 1, NGLLX
                  node_out = flatten_node(i, j, k)
                  kloc(1, node_out, dim_in, node_in) = force_x(i, j, k)
                  kloc(2, node_out, dim_in, node_in) = force_y(i, j, k)
                  kloc(3, node_out, dim_in, node_in) = force_z(i, j, k)
                enddo
              enddo
            enddo

            ux(p, q, r) = 0.0_rk
            uy(p, q, r) = 0.0_rk
            uz(p, q, r) = 0.0_rk
          enddo
        enddo
      enddo
    enddo
  end subroutine assemble_kloc

  subroutine compute_elemwise_Kxu(ispec, dummyx_loc, dummyy_loc, dummyz_loc, force_x, force_y, force_z, &
                                  is_nonlinear, compute_stress_and_strain, stress_loc, strain_loc)
    integer, intent(in) :: ispec
    real(rk), intent(in) :: dummyx_loc(NGLLX, NGLLY, NGLLZ)
    real(rk), intent(in) :: dummyy_loc(NGLLX, NGLLY, NGLLZ)
    real(rk), intent(in) :: dummyz_loc(NGLLX, NGLLY, NGLLZ)
    real(rk), intent(out) :: force_x(NGLLX, NGLLY, NGLLZ)
    real(rk), intent(out) :: force_y(NGLLX, NGLLY, NGLLZ)
    real(rk), intent(out) :: force_z(NGLLX, NGLLY, NGLLZ)
    real(rk), intent(out) :: stress_loc(NGLLX, NGLLY, NGLLZ, 6)
    real(rk), intent(out) :: strain_loc(NGLLX, NGLLY, NGLLZ, 6)
    logical, intent(in) :: is_nonlinear
    logical, intent(in) :: compute_stress_and_strain

    integer :: i, j, k
    real(rk) :: tempx1(NGLLX, NGLLY, NGLLZ), tempy1(NGLLX, NGLLY, NGLLZ), tempz1(NGLLX, NGLLY, NGLLZ)
    real(rk) :: tempx2(NGLLX, NGLLY, NGLLZ), tempy2(NGLLX, NGLLY, NGLLZ), tempz2(NGLLX, NGLLY, NGLLZ)
    real(rk) :: tempx3(NGLLX, NGLLY, NGLLZ), tempy3(NGLLX, NGLLY, NGLLZ), tempz3(NGLLX, NGLLY, NGLLZ)
    real(rk) :: duxdxl(NGLLX, NGLLY, NGLLZ), duxdyl(NGLLX, NGLLY, NGLLZ), duxdzl(NGLLX, NGLLY, NGLLZ)
    real(rk) :: duydxl(NGLLX, NGLLY, NGLLZ), duydyl(NGLLX, NGLLY, NGLLZ), duydzl(NGLLX, NGLLY, NGLLZ)
    real(rk) :: duzdxl(NGLLX, NGLLY, NGLLZ), duzdyl(NGLLX, NGLLY, NGLLZ), duzdzl(NGLLX, NGLLY, NGLLZ)
    real(rk) :: newtempx1(NGLLX, NGLLY, NGLLZ), newtempy1(NGLLX, NGLLY, NGLLZ), newtempz1(NGLLX, NGLLY, NGLLZ)
    real(rk) :: newtempx2(NGLLX, NGLLY, NGLLZ), newtempy2(NGLLX, NGLLY, NGLLZ), newtempz2(NGLLX, NGLLY, NGLLZ)
    real(rk) :: newtempx3(NGLLX, NGLLY, NGLLZ), newtempy3(NGLLX, NGLLY, NGLLZ), newtempz3(NGLLX, NGLLY, NGLLZ)
    real(rk) :: xixl, xiyl, xizl, etaxl, etayl, etazl, gammaxl, gammayl, gammazl, jacobianl
    real(rk) :: duxdyl_plus_duydxl, duzdxl_plus_duxdzl, duzdyl_plus_duydzl
    real(rk) :: sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_xz, sigma_yz, sigma_yx, sigma_zx, sigma_zy
    real(rk) :: lambdal, mul, lambdalplus2mul, kappal, fac1, fac2, fac3
    real(rk) :: strain_xx, strain_yy, strain_zz, strain_xy, strain_xz, strain_yz

    call mxm5_3comp_singleA(hprime_xxT, dummyx_loc, dummyy_loc, dummyz_loc, tempx1, tempy1, tempz1)
    call mxm5_3comp_3dmat_single(dummyx_loc, dummyy_loc, dummyz_loc, hprime_xxT, tempx2, tempy2, tempz2)
    call mxm5_3comp_singleB(dummyx_loc, dummyy_loc, dummyz_loc, hprime_xxT, tempx3, tempy3, tempz3)

    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          xixl = xixstore(i, j, k, ispec)
          xiyl = xiystore(i, j, k, ispec)
          xizl = xizstore(i, j, k, ispec)
          etaxl = etaxstore(i, j, k, ispec)
          etayl = etaystore(i, j, k, ispec)
          etazl = etazstore(i, j, k, ispec)
          gammaxl = gammaxstore(i, j, k, ispec)
          gammayl = gammaystore(i, j, k, ispec)
          gammazl = gammazstore(i, j, k, ispec)

          duxdxl(i, j, k) = xixl * tempx1(i, j, k) + etaxl * tempx2(i, j, k) + gammaxl * tempx3(i, j, k)
          duxdyl(i, j, k) = xiyl * tempx1(i, j, k) + etayl * tempx2(i, j, k) + gammayl * tempx3(i, j, k)
          duxdzl(i, j, k) = xizl * tempx1(i, j, k) + etazl * tempx2(i, j, k) + gammazl * tempx3(i, j, k)
          duydxl(i, j, k) = xixl * tempy1(i, j, k) + etaxl * tempy2(i, j, k) + gammaxl * tempy3(i, j, k)
          duydyl(i, j, k) = xiyl * tempy1(i, j, k) + etayl * tempy2(i, j, k) + gammayl * tempy3(i, j, k)
          duydzl(i, j, k) = xizl * tempy1(i, j, k) + etazl * tempy2(i, j, k) + gammazl * tempy3(i, j, k)
          duzdxl(i, j, k) = xixl * tempz1(i, j, k) + etaxl * tempz2(i, j, k) + gammaxl * tempz3(i, j, k)
          duzdyl(i, j, k) = xiyl * tempz1(i, j, k) + etayl * tempz2(i, j, k) + gammayl * tempz3(i, j, k)
          duzdzl(i, j, k) = xizl * tempz1(i, j, k) + etazl * tempz2(i, j, k) + gammazl * tempz3(i, j, k)
        enddo
      enddo
    enddo

    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          strain_xx = duxdxl(i, j, k)
          strain_yy = duydyl(i, j, k)
          strain_zz = duzdzl(i, j, k)
          strain_xy = 0.5_rk * (duxdyl(i, j, k) + duydxl(i, j, k))
          strain_xz = 0.5_rk * (duzdxl(i, j, k) + duxdzl(i, j, k))
          strain_yz = 0.5_rk * (duzdyl(i, j, k) + duydzl(i, j, k))

          if (is_nonlinear) then
            error stop 'nonlinear mode is not implemented in this standalone example'
          endif

          if (compute_stress_and_strain) then
            strain_loc(i, j, k, 1) = strain_xx
            strain_loc(i, j, k, 2) = strain_yy
            strain_loc(i, j, k, 3) = strain_zz
            strain_loc(i, j, k, 4) = strain_yz
            strain_loc(i, j, k, 5) = strain_xz
            strain_loc(i, j, k, 6) = strain_xy
          endif

          duxdyl_plus_duydxl = 2.0_rk * strain_xy
          duzdxl_plus_duxdzl = 2.0_rk * strain_xz
          duzdyl_plus_duydzl = 2.0_rk * strain_yz

          if (ANISOTROPY) then
            error stop 'anisotropy is not implemented in this standalone example'
          else
            kappal = kappastore(i, j, k, ispec)
            mul = mustore(i, j, k, ispec)
            lambdalplus2mul = kappal + FOUR_THIRDS * mul
            lambdal = lambdalplus2mul - 2.0_rk * mul

            sigma_xx = lambdalplus2mul * strain_xx + lambdal * (strain_yy + strain_zz)
            sigma_yy = lambdalplus2mul * strain_yy + lambdal * (strain_xx + strain_zz)
            sigma_zz = lambdalplus2mul * strain_zz + lambdal * (strain_xx + strain_yy)
            sigma_xy = mul * duxdyl_plus_duydxl
            sigma_xz = mul * duzdxl_plus_duxdzl
            sigma_yz = mul * duzdyl_plus_duydzl
          endif

          if (compute_stress_and_strain) then
            stress_loc(i, j, k, 1) = sigma_xx
            stress_loc(i, j, k, 2) = sigma_yy
            stress_loc(i, j, k, 3) = sigma_zz
            stress_loc(i, j, k, 4) = sigma_yz
            stress_loc(i, j, k, 5) = sigma_xz
            stress_loc(i, j, k, 6) = sigma_xy
          endif

          sigma_yx = sigma_xy
          sigma_zx = sigma_xz
          sigma_zy = sigma_yz

          xixl = xixstore(i, j, k, ispec)
          xiyl = xiystore(i, j, k, ispec)
          xizl = xizstore(i, j, k, ispec)
          etaxl = etaxstore(i, j, k, ispec)
          etayl = etaystore(i, j, k, ispec)
          etazl = etazstore(i, j, k, ispec)
          gammaxl = gammaxstore(i, j, k, ispec)
          gammayl = gammaystore(i, j, k, ispec)
          gammazl = gammazstore(i, j, k, ispec)
          jacobianl = jacobianstore(i, j, k, ispec)

          tempx1(i, j, k) = jacobianl * (sigma_xx * xixl + sigma_yx * xiyl + sigma_zx * xizl)
          tempy1(i, j, k) = jacobianl * (sigma_xy * xixl + sigma_yy * xiyl + sigma_zy * xizl)
          tempz1(i, j, k) = jacobianl * (sigma_xz * xixl + sigma_yz * xiyl + sigma_zz * xizl)
          tempx2(i, j, k) = jacobianl * (sigma_xx * etaxl + sigma_yx * etayl + sigma_zx * etazl)
          tempy2(i, j, k) = jacobianl * (sigma_xy * etaxl + sigma_yy * etayl + sigma_zy * etazl)
          tempz2(i, j, k) = jacobianl * (sigma_xz * etaxl + sigma_yz * etayl + sigma_zz * etazl)
          tempx3(i, j, k) = jacobianl * (sigma_xx * gammaxl + sigma_yx * gammayl + sigma_zx * gammazl)
          tempy3(i, j, k) = jacobianl * (sigma_xy * gammaxl + sigma_yy * gammayl + sigma_zy * gammazl)
          tempz3(i, j, k) = jacobianl * (sigma_xz * gammaxl + sigma_yz * gammayl + sigma_zz * gammazl)
        enddo
      enddo
    enddo

    call mxm5_3comp_singleA(hprimewgll_xx, tempx1, tempy1, tempz1, newtempx1, newtempy1, newtempz1)
    call mxm5_3comp_3dmat_single(tempx2, tempy2, tempz2, hprimewgll_xx, newtempx2, newtempy2, newtempz2)
    call mxm5_3comp_singleB(tempx3, tempy3, tempz3, hprimewgll_xx, newtempx3, newtempy3, newtempz3)

    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          if (ROTATION) error stop 'rotation is not implemented in this standalone example'

          fac1 = wgllwgll_yz_3D(i, j, k)
          fac2 = wgllwgll_xz_3D(i, j, k)
          fac3 = wgllwgll_xy_3D(i, j, k)
          force_x(i, j, k) = fac1 * newtempx1(i, j, k) + fac2 * newtempx2(i, j, k) + fac3 * newtempx3(i, j, k)
          force_y(i, j, k) = fac1 * newtempy1(i, j, k) + fac2 * newtempy2(i, j, k) + fac3 * newtempy3(i, j, k)
          force_z(i, j, k) = fac1 * newtempz1(i, j, k) + fac2 * newtempz2(i, j, k) + fac3 * newtempz3(i, j, k)
        enddo
      enddo
    enddo
  end subroutine compute_elemwise_Kxu

  subroutine mxm5_3comp_singleA(a, b1, b2, b3, c1, c2, c3)
    real(rk), intent(in) :: a(NGLLX, NGLLX)
    real(rk), intent(in) :: b1(NGLLX, NGLLY, NGLLZ), b2(NGLLX, NGLLY, NGLLZ), b3(NGLLX, NGLLY, NGLLZ)
    real(rk), intent(out) :: c1(NGLLX, NGLLY, NGLLZ), c2(NGLLX, NGLLY, NGLLZ), c3(NGLLX, NGLLY, NGLLZ)

    integer :: i, j, k, l

    real(rk) :: temp1, temp2, temp3,a0 
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          temp1 = 0.0_rk
          temp2 = 0.0_rk
          temp3 = 0.0_rk
          do l = 1, NGLLX
            a0 = a(l, i)
            temp1 = temp1 + a0 * b1(l, j, k)
            temp2 = temp2 + a0 * b2(l, j, k)
            temp3 = temp3 + a0 * b3(l, j, k)
          enddo
          c1(i, j, k) = temp1
          c2(i, j, k) = temp2
          c3(i, j, k) = temp3
        enddo
      enddo
    enddo
  end subroutine mxm5_3comp_singleA

  subroutine mxm5_3comp_3dmat_single(a1, a2, a3, b, c1, c2, c3)
    real(rk), intent(in) :: a1(NGLLX, NGLLY, NGLLZ), a2(NGLLX, NGLLY, NGLLZ), a3(NGLLX, NGLLY, NGLLZ)
    real(rk), intent(in) :: b(NGLLY, NGLLY)
    real(rk), intent(out) :: c1(NGLLX, NGLLY, NGLLZ), c2(NGLLX, NGLLY, NGLLZ), c3(NGLLX, NGLLY, NGLLZ)

    integer :: i, j, k,l
    real(rk) :: temp1, temp2, temp3,b0 


    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          temp1 = 0.0_rk
          temp2 = 0.0_rk
          temp3 = 0.0_rk
          do l = 1, NGLLY
            b0 = b(l, j)
            temp1 = temp1 + a1(i, l, k) * b0
            temp2 = temp2 + a2(i, l, k) * b0
            temp3 = temp3 + a3(i, l, k) * b0
          enddo
          c1(i, j, k) = temp1
          c2(i, j, k) = temp2
          c3(i, j, k) = temp3
        enddo
      enddo
    enddo
  end subroutine mxm5_3comp_3dmat_single

  subroutine mxm5_3comp_singleB(a1, a2, a3, b, c1, c2, c3)
    real(rk), intent(in) :: a1(NGLLX, NGLLY, NGLLZ), a2(NGLLX, NGLLY, NGLLZ), a3(NGLLX, NGLLY, NGLLZ)
    real(rk), intent(in) :: b(NGLLZ, NGLLZ)
    real(rk), intent(out) :: c1(NGLLX, NGLLY, NGLLZ), c2(NGLLX, NGLLY, NGLLZ), c3(NGLLX, NGLLY, NGLLZ)

    integer :: i, j, k,l
    real(rk) :: temp1, temp2, temp3,b0
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          temp1 = 0.0_rk
          temp2 = 0.0_rk
          temp3 = 0.0_rk
          do l = 1, NGLLZ
            b0 = b(l, k)
            temp1 = temp1 + a1(i, j, l) * b0
            temp2 = temp2 + a2(i, j, l) * b0
            temp3 = temp3 + a3(i, j, l) * b0
          enddo
          c1(i, j, k) = temp1
          c2(i, j, k) = temp2
          c3(i, j, k) = temp3
        enddo
      enddo
    enddo
  end subroutine mxm5_3comp_singleB

  subroutine apply_kloc(kloc, ux, uy, uz, force_x, force_y, force_z)
    real(rk), intent(in) :: kloc(NDIM, NGLLCUBE, NDIM, NGLLCUBE)
    real(rk), intent(in) :: ux(NGLLX, NGLLY, NGLLZ)
    real(rk), intent(in) :: uy(NGLLX, NGLLY, NGLLZ)
    real(rk), intent(in) :: uz(NGLLX, NGLLY, NGLLZ)
    real(rk), intent(out) :: force_x(NGLLX, NGLLY, NGLLZ)
    real(rk), intent(out) :: force_y(NGLLX, NGLLY, NGLLZ)
    real(rk), intent(out) :: force_z(NGLLX, NGLLY, NGLLZ)

    integer :: i, j, k, node_in, node_out
    real(rk) :: uvec(NDIM, NGLLCUBE)
    real(rk) :: fvec(NDIM, NGLLCUBE)

    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          node_in = flatten_node(i, j, k)
          uvec(1, node_in) = ux(i, j, k)
          uvec(2, node_in) = uy(i, j, k)
          uvec(3, node_in) = uz(i, j, k)
        enddo
      enddo
    enddo

    fvec = 0.0_rk
    do node_out = 1, NGLLCUBE
      do node_in = 1, NGLLCUBE
        fvec(:, node_out) = fvec(:, node_out) + matmul(kloc(:, node_out, :, node_in), uvec(:, node_in))
      enddo
    enddo

    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          node_out = flatten_node(i, j, k)
          force_x(i, j, k) = fvec(1, node_out)
          force_y(i, j, k) = fvec(2, node_out)
          force_z(i, j, k) = fvec(3, node_out)
        enddo
      enddo
    enddo
  end subroutine apply_kloc

  subroutine fill_rigid_translation(ux, uy, uz)
    real(rk), intent(out) :: ux(NGLLX, NGLLY, NGLLZ), uy(NGLLX, NGLLY, NGLLZ), uz(NGLLX, NGLLY, NGLLZ)

    ux = 1.25_rk
    uy = -0.50_rk
    uz = 0.75_rk
  end subroutine fill_rigid_translation

  subroutine fill_rigid_rotation_z(ux, uy, uz)
    real(rk), intent(out) :: ux(NGLLX, NGLLY, NGLLZ), uy(NGLLX, NGLLY, NGLLZ), uz(NGLLX, NGLLY, NGLLZ)

    ux = -ycoord
    uy = xcoord
    uz = 0.0_rk
  end subroutine fill_rigid_rotation_z

  subroutine fill_general_field(ux, uy, uz)
    real(rk), intent(out) :: ux(NGLLX, NGLLY, NGLLZ), uy(NGLLX, NGLLY, NGLLZ), uz(NGLLX, NGLLY, NGLLZ)

    ux = 0.30_rk * xcoord + 0.10_rk * ycoord * zcoord
    uy = -0.20_rk * ycoord + 0.05_rk * xcoord * zcoord
    uz = 0.40_rk * zcoord + 0.07_rk * xcoord * ycoord
  end subroutine fill_general_field

  subroutine run_rigid_mode_check(name, kloc, fill_displacement)
    character(len=*), intent(in) :: name
    real(rk), intent(in) :: kloc(NDIM, NGLLCUBE, NDIM, NGLLCUBE)
    procedure(displacement_filler) :: fill_displacement

    real(rk) :: ux(NGLLX, NGLLY, NGLLZ), uy(NGLLX, NGLLY, NGLLZ), uz(NGLLX, NGLLY, NGLLZ)
    real(rk) :: force_x(NGLLX, NGLLY, NGLLZ), force_y(NGLLX, NGLLY, NGLLZ), force_z(NGLLX, NGLLY, NGLLZ)
    real(rk) :: matrix_x(NGLLX, NGLLY, NGLLZ), matrix_y(NGLLX, NGLLY, NGLLZ), matrix_z(NGLLX, NGLLY, NGLLZ)
    real(rk) :: stress_loc(NGLLX, NGLLY, NGLLZ, 6), strain_loc(NGLLX, NGLLY, NGLLZ, 6)
    real(rk) :: max_force, max_matrix_force

    call fill_displacement(ux, uy, uz)
    call compute_elemwise_Kxu(1, ux, uy, uz, force_x, force_y, force_z, .false., .true., stress_loc, strain_loc)
    call apply_kloc(kloc, ux, uy, uz, matrix_x, matrix_y, matrix_z)

    max_force = max(maxabs3(force_x), max(maxabs3(force_y), maxabs3(force_z)))
    max_matrix_force = max(maxabs3(matrix_x), max(maxabs3(matrix_y), maxabs3(matrix_z)))

    write(*,'(A,1X,A,1X,ES12.5,1X,A,1X,ES12.5)') 'CHECK', trim(name), max_force, 'matrix', max_matrix_force
    if (max_force > TOL_ZERO .or. max_matrix_force > TOL_ZERO) then
      error stop 'rigid mode check failed'
    endif
  end subroutine run_rigid_mode_check

  subroutine run_consistency_check(kloc)
    real(rk), intent(in) :: kloc(NDIM, NGLLCUBE, NDIM, NGLLCUBE)

    real(rk) :: ux(NGLLX, NGLLY, NGLLZ), uy(NGLLX, NGLLY, NGLLZ), uz(NGLLX, NGLLY, NGLLZ)
    real(rk) :: direct_x(NGLLX, NGLLY, NGLLZ), direct_y(NGLLX, NGLLY, NGLLZ), direct_z(NGLLX, NGLLY, NGLLZ)
    real(rk) :: matrix_x(NGLLX, NGLLY, NGLLZ), matrix_y(NGLLX, NGLLY, NGLLZ), matrix_z(NGLLX, NGLLY, NGLLZ)
    real(rk) :: stress_loc(NGLLX, NGLLY, NGLLZ, 6), strain_loc(NGLLX, NGLLY, NGLLZ, 6)
    real(rk) :: max_diff

    call fill_general_field(ux, uy, uz)
    call compute_elemwise_Kxu(1, ux, uy, uz, direct_x, direct_y, direct_z, .false., .true., stress_loc, strain_loc)
    call apply_kloc(kloc, ux, uy, uz, matrix_x, matrix_y, matrix_z)

    max_diff = max(maxabs3(direct_x - matrix_x), max(maxabs3(direct_y - matrix_y), maxabs3(direct_z - matrix_z)))

    write(*,'(A,1X,A,1X,ES12.5)') 'CHECK', 'assembly_consistency', max_diff
    if (max_diff > TOL_MATCH) then
      error stop 'assembled matrix does not match direct kernel application'
    endif
  end subroutine run_consistency_check

  function maxabs3(field) result(value)
    real(rk), intent(in) :: field(NGLLX, NGLLY, NGLLZ)
    real(rk) :: value

    value = maxval(abs(field))
  end function maxabs3

  integer function flatten_node(i, j, k) result(node)
    integer, intent(in) :: i, j, k

    node = (k - 1) * NGLLY * NGLLX + (j - 1) * NGLLX + i
  end function flatten_node

end module one_element_static_support

program one_element_static_examples
  use one_element_static_support
  implicit none

  real(rk) :: tic,toc 

  real(rk) :: kloc(NDIM, NGLLCUBE, NDIM, NGLLCUBE)
  real(rk) :: kloc1(NDIM, NGLLCUBE, NDIM, NGLLCUBE)
  real(rk),allocatable :: grad_phi(:,:,:)

  allocate( grad_phi(NDIM, NGLLCUBE, NGLLCUBE) )

  call initialize_state()

  call cpu_time(tic)
  call assemble_kloc(kloc)
  call cpu_time(toc)
  print*, 'Time taken to assemble kloc:', toc - tic, 'seconds'

  call cpu_time(tic)
  call get_kloc_analytical(kloc1, grad_phi)
  call cpu_time(toc)
  print*, 'Time taken to compute analytical kloc:', toc - tic, 'seconds'

  print*, 'Max difference between assembled kloc and analytical kloc:', maxval(abs(kloc - kloc1))


  call run_rigid_mode_check('rigid_translation', kloc, fill_rigid_translation)
  call run_rigid_mode_check('rigid_rotation_z', kloc, fill_rigid_rotation_z)
  call run_consistency_check(kloc)

  write(*,'(A)') 'All one-element checks passed.'
end program one_element_static_examples