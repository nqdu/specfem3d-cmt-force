#include "mesh_constants_gpu.h"

// GPU port of the preconditioned CG solver in static_problem_impl() (static_module.f90),
// including compute_forces_static_gpu() (the element-wise K*u matrix-vector product, mirroring
// compute_forces_static()/compute_elemwise_Kxu() in static_module.f90). The nonlinear
// (finite-strain / second Piola-Kirchhoff) path is not ported -- static_problem_impl() only
// ever calls compute_forces_static() with is_nonlinear = .false., so compute_forces_static_gpu()
// hardcodes the linear-strain branch.
//
// MPI assembly of kdotu across partition boundaries mirrors how the dynamic solver's GPU path
// does it for d_accel (compute_forces_viscoelastic_calling_routine.F90 + assemble_MPI_vector_cuda.cu):
// gather only the interface-boundary values into a compact device buffer (cg_prepare_boundary_kdotu),
// copy just that (mp->size_mpi_buffer reals, not the whole NGLOB_AB array) to host, hand it to
// Fortran's assemble_MPI_vector_send_cuda() (send/recv on the packed buffers only, no full-array
// packing), then on the way back scatter the received buffer into kdotu (cg_assemble_boundary_kdotu).
// The prepare/assemble kernels used for d_accel (assemble_MPI_vector_cuda.cu) can't be reused
// directly since they're hardcoded to d_accel/d_b_accel and live in a separate device-link group
// from this file, so this file has its own copies operating on kdotu instead.

typedef struct {

  // double precision PCG vectors, shape(NDIM,NGLOB_AB)
  double *r, *p, *Ap, *z, *inv_pred, *u;

  // single precision (realw) buffers for the element-wise K*p matrix-vector product, shape(NDIM,NGLOB_AB)
  realw *kdotu, *p_cr, *Ap_cr, *z_cr;

  // external force vector, shape(NDIM,NGLOB_AB)
  realw *force_ext;

  // inverse multiplicity weight used in the parallel inner product, shape(NGLOB_AB)
  realw *inv_mult;

  // fixed (Dirichlet, vec = 0) boundary condition
  int num_fixed_bdry_faces;
  int *fixed_bdry_ispec;   // shape(num_fixed_bdry_faces)
  int *fixed_bdry_ijk;     // shape(NDIM,NGLL2,num_fixed_bdry_faces)

  // roller (remove normal component) boundary condition
  int num_roller_bdry_faces;
  int *roller_bdry_ispec;    // shape(num_roller_bdry_faces)
  int *roller_bdry_ijk;      // shape(NDIM,NGLL2,num_roller_bdry_faces)
  realw *roller_bdry_normal; // shape(NDIM,NGLL2,num_roller_bdry_faces)

  // stress/strain output of compute_forces_static_gpu(), voigt notation
  // (xx,yy,zz,yz,xz,xy), shape(NGLL3,NSPEC_AB,6) -- mirrors ssol%stress/ssol%strain
  realw *stress, *strain;

  int NGLOB_AB;
  int NSPEC_AB;
  int nspec_inner_elastic, nspec_outer_elastic;
  int myrank;

  // MPI vector assembly of kdotu. num_interfaces_ext_mesh/max_nibool_interfaces_ext_mesh/
  // d_nibool_interfaces_ext_mesh/d_ibool_interfaces_ext_mesh/size_mpi_buffer are read straight
  // off Mesh* at call time (already there for d_accel's own assembly) -- only the host-only
  // pieces (nibool as a host array, neighbor ranks, the packed buffers, the request handles)
  // need to be stashed here, as raw pointers into Fortran (specfem_par) arrays that live for
  // the whole run.
  int NPROC;
  int *nibool_interfaces_ext_mesh, *my_neighbors_ext_mesh; // host arrays
  realw *buffer_send_vector_ext_mesh, *buffer_recv_vector_ext_mesh; // host arrays (Fortran-owned)
  int *request_send_vector_ext_mesh, *request_recv_vector_ext_mesh; // host arrays (Fortran-owned)
  realw *d_send_buffer, *d_recv_buffer; // device-side compact boundary buffers, size mp->size_mpi_buffer

} CGSolver;


/* ----------------------------------------------------------------------------------------------- */

// Fortran subroutines called directly from C: assemble_MPI_vector_send_cuda() sends/receives the
// already-packed boundary buffers only (src/specfem3D/assemble_MPI_vector.f90); wait_req() waits
// on a single MPI request handle (src/shared/parallel.f90)

/* ----------------------------------------------------------------------------------------------- */

extern "C" void FC_FUNC_(assemble_mpi_vector_send_cuda,ASSEMBLE_MPI_VECTOR_SEND_CUDA)
    (int* NPROC,
     realw* buffer_send_vector_ext_mesh, realw* buffer_recv_vector_ext_mesh,
     int* num_interfaces_ext_mesh, int* max_nibool_interfaces_ext_mesh,
     int* nibool_interfaces_ext_mesh, int* my_neighbors_ext_mesh,
     int* request_send_vector_ext_mesh, int* request_recv_vector_ext_mesh);

extern "C" void FC_FUNC_(wait_req,WAIT_REQ)(int* req);


/* ----------------------------------------------------------------------------------------------- */

// small helpers missing from mesh_constants_wrapper.h (which only has int/realw/field variants)

/* ----------------------------------------------------------------------------------------------- */

static inline void gpuMalloc_double(void** d_array_addr_ptr,const size_t size){
#ifdef USE_CUDA
  if (run_cuda){
    print_CUDA_error_if_any(cudaMalloc((void**)d_array_addr_ptr,size*sizeof(double)),2601);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    print_HIP_error_if_any(hipMalloc((void**)d_array_addr_ptr,size*sizeof(double)),2601);
  }
#endif
}

static inline void gpuMemset_double(double* d_array, const size_t size, int value){
#ifdef USE_CUDA
  if (run_cuda){
    print_CUDA_error_if_any(cudaMemset(d_array,value,size*sizeof(double)),2602);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    print_HIP_error_if_any(hipMemset(d_array,value,size*sizeof(double)),2602);
  }
#endif
}

// overloads gpuFree(realw*)/gpuFree(int*) from mesh_constants_wrapper.h for double*
static inline void gpuFree(double* d_array){
#ifdef USE_CUDA
  if (run_cuda){ cudaFree(d_array); }
#endif
#ifdef USE_HIP
  if (run_hip){ hipFree(d_array); }
#endif
}

// dual CUDA/HIP kernel dispatch, to avoid repeating the same boilerplate for every launch below
#ifdef USE_CUDA
#define CG_LAUNCH_CUDA(kernel,grid,threads,stream,...) \
  if (run_cuda){ kernel<<<grid,threads,0,stream>>>(__VA_ARGS__); }
#else
#define CG_LAUNCH_CUDA(kernel,grid,threads,stream,...)
#endif
#ifdef USE_HIP
#define CG_LAUNCH_HIP(kernel,grid,threads,stream,...) \
  if (run_hip){ hipLaunchKernelGGL(kernel,dim3(grid),dim3(threads),0,stream,__VA_ARGS__); }
#else
#define CG_LAUNCH_HIP(kernel,grid,threads,stream,...)
#endif
#define CG_LAUNCH(kernel,grid,threads,stream,...) \
  do { CG_LAUNCH_CUDA(kernel,grid,threads,stream,__VA_ARGS__) CG_LAUNCH_HIP(kernel,grid,threads,stream,__VA_ARGS__) } while(0)

// 1D grid sizing for flat elementwise kernels operating on n scalars
static inline void cg_grid_1d(int n, dim3* grid, dim3* threads){
  int blocksize = BLOCKSIZE_TRANSFER;
  int size_padded = ((int)ceil(((double)n)/((double)blocksize)))*blocksize;
  int num_blocks_x,num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);
  *grid = dim3(num_blocks_x,num_blocks_y);
  *threads = dim3(blocksize,1,1);
}

// one block per boundary face, NGLL2 threads per block
static inline void cg_grid_faces(int nfaces, dim3* grid, dim3* threads){
  int num_blocks_x,num_blocks_y;
  get_blocks_xy(nfaces,&num_blocks_x,&num_blocks_y);
  *grid = dim3(num_blocks_x,num_blocks_y);
  *threads = dim3(NGLL2,1,1);
}


/* ----------------------------------------------------------------------------------------------- */

// elementwise vector kernels (operate on flat arrays of NDIM*NGLOB_AB doubles/realw)

/* ----------------------------------------------------------------------------------------------- */

__global__ void cg_set_const_kernel(double* vec, double val, int n){
  int i = threadIdx.x + (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x;
  if (i < n) vec[i] = val;
}

__global__ void cg_copy_kernel(double* dst, const double* __restrict__ src, int n){
  int i = threadIdx.x + (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x;
  if (i < n) dst[i] = src[i];
}

__global__ void cg_cast_d2f_kernel(const double* __restrict__ src, realw* dst, int n){
  int i = threadIdx.x + (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x;
  if (i < n) dst[i] = (realw) src[i];
}

__global__ void cg_cast_f2d_kernel(const realw* __restrict__ src, double* dst, int n){
  int i = threadIdx.x + (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x;
  if (i < n) dst[i] = (double) src[i];
}

// r = force_ext - kdotu (both realw), result stored as double
__global__ void cg_init_residual_kernel(const realw* __restrict__ force_ext, const realw* __restrict__ kdotu,
                                        double* r, int n){
  int i = threadIdx.x + (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x;
  if (i < n) r[i] = (double)force_ext[i] - (double)kdotu[i];
}

// y = y + alpha * x
__global__ void cg_axpy_kernel(double* y, double alpha, const double* __restrict__ x, int n){
  int i = threadIdx.x + (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x;
  if (i < n) y[i] += alpha * x[i];
}

// y = a .* b (elementwise)
__global__ void cg_mul_kernel(double* y, const double* __restrict__ a, const double* __restrict__ b, int n){
  int i = threadIdx.x + (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x;
  if (i < n) y[i] = a[i] * b[i];
}

// p = z + beta * p
__global__ void cg_xpay_kernel(double* p, const double* __restrict__ z, double beta, int n){
  int i = threadIdx.x + (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x;
  if (i < n) p[i] = z[i] + beta * p[i];
}


/* ----------------------------------------------------------------------------------------------- */

// reduction kernels for the parallel inner product (see parallel_inner_product() in static_module.f90)

/* ----------------------------------------------------------------------------------------------- */

__global__ void cg_max_abs_reduce_kernel(const double* __restrict__ vec, int n, double* blockmax){
  __shared__ double sdata[BLOCKSIZE_TRANSFER];
  unsigned int tid = threadIdx.x;
  unsigned int i = threadIdx.x + (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x;

  sdata[tid] = (i < n) ? fabs(vec[i]) : 0.0;
  __syncthreads();

  for (unsigned int s = blockDim.x/2; s > 0; s >>= 1){
    if (tid < s){
      if (sdata[tid] < sdata[tid+s]) sdata[tid] = sdata[tid+s];
    }
    __syncthreads();
  }
  if (tid == 0) blockmax[blockIdx.x+blockIdx.y*gridDim.x] = sdata[0];
}

// sum_node inv_mult(node) * dot(vec1(:,node)*inv_a, vec2(:,node)*inv_b)
__global__ void cg_weighted_dot_reduce_kernel(const double* __restrict__ vec1, const double* __restrict__ vec2,
                                              const realw* __restrict__ inv_mult, double inv_a, double inv_b,
                                              int nnode, double* blocksum){
  __shared__ double sdata[BLOCKSIZE_TRANSFER];
  unsigned int tid = threadIdx.x;
  unsigned int inode = threadIdx.x + (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x;

  double val = 0.0;
  if (inode < (unsigned int)nnode){
    double w = (double) inv_mult[inode];
    val = w * ( (vec1[inode*3]  *inv_a) * (vec2[inode*3]  *inv_b)
              + (vec1[inode*3+1]*inv_a) * (vec2[inode*3+1]*inv_b)
              + (vec1[inode*3+2]*inv_a) * (vec2[inode*3+2]*inv_b) );
  }
  sdata[tid] = val;
  __syncthreads();

  for (unsigned int s = blockDim.x/2; s > 0; s >>= 1){
    if (tid < s) sdata[tid] += sdata[tid+s];
    __syncthreads();
  }
  if (tid == 0) blocksum[blockIdx.x+blockIdx.y*gridDim.x] = sdata[0];
}


/* ----------------------------------------------------------------------------------------------- */

// boundary condition kernels (see enforce_fixed_bc() / enforce_roller_bc() in static_module.f90)

/* ----------------------------------------------------------------------------------------------- */

__global__ void cg_enforce_fixed_bc_kernel(double* vec, const int* __restrict__ d_ibool,
                                           const int* __restrict__ bdry_ispec, const int* __restrict__ bdry_ijk,
                                           int num_faces){
  int igll = threadIdx.x;
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  if (iface < num_faces){
    int ispec = bdry_ispec[iface]-1;

    int i = bdry_ijk[INDEX3(NDIM,NGLL2,0,igll,iface)]-1;
    int j = bdry_ijk[INDEX3(NDIM,NGLL2,1,igll,iface)]-1;
    int k = bdry_ijk[INDEX3(NDIM,NGLL2,2,igll,iface)]-1;
    int iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

    vec[iglob*3]   = 0.0;
    vec[iglob*3+1] = 0.0;
    vec[iglob*3+2] = 0.0;
  }
}

__global__ void cg_enforce_roller_bc_kernel(double* vec, const int* __restrict__ d_ibool,
                                            const int* __restrict__ bdry_ispec, const int* __restrict__ bdry_ijk,
                                            const realw* __restrict__ bdry_normal, int num_faces){
  int igll = threadIdx.x;
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  if (iface < num_faces){
    int ispec = bdry_ispec[iface]-1;

    int i = bdry_ijk[INDEX3(NDIM,NGLL2,0,igll,iface)]-1;
    int j = bdry_ijk[INDEX3(NDIM,NGLL2,1,igll,iface)]-1;
    int k = bdry_ijk[INDEX3(NDIM,NGLL2,2,igll,iface)]-1;
    int iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

    double nx = (double) bdry_normal[INDEX3(NDIM,NGLL2,0,igll,iface)];
    double ny = (double) bdry_normal[INDEX3(NDIM,NGLL2,1,igll,iface)];
    double nz = (double) bdry_normal[INDEX3(NDIM,NGLL2,2,igll,iface)];

    double vx = vec[iglob*3];
    double vy = vec[iglob*3+1];
    double vz = vec[iglob*3+2];
    double vn = vx*nx + vy*ny + vz*nz;

    // NOTE: nodes shared by several roller faces with different normals are raced on by
    // several blocks here. The CPU version avoids this by masking already-visited nodes
    // (first face wins); this GPU version does not reproduce that exact tie-break.
    vec[iglob*3]   = vx - vn*nx;
    vec[iglob*3+1] = vy - vn*ny;
    vec[iglob*3+2] = vz - vn*nz;
  }
}


/* ----------------------------------------------------------------------------------------------- */

// host-side helpers built on top of the kernels above

/* ----------------------------------------------------------------------------------------------- */

static void cg_enforce_fixed_bc(CGSolver* cg, Mesh* mp, double* vec){
  if (cg->num_fixed_bdry_faces == 0) return;
  dim3 grid,threads;
  cg_grid_faces(cg->num_fixed_bdry_faces,&grid,&threads);
  CG_LAUNCH(cg_enforce_fixed_bc_kernel,grid,threads,mp->compute_stream,
            vec,mp->d_ibool,cg->fixed_bdry_ispec,cg->fixed_bdry_ijk,cg->num_fixed_bdry_faces);
  GPU_ERROR_CHECKING("cg_enforce_fixed_bc_kernel");
}

static void cg_enforce_roller_bc(CGSolver* cg, Mesh* mp, double* vec){
  if (cg->num_roller_bdry_faces == 0) return;
  dim3 grid,threads;
  cg_grid_faces(cg->num_roller_bdry_faces,&grid,&threads);
  CG_LAUNCH(cg_enforce_roller_bc_kernel,grid,threads,mp->compute_stream,
            vec,mp->d_ibool,cg->roller_bdry_ispec,cg->roller_bdry_ijk,cg->roller_bdry_normal,cg->num_roller_bdry_faces);
  GPU_ERROR_CHECKING("cg_enforce_roller_bc_kernel");
}

// local (single block-reduction) max(|vec|) over the whole flat NDIM*NGLOB_AB array
static double cg_local_max_abs(const double* vec, int n, Mesh* mp,
                               dim3 grid, dim3 threads, int nblocks, double* d_partial, double* h_partial){
  CG_LAUNCH(cg_max_abs_reduce_kernel,grid,threads,mp->compute_stream,vec,n,d_partial);
  GPU_ERROR_CHECKING("cg_max_abs_reduce_kernel");
  gpuStreamSynchronize(mp->compute_stream);
  // note: d_partial holds `double`, mesh_constants_wrapper.h has no gpuMemcpy_tohost_double,
  // so we copy manually here (guarded the same way as the realw variants).
#ifdef USE_CUDA
  if (run_cuda) print_CUDA_error_if_any(cudaMemcpy(h_partial,d_partial,nblocks*sizeof(double),cudaMemcpyDeviceToHost),2603);
#endif
#ifdef USE_HIP
  if (run_hip) print_HIP_error_if_any(hipMemcpy(h_partial,d_partial,nblocks*sizeof(double),hipMemcpyDeviceToHost),2603);
#endif

  double m = h_partial[0];
  for (int b = 1; b < nblocks; b++) if (h_partial[b] > m) m = h_partial[b];
  return m;
}

// mirrors parallel_inner_product() in static_module.f90: scales by the max-abs of each vector
// before summing (avoids under/overflow), and reduces across MPI ranks.
static double cg_parallel_inner_product(CGSolver* cg, Mesh* mp, const double* vec1, const double* vec2){
  int n = NDIM * cg->NGLOB_AB;
  dim3 grid_flat,threads_flat;
  cg_grid_1d(n,&grid_flat,&threads_flat);
  int nblocks_flat = grid_flat.x*grid_flat.y;

  dim3 grid_node,threads_node;
  cg_grid_1d(cg->NGLOB_AB,&grid_node,&threads_node);
  int nblocks_node = grid_node.x*grid_node.y;

  int nblocks = nblocks_flat > nblocks_node ? nblocks_flat : nblocks_node;
  double* d_partial = NULL;
  double* h_partial = (double*) calloc(nblocks,sizeof(double));
  if (!h_partial) exit_on_error("Error allocating temporary host partial-reduction array");
  gpuMalloc_double((void**)&d_partial,nblocks);

  double a = cg_local_max_abs(vec1,n,mp,grid_flat,threads_flat,nblocks_flat,d_partial,h_partial);
  double b = cg_local_max_abs(vec2,n,mp,grid_flat,threads_flat,nblocks_flat,d_partial,h_partial);

#ifdef WITH_MPI
  double tmp;
  MPI_Allreduce(&a,&tmp,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD); a = tmp;
  MPI_Allreduce(&b,&tmp,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD); b = tmp;
#endif
  if (a < 1.0e-20) a = 1.0;
  if (b < 1.0e-20) b = 1.0;

  CG_LAUNCH(cg_weighted_dot_reduce_kernel,grid_node,threads_node,mp->compute_stream,
            vec1,vec2,cg->inv_mult,1.0/a,1.0/b,cg->NGLOB_AB,d_partial);
  GPU_ERROR_CHECKING("cg_weighted_dot_reduce_kernel");
  gpuStreamSynchronize(mp->compute_stream);
#ifdef USE_CUDA
  if (run_cuda) print_CUDA_error_if_any(cudaMemcpy(h_partial,d_partial,nblocks_node*sizeof(double),cudaMemcpyDeviceToHost),2604);
#endif
#ifdef USE_HIP
  if (run_hip) print_HIP_error_if_any(hipMemcpy(h_partial,d_partial,nblocks_node*sizeof(double),hipMemcpyDeviceToHost),2604);
#endif

  double local_sum = 0.0;
  for (int i = 0; i < nblocks_node; i++) local_sum += h_partial[i];

  double global_sum = local_sum;
#ifdef WITH_MPI
  MPI_Allreduce(&local_sum,&global_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif

  gpuFree(d_partial);
  free(h_partial);

  return global_sum * a * b;
}


/* ----------------------------------------------------------------------------------------------- */

// compute_forces_static_gpu(): element-wise stiffness application, mirrors
// compute_forces_static()/compute_forces_phase()/compute_elemwise_Kxu() in static_module.f90.
//
// The gradient/dot-product machinery below (sf_load_shared_memory_*, sf_sum_hprime*,
// sf_get_spatial_derivatives, sf_get_dot_product) is copied and trimmed from the standard
// elastic Kernel_2_noatt_iso_impl/Kernel_2_noatt_ani_impl kernels in
// gpu/kernels/Kernel_2_viscoelastic_impl.cu (same math, no attenuation/PML/gravity/textures,
// since compute_elemwise_Kxu() has none of those either). Only the nonlinear
// (is_nonlinear = .true., finite-strain / second Piola-Kirchhoff) branch of
// compute_elemwise_Kxu() is skipped, since static_problem_impl() never sets it.

/* ----------------------------------------------------------------------------------------------- */

__device__ __forceinline__ void sf_load_shared_memory_displ(int tx, int iglob, realw_const_p d_displ,
                                                             realw* sh_displx, realw* sh_disply, realw* sh_displz){
  sh_displx[tx] = d_displ[iglob*3];
  sh_disply[tx] = d_displ[iglob*3+1];
  sh_displz[tx] = d_displ[iglob*3+2];
}

__device__ __forceinline__ void sf_load_shared_memory_hprime(int tx, realw_const_p d_hprime_xx, realw* sh_hprime_xx){
  sh_hprime_xx[tx] = d_hprime_xx[tx];
}

__device__ __forceinline__ void sf_load_shared_memory_hprimewgll(int tx, realw_const_p d_hprimewgll_xx,
                                                                  realw* sh_hprimewgll_xx){
  sh_hprimewgll_xx[tx] = d_hprimewgll_xx[tx];
}

__device__ __forceinline__ void sf_sum_hprime_xi(int I,int J,int K,realw* tempxl,realw* tempyl,realw* tempzl,
                                                 realw* sh_tempx,realw* sh_tempy,realw* sh_tempz,realw* sh_hprime){
  realw sumx=0.f,sumy=0.f,sumz=0.f,fac;
  #pragma unroll
  for (int l=0;l<NGLLX;l++){
    fac = sh_hprime[l*NGLLX+I];
    sumx += sh_tempx[K*NGLL2+J*NGLLX+l] * fac;
    sumy += sh_tempy[K*NGLL2+J*NGLLX+l] * fac;
    sumz += sh_tempz[K*NGLL2+J*NGLLX+l] * fac;
  }
  *tempxl=sumx; *tempyl=sumy; *tempzl=sumz;
}

__device__ __forceinline__ void sf_sum_hprime_eta(int I,int J,int K,realw* tempxl,realw* tempyl,realw* tempzl,
                                                  realw* sh_tempx,realw* sh_tempy,realw* sh_tempz,realw* sh_hprime){
  realw sumx=0.f,sumy=0.f,sumz=0.f,fac;
  #pragma unroll
  for (int l=0;l<NGLLX;l++){
    fac = sh_hprime[l*NGLLX+J];
    sumx += sh_tempx[K*NGLL2+l*NGLLX+I] * fac;
    sumy += sh_tempy[K*NGLL2+l*NGLLX+I] * fac;
    sumz += sh_tempz[K*NGLL2+l*NGLLX+I] * fac;
  }
  *tempxl=sumx; *tempyl=sumy; *tempzl=sumz;
}

__device__ __forceinline__ void sf_sum_hprime_gamma(int I,int J,int K,realw* tempxl,realw* tempyl,realw* tempzl,
                                                    realw* sh_tempx,realw* sh_tempy,realw* sh_tempz,realw* sh_hprime){
  realw sumx=0.f,sumy=0.f,sumz=0.f,fac;
  #pragma unroll
  for (int l=0;l<NGLLX;l++){
    fac = sh_hprime[l*NGLLX+K];
    sumx += sh_tempx[l*NGLL2+J*NGLLX+I] * fac;
    sumy += sh_tempy[l*NGLL2+J*NGLLX+I] * fac;
    sumz += sh_tempz[l*NGLL2+J*NGLLX+I] * fac;
  }
  *tempxl=sumx; *tempyl=sumy; *tempzl=sumz;
}

__device__ __forceinline__ void sf_sum_hprimewgll_xi(int I,int J,int K,realw* tempxl,realw* tempyl,realw* tempzl,
                                                     realw* sh_tempx,realw* sh_tempy,realw* sh_tempz,realw* sh_hprimewgll){
  realw sumx=0.f,sumy=0.f,sumz=0.f,fac;
  #pragma unroll
  for (int l=0;l<NGLLX;l++){
    fac = sh_hprimewgll[I*NGLLX+l];
    sumx += sh_tempx[K*NGLL2+J*NGLLX+l] * fac;
    sumy += sh_tempy[K*NGLL2+J*NGLLX+l] * fac;
    sumz += sh_tempz[K*NGLL2+J*NGLLX+l] * fac;
  }
  *tempxl=sumx; *tempyl=sumy; *tempzl=sumz;
}

__device__ __forceinline__ void sf_sum_hprimewgll_eta(int I,int J,int K,realw* tempxl,realw* tempyl,realw* tempzl,
                                                      realw* sh_tempx,realw* sh_tempy,realw* sh_tempz,realw* sh_hprimewgll){
  realw sumx=0.f,sumy=0.f,sumz=0.f,fac;
  #pragma unroll
  for (int l=0;l<NGLLX;l++){
    fac = sh_hprimewgll[J*NGLLX+l];
    sumx += sh_tempx[K*NGLL2+l*NGLLX+I] * fac;
    sumy += sh_tempy[K*NGLL2+l*NGLLX+I] * fac;
    sumz += sh_tempz[K*NGLL2+l*NGLLX+I] * fac;
  }
  *tempxl=sumx; *tempyl=sumy; *tempzl=sumz;
}

__device__ __forceinline__ void sf_sum_hprimewgll_gamma(int I,int J,int K,realw* tempxl,realw* tempyl,realw* tempzl,
                                                        realw* sh_tempx,realw* sh_tempy,realw* sh_tempz,realw* sh_hprimewgll){
  realw sumx=0.f,sumy=0.f,sumz=0.f,fac;
  #pragma unroll
  for (int l=0;l<NGLLX;l++){
    fac = sh_hprimewgll[K*NGLLX+l];
    sumx += sh_tempx[l*NGLL2+J*NGLLX+I] * fac;
    sumy += sh_tempy[l*NGLL2+J*NGLLX+I] * fac;
    sumz += sh_tempz[l*NGLL2+J*NGLLX+I] * fac;
  }
  *tempxl=sumx; *tempyl=sumy; *tempzl=sumz;
}

__device__ __forceinline__ void sf_get_spatial_derivatives(realw* xixl,realw* xiyl,realw* xizl,
                                                            realw* etaxl,realw* etayl,realw* etazl,
                                                            realw* gammaxl,realw* gammayl,realw* gammazl,
                                                            realw* jacobianl,
                                                            int I,int J,int K,int tx,
                                                            realw* tempx1l,realw* tempy1l,realw* tempz1l,
                                                            realw* tempx2l,realw* tempy2l,realw* tempz2l,
                                                            realw* tempx3l,realw* tempy3l,realw* tempz3l,
                                                            realw* sh_tempx,realw* sh_tempy,realw* sh_tempz,
                                                            realw* sh_hprime_xx,
                                                            realw* duxdxl,realw* duxdyl,realw* duxdzl,
                                                            realw* duydxl,realw* duydyl,realw* duydzl,
                                                            realw* duzdxl,realw* duzdyl,realw* duzdzl,
                                                            realw_const_p d_xix,realw_const_p d_xiy,realw_const_p d_xiz,
                                                            realw_const_p d_etax,realw_const_p d_etay,realw_const_p d_etaz,
                                                            realw_const_p d_gammax,realw_const_p d_gammay,realw_const_p d_gammaz,
                                                            int ispec_irreg, realw xix_regular){
  if (ispec_irreg >= 0){
    int offset = ispec_irreg*NGLL3_PADDED + tx;
    *xixl = d_xix[offset]; *xiyl = d_xiy[offset]; *xizl = d_xiz[offset];
    *etaxl = d_etax[offset]; *etayl = d_etay[offset]; *etazl = d_etaz[offset];
    *gammaxl = d_gammax[offset]; *gammayl = d_gammay[offset]; *gammazl = d_gammaz[offset];

    *jacobianl = 1.f / ((*xixl)*((*etayl)*(*gammazl)-(*etazl)*(*gammayl))
                      -(*xiyl)*((*etaxl)*(*gammazl)-(*etazl)*(*gammaxl))
                      +(*xizl)*((*etaxl)*(*gammayl)-(*etayl)*(*gammaxl)));
  }

  sf_sum_hprime_xi(I,J,K,tempx1l,tempy1l,tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);
  sf_sum_hprime_eta(I,J,K,tempx2l,tempy2l,tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);
  sf_sum_hprime_gamma(I,J,K,tempx3l,tempy3l,tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx);

  __syncthreads();

  if (ispec_irreg >= 0){
    (*duxdxl) = (*xixl)*(*tempx1l) + (*etaxl)*(*tempx2l) + (*gammaxl)*(*tempx3l);
    (*duxdyl) = (*xiyl)*(*tempx1l) + (*etayl)*(*tempx2l) + (*gammayl)*(*tempx3l);
    (*duxdzl) = (*xizl)*(*tempx1l) + (*etazl)*(*tempx2l) + (*gammazl)*(*tempx3l);

    (*duydxl) = (*xixl)*(*tempy1l) + (*etaxl)*(*tempy2l) + (*gammaxl)*(*tempy3l);
    (*duydyl) = (*xiyl)*(*tempy1l) + (*etayl)*(*tempy2l) + (*gammayl)*(*tempy3l);
    (*duydzl) = (*xizl)*(*tempy1l) + (*etazl)*(*tempy2l) + (*gammazl)*(*tempy3l);

    (*duzdxl) = (*xixl)*(*tempz1l) + (*etaxl)*(*tempz2l) + (*gammaxl)*(*tempz3l);
    (*duzdyl) = (*xiyl)*(*tempz1l) + (*etayl)*(*tempz2l) + (*gammayl)*(*tempz3l);
    (*duzdzl) = (*xizl)*(*tempz1l) + (*etazl)*(*tempz2l) + (*gammazl)*(*tempz3l);
  } else {
    (*duxdxl) = xix_regular*(*tempx1l); (*duxdyl) = xix_regular*(*tempx2l); (*duxdzl) = xix_regular*(*tempx3l);
    (*duydxl) = xix_regular*(*tempy1l); (*duydyl) = xix_regular*(*tempy2l); (*duydzl) = xix_regular*(*tempy3l);
    (*duzdxl) = xix_regular*(*tempz1l); (*duzdyl) = xix_regular*(*tempz2l); (*duzdzl) = xix_regular*(*tempz3l);
  }
}

__device__ __forceinline__ void sf_get_dot_product(realw jacobianl,
                                                    realw sigma_xx,realw sigma_xy,realw sigma_yx,
                                                    realw sigma_xz,realw sigma_zx,realw sigma_yy,
                                                    realw sigma_yz,realw sigma_zy,realw sigma_zz,
                                                    realw Dxl,realw Dyl,realw Dzl,
                                                    realw* sh_tempx,realw* sh_tempy,realw* sh_tempz,
                                                    int tx, int ispec_irreg,realw xix_regular,realw jacobian_regular,
                                                    int component){
  if (threadIdx.x < NGLL3){
    if (ispec_irreg >= 0){
      sh_tempx[tx] = jacobianl * (sigma_xx*Dxl + sigma_yx*Dyl + sigma_zx*Dzl);
      sh_tempy[tx] = jacobianl * (sigma_xy*Dxl + sigma_yy*Dyl + sigma_zy*Dzl);
      sh_tempz[tx] = jacobianl * (sigma_xz*Dxl + sigma_yz*Dyl + sigma_zz*Dzl);
    } else if (component == 1){
      sh_tempx[tx] = jacobian_regular * (sigma_xx*xix_regular);
      sh_tempy[tx] = jacobian_regular * (sigma_xy*xix_regular);
      sh_tempz[tx] = jacobian_regular * (sigma_xz*xix_regular);
    } else if (component == 2){
      sh_tempx[tx] = jacobian_regular * (sigma_yx*xix_regular);
      sh_tempy[tx] = jacobian_regular * (sigma_yy*xix_regular);
      sh_tempz[tx] = jacobian_regular * (sigma_yz*xix_regular);
    } else {
      sh_tempx[tx] = jacobian_regular * (sigma_zx*xix_regular);
      sh_tempy[tx] = jacobian_regular * (sigma_zy*xix_regular);
      sh_tempz[tx] = jacobian_regular * (sigma_zz*xix_regular);
    }
  }
  __syncthreads();
}

// element-wise K*u application (isotropic or anisotropic) + rotation stiffness contribution
// + stress/strain output; mirrors compute_elemwise_Kxu() in static_module.f90 (linear-strain
// branch only). Scatters into d_kdotu via atomicAdd, exactly like compute_forces_phase() sums
// into kdotu(:,iglob) in Fortran.
__global__ void sf_compute_forces_static_kernel(
    const int nb_blocks_to_compute,
    const int* d_ibool,
    const int* d_phase_ispec_inner_elastic, const int num_phase_ispec_elastic, const int d_iphase,
    const int* d_irregular_element_number,
    realw_const_p d_displ,
    realw_p d_kdotu,
    realw_const_p d_xix,realw_const_p d_xiy,realw_const_p d_xiz,
    realw_const_p d_etax,realw_const_p d_etay,realw_const_p d_etaz,
    realw_const_p d_gammax,realw_const_p d_gammay,realw_const_p d_gammaz,
    const realw xix_regular, const realw jacobian_regular,
    realw_const_p d_hprime_xx, realw_const_p d_hprimewgll_xx,
    realw_const_p d_wgllwgll_xy, realw_const_p d_wgllwgll_xz, realw_const_p d_wgllwgll_yz,
    realw_const_p d_kappav, realw_const_p d_muv,
    const int ANISOTROPY,
    realw_const_p d_c11store,realw_const_p d_c12store,realw_const_p d_c13store,
    realw_const_p d_c14store,realw_const_p d_c15store,realw_const_p d_c16store,
    realw_const_p d_c22store,realw_const_p d_c23store,realw_const_p d_c24store,
    realw_const_p d_c25store,realw_const_p d_c26store,realw_const_p d_c33store,
    realw_const_p d_c34store,realw_const_p d_c35store,realw_const_p d_c36store,
    realw_const_p d_c44store,realw_const_p d_c45store,realw_const_p d_c46store,
    realw_const_p d_c55store,realw_const_p d_c56store,realw_const_p d_c66store,
    const int ROTATION,
    realw_const_p d_rhostore, realw_const_p d_wgll_cube,
    const realw omega_x, const realw omega_y, const realw omega_z,
    realw_p d_stress, realw_p d_strain, const int NSPEC_AB){

  int bx = blockIdx.y*gridDim.x + blockIdx.x;
  if (bx >= nb_blocks_to_compute) return;

  int tx = threadIdx.x;
  if (tx >= NGLL3) tx = NGLL3 - 1;

  int K = tx/NGLL2;
  int J = (tx-K*NGLL2)/NGLLX;
  int I = tx-K*NGLL2-J*NGLLX;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;
  realw fac1,fac2,fac3,lambdal,mul,lambdalplus2mul,kappal;
  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy;
  realw sum_terms1,sum_terms2,sum_terms3;
  realw c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66;
  realw ux0 = 0.f, uy0 = 0.f, uz0 = 0.f;

  __shared__ realw sh_tempx[NGLL3];
  __shared__ realw sh_tempy[NGLL3];
  __shared__ realw sh_tempz[NGLL3];
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

  int working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)] - 1;
  int ispec_irreg = d_irregular_element_number[working_element] - 1;
  int offset = working_element*NGLL3_PADDED + tx;
  int iglob = d_ibool[offset] - 1;

  if (threadIdx.x < NGLL3){
    sf_load_shared_memory_displ(tx,iglob,d_displ,sh_tempx,sh_tempy,sh_tempz);
    ux0 = sh_tempx[tx]; uy0 = sh_tempy[tx]; uz0 = sh_tempz[tx];
  }
  if (tx < NGLL2){
    sf_load_shared_memory_hprime(tx,d_hprime_xx,sh_hprime_xx);
    sf_load_shared_memory_hprimewgll(tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }

  __syncthreads();

  sf_get_spatial_derivatives(&xixl,&xiyl,&xizl,&etaxl,&etayl,&etazl,&gammaxl,&gammayl,&gammazl,&jacobianl,
                             I,J,K,tx,
                             &tempx1l,&tempy1l,&tempz1l,&tempx2l,&tempy2l,&tempz2l,&tempx3l,&tempy3l,&tempz3l,
                             sh_tempx,sh_tempy,sh_tempz,sh_hprime_xx,
                             &duxdxl,&duxdyl,&duxdzl,&duydxl,&duydyl,&duydzl,&duzdxl,&duzdyl,&duzdzl,
                             d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,
                             ispec_irreg,xix_regular);

  // jacobianl is only set above for irregular elements (see sf_get_spatial_derivatives)
  realw jac_l = (ispec_irreg >= 0) ? jacobianl : jacobian_regular;

  duxdxl_plus_duydyl = duxdxl+duydyl;
  duxdxl_plus_duzdzl = duxdxl+duzdzl;
  duydyl_plus_duzdzl = duydyl+duzdzl;
  duxdyl_plus_duydxl = duxdyl+duydxl;
  duzdxl_plus_duxdzl = duzdxl+duxdzl;
  duzdyl_plus_duydzl = duzdyl+duydzl;

  if (ANISOTROPY){
    c11=d_c11store[offset]; c12=d_c12store[offset]; c13=d_c13store[offset];
    c14=d_c14store[offset]; c15=d_c15store[offset]; c16=d_c16store[offset];
    c22=d_c22store[offset]; c23=d_c23store[offset]; c24=d_c24store[offset];
    c25=d_c25store[offset]; c26=d_c26store[offset]; c33=d_c33store[offset];
    c34=d_c34store[offset]; c35=d_c35store[offset]; c36=d_c36store[offset];
    c44=d_c44store[offset]; c45=d_c45store[offset]; c46=d_c46store[offset];
    c55=d_c55store[offset]; c56=d_c56store[offset]; c66=d_c66store[offset];

    sigma_xx = c11*duxdxl + c16*duxdyl_plus_duydxl + c12*duydyl + c15*duzdxl_plus_duxdzl + c14*duzdyl_plus_duydzl + c13*duzdzl;
    sigma_yy = c12*duxdxl + c26*duxdyl_plus_duydxl + c22*duydyl + c25*duzdxl_plus_duxdzl + c24*duzdyl_plus_duydzl + c23*duzdzl;
    sigma_zz = c13*duxdxl + c36*duxdyl_plus_duydxl + c23*duydyl + c35*duzdxl_plus_duxdzl + c34*duzdyl_plus_duydzl + c33*duzdzl;
    sigma_xy = c16*duxdxl + c66*duxdyl_plus_duydxl + c26*duydyl + c56*duzdxl_plus_duxdzl + c46*duzdyl_plus_duydzl + c36*duzdzl;
    sigma_xz = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl + c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl;
    sigma_yz = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl + c45*duzdxl_plus_duxdzl + c44*duzdyl_plus_duydzl + c34*duzdzl;
  } else {
    kappal = d_kappav[offset];
    mul = d_muv[offset];
    lambdalplus2mul = kappal + 1.33333333333333333333f * mul;
    lambdal = lambdalplus2mul - 2.0f * mul;

    sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
    sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
    sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;
    sigma_xy = mul*duxdyl_plus_duydxl;
    sigma_xz = mul*duzdxl_plus_duxdzl;
    sigma_yz = mul*duzdyl_plus_duydzl;
  }

  // stress/strain output, voigt notation (xx,yy,zz,yz,xz,xy), shape(NGLL3,NSPEC_AB,6)
  if (threadIdx.x < NGLL3){
    int idx = tx + NGLL3*working_element;
    int comp_stride = NGLL3*NSPEC_AB;

    d_strain[idx + comp_stride*0] = duxdxl;
    d_strain[idx + comp_stride*1] = duydyl;
    d_strain[idx + comp_stride*2] = duzdzl;
    d_strain[idx + comp_stride*3] = 0.5f*duzdyl_plus_duydzl;
    d_strain[idx + comp_stride*4] = 0.5f*duzdxl_plus_duxdzl;
    d_strain[idx + comp_stride*5] = 0.5f*duxdyl_plus_duydxl;

    d_stress[idx + comp_stride*0] = sigma_xx;
    d_stress[idx + comp_stride*1] = sigma_yy;
    d_stress[idx + comp_stride*2] = sigma_zz;
    d_stress[idx + comp_stride*3] = sigma_yz;
    d_stress[idx + comp_stride*4] = sigma_xz;
    d_stress[idx + comp_stride*5] = sigma_xy;
  }

  sigma_yx = sigma_xy;
  sigma_zx = sigma_xz;
  sigma_zy = sigma_yz;

  __syncthreads();
  sf_get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_yx,sigma_xz,sigma_zx,sigma_yy,sigma_yz,sigma_zy,sigma_zz,
                     xixl,xiyl,xizl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,1);
  sf_sum_hprimewgll_xi(I,J,K,&tempx1l,&tempy1l,&tempz1l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  __syncthreads();
  sf_get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_yx,sigma_xz,sigma_zx,sigma_yy,sigma_yz,sigma_zy,sigma_zz,
                     etaxl,etayl,etazl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,2);
  sf_sum_hprimewgll_eta(I,J,K,&tempx2l,&tempy2l,&tempz2l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  __syncthreads();
  sf_get_dot_product(jacobianl,sigma_xx,sigma_xy,sigma_yx,sigma_xz,sigma_zx,sigma_yy,sigma_yz,sigma_zy,sigma_zz,
                     gammaxl,gammayl,gammazl,sh_tempx,sh_tempy,sh_tempz,tx,ispec_irreg,xix_regular,jacobian_regular,3);
  sf_sum_hprimewgll_gamma(I,J,K,&tempx3l,&tempy3l,&tempz3l,sh_tempx,sh_tempy,sh_tempz,sh_hprimewgll_xx);

  fac1 = d_wgllwgll_yz[K*NGLLX+J];
  fac2 = d_wgllwgll_xz[K*NGLLX+I];
  fac3 = d_wgllwgll_xy[J*NGLLX+I];

  // note: no negation here (unlike the dynamic-solver accel kernels) -- this must match
  // force_x(i,j,k) = c11 = ... + fac1*newtempx1 + ... in compute_elemwise_Kxu(), which
  // accumulates +K*u directly (kdotu is used as "r = force_ext - kdotu" afterwards).
  sum_terms1 = fac1*tempx1l + fac2*tempx2l + fac3*tempx3l;
  sum_terms2 = fac1*tempy1l + fac2*tempy2l + fac3*tempy3l;
  sum_terms3 = fac1*tempz1l + fac2*tempz2l + fac3*tempz3l;

  // rotation stiffness contribution: + rho*jacobian*wgll_cube * (omega x (omega x u))
  // mirrors the two cross_product() calls + jacobianl scaling in compute_elemwise_Kxu()
  if (ROTATION){
    realw rho_jw = d_rhostore[offset] * jac_l * d_wgll_cube[tx];

    realw tx1 = omega_y*uz0 - omega_z*uy0;
    realw ty1 = omega_z*ux0 - omega_x*uz0;
    realw tz1 = omega_x*uy0 - omega_y*ux0;

    realw rx = omega_y*tz1 - omega_z*ty1;
    realw ry = omega_z*tx1 - omega_x*tz1;
    realw rz = omega_x*ty1 - omega_y*tx1;

    sum_terms1 += rho_jw * rx;
    sum_terms2 += rho_jw * ry;
    sum_terms3 += rho_jw * rz;
  }

  if (threadIdx.x < NGLL3){
    atomicAdd(&d_kdotu[iglob*3],   sum_terms1);
    atomicAdd(&d_kdotu[iglob*3+1], sum_terms2);
    atomicAdd(&d_kdotu[iglob*3+2], sum_terms3);
  }
}

static void sf_launch_phase(Mesh* mp, CGSolver* cg, int iphase, int num_elements,
                            realw* displ_in, realw* kdotu_out,
                            realw omega_x, realw omega_y, realw omega_z){
  if (num_elements == 0) return;

  int num_blocks_x,num_blocks_y;
  get_blocks_xy(num_elements,&num_blocks_x,&num_blocks_y);
  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLL3_PADDED,1,1);

  CG_LAUNCH(sf_compute_forces_static_kernel,grid,threads,mp->compute_stream,
            num_elements,
            mp->d_ibool,
            mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,iphase,
            mp->d_irregular_element_number,
            displ_in,
            kdotu_out,
            mp->d_xix,mp->d_xiy,mp->d_xiz,mp->d_etax,mp->d_etay,mp->d_etaz,mp->d_gammax,mp->d_gammay,mp->d_gammaz,
            mp->xix_regular,mp->jacobian_regular,
            mp->d_hprime_xx,mp->d_hprimewgll_xx,
            mp->d_wgllwgll_xy,mp->d_wgllwgll_xz,mp->d_wgllwgll_yz,
            mp->d_kappav,mp->d_muv,
            mp->ANISOTROPY,
            mp->d_c11store,mp->d_c12store,mp->d_c13store,mp->d_c14store,mp->d_c15store,mp->d_c16store,
            mp->d_c22store,mp->d_c23store,mp->d_c24store,mp->d_c25store,mp->d_c26store,mp->d_c33store,
            mp->d_c34store,mp->d_c35store,mp->d_c36store,mp->d_c44store,mp->d_c45store,mp->d_c46store,
            mp->d_c55store,mp->d_c56store,mp->d_c66store,
            mp->ROTATION,
            mp->d_rhostore,mp->d_wgll_cube,
            omega_x,omega_y,omega_z,
            cg->stress,cg->strain,cg->NSPEC_AB);
  GPU_ERROR_CHECKING("sf_compute_forces_static_kernel");
}

// gather/scatter kernels for kdotu's boundary-only MPI exchange (see compute_forces_static_gpu()
// below) -- copies of prepare_boundary_accel_on_device()/assemble_boundary_accel_on_device()
// (gpu/kernels/) operating on kdotu instead of d_accel; see the note above on why they can't
// just be called directly.

__global__ void cg_prepare_boundary_kdotu_kernel(realw* d_kdotu, realw* d_send_buffer,
                                                  const int num_interfaces_ext_mesh,
                                                  const int max_nibool_interfaces_ext_mesh,
                                                  const int* d_nibool_interfaces_ext_mesh,
                                                  const int* d_ibool_interfaces_ext_mesh){
  int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;
  for (int iinterface = 0; iinterface < num_interfaces_ext_mesh; iinterface++){
    if (id < d_nibool_interfaces_ext_mesh[iinterface]){
      int ientry = id + max_nibool_interfaces_ext_mesh*iinterface;
      int iglob = d_ibool_interfaces_ext_mesh[ientry] - 1;
      d_send_buffer[3*ientry]   = d_kdotu[3*iglob];
      d_send_buffer[3*ientry+1] = d_kdotu[3*iglob+1];
      d_send_buffer[3*ientry+2] = d_kdotu[3*iglob+2];
    }
  }
}

__global__ void cg_assemble_boundary_kdotu_kernel(realw* d_kdotu, realw* d_recv_buffer,
                                                   const int num_interfaces_ext_mesh,
                                                   const int max_nibool_interfaces_ext_mesh,
                                                   const int* d_nibool_interfaces_ext_mesh,
                                                   const int* d_ibool_interfaces_ext_mesh){
  int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;
  for (int iinterface = 0; iinterface < num_interfaces_ext_mesh; iinterface++){
    if (id < d_nibool_interfaces_ext_mesh[iinterface]){
      int ientry = id + max_nibool_interfaces_ext_mesh*iinterface;
      int iglob = d_ibool_interfaces_ext_mesh[ientry] - 1;
      atomicAdd(&d_kdotu[3*iglob],   d_recv_buffer[3*ientry]);
      atomicAdd(&d_kdotu[3*iglob+1], d_recv_buffer[3*ientry+1]);
      atomicAdd(&d_kdotu[3*iglob+2], d_recv_buffer[3*ientry+2]);
    }
  }
}

static void cg_boundary_grid(Mesh* mp, dim3* grid, dim3* threads){
  int blocksize = BLOCKSIZE_TRANSFER;
  int size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_ext_mesh)/((double)blocksize)))*blocksize;
  int num_blocks_x,num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);
  *grid = dim3(num_blocks_x,num_blocks_y);
  *threads = dim3(blocksize,1,1);
}

static void cg_prepare_boundary_kdotu(Mesh* mp, CGSolver* cg, realw* kdotu){
  dim3 grid,threads;
  cg_boundary_grid(mp,&grid,&threads);
  CG_LAUNCH(cg_prepare_boundary_kdotu_kernel,grid,threads,mp->compute_stream,
            kdotu,cg->d_send_buffer,mp->num_interfaces_ext_mesh,mp->max_nibool_interfaces_ext_mesh,
            mp->d_nibool_interfaces_ext_mesh,mp->d_ibool_interfaces_ext_mesh);
  GPU_ERROR_CHECKING("cg_prepare_boundary_kdotu_kernel");
}

static void cg_assemble_boundary_kdotu(Mesh* mp, CGSolver* cg, realw* kdotu){
  dim3 grid,threads;
  cg_boundary_grid(mp,&grid,&threads);
  CG_LAUNCH(cg_assemble_boundary_kdotu_kernel,grid,threads,mp->compute_stream,
            kdotu,cg->d_recv_buffer,mp->num_interfaces_ext_mesh,mp->max_nibool_interfaces_ext_mesh,
            mp->d_nibool_interfaces_ext_mesh,mp->d_ibool_interfaces_ext_mesh);
  GPU_ERROR_CHECKING("cg_assemble_boundary_kdotu_kernel");
}

// mirrors compute_forces_static() in static_module.f90: zeroes kdotu, then accumulates the
// outer- and inner-phase element contributions (elements chosen so that the outer phase covers
// every element touching an MPI interface, letting the interface contributions be sent while
// the inner phase is still computing -- see phase_ispec_inner_elastic in specfem_par).
static void compute_forces_static_gpu(Mesh* mp, CGSolver* cg, realw* displ_in, realw* kdotu_out,
                                      realw omega_x, realw omega_y, realw omega_z){
  int n = NDIM*cg->NGLOB_AB;
  gpuMemset_realw(kdotu_out,(size_t)n,0);

  sf_launch_phase(mp,cg,1,cg->nspec_outer_elastic,displ_in,kdotu_out,omega_x,omega_y,omega_z);

  // outer-phase elements are exactly the ones that touch an MPI interface, so kdotu_out is
  // already the complete local contribution at every interface node at this point (inner-phase
  // elements never touch interface nodes) -- safe to send now, same ordering as the CPU path.
  // only the interface-boundary values (mp->size_mpi_buffer reals) are copied to host, not the
  // whole NGLOB_AB array.
  if (cg->NPROC > 1 && mp->size_mpi_buffer > 0){
    cg_prepare_boundary_kdotu(mp,cg,kdotu_out);
    gpuStreamSynchronize(mp->compute_stream);
    gpuMemcpy_tohost_realw(cg->buffer_send_vector_ext_mesh,cg->d_send_buffer,mp->size_mpi_buffer);
    FC_FUNC_(assemble_mpi_vector_send_cuda,ASSEMBLE_MPI_VECTOR_SEND_CUDA)(
        &cg->NPROC,
        cg->buffer_send_vector_ext_mesh,cg->buffer_recv_vector_ext_mesh,
        &mp->num_interfaces_ext_mesh,&mp->max_nibool_interfaces_ext_mesh,
        cg->nibool_interfaces_ext_mesh,cg->my_neighbors_ext_mesh,
        cg->request_send_vector_ext_mesh,cg->request_recv_vector_ext_mesh);
  }

  sf_launch_phase(mp,cg,2,cg->nspec_inner_elastic,displ_in,kdotu_out,omega_x,omega_y,omega_z);

  if (cg->NPROC > 1 && mp->size_mpi_buffer > 0){
    for (int i = 0; i < mp->num_interfaces_ext_mesh; i++){
      FC_FUNC_(wait_req,WAIT_REQ)(&cg->request_recv_vector_ext_mesh[i]);
    }
    gpuMemcpy_todevice_realw(cg->d_recv_buffer,cg->buffer_recv_vector_ext_mesh,mp->size_mpi_buffer);
    cg_assemble_boundary_kdotu(mp,cg,kdotu_out);
    for (int i = 0; i < mp->num_interfaces_ext_mesh; i++){
      FC_FUNC_(wait_req,WAIT_REQ)(&cg->request_send_vector_ext_mesh[i]);
    }
  }
}


/* ----------------------------------------------------------------------------------------------- */

// allocates the CG solver structure and its device arrays

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(prepare_cg_solver_gpu,
              PREPARE_CG_SOLVER_GPU)(long* Container, long* Mesh_pointer, int* NGLOB_AB_f, int* NSPEC_AB_f, int* myrank_f,
                                     int* nspec_outer_elastic_f, int* nspec_inner_elastic_f,
                                     realw* force_ext, realw* inv_mult,
                                     int* num_fixed_bdry_faces_f, int* fixed_bdry_ispec, int* fixed_bdry_ijk,
                                     int* num_roller_bdry_faces_f, int* roller_bdry_ispec, int* roller_bdry_ijk,
                                     realw* roller_bdry_normal,
                                     int* NPROC_f,
                                     int* nibool_interfaces_ext_mesh, int* my_neighbors_ext_mesh,
                                     realw* buffer_send_vector_ext_mesh, realw* buffer_recv_vector_ext_mesh,
                                     int* request_send_vector_ext_mesh, int* request_recv_vector_ext_mesh){

  TRACE("prepare_cg_solver_gpu");

  Mesh* mp = (Mesh*)(*Mesh_pointer);

  // allocates structure
  CGSolver* cg = (CGSolver*) malloc(sizeof(CGSolver));
  if (! cg) exit_on_error("Error allocating CG solver pointer");

  // sets fortran pointer
  *Container = (long) cg;

  cg->NGLOB_AB = *NGLOB_AB_f;
  cg->NSPEC_AB = *NSPEC_AB_f;
  cg->myrank = *myrank_f;
  cg->nspec_outer_elastic = *nspec_outer_elastic_f;
  cg->nspec_inner_elastic = *nspec_inner_elastic_f;
  cg->num_fixed_bdry_faces = *num_fixed_bdry_faces_f;
  cg->num_roller_bdry_faces = *num_roller_bdry_faces_f;

  // MPI vector assembly: these Fortran arrays are allocated once for the whole run and stay
  // put, so it's safe to just hold onto the raw pointers instead of copying them. The interface
  // counts and device-side ibool/nibool arrays are already on Mesh* (shared with d_accel's own
  // assembly), so only the host-only pieces need to be stashed here.
  cg->NPROC = *NPROC_f;
  cg->nibool_interfaces_ext_mesh = nibool_interfaces_ext_mesh;
  cg->my_neighbors_ext_mesh = my_neighbors_ext_mesh;
  cg->buffer_send_vector_ext_mesh = buffer_send_vector_ext_mesh;
  cg->buffer_recv_vector_ext_mesh = buffer_recv_vector_ext_mesh;
  cg->request_send_vector_ext_mesh = request_send_vector_ext_mesh;
  cg->request_recv_vector_ext_mesh = request_recv_vector_ext_mesh;

  if (mp->size_mpi_buffer > 0){
    gpuMalloc_realw((void**)&cg->d_send_buffer,mp->size_mpi_buffer);
    gpuMalloc_realw((void**)&cg->d_recv_buffer,mp->size_mpi_buffer);
  } else {
    cg->d_send_buffer = NULL;
    cg->d_recv_buffer = NULL;
  }

  const int size = NDIM * cg->NGLOB_AB;

  // double precision PCG vectors (pure device scratch, no host data to copy)
  gpuMalloc_double((void**)&cg->r,size);
  gpuMalloc_double((void**)&cg->p,size);
  gpuMalloc_double((void**)&cg->Ap,size);
  gpuMalloc_double((void**)&cg->z,size);
  gpuMalloc_double((void**)&cg->inv_pred,size);
  gpuMalloc_double((void**)&cg->u,size);

  // preconditioner defaults to identity (no preconditioning implemented yet, matches Fortran)
  dim3 grid,threads;
  cg_grid_1d(size,&grid,&threads);
  CG_LAUNCH(cg_set_const_kernel,grid,threads,0,cg->inv_pred,1.0,size);
  GPU_ERROR_CHECKING("cg_set_const_kernel");

  // single precision buffers for the element-wise matrix-vector product
  gpuMalloc_realw((void**)&cg->kdotu,size);
  gpuMalloc_realw((void**)&cg->p_cr,size);
  gpuMalloc_realw((void**)&cg->Ap_cr,size);
  gpuMalloc_realw((void**)&cg->z_cr,size);

  // external force vector & inner-product weight
  gpuCreateCopy_todevice_realw((void**)&cg->force_ext,force_ext,size);
  gpuCreateCopy_todevice_realw((void**)&cg->inv_mult,inv_mult,cg->NGLOB_AB);

  // boundary condition data
  if (cg->num_fixed_bdry_faces > 0){
    gpuCreateCopy_todevice_int((void**)&cg->fixed_bdry_ispec,fixed_bdry_ispec,cg->num_fixed_bdry_faces);
    gpuCreateCopy_todevice_int((void**)&cg->fixed_bdry_ijk,fixed_bdry_ijk,NDIM*NGLL2*cg->num_fixed_bdry_faces);
  } else {
    cg->fixed_bdry_ispec = NULL;
    cg->fixed_bdry_ijk = NULL;
  }

  if (cg->num_roller_bdry_faces > 0){
    gpuCreateCopy_todevice_int((void**)&cg->roller_bdry_ispec,roller_bdry_ispec,cg->num_roller_bdry_faces);
    gpuCreateCopy_todevice_int((void**)&cg->roller_bdry_ijk,roller_bdry_ijk,NDIM*NGLL2*cg->num_roller_bdry_faces);
    gpuCreateCopy_todevice_realw((void**)&cg->roller_bdry_normal,roller_bdry_normal,NDIM*NGLL2*cg->num_roller_bdry_faces);
  } else {
    cg->roller_bdry_ispec = NULL;
    cg->roller_bdry_ijk = NULL;
    cg->roller_bdry_normal = NULL;
  }

  // stress/strain output of compute_forces_static_gpu(), shape(NGLL3,NSPEC_AB,6)
  gpuMalloc_realw((void**)&cg->stress,NGLL3*cg->NSPEC_AB*6);
  gpuMalloc_realw((void**)&cg->strain,NGLL3*cg->NSPEC_AB*6);
}


/* ----------------------------------------------------------------------------------------------- */

// runs the PCG solver, mirrors static_problem_impl() in static_module.f90
// (minus the ssol%USE_PETSC_AS_BACKEND / PETSc branch, which stays on the CPU side)

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(run_cg_solver_gpu,
              RUN_CG_SOLVER_GPU)(long* Container, long* Mesh_pointer, int* maxiter_f, realw* omega,
                                 double* rsinit_out, double* rsfinal_out, int* iter_out,
                                 realw* displ_out, realw* stress_out, realw* strain_out){

  TRACE("run_cg_solver_gpu");

  CGSolver* cg = (CGSolver*)(*Container);
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  const int n = NDIM * cg->NGLOB_AB;
  const int maxiter = *maxiter_f;
  const double tol = 1.0e-6;
  const realw omega_x = omega[0], omega_y = omega[1], omega_z = omega[2];

  dim3 grid,threads;
  cg_grid_1d(n,&grid,&threads);

  double alpha,beta,rsnew,rsinit,temp_sum,rho_old,rho_new;
  int iter = 0;

  // -----------------------------------------------------
  // 1. Initialization: u = 0, displ = 0
  // -----------------------------------------------------
  gpuMemset_double(cg->u,size_t(n),0);
  gpuMemset_realw(mp->d_displ,size_t(n),0);

  compute_forces_static_gpu(mp,cg,mp->d_displ,cg->kdotu,omega_x,omega_y,omega_z);
  CG_LAUNCH(cg_init_residual_kernel,grid,threads,mp->compute_stream,cg->force_ext,cg->kdotu,cg->r,n);
  GPU_ERROR_CHECKING("cg_init_residual_kernel");

  cg_enforce_fixed_bc(cg,mp,cg->r);
  cg_enforce_roller_bc(cg,mp,cg->r);

  rsinit = cg_parallel_inner_product(cg,mp,cg->r,cg->r);
  if (cg->myrank == 0){
    printf("\n----------------------------------------\n");
    printf("Starting PCG Solver for Static Problem (GPU)\n");
    printf("Initial residual norm: %e\n",sqrt(rsinit));
    printf("----------------------------------------\n\n");
  }
  if (rsinit < 1.0e-20) rsinit = 1.0;

  // -----------------------------------------------------
  // 2. Apply preconditioner (identity for now): z = M^-1 * r, p = z
  // -----------------------------------------------------
  CG_LAUNCH(cg_mul_kernel,grid,threads,mp->compute_stream,cg->z,cg->r,cg->inv_pred,n);
  GPU_ERROR_CHECKING("cg_mul_kernel");
  CG_LAUNCH(cg_copy_kernel,grid,threads,mp->compute_stream,cg->p,cg->z,n);
  GPU_ERROR_CHECKING("cg_copy_kernel");

  rho_old = cg_parallel_inner_product(cg,mp,cg->r,cg->z);

  // -----------------------------------------------------
  // 3. PCG loop
  // -----------------------------------------------------
  rsnew = rsinit;
  for (iter = 1; iter <= maxiter; iter++){

    cg_enforce_fixed_bc(cg,mp,cg->p);
    cg_enforce_roller_bc(cg,mp,cg->p);

    // Ap = K * p
    CG_LAUNCH(cg_cast_d2f_kernel,grid,threads,mp->compute_stream,cg->p,cg->p_cr,n);
    GPU_ERROR_CHECKING("cg_cast_d2f_kernel");
    compute_forces_static_gpu(mp,cg,cg->p_cr,cg->Ap_cr,omega_x,omega_y,omega_z);
    CG_LAUNCH(cg_cast_f2d_kernel,grid,threads,mp->compute_stream,cg->Ap_cr,cg->Ap,n);
    GPU_ERROR_CHECKING("cg_cast_f2d_kernel");

    temp_sum = cg_parallel_inner_product(cg,mp,cg->p,cg->Ap);
    if (temp_sum < 1.0e-30){
      if (cg->myrank == 0) printf("Warning: Curvature too small/negative, stopping. %d %e\n",iter,temp_sum);
      break;
    }
    alpha = rho_old / temp_sum;

    // u = u + alpha * p ; displ = u
    CG_LAUNCH(cg_axpy_kernel,grid,threads,mp->compute_stream,cg->u,alpha,cg->p,n);
    GPU_ERROR_CHECKING("cg_axpy_kernel (u)");
    CG_LAUNCH(cg_cast_d2f_kernel,grid,threads,mp->compute_stream,cg->u,mp->d_displ,n);
    GPU_ERROR_CHECKING("cg_cast_d2f_kernel (displ)");

    // r = r - alpha * Ap
    CG_LAUNCH(cg_axpy_kernel,grid,threads,mp->compute_stream,cg->r,-alpha,cg->Ap,n);
    GPU_ERROR_CHECKING("cg_axpy_kernel (r)");
    cg_enforce_fixed_bc(cg,mp,cg->r);
    cg_enforce_roller_bc(cg,mp,cg->r);

    rsnew = cg_parallel_inner_product(cg,mp,cg->r,cg->r);

    if (iter % 1000 == 0 || iter == 5 || iter == maxiter){
      // note: max_all_all_dp() in the Fortran version is a collective call, so this runs
      // on every rank -- only the printf below is gated to rank 0.
      int nblocks = grid.x*grid.y;
      double* d_tmp = NULL;
      double* h_tmp = (double*) calloc(nblocks,sizeof(double));
      if (!h_tmp) exit_on_error("Error allocating temporary host max-reduction array");
      gpuMalloc_double((void**)&d_tmp,nblocks);

      double max_u = cg_local_max_abs(cg->u,n,mp,grid,threads,nblocks,d_tmp,h_tmp);
#ifdef WITH_MPI
      double tmp;
      MPI_Allreduce(&max_u,&tmp,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      max_u = tmp;
#endif
      gpuFree(d_tmp);
      free(h_tmp);

      if (cg->myrank == 0)
        printf("Iter: %d  Rel Resid: %e  norm of Displ: %e\n",iter,sqrt(rsnew/rsinit),max_u);
    }

    if (sqrt(rsnew/rsinit) < tol) break;

    // z = r * inv_pred ; beta = (r.z)_new / (r.z)_old ; p = z + beta*p
    CG_LAUNCH(cg_mul_kernel,grid,threads,mp->compute_stream,cg->z,cg->r,cg->inv_pred,n);
    GPU_ERROR_CHECKING("cg_mul_kernel");
    rho_new = cg_parallel_inner_product(cg,mp,cg->r,cg->z);
    beta = rho_new / rho_old;
    CG_LAUNCH(cg_xpay_kernel,grid,threads,mp->compute_stream,cg->p,cg->z,beta,n);
    GPU_ERROR_CHECKING("cg_xpay_kernel");

    rho_old = rho_new;
  }

  // apply boundary conditions to u, then recompute stress/strain on the converged solution
  cg_enforce_fixed_bc(cg,mp,cg->u);
  cg_enforce_roller_bc(cg,mp,cg->u);
  CG_LAUNCH(cg_cast_d2f_kernel,grid,threads,mp->compute_stream,cg->u,mp->d_displ,n);
  GPU_ERROR_CHECKING("cg_cast_d2f_kernel (final displ)");
  compute_forces_static_gpu(mp,cg,mp->d_displ,cg->kdotu,omega_x,omega_y,omega_z);

  // copies the converged solution back to the host arrays passed in from Fortran
  // (displ, ssol%stress, ssol%strain all live on the host -- the CG loop above only
  // ever touches mp->d_displ / cg->stress / cg->strain on the device)
  gpuStreamSynchronize(mp->compute_stream);
  gpuMemcpy_tohost_realw(displ_out,mp->d_displ,n);
  gpuMemcpy_tohost_realw(stress_out,cg->stress,NGLL3*cg->NSPEC_AB*6);
  gpuMemcpy_tohost_realw(strain_out,cg->strain,NGLL3*cg->NSPEC_AB*6);

  *rsinit_out = sqrt(rsinit);
  *rsfinal_out = sqrt(rsnew/rsinit);
  *iter_out = iter;
}
