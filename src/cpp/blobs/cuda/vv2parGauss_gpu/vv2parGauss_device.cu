/*_________________________________________________________________________
* vv2parCC_device_066DP.cu - calculates the self-induced velocity field of the wake.
* Parallel version on GPU - CUDA code executed on device
*
* CUDA kernel function (executed on the device, called from the host) +
* CUDA block & thread functions (executed on the device, called from device)
* Manages data flow, launches and syncronize threads blocks
*
* DUWIND- Delft University Wind Energy Research Institute
* developer: Giuseppe Tescione
*
* Version: 0.6.6DP (alpha) - 20110824
* basic version with no loop unrolling, no wrap and no multithread bodies
* simple cut-off constant for desingularization
* double precision (for GPUs of computing capability 2.x)
*________________________________________________________________________*/

//Definition of double2 and double3 types
//typedef struct {
//double x, y;
//} double2;

//typedef struct {
//double x, y, z;
//} double3;



__constant__ int    blocksize_gpu;
__constant__ int    nParticles_gpu;
__constant__ int    nTargets_gpu;
__constant__ int    nParticleBlocks_gpu;
__constant__ int    nTargetBlocks_gpu;
__constant__ double  ksigmasqr_gpu;
__constant__ double  inv_2pi_gpu;
__constant__ double myeps_gpu;
/* constants (block dimension, number of particle, cut-off and 1/2pi) residing in
constant memory space, accessible from all threads within the grid and from the host.
Defined in host code*/

__device__ double2 vv2par_thread(double2 THR_velocities, double THR_xTarget, double THR_yTarget, double THR_xBlob, double THR_yBlob, double THR_wBlob)
/*THREAD FUNCTION - set of instructions performed parallely by each processor.
Calculates velocity induction on target particles by source particle.
Takes as input:
THR_UIND (2x double THR_UIND.x & THR_UIND.y) -> velocity induction of target particles
already computed by previous thread blocks to which adds the new induction;
THR_TARG (3x double THR_TARG.x & THR_TARG.y & THR_TARG.z) -> position (x, y) and
vorticity (z) of target particles.Position is needed to calculate induction
but vorticity is not used but kept to mantain data structure coeherency;
THR_SRC (3x double THR_SRC.x & THR_SRC.y & THR_SRC.z) -> position (x, y) and
vorticity (z) of source particles, needed to calculate induction.
Gives as output:
THR_UIND (2x double THR_UIND.x & THR_UIND.y) -> updated velocity induction of targets */

{
 //targets-particle distance, local variable [2 FLOPS]
// 	printf("Thread %d; xTarget %f; yTarget %f; xBlob %f; yBlob %f, wBlob %f\n -----------------------------------------\n",threadIdx.x,THR_xTarget,THR_yTarget,THR_xBlob,THR_yBlob,THR_wBlob);

 double2 RAD;
 RAD.x = THR_xTarget - THR_xBlob;
 RAD.y = THR_yTarget - THR_yBlob;

 //square of distance plus cut-off, local variable [4 FLOPS]
 double RADSQR = RAD.x * RAD.x + RAD.y * RAD.y + myeps_gpu;

 //vorticity/(2pi*sqr(rad)) [2 FLOPS]
 double S = THR_wBlob * inv_2pi_gpu / RADSQR;

 //update velocity induction [4 FLOPS]
 THR_velocities.x -= S * RAD.y * (1.0-exp(-RADSQR/(ksigmasqr_gpu)));
 THR_velocities.y += S * RAD.x * (1.0-exp(-RADSQR/(ksigmasqr_gpu)));

 return THR_velocities;
}

__device__ double2 vv2par_block(double BLK_xTarget, double BLK_yTarget, double2 BLK_velocities)
/*BLOCK FUNCTION - data & execution management for thread block
Evaluate induction in a pxp block
Takes as input:
BLK_TARG (3x double BLK_TARG.x & BLK_TARG.y & BLK_TARG.z) -> position (x, y)
and vorticity (z) of target. Passed unchanged to THREAD CODE as TARGET;
BLK_UIND (2x double BLK_UIND.x & BLK_UIND.y) -> velocity induction of target
particles. Passed unchanged to THREAD CODE as UIND.
Gives as output:
BLK_UIND (2x double BLK_UIND.x & BLK_UIND.y) -> updated velocity induction
of target. Received unchanged by THREAD CODE as UIND */

{
 extern __shared__ double BLK_blob [];  
 //extern __shared__ double BLK_yBlob []; 
 //extern __shared__ double BLK_wBlob []; 
 /* External variable residing in shared memory space of thread block,
 accessible from all threads within the block. Source particles data array
 (position (x & y) and vorticity (z)) common to the block.
 Size of the array is determined at launch time with instruction []     */

 //call thread function for every thread in block
 for ( int i = 0; i <blockDim.x; i++)
 {
   BLK_velocities = vv2par_thread(BLK_velocities, BLK_xTarget, BLK_yTarget, BLK_blob[i], BLK_blob[i+blocksize_gpu], BLK_blob[i+2*blocksize_gpu]);
   
 }

 return (BLK_velocities);
}

__global__ void vv2par_kernel(void *cxBlob_gpu_ondevice, void *cyBlob_gpu_ondevice, void *cwBlob_gpu_ondevice, void *cxTarget_gpu_ondevice, void *cyTarget_gpu_ondevice, void *cvx_gpu_ondevice, void *cvy_gpu_ondevice)
/*KERNEL FUNCTION - data & execution management for block grid
Kernel executed on the device, called from the host.
Manages memory passages from host to device and executes block function
Takes as input:
*ONDEV_POS and *ONDEV_IND -> pointers to global device memory for the
position and induction of particles                                     */
{
 extern __shared__ double BLK_blob [];  //see above
 //extern __shared__ double BLK_yBlob [];  //see above
 //extern __shared__ double BLK_wBlob [];  //see above

 //pointers passage
 double * KRN_xBlob = (double *)cxBlob_gpu_ondevice;
 double * KRN_yBlob = (double *)cyBlob_gpu_ondevice;
 double * KRN_wBlob = (double *)cwBlob_gpu_ondevice;
 double * KRN_xTarget = (double *)cxTarget_gpu_ondevice;
 double * KRN_yTarget = (double *)cyTarget_gpu_ondevice;
 double * KRN_vx = (double *)cvx_gpu_ondevice;
 double * KRN_vy = (double *)cvy_gpu_ondevice;

 //induction initialization
 double2 BLK_velocities;
 BLK_velocities.x = 0;
 BLK_velocities.y = 0;

 //target particles definition
 double BLK_xTarget;
 double BLK_yTarget;
 int NTHR = blockIdx.x * blockDim.x + threadIdx.x;
 BLK_xTarget = KRN_xTarget[NTHR];
 BLK_yTarget = KRN_yTarget[NTHR];
 //printf("Block %d; Thread %d :: Before the loop\n",blockIdx.x,threadIdx.x);
 int i, block;
 for (i = 0, block = 0; i < nParticles_gpu; i += blocksize_gpu, block++)//LOOP over blocks
 {
   //source particle definition (shared data)
   int id = block * blockDim.x + threadIdx.x;
   BLK_blob [threadIdx.x] = KRN_xBlob[id];
   BLK_blob [threadIdx.x + blocksize_gpu] = KRN_yBlob[id];
   BLK_blob [threadIdx.x + 2*blocksize_gpu] = KRN_wBlob[id];
	
	
   __syncthreads();
   // all shared memory locations are populated before starting computation

   BLK_velocities = vv2par_block(BLK_xTarget, BLK_yTarget, BLK_velocities); //block function call

   __syncthreads();
   //all threads within block finish computation before advancing next block
 }
 //save results in global memory
 double2 UIND = {BLK_velocities.x, BLK_velocities.y};
 KRN_vx[NTHR] = UIND.x;
 KRN_vy[NTHR] = UIND.y;
}
