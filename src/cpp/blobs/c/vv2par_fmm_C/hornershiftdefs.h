//Define addition operator for upward and downward shift, which element to store and if SSE2 should be used
//Both upward and downward shift are built from the same structure, 
//but with additions going in different directions. The element stored after each operation
//is different as well, but the load is identical. 

//Note that PERMUTE2 corresponds to UNROLLLEVEL1 etc, but operates on two elements (with one addition)
#ifdef __CUDACC__
#define ADDEUP(i1,i2) a[i2]+=a[i1];
#define SUBEUP(i1,i2) a[i2]-=a[i1];
#define ADDEDN(i1,i2) a[i1]+=a[i2];
#define LOADM1(i1,in,base,offset) a[i1]=in[(base)*M2PSMAXTHREADS+threadIdx.x];
#define STOREMUP(in,out,offset,i2,i3) in[(out)*M2PSMAXTHREADS+threadIdx.x]=a[i3]
#define STOREMDN(in,out,offset,i2,i3) in[(out)*M2PSMAXTHREADS+threadIdx.x]=a[i2]
#define LOADM1M2M(i1,in,base,offset) a[i1]=((double*)&localcoeff[base])[threadIdx.x&1];
#define STOREMUPM2M(in,out,offset,i2,i3) ((double*)&localcoeff[out])[threadIdx.x&1]=a[i3]
#define STOREMDNM2M(in,out,offset,i2,i3) ((double*)&localcoeff[out])[threadIdx.x&1]=a[i2]
#define LOADM1M2MT(i1,in,base,offset) a[i1]=wksp[(base)*blockDim.x+threadIdx.x];
#define STOREMUPM2MT(in,out,offset,i2,i3) wksp[(out)*blockDim.x+threadIdx.x]=a[i3]
#define STOREMDNM2MT(in,out,offset,i2,i3) wksp[(out)*blockDim.x+threadIdx.x]=a[i2]

#else

#ifdef SSE2INTRINSIC
#ifdef FULLCROSSORDER //alternates wksp1 and wksp2
#define ADDEUP2(i1,i2) a[i2]=_mm_add_pd(a[i2],a[i1]);a2[i2]=_mm_add_pd(a2[i2],a2[i1])
#define ADDEDN2(i1,i2) a[i1]=_mm_add_pd(a[i1],a[i2]);a2[i1]=_mm_add_pd(a2[i1],a2[i2])
#define LOADM2(i1,in,base,offset) a[i1]=_mm_load_pd(&wksp1[(base)*2]);a2[i1]=_mm_load_pd(&wksp2[(base)*2])
#define STOREMUP2(in,out,offset,i2,i3) _mm_store_pd(&wksp1[(out)*2],a[i3]);_mm_store_pd(&wksp2[(out)*2],a2[i3])
#define STOREMDN2(in,out,offset,i2,i3) _mm_store_pd(&wksp1[(out)*2],a[i2]);_mm_store_pd(&wksp2[(out)*2],a2[i2])
#endif

#define ADDEUP(i1,i2) a[i2]=_mm_add_pd(a[i2],a[i1])
#define ADDEDN(i1,i2) a[i1]=_mm_add_pd(a[i1],a[i2])
#define LOADM1(i1,in,base,offset) a[i1]=_mm_load_pd(&in[(base)*2])
//#define STOREM(in,out,offset,i2) _mm_store_pd(&in[(out)*2],a[i2])
#define STOREMUP(in,out,offset,i2,i3) _mm_store_pd(&in[(out)*2],a[i3])
#define STOREMDN(in,out,offset,i2,i3) _mm_store_pd(&in[(out)*2],a[i2])
#define READMULTIPLIER 1

#else /*SSE2INTRINSIC*/

#define READMULTIPLIER 2
#ifdef FULLCROSSORDER //alternates wksp1 and wksp2 and real and imag part
#define ADDEUP2(i1,i2) a[i2]+=a[i1];a2[i2]+=a2[i1];a3[i2]+=a3[i1];a4[i2]+=a4[i1]
#define ADDEDN2(i1,i2) a[i1]+=a[i2];a2[i1]+=a2[i2];a3[i1]+=a3[i2];a4[i1]+=a4[i2]
#define LOADM2(i1,in,base,offset) a[i1]=wksp1[(base)*READMULTIPLIER];    a2[i1]=wksp1[(base)*READMULTIPLIER+1];a3[i1]=wksp2[(base)*READMULTIPLIER];a4[i1]=wksp2[(base)*READMULTIPLIER+1]
#define STOREMUP2(in,out,offset,i2,i3) wksp1[(out)*READMULTIPLIER]=a[i3];wksp1[(out)*READMULTIPLIER+1]=a2[i3]; wksp2[(out)*READMULTIPLIER]=a3[i3]; wksp2[(out)*READMULTIPLIER+1]=a4[i3]
#define STOREMDN2(in,out,offset,i2,i3) wksp1[(out)*READMULTIPLIER]=a[i2];wksp1[(out)*READMULTIPLIER+1]=a2[i2]; wksp2[(out)*READMULTIPLIER]=a3[i2]; wksp2[(out)*READMULTIPLIER+1]=a4[i2]
#endif

#ifdef CROSSORDER
#define ADDEUP(i1,i2) a[i2]+=a[i1];a2[i2]+=a2[i1]
#define ADDEDN(i1,i2) a[i1]+=a[i2];a2[i1]+=a2[i2]
#define LOADM1(i1,in,base,offset) a[i1]=in[(base)*READMULTIPLIER];    a2[i1]=in[(base)*READMULTIPLIER+1]
//#define STOREM(in,out,offset,i2) in[(out)*READMULTIPLIER]=a[i2];     in[(out)*READMULTIPLIER+1]=a2[i2]
#define STOREMUP(in,out,offset,i2,i3) in[(out)*READMULTIPLIER]=a[i3];in[(out)*READMULTIPLIER+1]=a2[i3]
#define STOREMDN(in,out,offset,i2,i3) in[(out)*READMULTIPLIER]=a[i2];in[(out)*READMULTIPLIER+1]=a2[i2]
#else
#define ADDEUP(i1,i2) a[i2]+=a[i1]
#define ADDEDN(i1,i2) a[i1]+=a[i2]
#define LOADM1(i1,in,base,offset) a[i1]=in[(base)*READMULTIPLIER+offset]
//#define STOREM(in,out,offset,i2) in[(out)*READMULTIPLIER+offset]=a[i2]
#define STOREMUP(in,out,offset,i2,i3) in[(out)*READMULTIPLIER+offset]=a[i3]
#define STOREMDN(in,out,offset,i2,i3) in[(out)*READMULTIPLIER+offset]=a[i2]
#endif
#endif /*SSE2INTRINSIC*/
#endif /*__CUDACC__*/
//The permutations of the local variables, a should be a local variable of size equal to unrolllevel+1
#define PERMUTE2(i1,i2,ADDE) ADDE(i1,i2)
#define PERMUTE3(i1,i2,i3,ADDE) PERMUTE2(i1,i2,ADDE); ADDE(i2,i3)
#define PERMUTE4(i1,i2,i3,i4,ADDE) PERMUTE3(i1,i2,i3,ADDE); ADDE(i3,i4)
#define PERMUTE5(i1,i2,i3,i4,i5,ADDE) PERMUTE4(i1,i2,i3,i4,ADDE); ADDE(i4,i5)
#define PERMUTE6(i1,i2,i3,i4,i5,i6,ADDE) PERMUTE5(i1,i2,i3,i4,i5,ADDE); ADDE(i5,i6)
#define PERMUTE7(i1,i2,i3,i4,i5,i6,i7,ADDE) PERMUTE6(i1,i2,i3,i4,i5,i6,ADDE); ADDE(i6,i7)
#define PERMUTE8(i1,i2,i3,i4,i5,i6,i7,i8,ADDE) PERMUTE7(i1,i2,i3,i4,i5,i6,i7,ADDE); ADDE(i7,i8)

//Permutation and store last element to memory
#define PERMUTE2O(in,out,offset,STOREM,ADDE,i1,i2) PERMUTE2(i1,i2,ADDE); STOREM(in,out,offset,i1,i2)
#define PERMUTE3O(in,out,offset,STOREM,ADDE,i1,i2,i3) PERMUTE3(i1,i2,i3,ADDE); STOREM(in,out,offset,i2,i3)
#define PERMUTE4O(in,out,offset,STOREM,ADDE,i1,i2,i3,i4) PERMUTE4(i1,i2,i3,i4,ADDE); STOREM(in,out,offset,i3,i4)
#define PERMUTE5O(in,out,offset,STOREM,ADDE,i1,i2,i3,i4,i5) PERMUTE5(i1,i2,i3,i4,i5,ADDE); STOREM(in,out,offset,i4,i5)
#define PERMUTE6O(in,out,offset,STOREM,ADDE,i1,i2,i3,i4,i5,i6) PERMUTE6(i1,i2,i3,i4,i5,i6,ADDE); STOREM(in,out,offset,i5,i6)
#define PERMUTE7O(in,out,offset,STOREM,ADDE,i1,i2,i3,i4,i5,i6,i7) PERMUTE7(i1,i2,i3,i4,i5,i6,i7,ADDE); STOREM(in,out,offset,i6,i7)
#define PERMUTE8O(in,out,offset,STOREM,ADDE,i1,i2,i3,i4,i5,i6,i7,i8) PERMUTE8(i1,i2,i3,i4,i5,i6,i7,i8,ADDE); STOREM(in,out,offset,i7,i8)

//Permutation ans load one element as well. Used all times except first time, when all elements are loaded
#define PERMUTE2IO(in,base,offset,out,STOREM,LOADM,ADDE,i1,i2) LOADM(i1,in,base,offset);PERMUTE2O(in,out,offset,STOREM,ADDE,i1,i2)
#define PERMUTE3IO(in,base,offset,out,STOREM,LOADM,ADDE,i1,i2,i3) LOADM(i1,in,base,offset);PERMUTE3O(in,out,offset,STOREM,ADDE,i1,i2,i3)
#define PERMUTE4IO(in,base,offset,out,STOREM,LOADM,ADDE,i1,i2,i3,i4) LOADM(i1,in,base,offset);PERMUTE4O(in,out,offset,STOREM,ADDE,i1,i2,i3,i4)
#define PERMUTE5IO(in,base,offset,out,STOREM,LOADM,ADDE,i1,i2,i3,i4,i5) LOADM(i1,in,base,offset);PERMUTE5O(in,out,offset,STOREM,ADDE,i1,i2,i3,i4,i5)
#define PERMUTE6IO(in,base,offset,out,STOREM,LOADM,ADDE,i1,i2,i3,i4,i5,i6) LOADM(i1,in,base,offset);PERMUTE6O(in,out,offset,STOREM,ADDE,i1,i2,i3,i4,i5,i6)
#define PERMUTE7IO(in,base,offset,out,STOREM,LOADM,ADDE,i1,i2,i3,i4,i5,i6,i7) LOADM(i1,in,base,offset);PERMUTE7O(in,out,offset,STOREM,ADDE,i1,i2,i3,i4,i5,i6,i7)
#define PERMUTE8IO(in,base,offset,out,STOREM,LOADM,ADDE,i1,i2,i3,i4,i5,i6,i7,i8) LOADM(i1,in,base,offset);PERMUTE8O(in,out,offset,STOREM,ADDE,i1,i2,i3,i4,i5,i6,i7,i8)

//Different permutations
#define PERMUTE2I1(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE2IO(in,base,offset,out,STOREM,LOADM,ADDE,1,0)
#define PERMUTE2I2(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE2IO(in,base,offset,out,STOREM,LOADM,ADDE,0,1)

#define PERMUTE3I1(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE3IO(in,base,offset,out,STOREM,LOADM,ADDE,2,0,1)
#define PERMUTE3I2(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE3IO(in,base,offset,out,STOREM,LOADM,ADDE,1,2,0)
#define PERMUTE3I3(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE3IO(in,base,offset,out,STOREM,LOADM,ADDE,0,1,2)

#define PERMUTE4I1(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE4IO(in,base,offset,out,STOREM,LOADM,ADDE,3,0,1,2)
#define PERMUTE4I2(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE4IO(in,base,offset,out,STOREM,LOADM,ADDE,2,3,0,1)
#define PERMUTE4I3(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE4IO(in,base,offset,out,STOREM,LOADM,ADDE,1,2,3,0)
#define PERMUTE4I4(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE4IO(in,base,offset,out,STOREM,LOADM,ADDE,0,1,2,3)

#define PERMUTE5I1(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE5IO(in,base,offset,out,STOREM,LOADM,ADDE,4,0,1,2,3)
#define PERMUTE5I2(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE5IO(in,base,offset,out,STOREM,LOADM,ADDE,3,4,0,1,2)
#define PERMUTE5I3(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE5IO(in,base,offset,out,STOREM,LOADM,ADDE,2,3,4,0,1)
#define PERMUTE5I4(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE5IO(in,base,offset,out,STOREM,LOADM,ADDE,1,2,3,4,0)
#define PERMUTE5I5(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE5IO(in,base,offset,out,STOREM,LOADM,ADDE,0,1,2,3,4)

#define PERMUTE6I1(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE6IO(in,base,offset,out,STOREM,LOADM,ADDE,5,0,1,2,3,4)
#define PERMUTE6I2(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE6IO(in,base,offset,out,STOREM,LOADM,ADDE,4,5,0,1,2,3)
#define PERMUTE6I3(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE6IO(in,base,offset,out,STOREM,LOADM,ADDE,3,4,5,0,1,2)
#define PERMUTE6I4(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE6IO(in,base,offset,out,STOREM,LOADM,ADDE,2,3,4,5,0,1)
#define PERMUTE6I5(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE6IO(in,base,offset,out,STOREM,LOADM,ADDE,1,2,3,4,5,0)
#define PERMUTE6I6(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE6IO(in,base,offset,out,STOREM,LOADM,ADDE,0,1,2,3,4,5)

#define PERMUTE7I1(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE7IO(in,base,offset,out,STOREM,LOADM,ADDE,6,0,1,2,3,4,5)
#define PERMUTE7I2(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE7IO(in,base,offset,out,STOREM,LOADM,ADDE,5,6,0,1,2,3,4)
#define PERMUTE7I3(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE7IO(in,base,offset,out,STOREM,LOADM,ADDE,4,5,6,0,1,2,3)
#define PERMUTE7I4(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE7IO(in,base,offset,out,STOREM,LOADM,ADDE,3,4,5,6,0,1,2)
#define PERMUTE7I5(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE7IO(in,base,offset,out,STOREM,LOADM,ADDE,2,3,4,5,6,0,1)
#define PERMUTE7I6(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE7IO(in,base,offset,out,STOREM,LOADM,ADDE,1,2,3,4,5,6,0)
#define PERMUTE7I7(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE7IO(in,base,offset,out,STOREM,LOADM,ADDE,0,1,2,3,4,5,6)

#define PERMUTE8I1(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE8IO(in,base,offset,out,STOREM,LOADM,ADDE,7,0,1,2,3,4,5,6)
#define PERMUTE8I2(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE8IO(in,base,offset,out,STOREM,LOADM,ADDE,6,7,0,1,2,3,4,5)
#define PERMUTE8I3(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE8IO(in,base,offset,out,STOREM,LOADM,ADDE,5,6,7,0,1,2,3,4)
#define PERMUTE8I4(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE8IO(in,base,offset,out,STOREM,LOADM,ADDE,4,5,6,7,0,1,2,3)
#define PERMUTE8I5(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE8IO(in,base,offset,out,STOREM,LOADM,ADDE,3,4,5,6,7,0,1,2)
#define PERMUTE8I6(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE8IO(in,base,offset,out,STOREM,LOADM,ADDE,2,3,4,5,6,7,0,1)
#define PERMUTE8I7(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE8IO(in,base,offset,out,STOREM,LOADM,ADDE,1,2,3,4,5,6,7,0)
#define PERMUTE8I8(in,base,offset,out,STOREM,LOADM,ADDE) PERMUTE8IO(in,base,offset,out,STOREM,LOADM,ADDE,0,1,2,3,4,5,6,7)

#define READELEM1(in,base,offset,LOADM) LOADM(0,in,base,offset)
#define READELEM2(in,base,offset,LOADM) READELEM1(in,base,offset,LOADM); LOADM(1,in,base-1,offset)
#define READELEM3(in,base,offset,LOADM) READELEM2(in,base,offset,LOADM); LOADM(2,in,base-2,offset)
#define READELEM4(in,base,offset,LOADM) READELEM3(in,base,offset,LOADM); LOADM(3,in,base-3,offset)
#define READELEM5(in,base,offset,LOADM) READELEM4(in,base,offset,LOADM); LOADM(4,in,base-4,offset)
#define READELEM6(in,base,offset,LOADM) READELEM5(in,base,offset,LOADM); LOADM(5,in,base-5,offset)
#define READELEM7(in,base,offset,LOADM) READELEM6(in,base,offset,LOADM); LOADM(6,in,base-6,offset)
#define READELEM8(in,base,offset,LOADM) READELEM7(in,base,offset,LOADM); LOADM(7,in,base-7,offset)

//Initializes one level by loading all required elements to fill the local array, and performes the first permutation
#define INITS2(in,out,base,offset,STOREM,LOADM,ADDE) READELEM2(in,base,offset,LOADM); PERMUTE2O(in,out,offset,STOREM,ADDE,0,1)
#define INITS3(in,out,base,offset,STOREM,LOADM,ADDE) READELEM3(in,base,offset,LOADM); PERMUTE3O(in,out,offset,STOREM,ADDE,0,1,2)
#define INITS4(in,out,base,offset,STOREM,LOADM,ADDE) READELEM4(in,base,offset,LOADM); PERMUTE4O(in,out,offset,STOREM,ADDE,0,1,2,3)
#define INITS5(in,out,base,offset,STOREM,LOADM,ADDE) READELEM5(in,base,offset,LOADM); PERMUTE5O(in,out,offset,STOREM,ADDE,0,1,2,3,4)
#define INITS6(in,out,base,offset,STOREM,LOADM,ADDE) READELEM6(in,base,offset,LOADM); PERMUTE6O(in,out,offset,STOREM,ADDE,0,1,2,3,4,5)
#define INITS7(in,out,base,offset,STOREM,LOADM,ADDE) READELEM7(in,base,offset,LOADM); PERMUTE7O(in,out,offset,STOREM,ADDE,0,1,2,3,4,5,6)
#define INITS8(in,out,base,offset,STOREM,LOADM,ADDE) READELEM8(in,base,offset,LOADM); PERMUTE8O(in,out,offset,STOREM,ADDE,0,1,2,3,4,5,6,7)

//One block is when all elements are back to their original position. Always requires UNROLLLEVEL+1 shifts
#define BLOCK2S1(in,base,offset,STOREM,LOADM,ADDE,DN) PERMUTE2I1(in,base,offset,base-1+DN,STOREM,LOADM,ADDE)
#define BLOCK2(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK2S1(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE2I2(in,base+1,offset,base+DN,STOREM,LOADM,ADDE)

#define BLOCK3S1(in,base,offset,STOREM,LOADM,ADDE,DN) PERMUTE3I1(in,base,offset,base-2+DN,STOREM,LOADM,ADDE)
#define BLOCK3S2(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK3S1(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE3I2(in,base+1,offset,base-1+DN,STOREM,LOADM,ADDE)
#define BLOCK3(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK3S2(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE3I3(in,base+2,offset,base+DN,STOREM,LOADM,ADDE)

#define BLOCK4S1(in,base,offset,STOREM,LOADM,ADDE,DN) PERMUTE4I1(in,base,offset,base-3+DN,STOREM,LOADM,ADDE)
#define BLOCK4S2(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK4S1(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE4I2(in,base+1,offset,base-2+DN,STOREM,LOADM,ADDE)
#define BLOCK4S3(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK4S2(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE4I3(in,base+2,offset,base-1+DN,STOREM,LOADM,ADDE)
#define BLOCK4(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK4S3(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE4I4(in,base+3,offset,base+DN,STOREM,LOADM,ADDE)

#define BLOCK5S1(in,base,offset,STOREM,LOADM,ADDE,DN) PERMUTE5I1(in,base,offset,base-4+DN,STOREM,LOADM,ADDE)
#define BLOCK5S2(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK5S1(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE5I2(in,base+1,offset,base-3+DN,STOREM,LOADM,ADDE)
#define BLOCK5S3(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK5S2(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE5I3(in,base+2,offset,base-2+DN,STOREM,LOADM,ADDE)
#define BLOCK5S4(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK5S3(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE5I4(in,base+3,offset,base-1+DN,STOREM,LOADM,ADDE)
#define BLOCK5(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK5S4(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE5I5(in,base+4,offset,base+DN,STOREM,LOADM,ADDE)

#define BLOCK6S1(in,base,offset,STOREM,LOADM,ADDE,DN) PERMUTE6I1(in,base,offset,base-5+DN,STOREM,LOADM,ADDE)
#define BLOCK6S2(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK6S1(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE6I2(in,base+1,offset,base-4+DN,STOREM,LOADM,ADDE)
#define BLOCK6S3(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK6S2(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE6I3(in,base+2,offset,base-3+DN,STOREM,LOADM,ADDE)
#define BLOCK6S4(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK6S3(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE6I4(in,base+3,offset,base-2+DN,STOREM,LOADM,ADDE)
#define BLOCK6S5(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK6S4(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE6I5(in,base+4,offset,base-1+DN,STOREM,LOADM,ADDE)
#define BLOCK6(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK6S5(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE6I6(in,base+5,offset,base+DN,STOREM,LOADM,ADDE)

#define BLOCK7S1(in,base,offset,STOREM,LOADM,ADDE,DN) PERMUTE7I1(in,base,offset,base-6+DN,STOREM,LOADM,ADDE)
#define BLOCK7S2(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK7S1(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE7I2(in,base+1,offset,base-5+DN,STOREM,LOADM,ADDE)
#define BLOCK7S3(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK7S2(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE7I3(in,base+2,offset,base-4+DN,STOREM,LOADM,ADDE)
#define BLOCK7S4(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK7S3(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE7I4(in,base+3,offset,base-3+DN,STOREM,LOADM,ADDE)
#define BLOCK7S5(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK7S4(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE7I5(in,base+4,offset,base-2+DN,STOREM,LOADM,ADDE)
#define BLOCK7S6(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK7S5(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE7I6(in,base+5,offset,base-1+DN,STOREM,LOADM,ADDE)
#define BLOCK7(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK7S6(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE7I7(in,base+6,offset,base+DN,STOREM,LOADM,ADDE)

#define BLOCK8S1(in,base,offset,STOREM,LOADM,ADDE,DN) PERMUTE8I1(in,base,offset,base-7+DN,STOREM,LOADM,ADDE)
#define BLOCK8S2(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK8S1(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE8I2(in,base+1,offset,base-6+DN,STOREM,LOADM,ADDE)
#define BLOCK8S3(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK8S2(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE8I3(in,base+2,offset,base-5+DN,STOREM,LOADM,ADDE)
#define BLOCK8S4(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK8S3(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE8I4(in,base+3,offset,base-4+DN,STOREM,LOADM,ADDE)
#define BLOCK8S5(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK8S4(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE8I5(in,base+4,offset,base-3+DN,STOREM,LOADM,ADDE)
#define BLOCK8S6(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK8S5(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE8I6(in,base+5,offset,base-2+DN,STOREM,LOADM,ADDE)
#define BLOCK8S7(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK8S6(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE8I7(in,base+6,offset,base-1+DN,STOREM,LOADM,ADDE)
#define BLOCK8(in,base,offset,STOREM,LOADM,ADDE,DN) BLOCK8S7(in,base,offset,STOREM,LOADM,ADDE,DN); PERMUTE8I8(in,base+7,offset,base+DN,STOREM,LOADM,ADDE)

//Ending code to finish one level. Differs depending on where in the permutation the code is
#define BLOCK2END1BASE(in,base,offset,STOREM,ADDE,DN,O1,O2) PERMUTE2O(in,base-1+DN,offset,STOREM,ADDE,O1,1+O2)
#define BLOCK2END2BASE(in,base,offset,STOREM,ADDE,DN,O1,O2) PERMUTE2O(in,base-1+DN,offset,STOREM,ADDE,1+O1,O2)
#define BLOCK2END1(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK2END1BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET1)
#define BLOCK2END2(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK2END2BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1+OFFSET2,OFFSET1)

#define BLOCK3END1BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3) PERMUTE3O(in,base-2+DN,offset,STOREM,ADDE,O1,1+O2,2+O3); BLOCK2END1BASE(in,base,offset,STOREM,ADDE,DN,O1,O2)
#define BLOCK3END2BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3) PERMUTE3O(in,base-2+DN,offset,STOREM,ADDE,2+O1,O2,1+O3); BLOCK2END2BASE(in,base,offset,STOREM,ADDE,DN,O1+1,O2)
#define BLOCK3END3BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3) PERMUTE3O(in,base-2+DN,offset,STOREM,ADDE,1+O1,2+O2,O3); BLOCK2END1BASE(in,base,offset,STOREM,ADDE,DN,O1+1,O2+1)
#define BLOCK3END1(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK3END1BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET1,OFFSET1)
#define BLOCK3END2(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK3END2BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1+OFFSET2,OFFSET1,OFFSET1)
#define BLOCK3END3(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK3END3BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1)

#define BLOCK4END1BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4) PERMUTE4O(in,base-3+DN,offset,STOREM,ADDE,O1,1+O2,2+O3,3+O4); BLOCK3END1BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3)
#define BLOCK4END2BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4) PERMUTE4O(in,base-3+DN,offset,STOREM,ADDE,3+O1,O2,1+O3,2+O4); BLOCK3END2BASE(in,base,offset,STOREM,ADDE,DN,O1+1,O2,O3)
#define BLOCK4END3BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4) PERMUTE4O(in,base-3+DN,offset,STOREM,ADDE,2+O1,3+O2,O3,1+O4); BLOCK3END3BASE(in,base,offset,STOREM,ADDE,DN,O1+1,O2+1,O3)
#define BLOCK4END4BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4) PERMUTE4O(in,base-3+DN,offset,STOREM,ADDE,1+O1,2+O2,3+O3,O4); BLOCK3END1BASE(in,base,offset,STOREM,ADDE,DN,O1+1,O2+1,O3+1)
#define BLOCK4END1(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK4END1BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET1,OFFSET1,OFFSET1)
#define BLOCK4END2(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK4END2BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1+OFFSET2,OFFSET1,OFFSET1,OFFSET1)
#define BLOCK4END3(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK4END3BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1,OFFSET1)
#define BLOCK4END4(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK4END4BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1)

#define BLOCK5END1BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4,O5) PERMUTE5O(in,base-4+DN,offset,STOREM,ADDE,O1,1+O2,2+O3,3+O4,4+O5); BLOCK4END1BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4)
#define BLOCK5END2BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4,O5) PERMUTE5O(in,base-4+DN,offset,STOREM,ADDE,4+O1,O2,1+O3,2+O4,3+O5); BLOCK4END2BASE(in,base,offset,STOREM,ADDE,DN,O1+1,O2,O3,O4)
#define BLOCK5END3BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4,O5) PERMUTE5O(in,base-4+DN,offset,STOREM,ADDE,3+O1,4+O2,O3,1+O4,2+O5); BLOCK4END3BASE(in,base,offset,STOREM,ADDE,DN,O1+1,O2+1,O3,O4)
#define BLOCK5END4BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4,O5) PERMUTE5O(in,base-4+DN,offset,STOREM,ADDE,2+O1,3+O2,4+O3,O4,1+O5); BLOCK4END4BASE(in,base,offset,STOREM,ADDE,DN,O1+1,O2+1,O3+1,O4)
#define BLOCK5END5BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4,O5) PERMUTE5O(in,base-4+DN,offset,STOREM,ADDE,1+O1,2+O2,3+O3,4+O4,O5); BLOCK4END1BASE(in,base,offset,STOREM,ADDE,DN,O1+1,O2+1,O3+1,O4+1)
#define BLOCK5END1(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK5END1BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET1,OFFSET1,OFFSET1,OFFSET1) 
#define BLOCK5END2(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK5END2BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1+OFFSET2,OFFSET1,OFFSET1,OFFSET1,OFFSET1)
#define BLOCK5END3(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK5END3BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1,OFFSET1,OFFSET1)
#define BLOCK5END4(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK5END4BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1,OFFSET1)
#define BLOCK5END5(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK5END5BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1)

#define BLOCK6END1BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4,O5,O6) PERMUTE6O(in,base-5+DN,offset,STOREM,ADDE,O1,1+O2,2+O3,3+O4,4+O5,5+O6); BLOCK5END1BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4,O5)
#define BLOCK6END2BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4,O5,O6) PERMUTE6O(in,base-5+DN,offset,STOREM,ADDE,5+O1,O2,1+O3,2+O4,3+O5,4+O6); BLOCK5END2BASE(in,base,offset,STOREM,ADDE,DN,O1+1,O2,O3,O4,O5)
#define BLOCK6END3BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4,O5,O6) PERMUTE6O(in,base-5+DN,offset,STOREM,ADDE,4+O1,5+O2,O3,1+O4,2+O5,3+O6); BLOCK5END3BASE(in,base,offset,STOREM,ADDE,DN,O1+1,O2+1,O3,O4,O5)
#define BLOCK6END4BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4,O5,O6) PERMUTE6O(in,base-5+DN,offset,STOREM,ADDE,3+O1,4+O2,5+O3,O4,1+O5,2+O6); BLOCK5END4BASE(in,base,offset,STOREM,ADDE,DN,O1+1,O2+1,O3+1,O4,O5)
#define BLOCK6END5BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4,O5,O6) PERMUTE6O(in,base-5+DN,offset,STOREM,ADDE,2+O1,3+O2,4+O3,5+O4,O5,1+O6); BLOCK5END5BASE(in,base,offset,STOREM,ADDE,DN,O1+1,O2+1,O3+1,O4+1,O5)
#define BLOCK6END6BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4,O5,O6) PERMUTE6O(in,base-5+DN,offset,STOREM,ADDE,1+O1,2+O2,3+O3,4+O4,5+O5,O6); BLOCK5END1BASE(in,base,offset,STOREM,ADDE,DN,O1+1,O2+1,O3+1,O4+1,O5+1)
#define BLOCK6END1(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK6END1BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET1,OFFSET1,OFFSET1,OFFSET1,OFFSET1)
#define BLOCK6END2(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK6END2BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1+OFFSET2,OFFSET1,OFFSET1,OFFSET1,OFFSET1,OFFSET1)
#define BLOCK6END3(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK6END3BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1,OFFSET1,OFFSET1,OFFSET1)
#define BLOCK6END4(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK6END4BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1,OFFSET1,OFFSET1)
#define BLOCK6END5(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK6END5BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1,OFFSET1)
#define BLOCK6END6(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK6END6BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1)

#define BLOCK7END1BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4,O5,O6,O7) PERMUTE7O(in,base-6+DN,offset,STOREM,ADDE,O1,1+O2,2+O3,3+O4,4+O5,5+O6,6+O7); BLOCK6END1BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4,O5,O6)
#define BLOCK7END2BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4,O5,O6,O7) PERMUTE7O(in,base-6+DN,offset,STOREM,ADDE,6+O1,O2,1+O3,2+O4,3+O5,4+O6,5+O7); BLOCK6END2BASE(in,base,offset,STOREM,ADDE,DN,O1+1,O2,O3,O4,O5,O6)
#define BLOCK7END3BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4,O5,O6,O7) PERMUTE7O(in,base-6+DN,offset,STOREM,ADDE,5+O1,6+O2,O3,1+O4,2+O5,3+O6,4+O7); BLOCK6END3BASE(in,base,offset,STOREM,ADDE,DN,O1+1,O2+1,O3,O4,O5,O6)
#define BLOCK7END4BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4,O5,O6,O7) PERMUTE7O(in,base-6+DN,offset,STOREM,ADDE,4+O1,5+O2,6+O3,O4,1+O5,2+O6,3+O7); BLOCK6END4BASE(in,base,offset,STOREM,ADDE,DN,O1+1,O2+1,O3+1,O4,O5,O6)
#define BLOCK7END5BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4,O5,O6,O7) PERMUTE7O(in,base-6+DN,offset,STOREM,ADDE,3+O1,4+O2,5+O3,6+O4,O5,1+O6,2+O7); BLOCK6END5BASE(in,base,offset,STOREM,ADDE,DN,O1+1,O2+1,O3+1,O4+1,O5,O6)
#define BLOCK7END6BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4,O5,O6,O7) PERMUTE7O(in,base-6+DN,offset,STOREM,ADDE,2+O1,3+O2,4+O3,5+O4,6+O5,O6,1+O7); BLOCK6END6BASE(in,base,offset,STOREM,ADDE,DN,O1+1,O2+1,O3+1,O4+1,O5+1,O6)
#define BLOCK7END7BASE(in,base,offset,STOREM,ADDE,DN,O1,O2,O3,O4,O5,O6,O7) PERMUTE7O(in,base-6+DN,offset,STOREM,ADDE,1+O1,2+O2,3+O3,4+O4,5+O5,6+O6,O7); BLOCK6END1BASE(in,base,offset,STOREM,ADDE,DN,O1+1,O2+1,O3+1,O4+1,O5+1,O6+1)
#define BLOCK7END1(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK7END1BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET1,OFFSET1,OFFSET1,OFFSET1,OFFSET1,OFFSET1)
#define BLOCK7END2(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK7END2BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1+OFFSET2,OFFSET1,OFFSET1,OFFSET1,OFFSET1,OFFSET1,OFFSET1)
#define BLOCK7END3(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK7END3BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1,OFFSET1,OFFSET1,OFFSET1,OFFSET1)
#define BLOCK7END4(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK7END4BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1,OFFSET1,OFFSET1,OFFSET1)
#define BLOCK7END5(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK7END5BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1,OFFSET1,OFFSET1)
#define BLOCK7END6(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK7END6BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1,OFFSET1)
#define BLOCK7END7(in,base,offset,STOREM,ADDE,DN,OFFSET1,OFFSET2) BLOCK7END7BASE(in,base,offset,STOREM,ADDE,DN,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1+OFFSET2,OFFSET1)

//Basic code for inner loop. One version for each number of shifts to perform
#define INNERUNROLLBASE2(in,offset,STOREM,LOADM,ADDE,DN,startj) j = startj+1-DN;\
INITS2(in,j-1+DN,j,offset,STOREM,LOADM,ADDE);j++;\
for (; j <= pshift-7; j+=8) {\
  BLOCK2(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK2(in,j+2,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK2(in,j+4,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK2(in,j+6,offset,STOREM,LOADM,ADDE,DN);\
}\
if (j <= pshift-3) {\
  BLOCK2(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK2(in,j+2,offset,STOREM,LOADM,ADDE,DN);\
  j+=4;\
}\
if (j <= pshift-2) {\
  BLOCK2(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK2S1(in,j+2,offset,STOREM,LOADM,ADDE,DN);\
}\
else if (j <= pshift-1) {\
  BLOCK2(in,j,offset,STOREM,LOADM,ADDE,DN);\
}\
else if (j <= pshift) {\
  BLOCK2S1(in,j,offset,STOREM,LOADM,ADDE,DN);\
}
    
#define INNERUNROLLBASE3(in,offset,STOREM,LOADM,ADDE,DN,startj) j = startj+1-DN;\
INITS3(in,j-2+DN,j,offset,STOREM,LOADM,ADDE);j++;\
for (; j <= pshift-5; j+=6) {\
  BLOCK3(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK3(in,j+3,offset,STOREM,LOADM,ADDE,DN);\
}\
if (j <= pshift-2) {\
  BLOCK3(in,j,offset,STOREM,LOADM,ADDE,DN);\
  j+=3;\
}\
if (j <= pshift-1) {\
  BLOCK3S2(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK2END1(in,j+1,offset,STOREM,ADDE,DN,1,0);\
}\
else if (j <= pshift) {\
  BLOCK3S1(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK2END2(in,j,offset,STOREM,ADDE,DN,0,1);\
}\
else {\
  BLOCK2END1(in,j-1,offset,STOREM,ADDE,DN,0,0);\
}
        
#define INNERUNROLLBASE4(in,offset,STOREM,LOADM,ADDE,DN,startj) j = startj+1-DN;\
INITS4(in,j-3+DN,j,offset,STOREM,LOADM,ADDE);j++;\
for (; j <= pshift-7; j+=8) {\
  BLOCK4(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK4(in,j+4,offset,STOREM,LOADM,ADDE,DN);\
}\
if (j <= pshift-3) {\
  BLOCK4(in,j,offset,STOREM,LOADM,ADDE,DN);\
  j+=4;\
}\
if (j <= pshift-2) {\
  BLOCK4S3(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK3END1(in,j+2,offset,STOREM,ADDE,DN,1,0);\
}\
else if (j <= pshift-1) {\
  BLOCK4S2(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK3END3(in,j+1,offset,STOREM,ADDE,DN,0,1);\
}\
else if (j <= pshift) {\
  BLOCK4S1(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK3END2(in,j,offset,STOREM,ADDE,DN,0,1);\
}\
else {\
  BLOCK3END1(in,j-1,offset,STOREM,ADDE,DN,0,0);\
}

#define INNERUNROLLBASE5(in,offset,STOREM,LOADM,ADDE,DN,startj) j = startj+1-DN;\
INITS5(in,j-4+DN,j,offset,STOREM,LOADM,ADDE);j++;\
for (; j <= pshift-4; j+=5) {\
  BLOCK5(in,j,offset,STOREM,LOADM,ADDE,DN);\
}\
if (j <= pshift-3) {\
  BLOCK5S4(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK4END1(in,j+3,offset,STOREM,ADDE,DN,1,0);\
}\
  else if (j <= pshift-2) {\
  BLOCK5S3(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK4END4(in,j+2,offset,STOREM,ADDE,DN,0,1);\
}\
else if (j <= pshift-1) {\
  BLOCK5S2(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK4END3(in,j+1,offset,STOREM,ADDE,DN,0,1);\
}\
else if (j <= pshift) {\
  BLOCK5S1(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK4END2(in,j,offset,STOREM,ADDE,DN,0,1);\
}\
else {\
  BLOCK4END1(in,j-1,offset,STOREM,ADDE,DN,0,0);\
}

#define INNERUNROLLBASE6(in,offset,STOREM,LOADM,ADDE,DN,startj) j = startj+1-DN;\
INITS6(in,j-5+DN,j,offset,STOREM,LOADM,ADDE);j++;\
for (; j <= pshift-5; j+=6) {\
  BLOCK6(in,j,offset,STOREM,LOADM,ADDE,DN);\
}\
if (j <= pshift-4) {\
  BLOCK6S5(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK5END1(in,j+4,offset,STOREM,ADDE,DN,1,0);\
}\
else if (j <= pshift-3) {\
  BLOCK6S4(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK5END5(in,j+3,offset,STOREM,ADDE,DN,0,1);\
}\
  else if (j <= pshift-2) {\
  BLOCK6S3(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK5END4(in,j+2,offset,STOREM,ADDE,DN,0,1);\
}\
else if (j <= pshift-1) {\
  BLOCK6S2(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK5END3(in,j+1,offset,STOREM,ADDE,DN,0,1);\
}\
else if (j <= pshift) {\
  BLOCK6S1(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK5END2(in,j,offset,STOREM,ADDE,DN,0,1);\
}\
else {\
  BLOCK5END1(in,j-1,offset,STOREM,ADDE,DN,0,0);\
}

#define INNERUNROLLBASE7(in,offset,STOREM,LOADM,ADDE,DN,startj) j = startj+1-DN;\
INITS7(in,j-6+DN,j,offset,STOREM,LOADM,ADDE);j++;\
for (; j <= pshift-6; j+=7) {\
  BLOCK7(in,j,offset,STOREM,LOADM,ADDE,DN);\
}\
if (j <= pshift-5) {\
  BLOCK7S6(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK6END1(in,j+5,offset,STOREM,ADDE,DN,1,0);\
}\
else if (j <= pshift-4) {\
  BLOCK7S5(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK6END6(in,j+4,offset,STOREM,ADDE,DN,0,1);\
}\
else if (j <= pshift-3) {\
  BLOCK7S4(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK6END5(in,j+3,offset,STOREM,ADDE,DN,0,1);\
}\
  else if (j <= pshift-2) {\
  BLOCK7S3(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK6END4(in,j+2,offset,STOREM,ADDE,DN,0,1);\
}\
else if (j <= pshift-1) {\
  BLOCK7S2(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK6END3(in,j+1,offset,STOREM,ADDE,DN,0,1);\
}\
else if (j <= pshift) {\
  BLOCK7S1(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK6END2(in,j,offset,STOREM,ADDE,DN,0,1);\
}\
else {\
  BLOCK6END1(in,j-1,offset,STOREM,ADDE,DN,0,0);\
}

#define INNERUNROLLBASE8(in,offset,STOREM,LOADM,ADDE,DN,startj) j = startj+1-DN;\
/*mexPrintf("start pshift=%d j=%d,k=%d\n",pshift,j,k);*/\
INITS8(in,j-7+DN,j,offset,STOREM,LOADM,ADDE);j++;\
for (; j <= pshift-7; j+=8) {\
  BLOCK8(in,j,offset,STOREM,LOADM,ADDE,DN);\
  /*mexPrintf("loop pshift=%d j=%d,k=%d\n",pshift,j,k);*/\
}\
if (j <= pshift-6) {\
  BLOCK8S7(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK7END1(in,j+6,offset,STOREM,ADDE,DN,1,0);\
  /*mexPrintf("j <= pshift-6 pshift=%d j=%d,k=%d\n",pshift,j,k);*/\
}\
else if (j <= pshift-5) {\
  BLOCK8S6(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK7END7(in,j+5,offset,STOREM,ADDE,DN,0,1);\
  /*mexPrintf("j <= pshift-5 pshift=%d j=%d,k=%d\n",pshift,j,k);*/\
}\
else if (j <= pshift-4) {\
  BLOCK8S5(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK7END6(in,j+4,offset,STOREM,ADDE,DN,0,1);\
  /*mexPrintf("j <= pshift-4 pshift=%d j=%d,k=%d\n",pshift,j,k);*/\
}\
else if (j <= pshift-3) {\
  BLOCK8S4(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK7END5(in,j+3,offset,STOREM,ADDE,DN,0,1);\
  /*mexPrintf("j <= pshift-3 pshift=%d j=%d,k=%d\n",pshift,j,k);*/\
}\
else if (j <= pshift-2) {\
  BLOCK8S3(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK7END4(in,j+2,offset,STOREM,ADDE,DN,0,1);\
  /*mexPrintf("j <= pshift-2 pshift=%d j=%d,k=%d\n",pshift,j,k);*/\
}\
else if (j <= pshift-1) {\
  BLOCK8S2(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK7END3(in,j+1,offset,STOREM,ADDE,DN,0,1);\
  /*mexPrintf("j <= pshift-1 pshift=%d j=%d,k=%d\n",pshift,j,k);*/\
}\
else if (j <= pshift) {\
  BLOCK8S1(in,j,offset,STOREM,LOADM,ADDE,DN);\
  BLOCK7END2(in,j,offset,STOREM,ADDE,DN,0,1);\
  /*mexPrintf("j <= pshift pshift=%d j=%d,k=%d\n",pshift,j,k);*/\
}\
else {\
  BLOCK7END1(in,j-1,offset,STOREM,ADDE,DN,0,0);\
  /*mexPrintf("else pshift=%d j=%d,k=%d\n",pshift,j,k);*/\
}
