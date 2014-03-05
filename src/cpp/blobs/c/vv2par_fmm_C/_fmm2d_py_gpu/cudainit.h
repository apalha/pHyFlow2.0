#define initthreadperblock 64
#define initnrcoeff 4
#define initoffset (initthreadperblock+1)
#define INITUPDATECOEFF for(k=1;k<8;k++) {/*make an intermediate step to keep better load balance between the threads, 8 points in coeff matrix for each thread otherwise*/\
  coeffmatrix[dest3]+=coeffmatrix[dest3+k];\
}\
__syncthreads();\
if(threadIdx.x<(blockDim.x>>1)) { /*a second intermediate step*/\
  coeffmatrix[dest2]+=coeffmatrix[dest2+8];\
}\
if(threadIdx.x<(blockDim.x>>2)) { /*a second intermediate step*/\
  coeffmatrix[dest]+=coeffmatrix[dest+16];\
}\
__syncthreads();\
if(threadIdx.x<(blockDim.x>>3)) { /*write results*/\
  coeffmatrix[threadIdx.x*initoffset]+=coeffmatrix[threadIdx.x*initoffset+32];\
  if(isfinite(coeffmatrix[threadIdx.x*initoffset]))\
    coeff[2*j+threadIdx.x]-=coeffmatrix[threadIdx.x*initoffset];\
}

#define INITUPDATECOEFFWITHCHECK if(threadIdx.x<innercount*16+16) { /*make an intermediate step to keep better load balance between the threads, 8 points in coeff matrix for each thread otherwise*/\
  for(k=1;k<8;k++) {/*make an intermediate step to keep better load balance between the threads, 8 points in coeff matrix for each thread otherwise*/\
    coeffmatrix[dest3]+=coeffmatrix[dest3+k];\
  }\
}\
__syncthreads();\
if(threadIdx.x<innercount*8+8) {  /*a second intermediate step*/\
  coeffmatrix[dest2]+=coeffmatrix[dest2+8];\
}\
__syncthreads();\
if(threadIdx.x<innercount*4+4) {  /*a second intermediate step*/\
  coeffmatrix[dest]+=coeffmatrix[dest+16];\
}\
__syncthreads();\
if(threadIdx.x<innercount*2+2) { /*write results*/\
  coeffmatrix[threadIdx.x*initoffset]+=coeffmatrix[threadIdx.x*initoffset+32];\
  if(isfinite(coeffmatrix[threadIdx.x*initoffset]))\
    coeff[2*j+threadIdx.x]-=coeffmatrix[threadIdx.x*initoffset];\
}

//the ordering is 0,2,1,3,4,6,5,7 etc from coeffmatrix to coeff
#define INITUPDATECOEFF2 for(k=1;k<8;k++) {/*make an intermediate step to keep better load balance between the threads, 8 points in coeff matrix for each thread otherwise*/\
  coeffmatrix[dest3]+=coeffmatrix[dest3+k];\
}\
__syncthreads();\
if(threadIdx.x<(blockDim.x>>1)) { /*a second intermediate step*/\
  coeffmatrix[dest2]+=coeffmatrix[dest2+8];\
}\
if(threadIdx.x<(blockDim.x>>2)) { /*a second intermediate step*/\
  coeffmatrix[dest]+=coeffmatrix[dest+16];\
  if(isfinite(coeffmatrix[dest]))\
    coeff[2*j+threadIdx.x+(threadIdx.x&1)-((threadIdx.x&2)>>1)]-=coeffmatrix[dest];\
}

#define INITUPDATECOEFFWITHCHECK2 if(threadIdx.x<innercount*8+16) { /*make an intermediate step to keep better load balance between the threads, 8 points in coeff matrix for each thread otherwise*/\
  for(k=1;k<8;k++) {/*make an intermediate step to keep better load balance between the threads, 8 points in coeff matrix for each thread otherwise*/\
    coeffmatrix[dest3]+=coeffmatrix[dest3+k];\
  }\
}\
__syncthreads();\
if(threadIdx.x<(blockDim.x>>1)) { /*a second intermediate step*/\
  coeffmatrix[dest2]+=coeffmatrix[dest2+8];\
}\
__syncthreads();\
if(threadIdx.x<(blockDim.x>>2)) { /*a second intermediate step*/\
  coeffmatrix[dest]+=coeffmatrix[dest+16];\
  k=2*j+threadIdx.x+(threadIdx.x&1)-((threadIdx.x&2)>>1);\
  if(isfinite(coeffmatrix[dest])&&k<(p+1)*2)\
    coeff[k]-=coeffmatrix[dest];\
}


//Initiates the coefficients at the lowers level from the potential points. 
//Different threads take different points and write results to a temporary array, 
//which is later summed to the coefficient matrix. This is necessary, since all points
//has to write their coefficients to the same array
extern __shared__ double coeff[];
__global__ void cuda_mpexp_init(int p, double *coeff1, double *coeff2, int Nf,int fetchbase,int *j2cptr,int *kcptr,int *ir,int complexpoint,int pot DEBUGVECTORSTRING)
{
  __shared__ double coeffmatrix[initoffset*initnrcoeff*2];
  __shared__ int count,loopmax,m,m1,base,innercount,distbox; /*save space, same for all threads*/
  __shared__ SORT_DCMPLX z0;
  int boxnr;
  int dest,dest2,dest3;
  double re,im,zre0,zre,zim;
  int i,j,k,multiple;
  for(boxnr=blockIdx.x;boxnr<Nf;boxnr+=gridDim.x) { //The loop over the boxes
    if(threadIdx.x==0) {
      m=j2cptr[boxnr];//distant box interactions
      m1=kcptr[boxnr];
      count=tex1Dfetch(tixptr, boxnr+1); /*right now, it is possible that two boxes has same evalbegin if one box has 0 elements*/
      loopmax=count-(count+blockDim.x-1-(tex1Dfetch(tixptr,boxnr)%blockDim.x))%blockDim.x-1+blockDim.x;//All threads need to loop the same number of times due to __syncthreads()
      base=(boxnr+fetchbase)*(p+1)*sizeof(dcmplx)/sizeof(double);
    }
    if(threadIdx.x<sizeof(SORT_DCMPLX)/sizeof(float)) { //get the box position
      ((float*)&z0)[threadIdx.x] =tex1Dfetch(tz0,(boxnr+fetchbase)*sizeof(SORT_DCMPLX)/sizeof(float)+threadIdx.x);
    }
    __syncthreads();
    if(p<127) { //start by initiating the coefficients in the current box
      for(j=threadIdx.x;j<(p+1)*2;j+=blockDim.x) {
        coeff[j]=0;      
      }
      dest=(threadIdx.x>>1)*initoffset+(threadIdx.x&1)*32; //save this value, since it is used frequently
      dest2=(threadIdx.x>>2)*initoffset+(threadIdx.x&3)*16; //save this value, since it is used frequently
      dest3=(threadIdx.x>>3)*initoffset+(threadIdx.x&7)*8; //save this value, since it is used frequently
//             //__syncthreads(); /*should not be necessary*/
      for(i=threadIdx.x+tex1Dfetch(tixptr, boxnr);i<loopmax;i+=blockDim.x) { /*if more points than threads*/
        multiple=0;
        if(count-i+threadIdx.x<blockDim.x/2) {
          multiple=1;
          if(threadIdx.x>=blockDim.x/2)
            i-=blockDim.x/2;
        }
        if(i<count) { //if point belongs to box, read the input for this point and save in local memory
          TEXTUREFETCH(re, tzr, i);
          TEXTUREFETCH(im, tzi, i);
          TEXTUREFETCH(zre, tmr, i);
          re -= creal(z0);
          im -= cimag(z0);
          if(complexpoint) {
            TEXTUREFETCH(zim, tmi, i);
          }
          else {
            zim = 0;
          }
        }
        for(k=0;k<8;k++)
          coeffmatrix[initoffset*k+threadIdx.x]=0; //init entire matrix to 0 to allow for additions with extra elements
        if(multiple) {
          if(threadIdx.x>=blockDim.x/2) {
            zre0 = zre*re-zim*im;
            zim = zim*re+zre*im;
            zre = zre0;
          }
          zre0=re*re-im*im;
          im=2*re*im;
          re=zre0;
//           debugvector[threadIdx.x+blockIdx.x*gridDim.x]=1;
          for(j=pot;j<=p-initnrcoeff*2;j+=initnrcoeff*2) { //due to limited memory, do this in loops
            //first, make one pass calculating the coffient values and store them in a large matrix, then sum the matrix. This is because all points contribute to the same coefficien values
            if(i<count) {
              if(pot==0) { //special case for pot==0
                //The same algorithm as for the CPU case, but temporary save the results in coeffmatrix instead
                if(j==0) {
                  if(threadIdx.x<blockDim.x/2) {
                    coeffmatrix[threadIdx.x]=zre;
                    coeffmatrix[initoffset+threadIdx.x]=zim;
                  }
                  else {
                    coeffmatrix[threadIdx.x]=-zre;
                    coeffmatrix[initoffset+threadIdx.x]=-zim;
                  }
                  zre0 = zre*re-zim*im;
                  zim = zim*re+zre*im;
                  zre = zre0;
                  for(k=1;k<initnrcoeff;k++) {
                    coeffmatrix[k*2*initoffset+threadIdx.x]=-zre/(j+k*2+(threadIdx.x>=blockDim.x/2));
                    coeffmatrix[(k*2+1)*initoffset+threadIdx.x]=-zim/(j+k*2+(threadIdx.x>=blockDim.x/2));
                    zre0 = zre*re-zim*im;
                    zim = zim*re+zre*im;
                    zre = zre0;
                  }
                }
                else {
                  for(k=0;k<initnrcoeff;k++) {
                    coeffmatrix[k*2*initoffset+threadIdx.x]=-zre/(j+k*2+(threadIdx.x>=blockDim.x/2));
                    coeffmatrix[(k*2+1)*initoffset+threadIdx.x]=-zim/(j+k*2+(threadIdx.x>=blockDim.x/2));
                    zre0 = zre*re-zim*im;
                    zim = zim*re+zre*im;
                    zre = zre0;
                  }
                }
              }
              else { //pot==1

                //pass 1, calculate coeff matrix
                for(k=0;k<initnrcoeff;k++) {
                  coeffmatrix[k*2*initoffset+threadIdx.x]=zre;
                  coeffmatrix[(k*2+1)*initoffset+threadIdx.x]=zim;
                  zre0 = zre*re-zim*im;
                  zim = zim*re+zre*im;
                  zre = zre0;
                }
              }
            }
            __syncthreads();

            //make the summation by letting different threads take different points in the coeff array instead
            INITUPDATECOEFF2
            __syncthreads();
          }

        }
        else {
          for(j=pot;j<=p-initnrcoeff;j+=initnrcoeff) { //due to limited memory, do this in loops
            //first, make one pass calculating the coffient values and store them in a large matrix, then sum the matrix. This is because all points contribute to the same coefficien values
            if(i<count) {
              if(pot==0) { //special case for pot==0
                //The same algorithm as for the CPU case, but temporary save the results in coeffmatrix instead
                if(j==0) {
                  coeffmatrix[threadIdx.x]=zre;
                  coeffmatrix[initoffset+threadIdx.x]=zim;
                  zre0 = zre*re-zim*im;
                  zim = zim*re+zre*im;
                  zre = zre0;
                  coeffmatrix[2*initoffset+threadIdx.x]=-zre;
                  coeffmatrix[(1*2+1)*initoffset+threadIdx.x]=-zim;
                  zre0 = zre*re-zim*im;
                  zim = zim*re+zre*im;
                  zre = zre0;
                  for(k=2;k<initnrcoeff;k++) {
                    coeffmatrix[k*2*initoffset+threadIdx.x]=-zre/(j+k);
                    coeffmatrix[(k*2+1)*initoffset+threadIdx.x]=-zim/(j+k);
                    zre0 = zre*re-zim*im;
                    zim = zim*re+zre*im;
                    zre = zre0;
                  }
                }
                else {
                  for(k=0;k<initnrcoeff;k++) {
                    coeffmatrix[k*2*initoffset+threadIdx.x]=-zre/(j+k);
                    coeffmatrix[(k*2+1)*initoffset+threadIdx.x]=-zim/(j+k);
                    zre0 = zre*re-zim*im;
                    zim = zim*re+zre*im;
                    zre = zre0;
                  }
                }
              }
              else { //pot==1

                //pass 1, calculate coeff matrix
                for(k=0;k<initnrcoeff;k++) {
                  coeffmatrix[k*2*initoffset+threadIdx.x]=zre;
                  coeffmatrix[(k*2+1)*initoffset+threadIdx.x]=zim;
                  zre0 = zre*re-zim*im;
                  zim = zim*re+zre*im;
                  zre = zre0;
                }
              }
            }
            __syncthreads();

            //make the summation by letting different threads take different points in the coeff array instead
            INITUPDATECOEFF
            __syncthreads();
          }
        }
        if(threadIdx.x==0) {
          innercount=p-j;
        }
        
        //one additional pass for the rest of the values. Almost same as above, but the matrix will not be full
        //This is for performance issues that the full loops are handled separately. This is basically the same
        //code as above, but with additional if-statements
        __syncthreads();
        if(multiple) {
          if(i<count) {
            if(pot==0) {
              //The same algorithm as for the CPU case, but temporary save the results in coeffmatrix instead
              if(j==0) {
                if(threadIdx.x<blockDim.x/2) {
                  coeffmatrix[threadIdx.x]=zre;
                  coeffmatrix[initoffset+threadIdx.x]=zim;
                }
                else if(innercount>=1){
                  coeffmatrix[threadIdx.x]=-zre;
                  coeffmatrix[initoffset+threadIdx.x]=-zim;
                }
                zre0 = zre*re-zim*im;
                zim = zim*re+zre*im;
                zre = zre0;
                for(k=2;k<innercount-1;k+=2) {
                  coeffmatrix[k*initoffset+threadIdx.x]=-zre/(j+k+(threadIdx.x>=blockDim.x/2));
                  coeffmatrix[(k+1)*initoffset+threadIdx.x]=-zim/(j+k+(threadIdx.x>=blockDim.x/2));
                  zre0 = zre*re-zim*im;
                  zim = zim*re+zre*im;
                  zre = zre0;
                }
//                 if(innercount>=4) {
                  coeffmatrix[k*initoffset+threadIdx.x]=-zre/(j+k+(threadIdx.x>=blockDim.x/2));
                  coeffmatrix[(k+1)*initoffset+threadIdx.x]=-zim/(j+k+(threadIdx.x>=blockDim.x/2));
//                 }
              }
              else {
                for(k=0;k<innercount-1;k+=2) {
                  coeffmatrix[k*initoffset+threadIdx.x]=-zre/(j+k+(threadIdx.x>=blockDim.x/2));
                  coeffmatrix[(k+1)*initoffset+threadIdx.x]=-zim/(j+k+(threadIdx.x>=blockDim.x/2));
                  zre0 = zre*re-zim*im;
                  zim = zim*re+zre*im;
                  zre = zre0;
                }
                coeffmatrix[k*initoffset+threadIdx.x]=-zre/(j+k+(threadIdx.x>=blockDim.x/2));
                coeffmatrix[(k+1)*initoffset+threadIdx.x]=-zim/(j+k+(threadIdx.x>=blockDim.x/2));
              }
            }
            else {
  //                     pass 1, calculate coeff matrix
              for(k=0;k<innercount-1;k+=2) {
                coeffmatrix[k*initoffset+threadIdx.x]=zre;
                coeffmatrix[(k+1)*initoffset+threadIdx.x]=zim;
                zre0 = zre*re-zim*im;
                zim = zim*re+zre*im;
                zre = zre0;
              }
              coeffmatrix[k*initoffset+threadIdx.x]=zre;
              coeffmatrix[(k+1)*initoffset+threadIdx.x]=zim;
            }
          }
          __syncthreads();
          //make the summation by letting different threads take different points in the coeff array instead
          INITUPDATECOEFFWITHCHECK2
          __syncthreads();
        }
        else {
          if(i<count) {
            if(pot==0) {
              //The same algorithm as for the CPU case, but temporary save the results in coeffmatrix instead
              if(j==0) {
                coeffmatrix[threadIdx.x]=zre;
                coeffmatrix[initoffset+threadIdx.x]=zim;
                zre0 = zre*re-zim*im;
                zim = zim*re+zre*im;
                zre = zre0;
                if(innercount>=1) {
                  coeffmatrix[2*initoffset+threadIdx.x]=-zre;
                  coeffmatrix[(1*2+1)*initoffset+threadIdx.x]=-zim;
                  zre0 = zre*re-zim*im;
                  zim = zim*re+zre*im;
                  zre = zre0;
                }
                for(k=2;k<innercount;k++) {
                  coeffmatrix[k*2*initoffset+threadIdx.x]=-zre/(j+k);
                  coeffmatrix[(k*2+1)*initoffset+threadIdx.x]=-zim/(j+k);
                  zre0 = zre*re-zim*im;
                  zim = zim*re+zre*im;
                  zre = zre0;
                }
                if(innercount>=2) {
                  coeffmatrix[k*2*initoffset+threadIdx.x]=-zre/(j+k);
                  coeffmatrix[(k*2+1)*initoffset+threadIdx.x]=-zim/(j+k);
                }
              }
              else {
                for(k=0;k<innercount;k++) {
                  coeffmatrix[k*2*initoffset+threadIdx.x]=-zre/(j+k);
                  coeffmatrix[(k*2+1)*initoffset+threadIdx.x]=-zim/(j+k);
                  zre0 = zre*re-zim*im;
                  zim = zim*re+zre*im;
                  zre = zre0;
                }
                coeffmatrix[k*2*initoffset+threadIdx.x]=-zre/(j+k);
                coeffmatrix[(k*2+1)*initoffset+threadIdx.x]=-zim/(j+k);
              }
            }
            else {
  //                     pass 1, calculate coeff matrix
              for(k=0;k<innercount;k++) {
                coeffmatrix[k*2*initoffset+threadIdx.x]=zre;
                coeffmatrix[(k*2+1)*initoffset+threadIdx.x]=zim;
                zre0 = zre*re-zim*im;
                zim = zim*re+zre*im;
                zre = zre0;
              }
              coeffmatrix[k*2*initoffset+threadIdx.x]=zre;
              coeffmatrix[(k*2+1)*initoffset+threadIdx.x]=zim;
            }
          }
          __syncthreads();
          //make the summation by letting different threads take different points in the coeff array instead
          INITUPDATECOEFFWITHCHECK
          __syncthreads();
        }
      }
      for(k=threadIdx.x;k<(p+1)*2;k+=blockDim.x) {
        coeff1[base+k]=coeff[k];
      }
      __syncthreads();
      for(j=threadIdx.x;j<(p+1)*2;j+=blockDim.x) {
        coeff[j] =0;
      }
      
      //Now, take care of the distant interactions in the same way. Works according to the same principle as above, but with a different formula for the coefficients
      for(;m<m1;) {
        __syncthreads(); //if the continue statement is true
        if(threadIdx.x==0) { //update the loop parameters for the new box
          distbox=ir[m];
          if(distbox>0) {
            count=tex1Dfetch(tixptr, distbox+1); /*right now, it is possible that two boxes has same evalbegin if one box has 0 elements*/
            loopmax=count-(count+blockDim.x-1-(tex1Dfetch(tixptr, distbox)%blockDim.x))%blockDim.x-1+blockDim.x;
          }
        }
        __syncthreads();
        if(distbox<0) {
          if(threadIdx.x==0)
            m++;
          __syncthreads();
          continue;
        }
        for(i=threadIdx.x+tex1Dfetch(tixptr, distbox);i<loopmax;i+=blockDim.x) { /*if more points than threads*/
          for(k=0;k<8;k++)
            coeffmatrix[initoffset*k+threadIdx.x]=0; //init entire matrix to 0 to allow for additions with extra elements
          if(count-i+threadIdx.x<blockDim.x/2) {
            if(threadIdx.x>=blockDim.x/2)
              i-=blockDim.x/2;
            else if(i<count) {
              TEXTUREFETCH(zre, tzr, i);
              TEXTUREFETCH(zim, tzi, i);
              zre-=creal(z0);
              zim-=cimag(z0);
              zre0=zre*zre+zim*zim;
              if(zre0==0.0) {
                if(zre!=0.0) {
                  re=1/(zre+(zim/zre)*zim);
                  im=-(re/zre)*zim;
                }
                else {
                  re=0;
                  im=-1/zim;
                }
              }
              else {
                re=zre/zre0;
                im=-zim/zre0;
              }
              if(complexpoint) {
                TEXTUREFETCH(zre, tmr, i);
                TEXTUREFETCH(zim, tmi, i);
                if(pot==0) {
                  zre0=0.5*log(re*re+im*im);
                  coeffmatrix[threadIdx.x]=zre*zre0;
                  coeffmatrix[initoffset+threadIdx.x]=zim*zre0;
                  zre0=atan2(-im, -re);
                  coeffmatrix[threadIdx.x]-=zim*zre0;
                  coeffmatrix[initoffset+threadIdx.x]+=zre*zre0;
                  
                  zre0=zre*re-zim*im;
                  zim=zre*im+zim*re;
                  zre=zre0;
                  coeffmatrix[threadIdx.x+blockDim.x/2]=0.5*zre;
                  coeffmatrix[initoffset+threadIdx.x+blockDim.x/2]=0.5*zim;
                }
                zre0=zre*re-zim*im;
                zim=zre*im+zim*re;
                zre=zre0;
                
              }
              else {
                TEXTUREFETCH(zre, tmr, i);
                if(pot==0) {
                  zre0=0.5*log(re*re+im*im);
                  zim=atan2(-im, -re);
                  coeffmatrix[threadIdx.x]=zre*zre0;
                  coeffmatrix[initoffset+threadIdx.x]=zim*zre0;
                  
                  zim=zre*im;
                  zre*=re;
                  coeffmatrix[threadIdx.x+blockDim.x/2]=zre;
                  coeffmatrix[initoffset+threadIdx.x+blockDim.x/2]=zim;
                  zre0=zre*re-zim*im;
                  zim=zre*im+zim*re;
                  zre=zre0;
                }
                else {
                  zim=zre*im;
                  zre*=re;
                }
              }
              coeffmatrix[initoffset*2+threadIdx.x]=re;
              coeffmatrix[initoffset*3+threadIdx.x]=im;
              coeffmatrix[initoffset*4+threadIdx.x]=zre;
              coeffmatrix[initoffset*5+threadIdx.x]=zim;
            }
            __syncthreads();
            //due to the complex operations in this case, only let half of the threads perform them
            //and copy the results.
            if(i<count&&threadIdx.x>=blockDim.x/2){
              re=coeffmatrix[initoffset*2+threadIdx.x-blockDim.x/2];
              im=coeffmatrix[initoffset*3+threadIdx.x-blockDim.x/2];
              zre=coeffmatrix[initoffset*4+threadIdx.x-blockDim.x/2];
              zim=coeffmatrix[initoffset*5+threadIdx.x-blockDim.x/2];
              zre0=zre*re-zim*im;
              zim=zre*im+zim*re;
              zre=zre0;
            }
            zre0=re*re-im*im;
            im=2*re*im;
            re=zre0;
            __syncthreads();
            for(j=0;j<=p-initnrcoeff*2;j+=initnrcoeff*2) { //take initnrcoeff coefficients each loop. High values gives less idle threads, but used more memory (which slows down the code)
              if(i<count) {
                //same algorithm as for CPU
                if(pot==0) {
                  for(k=(j==0);k<initnrcoeff;k++) {
                    coeffmatrix[k*2*initoffset+threadIdx.x]=zre/(j+k*2+(threadIdx.x>=blockDim.x/2));
                    coeffmatrix[(k*2+1)*initoffset+threadIdx.x]=zim/(j+k*2+(threadIdx.x>=blockDim.x/2));
                    zre0 = zre*re-zim*im;
                    zim = zim*re+zre*im;
                    zre = zre0;
                  }
                }
                else {

                  //pass 1, calculate coeff matrix
                  for(k=0;k<initnrcoeff;k++) {
                    coeffmatrix[k*2*initoffset+threadIdx.x]=zre;
                    coeffmatrix[(k*2+1)*initoffset+threadIdx.x]=zim;
                    zre0 = zre*re-zim*im;
                    zim = zim*re+zre*im;
                    zre = zre0;

                  }
                }
              }
              __syncthreads();

              //The same summation as for the local coefficients
              INITUPDATECOEFF2
              __syncthreads();
            }
            if(threadIdx.x==0) {
              innercount=p-j;
            }

            //same as before, but without a full matrix. See above for comments on the algorithm
            __syncthreads();
            if(i<count) {
              if(pot==0) {
  //                     pass 1, calculate coeff matrix
                for(k=(j==0)*2;k<innercount-1;k+=2) {
                  coeffmatrix[k*initoffset+threadIdx.x]=zre/(j+k+(threadIdx.x>=blockDim.x/2));
                  coeffmatrix[(k+1)*initoffset+threadIdx.x]=zim/(j+k+(threadIdx.x>=blockDim.x/2));
                  zre0 = zre*re-zim*im;
                  zim = zim*re+zre*im;
                  zre = zre0;

                }
                coeffmatrix[k*initoffset+threadIdx.x]=zre/(j+k+(threadIdx.x>=blockDim.x/2));
                coeffmatrix[(k+1)*initoffset+threadIdx.x]=zim/(j+k+(threadIdx.x>=blockDim.x/2));
              }
              else {
                for(k=0;k<innercount-1;k+=2) {
                  coeffmatrix[k*initoffset+threadIdx.x]=zre;
                  coeffmatrix[(k+1)*initoffset+threadIdx.x]=zim;
                  zre0 = zre*re-zim*im;
                  zim = zim*re+zre*im;
                  zre = zre0;

                }
                coeffmatrix[k*initoffset+threadIdx.x]=zre;
                coeffmatrix[(k+1)*initoffset+threadIdx.x]=zim;
              }
            }
            __syncthreads();
            //The same summation as for the local coefficients
            INITUPDATECOEFFWITHCHECK2
            __syncthreads(); 
          }
          else {
            if(i<count) { //if point belongs to the box, read input, check CPU version for how the algorithm works, it is the same, but the temporary coeffmatrix is used here as well
              TEXTUREFETCH(zre, tzr, i);
              TEXTUREFETCH(zim, tzi, i);
              zre-=creal(z0);
              zim-=cimag(z0);
              zre0=zre*zre+zim*zim;
              if(zre0==0.0) {
                if(zre!=0.0) {
                  re=1/(zre+(zim/zre)*zim);
                  im=-(re/zre)*zim;
                }
                else {
                  re=0;
                  im=-1/zim;
                }
              }
              else {
                re=zre/zre0;
                im=-zim/zre0;
              }
              if(complexpoint) {
                TEXTUREFETCH(zre, tmr, i);
                TEXTUREFETCH(zim, tmi, i);
                if(pot==0) {
                  zre0=0.5*log(re*re+im*im);
                  coeffmatrix[threadIdx.x]=zre*zre0;
                  coeffmatrix[initoffset+threadIdx.x]=zim*zre0;
                  zre0=atan2(-im, -re);
                  coeffmatrix[threadIdx.x]-=zim*zre0;
                  coeffmatrix[initoffset+threadIdx.x]+=zre*zre0;
                }
                zre0=zre*re-zim*im;
                zim=zre*im+zim*re;
                zre=zre0;
              }
              else {
                TEXTUREFETCH(zre, tmr, i);
                if(pot==0) {
                  zre0=0.5*log(re*re+im*im);
                  zim=atan2(-im, -re);
                  coeffmatrix[threadIdx.x]=zre*zre0;
                  coeffmatrix[initoffset+threadIdx.x]=zim*zre0;
                }
                zim=zre*im;
                zre*=re;
              }

            }

            for(j=0;j<=p-initnrcoeff;j+=initnrcoeff) { //take initnrcoeff coefficients each loop. High values gives less idle threads, but used more memory (which slows down the code)
              if(i<count) {
                //same algorithm as for CPU
                if(pot==0) {
                  for(k=(j==0);k<initnrcoeff;k++) {
                    coeffmatrix[k*2*initoffset+threadIdx.x]=zre/(j+k);
                    coeffmatrix[(k*2+1)*initoffset+threadIdx.x]=zim/(j+k);
                    zre0 = zre*re-zim*im;
                    zim = zim*re+zre*im;
                    zre = zre0;
                  }
                }
                else {

                  //pass 1, calculate coeff matrix
                  for(k=0;k<initnrcoeff;k++) {
                    coeffmatrix[k*2*initoffset+threadIdx.x]=zre;
                    coeffmatrix[(k*2+1)*initoffset+threadIdx.x]=zim;
                    zre0 = zre*re-zim*im;
                    zim = zim*re+zre*im;
                    zre = zre0;

                  }
                }
              }
              __syncthreads();

              //The same summation as for the local coefficients
              INITUPDATECOEFF
              __syncthreads();
            }
            if(threadIdx.x==0) {
              innercount=p-j;
            }

            //same as before, but without a full matrix. See above for comments on the algorithm
            __syncthreads();
            if(i<count) {
              if(pot==0) {
  //                     pass 1, calculate coeff matrix
                for(k=(j==0);k<innercount;k++) {
                  coeffmatrix[k*2*initoffset+threadIdx.x]=zre/(j+k);
                  coeffmatrix[(k*2+1)*initoffset+threadIdx.x]=zim/(j+k);
                  zre0 = zre*re-zim*im;
                  zim = zim*re+zre*im;
                  zre = zre0;

                }
                coeffmatrix[k*2*initoffset+threadIdx.x]=zre/(j+k);
                coeffmatrix[(k*2+1)*initoffset+threadIdx.x]=zim/(j+k);
              }
              else {
                for(k=0;k<innercount;k++) {
                  coeffmatrix[k*2*initoffset+threadIdx.x]=zre;
                  coeffmatrix[(k*2+1)*initoffset+threadIdx.x]=zim;
                  zre0 = zre*re-zim*im;
                  zim = zim*re+zre*im;
                  zre = zre0;

                }
                coeffmatrix[k*2*initoffset+threadIdx.x]=zre;
                coeffmatrix[(k*2+1)*initoffset+threadIdx.x]=zim;
              }
            }
            __syncthreads();
            //The same summation as for the local coefficients
            INITUPDATECOEFFWITHCHECK
            __syncthreads();
          }
        }
        if(threadIdx.x==0)
          m++;
        __syncthreads();
      }
      for(k=threadIdx.x;k<(p+1)*sizeof(dcmplx)/sizeof(double);k+=blockDim.x) {
        coeff2[base+k]=-coeff[k];
      }
        __syncthreads();
    }
  }
}