__shared__ dcmplx rpow[M2PSMAXTHREADS];
#ifdef LOGM2PS
__shared__ dcmplx Thisthat0[M2PSMAXTHREADS/2];
#endif
double re,im,pre,pim,resqr,imsqr;
int shifttmp;
int j,k;
__shared__ int currentlevel,offset,m,mend,loopend,boxshift,loopcount,oldboxnr,boxnr;
__shared__ int *jcptr,*kcptr,*ir;
__shared__ double* result;
__shared__ SORT_DCMPLX zlocal[2];
double a[UNROLLLEVEL+1];
if(threadIdx.x==0){
  result=wksp+(pshift+1)*M2PSMAXTHREADS;
  currentlevel=startlevel;
  offset=0;
  for(j=0;j<currentlevel;j++)
    offset+=(1<<(j<<1));
  boxnr=blockIdx.x;
  while(boxnr>=(1<<(currentlevel<<1))) {
    boxnr-=(1<<(currentlevel<<1));
    offset+=(1<<(currentlevel<<1));
    currentlevel++;
  }
  oldboxnr=boxnr+offset;
  if(currentlevel<=nlevel) {
    jcptr=jcptr2[currentlevel];
    kcptr=kcptr2[currentlevel];
    ir=ir2[currentlevel];
  }
  m=1;
  mend=0;
}
__syncthreads();
for(k=threadIdx.x;k<=pshift*2+1;k+=blockDim.x)
  result[k]=0;
__syncthreads();
//choose working box
for(;;) {
  if(currentlevel>nlevel)
    break;
  if(m>=mend) { //operate on a new box?
    if(threadIdx.x==0) {
      m=kcptr[boxnr];
      mend=jcptr[boxnr+1];
      loopend=imin(mend, m+blockDim.x/2);
      boxshift=loopend;
      loopcount=loopend-m;
      boxshift=loopcount;
    }
    __syncthreads();
    //load box position
    if(threadIdx.x<sizeof(SORT_DCMPLX)/sizeof(float)) {
      ((float*)&zlocal)[threadIdx.x] =tex1Dfetch(tz0,(boxnr+offset)*sizeof(SORT_DCMPLX)/sizeof(float)+threadIdx.x);
    }
  }
  
  
  if((threadIdx.x>>1)<loopcount) {
    k=ir[m+(threadIdx.x/2)]+offset;
    wksp[threadIdx.x]=((SORT_REAL*)&z0[k])[threadIdx.x&1];
#ifdef LOGM2PS
    ((double*)Thisthat0)[threadIdx.x]=((double*)&coeff1[k*(pshift+1)])[threadIdx.x&1];
#endif
  }
  __syncthreads();
  if(mend-m<(blockDim.x>>1)) { //operate on two boxes if space is available
    if(threadIdx.x==0) {
      boxnr+=gridDim.x;
      while(boxnr>=(1<<(currentlevel<<1))) {
        boxnr-=(1<<(currentlevel<<1));
        offset+=(1<<(currentlevel<<1));
        currentlevel++;
      }
      if(currentlevel<=nlevel) {
        jcptr=jcptr2[currentlevel];
        kcptr=kcptr2[currentlevel];
        ir=ir2[currentlevel];
        m=kcptr[boxnr];
        mend=jcptr[boxnr+1];
        loopend=imin(mend, m+blockDim.x/2-boxshift);
        loopcount+=loopend-m;
      }
    }
    __syncthreads();
    if(currentlevel<=nlevel) {
      if(threadIdx.x<2)
        ((SORT_REAL*)&zlocal[1])[threadIdx.x]=((SORT_REAL*)&z0[(boxnr+offset)])[threadIdx.x];
      if((threadIdx.x>>1)>=boxshift&&(threadIdx.x>>1)<loopcount) {
        k=ir[m+(threadIdx.x/2)-boxshift]+offset;
        wksp[threadIdx.x]=((SORT_REAL*)&z0[k])[threadIdx.x&1];
#ifdef LOGM2PS
        ((double*)Thisthat0)[threadIdx.x]=((double*)&coeff1[k*(pshift+1)])[threadIdx.x&1];
#endif
      }
    }
  }
  __syncthreads();
  if(threadIdx.x<loopcount) {
    pre=creal(zlocal[threadIdx.x>=boxshift])-wksp[threadIdx.x*2];
    pim=cimag(zlocal[threadIdx.x>=boxshift])-wksp[threadIdx.x*2+1];
    //correct division
    //code is the same as the commented code below, but for cuda, it is better to keep the warp together during the division
    int testres=fabs(pre)<fabs(pim);
    if(testres) {
      double retmp=pre;
      pre=pim;
      pim=retmp;
    }
    im=pim/pre;
    re=pre+im*pim;
    double tmp2=1.0/re;
    double tmp3=im*tmp2;
    if(testres) {
      COMPLEXASSIGN(rpow[threadIdx.x*2], tmp3, -tmp2);
    }
    else {
       COMPLEXASSIGN(rpow[threadIdx.x*2], tmp2, -tmp3);
    }
//     if(fabs(pre)>fabs(pim)) {
//       im=pim/pre;
//       re=pre+im*pim;
//       COMPLEXASSIGN(rpow[threadIdx.x*2], 1/re, -im*(1/re));
//     }
//     else {
//       im=pre/pim;
//       re=pim+im*pre;
//       COMPLEXASSIGN(rpow[threadIdx.x*2], im*(1/re), -1/re);
//     }
    //calculate sqr
    COMPLEXASSIGN(rpow[threadIdx.x*2+1], creal(rpow[threadIdx.x*2])*creal(rpow[threadIdx.x*2])-cimag(rpow[threadIdx.x*2])*cimag(rpow[threadIdx.x*2]), 2*creal(rpow[threadIdx.x*2])*cimag(rpow[threadIdx.x*2]));
  }
  __syncthreads();
  wksp[threadIdx.x]=0;
  __syncthreads();

  if((threadIdx.x>>1)<loopcount&&(threadIdx.x&1)<=pshift-1) {
    dcmplx pz;
    re=creal(rpow[threadIdx.x]); //scaling value
    im=cimag(rpow[threadIdx.x]);
    fetch_double2((double2*)&pz,tcoeff1d,k*(pshift+1)+(1+(threadIdx.x&1)));
    shifttmp=M2PSMAXTHREADS*((threadIdx.x&1)+1)+(threadIdx.x&(~1));
    wksp[shifttmp]=creal(pz)*re-cimag(pz)*im; //scale the coeffient
    wksp[shifttmp+1] = creal(pz)*im+cimag(pz)*re;
    resqr=creal(rpow[(threadIdx.x|1)]);
    imsqr=cimag(rpow[(threadIdx.x|1)]);
  }
  __syncthreads(); //not necessary, but on cuda 3.2 with gtx480, the code is actually slightly (1-2%) faster with this syncthreads enabled

  //loop is manually unrolled
  if((threadIdx.x>>1)<loopcount) {
    for(j=3;j<=pshift-7;j+=8) {
      dcmplx pz;
      shifttmp=M2PSMAXTHREADS*(j+(threadIdx.x&1))+(threadIdx.x&(~1));
      pre=re*resqr-im*imsqr; //calculate next power of r
      im=re*imsqr+im*resqr;
      re=pre;
      fetch_double2((double2*)&pz,tcoeff1d,k*(pshift+1)+(j+(threadIdx.x&1)));
      wksp[shifttmp]=creal(pz)*re-cimag(pz)*im; //scale the coeffient
      wksp[shifttmp+1] = creal(pz)*im+cimag(pz)*re;
      
      shifttmp+=2*M2PSMAXTHREADS;
      pre=re*resqr-im*imsqr; //calculate next power of r
      im=re*imsqr+im*resqr;
      re=pre;
      fetch_double2((double2*)&pz,tcoeff1d,k*(pshift+1)+(j+2+(threadIdx.x&1)));
      wksp[shifttmp]=creal(pz)*re-cimag(pz)*im; //scale the coeffient
      wksp[shifttmp+1] = creal(pz)*im+cimag(pz)*re;
      
      shifttmp+=2*M2PSMAXTHREADS;
      pre=re*resqr-im*imsqr; //calculate next power of r
      im=re*imsqr+im*resqr;
      re=pre;
      fetch_double2((double2*)&pz,tcoeff1d,k*(pshift+1)+(j+4+(threadIdx.x&1)));
      wksp[shifttmp]=creal(pz)*re-cimag(pz)*im; //scale the coeffient
      wksp[shifttmp+1] = creal(pz)*im+cimag(pz)*re;
      
      shifttmp+=2*M2PSMAXTHREADS;
      pre=re*resqr-im*imsqr; //calculate next power of r
      im=re*imsqr+im*resqr;
      re=pre;
      fetch_double2((double2*)&pz,tcoeff1d,k*(pshift+1)+(j+6+(threadIdx.x&1)));
      wksp[shifttmp]=creal(pz)*re-cimag(pz)*im; //scale the coeffient
      wksp[shifttmp+1] = creal(pz)*im+cimag(pz)*re;
    }
    //second manual unroll
    for(/*j=3*/;j<=pshift-3;j+=4) {
      dcmplx pz;
      shifttmp=M2PSMAXTHREADS*(j+(threadIdx.x&1))+(threadIdx.x&(~1));
      pre=re*resqr-im*imsqr; //calculate next power of r
      im=re*imsqr+im*resqr;
      re=pre;
      fetch_double2((double2*)&pz,tcoeff1d,k*(pshift+1)+(j+(threadIdx.x&1)));
      wksp[shifttmp]=creal(pz)*re-cimag(pz)*im; //scale the coeffient
      wksp[shifttmp+1] = creal(pz)*im+cimag(pz)*re;
      
      shifttmp+=2*M2PSMAXTHREADS;
      pre=re*resqr-im*imsqr; //calculate next power of r
      im=re*imsqr+im*resqr;
      re=pre;
      fetch_double2((double2*)&pz,tcoeff1d,k*(pshift+1)+(j+2+(threadIdx.x&1)));
      wksp[shifttmp]=creal(pz)*re-cimag(pz)*im; //scale the coeffient
      wksp[shifttmp+1] = creal(pz)*im+cimag(pz)*re;
    }
    for(;j<=pshift;j+=2) {
      dcmplx pz;
      shifttmp=M2PSMAXTHREADS*(j+(threadIdx.x&1))+(threadIdx.x&(~1));
      if((j+(threadIdx.x&1))<=pshift) {
        pre=re*resqr-im*imsqr; //calculate next power of r
        im=re*imsqr+im*resqr;
        re=pre;
        fetch_double2((double2*)&pz,tcoeff1d,k*(pshift+1)+(j+(threadIdx.x&1)));
        wksp[shifttmp]=creal(pz)*re-cimag(pz)*im; //scale the coeffient
        wksp[shifttmp+1] = creal(pz)*im+cimag(pz)*re;
      }
    }
//     __syncthreads(); //should not be necessary
  }
  __syncthreads();
  
  
  //perform horner scheme interaction
      // (1) "upward" phase
  if((threadIdx.x>>1)<loopcount) {
    k=1;
    #if UNROLLLEVEL>=7 //loop is performed on different levels depending on UNROLLLEVEL
    for (; k <= pshift-6; k+=7) {
      INNERUNROLLBASE8(wksp,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k) //here, imaginary part is taken directly after real part for each operation
    }
    if (k <= pshift-5)
    #elif UNROLLLEVEL>=6 //The pattern repeats itself for each unroll level
    for (;k <= pshift-5;)
    #endif
    #if UNROLLLEVEL>=6
    {
      INNERUNROLLBASE7(wksp,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
      k+=6;
    }
    if (k <= pshift-4)
    #elif UNROLLLEVEL>=5
    for (;k <= pshift-4;)
    #endif
    #if UNROLLLEVEL>=5
    {
      INNERUNROLLBASE6(wksp,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
      k+=5;
    }
    if (k <= pshift-3)
    #elif UNROLLLEVEL>=4
    for (;k <= pshift-3;)
    #endif
    #if UNROLLLEVEL>=4
    {
      INNERUNROLLBASE5(wksp,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
      k+=4;
    }
    if (k <= pshift-2)
    #elif UNROLLLEVEL>=3
    for (;k <= pshift-2;)
    #endif
    #if UNROLLLEVEL>=3
    {
      INNERUNROLLBASE4(wksp,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
      k+=3;
    }
    if (k <= pshift-1)
    #elif UNROLLLEVEL>=2
    for (;k <= pshift-1;)
    #endif
    #if UNROLLLEVEL>=2
    {
      INNERUNROLLBASE3(wksp,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
      k+=2;
    }
    if (k <= pshift)
    #else
    for (;k <= pshift;k++)
    #endif
    {
      INNERUNROLLBASE2(wksp,0,STOREMUP,LOADM1,ADDEUP,0,pshift-k)
    }
    
    //downward
    k=pshift;
    #if UNROLLLEVEL>=7
    for (; k>7; k-=7) {
      INNERUNROLLBASE8(wksp,0,STOREMDN,LOADM1,ADDEDN,1,k)
    }
    if (k>6)
    #elif UNROLLLEVEL>=6
    for (;k>6;)
    #endif
    #if UNROLLLEVEL>=6
    {
      INNERUNROLLBASE7(wksp,0,STOREMDN,LOADM1,ADDEDN,1,k)
      k-=6;
    }
    if (k>5)
    #elif UNROLLLEVEL>=5
    for (;k>5;)
    #endif
    #if UNROLLLEVEL>=5
    {
      INNERUNROLLBASE6(wksp,0,STOREMDN,LOADM1,ADDEDN,1,k)
      k-=5;
    }
    if (k>4)
    #elif UNROLLLEVEL>=4
    for (;k>4;)
    #endif
    #if UNROLLLEVEL>=4
    {
      INNERUNROLLBASE5(wksp,0,STOREMDN,LOADM1,ADDEDN,1,k)
      k-=4;
    }
    if (k>3)
    #elif UNROLLLEVEL>=3
    for (;k>3;)
    #endif
    #if UNROLLLEVEL>=3
    {
      INNERUNROLLBASE4(wksp,0,STOREMDN,LOADM1,ADDEDN,1,k)
      k-=3;
    }
    if (k>2)
    #elif UNROLLLEVEL>=2
    for (;k>2;)
    #endif
    #if UNROLLLEVEL>=2
    {
      INNERUNROLLBASE3(wksp,0,STOREMDN,LOADM1,ADDEDN,1,k)
      k-=2;
    }
    if (k > 1)
    #else
    for (;k > 1;k--)
    #endif
    {
      INNERUNROLLBASE2(wksp,0,STOREMDN,LOADM1,ADDEDN,1,k)
    }
    //this code would give the same results s the code above, but not as fast
//     for (k = 1; k <= pshift; k++)
// //       for (j = pshift-k; j < pshift; j++) {
// //         wksp[j*M2PSMAXTHREADS+threadIdx.x] += wksp[(j+1)*M2PSMAXTHREADS+threadIdx.x];
// //       }
    // (2) "downward" phase
//     for (k = pshift; k > 1; k--)
// //       for (j = k; j <= pshift; j++) {
// //         wksp[j*M2PSMAXTHREADS+threadIdx.x] += wksp[(j-1)*M2PSMAXTHREADS+threadIdx.x];
// //       }

  }
  __syncthreads();
  //perform post-scaling
  if((threadIdx.x>>1)<loopcount) {
    //the first part already has r values calculated, take these as a separate step
    shifttmp=M2PSMAXTHREADS*(threadIdx.x&1)+(threadIdx.x&(~1));
#ifdef LOGM2PS
    if((threadIdx.x&1)==0) { //coefficient 0 has additional terms for logarithmic interaction
      re=0.5*log(creal(rpow[(threadIdx.x&(~1))])*creal(rpow[(threadIdx.x&(~1))])+cimag(rpow[(threadIdx.x&(~1))])*cimag(rpow[(threadIdx.x&(~1))]));
      im=atan2(creal(rpow[(threadIdx.x&(~1))]), cimag(rpow[(threadIdx.x&(~1))]));
      wksp[shifttmp]-=re*creal(Thisthat0[(threadIdx.x>>1)])+cimag(Thisthat0[(threadIdx.x>>1)])*im;
      wksp[shifttmp+1]+=creal(Thisthat0[(threadIdx.x>>2)*2])*im+cimag(Thisthat0[(threadIdx.x>>2)*2])*re;
      re=-1;
      im=0;
    }
    else if((threadIdx.x&1)<=pshift){ //scale the other three elements
      re=creal(rpow[threadIdx.x-1]); //set scaling term starting value
      im=cimag(rpow[threadIdx.x-1]);
      wksp[shifttmp]-=creal(Thisthat0[(threadIdx.x>>1)]); //apply the additional term from the logarithm
      wksp[shifttmp+1]-=cimag(Thisthat0[(threadIdx.x>>1)]);
      //perform the scaling
      pre=-(wksp[shifttmp]*re-wksp[shifttmp+1]*im)/**(1-(signed int)((threadIdx.x&1)<<1))*/;
      wksp[shifttmp+1]=-(wksp[shifttmp]*im+wksp[shifttmp+1]*re)/**(1-(signed int)((threadIdx.x&1)<<1))*/;
      wksp[shifttmp]=pre;
    }
#else
    if((threadIdx.x&1)==0) { //in this phase, the r vector should start from 0, while the saved one starts from 1, add the first unity element
      re=-1; //set scaling term starting value
      im=0;
    }
    else {
      re=creal(rpow[threadIdx.x-1]); //set scaling term starting value
      im=cimag(rpow[threadIdx.x-1]);
    }
    if((threadIdx.x&1)<=pshift){ //perform the scaling
      pre=-(wksp[shifttmp]*re-wksp[shifttmp+1]*im)/**(1-(signed int)((threadIdx.x&1)<<1))*/;
      wksp[shifttmp+1]=-(wksp[shifttmp]*im+wksp[shifttmp+1]*re)/**(1-(signed int)((threadIdx.x&1)<<1))*/;
      wksp[shifttmp]=pre;
    }
#endif
    for(j=(threadIdx.x&1)+2;j<=pshift;j+=2) { //scale the rest of the elements
      pre=re*resqr-im*imsqr; //calculate the next term in the scaling by multiplying with 1/r^4
      im=re*imsqr+im*resqr;
      re=pre;
      shifttmp=M2PSMAXTHREADS*j+(threadIdx.x&(~1));

#ifdef LOGM2PS
      wksp[shifttmp]-=creal(Thisthat0[(threadIdx.x>>1)])/j; //apply the additional term from the logarithm
      wksp[shifttmp+1]-=cimag(Thisthat0[(threadIdx.x>>1)])/j;
#endif
      //perform the scaling
      pre=-(wksp[shifttmp]*re-wksp[shifttmp+1]*im)/**(1-(signed int)((j&1)<<1))*/;
      wksp[shifttmp+1]=-(wksp[shifttmp]*im+wksp[shifttmp+1]*re)/**(1-(signed int)((j&1)<<1))*/;
      wksp[shifttmp]=pre;
    }
  }
  __syncthreads();
  //sum the results
  /*note: the complicated scheme here is to avoid bank conflicts, let all summations start with an offset of 1*/
#if M2PSMAXTHREADS==8 || M2PSMAXTHREADS==16 || M2PSMAXTHREADS==32 || M2PSMAXTHREADS==64 || M2PSMAXTHREADS==128 || M2PSMAXTHREADS==256
  if(boxshift==M2PSMAXTHREADS/2) { //multiple of 32/64 etc, modulo operator % can be changed to &
    for(k=threadIdx.x;k<=pshift*2+1;k+=blockDim.x) {
      for(j=0,shifttmp=(threadIdx.x&(M2PSMAXTHREADS/2-1));j<M2PSMAXTHREADS/2;j++,shifttmp++,shifttmp&=(M2PSMAXTHREADS/2-1)) {
        result[k]+=wksp[(k>>1)*M2PSMAXTHREADS+(shifttmp<<1)+(threadIdx.x&1)];
      }
    }
  }
  else {
    for(k=threadIdx.x;k<=pshift*2+1;k+=blockDim.x) {
      for(j=0,shifttmp=threadIdx.x%boxshift;j<boxshift;j++,shifttmp++,shifttmp=shifttmp==boxshift?0:shifttmp){
        result[k]+=wksp[(k>>1)*M2PSMAXTHREADS+(shifttmp<<1)+(threadIdx.x&1)];
      }
    }
  }
#else
  for(k=threadIdx.x;k<=pshift*2+1;k+=blockDim.x) {
    for(j=0,shifttmp=threadIdx.x%boxshift;j<boxshift;j++,shifttmp++,shifttmp=shifttmp==boxshift?0:shifttmp){
      result[k]+=wksp[(k>>1)*M2PSMAXTHREADS+(shifttmp<<1)+(threadIdx.x&1)];
    }
  }
#endif
  __syncthreads();
  //write results to memory
  if(loopend==mend||boxshift<blockDim.x/2) {
    for(k=threadIdx.x;k<=pshift*2+1;k+=blockDim.x) {
      ((double*)&coeff2[oldboxnr*(pshift+1)])[k]+=result[k];
      result[k]=0;
    }
    __syncthreads();
  }
  
  //two boxes involved in split
  if(loopcount!=boxshift) {
    /*note: using the complicated scheme above did actually reduce performane in this case*/
    for(k=threadIdx.x;k<=pshift*2+1;k+=blockDim.x) {
      for(j=boxshift;j<loopcount;j++){
        result[k]+=wksp[(k>>1)*M2PSMAXTHREADS+(j<<1)+(threadIdx.x&1)];
      }
    }
    __syncthreads();
    if(loopend==mend) { //second box finished as well
      for(k=threadIdx.x;k<=pshift*2+1;k+=blockDim.x) {
        ((double*)&coeff2[(boxnr+offset)*(pshift+1)])[k]+=result[k];
        result[k]=0;
      }
    }
    else {
      if(threadIdx.x==0) {
        oldboxnr=boxnr+offset;
        COMPLEXASSIGN(zlocal[0],creal(zlocal[1]),cimag(zlocal[1]));
      }
    }
    __syncthreads();
  }
  
  if(threadIdx.x==0) {
    m=loopend;
    loopend=imin(mend, m+blockDim.x/2);
    loopcount=loopend-m;
    boxshift=loopcount;
  }
  __syncthreads();
  if(m>=mend){
    //update needs to be performed her for the if-statement at the beginning of the for-loop to work
    if(threadIdx.x==0) {
      boxnr+=gridDim.x;
      while(boxnr>=(1<<(currentlevel<<1))) {
        boxnr-=(1<<(currentlevel<<1));
        offset+=(1<<(currentlevel<<1));
        currentlevel++;
      }
      oldboxnr=boxnr+offset;
      if(currentlevel<=nlevel) {
        jcptr=jcptr2[currentlevel];
        kcptr=kcptr2[currentlevel];
        ir=ir2[currentlevel];
      }
    }
    __syncthreads();
  }
}
