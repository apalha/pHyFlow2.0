#if defined(CUDASUPPORT) && defined(CUDASORT)

//starting values for indices, [1,2,3,.....,count]
__global__ void initiateindices(int* indices,int count) 
{
  for(int i=blockDim.x * blockIdx.x + threadIdx.x;i<count;i+=blockDim.x*gridDim.x) {
    indices[i]=i;
  }
}

//Currently, the debug code is being kept here, in case of future modifications of the sorting. Kept outside the main code, to make main code more readable. To be removed later
#define DEBUGSINGLEBLOCKPARTITIONLINE1 debugvector[32+20*lc+splitindex*1000]=lcount+localsplitstart-splitstart;\
debugvector[33+40*lc+splitindex*1000]=splitend-localsplitend+rcount;\
debugvector[35+40*lc+splitindex*1000]=lcount;
                
#define DEBUGSINGLEBLOCKPARTITIONLINE2  if(threadIdx.x==0) {\
  debugvector[29+40*lc+splitindex*1000]=oldsplitstart;\
  debugvector[30+40*lc+splitindex*1000]=oldsplitend;\
}
#define DEBUGSINGLEBLOCKPARTITIONLINE3 if(threadIdx.x==0) {\
  debugvector[20+40*lc+splitindex*1000]=point1;\
  debugvector[21+40*lc+splitindex*1000]=point2;\
  debugvector[22+40*lc+splitindex*1000]=point3;\
  debugvector[23+40*lc+splitindex*1000]=localsplitstart;\
  debugvector[24+40*lc+splitindex*1000]=localsplitend;\
  if(newactive) {\
    debugvector[25+40*lc+splitindex*1000]=activenewpositions[localsplitstart];\
    debugvector[26+40*lc+splitindex*1000]=activenewpositions[(localsplitstart+localsplitend)/2];\
    debugvector[27+40*lc+splitindex*1000]=activenewpositions[localsplitend-1];\
  }\
  else {\
    debugvector[25+40*lc+splitindex*1000]=activepositions[localsplitstart];\
    debugvector[26+40*lc+splitindex*1000]=activepositions[(localsplitstart+localsplitend)/2];\
    debugvector[27+40*lc+splitindex*1000]=activepositions[localsplitend-1];\
  }\
  debugvector[28+40*lc+splitindex*1000]=splitpoint;\
  debugvector[31+40*lc+splitindex*1000]=newactive;\
  debugvector[34+40*lc+splitindex*1000]=count;\
  debugvector[36+40*lc+splitindex*1000]=(localsplitstart+localsplitend)/2;\
  debugvector[37+40*lc+splitindex*1000]=1;\
  debugvector[38+40*lc+splitindex*1000]=threesplit;\
  lc++;\
}
#define DEBUGSINGLEBLOCKPARTITIONLINE4 debugvector[32+20*lc+splitindex*1000]=addpoint;\
debugvector[33+40*lc+splitindex*1000]=index;\
debugvector[38+40*lc+splitindex*1000]=localleftcount[nractivethreads-1];
#define DEBUGSINGLEBLOCKPARTITIONLINE5 if(threadIdx.x==0) {\
  debugvector[20+40*lc+splitindex*1000]=point1;\
  debugvector[21+40*lc+splitindex*1000]=point2;\
  debugvector[22+40*lc+splitindex*1000]=point3;\
  debugvector[23+40*lc+splitindex*1000]=innersplitstart;\
  debugvector[24+40*lc+splitindex*1000]=innersplitend;\
  if(innernewactive) {\
    debugvector[25+40*lc+splitindex*1000]=positioncache[indexcache2[innersplitstart]];\
    debugvector[26+40*lc+splitindex*1000]=positioncache[indexcache2[(innersplitstart+innersplitend)/2]];\
    debugvector[27+40*lc+splitindex*1000]=positioncache[indexcache2[innersplitend-1]];\
  }\
  else {\
    debugvector[25+40*lc+splitindex*1000]=positioncache[indexcache[innersplitstart]];\
    debugvector[26+40*lc+splitindex*1000]=positioncache[indexcache[(innersplitstart+innersplitend)/2]];\
    debugvector[27+40*lc+splitindex*1000]=positioncache[indexcache[innersplitend-1]];\
  }\
  debugvector[28+40*lc+splitindex*1000]=splitpoint;\
  debugvector[31+40*lc+splitindex*1000]=innernewactive;\
  debugvector[34+40*lc+splitindex*1000]=count;\
  debugvector[36+40*lc+splitindex*1000]=(innersplitstart+innersplitend)/2;\
  debugvector[37+40*lc+splitindex*1000]=1;\
  lc++;\
}

//This function performs partitioning using only one block, which allows it to
//complete in one single call. The function consists of two splitting method
//The first one works in all cases, and reglularly reads from memory. The second
//caches the entire array, and only works from shared memory.
#define MEDIANCOUNT 32
extern __shared__ int indexcache[];
__global__ void singleblockpartition(SORT_REAL* positions,SORT_REAL *positions2,int* indices,int* newindices,SORT_REAL *splitpoints,int *llimits,int* rlimits,int *newllimits,int* newrlimits,int* basellimits,int* baserlimits,int splitcount,SORT_REAL* newpositions,SORT_REAL *newpositions2,int* ysplit,int resultsinoriginal SORTLIMITSTRING DEBUGVECTORSTRING) 
//positions is input positions for x values
//positions2 is input positions for y values
//indices will be the permutation matrix
//newindices will contain the output indices (unless resultsinoriginal is set)
//splitpoints contains the values for the first split, will return the new splitpoints
//llimits is left limits for the partitions
//rlimits is right limtis for the partitions
//newllimits will be the new left limits after the partitioning
//newrlimits will be the new right limits after the partitioning
//basellimits,baserlimits are the limits for the start of the partitioning. If only singleblockpartition is used, set these to the same as llimits and rlimits, but if used in combination with multiblockpartition, these are the original starting positions
//splitcount is the number of splits
//newpositions,newpositions2 are the permuted values for positions,positions2 (unless resultsinoriginal is set)
//ysplit will determine if split should be performed using positions or positions2, true means positions2
//resultsinoriginal use positions,positions2,indices as output values instead
//debugvector debug purpose only

//Note that even though results in original is set, newpositions etc. are needed, since the algorithm will alternate the current permutation between these two
{
  __shared__ int localleftcount[threadsperblock],localmiddlecount[threadsperblock];
  __shared__ int threadcount,blockaddstartleft,blockaddstartright,blockaddstartmiddle,count,splitstart,splitend,nractivethreads,lcount,rcount,newactive,innernewactive,localsplitstart,localsplitend,oldsplitstart,oldsplitend,innersplitstart,innersplitend,innersplitcount,mask,threesplit,fullthreesplit;
#ifdef DEBUGSINGLEBLOCKPARTITION
  __shared__ int lc;
  __shared__ SORT_REAL point1,point2,point3;
#endif
#ifdef MEDIAN_OF_32
  __shared__ SORT_REAL point[MEDIANCOUNT];
  __shared__ int pointcount[MEDIANCOUNT],sortedpoint[MEDIANCOUNT];
#else
  __shared__ SORT_REAL point1,point2,point3;
#endif
  __shared__ SORT_REAL *activepositions,*activenewpositions,*passivepositions,*passivenewpositions,splitpoint;
#ifdef SORTLIMIT
  __shared__ SORT_REAL leftlimitvalue,rightlimitvalue;
#endif
  int start,comparisonpoint,comparisonvalue,i,localleft,localright,addpoint,index,loopend,loopend2,defaultright;
//   __shared__ int indexcache[INDEXCACHELENGTH*threadsperblock]; //for cached case, shift the indices instead, since they require less storage space, and shared memory should be fast enough without being ordered
//   __shared__ int indexcache2[INDEXCACHELENGTH*threadsperblock];
  __shared__ int leftside,rightside;
//   __shared__ SORT_REAL positioncache[INDEXCACHELENGTH*threadsperblock],splitpoint;
  __shared__ int *indexcache2;
  __shared__ SORT_REAL *positioncache;
  if(threadIdx.x==0) {
    indexcache2=indexcache+INDEXCACHELENGTH*blockDim.x;
    positioncache=(SORT_REAL*)(indexcache+2*INDEXCACHELENGTH*blockDim.x);
  }
  for(int splitindex=blockIdx.x;splitindex<splitcount;splitindex+=gridDim.x) { //loop through all splits that this block has to do
    if(threadIdx.x==0) { //setup input values
      mask=1; //for threesplit, Different markings for different results
      threesplit=0; //only use threesplit when necessary, threesplit separates equal elements in a separate box (like quicksort3)
      fullthreesplit=0;
      defaultright=0;
      splitpoint=splitpoints[splitindex];
      oldsplitstart=localsplitstart=llimits[splitindex]; //it is necessary to remember the split at least one extra previous step in order to transfer the outside region properly
      oldsplitend=localsplitend=rlimits[splitindex];
#ifdef SORTLIMIT
      leftlimitvalue=leftlimitvalues[splitindex];
      rightlimitvalue=rightlimitvalues[splitindex];
#endif
      if(basellimits!=NULL) //to operate with multiblockpartition, it has to be able to determine the midpoint
        splitstart=basellimits[splitindex]; //If it takes over in the middle of another split, baselimits and limits may not be the same
      else
        splitstart=localsplitstart;
      if(baserlimits!=NULL)
        splitend=baserlimits[splitindex];
      else
        splitend=localsplitend;
      //these values can be set right away
      newllimits[splitindex<<1]=splitstart;
      newrlimits[(splitindex<<1)+1]=splitend;
      if(splitindex==splitcount-1) //last element, make it possible to get all indices from newllimits at the end, making rlimits useless after the split. This will make it like the CPU version of the FMM code
        newllimits[(splitindex<<1)+2]=splitend;
#ifdef SORTLIMIT
      if(isinf(leftlimitvalue))
        splitend+=(int)((splitend-splitstart-1)*disttarget); //note: splitend not part of interval, therefore -1
      if(isinf(rightlimitvalue))
        splitstart-=(int)((splitend-splitstart)*disttarget);
#endif
      count=localsplitend-localsplitstart;
      nractivethreads=blockDim.x;
      blockaddstartleft=localsplitstart;
      blockaddstartright=localsplitend;
      lcount=0;
      rcount=0;
      newactive=0;
      if(count>0) {
        while(count<nractivethreads) //always more points than threads, so each thread has atleast one element to split
          nractivethreads>>=1;
      }
      threadcount=(count+nractivethreads-1)/nractivethreads;
      if(ysplit[splitindex]) { //choose split direction (x/y)
        activepositions=positions2;
        activenewpositions=newpositions2;
        passivepositions=positions;
        passivenewpositions=newpositions;
      }
      else {
        activepositions=positions;
        activenewpositions=newpositions;
        passivepositions=positions2;
        passivenewpositions=newpositions2;
      }
      if(count==1)
        splitpoint=activepositions[localsplitstart]; //special case for count==1, since it will not set splitpoint properly otherwise cause no iteration is perform
#ifdef DEBUGSINGLEBLOCKPARTITION
    lc=0;
#endif
#ifdef SORTLIMIT
      if(splitpoint<leftlimitvalue)
        splitpoint=leftlimitvalue;
      if(splitpoint>rightlimitvalue)
        splitpoint=rightlimitvalue;
#endif
    }
    __syncthreads();
    //As long as the number of elements to split is larger than the internal cache, run this code, which loops through the elements, and writes back to global memory all the time
    while(count>INDEXCACHELENGTH*blockDim.x) { //when all elements cannot be cached
      for(i=0;i<threadcount;i+=INDEXCACHELENGTH) { //fill cache and split, then refill etc
        if(threadIdx.x<nractivethreads) {
          loopend=loopend2=imin(threadcount-i, INDEXCACHELENGTH); //number of elements to handle for each thread
          start=localsplitstart+threadIdx.x*loopend+i*nractivethreads;
          if(start+loopend>localsplitstart+count) //make sure it does not get elements outside the split
            loopend=localsplitstart+count-start;
          if(newactive) { //depending on if elements are in newindices etc or indices
            for(addpoint=0;addpoint<loopend;addpoint++) {
              indexcache[threadIdx.x*INDEXCACHELENGTH+addpoint]=newindices[addpoint+start];
              positioncache[threadIdx.x*INDEXCACHELENGTH+addpoint]=activenewpositions[addpoint+start];
            }
          }
          else {
            for(addpoint=0;addpoint<loopend;addpoint++) {
              indexcache[threadIdx.x*INDEXCACHELENGTH+addpoint]=indices[addpoint+start];
              positioncache[threadIdx.x*INDEXCACHELENGTH+addpoint]=activepositions[addpoint+start];
            }
          }
          localleftcount[threadIdx.x]=0;
          comparisonvalue=0;
          comparisonpoint=1;
          
          //make first loop and save results
          if(threesplit) { //more complicated if the split uses three equal as an own group
            localmiddlecount[threadIdx.x]=0;
            for(addpoint=0;addpoint<loopend;addpoint++, comparisonpoint<<=2) {
              if(positioncache[threadIdx.x*INDEXCACHELENGTH+addpoint]<splitpoint) { //0 for right side,1 for left side, 2 for middle
                localleftcount[threadIdx.x]++;
                comparisonvalue^=comparisonpoint;
              }
              else if(positioncache[threadIdx.x*INDEXCACHELENGTH+addpoint]==splitpoint){
                localmiddlecount[threadIdx.x]++;
                comparisonvalue^=(comparisonpoint<<1);
              }
            }
          }
          else {
            for(addpoint=0;addpoint<loopend;addpoint++, comparisonpoint<<=2) { //only 1 for left side, increase with 2 to interoperate with threesplit
              if(positioncache[threadIdx.x*INDEXCACHELENGTH+addpoint]<=splitpoint) {
                localleftcount[threadIdx.x]++;
                comparisonvalue^=comparisonpoint;
              }
            }
          }
        }
        __syncthreads();
        index=1;
        addpoint=1;
        
        //calculate culmulative sums, the number of elements on the right side can be calculated if the other ones are known
        //culmulative sums are necessary for the threads to determine where to place its own elements
        if(threesplit) {
          for(localleft=nractivethreads>>1;localleft>=1;localleft>>=1, index<<=1, addpoint++) {
            if(threadIdx.x<localleft) {
              localleftcount[(threadIdx.x<<addpoint)+(index<<1)-1]+=localleftcount[(threadIdx.x<<addpoint)+(index)-1];
              localmiddlecount[(threadIdx.x<<addpoint)+(index<<1)-1]+=localmiddlecount[(threadIdx.x<<addpoint)+(index)-1];
            }
            __syncthreads();
          }
          addpoint-=2;
          index=(nractivethreads>>2);
          for(localleft=2;index>=1;localleft<<=1, index>>=1, addpoint--) {
            if(threadIdx.x<localleft-1) {
              localleftcount[((threadIdx.x+1)<<addpoint)+index-1]+=localleftcount[((threadIdx.x+1)<<addpoint)-1];
              localmiddlecount[((threadIdx.x+1)<<addpoint)+index-1]+=localmiddlecount[((threadIdx.x+1)<<addpoint)-1];
            }
            __syncthreads();
          }
        }
        else {
          for(localleft=nractivethreads>>1;localleft>=1;localleft>>=1, index<<=1, addpoint++) {
            if(threadIdx.x<localleft)
              localleftcount[(threadIdx.x<<addpoint)+(index<<1)-1]+=localleftcount[(threadIdx.x<<addpoint)+(index)-1];
            __syncthreads();
          }
          addpoint-=2;
          index=(nractivethreads>>2);
          for(localleft=2;index>=1;localleft<<=1, index>>=1, addpoint--) {
            if(threadIdx.x<localleft-1)
              localleftcount[((threadIdx.x+1)<<addpoint)+index-1]+=localleftcount[((threadIdx.x+1)<<addpoint)-1];
            __syncthreads();
          }
        }
        if(threadIdx.x==0) { //set up the redistribution of elements
          addpoint=imin(count-i*nractivethreads, INDEXCACHELENGTH*nractivethreads);
          localleft=blockaddstartleft=lcount+localsplitstart;
          if(threesplit) {
            if(lcount+2*localleftcount[nractivethreads-1]+defaultright>rcount+addpoint-localmiddlecount[nractivethreads-1]) { //give the middle elements to the smallest count. Should quite evenly split them, even though all elements are equal. Could be done smarter by checking the best interval, but there is probably no use in optimizing this part, since no reasonable distribution should look like this, and all that is important is that the code works.
              lcount+=localleftcount[nractivethreads-1];
              rcount+=(addpoint-localleftcount[nractivethreads-1]);
              mask=1;
            }
            else { //mask is used to enable distribution of middle elements to either side
              lcount+=localleftcount[nractivethreads-1]+localmiddlecount[nractivethreads-1];
              rcount+=(addpoint-localleftcount[nractivethreads-1]-localmiddlecount[nractivethreads-1]);
              mask=3;
            }
          }
          else {
            lcount+=localleftcount[nractivethreads-1];
            rcount+=(addpoint-localleftcount[nractivethreads-1]);
          }
          localright=blockaddstartright=localsplitend-rcount;
        }
        
        __syncthreads();
        
        //determine the starting point for each thread to add elements
        if(threadIdx.x<nractivethreads) {
          if(threadIdx.x!=0) {
            localleft=blockaddstartleft+localleftcount[threadIdx.x-1];
            localright=blockaddstartright+loopend2*threadIdx.x-localleftcount[threadIdx.x-1]; //some localright will be outside interval, but these should not have any elements to add
            if(threesplit) {
              if(mask==3) { //give the middle elements to the proper side. No point putting them in the middle, since only the median is desired, not a full sorted array
                localright-=localmiddlecount[threadIdx.x-1];
                localleft+=localmiddlecount[threadIdx.x-1];
              }
            }
          }
          //distribute the elements. Note that fullthreesplit does not exist here, since it is only necessary at the end. The current code will always distribute the elements in a way that the count will be small enough to move to next step
          if(newactive) {
            for(index=0, comparisonpoint=mask;index<loopend;index++, comparisonpoint<<=2) {
              if(comparisonvalue&comparisonpoint)
                addpoint=localleft++;
              else
                addpoint=localright++;
              indices[addpoint]=indexcache[threadIdx.x*INDEXCACHELENGTH+index]; //the idea of addpoint is to keep the warp intact at the memory access point
              activepositions[addpoint]=positioncache[threadIdx.x*INDEXCACHELENGTH+index];
              passivepositions[addpoint]=passivenewpositions[start+index];
            }
          }
          else {
            for(index=0, comparisonpoint=mask;index<loopend;index++, comparisonpoint<<=2) {
              if(comparisonvalue&comparisonpoint)
                addpoint=localleft++;
              else
                addpoint=localright++;
              newindices[addpoint]=indexcache[threadIdx.x*INDEXCACHELENGTH+index]; //the idea of addpoint is to keep the warp intact at the memory access point
              activenewpositions[addpoint]=positioncache[threadIdx.x*INDEXCACHELENGTH+index];
              passivenewpositions[addpoint]=passivepositions[start+index];
            }
          }
        }
      }
      __syncthreads();
      if(threadIdx.x==0) { //split is done, set up the next one
        newactive^=1;
#ifdef DEBUGSINGLEBLOCKPARTITION
        DEBUGSINGLEBLOCKPARTITIONLINE1
#endif
        oldsplitend=localsplitend;
        oldsplitstart=localsplitstart;
        leftside=lcount+localsplitstart-splitstart;
        rightside=splitend-localsplitend+rcount;
        if(leftside>rightside) { //more elements on the left side
          localsplitend=localsplitend-rcount;
#ifdef SORTLIMIT
          if(splitpoint==leftlimitvalue) {
            if(leftside*0.5*(1-distlimit)>rightside*(1-0.5*(1-distlimit))) { //check if enough points are on each side
              leftlimitvalue=NEGINF; //if not, ignore geometric limit and move center position for the split
              splitend+=(int)((splitend-splitstart-1)*disttarget); //note -1 necessary here as splitend is not apart of the interval. (int) necessary to prevent roundoff errors
            }
            else {
              localsplitstart=localsplitstart+lcount;
            }
          }
#endif
          if(rcount==0) { //nothing happend in the split, move to threesplit mode and continue
            if(threesplit==1)
              defaultright=1;
            threesplit=1;
          }
        }
        else if(leftside==rightside) {
          localsplitend=localsplitend-rcount;
          localsplitstart=localsplitstart+lcount;
        }
        else {
          localsplitstart=localsplitstart+lcount;
#ifdef SORTLIMIT
          if(splitpoint==rightlimitvalue) {
            if(rightside*0.5*(1-distlimit)>leftside*(1-0.5*(1-distlimit))) {
              rightlimitvalue=INF;
              splitstart-=(int)((splitend-splitstart)*disttarget);
            }
            else
              localsplitend=localsplitend-rcount;
          }
#endif
          if(lcount==0) { //nothing happend in the split, move to threesplit mode and continue
            if(threesplit==1)
              defaultright=0;
            threesplit=1;
          }
        }
        count=localsplitend-localsplitstart;
        if(count>0) {
          while(count<nractivethreads) //not really necessary
            nractivethreads>>=1;
        }
        threadcount=(count+nractivethreads-1)/nractivethreads;
#ifdef MEDIAN_OF_32
      }
      if(localsplitend-localsplitstart>MEDIANCOUNT) {
        if(threadIdx.x<MEDIANCOUNT)
          point[threadIdx.x]=activepositions[localsplitstart+(localsplitend-localsplitstart-1)*threadIdx.x/(MEDIANCOUNT-1)];
        __syncthreads();
        if(threadIdx.x<MEDIANCOUNT) {
          pointcount[threadIdx.x]=0;
          for(int k=0;k<MEDIANCOUNT;k++) {
            if(point[threadIdx.x]>point[k]||point[threadIdx.x]==point[k]&&threadIdx.x>k)
              pointcount[threadIdx.x]++;
          }
        }
        __syncthreads();
        if(threadIdx.x<MEDIANCOUNT)
          sortedpoint[pointcount[threadIdx.x]]=threadIdx.x;
        __syncthreads();
        if(threadIdx.x==0) {
          splitpoint=point[sortedpoint[2+(int)((SORT_REAL)(((splitend+splitstart)>>1)-localsplitstart)*28/(SORT_REAL)(localsplitend-localsplitstart))]];
        }
      }
      else if(count!=0)
        if(threadIdx.x==0)
          splitpoint=activepositions[(localsplitstart+localsplitend)/2];
#ifdef SORTLIMIT
      if(splitpoint<leftlimitvalue)
        splitpoint=leftlimitvalue;
      if(splitpoint>rightlimitvalue)
        splitpoint=rightlimitvalue;
#endif
      if(threadIdx.x==0) {
        lcount=0;
        rcount=0;
      }
#else /*MEDAIN_OF_32*/
        if(localsplitend-localsplitstart>3) { //use median of three if enough elements are available
          if(newactive) {
            point1=activenewpositions[localsplitstart];
            point2=activenewpositions[(localsplitstart+localsplitend)/2];
            point3=activenewpositions[localsplitend-1];
          }
          else {
            point1=activepositions[localsplitstart];
            point2=activepositions[(localsplitstart+localsplitend)/2];
            point3=activepositions[localsplitend-1];
          }
          if(point1>point2) {
            SWAP(point1, point2, splitpoint);
          }
          if(point2>point3) {
            point2=point3;
          }
          if(point1>point2)
            splitpoint=point1;
          else
            splitpoint=point2;
        }
        else  //if not enough elements are available, split with middle element
          splitpoint= activepositions[(localsplitstart+localsplitend)/2];
#ifdef SORTLIMIT
        if(splitpoint<leftlimitvalue)
          splitpoint=leftlimitvalue;
        if(splitpoint>rightlimitvalue)
          splitpoint=rightlimitvalue;
#endif
        lcount=0;
        rcount=0;
      }
#endif
      __syncthreads();
      //copy elements not to be shifted directly to the other array to keep both synchronized at the end. This is only necessary for the output array
      if(!resultsinoriginal&&!newactive) { //only this direction useful, copy the elements not involved in next split to the output vector. The other one does not have to be valid
        for(addpoint=oldsplitstart+threadIdx.x;addpoint<localsplitstart;addpoint+=blockDim.x) {
          newindices[addpoint]=indices[addpoint];
          activenewpositions[addpoint]=activepositions[addpoint];
          passivenewpositions[addpoint]=passivepositions[addpoint];
        }
        for(addpoint=localsplitend+threadIdx.x;addpoint<oldsplitend;addpoint+=blockDim.x) {
          newindices[addpoint]=indices[addpoint];
          activenewpositions[addpoint]=activepositions[addpoint];
          passivenewpositions[addpoint]=passivepositions[addpoint];
        }
#ifdef DEBUGSINGLEBLOCKPARTITION
        DEBUGSINGLEBLOCKPARTITIONLINE2
#endif
      }
      if(resultsinoriginal&&newactive) { //only this direction useful
        for(addpoint=oldsplitstart+threadIdx.x;addpoint<localsplitstart;addpoint+=blockDim.x) {
          indices[addpoint]=newindices[addpoint];
          activepositions[addpoint]=activenewpositions[addpoint];
          passivepositions[addpoint]=passivenewpositions[addpoint];
        }
        for(addpoint=localsplitend+threadIdx.x;addpoint<oldsplitend;addpoint+=blockDim.x) {
          indices[addpoint]=newindices[addpoint];
          activepositions[addpoint]=activenewpositions[addpoint];
          passivepositions[addpoint]=passivenewpositions[addpoint];
        }
#ifdef DEBUGSINGLEBLOCKPARTITION
        DEBUGSINGLEBLOCKPARTITIONLINE2
#endif
      }
#ifdef DEBUGSINGLEBLOCKPARTITION
      DEBUGSINGLEBLOCKPARTITIONLINE3
#endif
      __syncthreads();
    }

    //Now, all elements will fit in the cache. To reduce memory reads, run the rest of the sort from cache, and write to memory aftewards
    //Cache all elements for the final stage.
    if(newactive) { //check which of the two arrays, that contains the most recent data, and load from that array
      for(addpoint=threadIdx.x;addpoint<localsplitend-localsplitstart;addpoint+=blockDim.x) {
        positioncache[addpoint]=activenewpositions[addpoint+localsplitstart];
        indexcache[addpoint]=addpoint; //create a permutation array for the cache
      }
    }
    else {
      for(addpoint=threadIdx.x;addpoint<localsplitend-localsplitstart;addpoint+=blockDim.x) {
        positioncache[addpoint]=activepositions[addpoint+localsplitstart];
        indexcache[addpoint]=addpoint; //create a permutation array for the cache
      }
    }
    if(threadIdx.x==0) {
      oldsplitstart=innersplitstart=0;
      innersplitcount=oldsplitend=innersplitend=count;
      innernewactive=0;
      
    }
    __syncthreads();

    //complete the split with all elements cached, and only split in cache. Otherwise, the method is basically as the other one
    while(count>1) {
      if(threadIdx.x<nractivethreads) {
        loopend=threadcount;
        start=innersplitstart+threadIdx.x*threadcount;
        if(start+loopend>innersplitstart+count)
          loopend=innersplitstart+count-start; //setup limits of all threads
        localleftcount[threadIdx.x]=0;
        comparisonvalue=0;
        comparisonpoint=1;
        if(threesplit) { //same as above. Here, the positions are not moved around to save space in shared memory. Since shared memory access does not require ordering, irregular accesses should be ok
          localmiddlecount[threadIdx.x]=0;
          if(innernewactive) {
            for(addpoint=0;addpoint<loopend;addpoint++, comparisonpoint<<=2) {
              if(positioncache[indexcache2[start+addpoint]]<splitpoint) {
                localleftcount[threadIdx.x]++;
                comparisonvalue^=comparisonpoint;
              }
              else if(positioncache[indexcache2[start+addpoint]]==splitpoint) {
                localmiddlecount[threadIdx.x]++;
                comparisonvalue^=(comparisonpoint<<1);
              }
            }
          }
          else {
            for(addpoint=0;addpoint<loopend;addpoint++, comparisonpoint<<=2) {
              if(positioncache[indexcache[start+addpoint]]<splitpoint) {
                localleftcount[threadIdx.x]++;
                comparisonvalue^=comparisonpoint;
              }
              else if(positioncache[indexcache[start+addpoint]]==splitpoint) {
                localmiddlecount[threadIdx.x]++;
                comparisonvalue^=(comparisonpoint<<1);
              }
            }
          }
        }
        else { //!threesplit
          if(innernewactive) {
            for(addpoint=0;addpoint<loopend;addpoint++, comparisonpoint<<=2) {
              if(positioncache[indexcache2[start+addpoint]]<=splitpoint) {
                localleftcount[threadIdx.x]++;
                comparisonvalue^=comparisonpoint;
              }
            }
          }
          else {
            for(addpoint=0;addpoint<loopend;addpoint++, comparisonpoint<<=2) {
              if(positioncache[indexcache[start+addpoint]]<=splitpoint) {
                localleftcount[threadIdx.x]++;
                comparisonvalue^=comparisonpoint;
              }
            }
          }
        }
      }
      
      __syncthreads();
      index=1;
      addpoint=1;
      
      //calculate culmulative sums. In this case, it may actually be necessary to split with all three cases, if the middle elements end up in the middle
      if(threesplit) {
        for(comparisonpoint=nractivethreads>>1;comparisonpoint>=1;comparisonpoint>>=1, index<<=1, addpoint++) {
          if(threadIdx.x<comparisonpoint) {
            localleftcount[(threadIdx.x<<addpoint)+(index<<1)-1]+=localleftcount[(threadIdx.x<<addpoint)+(index)-1];
            localmiddlecount[(threadIdx.x<<addpoint)+(index<<1)-1]+=localmiddlecount[(threadIdx.x<<addpoint)+(index)-1];
          }
          __syncthreads();
        }
        addpoint-=2;
        index=(nractivethreads>>2);
        for(comparisonpoint=2;index>=1;comparisonpoint<<=1, index>>=1, addpoint--) {
          if(threadIdx.x<comparisonpoint-1) {
            localleftcount[((threadIdx.x+1)<<addpoint)+index-1]+=localleftcount[((threadIdx.x+1)<<addpoint)-1];
            localmiddlecount[((threadIdx.x+1)<<addpoint)+index-1]+=localmiddlecount[((threadIdx.x+1)<<addpoint)-1];
          }
          __syncthreads();
        }
        if(threadIdx.x==0) { //setup the three split
          localleft=blockaddstartleft=innersplitstart;
          addpoint=innersplitstart+localsplitstart-splitstart+localleftcount[nractivethreads-1];//nr of elements to the left
          index=splitend-innersplitend-localsplitend+count+innersplitcount-localleftcount[nractivethreads-1];//nr of elements to the right
          if(addpoint>=index) { //give middle elements to right
            localright=blockaddstartright=innersplitstart+localleftcount[nractivethreads-1];
            mask=1;
            fullthreesplit=0;
          }
          else if(index-(localmiddlecount[nractivethreads-1]<<1)>addpoint) { //give middle elements to left
            localright=blockaddstartright=innersplitstart+localleftcount[nractivethreads-1]+localmiddlecount[nractivethreads-1];
            mask=3;
            fullthreesplit=0;
          }
          else { //the case where the middle elements actually ends up in the middle, split in three this time, and the split is done. This is new for this final loop
            fullthreesplit=1;
            index=blockaddstartmiddle=innersplitstart+localleftcount[nractivethreads-1];
            localright=blockaddstartright=blockaddstartmiddle+localmiddlecount[nractivethreads-1];
          }
        }
      }
      else { //normal split, calculate culmulative sum
        for(comparisonpoint=nractivethreads>>1;comparisonpoint>=1;comparisonpoint>>=1, index<<=1, addpoint++) {
          if(threadIdx.x<comparisonpoint)
            localleftcount[(threadIdx.x<<addpoint)+(index<<1)-1]+=localleftcount[(threadIdx.x<<addpoint)+(index)-1];
          __syncthreads();
        }
        addpoint-=2;
        index=(nractivethreads>>2);
        for(comparisonpoint=2;index>=1;comparisonpoint<<=1, index>>=1, addpoint--) {
          if(threadIdx.x<comparisonpoint-1)
            localleftcount[((threadIdx.x+1)<<addpoint)+index-1]+=localleftcount[((threadIdx.x+1)<<addpoint)-1];
          __syncthreads();
        }
        if(threadIdx.x==0) {
          localleft=blockaddstartleft=innersplitstart;
          localright=blockaddstartright=innersplitstart+localleftcount[nractivethreads-1];
        }
      }
      __syncthreads();
      if(threadIdx.x<nractivethreads) { //set starting point for each thread in the reordering
        if(threadIdx.x!=0) {
          localleft=blockaddstartleft+localleftcount[threadIdx.x-1];
          localright=blockaddstartright+threadcount*threadIdx.x-localleftcount[threadIdx.x-1]; //some localright will be outside interval, but these should not have any elements to add
          if(threesplit) {
            if(fullthreesplit) {
              index=blockaddstartmiddle+localmiddlecount[threadIdx.x-1];
              localright-=localmiddlecount[threadIdx.x-1];
            }
            else if(mask==3) {
              localright-=localmiddlecount[threadIdx.x-1];
              localleft+=localmiddlecount[threadIdx.x-1];
            }
          }
        }
        
        
        //reorder the elements
        if(fullthreesplit) { //the last split. The equal elements will end up in the middle, therefore, they can't be given to any of the sides
          if(innernewactive) {
            for(addpoint=0, comparisonpoint=1;addpoint<loopend;addpoint++, comparisonpoint<<=2) {
              if(comparisonvalue&comparisonpoint)
                indexcache[localleft++]=indexcache2[start+addpoint];
              else if(comparisonvalue&(comparisonpoint<<1))
                indexcache[index++]=indexcache2[start+addpoint];
              else
                indexcache[localright++]=indexcache2[start+addpoint];
            }
          }
          else {
            for(addpoint=0, comparisonpoint=1;addpoint<loopend;addpoint++, comparisonpoint<<=2) {
              if(comparisonvalue&comparisonpoint)
                indexcache2[localleft++]=indexcache[start+addpoint];
              else if(comparisonvalue&(comparisonpoint<<1))
                indexcache2[index++]=indexcache[start+addpoint];
              else
                indexcache2[localright++]=indexcache[start+addpoint];
            }
          }
        }
        else {  //normal split (not fullthreesplit)
          if(innernewactive) {
            for(addpoint=0, comparisonpoint=mask;addpoint<loopend;addpoint++, comparisonpoint<<=2) {
              if(comparisonvalue&comparisonpoint)
                indexcache[localleft++]=indexcache2[start+addpoint];
              else
                indexcache[localright++]=indexcache2[start+addpoint];
            }
          }
          else {
            for(addpoint=0, comparisonpoint=mask;addpoint<loopend;addpoint++, comparisonpoint<<=2) {
              if(comparisonvalue&comparisonpoint)
                indexcache2[localleft++]=indexcache[start+addpoint];
              else
                indexcache2[localright++]=indexcache[start+addpoint];
            }
          }
        }
      }
      __syncthreads();
      if(threadIdx.x==0) {  //split done, set up the next one
        innernewactive^=1;
        if(fullthreesplit) { //if in fullthreesplit mode, the partitioning is complete, since the middle element actually ended up in the middle
          if(!innernewactive) { //probably quite a useless if statement, since evaluating the expressions every time should not be that expensive
            oldsplitstart=innersplitstart;
            oldsplitend=innersplitend;
          }
          innersplitstart=innersplitend=((splitend-splitstart)>>1)-(localsplitstart-splitstart); //set the final endpoints of the split to the middle
          count=0;
        }
        else { //check which side to split next time
          leftside=innersplitstart+localsplitstart-splitstart+localleftcount[nractivethreads-1];
          rightside=splitend-innersplitend-localsplitend+count+innersplitcount-localleftcount[nractivethreads-1];
          if(threesplit&&mask==3) {
            leftside+=localmiddlecount[nractivethreads-1];
            rightside-=localmiddlecount[nractivethreads-1];
          }
#ifdef DEBUGSINGLEBLOCKPARTITION
          DEBUGSINGLEBLOCKPARTITIONLINE4
#endif
          if(!innernewactive) { //probably quite a useless if statement, since evaluating the expressions every time should not be that expensive
            oldsplitstart=innersplitstart;
            oldsplitend=innersplitend;
          }
          if(leftside!=rightside) {
            if(leftside>rightside) { //more elements on the left side
              innersplitend+=-count+localleftcount[nractivethreads-1];
              if(threesplit&&mask==3) {
                innersplitend+=localmiddlecount[nractivethreads-1];
              }
              if(count-localleftcount[nractivethreads-1]==0) //safety to prevent infinite loops, threesplit should always work, compared to the normal one, which can fail with several equal elements
                threesplit=1;
#ifdef SORTLIMIT
              if(splitpoint==leftlimitvalue) {
                if(leftside*0.5*(1-distlimit)>rightside*(1-0.5*(1-distlimit))) { //check if enough points are on each side
                  leftlimitvalue=NEGINF; //if not, ignore geometric limit and move center position for the split
                  splitend+=(int)((splitend-splitstart-1)*disttarget);
                }
                else
                  innersplitstart=innersplitend;
              }
#endif
            }
            else {
              innersplitstart+=localleftcount[nractivethreads-1];
              if(threesplit&&mask==3) {
                innersplitstart+=localmiddlecount[nractivethreads-1];
              }
              if(localleftcount[nractivethreads-1]==0) //safety, same as above, if no elements were moved, continue in threesplit mode, which always moves atleast one element
                threesplit=1;
#ifdef SORTLIMIT
              if(splitpoint==rightlimitvalue) {
                if(rightside*0.5*(1-distlimit)>leftside*(1-0.5*(1-distlimit))) {
                  rightlimitvalue=INF;
                  splitstart-=(int)((splitend-splitstart)*disttarget);
                }
                else
                  innersplitend=innersplitstart;
              }
#endif
            }
            count=innersplitend-innersplitstart;
            if(count>0) { //if count==0, the partitioning is done. Otherwise, choose a new splitpoint
              while(count<nractivethreads)  //reduce nr of active threads if number of elements left is small
                nractivethreads>>=1;
              threadcount=(count+nractivethreads-1)/nractivethreads;
#ifdef MEDIAN_OF_32
            }
          }
        }
      }
      __syncthreads();
      if(!fullthreesplit&&leftside!=rightside&&count>0) {
        if(count>MEDIANCOUNT) {
          if(threadIdx.x<MEDIANCOUNT) {
            if(innernewactive)
              point[threadIdx.x]=positioncache[indexcache2[innersplitstart+(innersplitend-innersplitstart-1)*threadIdx.x/(MEDIANCOUNT-1)]];
            else
              point[threadIdx.x]=positioncache[indexcache[innersplitstart+(innersplitend-innersplitstart-1)*threadIdx.x/(MEDIANCOUNT-1)]];
          }
          __syncthreads();
          if(threadIdx.x<MEDIANCOUNT) {
            pointcount[threadIdx.x]=0;
            for(int k=0;k<MEDIANCOUNT;k++) {
              if(point[threadIdx.x]>point[k]||point[threadIdx.x]==point[k]&&threadIdx.x>k)
                pointcount[threadIdx.x]++;
            }
          }
          __syncthreads();
          if(threadIdx.x<MEDIANCOUNT)
            sortedpoint[pointcount[threadIdx.x]]=threadIdx.x;
          __syncthreads();
          if(threadIdx.x==0) {
            splitpoint=point[sortedpoint[2+(int)((SORT_REAL)(((splitend+splitstart)>>1)-localsplitstart-innersplitstart)*28/(SORT_REAL)count)]];
          }
#ifdef SORTLIMIT
          if(splitpoint<leftlimitvalue)
            splitpoint=leftlimitvalue;
          if(splitpoint>rightlimitvalue)
            splitpoint=rightlimitvalue;
#endif
        }
        else { //less than 32 elements, finalize split
          if(threadIdx.x<count) {
            if(innernewactive)
              point[threadIdx.x]=positioncache[indexcache2[innersplitstart+threadIdx.x]];
            else
              point[threadIdx.x]=positioncache[indexcache[innersplitstart+threadIdx.x]];
          }
          __syncthreads();
          if(threadIdx.x<count) { //counting sort
            pointcount[threadIdx.x]=0;
            for(int k=0;k<count;k++) {
              if(point[threadIdx.x]>point[k]||point[threadIdx.x]==point[k]&&threadIdx.x>k)
                pointcount[threadIdx.x]++;
            }
          }
          __syncthreads();
          if(threadIdx.x<count) //permutation index
            sortedpoint[pointcount[threadIdx.x]]=threadIdx.x;
          __syncthreads();
          if(threadIdx.x<count){ //permuted array
            if(innernewactive)
              pointcount[threadIdx.x]=indexcache2[innersplitstart+threadIdx.x];
            else
              pointcount[threadIdx.x]=indexcache[innersplitstart+threadIdx.x];
          }
          __syncthreads();
          if(threadIdx.x<count) { //write sorted list back
            if(innernewactive)
              indexcache2[innersplitstart+threadIdx.x]=pointcount[sortedpoint[threadIdx.x]];
            else
              indexcache[innersplitstart+threadIdx.x]=pointcount[sortedpoint[threadIdx.x]];
          }
          __syncthreads();
          if(threadIdx.x==0) {
#ifdef SORTLIMIT
            mask=innersplitstart; //just reusing variable mask
#endif
            innersplitstart=((splitend+splitstart)>>1)-localsplitstart;
            if(innernewactive)
              splitpoint=positioncache[indexcache2[innersplitstart]];
            else
              splitpoint=positioncache[indexcache[innersplitstart]];
#ifdef SORTLIMIT //check if splitpoint is outside limits, and move point in this case
//note that distlimit is currently not implemented here, but since it can at most be 32 elements on the wrong side, this is acceptable as it would take longer time to check it as well
            while(splitpoint<leftlimitvalue&&innersplitstart<innersplitend-1){
              innersplitstart++;
              if(innernewactive)
                splitpoint=positioncache[indexcache2[innersplitstart]];
              else
                splitpoint=positioncache[indexcache[innersplitstart]];
            }
            while(splitpoint>rightlimitvalue&&innersplitstart>mask) {
              innersplitstart--;
              if(innernewactive)
                splitpoint=positioncache[indexcache2[innersplitstart]];
              else
                splitpoint=positioncache[indexcache[innersplitstart]];
            }
#endif
            innersplitend=innersplitstart;
            count=0;
          }
        }
      }
      else if(leftside==rightside&&threadIdx.x==0){ //done
        innersplitstart+=localleftcount[nractivethreads-1];
        innersplitend=innersplitstart; //should give the same as line above, but it is a more safe solution to make them equal
        count=0;
      }
          
          
              #else
              if(count>3) { //if more than three elements left, use median of 3 for new splitpoint
                if(innernewactive) {
                  point1=positioncache[indexcache2[innersplitstart]];
                  point2=positioncache[indexcache2[(innersplitstart+innersplitend)/2]];
                  point3=positioncache[indexcache2[innersplitend-1]];
                }
                else {
                  point1=positioncache[indexcache[innersplitstart]];
                  point2=positioncache[indexcache[(innersplitstart+innersplitend)/2]];
                  point3=positioncache[indexcache[innersplitend-1]];
                }
                if(point1>point2) {
                  SWAP(point1, point2, splitpoint);
                }
                if(point2>point3) {
                  point2=point3;
                }
                if(point1>point2)
                  splitpoint=point1;
                else
                  splitpoint=point2;
              }
              else { //less or equal than 3 elements left. Complete the partitioning here. One case for each number of elements left
                if(count==3) { //complete the sort on one thread
                  if(innernewactive) {
                    addpoint=indexcache2[innersplitstart];
                    index=indexcache2[innersplitstart+1];
                    localleft=indexcache2[innersplitstart+2];
                  }
                  else {
                    addpoint=indexcache[innersplitstart];
                    index=indexcache[innersplitstart+1];
                    localleft=indexcache[innersplitstart+2];
                  }
                  if(positioncache[addpoint]>positioncache[index]) {
                    SWAP(addpoint, index, localright);
                  }
                  if(positioncache[index]>positioncache[localleft]) {
                    SWAP(index, localleft, localright);
                  }
                  if(positioncache[addpoint]>positioncache[index]) {
                    SWAP(addpoint, index, localright);
                  }
                  if(innernewactive) {
                    indexcache2[innersplitstart]=addpoint;
                    indexcache2[innersplitstart+1]=index;
                    indexcache2[innersplitstart+2]=localleft;
                  }
                  else {
                    indexcache[innersplitstart]=addpoint;
                    indexcache[innersplitstart+1]=index;
                    indexcache[innersplitstart+2]=localleft;
                  }
                  addpoint=innersplitstart+localsplitstart-splitstart-(splitend-innersplitend-localsplitend+innersplitcount);
                  innersplitstart+=(3-addpoint)/2;
                  innersplitend=innersplitstart;
                  count=0;
                  
                }
                else if(count==2) {  //one case for each number of elements
                  if(innernewactive) {
                    addpoint=indexcache2[innersplitstart];
                    index=indexcache2[innersplitstart+1];
                    if(positioncache[addpoint]>positioncache[index]) {
                      SWAP(addpoint, index, localright);
                    }
                    indexcache2[innersplitstart]=addpoint;
                    indexcache2[innersplitstart+1]=index;
                  }
                  else {
                    addpoint=indexcache[innersplitstart];
                    index=indexcache[innersplitstart+1];
                    if(positioncache[addpoint]>positioncache[index]) {
                      SWAP(addpoint, index, localright);
                    }
                    indexcache[innersplitstart]=addpoint;
                    indexcache[innersplitstart+1]=index;
                  }
                  addpoint=innersplitstart+localsplitstart-splitstart-(splitend-innersplitend-localsplitend+innersplitcount);
                  innersplitstart+=(2-addpoint)/2;  //think it always is one on each side, which could be implemented directly
                  innersplitend=innersplitstart;
                  count=0;
                }
                else if(count==1) {
                  addpoint=innersplitstart+localsplitstart-splitstart-(splitend-innersplitend-localsplitend+innersplitcount);
                  innersplitstart+=(1-addpoint)/2;
                  innersplitend=innersplitstart;
                  count=0;
                }
                if(innernewactive) {
                  splitpoint=positioncache[indexcache2[innersplitstart]];
                }
                else {
                  splitpoint=positioncache[indexcache[innersplitstart]];
                }
#ifdef SORTLIMIT
                if(splitpoint<leftlimitvalue)
                  splitpoint=leftlimitvalue;
                if(splitpoint>rightlimitvalue)
                  splitpoint=rightlimitvalue;
#endif
              }
            }
          }
          else { //done
            innersplitstart+=localleftcount[nractivethreads-1];
            innersplitend=innersplitstart; //should give the same as line above, but it is a more safe solution to make them equal
            count=0;
          }
        }
      }
      #endif /*MEDIAN_OF_32*/

      __syncthreads();
      //copy elements not to be shifted directly to the other array to keep both synchronized at the end
#ifdef DEBUGSINGLEBLOCKPARTITION
      DEBUGSINGLEBLOCKPARTITIONLINE2
#endif
      if(!innernewactive) { //only this direction useful, copy elements to right vector (always indexcache2 in the inner loop)
        for(addpoint=oldsplitstart+threadIdx.x;addpoint<innersplitstart;addpoint+=blockDim.x) {
          indexcache2[addpoint]=indexcache[addpoint];
        }
        for(addpoint=innersplitend+threadIdx.x;addpoint<oldsplitend;addpoint+=blockDim.x) {
          indexcache2[addpoint]=indexcache[addpoint];
        }
      }
#ifdef DEBUGSINGLEBLOCKPARTITION
      DEBUGSINGLEBLOCKPARTITIONLINE5
#endif
      __syncthreads();
    }

    //here, the main loop is done
    if(!innernewactive) {  //is this necessary?
      for(addpoint=threadIdx.x+innersplitstart;addpoint<innersplitend;addpoint++)
        indexcache2[addpoint]=indexcache[addpoint];
    }
    __syncthreads();
    
    //write results back to memory. Start by caching newindices
    if(newactive) { //which element did it end up in, before the cached code started?
      for(addpoint=threadIdx.x;addpoint<innersplitcount;addpoint+=blockDim.x)
        indexcache[addpoint]=newindices[localsplitstart+addpoint];
    }
    else {
      for(addpoint=threadIdx.x;addpoint<innersplitcount;addpoint+=blockDim.x)
        indexcache[addpoint]=indices[localsplitstart+addpoint];
    }
    __syncthreads();
    
    //now it is written back
    if(resultsinoriginal) {
      for(addpoint=threadIdx.x;addpoint<innersplitcount;addpoint+=blockDim.x) {
        activepositions[localsplitstart+addpoint]=positioncache[indexcache2[addpoint]];
        indices[localsplitstart+addpoint]=indexcache[indexcache2[addpoint]];
      }
    }
    else {
      for(addpoint=threadIdx.x;addpoint<innersplitcount;addpoint+=blockDim.x) {
        activenewpositions[localsplitstart+addpoint]=positioncache[indexcache2[addpoint]];
        newindices[localsplitstart+addpoint]=indexcache[indexcache2[addpoint]];
      }
    }
    __syncthreads();
    //since it is no use caching the passive positions, reuse positioncache for these as well. This reordering is done directly from global to global memory
    if(newactive) {
      for(addpoint=threadIdx.x;addpoint<innersplitcount;addpoint+=blockDim.x)
        positioncache[addpoint]=passivenewpositions[localsplitstart+addpoint];
    }
    else {
      for(addpoint=threadIdx.x;addpoint<innersplitcount;addpoint+=blockDim.x)
        positioncache[addpoint]=passivepositions[localsplitstart+addpoint];
    }
    __syncthreads();
    if(resultsinoriginal) {
      for(addpoint=threadIdx.x;addpoint<innersplitcount;addpoint+=blockDim.x) {
        passivepositions[localsplitstart+addpoint]=positioncache[indexcache2[addpoint]];
      }
    }
    else {
      for(addpoint=threadIdx.x;addpoint<innersplitcount;addpoint+=blockDim.x) {
        passivenewpositions[localsplitstart+addpoint]=positioncache[indexcache2[addpoint]];
      }
    }
    __syncthreads();
    if(threadIdx.x==0) { //set output values
      splitpoints[splitindex]=splitpoint;
//       newllimits[splitindex<<1]=splitstart;
      newrlimits[splitindex<<1]=localsplitstart+innersplitstart;
      newllimits[(splitindex<<1)+1]=localsplitstart+innersplitstart;
//       newrlimits[(splitindex<<1)+1]=splitend;
//       if(splitindex==splitcount-1) //last element, make it possible to get all indices from newllimits at the end, making rlimits useless after the split. This will make it like the CPU version of the FMM code
//         newllimits[(splitindex<<1)+2]=splitend;
    }
  }
}

//this function performs a single split using multiple blocks per split. Works the same way as singleblockpartition in the splits, but does not have the inner loop to split several times
//do not call with more splits than blocks. 
__global__ void partitionsplit(SORT_REAL* positions, SORT_REAL *positions2, int* indices, int* newindices, SORT_REAL *splitpoints, int *llimits, int* rlimits, int *lrcount, int splitcount, SORT_REAL* newpositions, SORT_REAL *newpositions2, int* ysplit, int* threesplitvector DEBUGVECTORSTRING) 
//positions is input positions for x values
//positions2 is input positions for y values
//indices will be the permutation matrix
//newindices will contain the output indices
//splitpoints contains the values for the split
//llimits is left limits for the partitions
//rlimits is right limtis for the partitions
//lrcount will contain the number of elements on each side in the split
//splitcount is the number of splits
//newpositions,newpositions2 are the permuted values for positions,positions2
//ysplit will determine if split should be performed using positions or positions2, true means positions2
//threesplitvector will determine if threesplitmode should be used
//debugvector debug purpose only
{
  __shared__ int localleftcount[threadsperblock], localmiddlecount[threadsperblock];
  __shared__ int threadcount, blockaddstartleft, blockaddstartright, count, splitstart, splitend, nractivethreads, lcount, rcount, mask, threesplit,splitindex,blockindex,blockcount,splitblockcount,blockstart;
  __shared__ SORT_REAL *activepositions, *activenewpositions, *passivepositions, *passivenewpositions,splitpoint;
  int start, comparisonpoint, comparisonvalue, i, localleft, localright, addpoint, index, loopend, loopend2;
  int indexcache[INDEXCACHELENGTH]; //for cached case, shift the indices instead, since they require less storage space, and shared memory should be fast enough without being ordered
  SORT_REAL positioncache[INDEXCACHELENGTH]; //could probably be used as a local variable as well
  if(threadIdx.x==0) { //setup
    splitindex=blockIdx.x*splitcount/gridDim.x;
#ifdef EVENSPLIT //if all splits have the same number of blocks.
    splitblockcount=gridDim.x/splitcount;
    blockindex=blockIdx.x%splitcount;
#else
    blockindex=blockIdx.x*splitcount%gridDim.x;
    splitblockcount=1+blockindex/splitcount+(gridDim.x-blockindex-1)/splitcount; //number of blocks in this split
    blockindex/=splitcount;
#endif
    splitpoint=splitpoints[splitindex];
    splitstart=llimits[splitindex];
    splitend=rlimits[splitindex];
    count=splitend-splitstart;
    if(splitblockcount*blockDim.x>count) { //fill as many blocks as possible with atleast one elements per thread, let the rest be idle
      splitblockcount=count/blockDim.x;
      if(splitblockcount==0)
        splitblockcount=1;
      if(blockindex>=splitblockcount)
        blockcount=0;
      else
        blockcount=(count+splitblockcount-1)/splitblockcount;
    }
    else
      blockcount=(count+splitblockcount-1)/splitblockcount;
    nractivethreads=blockDim.x;

    blockstart=splitstart+blockcount*blockindex;
    if(blockstart+blockcount>splitend) //make sure to stay within the limits
      blockcount=splitend-blockstart;
    if(blockcount>0) {
      while(blockcount<nractivethreads)
        nractivethreads>>=1;
    }
    mask=1;
    threadcount=(blockcount+nractivethreads-1)/nractivethreads;
    if(threesplitvector!=NULL)
      threesplit=threesplitvector[splitindex];
    lcount=0;
    rcount=0;
    if(ysplit[splitindex]) { //split direction. This is necessary, since different splits can be in different directions
      activepositions=positions2;
      activenewpositions=newpositions2;
      passivepositions=positions;
      passivenewpositions=newpositions;
    }
    else {
      activepositions=positions;
      activenewpositions=newpositions;
      passivepositions=positions2;
      passivenewpositions=newpositions2;
    }
  }
  __syncthreads();

  //now, start the splitting
  for(i=0;i<threadcount;i+=INDEXCACHELENGTH) { //fill cache and split cache
    if(threadIdx.x<nractivethreads) {
      loopend=loopend2=imin(threadcount-i, INDEXCACHELENGTH);
      start=blockstart+threadIdx.x*loopend+i*nractivethreads;
      if(start+loopend>blockstart+blockcount)
        loopend=blockstart+blockcount-start;
      for(addpoint=0;addpoint<loopend;addpoint++) { //load into cache
        indexcache[addpoint]=indices[addpoint+start];
        positioncache[addpoint]=activepositions[addpoint+start];
      }
      
      //pass 1, check if values are larger or smaller than splitpoint, and count the elements
      localleftcount[threadIdx.x]=0;
      comparisonvalue=0;
      comparisonpoint=1;
      if(threesplit) {
        localmiddlecount[threadIdx.x]=0;
        for(addpoint=0;addpoint<loopend;addpoint++, comparisonpoint<<=2) {
          if(positioncache[addpoint]<splitpoint) {
            localleftcount[threadIdx.x]++;
            comparisonvalue^=comparisonpoint;
          }
          else if(positioncache[addpoint]==splitpoint){
            localmiddlecount[threadIdx.x]++;
            comparisonvalue^=(comparisonpoint<<1);
          }
        }
      }
      else {
        for(addpoint=0;addpoint<loopend;addpoint++, comparisonpoint<<=2) {
          if(positioncache[addpoint]<=splitpoint) {
            localleftcount[threadIdx.x]++;
            comparisonvalue^=comparisonpoint;
          }
        }
      }
    }
    __syncthreads();
    
    //calculate culmulative sum for the block internally
    index=1;
    addpoint=1;
    if(threesplit) { //two culmulative sums has to be calculated here
      for(localleft=nractivethreads>>1;localleft>=1;localleft>>=1, index<<=1, addpoint++) {
        if(threadIdx.x<localleft) {
          localleftcount[(threadIdx.x<<addpoint)+(index<<1)-1]+=localleftcount[(threadIdx.x<<addpoint)+(index)-1];
          localmiddlecount[(threadIdx.x<<addpoint)+(index<<1)-1]+=localmiddlecount[(threadIdx.x<<addpoint)+(index)-1];
        }
        __syncthreads();
      }
      addpoint-=2;
      index=(nractivethreads>>2);
      for(localleft=2;index>=1;localleft<<=1, index>>=1, addpoint--) {
        if(threadIdx.x<localleft-1) {
          localleftcount[((threadIdx.x+1)<<addpoint)+index-1]+=localleftcount[((threadIdx.x+1)<<addpoint)-1];
          localmiddlecount[((threadIdx.x+1)<<addpoint)+index-1]+=localmiddlecount[((threadIdx.x+1)<<addpoint)-1];
        }
        __syncthreads();
      }
    }
    else {
      for(localleft=nractivethreads>>1;localleft>=1;localleft>>=1, index<<=1, addpoint++) {
        if(threadIdx.x<localleft)
          localleftcount[(threadIdx.x<<addpoint)+(index<<1)-1]+=localleftcount[(threadIdx.x<<addpoint)+(index)-1];
        __syncthreads();
      }
      addpoint-=2;
      index=(nractivethreads>>2);
      for(localleft=2;index>=1;localleft<<=1, index>>=1, addpoint--) {
        if(threadIdx.x<localleft-1)
          localleftcount[((threadIdx.x+1)<<addpoint)+index-1]+=localleftcount[((threadIdx.x+1)<<addpoint)-1];
        __syncthreads();
      }
    }
    
    //now the internal culmulative sum has been calculated. Continue and make a global one over all blocks as well.
    //To avoid unnecessary CPU syncronization barriers, use atomic operations instead.
    //This causes the code to be not-deterministic, sicne the order of the atomic operations between the blocks ins unknown,
    //but should be faster than making a separate kernel call form the CPU for this purpose
    if(threadIdx.x==0) {
      addpoint=imin(blockcount-i*nractivethreads, INDEXCACHELENGTH*nractivethreads);
      if(threesplit) { //note that this part is not fully thread safe, since other blocks can update lrcound after the if statement. However, the ill effect is that the middle elements will be given to the wrong side, but it should still be a valid split. It may however cause that no elements gets moved if all are equal. In this case, move to singlethreadpartition and complete the work
        if(lrcount[splitindex*2]+localleftcount[nractivethreads-1]>lrcount[splitindex*2+1]+(addpoint-localleftcount[nractivethreads-1])) { //give the middle elements to the smallest count. Should quite evenly split them, even though all elements are equal. Could be done smarter by checking the best interval, but there is probably no use in optimizing this part, since no reasonable distribution should look like this, and all that is important is that the code works.
          lcount=atomicAdd(lrcount+splitindex*2, localleftcount[nractivethreads-1]); //the global culmulative sum. Note that this makes the code non-deterministic, since the order of the blocks varies
          rcount=atomicAdd(lrcount+splitindex*2+1, addpoint-localleftcount[nractivethreads-1]);
          mask=1;
          localright=blockaddstartright=splitend-rcount-(addpoint-localleftcount[nractivethreads-1]);
          localleft=blockaddstartleft=lcount+splitstart;
        }
        else {
          lcount=atomicAdd(lrcount+splitindex*2, localleftcount[nractivethreads-1]+localmiddlecount[nractivethreads-1]);
          rcount=atomicAdd(lrcount+splitindex*2+1, addpoint-localleftcount[nractivethreads-1]-localmiddlecount[nractivethreads-1]);
          mask=3;
          localright=blockaddstartright=splitend-rcount-(addpoint-localleftcount[nractivethreads-1]-localmiddlecount[nractivethreads-1]);
          localleft=blockaddstartleft=lcount+splitstart;
        }
      }
      else {
        
        //global culmulativ sum
        lcount=atomicAdd(lrcount+splitindex*2, localleftcount[nractivethreads-1]);
        rcount=atomicAdd(lrcount+splitindex*2+1, addpoint-localleftcount[nractivethreads-1]);
        localright=blockaddstartright=splitend-rcount-(addpoint-localleftcount[nractivethreads-1]);
        localleft=blockaddstartleft=lcount+splitstart;
      }
    }
    
    __syncthreads();
    if(threadIdx.x<nractivethreads) {
      if(threadIdx.x!=0) { //determine starting point for each thread (threadIdx.x==0 is already set up by previous if-statement)
        localleft=blockaddstartleft+localleftcount[threadIdx.x-1];
        localright=blockaddstartright+loopend2*threadIdx.x-localleftcount[threadIdx.x-1]; //some localright will be outside interval, but these should not have any elements to add
        if(threesplit) {
          if(mask==3) {
            localright-=localmiddlecount[threadIdx.x-1];
            localleft+=localmiddlecount[threadIdx.x-1];
          }
        }
      }
      //reorder the elements
      for(index=0, comparisonpoint=mask;index<loopend;index++, comparisonpoint<<=2) {
        if(comparisonvalue&comparisonpoint)
          addpoint=localleft++;
        else
          addpoint=localright++;
        newindices[addpoint]=indexcache[index]; //the idea of addpoint is to keep the warp intact at the memory access point
        activenewpositions[addpoint]=positioncache[index];
        passivenewpositions[addpoint]=passivepositions[start+index];
      }
    }
  }
}


//this function is the same as the multi-block version, except that it does run on one block, which simplifies slightly, and it does not support threesplit (it is not used in partitioning, only for splits of evaluation points and if eta is used)
//this is the simplest of the split implementations
__global__ void partitionsplitsinglethread(SORT_REAL* positions, SORT_REAL *positions2, int* indices, int* newindices, SORT_REAL *splitpoints, int *llimits, int* rlimits, int *lrcount, int splitcount, SORT_REAL* newpositions, SORT_REAL *newpositions2, int* ysplit DEBUGVECTORSTRING) {
  __shared__ int localleftcount[threadsperblock];
  __shared__ int threadcount, blockaddstartleft, blockaddstartright, count, splitstart, splitend, nractivethreads, lcount, rcount;
  __shared__ SORT_REAL *activepositions, *activenewpositions, *passivepositions, *passivenewpositions, splitpoint;
  int start, comparisonpoint, comparisonvalue, i, localleft, localright, addpoint, index, loopend, loopend2;
  int indexcache[INDEXCACHELENGTH]; //for cached case, shift the indices instead, since they require less storage space, and shared memory should be fast enough without being ordered
  SORT_REAL positioncache[INDEXCACHELENGTH]; //could probably be used as a local variable as well
  for(int splitindex=blockIdx.x;splitindex<splitcount;splitindex+=gridDim.x) { //handle more splits than blocks
    if(threadIdx.x==0) { //simplified setup
      splitpoint=splitpoints[splitindex];
      splitstart=llimits[splitindex];
      splitend=rlimits[splitindex];
      count=splitend-splitstart;
      nractivethreads=blockDim.x;
      if(count>0) {
        while(count<nractivethreads)
          nractivethreads>>=1;
      }
      threadcount=(count+nractivethreads-1)/nractivethreads;
      lcount=0;
      rcount=0;
      if(ysplit[splitindex]) { //choose direciton
        activepositions=positions2;
        activenewpositions=newpositions2;
        passivepositions=positions;
        passivenewpositions=newpositions;
      }
      else {
        activepositions=positions;
        activenewpositions=newpositions;
        passivepositions=positions2;
        passivenewpositions=newpositions2;
      }
    }
    __syncthreads();
    //start the split
    for(i=0;i<threadcount;i+=INDEXCACHELENGTH) {  //fill cache etc...
      if(threadIdx.x<nractivethreads) {
        loopend=loopend2=imin(threadcount-i, INDEXCACHELENGTH);
        start=splitstart+threadIdx.x*loopend+i*nractivethreads;
        if(start+loopend>splitstart+count)
          loopend=splitstart+count-start;
        for(addpoint=0;addpoint<loopend;addpoint++) {
          indexcache[addpoint]=indices[addpoint+start];
          positioncache[addpoint]=activepositions[addpoint+start];
        }
        localleftcount[threadIdx.x]=0;
        comparisonvalue=0;
        comparisonpoint=1;
        for(addpoint=0;addpoint<loopend;addpoint++, comparisonpoint<<=2) { //pass 1, count the elements for each thread
          if(positioncache[addpoint]<=splitpoint) {
            localleftcount[threadIdx.x]++;
            comparisonvalue^=comparisonpoint;
          }
        }
      }
      __syncthreads();
      //culmulative sum, to calculate where each thread should start adding elements
      index=1;
      addpoint=1;
      for(localleft=nractivethreads>>1;localleft>=1;localleft>>=1, index<<=1, addpoint++) {
        if(threadIdx.x<localleft)
          localleftcount[(threadIdx.x<<addpoint)+(index<<1)-1]+=localleftcount[(threadIdx.x<<addpoint)+(index)-1];
        __syncthreads();
      }
      addpoint-=2;
      index=(nractivethreads>>2);
      for(localleft=2;index>=1;localleft<<=1, index>>=1, addpoint--) {
        if(threadIdx.x<localleft-1)
          localleftcount[((threadIdx.x+1)<<addpoint)+index-1]+=localleftcount[((threadIdx.x+1)<<addpoint)-1];
        __syncthreads();
      }
      
      //set up the redistribution
      if(threadIdx.x==0) {
        addpoint=imin(count-i*nractivethreads, INDEXCACHELENGTH*nractivethreads);
        
        localleft=blockaddstartleft=lcount+splitstart;
        lcount+=localleftcount[nractivethreads-1]; //simplified calculation of lcount when atomicAdd isn't necessary
        rcount+=addpoint-localleftcount[nractivethreads-1];
        localright=blockaddstartright=splitend-rcount;
      }
      
      __syncthreads();
      if(threadIdx.x<nractivethreads) { //determine starting point of each thread
        if(threadIdx.x!=0) { //threadIdx.x==0 is set up in previous if-statement, starting points for the threads
          localleft=blockaddstartleft+localleftcount[threadIdx.x-1];
          localright=blockaddstartright+loopend2*threadIdx.x-localleftcount[threadIdx.x-1]; //some localright will be outside interval, but these should not have any elements to add
        }
        
        //reorder the elements
        for(index=0, comparisonpoint=1;index<loopend;index++, comparisonpoint<<=2) {
          if(comparisonvalue&comparisonpoint)
            addpoint=localleft++;
          else
            addpoint=localright++;
          newindices[addpoint]=indexcache[index]; //the idea of addpoint is to keep the warp intact at the memory access point
          activenewpositions[addpoint]=positioncache[index];
          passivenewpositions[addpoint]=passivepositions[start+index];
        }
      }
    }
    if(threadIdx.x==0) { //to keep output consistent with multithreaded version
      lrcount[splitindex*2]=lcount;
      lrcount[splitindex*2+1]=rcount;
    }
    __syncthreads();// to avoid changing the pointers activenewpositions etc
  }
}

//this function is used together with partition split. Makes the preparation step otherwise done inside singleblockpartition
//SORTLIMITS not properly implemented in this one, preparesplit32 should be used instead
__global__ void preparesplit(SORT_REAL* positions, SORT_REAL *positions2,int* indices, SORT_REAL *splitpoints, int *llimits, int* rlimits, int* newllimits,int* newrlimits,int* basellimits,int* baserlimits, int *lrcount, int splitcount, int* ysplit, int* threesplitvector,int * outputvector SORTLIMITSTRING DEBUGVECTORSTRING) 
{
  int llimit, rlimit, lcount, rcount, basellimit, baserlimit;
  SORT_REAL point1, point2, point3, splitpoint;
  SORT_REAL *activepositions;
  __shared__ int  maxelements[PREPARESPLITTHREADCOUNT];
#ifdef SORTLIMIT
  SORT_REAL leftlimitvalue,rightlimitvalue;
#endif
  maxelements[threadIdx.x]=0;
  for(int i=threadIdx.x+blockIdx.x*blockDim.x;i<splitcount;i+=blockDim.x*gridDim.x) {
    llimit=llimits[i];
    rlimit=rlimits[i]; //should really be llimit+lcount+rcount
    if(llimit!=rlimit) { //if there is any split to perform
      if(ysplit[i])
        activepositions=positions2;
      else
        activepositions=positions;
      basellimit=basellimits[i];
      baserlimit=baserlimits[i];
#ifdef SORTLIMIT
      leftlimitvalue=leftlimitvalues[i];
      rightlimitvalue=rightlimitvalues[i];
#endif
      lcount=lrcount[2*i];
      rcount=lrcount[2*i+1];
      if(lcount==0||rcount==0) { //if there were no split, this is a safety to prevent deadlock
        if(threesplitvector[i]) //if it were in three split mode, this algorithm does not work, move to singlethreadpartition instead, which should always work
          outputvector[1]=1; //this should move the algorithm to three split mode
        else
          threesplitvector[i]=1; //attempt to run three split instead
      }
      if(llimit-basellimit+lcount>baserlimit-rlimit+rcount) //choose side to make next split
        rlimit-=rcount;
      else
        llimit+=lcount;
      newllimits[i]=llimit;
      newrlimits[i]=rlimit;
      if(rlimit-llimit>=3) { //complete the split if less or equal than 3 elements
        point1=activepositions[llimit];
        point2=activepositions[(llimit+rlimit)/2];
        point3=activepositions[rlimit-1];
        if(point1>point2) {
          SWAP(point1, point2, splitpoint);
        }
        if(point2>point3) {
          point2=point3;
        }
        if(point1>point2)
          splitpoint=point1;
        else
          splitpoint=point2;
        splitpoints[i]=splitpoint;
      }
      else if(rlimit-llimit==2) { //should be quite rare
        point1=activepositions[llimit];
        point2=activepositions[llimit+1];
        if(point1>point2) {
          activepositions[llimit]=point2;
          activepositions[llimit+1]=point1;
          if(activepositions==positions) {
            SWAP(positions2[llimit], positions2[llimit+1], splitpoint);
          }
          else {
            SWAP(positions[llimit], positions[llimit+1], splitpoint);
          }
          rcount=indices[llimit];
          indices[llimit]=indices[llimit+1];
          indices[llimit+1]=rcount;
        }
        rlimit=llimit-basellimit-(baserlimit-rlimit);
        llimit+=(2-rlimit)/2;
        rlimit=llimit;
        newllimits[i]=llimit;
        newrlimits[i]=rlimit;
        splitpoints[i]=activepositions[llimit];
        
      }
      else {
        rlimit=llimit-basellimit-(baserlimit-rlimit);
        llimit+=(1-rlimit)/2;
        rlimit=llimit;
        splitpoints[i]=activepositions[llimit];
        newllimits[i]=llimit;
        newrlimits[i]=rlimit;
      }
    }
    else { //no split. Just copy the limits
      newllimits[i]=llimits[i];
      newrlimits[i]=rlimits[i];
    }
#ifdef SORTLIMIT
    if(splitpoint<leftlimitvalue)
      splitpoint=leftlimitvalue;
    if(splitpoint>rightlimitvalue)
      splitpoint=rightlimitvalue;
#endif
    if(rlimit-llimit>maxelements[threadIdx.x]) //determine the number of elements in the split
      maxelements[threadIdx.x]=rlimit-llimit;
  }
  __syncthreads();
  
  //the maximum for the entire block
  for(basellimit=(blockDim.x>>1);basellimit>0;basellimit>>=1) {
    if(threadIdx.x<basellimit) {
      if(maxelements[threadIdx.x]<maxelements[threadIdx.x+basellimit])
        maxelements[threadIdx.x]=maxelements[threadIdx.x+basellimit];
    }
    __syncthreads();
  }
  if(threadIdx.x==0)
    atomicMax(outputvector, maxelements[threadIdx.x]); //the maximum for all blocks
}

//Similar to preparesplit, but uses median of 32 instead
__global__ void preparesplit32(SORT_REAL* positions, SORT_REAL *positions2,int* indices, SORT_REAL *splitpoints, int *llimits, int* rlimits, int* newllimits,int* newrlimits,int* basellimits,int* baserlimits, int *lrcount, int splitcount, int* ysplit, int* threesplitvector,int * outputvector SORTLIMITSTRING DEBUGVECTORSTRING) 
{
  __shared__ int llimit, rlimit, lcount, rcount, basellimit, baserlimit;
  __shared__ SORT_REAL point[PREPARESPLIT32THREADCOUNT];
  __shared__ SORT_REAL *activepositions;
  __shared__ int  maxelements,pointcount[PREPARESPLIT32THREADCOUNT],sortedpoint[PREPARESPLIT32THREADCOUNT]/*[PREPARESPLITTHREADCOUNT]*/;
  int rightside,leftside;
#ifdef SORTLIMIT
  __shared__ SORT_REAL leftlimitvalue,rightlimitvalue,splitpoint;
  __shared__ int tmp;
#endif
  if(threadIdx.x==0)
    maxelements=0;
  for(int i=blockIdx.x;i<splitcount;i+=gridDim.x) {
    if(threadIdx.x==0) {
      llimit=llimits[i];
      rlimit=rlimits[i]; //should really be llimit+lcount+rcount
      if(llimit!=rlimit) { //if there is any split to perform
        if(ysplit[i])
          activepositions=positions2;
        else
          activepositions=positions;
        basellimit=basellimits[i];
        baserlimit=baserlimits[i];
#ifdef SORTLIMIT
        leftlimitvalue=leftlimitvalues[i];
        rightlimitvalue=rightlimitvalues[i];
        splitpoint=splitpoints[i];
        if(isinf(leftlimitvalue)) //if previous step triggered this
          baserlimit+=(int)((baserlimit-basellimit-1)*disttarget);
        if(isinf(rightlimitvalue))
          basellimit-=(int)((baserlimit-basellimit)*disttarget);
#endif
        lcount=lrcount[2*i];
        rcount=lrcount[2*i+1];
        if(lcount==0||rcount==0) { //if there were no split, this is a safety to prevent deadlock
          if(threesplitvector[i]) //if it were in three split mode, this algorithm does not work, move to singlethreadpartition instead, which should always work
            outputvector[1]=1; //this should move the algorithm to three split mode
          else
            threesplitvector[i]=1; //attempt to run three split instead
        }
        leftside=llimit-basellimit+lcount;
        rightside=baserlimit-rlimit+rcount;
        if(leftside>rightside) {//choose side to make next split
          rlimit-=rcount;
#ifdef SORTLIMIT
          if(splitpoint==leftlimitvalue) {
            if(leftside*0.5*(1-distlimit)>rightside*(1-0.5*(1-distlimit))) { //check if enough points are on each side
              leftlimitvalue=leftlimitvalues[i]=NEGINF; //if not, ignore geometric limit and move center position for the split
              baserlimit+=(int)((baserlimit-basellimit-1)*disttarget);
            }
            else
              llimit+=lcount; //should be same ar llimit=rlimit;
          }
#endif
        }
        else {
          llimit+=lcount;
#ifdef SORTLIMIT
          if(splitpoint==rightlimitvalue) {
            if(rightside*0.5*(1-distlimit)>leftside*(1-0.5*(1-distlimit))) {
              rightlimitvalue=rightlimitvalues[i]=INF;
              basellimit-=(int)((baserlimit-basellimit)*disttarget);
            }
            else
              rlimit-=rcount;
          }
#endif
        }
        newllimits[i]=llimit;
        newrlimits[i]=rlimit;
      }
    }
    __syncthreads();
    if(rlimit-llimit>PREPARESPLIT32THREADCOUNT) {
      point[threadIdx.x]=activepositions[llimit+(rlimit-llimit-1)*threadIdx.x/(blockDim.x-1)];
      __syncthreads();
      pointcount[threadIdx.x]=0;
      for(int k=0;k<PREPARESPLIT32THREADCOUNT;k++) {
        if(point[threadIdx.x]>point[k]||point[threadIdx.x]==point[k]&&threadIdx.x>k)
          pointcount[threadIdx.x]++;
      }
      __syncthreads();
      sortedpoint[pointcount[threadIdx.x]]=threadIdx.x;
       __syncthreads();
        
      if(threadIdx.x==0) {
        splitpoints[i]=point[sortedpoint[2+(int)((SORT_REAL)(((baserlimit+basellimit)>>1)-llimit)*28/(SORT_REAL)(rlimit-llimit))]];
#ifdef SORTLIMIT
        if(splitpoints[i]<leftlimitvalue)
          splitpoints[i]=leftlimitvalue;
        if(splitpoints[i]>rightlimitvalue)
          splitpoints[i]=rightlimitvalue;
#endif
      }
    }
    else {//less than 32 elements, finalize split
      if(rlimit!=llimit) {
        if(threadIdx.x<rlimit-llimit) {
          point[threadIdx.x]=activepositions[llimit+threadIdx.x];
        }
        
        __syncthreads();
        pointcount[threadIdx.x]=0;
        if(threadIdx.x<rlimit-llimit) {//counting sort
          for(int k=0;k<rlimit-llimit;k++) {
            if(point[threadIdx.x]>point[k]||point[threadIdx.x]==point[k]&&threadIdx.x>k)
              pointcount[threadIdx.x]++;
          }
        }
        __syncthreads();
        if(threadIdx.x<rlimit-llimit) { //permutation index
          sortedpoint[pointcount[threadIdx.x]]=threadIdx.x;
        }
        __syncthreads();
        if(threadIdx.x<rlimit-llimit) { //permuted array
          activepositions[llimit+threadIdx.x]=point[sortedpoint[threadIdx.x]];
        }
        __syncthreads();
        if(activepositions==positions) {//read points
          if(threadIdx.x<rlimit-llimit) {
            point[threadIdx.x]=positions2[llimit+threadIdx.x];
          }
          __syncthreads();
          if(threadIdx.x<rlimit-llimit) {//write sorted list back
            positions2[llimit+threadIdx.x]=point[sortedpoint[threadIdx.x]];
          }
        }
        else { //positions instead of positions2
          if(threadIdx.x<rlimit-llimit) {
            point[threadIdx.x]=positions[llimit+threadIdx.x];
          }
          __syncthreads();
          if(threadIdx.x<rlimit-llimit) {
            positions[llimit+threadIdx.x]=point[sortedpoint[threadIdx.x]];
          }
        }
        if(threadIdx.x<rlimit-llimit) {
          pointcount[threadIdx.x]=indices[llimit+threadIdx.x]; //reuse of pointcount for temporary storage
        }
        __syncthreads();
        if(threadIdx.x<rlimit-llimit) {
          indices[llimit+threadIdx.x]=pointcount[sortedpoint[threadIdx.x]];
        }
        if(threadIdx.x==0) {
          splitpoints[i]=activepositions[((baserlimit+basellimit)>>1)];
#ifdef SORTLIMIT
          tmp=((baserlimit+basellimit)>>1);
          while(splitpoints[i]<leftlimitvalue&&tmp<rlimit-1){
            tmp++;
            splitpoints[i]=activepositions[tmp];
          }
          while(splitpoints[i]>rightlimitvalue&&tmp>llimit) {
            tmp--;
            splitpoints[i]=activepositions[tmp];
          }
          rlimit=llimit=tmp;
#else
          rlimit=llimit=((baserlimit+basellimit)>>1);
#endif
          newllimits[i]=llimit;
          newrlimits[i]=rlimit;
        }
      }
      else { //no split. Just copy the limits
        newllimits[i]=llimit;
        newrlimits[i]=rlimit;
      }
    }
    __syncthreads();
    if(rlimit-llimit>maxelements&&threadIdx.x==0) //determine the number of elements in the split
      maxelements=rlimit-llimit;
    __syncthreads();
  }
  __syncthreads();
  
  if(threadIdx.x==0)
    atomicMax(outputvector, maxelements); //the maximum for all blocks
}

//this is a copy function to move elements between positions and newpositions etc. This to just copy one half of the elements and split the other one. Will copy elements moved outside the split region
__global__ void splitoutsidecopy(int* llimits,int* rlimits,int* newllimits,int* newrlimits,SORT_REAL* positions,SORT_REAL *positions2,int* indices,SORT_REAL* newpositions,SORT_REAL *newpositions2,int* newindices,int splitcount DEBUGVECTORSTRING)
{
  __shared__ int splitindex, splitblockcount, blockindex, splitstart, splitend, count, newsplitblockcount, blockcount, blockstart;
  if(threadIdx.x==0) {  //setup, copy left side first
    splitindex=blockIdx.x*splitcount/gridDim.x;
#ifdef EVENSPLIT //if all splits have the same number of blocks.
    splitblockcount=gridDim.x/splitcount;
    blockindex=blockIdx.x%splitcount;
#else
    blockindex=blockIdx.x*splitcount%gridDim.x;
    splitblockcount=1+blockindex/splitcount+(gridDim.x-blockindex-1)/splitcount; //number of blocks in this split
    blockindex/=splitcount;
#endif
    splitstart=llimits[splitindex]; //load the limits of both previous steps
    splitend=newllimits[splitindex];
    count=splitend-splitstart;
    if(splitblockcount*blockDim.x>count) { //if uncecessary number of blocks for this split, atleast all threads should have one element
      newsplitblockcount=count/blockDim.x; //decrease number of active blocks
      if(newsplitblockcount==0)
        newsplitblockcount=1;
      if(blockindex>=newsplitblockcount)
        blockcount=0;
      else
        blockcount=(count+newsplitblockcount-1)/newsplitblockcount;
    }
    else
      blockcount=(count+splitblockcount-1)/splitblockcount; //number of elements in this block
    blockstart=splitstart+blockcount*blockindex;
    if(blockstart+blockcount>splitend) //set limits
      blockcount=splitend-blockstart;
  }
  __syncthreads();
  //make the copy of all elements that belonged to the old split limits, but not the new ones
  for(int i=blockstart+threadIdx.x;i<blockstart+blockcount;i+=blockDim.x) {
    newpositions[i]=positions[i];
    newpositions2[i]=positions2[i];
    newindices[i]=indices[i];
  }
  __syncthreads();
  if(threadIdx.x==0) { //setup again for right side. In practice, only one of the left or right side should have elements
    splitstart=newrlimits[splitindex];
    splitend=rlimits[splitindex];
    count=splitend-splitstart;
    if(splitblockcount*blockDim.x>count) { //same setup as above
      newsplitblockcount=count/blockDim.x;
      if(newsplitblockcount==0)
        newsplitblockcount=1;
      if(blockindex>=newsplitblockcount)
        blockcount=0;
      else
        blockcount=(count+newsplitblockcount-1)/newsplitblockcount;
    }
    else
      blockcount=(count+splitblockcount-1)/splitblockcount;
    blockstart=splitstart+blockcount*blockindex;
    if(blockstart+blockcount>splitend)
      blockcount=splitend-blockstart;
  }
  __syncthreads(); //the copy of the other side
  for(int i=blockstart+threadIdx.x;i<blockstart+blockcount;i+=blockDim.x) {
    newpositions[i]=positions[i];
    newpositions2[i]=positions2[i];
    newindices[i]=indices[i];
  }
}
/*------------------------------------------------------------------------*/

//determines z and d for a box. This version uses one block for each box.
extern __shared__ SORT_REAL maxxvalues[];
__global__ void findzd(int* llimits,int* rlimits,SORT_REAL *xpositions,SORT_REAL *ypositions,int* llimitsNE,int* rlimitsNE,SORT_REAL *xpositionsNE,SORT_REAL *ypositionsNE,SORT_DCMPLX *z,SORT_DCMPLX *d,SORT_REAL *dabs,int count)
{
  __shared__ SORT_REAL *maxyvalues,*minxvalues,*minyvalues,xbase,ybase;
  __shared__ int blockstart, blockend, nractivethreads, blockstartNE, blockendNE, blockcount, blockcountNE;
  int i;
  SORT_REAL value;
  #ifdef RADIALSHRINK
  SORT_REAL value2;
  #endif
  if(threadIdx.x==0) {
    maxyvalues=maxxvalues+blockDim.x;
    minxvalues=maxxvalues+2*blockDim.x;
    minyvalues=maxxvalues+3*blockDim.x;
  }
  for(int index=blockIdx.x;index<count;index+=gridDim.x) { //if more boxes than blocks
    if(threadIdx.x==0) { //setup
      blockstart=llimits[index]; //limits
      blockend=rlimits[index];
      blockcount=blockend-blockstart;
      nractivethreads=blockDim.x;
      if(xpositionsNE!=NULL) { //if evaluation points exists, set up these as well
        blockstartNE=llimitsNE[index];
        blockendNE=rlimitsNE[index];
        blockcountNE=blockendNE-blockstartNE;
        if(blockcount>blockcountNE) {
          if(blockcount>0) {
            while(nractivethreads>blockcount)
              nractivethreads>>=1;
          }
        }
        else {
          if(blockcountNE>0) {
            while(nractivethreads>blockcountNE)
              nractivethreads>>=1;
          }
        }
      }
      else {
        blockcountNE=0;
        if(blockcount>0) {
          while(nractivethreads>blockcount)
            nractivethreads>>=1;
        }
      }
    }
    __syncthreads();
    
    //start by calculating maximum values for each thread
    if(blockcount>0) {
      if(threadIdx.x<nractivethreads) {
        if(threadIdx.x<blockcount) { //initial value
          maxxvalues[threadIdx.x]=minxvalues[threadIdx.x]=xpositions[blockstart+threadIdx.x];
          maxyvalues[threadIdx.x]=minyvalues[threadIdx.x]=ypositions[blockstart+threadIdx.x];
        }
        //determine max/min
        for(i=blockstart+nractivethreads+threadIdx.x;i<blockend;i+=nractivethreads) {
          value=xpositions[i];
          if(value>maxxvalues[threadIdx.x])
            maxxvalues[threadIdx.x]=value;
          if(value<minxvalues[threadIdx.x])
            minxvalues[threadIdx.x]=value;
          value=ypositions[i];
          if(value>maxyvalues[threadIdx.x])
            maxyvalues[threadIdx.x]=value;
          if(value<minyvalues[threadIdx.x])
            minyvalues[threadIdx.x]=value;
        }
      }
    }
    if(xpositionsNE!=NULL) { //include evaluation points for each thread
      if(threadIdx.x<nractivethreads) {
        //if value has not been set during first step
        if(threadIdx.x>=blockcount&&threadIdx.x<blockcountNE) { //is it possible for the second condition to be false if the first one is true???
          maxxvalues[threadIdx.x]=minxvalues[threadIdx.x]=xpositionsNE[blockstartNE+threadIdx.x];
          maxyvalues[threadIdx.x]=minyvalues[threadIdx.x]=ypositionsNE[blockstartNE+threadIdx.x];
        }
        //determine max/min
        for(i=blockstartNE+threadIdx.x;i<blockendNE;i+=nractivethreads) {
          value=xpositionsNE[i];
          if(value>maxxvalues[threadIdx.x])
            maxxvalues[threadIdx.x]=value;
          if(value<minxvalues[threadIdx.x])
            minxvalues[threadIdx.x]=value;
          value=ypositionsNE[i];
          if(value>maxyvalues[threadIdx.x])
            maxyvalues[threadIdx.x]=value;
          if(value<minyvalues[threadIdx.x])
            minyvalues[threadIdx.x]=value;
        }
      }
    }
    //if there were any elements. Calculate max/min for the block
    if(blockcount>0||blockcountNE>0) {
      __syncthreads();
      for(i=(nractivethreads>>1);i>0;i>>=1) {
        if(threadIdx.x<i) {
          if(maxxvalues[threadIdx.x+i]>maxxvalues[threadIdx.x])
            maxxvalues[threadIdx.x]=maxxvalues[threadIdx.x+i];
          if(minxvalues[threadIdx.x+i]<minxvalues[threadIdx.x])
            minxvalues[threadIdx.x]=minxvalues[threadIdx.x+i];
          if(maxyvalues[threadIdx.x+i]>maxyvalues[threadIdx.x])
            maxyvalues[threadIdx.x]=maxyvalues[threadIdx.x+i];
          if(minyvalues[threadIdx.x+i]<minyvalues[threadIdx.x])
            minyvalues[threadIdx.x]=minyvalues[threadIdx.x+i];
        }
        __syncthreads();
      }
      if(threadIdx.x==0) { //write results
        #ifndef RADIALSHRINK
        xbase=(maxxvalues[0]-minxvalues[0])*0.5;
        ybase=(maxyvalues[0]-minyvalues[0])*0.5;
        COMPLEXASSIGN(d[index], xbase, ybase);
        dabs[index]=hypot(xbase,ybase);
        #else
        COMPLEXASSIGN(d[index], (maxxvalues[0]-minxvalues[0])*0.5, (maxyvalues[0]-minyvalues[0])*0.5);
        #endif
        xbase=(minxvalues[0]+maxxvalues[0])*0.5;
        ybase=(minyvalues[0]+maxyvalues[0])*0.5;
        COMPLEXASSIGN(z[index], xbase, ybase);
      }
      #ifdef RADIALSHRINK
      __syncthreads();
      if(blockcount+blockcountNE<RADIALSHRINKLIMIT) {
        if(blockcount>0) {
          if(threadIdx.x<nractivethreads) {
            if(threadIdx.x<blockcount) { //initial value
              value=xbase-xpositions[blockstart+threadIdx.x];
              value2=ybase-ypositions[blockstart+threadIdx.x];
              maxxvalues[threadIdx.x]=value*value+value2*value2;
            }
            //determine max/min
            for(i=blockstart+nractivethreads+threadIdx.x;i<blockend;i+=nractivethreads) {
              value=xbase-xpositions[i];
              value2=ybase-ypositions[i];
              value=value*value+value2*value2;
              if(value>maxxvalues[threadIdx.x])
                maxxvalues[threadIdx.x]=value;
            }
          }
        }
        if(xpositionsNE!=NULL) { //include evaluation points for each thread
          if(threadIdx.x<nractivethreads) {
            //if value has not been set during first step
            if(threadIdx.x>=blockcount&&threadIdx.x<blockcountNE) { //is it possible for the second condition to be false if the first one is true???
              value=xbase-xpositionsNE[blockstartNE+threadIdx.x];
              value2=ybase-ypositionsNE[blockstartNE+threadIdx.x];
              maxxvalues[threadIdx.x]=value*value+value2*value2;
            }
            //determine max/min
            for(i=blockstartNE+threadIdx.x;i<blockendNE;i+=nractivethreads) {
              value=xbase-xpositionsNE[i];
              value2=ybase-ypositionsNE[i];
              value=value*value+value2*value2;
              if(value>maxxvalues[threadIdx.x])
                maxxvalues[threadIdx.x]=value;
            }
          }
        }
        __syncthreads();
        for(i=(nractivethreads>>1);i>0;i>>=1) {
          if(threadIdx.x<i) {
            if(maxxvalues[threadIdx.x+i]>maxxvalues[threadIdx.x])
              maxxvalues[threadIdx.x]=maxxvalues[threadIdx.x+i];
          }
          __syncthreads();
        }
        if(threadIdx.x==0) { //write results
          dabs[index]=sqrt(maxxvalues[0]);
        }
      }
      else
        if(threadIdx.x==0)
          dabs[index]=hypot((maxxvalues[0]-minxvalues[0])*0.5,(maxyvalues[0]-minyvalues[0])*0.5);
      __syncthreads();
      #endif /*RADIALSHRINK*/
    }
    else if(threadIdx.x==0) { //default if no elements
      COMPLEXASSIGN(z[index], 0, 0);
      COMPLEXASSIGN(d[index], 0, 0);
      dabs[index]=0;
    }
  }
}

//this is the second part of the multi-block version of findzd. Has to be implemented as a two pass version, since cuda does not support atomicMax for doubles in the current version
//This function assumes that the variables maxxposition etc. has been set by a previous cuda call (those values contains maximum valus for each block), now combine the results from each block
__global__ void findzdmultistep2(int* llimits,int* rlimits,SORT_REAL *maxxpositions,SORT_REAL *maxypositions,SORT_REAL *minxpositions,SORT_REAL *minypositions,SORT_DCMPLX *z,SORT_DCMPLX *d,SORT_REAL *dabs,int count)
{
  __shared__ SORT_REAL maxxvalues[ZDMAXTHREADS], maxyvalues[ZDMAXTHREADS], minxvalues[ZDMAXTHREADS], minyvalues[ZDMAXTHREADS];
  __shared__ int blockstart, blockend, nractivethreads;
  int i;
  SORT_REAL value;
  for(int index=blockIdx.x;index<count;index+=gridDim.x) {
    if(threadIdx.x==0) {
      blockstart=llimits[index]; //note: this is not the traditional llimits, it is a special one, that contains how to combine the maxxpositions etc. arrays
      blockend=rlimits[index];
      nractivethreads=blockDim.x;
      if(blockend-blockstart>0)
        while(nractivethreads>blockend-blockstart) //atleast one point/thread
          nractivethreads>>=1;
    }
    __syncthreads();
    if(blockend-blockstart>0) {
      if(threadIdx.x<nractivethreads) { //set initial values and make the loop
        maxxvalues[threadIdx.x]=maxxpositions[blockstart+threadIdx.x];
        minxvalues[threadIdx.x]=minxpositions[blockstart+threadIdx.x];
        maxyvalues[threadIdx.x]=maxypositions[blockstart+threadIdx.x];
        minyvalues[threadIdx.x]=minypositions[blockstart+threadIdx.x];
        for(i=blockstart+nractivethreads+threadIdx.x;i<blockend;i+=nractivethreads) {
          value=maxxpositions[i];
          if(value>maxxvalues[threadIdx.x])
            maxxvalues[threadIdx.x]=value;
          value=minxpositions[i];
          if(value<minxvalues[threadIdx.x])
            minxvalues[threadIdx.x]=value;
          value=maxypositions[i];
          if(value>maxyvalues[threadIdx.x])
            maxyvalues[threadIdx.x]=value;
          value=minypositions[i];
          if(value<minyvalues[threadIdx.x])
            minyvalues[threadIdx.x]=value;
        }
      }
      __syncthreads();
      //sum for the block
      for(i=(nractivethreads>>1);i>0;i>>=1) {
        if(threadIdx.x<i) {
          if(maxxvalues[threadIdx.x+i]>maxxvalues[threadIdx.x])
            maxxvalues[threadIdx.x]=maxxvalues[threadIdx.x+i];
          if(minxvalues[threadIdx.x+i]<minxvalues[threadIdx.x])
            minxvalues[threadIdx.x]=minxvalues[threadIdx.x+i];
          if(maxyvalues[threadIdx.x+i]>maxyvalues[threadIdx.x])
            maxyvalues[threadIdx.x]=maxyvalues[threadIdx.x+i];
          if(minyvalues[threadIdx.x+i]<minyvalues[threadIdx.x])
            minyvalues[threadIdx.x]=minyvalues[threadIdx.x+i];
        }
        __syncthreads();
      }
      if(threadIdx.x==0) { //write results
        COMPLEXASSIGN(z[index], (minxvalues[0]+maxxvalues[0])*0.5, (minyvalues[0]+maxyvalues[0])*0.5);
        COMPLEXASSIGN(d[index], (maxxvalues[0]-minxvalues[0])*0.5, (maxyvalues[0]-minyvalues[0])*0.5);
        dabs[index]=hypot((maxxvalues[0]-minxvalues[0])*0.5,(maxyvalues[0]-minyvalues[0])*0.5);
      }
    }
    else { //default
      COMPLEXASSIGN(z[index], 0, 0);
      COMPLEXASSIGN(d[index], 0, 0);
    }
    
  }
}

//first part of the multiblock version of findzd. Each block sums a part of the box, and writes its value to the new vectors maxspositions etc, to be used by step 2 to complete the max/min calculation
__global__ void findzdmulti(int* llimits,int* rlimits,SORT_REAL *xpositions,SORT_REAL *ypositions,int* llimitsNE,int* rlimitsNE,SORT_REAL *xpositionsNE,SORT_REAL *ypositionsNE,int* newllimits,int* newrlimits,SORT_REAL *maxxpositions,SORT_REAL *maxypositions,SORT_REAL *minxpositions,SORT_REAL *minypositions,int count,int nancheck,int* nanoutput DEBUGVECTORSTRING)
{
  __shared__ SORT_REAL maxxvalues[ZDMAXTHREADS], maxyvalues[ZDMAXTHREADS], minxvalues[ZDMAXTHREADS], minyvalues[ZDMAXTHREADS];
  __shared__ int blockstart, blockend, nractivethreads, index, blockcount, blockindex, sumcount, splitblockcount;
  __shared__ int blockstartNE, blockendNE, blockcountNE, splitblockcountNE;
  int i;
  SORT_REAL value;
  
  if(threadIdx.x==0) { //normal setup for the multiblock codes
    index=blockIdx.x*count/gridDim.x;
#ifdef EVENSPLIT //if all splits have the same number of blocks.
    splitblockcount=gridDim.x/count;
    blockindex=blockIdx.x%count;
#else
    blockindex=blockIdx.x*count%gridDim.x;
    splitblockcount=1+blockindex/count+(gridDim.x-blockindex-1)/count; //number of blocks in this split, remember that this is integer divisions, two divisions necessary
    blockindex/=count;
#endif
    blockstart=llimits[index];
    blockend=rlimits[index];
    sumcount=blockend-blockstart;
    splitblockcountNE=splitblockcount;
    if(splitblockcount*blockDim.x>sumcount) { //all threads should have atleast one elements, if number of blocks is too high, reduce it
      splitblockcount=sumcount/blockDim.x;
      if(splitblockcount==0)
        splitblockcount=1;
      if(blockindex>=splitblockcount)
        blockcount=0;
      else
        blockcount=(sumcount+splitblockcount-1)/splitblockcount;
    }
    else
      blockcount=(sumcount+splitblockcount-1)/splitblockcount;
    
    blockstart+=blockcount*blockindex;
    if(blockstart+blockcount>blockend) //correct the limits in each thread
      blockcount=blockend-blockstart;
    blockend=blockstart+blockcount;
    if(xpositionsNE!=NULL) { //if evaluation points, set up these in the same way
      blockstartNE=llimitsNE[index];
      blockendNE=rlimitsNE[index];
      sumcount=blockendNE-blockstartNE;
      if(splitblockcountNE*blockDim.x>sumcount) {
        splitblockcountNE=sumcount/blockDim.x;
        if(splitblockcountNE==0)
          splitblockcountNE=1;
        if(blockindex>=splitblockcountNE)
          blockcountNE=0;
        else
          blockcountNE=(sumcount+splitblockcountNE-1)/splitblockcountNE;
      }
      else
        blockcountNE=(sumcount+splitblockcountNE-1)/splitblockcountNE;
      blockstartNE+=blockcountNE*blockindex;
      if(blockstartNE+blockcountNE>blockendNE)
        blockcountNE=blockendNE-blockstartNE;
      blockendNE=blockstartNE+blockcountNE;
      nractivethreads=blockDim.x;
      sumcount=imax(blockcount, blockcountNE); //let the largest one determine number of threads to use
      if(sumcount>0)
        while(nractivethreads>sumcount)
          nractivethreads>>=1;
      if(blockindex==0) { //set the output limits for step 2 to use
        newllimits[index]=blockIdx.x;
        newrlimits[index]=blockIdx.x+imax(splitblockcount, splitblockcountNE);
      }
    }
    else {
      blockcountNE=0;
      if(blockindex==0) {//set the output limits for step 2 to use
        newllimits[index]=blockIdx.x;
        newrlimits[index]=blockIdx.x+splitblockcount;
      }
      nractivethreads=blockDim.x;
      if(blockcount>0)
        while(nractivethreads>blockcount)
          nractivethreads>>=1;
    }
  }
  __syncthreads();
  if(blockcount>0) { //if there are any elements
    if(threadIdx.x<nractivethreads) {
      if(threadIdx.x<blockcount) { //init values
        maxxvalues[threadIdx.x]=minxvalues[threadIdx.x]=xpositions[blockstart+threadIdx.x];
        maxyvalues[threadIdx.x]=minyvalues[threadIdx.x]=ypositions[blockstart+threadIdx.x];
      }
      
      //loop and determin max/min
      for(i=blockstart+nractivethreads+threadIdx.x;i<blockend;i+=nractivethreads) {
        value=xpositions[i];
        if(value>maxxvalues[threadIdx.x])
          maxxvalues[threadIdx.x]=value;
        if(value<minxvalues[threadIdx.x])
          minxvalues[threadIdx.x]=value;
        if(nancheck&&isnan(value))
          maxxvalues[threadIdx.x]=value;
        value=ypositions[i];
        if(value>maxyvalues[threadIdx.x])
          maxyvalues[threadIdx.x]=value;
        if(value<minyvalues[threadIdx.x])
          minyvalues[threadIdx.x]=value;
        if(nancheck&&isnan(value))
          maxxvalues[threadIdx.x]=value;
        
      }
    }
  }
  if(xpositionsNE!=NULL) { //evaluation points
    if(blockcountNE>0) {
      if(threadIdx.x<nractivethreads) {
        if(threadIdx.x>=blockcount&&threadIdx.x<blockcountNE) { //if no value has been set, set one
          maxxvalues[threadIdx.x]=minxvalues[threadIdx.x]=xpositionsNE[blockstartNE+threadIdx.x];
          maxyvalues[threadIdx.x]=minyvalues[threadIdx.x]=ypositionsNE[blockstartNE+threadIdx.x];
        }
        //loop and determin max/min
        for(i=blockstartNE+threadIdx.x;i<blockendNE;i+=nractivethreads) {
          value=xpositionsNE[i];
          if(value>maxxvalues[threadIdx.x])
            maxxvalues[threadIdx.x]=value;
          if(value<minxvalues[threadIdx.x])
            minxvalues[threadIdx.x]=value;
          if(nancheck&&isnan(value))
            maxxvalues[threadIdx.x]=value;
          value=ypositionsNE[i];
          if(value>maxyvalues[threadIdx.x])
            maxyvalues[threadIdx.x]=value;
          if(value<minyvalues[threadIdx.x])
            minyvalues[threadIdx.x]=value;
          if(nancheck&&isnan(value))
            maxxvalues[threadIdx.x]=value;
        }
      }
    }
  }
  if(blockcount>0||blockcountNE>0) { //if any elements
    __syncthreads();
    //calculate max/min for the block
    for(i=(nractivethreads>>1);i>0;i>>=1) {
      if(threadIdx.x<i) {
        if(nancheck) {
          if(isnan(maxxvalues[threadIdx.x+i]))
            maxxvalues[threadIdx.x]=maxxvalues[threadIdx.x+i];
          else if(!isnan(maxxvalues[threadIdx.x])) {
            if(maxxvalues[threadIdx.x+i]>maxxvalues[threadIdx.x])
              maxxvalues[threadIdx.x]=maxxvalues[threadIdx.x+i];
          }
        }
        else
          if(maxxvalues[threadIdx.x+i]>maxxvalues[threadIdx.x])
            maxxvalues[threadIdx.x]=maxxvalues[threadIdx.x+i];
        if(minxvalues[threadIdx.x+i]<minxvalues[threadIdx.x])
          minxvalues[threadIdx.x]=minxvalues[threadIdx.x+i];
        if(maxyvalues[threadIdx.x+i]>maxyvalues[threadIdx.x])
          maxyvalues[threadIdx.x]=maxyvalues[threadIdx.x+i];
        if(minyvalues[threadIdx.x+i]<minyvalues[threadIdx.x])
          minyvalues[threadIdx.x]=minyvalues[threadIdx.x+i];
        if(nancheck&&isnan(maxxvalues[threadIdx.x+i]))
          maxxvalues[threadIdx.x]=maxxvalues[threadIdx.x+i];
      }
      __syncthreads();
    }
    if(threadIdx.x==0) { //write output
      maxxpositions[blockIdx.x]=maxxvalues[0];
      maxypositions[blockIdx.x]=maxyvalues[0];
      minxpositions[blockIdx.x]=minxvalues[0];
      minypositions[blockIdx.x]=minyvalues[0];
      if(nancheck&&isnan(maxxpositions[blockIdx.x]))
        nanoutput[0]=1;
    }
  }
}

//this function uses old values for z and d, and the splitpoint to approximate the value for the new z and d. Much faster than to calculate it with findzd etc. for large boxes
__global__ void approximatezd(SORT_DCMPLX *oldz,SORT_DCMPLX *oldd,SORT_REAL *splitpoints,SORT_DCMPLX *newz,SORT_DCMPLX *newd,int* ysplit,int count)
{
  SORT_DCMPLX zold, dold, znew1, dnew1, znew2, dnew2; //save as local variables, minimize access to global memory, and make it outside if-statements
  SORT_REAL splitpoint;
  int ysplitlocal;
  for(int i=threadIdx.x+blockIdx.x*blockDim.x;i<count;i+=blockDim.x*gridDim.x) {
    COMPLEXASSIGN(zold, creal(oldz[i]), cimag(oldz[i]));
    COMPLEXASSIGN(dold, creal(oldd[i]), cimag(oldd[i]));
    splitpoint=splitpoints[i];
    ysplitlocal=ysplit[i];
    if(ysplitlocal) { //depending on if last split was a x or y split
      COMPLEXASSIGN(znew1, creal(zold), 0.5*(cimag(zold)-cimag(dold)+splitpoint));
      COMPLEXASSIGN(dnew1, creal(dold), cimag(znew1)-(cimag(zold)-cimag(dold)));
      COMPLEXASSIGN(znew2, creal(zold), 0.5*(cimag(zold)+cimag(dold)+splitpoint));
      COMPLEXASSIGN(dnew2, creal(dold), (cimag(zold)+cimag(dold))-cimag(znew2));
    }
    else {
      COMPLEXASSIGN(znew1, 0.5*(creal(zold)-creal(dold)+splitpoint), cimag(zold));
      COMPLEXASSIGN(dnew1, creal(znew1)-(creal(zold)-creal(dold)), cimag(dold));
      COMPLEXASSIGN(znew2, 0.5*(creal(zold)+creal(dold)+splitpoint), cimag(zold));
      COMPLEXASSIGN(dnew2, (creal(zold)+creal(dold))-creal(znew2), cimag(dold));
    }
    //write results. To keep warp together at the memory access point, it is written to a temporary variable first (speed difference not checked)
    COMPLEXASSIGN(newz[2*i], creal(znew1), cimag(znew1));
    COMPLEXASSIGN(newz[2*i+1], creal(znew2), cimag(znew2));
    COMPLEXASSIGN(newd[2*i], creal(dnew1), cimag(dnew1));
    COMPLEXASSIGN(newd[2*i+1], creal(dnew2), cimag(dnew2));
  }
}

//for the evaluation points where only a split is performed. Since the split functions do not update the limits themselves, compared to singleblockpartition, do it with this function
__global__ void setNElimits(int* llimits,int* rlimits,int* newllimits, int* newrlimits,int* lrcount,int count)
{
  int localllimits, localrlimits, lcount;
  for(int i=threadIdx.x+blockIdx.x*blockDim.x;i<count;i+=blockDim.x*gridDim.x) {
    localllimits=llimits[i];
    localrlimits=rlimits[i];
    lcount=lrcount[2*i];
    newllimits[2*i]=localllimits;
    newllimits[2*i+1]=localllimits+lcount;
    newrlimits[2*i]=localllimits+lcount;
    newrlimits[2*i+1]=localrlimits;
    if(i==count-1)
      newllimits[2*i+2]=localrlimits;
  }
}


//prepare the partitioning. Determine if the split should be performed in the x or y direction
__global__ void setuppartition(SORT_DCMPLX *z,SORT_DCMPLX *d,SORT_REAL* splitpoints,int* ysplit,int prefery,int count SORTLIMITSTRING2)
{
  SORT_DCMPLX zlocal, dlocal;
  int ylocal;
  SORT_REAL splitpointlocal,leftlimitlocal,rightlimitlocal;
  for(int i=threadIdx.x+blockIdx.x*blockDim.x;i<count;i+=blockDim.x*gridDim.x) {
    COMPLEXASSIGN(zlocal, creal(z[i]), cimag(z[i]));
    COMPLEXASSIGN(dlocal, creal(d[i]), cimag(d[i]));
    if(prefery) { //depending on if y-direction is the desired direction or not
      if (creal(dlocal) > meshratio*cimag(dlocal))
        ylocal=0;
      else
        ylocal=1;
    }
    else {
      if (meshratio*creal(dlocal) > cimag(dlocal))
        ylocal=0;
      else
        ylocal=1;
    }
#ifdef SORTLIMIT
    if(ylocal) {
      leftlimitlocal=cimag(zlocal)-cimag(dlocal)*sortlimit;
      rightlimitlocal=cimag(zlocal)+cimag(dlocal)*sortlimit;
    }
    else {
      leftlimitlocal=creal(zlocal)-creal(dlocal)*sortlimit;
      rightlimitlocal=creal(zlocal)+creal(dlocal)*sortlimit;
    }
    leftlimitvalues[i]=leftlimitlocal;
    rightlimitvalues[i]=rightlimitlocal;
//debug
//    if(ylocal) {
//      leftlimitlocal=cimag(zlocal)-cimag(dlocal)*2;
//      rightlimitlocal=cimag(zlocal)+cimag(dlocal)*2;
//    }
//    else {
//      leftlimitlocal=creal(zlocal)-creal(dlocal)*2;
//      rightlimitlocal=creal(zlocal)+creal(dlocal)*2;
//    }
//    leftlimitvaluesdebug[i]=leftlimitlocal;
//    rightlimitvaluesdebug[i]=rightlimitlocal;
//end debug
#endif
    ysplit[i]=ylocal;
    if(ylocal)
      splitpointlocal=cimag(zlocal);
    else
      splitpointlocal=creal(zlocal);
    splitpoints[i]=splitpointlocal; //set the splitpoint for the first split
  }
}

//if eta are to be used. determine which side that requires an additional split
__global__ void setupetasplit(SORT_DCMPLX *z,SORT_REAL* splitpoints,int* ysplit,SORT_REAL eta,int count,int* splitside,int* llimits,int* rlimits,int* newllimits,int* newrlimits)
{
  SORT_REAL oldsplitpointlocal;
  int ylocal, limitindex;
  SORT_REAL splitpointlocal;
  for(int i=threadIdx.x+blockIdx.x*blockDim.x;i<count;i+=blockDim.x*gridDim.x) {
    ylocal=ysplit[i];
    if(ylocal)
      oldsplitpointlocal=cimag(z[i]);
    else
      oldsplitpointlocal=creal(z[i]);
    splitpointlocal=splitpoints[i];
    if(oldsplitpointlocal>splitpointlocal) //which value is largest, the midpoint or the splitpoint? This tells which side that should be splitted
      limitindex=2*i+1;
    else
      limitindex=2*i;
    newllimits[i]=llimits[limitindex];
    newrlimits[i]=rlimits[limitindex];
    splitside[i]=limitindex&1;
    splitpointlocal=splitpointlocal*eta+(1-eta)*oldsplitpointlocal;
    splitpoints[i]=splitpointlocal;
  }
}

//update the limits after the eta split
__global__ void correctetalimits(int *llimits,int* rlimits,int* splitside,int* lrcount,int count)
{
  int splitsidelocal, midpoint;
  for(int i=threadIdx.x+blockIdx.x*blockDim.x;i<count;i+=blockDim.x*gridDim.x) {
    splitsidelocal=splitside[i];
    midpoint=llimits[2*i+1];
    if(splitsidelocal)
      midpoint+=lrcount[2*i];
    else
      midpoint-=lrcount[2*i+1];
    llimits[2*i+1]=rlimits[2*i]=midpoint;
  }
}

//since only half of the elements were moved in the eta split, move these ones back to the original vector
__global__ void splitetacopymultithread(int* llimits,int* rlimits,SORT_REAL* positions,SORT_REAL *positions2,int* indices,SORT_REAL* newpositions,SORT_REAL *newpositions2,int* newindices,int splitcount)
{
  __shared__ int splitindex, splitblockcount, blockindex, splitstart, splitend, count, newsplitblockcount, blockcount, blockstart;
  if(threadIdx.x==0) { //traditional setup
    splitindex=blockIdx.x*splitcount/gridDim.x;
#ifdef EVENSPLIT //if all splits have the same number of blocks.
    splitblockcount=gridDim.x/splitcount;
    blockindex=blockIdx.x%splitcount;
#else
    blockindex=blockIdx.x*splitcount%gridDim.x;
    splitblockcount=1+blockindex/splitcount+(gridDim.x-blockindex-1)/splitcount; //number of blocks in this split
    blockindex/=splitcount;
#endif
    splitstart=llimits[splitindex];
    splitend=rlimits[splitindex];
    count=splitend-splitstart;
    if(splitblockcount*blockDim.x>count) {
      newsplitblockcount=count/blockDim.x;
      if(newsplitblockcount==0)
        newsplitblockcount=1;
      if(blockindex>=newsplitblockcount)
        blockcount=0;
      else
        blockcount=(count+newsplitblockcount-1)/newsplitblockcount;
    }
    else
      blockcount=(count+splitblockcount-1)/splitblockcount;
    blockstart=splitstart+blockcount*blockindex;
    if(blockstart+blockcount>splitend)
      blockcount=splitend-blockstart;
  }
  __syncthreads();
  
  //make the copy
  for(int i=blockstart+threadIdx.x;i<blockstart+blockcount;i+=blockDim.x) {
    newpositions[i]=positions[i];
    newpositions2[i]=positions2[i];
    newindices[i]=indices[i];
  }
}

//same as above, but only one block per copy
__global__ void splitetacopysinglethread(int* llimits,int* rlimits,SORT_REAL* positions,SORT_REAL *positions2,int* indices,SORT_REAL* newpositions,SORT_REAL *newpositions2,int* newindices,int splitcount)
{
  __shared__ int splitstart, splitend, count;
  for(int i=blockIdx.x;i<splitcount;i+=gridDim.x) {
    if(threadIdx.x==0) {
      splitstart=llimits[i];
      splitend=rlimits[i];
      count=splitend-splitstart;
    }
    __syncthreads();
    for(int i=splitstart+threadIdx.x;i<splitstart+count;i+=blockDim.x) {
      newpositions[i]=positions[i];
      newpositions2[i]=positions2[i];
      newindices[i]=indices[i];
    }
  }
}
#endif /*defined(CUDASUPPORT) && defined(CUDASORT)*/

#ifdef CUDASUPPORT

//the cuda verion of the theta test. Should be almost identical to the CPU version, except the input
__device__ int mpexp_theta_cuda(SORT_DCMPLX z0,SORT_REAL Thisrad,SORT_DCMPLX z1,SORT_REAL thatrad,SORT_REAL cutoff)
{
    
    /* note: criterion uses in fact a strict inequality to get rid of
     the singular case when both boxes are the same points (This can
     happen!) */
  SORT_REAL d = sqrt((creal(z0)-creal(z1))*(creal(z0)-creal(z1)) + (cimag(z0)-cimag(z1))*(cimag(z0)-cimag(z1)));
  if(d-Thisrad-thatrad<cutoff)
    return false;
  if (Thisrad >= thatrad)
    return Thisrad+theta*thatrad < theta*d;
  else
    return thatrad+theta*Thisrad < theta*d;
}

//cuda version of second theta test. Should also almost be identical to CPU version
__device__ int mpexp_theta2_cuda(SORT_DCMPLX z0,SORT_REAL Thisrad,SORT_DCMPLX z1,SORT_REAL thatrad,SORT_REAL cutoff)
/* 2nd theta criterion: if mpexp_theta() above is false, is it due to
 * only one of the boxes? Then the multipoles of the larger box can be
 * shifted into a polynomial expansion in the smaller, and the
 * multipole expansion in the smaller can be evaluated at each point
 * in the larger box. Returns -1 if This is the smaller box (cluster
 * to particle) and 1 otherwise (particle to cluster). */
{
  SORT_REAL d = sqrt((creal(z0)-creal(z1))*(creal(z0)-creal(z1)) + (cimag(z0)-cimag(z1))*(cimag(z0)-cimag(z1)));
  if(d-Thisrad-thatrad<cutoff)
    return 0;
  if (Thisrad >= thatrad)
    return thatrad+theta*Thisrad < theta*d ? -1 : 0;
  else
    return Thisrad+theta*thatrad < theta*d ? 1 : 0;
}

/*-----------------------------------------------------------------------*/
__global__ void cudacreateconnectivity(int *jcptr,int *kcptr,int *ir,
                                       int *oldjcptr,int *oldkcptr,int *oldir,
                                       int count,int maxm2p,
                                       SORT_DCMPLX *z,SORT_REAL *dabs,
                                       SORT_REAL cutoff,int lastlevel,
                                       int *outputvector DEBUGVECTORSTRING)
// Create connectivity matrices. Basically the same as the CPU version.
{
  __shared__ int nc[MAXCONNECTIVITYTHREADS],maxm2plocal[MAXCONNECTIVITYTHREADS];
  int fwd,bwd,i,k,d,dir1,dir2;

  for (int j = threadIdx.x+blockIdx.x*blockDim.x; j < count;
      j += blockDim.x*gridDim.x) {
    maxm2plocal[threadIdx.x] = maxm2p;
    fwd = dir2 = jcptr[j];
    bwd = dir1 = jcptr[j+1];
    
    if(fwd==bwd) { //special case if no evaluation points in box
      nc[threadIdx.x]=0;
      maxm2plocal[threadIdx.x]=0;
      kcptr[j] = fwd;
      if(lastlevel)
        kcptr[j+count] = bwd;
    }
    else {

      // parents (near) connections (includes the parent itself)
      for (i = oldjcptr[j >> 2]; i < oldkcptr[j >> 2]; i++) {
        /* the children of the near connections of box i become either
           near connections or far connections of box j */
        k = oldir[i] << 2; // 0-child of connection
        for (d = 0; d < 4; d++, k++) {
          if (!mpexp_theta_cuda(z[j],dabs[j],z[k],dabs[k],cutoff))
            ir[fwd++] = k; // near connection
          else
            ir[--bwd] = k; // far connection
        }
      }
      kcptr[j] = bwd;
      /* this splits the column into near connections
         [jcptr[j],kcptr[j]) and far connections
         [kcptr[j],jcptr[j+1]) */
      nc[threadIdx.x] = bwd-jcptr[j];
      if (dir1-bwd > maxm2plocal[threadIdx.x])
        maxm2plocal[threadIdx.x] = dir1-bwd;

      // last level: sort the near connections one step further
      if (lastlevel) {
        fwd = dir2-1;
        // since the format uses the sign of the index to denote the
        // type of the near connection, box 0 must be skipped
        if (j == 0)
          bwd--; // skip box 0
        else {
          // again, skip box 0 (must be first since column is sorted)
          if (ir[fwd+1] == 0)
            fwd++;
          for ( ; ; ) {
            do
              fwd++;
            while (fwd < bwd &&
                   !(dir1 = mpexp_theta2_cuda(z[j],dabs[j],
                                              z[ir[fwd]],dabs[ir[fwd]],
                                              cutoff)));
            do {
              bwd--;
              if (fwd <= bwd &&
                  (dir2 = mpexp_theta2_cuda(z[j],dabs[j],
                                            z[ir[bwd]],dabs[ir[bwd]],cutoff)))
                ir[bwd] *= dir2;
              else
                break;
            } while (1);

            if (fwd >= bwd) break;
            d = ir[fwd];
            ir[fwd] = ir[bwd];
            ir[bwd] = d*dir1;
          }
        }
        kcptr[j+count] = bwd+1;
        /* this splits the near connections into 'truly near'
           [jcptr[j],kcptr[j+N]) and 'less near' [kcptr[j+N],kcptr[j])
           connections */
      }
    }
    __syncthreads();
 
    // calculate the block sum and max and save these to the outputvector
    for (d = blockDim.x >> 1; d > 0; d >>= 1) {
      if (threadIdx.x < d && d < count-j) {
        if (maxm2plocal[threadIdx.x+d] > maxm2plocal[threadIdx.x])
          maxm2plocal[threadIdx.x] = maxm2plocal[threadIdx.x+d];
        nc[threadIdx.x] += nc[threadIdx.x+d];
      }
      __syncthreads();
    }
    if (threadIdx.x == 0) {
      atomicMax(outputvector,maxm2plocal[0]);
      atomicAdd(outputvector+1,nc[0]);
    }
  }
}
/*-----------------------------------------------------------------------*/
//calculates the absolute value of d for all boxes. This is to avoid calculating it several times for each box in the theta criterion
__global__ void calculatedabs(const SORT_DCMPLX *d,SORT_REAL *dabs,int count)
{
  for(int j = threadIdx.x+blockIdx.x*blockDim.x; j < count;
      j += blockDim.x*gridDim.x) {
#ifdef FLOATSORT
    dabs[j] = hypotf(creal(d[j]),cimag(d[j]));
#else
    dabs[j] = hypot(creal(d[j]),cimag(d[j]));
#endif
  }
}
/*-----------------------------------------------------------------------*/
//converts position data to float for FLOATSORT
__global__ void converttofloat(float *out1,float *out2,const double* in1,const double *in2,int N)
{
  for(int j = threadIdx.x+blockIdx.x*blockDim.x; j < N;
      j += blockDim.x*gridDim.x) {
    out1[j]=(float)in1[j];
    out2[j]=(float)in2[j];
  }
}
/*------------------------------------------------------------------------*/
__global__ void cumsuminitevalonly(int *oldjcptr,int* oldkcptr, int* jcptr, size_t count,size_t blockcount,int* jxptr,int evalshift)
{
  __shared__ double localarray[threadsperblock];
  for(size_t blockindex=blockIdx.x;blockindex<blockcount;blockindex+=gridDim.x) {
    size_t index=blockindex*blockDim.x+threadIdx.x;
    
    //calculate starting vector
    if(index<count) {
      if(jxptr[(index+1)<<evalshift]==jxptr[index<<evalshift]) //no evaluation points, no interactions necessary
        localarray[threadIdx.x]=0;
      else
        localarray[threadIdx.x]=4*(oldkcptr[index>>2]-oldjcptr[index>>2]);
    }
    else
      localarray[threadIdx.x]=0;
    __syncthreads();
    
    //calculate cumulative sum within block
    int i=1,a=1;
    for(int c=blockDim.x>>1;c>=1;c>>=1, i<<=1, a++) {
      if(threadIdx.x<c)
        localarray[(threadIdx.x<<a)+(i<<1)-1]+=localarray[(threadIdx.x<<a)+i-1];
      __syncthreads();
    }
    i=(blockDim.x>>2);
    a-=2;
    for(int c=2;i>=1;c<<=1, i>>=1, a--) {
      if(threadIdx.x<c-1)
        localarray[((threadIdx.x+1)<<a)+i-1]+=localarray[((threadIdx.x+1)<<a)-1];
      __syncthreads();
    }
    if(index<count)
      jcptr[index]=localarray[threadIdx.x];
  }
}
/*------------------------------------------------------------------------*/
__global__ void cumsuminit(int *oldjcptr,int* oldkcptr, int* jcptr, size_t count,size_t blockcount,int* ixptr,int* jxptr,int evalshift)
{
  __shared__ double localarray[threadsperblock];
  for(size_t blockindex=blockIdx.x;blockindex<blockcount;blockindex+=gridDim.x) {
    size_t index=blockindex*blockDim.x+threadIdx.x;
    
    //calculate starting vector
    if(index<count)
      if(ixptr[(index+1)<<evalshift]==ixptr[index<<evalshift]) { //if box is empty, no connections is necessary
        if(jxptr!=NULL&&jxptr[(index+1)<<evalshift]!=jxptr[index<<evalshift]) //if it has evaluation points, it is not empty
          localarray[threadIdx.x]=4*(oldkcptr[index>>2]-oldjcptr[index>>2]);
        else
          localarray[threadIdx.x]=0;
      }
      else
        localarray[threadIdx.x]=4*(oldkcptr[index>>2]-oldjcptr[index>>2]);
    else
      localarray[threadIdx.x]=0;
    __syncthreads();
    
    //calculate cumulative sum within block
    int i=1,a=1;
    for(int c=blockDim.x>>1;c>=1;c>>=1, i<<=1, a++) {
      if(threadIdx.x<c)
        localarray[(threadIdx.x<<a)+(i<<1)-1]+=localarray[(threadIdx.x<<a)+i-1];
      __syncthreads();
    }
    i=(blockDim.x>>2);
    a-=2;
    for(int c=2;i>=1;c<<=1, i>>=1, a--) {
      if(threadIdx.x<c-1)
        localarray[((threadIdx.x+1)<<a)+i-1]+=localarray[((threadIdx.x+1)<<a)-1];
      __syncthreads();
    }
    if(index<count)
      jcptr[index]=localarray[threadIdx.x];
  }
}
/*------------------------------------------------------------------------*/
__global__ void cumsumpass1(int *carray, size_t count,size_t blockcount,size_t shift)
{
  __shared__ int localarray[threadsperblock];
  for(size_t blockindex=blockIdx.x;blockindex<blockcount;blockindex+=gridDim.x) {
    size_t index=(blockindex*blockDim.x+threadIdx.x)<<shift;
    
    //copy data to shared memory
    if(index<count)
      localarray[threadIdx.x]=carray[index];
    else
      localarray[threadIdx.x]=0;
    __syncthreads();
    
    //calculate cumulative sum within block
    int i=1,a=1;
    for(int c=blockDim.x>>1;c>=1;c>>=1, i<<=1, a++) {
      if(threadIdx.x<c)
        localarray[(threadIdx.x<<a)+(i<<1)-1]+=localarray[(threadIdx.x<<a)+i-1];
      __syncthreads();
    }
    i=(blockDim.x>>2);
    a-=2;
    for(int c=2;i>=1;c<<=1, i>>=1, a--) {
      if(threadIdx.x<c-1)
        localarray[((threadIdx.x+1)<<a)+i-1]+=localarray[((threadIdx.x+1)<<a)-1];
      __syncthreads();
    }
    if(index<count)
      carray[index]=localarray[threadIdx.x];
  }
}
/*------------------------------------------------------------------------*/
__global__ void cumsumpass2(int *carray, size_t count,size_t blockcount,size_t shift/*,double *debugvector*/)
{
  for(size_t blockindex=blockIdx.x+1;blockindex<blockcount;blockindex+=gridDim.x) {
    size_t index=(blockindex*blockDim.x+threadIdx.x)<<shift;
    if(index<count&&threadIdx.x!=blockDim.x-1)
      carray[index]+=carray[(blockindex*blockDim.x-1)<<shift];
//     if(shift==8&&index<count)
//       debugvector[index]=carray[(blockindex*blockDim.x-1)<<shift];
  }
}
/*-----------------------------------------------------------------------*/
#endif /* CUDASUPPORT */