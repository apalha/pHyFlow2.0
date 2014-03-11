
#ifndef UNROLLLEVEL
#define UNROLLLEVEL 7
#endif
#include "hornershiftdefs.h"
#include "expint.h"

//kernels for direct evaluation. These will be used several times, therefore use define to avoid repeating code
#define DIRACSYNCKERNEL(xvelocity,yvelocity,circ,icirc) sqrdist=xdist*xdist+ydist*ydist;\
if(sqrdist!=0.0) {\
  sqrdist=circ*(1.0/(sqrdist));\
  xvelocity+=sqrdist*xdist;\
  yvelocity-=sqrdist*ydist;\
}
#define IDIRACSYNCKERNEL(xvelocity,yvelocity,circ,icirc) sqrdist=xdist*xdist+ydist*ydist;\
if(sqrdist!=0.0) {\
  sqrdist=1/(sqrdist);\
  xvelocity+=sqrdist*(circ*xdist+icirc*ydist);\
  yvelocity-=sqrdist*(circ*ydist-icirc*xdist);\
}
#define RANKINESYNCKERNEL(xvelocity,yvelocity,circ,icirc) sqrdist=1/(xdist*xdist+ydist*ydist);\
if (sqrdist<shape){/*since both contains the 1/X, the comparison has to change sign*/\
  xvelocity+=circ*sqrdist*xdist;\
  yvelocity-=circ*sqrdist*ydist;\
}\
else if(sqrdist!=0){\
  xvelocity+=circ*shape*xdist;\
  yvelocity-=circ*shape*ydist;\
}
#define IRANKINESYNCKERNEL(xvelocity,yvelocity,circ,icirc) sqrdist=1/(xdist*xdist+ydist*ydist);\
if (sqrdist<shape){/*since both contains the 1/X, the comparison has to change sign*/\
  xvelocity+=sqrdist*(circ*xdist+icirc*ydist);\
  yvelocity-=sqrdist*(circ*ydist-icirc*xdist);\
}\
else if(sqrdist!=0){\
  xvelocity+=shape*(circ*xdist+icirc*ydist);\
  yvelocity-=shape*(circ*ydist-icirc*xdist);\
}
#define SCULLYSYNCKERNEL(xvelocity,yvelocity,circ,icirc) sqrdist=xdist*xdist+ydist*ydist;\
if (sqrdist>cutoffsqr) {/*since both contains the 1/X, the comparison has to change sign**/\
  sqrdist=circ/sqrdist;\
  xvelocity+=sqrdist*xdist;\
  yvelocity-=sqrdist*ydist;\
}\
else if(sqrdist!=0) {\
  sqrdist=circ/(sqrdist+shape)*scale;\
  xvelocity+=sqrdist*xdist;\
  yvelocity-=sqrdist*ydist;\
}
#define ISCULLYSYNCKERNEL(xvelocity,yvelocity,circ,icirc) sqrdist=xdist*xdist+ydist*ydist;\
if (sqrdist>cutoffsqr) {/*since both contains the 1/X, the comparison has to change sign**/\
  sqrdist=1/sqrdist;\
  xvelocity+=sqrdist*(circ*xdist+icirc*ydist);\
  yvelocity-=sqrdist*(circ*ydist-icirc*xdist);\
}\
else if(sqrdist!=0) {\
  sqrdist=1/(sqrdist+shape)*scale;\
  xvelocity+=sqrdist*(circ*xdist+icirc*ydist);\
  yvelocity-=sqrdist*(circ*ydist-icirc*xdist);\
}
#define OSEENSYNCKERNEL(xvelocity,yvelocity,circ,icirc) sqrdist=xdist*xdist+ydist*ydist;\
if (sqrdist>cutoffsqr) {/*since both contains the 1/X, the comparison has to change sign*/\
  sqrdist=circ/sqrdist;\
  xvelocity+=sqrdist*xdist;\
  yvelocity-=sqrdist*ydist;\
}\
else if(sqrdist!=0) {  /*can still fail if cutoffsqr==0*/\
  sqrdist=-circ/sqrdist*scale*expm1(-shape*sqrdist);\
  xvelocity+=sqrdist*xdist;\
  yvelocity-=sqrdist*ydist;\
}//since two version should exist, one that takes data from evaluation points and one from potential points, make a prototype
#define IOSEENSYNCKERNEL(xvelocity,yvelocity,circ,icirc) sqrdist=xdist*xdist+ydist*ydist;\
if (sqrdist>cutoffsqr) {/*since both contains the 1/X, the comparison has to change sign*/\
  sqrdist=1/sqrdist;\
  xvelocity+=sqrdist*(circ*xdist+icirc*ydist);\
  yvelocity-=sqrdist*(circ*ydist-icirc*xdist);\
}\
else if(sqrdist!=0) {  /*can still fail if cutoffsqr==0*/\
  sqrdist=-1/sqrdist*scale*expm1(-shape*sqrdist);\
  xvelocity+=sqrdist*(circ*xdist+icirc*ydist);\
  yvelocity-=sqrdist*(circ*ydist-icirc*xdist);\
}
#define DIRACLOGSYNCKERNEL(xvelocity,yvelocity,circ,icirc) sqrdist=xdist*xdist+ydist*ydist;\
if(sqrdist!=0.0) {\
  xvelocity-=0.5*circ*log(sqrdist);\
}
#define IDIRACLOGSYNCKERNEL(xvelocity,yvelocity,circ,icirc) sqrdist=xdist*xdist+ydist*ydist;\
if(sqrdist!=0.0) {\
  sqrdist=0.5*log(sqrdist);\
  xvelocity-=circ*sqrdist;\
  yvelocity-=icirc*sqrdist;\
}
#define RANKINELOGSYNCKERNEL(xvelocity,yvelocity,circ,icirc) sqrdist=xdist*xdist+ydist*ydist;\
if (sqrdist >= cutoffsqr)\
  sqrdist = 0.5*log(sqrdist);\
else\
  sqrdist = 0.5*(sqrdist/cutoffsqr)+shape;\
xvelocity-=circ*sqrdist;
#define IRANKINELOGSYNCKERNEL(xvelocity,yvelocity,circ,icirc) sqrdist=xdist*xdist+ydist*ydist;\
if (sqrdist >= cutoffsqr)\
  sqrdist = 0.5*log(sqrdist);\
else\
  sqrdist = 0.5*(sqrdist/cutoffsqr)+shape;\
xvelocity-=circ*sqrdist;\
yvelocity-=icirc*sqrdist;
#define SCULLYLOGSYNCKERNEL(xvelocity,yvelocity,circ,icirc) sqrdist=xdist*xdist+ydist*ydist;\
if (sqrdist >= cutoffsqr)\
  sqrdist = 0.5*log(sqrdist);\
else\
  sqrdist = 0.5*log(shape+sqrdist)*scale;\
xvelocity-=circ*sqrdist;
#define ISCULLYLOGSYNCKERNEL(xvelocity,yvelocity,circ,icirc) sqrdist=xdist*xdist+ydist*ydist;\
if (sqrdist >= cutoffsqr)\
  sqrdist = 0.5*log(sqrdist);\
else\
  sqrdist = 0.5*log(shape+sqrdist)*scale;\
xvelocity-=circ*sqrdist;\
yvelocity-=icirc*sqrdist;
#define OSEENLOGSYNCKERNEL(xvelocity,yvelocity,circ,icirc) sqrdist=xdist*xdist+ydist*ydist;\
if (sqrdist >= cutoffsqr)\
  sqrdist = 0.5*log(sqrdist);\
else if(sqrdist!=0.0)\
  sqrdist = (0.5*log(sqrdist)+0.5*expint(shape*sqrdist))*scale;\
else\
  sqrdist = -0.5*(log(shape)+gamma_)*scale;\
xvelocity-=circ*sqrdist;
#define IOSEENLOGSYNCKERNEL(xvelocity,yvelocity,circ,icirc) sqrdist=xdist*xdist+ydist*ydist;\
if (sqrdist >= cutoffsqr)\
  sqrdist = 0.5*log(sqrdist);\
else if(sqrdist!=0.0)\
  sqrdist = (0.5*log(sqrdist)+0.5*expint(shape*sqrdist))*scale;\
else\
  sqrdist = -0.5*(log(shape)+gamma_)*scale;\
xvelocity-=circ*sqrdist;\
yvelocity-=icirc*sqrdist;

//depending on complex mass, complex answer etc, some changes has to be made to the kernel code. These strings are the differences
#define DUMMYKERNELSTRING
#define IKERNELSTRING1 __shared__ double igammav[MAXTHREADSEVALUATION];
// #define IKERNELSTRING1DIRECT __shared__ double igammav[DIRECTMAXTHREADS];
#define IKERNELSTRING1DIRECT __shared__ double *igammav;
#define IKERNELSTRING1DIRECT2 igammav=gammav+blockDim.x;
#define IKERNELSTRING2 TEXTUREFETCH(igammav[threadIdx.x],tmi,k+threadIdx.x);
#define COMPLEXANSWER1 double yvelocity;
#define COMPLEXANSWER2 yvelocity=0;
#define COMPLEXANSWER3 yvelocities[i]+=yvelocity;
#define COMPLEXANSWER3DIRECT yvelocities[i]=yvelocity;

#define IKERNELSTRING2_2 TEXTUREFETCH(igammav[threadIdx.x+addstart],tmi,k+threadIdx.x);
#define COMPLEXANSWER1_2 __shared__ double yvelocity[MAXTHREADSEVALUATION],yvelocitytmp[MAXTHREADSEVALUATION];
#define COMPLEXANSWER2_2 yvelocity[threadIdx.x]=0;
#define COMPLEXANSWER3_2 yvelocities[i]+=yvelocity[threadIdx.x];
#define COMPLEXANSWER4_2 yvelocitytmp[threadIdx.x]=0;
#define COMPLEXANSWER5_2 yvelocity[threadIdx.x]+=yvelocitytmp[threadIdx.x<<1]+yvelocitytmp[(threadIdx.x<<1)+1];
#define COMPLEXANSWER6_2 yvelocity[threadIdx.x+ilocal]+=yvelocitytmp[threadIdx.x<<2]+yvelocitytmp[(threadIdx.x<<2)+1]+yvelocitytmp[(threadIdx.x<<2)+2]+yvelocitytmp[(threadIdx.x<<2)+3];
//kernel code prototype
#define SYNCNRONIZEDKERNEL(kernel,tevalr,tevali,txptr,ISTRING1,ISTRING2,CANSWER1,CANSWER2,CANSWER3,CANSWER4,CANSWER5,CANSWER6) int i,j,k,boxnr,m;\
  double sqrdist,xdist,ydist,xvelocity;\
  CANSWER1\
  __shared__ double vxpositions[MAXTHREADSEVALUATION];\
  __shared__ double vypositions[MAXTHREADSEVALUATION];\
  __shared__ double gammav[MAXTHREADSEVALUATION];\
  ISTRING1\
  __shared__ double xpositions[MAXTHREADSEVALUATION];/*makes code approx 1% slower for heavy load, but for some reason, doubles has to be fetched into an array*/\
  __shared__ double ypositions[MAXTHREADSEVALUATION];\
  __shared__ int count,loopmax,vcount,begin,end,loopstart,mstart,mend;\
  for(boxnr=blockIdx.x;boxnr<Nf;boxnr+=gridDim.x) {/*fist loop is over the boxes, one block take one box*/\
    if(threadIdx.x==0) {\
      count=tex1Dfetch(txptr,boxnr+1); /*right now, it is possible that two boxes has same evalbegin if one box has 0 elements*/\
      loopstart=tex1Dfetch(txptr,boxnr);\
      loopmax=count-(count+blockDim.x-1-(loopstart%blockDim.x))%blockDim.x-1+blockDim.x;/*This one makes sure all threads loop the same number of times, necessary for __syncthreads()*/\
      mstart=jcptr[boxnr];\
      mend=j2cptr[boxnr];\
    }\
    __syncthreads();\
    for(i=threadIdx.x+loopstart;i<loopmax;i+=blockDim.x) { /*if more points than threads*/\
      if(i<count) {/*Always check if the point belongs to the box*/\
        xvelocity=0; /*use local variables as much as possible to help the optimizer understand that they are constant*/\
        /*yvelocity=0;*/CANSWER2\
        TEXTUREFETCH(xpositions[threadIdx.x],tevalr,i);\
        TEXTUREFETCH(ypositions[threadIdx.x],tevali,i);\
      }\
      for(m=mstart;m<mend;m++) {/*Loop over all other boxes*/\
        if(threadIdx.x==0) {/*Set up the loop limits*/\
          begin=tex1Dfetch(tixptr,/*tex1Dfetch(tpotbegin, m)*/ir[m]);\
          end=tex1Dfetch(tixptr,/*tex1Dfetch(tpotbegin, m)*/ir[m]+1);\
          vcount=(end-begin)%blockDim.x;\
        }\
        __syncthreads();\
        for(k=begin;k+blockDim.x<=end;k+=blockDim.x) {/*Loop over the elements in the external box, this step is for when there are more than blockDim.x elements left*/\
          TEXTUREFETCH(vxpositions[threadIdx.x],tzr,k+threadIdx.x);\
          TEXTUREFETCH(vypositions[threadIdx.x],tzi,k+threadIdx.x);\
          TEXTUREFETCH(gammav[threadIdx.x],tmr,k+threadIdx.x);\
          ISTRING2\
          __syncthreads();\
          if(i<count) {\
            for(j=0;j<blockDim.x;j++) {/*the interaction loop*/\
              xdist=vxpositions[j]-xpositions[threadIdx.x];\
              ydist=vypositions[j]-ypositions[threadIdx.x];\
              kernel(xvelocity,yvelocity,circ,icirc);\
            }\
          }\
          __syncthreads();\
        }\
        if(k+threadIdx.x<end) {/*The remaining points*/\
          TEXTUREFETCH(vxpositions[threadIdx.x],tzr,k+threadIdx.x);\
          TEXTUREFETCH(vypositions[threadIdx.x],tzi,k+threadIdx.x);\
          TEXTUREFETCH(gammav[threadIdx.x],tmr,k+threadIdx.x);\
          ISTRING2\
        }\
        __syncthreads();\
        if(i<count) {\
          for(j=0;j<vcount;j++) {/*the interaction loop*/\
            xdist=vxpositions[j]-xpositions[threadIdx.x];\
            ydist=vypositions[j]-ypositions[threadIdx.x];\
            kernel(xvelocity,yvelocity,circ,icirc);\
          }\
        }\
        __syncthreads();\
      }\
      if(i<count) {/*write results*/\
        xvelocities[i]+=xvelocity;\
        /*yvelocities[i]+=yvelocity*/CANSWER3\
      }\
    }\
    __syncthreads();\
  }

//#define BOXLOOPSIZE 32 //cache this many boxes indices. Too high number results in more memory use. Do not use higher number than number of threads
#define SOURCEMULTIPLIER 1 //most likely more source points than evaluation points, cache more source points to increase efficiency by reducing loop overhead

/*This code caches a full line of points before any kernel evaluation,
 *and if number of evaluation points is smaller than number of threads
 *the code will use more threads for each point. At the moment, a maximum of
 *4 threads for each point is used, since it was observed that the evaluation
 *of high values of ndirect slightly increased when increasing the number of splits*/
#define SYNCNRONIZEDKERNEL2(kernel,tevalr,tevali,txptr,ISTRING1,ISTRING2,CANSWER1,CANSWER2,CANSWER3,CANSWER4,CANSWER5,CANSWER6) int i,j,k,boxnr,m,addstart,ilocal;\
  double sqrdist,xdist,ydist;\
  __shared__ double xvelocity[MAXTHREADSEVALUATION];\
  CANSWER1\
  __shared__ double xvelocitytmp[MAXTHREADSEVALUATION];\
  __shared__ double vxpositions[MAXTHREADSEVALUATION*SOURCEMULTIPLIER];\
  __shared__ double vypositions[MAXTHREADSEVALUATION*SOURCEMULTIPLIER];\
  __shared__ double gammav[MAXTHREADSEVALUATION*SOURCEMULTIPLIER];\
  ISTRING1\
  __shared__ double xpositions[MAXTHREADSEVALUATION];/*makes code approx 1% slower for heavy load, but for some reason, doubles has to be fetched into an array*/\
  __shared__ double ypositions[MAXTHREADSEVALUATION];\
  __shared__ int count,loopmax,begin,end,loopstart,mstart,mend;\
  for(boxnr=blockIdx.x;boxnr<Nf;boxnr+=gridDim.x) {\
    if(threadIdx.x==0) {\
      count=tex1Dfetch(txptr,boxnr+1); /*right now, it is possible that two boxes has same evalbegin if one box has 0 elements*/\
      loopstart=tex1Dfetch(txptr,boxnr);\
      loopmax=count-(count+blockDim.x-1-(loopstart%blockDim.x))%blockDim.x-1+blockDim.x;\
      mstart=jcptr[boxnr];\
      mend=j2cptr[boxnr];\
    }\
    __syncthreads();\
    for(i=threadIdx.x+loopstart;i<loopmax;i+=blockDim.x) { /*if more points than threads*/\
      if(i<count) {\
        xvelocity[threadIdx.x]=0; /*use local variables as much as possible to help the optimizer understand that they are constant*/\
        CANSWER2\
        TEXTUREFETCH(xpositions[threadIdx.x],tevalr,i);\
        TEXTUREFETCH(ypositions[threadIdx.x],tevali,i);\
      }\
      addstart=0;\
      for(m=mstart;m<mend;m++) {/*Loop through all interaction boxes*/\
        __syncthreads();\
        if(threadIdx.x==0) {/*cache start and end position*/\
          k=ir[m];\
          begin=tex1Dfetch(tixptr,k);\
          end=tex1Dfetch(tixptr,k+1);\
        }\
        __syncthreads();\
        /*if(end-begin<MAXTHREADSEVALUATION*SOURCEMULTIPLIER-addstart&&tex1Dfetch(tboxstart, boxnr+1)!=m+1) {/*whole box can be cached. This code could have been faster, but did not give any improvements*/\
          /*for(k=begin;k+blockDim.x<=end;k+=blockDim.x,addstart+=blockDim.x) {*/\
             /*TEXTUREFETCH(vxpositions[addstart+threadIdx.x],tzr,k+threadIdx.x);*/\
             /*TEXTUREFETCH(vypositions[addstart+threadIdx.x],tzi,k+threadIdx.x);*/\
             /*TEXTUREFETCH(gammav[addstart+threadIdx.x],tmr,k+threadIdx.x);*/\
             /*ISTRING2 /*Fix this one*/\
          /*}*/\
          /*if(k+threadIdx.x<end) {*/\
            /*TEXTUREFETCH(vxpositions[addstart+threadIdx.x],tzr,k+threadIdx.x);*/\
             /*TEXTUREFETCH(vypositions[addstart+threadIdx.x],tzi,k+threadIdx.x);*/\
             /*TEXTUREFETCH(gammav[addstart+threadIdx.x],tmr,k+threadIdx.x);*/\
             /*ISTRING2 /*Fix this one*/\
          /*}*/\
          /*addstart+=end-k;*/\
          /*continue;*/\
        /*}*/\
        for(k=begin;k<=end/*&&addstart<MAXTHREADSEVALUATION*SOURCEMULTIPLIER*/;k+=blockDim.x) {/*Start and fill the array. End condition is when box is done. Not that code has to enter the loop even if no points, because it may be the last loop, and summation is inside loop*/\
          if(end-k>=MAXTHREADSEVALUATION*SOURCEMULTIPLIER-addstart&&MAXTHREADSEVALUATION*SOURCEMULTIPLIER-addstart<=blockDim.x) {/*full array*/\
            if(addstart+threadIdx.x<MAXTHREADSEVALUATION*SOURCEMULTIPLIER) {/*If it will fill the whole array, fill it and perform the interaction*/\
              TEXTUREFETCH(vxpositions[addstart+threadIdx.x],tzr,k+threadIdx.x);\
              TEXTUREFETCH(vypositions[addstart+threadIdx.x],tzi,k+threadIdx.x);\
              TEXTUREFETCH(gammav[addstart+threadIdx.x],tmr,k+threadIdx.x);\
              ISTRING2\
            }\
            __syncthreads();\
            ilocal=0;/*remembers how many interactions that has been doen*/\
            if(count-(i-threadIdx.x)>=blockDim.x) {/*if more calculation points than threads, one thread per evaluation point*/\
              for(j=0;j<MAXTHREADSEVALUATION*SOURCEMULTIPLIER;j++) {\
                xdist=vxpositions[j]-xpositions[threadIdx.x];\
                ydist=vypositions[j]-ypositions[threadIdx.x];\
                kernel(xvelocity[threadIdx.x],yvelocity[threadIdx.x],gammav[j],igammav[j]);\
              }\
            }\
            else { /*else, use several threads for each point to get better balance*/\
              if(count-(i-threadIdx.x)>=(blockDim.x>>1)){ /*more than half of the points active?*/\
                xvelocitytmp[threadIdx.x]=0;\
                CANSWER4/*yvelocitytmp[threadIdx.x]=0*/;\
                for(j=0;j<MAXTHREADSEVALUATION*SOURCEMULTIPLIER;j+=2) {/*calculate temporary points*/\
                  xdist=vxpositions[j+(threadIdx.x&1)]-xpositions[(threadIdx.x>>1)];\
                  ydist=vypositions[j+(threadIdx.x&1)]-ypositions[(threadIdx.x>>1)];\
                  kernel(xvelocitytmp[threadIdx.x],yvelocitytmp[threadIdx.x],gammav[j+(threadIdx.x&1)],igammav[j+(threadIdx.x&1)]);\
                }\
                __syncthreads();\
                if(threadIdx.x<(blockDim.x>>1)) {/*sum, note that ilocal only can be zero here*/\
                  xvelocity[threadIdx.x]+=xvelocitytmp[threadIdx.x<<1]+xvelocitytmp[(threadIdx.x<<1)+1];\
                  CANSWER5/*yvelocity[threadIdx.x]+=yvelocitytmp[threadIdx.x<<1]+yvelocitytmp[(threadIdx.x<<1)+1];*/\
                }\
                __syncthreads();\
                ilocal+=(blockDim.x>>1);\
              }\
              if(count-(i-threadIdx.x)-ilocal>=(blockDim.x>>2)){/*if remaining points are more than 1/4 of the number of threads*/\
                xvelocitytmp[threadIdx.x]=0;\
                CANSWER4/*yvelocitytmp[threadIdx.x]=0*/;\
                for(j=0;j<MAXTHREADSEVALUATION*SOURCEMULTIPLIER;j+=4) {/*calculate temporary points*/\
                  xdist=vxpositions[j+(threadIdx.x&3)]-xpositions[(threadIdx.x>>2)+ilocal];\
                  ydist=vypositions[j+(threadIdx.x&3)]-ypositions[(threadIdx.x>>2)+ilocal];\
                  kernel(xvelocitytmp[threadIdx.x],yvelocitytmp[threadIdx.x],gammav[j+(threadIdx.x&3)],igammav[j+(threadIdx.x&3)]);\
                }\
                __syncthreads();\
                if(threadIdx.x<(blockDim.x>>2)) {/*sum*/\
                  xvelocity[threadIdx.x+ilocal]+=xvelocitytmp[threadIdx.x<<2]+xvelocitytmp[(threadIdx.x<<2)+1]+xvelocitytmp[(threadIdx.x<<2)+2]+xvelocitytmp[(threadIdx.x<<2)+3];\
                  CANSWER6/*yvelocity[threadIdx.x+ilocal]+=yvelocitytmp[threadIdx.x<<2]+yvelocitytmp[(threadIdx.x<<2)+1]+yvelocitytmp[(threadIdx.x<<2)+2]+yvelocitytmp[(threadIdx.x<<2)+3];*/\
                }\
                __syncthreads();\
                ilocal+=(blockDim.x>>2);\
              }\
              if((threadIdx.x>>2)+ilocal+(i-threadIdx.x)<count) {/*take care of the remaining points. More additional steps could be performed, but for high ndirect, adding more steps seems to decrease performance slightly*/\
                xvelocitytmp[threadIdx.x]=0;\
                CANSWER4/*yvelocitytmp[threadIdx.x]=0*/;\
                for(j=0;j<MAXTHREADSEVALUATION*SOURCEMULTIPLIER;j+=4) {\
                  xdist=vxpositions[j+(threadIdx.x&3)]-xpositions[(threadIdx.x>>2)+ilocal];\
                  ydist=vypositions[j+(threadIdx.x&3)]-ypositions[(threadIdx.x>>2)+ilocal];\
                  kernel(xvelocitytmp[threadIdx.x],yvelocitytmp[threadIdx.x],gammav[j+(threadIdx.x&3)],igammav[j+(threadIdx.x&3)]);\
                }\
              }\
              __syncthreads();\
              if(threadIdx.x+ilocal+(i-threadIdx.x)<count) {/*sum*/\
                xvelocity[threadIdx.x+ilocal]+=xvelocitytmp[threadIdx.x<<2]+xvelocitytmp[(threadIdx.x<<2)+1]+xvelocitytmp[(threadIdx.x<<2)+2]+xvelocitytmp[(threadIdx.x<<2)+3];\
                CANSWER6/*yvelocity[threadIdx.x+ilocal]+=yvelocitytmp[threadIdx.x<<2]+yvelocitytmp[(threadIdx.x<<2)+1]+yvelocitytmp[(threadIdx.x<<2)+2]+yvelocitytmp[(threadIdx.x<<2)+3];*/\
              }\
            }\
            __syncthreads();\
            k-=blockDim.x-(MAXTHREADSEVALUATION*SOURCEMULTIPLIER-addstart);/*All blockDim.x points was not read, correct the value for the next loop*/\
            addstart=0;\
            __syncthreads();\
            continue;\
          }\
          else if(end-k<=blockDim.x&&/*tex1Dfetch(tboxstart, boxnr+1)*/mend==m+1) {/*time to compute array due to no more available interaction points. Here, the list will not be full*/\
            if(k+threadIdx.x<end) {\
              TEXTUREFETCH(vxpositions[addstart+threadIdx.x],tzr,k+threadIdx.x);\
              TEXTUREFETCH(vypositions[addstart+threadIdx.x],tzi,k+threadIdx.x);\
              TEXTUREFETCH(gammav[addstart+threadIdx.x],tmr,k+threadIdx.x);\
              ISTRING2\
            }\
            addstart+=end-k;\
            __syncthreads();/*The remaining part is the same as above, but with the loops only running to addstart*/\
            if(count-(i-threadIdx.x)>=blockDim.x) {/*if more calculation points than threads*/\
              for(j=0;j<addstart;j++) {\
                xdist=vxpositions[j]-xpositions[threadIdx.x];\
                ydist=vypositions[j]-ypositions[threadIdx.x];\
                kernel(xvelocity[threadIdx.x],yvelocity[threadIdx.x],gammav[j],igammav[j]);\
              }\
            }\
            else { /*else, use several threads for each point to get better balance*/\
              ilocal=0;\
              if(count-(i-threadIdx.x)>=(blockDim.x>>1)){ /*more than half of the points active?*/\
                xvelocitytmp[threadIdx.x]=0;\
                CANSWER4/*yvelocitytmp[threadIdx.x]=0*/;\
                for(j=0;j+(threadIdx.x&1)<addstart;j+=2) {/*calculate temporary points*/\
                  xdist=vxpositions[j+(threadIdx.x&1)]-xpositions[(threadIdx.x>>1)];\
                  ydist=vypositions[j+(threadIdx.x&1)]-ypositions[(threadIdx.x>>1)];\
                  kernel(xvelocitytmp[threadIdx.x],yvelocitytmp[threadIdx.x],gammav[j+(threadIdx.x&1)],igammav[j+(threadIdx.x&1)]);\
                }\
                __syncthreads();\
                if(threadIdx.x<(blockDim.x>>1)) {/*sum*/\
                  xvelocity[threadIdx.x]+=xvelocitytmp[threadIdx.x<<1]+xvelocitytmp[(threadIdx.x<<1)+1];\
                  CANSWER5/*yvelocity[threadIdx.x]+=yvelocitytmp[threadIdx.x<<1]+yvelocitytmp[(threadIdx.x<<1)+1];*/\
                }\
                __syncthreads();\
                ilocal+=(blockDim.x>>1);\
              }\
              if(count-(i-threadIdx.x)-ilocal>=(blockDim.x>>2)){/*if remaining points are more than 1/4 of the number of threads*/\
                xvelocitytmp[threadIdx.x]=0;\
                CANSWER4/*yvelocitytmp[threadIdx.x]=0*/;\
                for(j=0;j+(threadIdx.x&3)<addstart;j+=4) {\
                  xdist=vxpositions[j+(threadIdx.x&3)]-xpositions[(threadIdx.x>>2)+ilocal];\
                  ydist=vypositions[j+(threadIdx.x&3)]-ypositions[(threadIdx.x>>2)+ilocal];\
                  kernel(xvelocitytmp[threadIdx.x],yvelocitytmp[threadIdx.x],gammav[j+(threadIdx.x&3)],igammav[j+(threadIdx.x&3)]);\
                }\
                __syncthreads();\
                if(threadIdx.x<(blockDim.x>>2)) {/*sum*/\
                  xvelocity[threadIdx.x+ilocal]+=xvelocitytmp[threadIdx.x<<2]+xvelocitytmp[(threadIdx.x<<2)+1]+xvelocitytmp[(threadIdx.x<<2)+2]+xvelocitytmp[(threadIdx.x<<2)+3];\
                  CANSWER6/*yvelocity[threadIdx.x+ilocal]+=yvelocitytmp[threadIdx.x<<2]+yvelocitytmp[(threadIdx.x<<2)+1]+yvelocitytmp[(threadIdx.x<<2)+2]+yvelocitytmp[(threadIdx.x<<2)+3];*/\
                }\
                __syncthreads();\
                ilocal+=(blockDim.x>>2);\
              }\
              if((threadIdx.x>>2)+ilocal+(i-threadIdx.x)<count) {/*take care of the remaining points. More additional steps could be performed, but for high ndirect, adding more steps seems to decrease performance slightly*/\
                xvelocitytmp[threadIdx.x]=0;\
                CANSWER4/*yvelocitytmp[threadIdx.x]=0*/;\
                for(j=0;j+(threadIdx.x&3)<addstart;j+=4) {\
                  xdist=vxpositions[j+(threadIdx.x&3)]-xpositions[(threadIdx.x>>2)+ilocal];\
                  ydist=vypositions[j+(threadIdx.x&3)]-ypositions[(threadIdx.x>>2)+ilocal];\
                  kernel(xvelocitytmp[threadIdx.x],yvelocitytmp[threadIdx.x],gammav[j+(threadIdx.x&3)],igammav[j+(threadIdx.x&3)]);\
                }\
              }\
              __syncthreads();\
              if(threadIdx.x+ilocal+(i-threadIdx.x)<count) {\
                xvelocity[threadIdx.x+ilocal]+=xvelocitytmp[threadIdx.x<<2]+xvelocitytmp[(threadIdx.x<<2)+1]+xvelocitytmp[(threadIdx.x<<2)+2]+xvelocitytmp[(threadIdx.x<<2)+3];\
                CANSWER6/*yvelocity[threadIdx.x+ilocal]+=yvelocitytmp[threadIdx.x<<2]+yvelocitytmp[(threadIdx.x<<2)+1]+yvelocitytmp[(threadIdx.x<<2)+2]+yvelocitytmp[(threadIdx.x<<2)+3];*/\
              }\
            }\
            addstart=0;\
            __syncthreads();\
            continue;\
          }\
          if(k+threadIdx.x<end) {/*If the list is not full, add the points and continue to next loop*/\
            TEXTUREFETCH(vxpositions[addstart+threadIdx.x],tzr,k+threadIdx.x);\
            TEXTUREFETCH(vypositions[addstart+threadIdx.x],tzi,k+threadIdx.x);\
            TEXTUREFETCH(gammav[addstart+threadIdx.x],tmr,k+threadIdx.x);\
            ISTRING2\
          }\
          if(k+blockDim.x<end)/*update addstart with the number of points added this loop*/\
            addstart+=blockDim.x;/*for sourcemultiplier==1, this is really unnecessary (but ifdefs cannot be used inside defines)*/\
          else\
            addstart+=end-k;\
        }\
      }\
      if(i<count) {\
        xvelocities[i]+=xvelocity[threadIdx.x];/*add the velocity*/\
        CANSWER3\
      }\
    }\
    __syncthreads();\
  }
    //for direct summation (tol==0) some simplifications can be made, since there are no mesh
extern __shared__ double vxpositions[];
#define SYNCNRONIZEDKERNELDIRECTSUM(kernel,tevalr,tevali,evalcount,ISTRING1,ISTRING2,ISTRING3,CANSWER1,CANSWER2,CANSWER3) int i,j,k;\
  double sqrdist,xdist,ydist,xvelocity;\
  CANSWER1\
  /*__shared__ double vxpositions[DIRECTMAXTHREADS];*/\
  __shared__ double *vypositions;/*__shared__ double vypositions[DIRECTMAXTHREADS]*/;\
  __shared__ double *gammav; /*__shared__ double gammav[DIRECTMAXTHREADS];*/\
  ISTRING1\
  double xposition;\
  double yposition;\
  __shared__ int loopmax,vcount;\
  if(threadIdx.x==0) {\
    vypositions=vxpositions+blockDim.x;\
    gammav=vypositions+blockDim.x;\
    ISTRING3\
    loopmax=evalcount+blockDim.x-(evalcount+blockDim.x-1)%blockDim.x-1;\
    vcount=count%blockDim.x;\
  }\
  __syncthreads();\
  for(i=threadIdx.x+blockIdx.x*blockDim.x;i<loopmax;i+=blockDim.x*gridDim.x) { /*if more points than threads*/\
    if(i<evalcount) {\
      xvelocity=0; /*use local variables as much as possible to help the optimizer understand that they are constant*/\
      /*yvelocity=0;*/CANSWER2\
      TEXTUREFETCH(xposition,tevalr,i);\
      TEXTUREFETCH(yposition,tevali,i);\
    }\
    for(k=0;k+blockDim.x<=count;k+=blockDim.x) {\
      TEXTUREFETCH(vxpositions[threadIdx.x],tzr,k+threadIdx.x);\
      TEXTUREFETCH(vypositions[threadIdx.x],tzi,k+threadIdx.x);\
      TEXTUREFETCH(gammav[threadIdx.x],tmr,k+threadIdx.x);\
      ISTRING2\
      __syncthreads();\
      if(i<evalcount) {\
        for(j=0;j<blockDim.x;j++) {\
          xdist=vxpositions[j]-xposition;\
          ydist=vypositions[j]-yposition;\
          kernel(xvelocity,yvelocity,circ,icirc);\
        }\
      }\
      __syncthreads();\
    }\
    if(k+threadIdx.x<count) {\
      TEXTUREFETCH(vxpositions[threadIdx.x],tzr,k+threadIdx.x);\
      TEXTUREFETCH(vypositions[threadIdx.x],tzi,k+threadIdx.x);\
      TEXTUREFETCH(gammav[threadIdx.x],tmr,k+threadIdx.x);\
      ISTRING2\
    }\
    __syncthreads();\
    if(i<evalcount) {\
      for(j=0;j<vcount;j++) {\
        xdist=vxpositions[j]-xposition;\
        ydist=vypositions[j]-yposition;\
        kernel(xvelocity,yvelocity,circ,icirc);\
      }\
      xvelocities[i]=xvelocity;\
      CANSWER3\
    }\
    __syncthreads();\
  }
#define circ gammav[j]
#define icirc igammav[j]

//choose between the old and new SYNCRONIZEDKERNEL by changing the defines
#ifdef SYNCKERNELOLD
#define IKERNELSTR2 IKERNELSTRING2
#define CANS1 COMPLEXANSWER1
#define CANS2 COMPLEXANSWER2
#define CANS3 COMPLEXANSWER3
#define CANS4
#define CANS5
#define CANS6
#define SYNCKERNEL SYNCNRONIZEDKERNEL

#else
#define IKERNELSTR2 IKERNELSTRING2_2
#define CANS1 COMPLEXANSWER1_2
#define CANS2 COMPLEXANSWER2_2
#define CANS3 COMPLEXANSWER3_2
#define CANS4 COMPLEXANSWER4_2
#define CANS5 COMPLEXANSWER5_2
#define CANS6 COMPLEXANSWER6_2
#define SYNCKERNEL SYNCNRONIZEDKERNEL2
#endif
//create all the functions
//there are one version that reads data from evaluation points, and one that reads from potential point.
__global__ void pot1_nomi_dirac_synchronized(double* xvelocities,double* yvelocities,int Nf,int* jcptr,int* j2cptr,int* ir DEBUGVECTORSTRING)
{
  SYNCKERNEL(DIRACSYNCKERNEL,ter,tei,tjxptr,DUMMYKERNELSTRING,DUMMYKERNELSTRING,CANS1,CANS2,CANS3,CANS4,CANS5,CANS6);
}
__global__ void pot1_mi_dirac_synchronized(double* xvelocities,double* yvelocities,int Nf,int* jcptr,int* j2cptr,int* ir DEBUGVECTORSTRING)
{
  SYNCKERNEL(IDIRACSYNCKERNEL,ter,tei,tjxptr,IKERNELSTRING1,IKERNELSTR2,CANS1,CANS2,CANS3,CANS4,CANS5,CANS6);
}
__global__ void pot1_nomi_rankine_synchronized(double* xvelocities,double* yvelocities,int Nf,int* jcptr,int* j2cptr,int* ir,double shape DEBUGVECTORSTRING)
{
  SYNCKERNEL(RANKINESYNCKERNEL,ter,tei,tjxptr,DUMMYKERNELSTRING,DUMMYKERNELSTRING,CANS1,CANS2,CANS3,CANS4,CANS5,CANS6);
}
__global__ void pot1_mi_rankine_synchronized(double* xvelocities,double* yvelocities,int Nf,int* jcptr,int* j2cptr,int* ir,double shape DEBUGVECTORSTRING)
{
  SYNCKERNEL(IRANKINESYNCKERNEL,ter,tei,tjxptr,IKERNELSTRING1,IKERNELSTR2,CANS1,CANS2,CANS3,CANS4,CANS5,CANS6);
}
__global__ void pot1_nomi_scully_synchronized(double* xvelocities,double* yvelocities,int Nf,int* jcptr,int* j2cptr,int* ir,double shape,double scale,double cutoffsqr DEBUGVECTORSTRING)
{
  SYNCKERNEL(SCULLYSYNCKERNEL,ter,tei,tjxptr,DUMMYKERNELSTRING,DUMMYKERNELSTRING,CANS1,CANS2,CANS3,CANS4,CANS5,CANS6);
}
__global__ void pot1_mi_scully_synchronized(double* xvelocities,double* yvelocities,int Nf,int* jcptr,int* j2cptr,int* ir,double shape,double scale,double cutoffsqr DEBUGVECTORSTRING)
{
  SYNCKERNEL(ISCULLYSYNCKERNEL,ter,tei,tjxptr,IKERNELSTRING1,IKERNELSTR2,CANS1,CANS2,CANS3,CANS4,CANS5,CANS6);
}
__global__ void pot1_nomi_oseen_synchronized(double* xvelocities,double* yvelocities,int Nf,int* jcptr,int* j2cptr,int* ir,double shape,double scale,double cutoffsqr DEBUGVECTORSTRING)
{
  SYNCKERNEL(OSEENSYNCKERNEL,ter,tei,tjxptr,DUMMYKERNELSTRING,DUMMYKERNELSTRING,CANS1,CANS2,CANS3,CANS4,CANS5,CANS6);
}
__global__ void pot1_mi_oseen_synchronized(double* xvelocities,double* yvelocities,int Nf,int* jcptr,int* j2cptr,int* ir,double shape,double scale,double cutoffsqr DEBUGVECTORSTRING)
{
  SYNCKERNEL(IOSEENSYNCKERNEL,ter,tei,tjxptr,IKERNELSTRING1,IKERNELSTR2,CANS1,CANS2,CANS3,CANS4,CANS5,CANS6);
}

//log kernels
__global__ void pot0_nomi_dirac_synchronized(double* xvelocities,int Nf,int* jcptr,int* j2cptr,int* ir DEBUGVECTORSTRING)
{
  SYNCKERNEL(DIRACLOGSYNCKERNEL,ter,tei,tjxptr,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING);
}
__global__ void pot0_mi_dirac_synchronized(double* xvelocities,double* yvelocities,int Nf,int* jcptr,int* j2cptr,int* ir DEBUGVECTORSTRING)
{
  SYNCKERNEL(IDIRACLOGSYNCKERNEL,ter,tei,tjxptr,IKERNELSTRING1,IKERNELSTR2,CANS1,CANS2,CANS3,CANS4,CANS5,CANS6);
}
__global__ void pot0_nomi_rankine_synchronized(double* xvelocities,int Nf,int* jcptr,int* j2cptr,int* ir,double shape,double cutoffsqr DEBUGVECTORSTRING)
{
  SYNCKERNEL(RANKINELOGSYNCKERNEL,ter,tei,tjxptr,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING);
}
__global__ void pot0_mi_rankine_synchronized(double* xvelocities,double* yvelocities,int Nf,int* jcptr,int* j2cptr,int* ir,double shape,double cutoffsqr DEBUGVECTORSTRING)
{
  SYNCKERNEL(IRANKINELOGSYNCKERNEL,ter,tei,tjxptr,IKERNELSTRING1,IKERNELSTR2,CANS1,CANS2,CANS3,CANS4,CANS5,CANS6);
}
__global__ void pot0_nomi_scully_synchronized(double* xvelocities,int Nf,int* jcptr,int* j2cptr,int* ir,double shape,double scale,double cutoffsqr DEBUGVECTORSTRING)
{
  SYNCKERNEL(SCULLYLOGSYNCKERNEL,ter,tei,tjxptr,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING);
}
__global__ void pot0_mi_scully_synchronized(double* xvelocities,double* yvelocities,int Nf,int* jcptr,int* j2cptr,int* ir,double shape,double scale,double cutoffsqr DEBUGVECTORSTRING)
{
  SYNCKERNEL(ISCULLYLOGSYNCKERNEL,ter,tei,tjxptr,IKERNELSTRING1,IKERNELSTR2,CANS1,CANS2,CANS3,CANS4,CANS5,CANS6);
}
__global__ void pot0_nomi_oseen_synchronized(double* xvelocities,int Nf,int* jcptr,int* j2cptr,int* ir,double shape,double scale,double cutoffsqr DEBUGVECTORSTRING)
{
  SYNCKERNEL(OSEENLOGSYNCKERNEL,ter,tei,tjxptr,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING);
}
__global__ void pot0_mi_oseen_synchronized(double* xvelocities,double* yvelocities,int Nf,int* jcptr,int* j2cptr,int* ir,double shape,double scale,double cutoffsqr DEBUGVECTORSTRING)
{
  SYNCKERNEL(IOSEENLOGSYNCKERNEL,ter,tei,tjxptr,IKERNELSTRING1,IKERNELSTR2,CANS1,CANS2,CANS3,CANS4,CANS5,CANS6);
}

//direct eval kernels
__global__ void direct_pot1_nomi_dirac_synchronized(double* xvelocities,double* yvelocities,int count,int evalcount)
{
  SYNCNRONIZEDKERNELDIRECTSUM(DIRACSYNCKERNEL,ter,tei,evalcount,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,COMPLEXANSWER1,COMPLEXANSWER2,COMPLEXANSWER3DIRECT);
}
__global__ void direct_pot1_mi_dirac_synchronized(double* xvelocities,double* yvelocities,int count,int evalcount)
{
  SYNCNRONIZEDKERNELDIRECTSUM(IDIRACSYNCKERNEL,ter,tei,evalcount,IKERNELSTRING1DIRECT,IKERNELSTRING2,IKERNELSTRING1DIRECT2,COMPLEXANSWER1,COMPLEXANSWER2,COMPLEXANSWER3DIRECT);
}
__global__ void direct_pot1_nomi_rankine_synchronized(double* xvelocities,double* yvelocities,int count,int evalcount,double shape)
{
  SYNCNRONIZEDKERNELDIRECTSUM(RANKINESYNCKERNEL,ter,tei,evalcount,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,COMPLEXANSWER1,COMPLEXANSWER2,COMPLEXANSWER3DIRECT);
}
__global__ void direct_pot1_mi_rankine_synchronized(double* xvelocities,double* yvelocities,int count,int evalcount,double shape)
{
  SYNCNRONIZEDKERNELDIRECTSUM(IRANKINESYNCKERNEL,ter,tei,evalcount,IKERNELSTRING1DIRECT,IKERNELSTRING2,IKERNELSTRING1DIRECT2,COMPLEXANSWER1,COMPLEXANSWER2,COMPLEXANSWER3DIRECT);
}
__global__ void direct_pot1_nomi_scully_synchronized(double* xvelocities,double* yvelocities,int count,int evalcount,double shape,double scale,double cutoffsqr)
{
  SYNCNRONIZEDKERNELDIRECTSUM(SCULLYSYNCKERNEL,ter,tei,evalcount,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,COMPLEXANSWER1,COMPLEXANSWER2,COMPLEXANSWER3DIRECT);
}
__global__ void direct_pot1_mi_scully_synchronized(double* xvelocities,double* yvelocities,int count,int evalcount,double shape,double scale,double cutoffsqr)
{
  SYNCNRONIZEDKERNELDIRECTSUM(ISCULLYSYNCKERNEL,ter,tei,evalcount,IKERNELSTRING1DIRECT,IKERNELSTRING2,IKERNELSTRING1DIRECT2,COMPLEXANSWER1,COMPLEXANSWER2,COMPLEXANSWER3DIRECT);
}
__global__ void direct_pot1_nomi_oseen_synchronized(double* xvelocities,double* yvelocities,int count,int evalcount,double shape,double scale,double cutoffsqr)
{
  SYNCNRONIZEDKERNELDIRECTSUM(OSEENSYNCKERNEL,ter,tei,evalcount,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,COMPLEXANSWER1,COMPLEXANSWER2,COMPLEXANSWER3DIRECT);
}
__global__ void direct_pot1_mi_oseen_synchronized(double* xvelocities,double* yvelocities,int count,int evalcount,double shape,double scale,double cutoffsqr)
{
  SYNCNRONIZEDKERNELDIRECTSUM(IOSEENSYNCKERNEL,ter,tei,evalcount,IKERNELSTRING1DIRECT,IKERNELSTRING2,IKERNELSTRING1DIRECT2,COMPLEXANSWER1,COMPLEXANSWER2,COMPLEXANSWER3DIRECT);
}
//log kernels
__global__ void direct_pot0_nomi_dirac_synchronized(double* xvelocities,int count,int evalcount)
{
  SYNCNRONIZEDKERNELDIRECTSUM(DIRACLOGSYNCKERNEL,ter,tei,evalcount,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING);
}
__global__ void direct_pot0_mi_dirac_synchronized(double* xvelocities,double* yvelocities,int count,int evalcount)
{
  SYNCNRONIZEDKERNELDIRECTSUM(IDIRACLOGSYNCKERNEL,ter,tei,evalcount,IKERNELSTRING1DIRECT,IKERNELSTRING2,IKERNELSTRING1DIRECT2,COMPLEXANSWER1,COMPLEXANSWER2,COMPLEXANSWER3DIRECT);
}
__global__ void direct_pot0_nomi_rankine_synchronized(double* xvelocities,int count,int evalcount,double shape,double cutoffsqr)
{
  SYNCNRONIZEDKERNELDIRECTSUM(RANKINELOGSYNCKERNEL,ter,tei,evalcount,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING);
}
__global__ void direct_pot0_mi_rankine_synchronized(double* xvelocities,double* yvelocities,int count,int evalcount,double shape,double cutoffsqr)
{
  SYNCNRONIZEDKERNELDIRECTSUM(IRANKINELOGSYNCKERNEL,ter,tei,evalcount,IKERNELSTRING1DIRECT,IKERNELSTRING2,IKERNELSTRING1DIRECT2,COMPLEXANSWER1,COMPLEXANSWER2,COMPLEXANSWER3DIRECT);
}
__global__ void direct_pot0_nomi_scully_synchronized(double* xvelocities,int count,int evalcount,double shape,double scale,double cutoffsqr)
{
  SYNCNRONIZEDKERNELDIRECTSUM(SCULLYLOGSYNCKERNEL,ter,tei,evalcount,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING);
}
__global__ void direct_pot0_mi_scully_synchronized(double* xvelocities,double* yvelocities,int count,int evalcount,double shape,double scale,double cutoffsqr)
{
  SYNCNRONIZEDKERNELDIRECTSUM(ISCULLYLOGSYNCKERNEL,ter,tei,evalcount,IKERNELSTRING1DIRECT,IKERNELSTRING2,IKERNELSTRING1DIRECT2,COMPLEXANSWER1,COMPLEXANSWER2,COMPLEXANSWER3DIRECT);
}
__global__ void direct_pot0_nomi_oseen_synchronized(double* xvelocities,int count,int evalcount,double shape,double scale,double cutoffsqr)
{
  SYNCNRONIZEDKERNELDIRECTSUM(OSEENLOGSYNCKERNEL,ter,tei,evalcount,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING,DUMMYKERNELSTRING);
}
__global__ void direct_pot0_mi_oseen_synchronized(double* xvelocities,double* yvelocities,int count,int evalcount,double shape,double scale,double cutoffsqr)
{
  SYNCNRONIZEDKERNELDIRECTSUM(IOSEENLOGSYNCKERNEL,ter,tei,evalcount,IKERNELSTRING1DIRECT,IKERNELSTRING2,IKERNELSTRING1DIRECT2,COMPLEXANSWER1,COMPLEXANSWER2,COMPLEXANSWER3DIRECT);
}
#undef circ


#define tevalr ter
#define tevali tei
#define DEBUGCHECKBOX 250

//multipole evaluations on gpu
#ifdef MULTIPOLEEVALCUDA
__global__ void cuda_mpexp_eval(int p, double *xvelocities, double *yvelocities, int Nf,int fetchbase,int *j2cptr,int *kcptr,int *ir,int pot DEBUGVECTORSTRING)
        /* Evaluation of near field at coordinates (zr,zi). */
{
  __shared__ dcmplx coeff[128]; //max 128 coefficients supported here
  __shared__ int count, loopmax, m1, base, distbase, k; /*save space, same for all threads*/
  __shared__ double xpositions[MAXTHREADSEVALUATION];/*makes code approx 1% slower for heavy load, but for some reason, doubles has to be fetched into an array*/
  __shared__ double ypositions[MAXTHREADSEVALUATION];
  __shared__ SORT_DCMPLX z0, dz0;
  double xvelocity, yvelocity;
  double re, im, pre, pim, pre0;
  int i, j, boxnr, m;
  for(boxnr=blockIdx.x;boxnr<Nf;boxnr+=gridDim.x) { //if more boxes than number of blocks, loop the code
    if(threadIdx.x==0) {
      m1=kcptr[boxnr];//distant box interactions
      count=tex1Dfetch(tjxptr, boxnr+1); /*right now, it is possible that two boxes has same evalbegin if one box has 0 elements*/
      loopmax=count-(count+blockDim.x-1-(tex1Dfetch(tjxptr, boxnr)%blockDim.x))%blockDim.x-1+blockDim.x; //loop end, this is to make sure all threads run the same number of loops, since __syncthreads() has to be called by all threads at the same time
      base=(boxnr+fetchbase)*(p+1); //to read coefficients and z0 positions, find the point in the array
    }
    
    if(threadIdx.x<sizeof(SORT_DCMPLX)/sizeof(float)) { //the box position
      ((float*)&z0)[threadIdx.x] =tex1Dfetch(tz0, (boxnr+fetchbase)*sizeof(SORT_DCMPLX)/sizeof(float)+threadIdx.x);
    }
    __syncthreads();
    //start by performing the interaction within the box
    if(p<127) { //should perhaps be corrected to handle more than this many coefficients
      for(j=threadIdx.x;j<p+1;j+=blockDim.x) { //read coefficients to shared memory
        fetch_double2((double2*)&coeff[j],tcoeff2d,base+j);
      }
      //__syncthreads(); /*should not be necessary*/
      for(i=threadIdx.x+tex1Dfetch(tjxptr, boxnr);i<loopmax;i+=blockDim.x) { /*if more points than threads*/
        if(i<count) { //get local positions
          TEXTUREFETCH(xpositions[threadIdx.x], tevalr, i);
          TEXTUREFETCH(ypositions[threadIdx.x], tevali, i);
        }
        __syncthreads();
        if(i<count) { //if the thread has a point to evaluate
          //use same algorithm as CPU code to perform multipole evaluation in box
          re=xpositions[threadIdx.x]-creal(z0);
          im=ypositions[threadIdx.x]-cimag(z0);
          pre=creal(coeff[p]);
          pim=cimag(coeff[p]);
          
          for (j = p-1; j >= 0; j--) {
            pre0 = re*pre-im*pim+creal(coeff[j]);
            pim = re*pim+im*pre+cimag(coeff[j]);
            pre = pre0;
          }
          xvelocity=pre;
          yvelocity=pim;
        }
        
        //now, continue and perform all interactions with other boxes according to the lists
        for(m=j2cptr[boxnr];m<m1;m++) {
          __syncthreads(); //if the continue statement is true
          if(threadIdx.x==0){ //the base value for the corresponding box
            k = -ir[m];
            if(k>0)
              distbase=(k+fetchbase)*(p+1);
            
          }
          __syncthreads();
          if(k<0) {
            if(m+1==m1&&i+blockDim.x<loopmax) { //if final step, reload coefficient matrix with previous values if necessary
              for(j=threadIdx.x;j<p+1;j+=blockDim.x) { //read coefficients to shared memory
                fetch_double2((double2*)&coeff[j],tcoeff2d,base+j);
              }
            }
            continue;
          }
          for(j=threadIdx.x;j<p+1;j+=blockDim.x) { //read coefficients to shared memory
            fetch_double2((double2*)&coeff[j],tcoeff1d,distbase+j);
          }
          if(threadIdx.x<sizeof(SORT_DCMPLX)/sizeof(float)) { //position of corresponding box
            ((float*)&dz0)[threadIdx.x] =tex1Dfetch(tz0, (k+fetchbase)*sizeof(SORT_DCMPLX)/sizeof(float)+threadIdx.x);
          }
          __syncthreads();
          if(i<count) {//if the thread has a point to evaluate
            //use same algorithm as CPU code to perform multipole evaluation between boxes
            pre=xpositions[threadIdx.x]-creal(dz0);
            pim=ypositions[threadIdx.x]-cimag(dz0);
            pre0=pre*pre+pim*pim;
            re=pre/pre0;
            im=-pim/pre0;
            pre=creal(coeff[p]);
            pim=cimag(coeff[p]);
            for (j = p-1; j > 0; j--) {
              pre0 = re*pre-im*pim+creal(coeff[j]);
              pim = re*pim+im*pre+cimag(coeff[j]);
              pre = pre0;
            } {
              pre0 = re*pre-im*pim;
              pim = re*pim+im*pre;
              pre = pre0;
              if (pot == 0) {
                pre0=0.5*log(re*re+im*im);
                im=atan2(im, re);
                pre -= creal(coeff[0])*pre0-cimag(coeff[0])*im;
                pim -= creal(coeff[0])*im+cimag(coeff[0])*pre0;
              }
            }
            xvelocity+=pre;
            yvelocity+=pim;
            
          }
          __syncthreads();
          if(m+1==m1&&i+blockDim.x<loopmax) { //if final step, reload coefficient matrix with previous values if necessary
            for(j=threadIdx.x;j<p+1;j+=blockDim.x) { //read coefficients to shared memory
              fetch_double2((double2*)&coeff[j],tcoeff2d,base+j);
            }
          }
        }
        __syncthreads();
        if(i<count) { //write velocity to global memory
          xvelocities[i]+=xvelocity;
          if(yvelocities!=NULL)
            yvelocities[i]+=yvelocity;
        }
      }
      __syncthreads();
    }
  }
}
#endif
#undef tevalr
#undef tevali

#ifdef MULTIPOLEINITCUDA
#include "cudainit.h"

extern __shared__ double wksp[];
extern __shared__ dcmplx localcoeffs[];
//performs m2m interaction using horner based scheme
__global__ void shift_m2m_cuda_horner(int parentlevel,int parentstart,int pshift,
               dcmplx *This,const SORT_DCMPLX *z,int potshift  DEBUGVECTORSTRING)
{
//     dcmplx localcoeff[PMAX];
  dcmplx *localcoeff;
  dcmplx tmp;
  SORT_DCMPLX z0;
  __shared__ SORT_DCMPLX z0base[(SHIFT_M2M_MAXTHREADS+3)/4];
  int i, j;
  for(int boxnr=blockIdx.x*blockDim.x/4;boxnr<(1<<(parentlevel<<1));boxnr+=gridDim.x*blockDim.x/4) { //loop over all boxes. One block takes care of blockDim.x/4 boxes
    localcoeff=(dcmplx*)&localcoeffs[threadIdx.x*(pshift+1)];
    if(threadIdx.x<blockDim.x/4&&threadIdx.x<(1<<(parentlevel<<1))) { /*only works if sizeof(dcmplx)/sizeof(float)<=4, which is the case for both double and float*/
      z0base[threadIdx.x]=z[boxnr+parentstart+threadIdx.x];
    }
    __syncthreads();
    if(boxnr+(threadIdx.x>>2)<(1<<(parentlevel<<1))) { /*safety check. For standard parameters, this should only be an issue for parentlevel = 1 - 3, then, it should be a multiple of blockDim/4*/
      z0=z[(boxnr<<2)+parentstart+(1<<(parentlevel<<1))+threadIdx.x];
      //calculate distances
      COMPLEXASSIGN(z0, creal(z0) - creal(z0base[threadIdx.x>>2]), cimag(z0) - cimag(z0base[threadIdx.x>>2])); //be careful with this use of COMPLEXASSIGN, remember that it is a macro
      j=((boxnr<<2)+parentstart+(1<<(parentlevel<<1))+threadIdx.x)*(pshift+1);
      for(i=0;i<pshift+1;i++) //read coefficient parameters
        COMPLEXASSIGN(localcoeff[i], creal(This[j+i]), cimag(This[j+i]));
      for(j=pshift;j>=2;j--) { //perform the Horner scheme
        for(i=j;i<=pshift;i++) {
          COMPLEXADD(localcoeff[i], creal(z0)*creal(localcoeff[i-1])-cimag(z0)*cimag(localcoeff[i-1]), creal(z0)*cimag(localcoeff[i-1])+cimag(z0)*creal(localcoeff[i-1]));
        }
      }
    }
    if(potshift==0) {
//       tmp=z0;
      COMPLEXASSIGN(tmp,creal(z0),cimag(z0));
      COMPLEXSUB(localcoeff[1], (creal(localcoeff[0])*creal(tmp)-cimag(localcoeff[0])*cimag(tmp)), (creal(localcoeff[0])*cimag(tmp)+cimag(localcoeff[0])*creal(tmp)));
      for(i=2;i<=pshift;i++) {
        double dtmp=creal(z0)*cimag(tmp)+cimag(z0)*creal(tmp); //depending on the definition of COMPLEXASSIGN, a poorly designed macro will overwrite tmp at first step, and use new value in second stage (normally not the case, but this is to prevent future problems)
        COMPLEXASSIGN(tmp, creal(tmp)*creal(z0)-cimag(tmp)*cimag(z0), dtmp);
        COMPLEXSUB(localcoeff[i], (creal(localcoeff[0])*creal(tmp)-cimag(localcoeff[0])*cimag(tmp))*(1.0/i), (creal(localcoeff[0])*cimag(tmp)+cimag(localcoeff[0])*creal(tmp))*(1.0/i));
      }
    }
    __syncthreads();
    //calculate the sum of the for contributing boxes
    if(boxnr+(threadIdx.x>>2)<(1<<(parentlevel<<1))) { /*safety check. For standard parameters, this should only be an issue for parentlevel=1 and 2, then, it should be a multiple of blockDim/4*/
      //write results, here, each thread has its own box to write results to
      localcoeff=(dcmplx*)&localcoeffs[(threadIdx.x&(~1))*(pshift+1)];
      for(i=(threadIdx.x&1);i<=pshift;i+=2) {
        COMPLEXADD(localcoeff[i], creal(localcoeff[i+pshift+1]), cimag(localcoeff[i+pshift+1]));
      }
    }
    __syncthreads();
    if(boxnr+(threadIdx.x>>2)<(1<<(parentlevel<<1))) { /*safety check. For standard parameters, this should only be an issue for parentlevel=1 and 2, then, it should be a multiple of blockDim/4*/
      //write results, here, each thread has its own box to write results to
      localcoeff=(dcmplx*)&localcoeffs[(threadIdx.x&(~3))*(pshift+1)];
      j=(boxnr+(threadIdx.x>>2)+parentstart)*(pshift+1);
      for(i=(threadIdx.x&3);i<=pshift;i+=4) {
        COMPLEXADD(localcoeff[i], creal(localcoeff[i+2*(pshift+1)]), cimag(localcoeff[i+2*(pshift+1)]));
        COMPLEXASSIGN(This[j+i], creal(localcoeff[i]), cimag(localcoeff[i]));
      }
    }
  }
}


//performs m2m interaction using horner based scheme, with pre and post scaling
__global__ void shift_m2m_cuda_horner_scaled(int parentlevel,int parentstart,int pshift,
               dcmplx *This,const SORT_DCMPLX *z,int potshift,double rminshift  DEBUGVECTORSTRING)
{
//     dcmplx localcoeff[PMAX];
  dcmplx *localcoeff;
  dcmplx tmp, invr;
  double rtmp, itmp;
  SORT_DCMPLX z0;
  __shared__ SORT_DCMPLX z0base[(SHIFT_M2M_MAXTHREADS+3)/4];
  int i, j, k;
  for(int boxnr=blockIdx.x*blockDim.x/8;boxnr<(1<<(parentlevel<<1));boxnr+=gridDim.x*blockDim.x/8) { //loop over all boxes. One block takes care of blockDim.x/4 boxes
    localcoeff=(dcmplx*)&localcoeffs[(threadIdx.x>>1)*(pshift+1)]; //two threads for each set of coefficients
    if(threadIdx.x<blockDim.x/8&&threadIdx.x<(1<<(parentlevel<<1))) { /*only works if sizeof(dcmplx)/sizeof(float)<=4, which is the case for both double and float*/
      z0base[threadIdx.x]=z[boxnr+parentstart+threadIdx.x];
    }
    __syncthreads();
    if(boxnr+(threadIdx.x>>3)<(1<<(parentlevel<<1))) { /*safety check. For standard parameters, this should only be an issue for parentlevel = 1 - 3, then, it should be a multiple of blockDim/4*/
      z0=z[(boxnr<<2)+parentstart+(1<<(parentlevel<<1))+threadIdx.x/2];
      //calculate distances
      COMPLEXASSIGN(z0, creal(z0) - creal(z0base[threadIdx.x>>3]), cimag(z0) - cimag(z0base[threadIdx.x>>3])); //be careful with this use of COMPLEXASSIGN, remember that it is a macro
      k=((boxnr<<2)+parentstart+(1<<(parentlevel<<1))+threadIdx.x/2)*(pshift+1);
      for(i=(threadIdx.x&1);i<pshift+1;i+=2) //read coefficient parameters
        COMPLEXASSIGN(localcoeff[i], creal(This[k+i]), cimag(This[k+i]));
      //safe calculation of 1/z0
      if(fabs(creal(z0))>fabs(cimag(z0))) {
        itmp=cimag(z0)*(1.0/creal(z0));
        rtmp=creal(z0)+itmp*cimag(z0);
        COMPLEXASSIGN(tmp, 1.0/rtmp, -itmp*(1.0/rtmp));
      }
      else {
        itmp=creal(z0)*(1.0/cimag(z0));
        rtmp=cimag(z0)+itmp*creal(z0);
        COMPLEXASSIGN(tmp, itmp*(1.0/rtmp), -1.0/rtmp);
      }
      COMPLEXASSIGN(invr, creal(tmp)*creal(tmp)-cimag(tmp)*cimag(tmp), creal(tmp)*cimag(tmp)*2);
      if(fabs(creal(z0))<rminshift&&fabs(cimag(z0))<rminshift) { //if overflow is possible, run safe code, i.e. traditional horner shift
        if(threadIdx.x&1) {
          for(k=pshift;k>=2;k--) { //perform the Horner scheme
            for(i=k;i<=pshift;i++) {
              COMPLEXADD(localcoeff[i], creal(z0)*creal(localcoeff[i-1])-cimag(z0)*cimag(localcoeff[i-1]), creal(z0)*cimag(localcoeff[i-1])+cimag(z0)*creal(localcoeff[i-1]));
            }
          }
          if(potshift==0) {
//             tmp=z0;
            COMPLEXASSIGN(tmp,creal(z0),cimag(z0));
            COMPLEXSUB(localcoeff[1], (creal(localcoeff[0])*creal(tmp)-cimag(localcoeff[0])*cimag(tmp)), (creal(localcoeff[0])*cimag(tmp)+cimag(localcoeff[0])*creal(tmp)));
            for(i=2;i<=pshift;i++) {
              double dtmp=creal(z0)*cimag(tmp)+cimag(z0)*creal(tmp); //depending on the definition of COMPLEXASSIGN, a poorly designed macro will overwrite tmp at first step, and use new value in second stage (normally not the case, but this is to prevent future problems)
              COMPLEXASSIGN(tmp, creal(tmp)*creal(z0)-cimag(tmp)*cimag(z0), dtmp);
              COMPLEXSUB(localcoeff[i], (creal(localcoeff[0])*creal(tmp)-cimag(localcoeff[0])*cimag(tmp))*(1.0/i), (creal(localcoeff[0])*cimag(tmp)+cimag(localcoeff[0])*creal(tmp))*(1.0/i));
            }
          }
        }
      }
      else {
        //perform scaling
        if((threadIdx.x&1)==1) { //set up starting values for the two threads working on the same set of data
          itmp=creal(localcoeff[1])*cimag(tmp)+cimag(localcoeff[1])*creal(tmp);
          COMPLEXASSIGN(localcoeff[1], creal(localcoeff[1])*creal(tmp)-cimag(localcoeff[1])*cimag(tmp), itmp);
        }
        else {
          COMPLEXASSIGN(tmp, 1, 0);
        }
        for(k=2+(threadIdx.x&1);k<=pshift;k+=2) {
          itmp=creal(invr)*cimag(tmp)+cimag(invr)*creal(tmp);
          COMPLEXASSIGN(tmp, creal(invr)*creal(tmp)-cimag(invr)*cimag(tmp), itmp);
          itmp=creal(localcoeff[k])*cimag(tmp)+cimag(localcoeff[k])*creal(tmp);
          COMPLEXASSIGN(localcoeff[k], creal(localcoeff[k])*creal(tmp)-cimag(localcoeff[k])*cimag(tmp), itmp);
        }
//         for(k=pshift;k>=2;k--) { //perform the Horner scheme
//           for(i=k;i<=pshift;i++) {
// //                    COMPLEXADD(localcoeff[i], creal(localcoeff[i-1]), cimag(localcoeff[i-1]));
//             ((double*)&localcoeff[i])[threadIdx.x&1]+=((double*)&localcoeff[i-1])[threadIdx.x&1];
//           }
//         }
        k=pshift;
        double a[UNROLLLEVEL+1];
        #if UNROLLLEVEL>=7
        for (; k>7; k-=7) {
          INNERUNROLLBASE8(wksp,0,STOREMDNM2M,LOADM1M2M,ADDEDN,1,k)
        }
        if (k>6)
        #elif UNROLLLEVEL>=6
        for (;k>6;)
        #endif
        #if UNROLLLEVEL>=6
        {
          INNERUNROLLBASE7(wksp,0,STOREMDNM2M,LOADM1M2M,ADDEDN,1,k)
          k-=6;
        }
        if (k>5)
        #elif UNROLLLEVEL>=5
        for (;k>5;)
        #endif
        #if UNROLLLEVEL>=5
        {
          INNERUNROLLBASE6(wksp,0,STOREMDNM2M,LOADM1M2M,ADDEDN,1,k)
          k-=5;
        }
        if (k>4)
        #elif UNROLLLEVEL>=4
        for (;k>4;)
        #endif
        #if UNROLLLEVEL>=4
        {
          INNERUNROLLBASE5(wksp,0,STOREMDNM2M,LOADM1M2M,ADDEDN,1,k)
          k-=4;
        }
        if (k>3)
        #elif UNROLLLEVEL>=3
        for (;k>3;)
        #endif
        #if UNROLLLEVEL>=3
        {
          INNERUNROLLBASE4(wksp,0,STOREMDNM2M,LOADM1M2M,ADDEDN,1,k)
          k-=3;
        }
        if (k>2)
        #elif UNROLLLEVEL>=2
        for (;k>2;)
        #endif
        #if UNROLLLEVEL>=2
        {
          INNERUNROLLBASE3(wksp,0,STOREMDNM2M,LOADM1M2M,ADDEDN,1,k)
          k-=2;
        }
        if (k > 1)
        #else
        for (;k > 1;k--)
        #endif
        {
          INNERUNROLLBASE2(wksp,0,STOREMDNM2M,LOADM1M2M,ADDEDN,1,k)
        }
        
        if(potshift==0) { //depending on logarithmic potential or not, an additional term may be necesarry
          //perform scaling
//           tmp=z0;
          COMPLEXASSIGN(tmp,creal(z0),cimag(z0));
          COMPLEXASSIGN(invr, creal(tmp)*creal(tmp)-cimag(tmp)*cimag(tmp), creal(tmp)*cimag(tmp)*2);
          if((threadIdx.x&1)==1) {
            COMPLEXSUB(localcoeff[1], creal(localcoeff[0]), cimag(localcoeff[0]));
            itmp=(creal(localcoeff[1])*cimag(tmp)+cimag(localcoeff[1])*creal(tmp));
            COMPLEXASSIGN(localcoeff[1], (creal(localcoeff[1])*creal(tmp)-cimag(localcoeff[1])*cimag(tmp)), itmp);
          }
          else {
            COMPLEXASSIGN(tmp, 1, 0);
          }
          for(k=2+(threadIdx.x&1);k<=pshift;k+=2) {
            itmp=creal(invr)*cimag(tmp)+cimag(invr)*creal(tmp); //depending on the definition of COMPLEXASSIGN, a poorly designed macro will overwrite tmp at first step, and use new value in second stage (normally not the case, but this is to prevent future problems)
            COMPLEXASSIGN(tmp, creal(tmp)*creal(invr)-cimag(tmp)*cimag(invr), itmp);
            COMPLEXSUB(localcoeff[k], creal(localcoeff[0])*(1.0/k), cimag(localcoeff[0])*(1.0/k));
            itmp=creal(localcoeff[k])*cimag(tmp)+cimag(localcoeff[k])*creal(tmp);
            COMPLEXASSIGN(localcoeff[k], creal(localcoeff[k])*creal(tmp)-cimag(localcoeff[k])*cimag(tmp), itmp);
          }
        }
        else {
          //perform scaling
//           tmp=z0;
          COMPLEXASSIGN(tmp,creal(z0),cimag(z0));
          COMPLEXASSIGN(invr, creal(tmp)*creal(tmp)-cimag(tmp)*cimag(tmp), creal(tmp)*cimag(tmp)*2);
          if((threadIdx.x&1)==1) {
            itmp=(creal(localcoeff[1])*cimag(tmp)+cimag(localcoeff[1])*creal(tmp));
            COMPLEXASSIGN(localcoeff[1], (creal(localcoeff[1])*creal(tmp)-cimag(localcoeff[1])*cimag(tmp)), itmp);
          }
          else {
            COMPLEXASSIGN(tmp, 1, 0);
          }
          for(k=2+(threadIdx.x&1);k<=pshift;k+=2) {
            itmp=creal(invr)*cimag(tmp)+cimag(invr)*creal(tmp); //depending on the definition of COMPLEXASSIGN, a poorly designed macro will overwrite tmp at first step, and use new value in second stage (normally not the case, but this is to prevent future problems)
            COMPLEXASSIGN(tmp, creal(tmp)*creal(invr)-cimag(tmp)*cimag(invr), itmp);
            itmp=creal(localcoeff[k])*cimag(tmp)+cimag(localcoeff[k])*creal(tmp);
            COMPLEXASSIGN(localcoeff[k], creal(localcoeff[k])*creal(tmp)-cimag(localcoeff[k])*cimag(tmp), itmp);
          }
        }
      }
    }
            
    __syncthreads();
    //calculate the sum of the for contributing boxes
    if(boxnr+(threadIdx.x>>3)<(1<<(parentlevel<<1))) { /*safety check. For standard parameters, this should only be an issue for parentlevel=1 and 2, then, it should be a multiple of blockDim/4*/
      //write results, here, each thread has its own box to write results to
      localcoeff=(dcmplx*)&localcoeffs[((threadIdx.x>>1)&(~1))*(pshift+1)];
      for(i=(threadIdx.x&3);i<=pshift;i+=4) {
        COMPLEXADD(localcoeff[i], creal(localcoeff[i+pshift+1]), cimag(localcoeff[i+pshift+1]));
      }
    }
    __syncthreads();
    //second part of the sum, and write results back to array
    if(boxnr+(threadIdx.x>>3)<(1<<(parentlevel<<1))) { /*safety check. For standard parameters, this should only be an issue for parentlevel=1 and 2, then, it should be a multiple of blockDim/4*/
      //write results, here, each thread has its own box to write results to
      localcoeff=(dcmplx*)&localcoeffs[((threadIdx.x>>1)&(~3))*(pshift+1)];
      k=(boxnr+(threadIdx.x>>3)+parentstart)*(pshift+1);
      for(i=(threadIdx.x&7);i<=pshift;i+=8) {
        COMPLEXADD(localcoeff[i], creal(localcoeff[i+2*(pshift+1)]), cimag(localcoeff[i+2*(pshift+1)]));
        COMPLEXASSIGN(This[k+i], creal(localcoeff[i]), cimag(localcoeff[i]));
      }
    }
  }
}

//performs m2m interaction using horner based scheme, with pre and post scaling
__global__ void shift_m2m_cuda_horner_scaled_t(int parentlevel,int parentstart,int pshift,
               dcmplx *This,const SORT_DCMPLX *z,int potshift,double rminshift  DEBUGVECTORSTRING)
{
  dcmplx tmp, invr;
  double rtmp, itmp;
  SORT_DCMPLX z0;
  __shared__ SORT_DCMPLX z0base[(SHIFT_M2M_MAXTHREADS+3)/4];
  int i, j, k;
  for(int boxnr=blockIdx.x*blockDim.x/8;boxnr<(1<<(parentlevel<<1));boxnr+=gridDim.x*blockDim.x/8) { //loop over all boxes. One block takes care of blockDim.x/4 boxes
//     localcoeff=(dcmplx*)&localcoeffs[(threadIdx.x>>1)*(pshift+1)]; //two threads for each set of coefficients
    if(threadIdx.x<blockDim.x/8&&threadIdx.x<(1<<(parentlevel<<1))) { /*only works if sizeof(dcmplx)/sizeof(float)<=4, which is the case for both double and float*/
      z0base[threadIdx.x]=z[boxnr+parentstart+threadIdx.x];
    }
    __syncthreads();
    if(boxnr+(threadIdx.x>>3)<(1<<(parentlevel<<1))) { /*safety check. For standard parameters, this should only be an issue for parentlevel = 1 - 3, then, it should be a multiple of blockDim/4*/
      z0=z[(boxnr<<2)+parentstart+(1<<(parentlevel<<1))+threadIdx.x/2];
      //calculate distances
      COMPLEXASSIGN(z0, creal(z0) - creal(z0base[threadIdx.x>>3]), cimag(z0) - cimag(z0base[threadIdx.x>>3])); //be careful with this use of COMPLEXASSIGN, remember that it is a macro
      k=((boxnr<<2)+parentstart+(1<<(parentlevel<<1))+threadIdx.x/2)*(pshift+1);
      for(i=(threadIdx.x&1);i<pshift+1;i+=2) {//read coefficient parameters
        fetch_double2((double2*)&wksp[i*blockDim.x+(threadIdx.x&(~1))],tcoeff1d,k+i);
      }
      //safe calculation of 1/z0
//       double div1=creal(z0);
//       double div2=cimag(z0);
//       int comparison=fabs(div1)<fabs(div2);
//       if(comparison) {
//         itmp=div2;
//         div2=div1;
//         div1=itmp;
//       }
//       itmp=div2*(1.0/div1);
//       rtmp=div1+itmp*div2;
//       div1=1.0/rtmp;
//       div2=itmp*div1;
//       if(comparison) {
//         COMPLEXASSIGN(tmp,div2,-div1);
//       }
//       else {
//         COMPLEXASSIGN(tmp,div1,-div2);
//       }
      if(fabs(creal(z0))>fabs(cimag(z0))) {
        itmp=cimag(z0)*(1.0/creal(z0));
        rtmp=creal(z0)+itmp*cimag(z0);
        COMPLEXASSIGN(tmp, 1.0/rtmp, -itmp*(1.0/rtmp));
      }
      else {
        itmp=creal(z0)*(1.0/cimag(z0));
        rtmp=cimag(z0)+itmp*creal(z0);
        COMPLEXASSIGN(tmp, itmp*(1.0/rtmp), -1.0/rtmp);
      }
      COMPLEXASSIGN(invr, creal(tmp)*creal(tmp)-cimag(tmp)*cimag(tmp), creal(tmp)*cimag(tmp)*2);
      if(fabs(creal(z0))<rminshift&&fabs(cimag(z0))<rminshift) { //if overflow is possible, run safe code, i.e. traditional horner shift
        if(threadIdx.x&1) {
          for(k=pshift;k>=2;k--) { //perform the Horner scheme
            for(i=k;i<=pshift;i++) {

              wksp[i*blockDim.x+threadIdx.x-1]+=creal(z0)*wksp[(i-1)*blockDim.x+threadIdx.x-1]-cimag(z0)*wksp[(i-1)*blockDim.x+threadIdx.x];
              wksp[i*blockDim.x+threadIdx.x]+=creal(z0)*wksp[(i-1)*blockDim.x+threadIdx.x]+cimag(z0)*wksp[(i-1)*blockDim.x+threadIdx.x-1];
            }
          }
          if(potshift==0) {
//             tmp=z0;
            COMPLEXASSIGN(tmp,creal(z0),cimag(z0));
            wksp[blockDim.x+threadIdx.x-1]-=wksp[threadIdx.x-1]*creal(tmp)-wksp[threadIdx.x]*cimag(tmp);
            wksp[blockDim.x+threadIdx.x]-=wksp[threadIdx.x-1]*cimag(tmp)+wksp[threadIdx.x]*creal(tmp);
            for(i=2;i<=pshift;i++) {
              double dtmp=creal(z0)*cimag(tmp)+cimag(z0)*creal(tmp); //depending on the definition of COMPLEXASSIGN, a poorly designed macro will overwrite tmp at first step, and use new value in second stage (normally not the case, but this is to prevent future problems)
              COMPLEXASSIGN(tmp, creal(tmp)*creal(z0)-cimag(tmp)*cimag(z0), dtmp);
              wksp[i*blockDim.x+threadIdx.x-1]-=wksp[threadIdx.x-1]*creal(tmp)-wksp[threadIdx.x]*cimag(tmp)*(1.0/i);
              wksp[i*blockDim.x+threadIdx.x]-=wksp[threadIdx.x-1]*cimag(tmp)+wksp[threadIdx.x]*creal(tmp)*(1.0/i);
            }
          }
        }
      }
      else {
        //perform scaling
        if((threadIdx.x&1)==1) { //set up starting values for the two threads working on the same set of data
          itmp=wksp[blockDim.x+threadIdx.x-1]*cimag(tmp)+wksp[blockDim.x+threadIdx.x]*creal(tmp);
          wksp[blockDim.x+threadIdx.x-1]=wksp[blockDim.x+threadIdx.x-1]*creal(tmp)-wksp[blockDim.x+threadIdx.x]*cimag(tmp);
          wksp[blockDim.x+threadIdx.x]=itmp;
        }
        else {
          COMPLEXASSIGN(tmp, 1, 0);
        }
        for(k=2+(threadIdx.x&1);k<=pshift;k+=2) {
          itmp=creal(invr)*cimag(tmp)+cimag(invr)*creal(tmp);
          COMPLEXASSIGN(tmp, creal(invr)*creal(tmp)-cimag(invr)*cimag(tmp), itmp);
          itmp=wksp[k*blockDim.x+(threadIdx.x&(~1))]*cimag(tmp)+wksp[k*blockDim.x+(threadIdx.x&(~1))+1]*creal(tmp);
          wksp[k*blockDim.x+(threadIdx.x&(~1))]=wksp[k*blockDim.x+(threadIdx.x&(~1))]*creal(tmp)-wksp[k*blockDim.x+(threadIdx.x&(~1))+1]*cimag(tmp);
          wksp[k*blockDim.x+(threadIdx.x&(~1))+1]=itmp;
        }
//         for(k=pshift;k>=2;k--) { //perform the Horner scheme
//           for(i=k;i<=pshift;i++) {
// //                    COMPLEXADD(localcoeff[i], creal(localcoeff[i-1]), cimag(localcoeff[i-1]));
//             ((double*)&localcoeff[i])[threadIdx.x&1]+=((double*)&localcoeff[i-1])[threadIdx.x&1];
//           }
//         }
        k=pshift;
        double a[UNROLLLEVEL+1];
        #if UNROLLLEVEL>=7
        for (; k>7; k-=7) {
          INNERUNROLLBASE8(wksp,0,STOREMDNM2MT,LOADM1M2MT,ADDEDN,1,k)
        }
        if (k>6)
        #elif UNROLLLEVEL>=6
        for (;k>6;)
        #endif
        #if UNROLLLEVEL>=6
        {
          INNERUNROLLBASE7(wksp,0,STOREMDNM2MT,LOADM1M2MT,ADDEDN,1,k)
          k-=6;
        }
        if (k>5)
        #elif UNROLLLEVEL>=5
        for (;k>5;)
        #endif
        #if UNROLLLEVEL>=5
        {
          INNERUNROLLBASE6(wksp,0,STOREMDNM2MT,LOADM1M2MT,ADDEDN,1,k)
          k-=5;
        }
        if (k>4)
        #elif UNROLLLEVEL>=4
        for (;k>4;)
        #endif
        #if UNROLLLEVEL>=4
        {
          INNERUNROLLBASE5(wksp,0,STOREMDNM2MT,LOADM1M2MT,ADDEDN,1,k)
          k-=4;
        }
        if (k>3)
        #elif UNROLLLEVEL>=3
        for (;k>3;)
        #endif
        #if UNROLLLEVEL>=3
        {
          INNERUNROLLBASE4(wksp,0,STOREMDNM2MT,LOADM1M2MT,ADDEDN,1,k)
          k-=3;
        }
        if (k>2)
        #elif UNROLLLEVEL>=2
        for (;k>2;)
        #endif
        #if UNROLLLEVEL>=2
        {
          INNERUNROLLBASE3(wksp,0,STOREMDNM2MT,LOADM1M2MT,ADDEDN,1,k)
          k-=2;
        }
        if (k > 1)
        #else
        for (;k > 1;k--)
        #endif
        {
          INNERUNROLLBASE2(wksp,0,STOREMDNM2MT,LOADM1M2MT,ADDEDN,1,k)
        }
        
        if(potshift==0) { //depending on logarithmic potential or not, an additional term may be necesarry
          //perform scaling
//           tmp=z0;
          COMPLEXASSIGN(tmp,creal(z0),cimag(z0));
          COMPLEXASSIGN(invr, creal(tmp)*creal(tmp)-cimag(tmp)*cimag(tmp), creal(tmp)*cimag(tmp)*2);
          if((threadIdx.x&1)==1) {
            wksp[blockDim.x+threadIdx.x-1]-=wksp[threadIdx.x-1];
            wksp[blockDim.x+threadIdx.x]-=wksp[threadIdx.x];
            itmp=wksp[blockDim.x+threadIdx.x-1]*cimag(tmp)+wksp[blockDim.x+threadIdx.x]*creal(tmp);
            wksp[blockDim.x+threadIdx.x-1]=wksp[blockDim.x+threadIdx.x-1]*creal(tmp)-wksp[blockDim.x+threadIdx.x]*cimag(tmp);
            wksp[blockDim.x+threadIdx.x]=itmp;
          }
          else {
            COMPLEXASSIGN(tmp, 1, 0);
          }
          for(k=2+(threadIdx.x&1);k<=pshift;k+=2) {
            itmp=creal(invr)*cimag(tmp)+cimag(invr)*creal(tmp); //depending on the definition of COMPLEXASSIGN, a poorly designed macro will overwrite tmp at first step, and use new value in second stage (normally not the case, but this is to prevent future problems)
            COMPLEXASSIGN(tmp, creal(tmp)*creal(invr)-cimag(tmp)*cimag(invr), itmp);
            wksp[k*blockDim.x+(threadIdx.x&(~1))]-=wksp[(threadIdx.x&(~1))]*(1.0/k);
            wksp[k*blockDim.x+(threadIdx.x&(~1))+1]-=wksp[(threadIdx.x&(~1))+1]*(1.0/k);
            itmp=wksp[k*blockDim.x+(threadIdx.x&(~1))]*cimag(tmp)+wksp[k*blockDim.x+(threadIdx.x&(~1))+1]*creal(tmp);
            wksp[k*blockDim.x+(threadIdx.x&(~1))]=wksp[k*blockDim.x+(threadIdx.x&(~1))]*creal(tmp)-wksp[k*blockDim.x+(threadIdx.x&(~1))+1]*cimag(tmp);
            wksp[k*blockDim.x+(threadIdx.x&(~1))+1]=itmp;
          }
        }
        else {
          //perform scaling
//           tmp=z0;
          COMPLEXASSIGN(tmp,creal(z0),cimag(z0));
          COMPLEXASSIGN(invr, creal(tmp)*creal(tmp)-cimag(tmp)*cimag(tmp), creal(tmp)*cimag(tmp)*2);
          if((threadIdx.x&1)==1) {
            itmp=wksp[blockDim.x+(threadIdx.x&(~1))]*cimag(tmp)+wksp[blockDim.x+(threadIdx.x&(~1))+1]*creal(tmp);
            wksp[blockDim.x+(threadIdx.x&(~1))]=wksp[blockDim.x+(threadIdx.x&(~1))]*creal(tmp)-wksp[blockDim.x+(threadIdx.x&(~1))+1]*cimag(tmp);
            wksp[blockDim.x+(threadIdx.x&(~1))+1]=itmp;
          }
          else {
            COMPLEXASSIGN(tmp, 1, 0);
          }
          for(k=2+(threadIdx.x&1);k<=pshift;k+=2) {
            itmp=creal(invr)*cimag(tmp)+cimag(invr)*creal(tmp); //depending on the definition of COMPLEXASSIGN, a poorly designed macro will overwrite tmp at first step, and use new value in second stage (normally not the case, but this is to prevent future problems)
            COMPLEXASSIGN(tmp, creal(tmp)*creal(invr)-cimag(tmp)*cimag(invr), itmp);
            itmp=wksp[k*blockDim.x+(threadIdx.x&(~1))]*cimag(tmp)+wksp[k*blockDim.x+(threadIdx.x&(~1))+1]*creal(tmp);
            wksp[k*blockDim.x+(threadIdx.x&(~1))]=wksp[k*blockDim.x+(threadIdx.x&(~1))]*creal(tmp)-wksp[k*blockDim.x+(threadIdx.x&(~1))+1]*cimag(tmp);
            wksp[k*blockDim.x+(threadIdx.x&(~1))+1]=itmp;
          }
        }
      }
    }
    __syncthreads();
    
    //calculate the sum of the for contributing boxes
    if(boxnr+(threadIdx.x>>3)<(1<<(parentlevel<<1))) { /*safety check. For standard parameters, this should only be an issue for parentlevel=1 and 2, then, it should be a multiple of blockDim/4*/
      //write results, here, each thread has its own box to write results to
      for(i=(threadIdx.x&2)/2;i<=pshift;i+=2) {
        wksp[i*blockDim.x+(threadIdx.x&(~2))]+=wksp[i*blockDim.x+(threadIdx.x&(~2))+2];
      }
    }
    __syncthreads();
    //second part of the sum, and write results back to array
    if(boxnr+(threadIdx.x>>3)<(1<<(parentlevel<<1))) { /*safety check. For standard parameters, this should only be an issue for parentlevel=1 and 2, then, it should be a multiple of blockDim/4*/
      //write results, here, each thread has its own box to write results to
      k=(boxnr+(threadIdx.x>>3)+parentstart)*(pshift+1);
      for(i=(threadIdx.x&6)/2;i<=pshift;i+=4) {
        ((double*)&This[k+i])[threadIdx.x&1]=wksp[i*blockDim.x+(threadIdx.x&(~6))]+wksp[i*blockDim.x+(threadIdx.x&(~6))+4];
      }
    }
  }
}
#endif //MULTIPOLEINITCUDA


#ifdef MULTIPOLESHIFTCUDA

//number of threads should be a multiple of 4 (4 threads per interaction)

//corresponds to shift_p2p in the CPU code, Moves potential expansions by using the Horner scheme
__global__ void shift_p2p_cuda(int parentlevel,int parentstart,int pshift,
               dcmplx *This  DEBUGVECTORSTRING)
{
  dcmplx *localcoeff=(dcmplx*)&localcoeffs[threadIdx.x*(pshift+1)];
  __shared__ SORT_DCMPLX z0[P2PMAXTHREADS];
  __shared__ SORT_DCMPLX z0base[(P2PMAXTHREADS+3)/4];
  int i, j;
  for(int boxnr=blockIdx.x*blockDim.x/4;boxnr<(1<<(parentlevel<<1));boxnr+=gridDim.x*blockDim.x/4) { //loop over all boxes. One block takes care of blockDim.x/4 boxes
    if(threadIdx.x<blockDim.x/4*sizeof(SORT_DCMPLX)/sizeof(float)) { /*only works if sizeof(dcmplx)/sizeof(float)<=4, which is the case for both double and float*/
      ((float*)z0base)[threadIdx.x] =tex1Dfetch(tz0, (boxnr+parentstart)*sizeof(SORT_DCMPLX)/sizeof(float)+threadIdx.x); /*not really a problem if values outside the list would be fetched, since they are 0*/
    }
    __syncthreads();
    if(boxnr+(threadIdx.x>>2)<(1<<(parentlevel<<1))) { /*safety check. For standard parameters, this should only be an issue for parentlevel=1 and 2, then, it should be a multiple of blockDim/4*/
      for(i=0;i<sizeof(dcmplx)/sizeof(float);i++) //read position parameters
        ((float*)&z0[threadIdx.x])[i] =tex1Dfetch(tz0, ((boxnr<<2)+parentstart+(1<<(parentlevel<<1))+threadIdx.x)*sizeof(dcmplx)/sizeof(float)+i);
      //calculate distances
      COMPLEXASSIGN(z0[threadIdx.x], creal(z0base[threadIdx.x>>2])-creal(z0[threadIdx.x]), cimag(z0base[threadIdx.x>>2])-cimag(z0[threadIdx.x]));
      j=(boxnr+(threadIdx.x>>2)+parentstart)*(pshift+1);
      for(i=0;i<pshift+1;i++) //read coefficient parameters
        COMPLEXASSIGN(localcoeff[i], creal(This[j+i]), cimag(This[j+i]));
      for(j=0;j<=pshift;j++) { //perform the Horner scheme
        for(i=pshift-j;i<pshift;i++) {
          COMPLEXSUB(localcoeff[i], creal(z0[threadIdx.x])*creal(localcoeff[i+1])-cimag(z0[threadIdx.x])*cimag(localcoeff[i+1]), creal(z0[threadIdx.x])*cimag(localcoeff[i+1])+cimag(z0[threadIdx.x])*creal(localcoeff[i+1]));
        }
      }
      //write results, here, each thread has its own box to write results to
      j=((boxnr<<2)+parentstart+(1<<(parentlevel<<1))+threadIdx.x)*(pshift+1);
      for(i=0;i<=pshift;i++) {
        COMPLEXADD(This[j+i], creal(localcoeff[i]), cimag(localcoeff[i]));
      }
    }
  }
}
//corresponds to shift_p2p in the CPU code, Moves potential expansions by using the Horner scheme with pre and post scaling
__global__ void shift_p2p_cuda_horner_scaled(int parentlevel,int parentstart,int pshift,
               dcmplx *This,const SORT_DCMPLX *z,double rminshift  DEBUGVECTORSTRING)
{
  dcmplx *localcoeff;
  dcmplx tmp, invr;
  double rtmp, itmp;
  SORT_DCMPLX z0;
  __shared__ SORT_DCMPLX z0base[(SHIFT_M2M_MAXTHREADS+3)/4];
  int i, j, k;
  for(int boxnr=blockIdx.x*blockDim.x/8;boxnr<(1<<(parentlevel<<1));boxnr+=gridDim.x*blockDim.x/8) { //loop over all boxes. One block takes care of blockDim.x/2 boxes
    localcoeff=(dcmplx*)&localcoeffs[(threadIdx.x>>1)*(pshift+1)];
    if(threadIdx.x<blockDim.x/8&&threadIdx.x<(1<<(parentlevel<<1))) {
      z0base[threadIdx.x]=z[boxnr+parentstart+threadIdx.x];
    }
    __syncthreads();
    if(boxnr+(threadIdx.x>>3)<(1<<(parentlevel<<1))) { /*safety check. For standard parameters, this should only be an issue for parentlevel = 1 - 3, then, it should be a multiple of blockDim/4*/
      z0=z[(boxnr<<2)+parentstart+(1<<(parentlevel<<1))+threadIdx.x/2];
      //calculate distances
      COMPLEXASSIGN(z0, creal(z0base[threadIdx.x>>3]) - creal(z0), cimag(z0base[threadIdx.x>>3]) - cimag(z0)); //be careful with this use of COMPLEXASSIGN, remember that it is a macro
      k=(boxnr+(threadIdx.x>>3)+parentstart)*(pshift+1);
      for(i=(threadIdx.x&1);i<pshift+1;i+=2) //read coefficient parameters
        COMPLEXASSIGN(localcoeff[i], creal(This[k+i]), cimag(This[k+i]));
//       tmp=z0;
      COMPLEXASSIGN(tmp,creal(z0),cimag(z0));
      COMPLEXASSIGN(invr, creal(tmp)*creal(tmp)-cimag(tmp)*cimag(tmp), creal(tmp)*cimag(tmp)*2);
      if(fabs(creal(z0))<rminshift&&fabs(cimag(z0))<rminshift) { //if overflow is possible, run safe code i.e. traditional horner shift
        if(threadIdx.x&1) {
          for(k=0;k<=pshift;k++) { //perform the Horner scheme
            for(i=pshift-k;i<pshift;i++) {
              COMPLEXSUB(localcoeff[i], creal(z0)*creal(localcoeff[i+1])-cimag(z0)*cimag(localcoeff[i+1]), creal(z0)*cimag(localcoeff[i+1])+cimag(z0)*creal(localcoeff[i+1]));
            }
          }
        }
      }
      else { //use scaling
        if((threadIdx.x&1)==1) { //initiate one thread to z0^2, the other to z0
          tmp=invr;
        }
        
        //perform the scaling
        itmp=creal(localcoeff[(threadIdx.x&1)])*cimag(tmp)+cimag(localcoeff[(threadIdx.x&1)])*creal(tmp);
        COMPLEXASSIGN(localcoeff[(threadIdx.x&1)], creal(localcoeff[(threadIdx.x&1)])*creal(tmp)-cimag(localcoeff[(threadIdx.x&1)])*cimag(tmp), itmp);
        for(k=2+(threadIdx.x&1);k<=pshift;k+=2) {
          itmp=creal(invr)*cimag(tmp)+cimag(invr)*creal(tmp);
          COMPLEXASSIGN(tmp, creal(invr)*creal(tmp)-cimag(invr)*cimag(tmp), itmp);
          itmp=creal(localcoeff[k])*cimag(tmp)+cimag(localcoeff[k])*creal(tmp);
          COMPLEXASSIGN(localcoeff[k], creal(localcoeff[k])*creal(tmp)-cimag(localcoeff[k])*cimag(tmp), itmp);
        }
//         for(k=0;k<=pshift;k++) { //perform the Horner scheme
//           for(i=pshift-k;i<pshift;i++) {
//             ((double*)&localcoeff[i])[threadIdx.x&1]-=((double*)&localcoeff[i+1])[threadIdx.x&1];
//           }
//         }
        k=1;
        double a[UNROLLLEVEL+1];
        #if UNROLLLEVEL>=7 //loop is performed on different levels depending on UNROLLLEVEL
        for (; k <= pshift-6; k+=7) {
          INNERUNROLLBASE8(wksp,0,STOREMUPM2M,LOADM1M2M,SUBEUP,0,pshift-k) //here, imaginary part is taken directly after real part for each operation
        }
        if (k <= pshift-5)
        #elif UNROLLLEVEL>=6 //The pattern repeats itself for each unroll level
        for (;k <= pshift-5;)
        #endif
        #if UNROLLLEVEL>=6
        {
          INNERUNROLLBASE7(wksp,0,STOREMUPM2M,LOADM1M2M,SUBEUP,0,pshift-k)
          k+=6;
        }
        if (k <= pshift-4)
        #elif UNROLLLEVEL>=5
        for (;k <= pshift-4;)
        #endif
        #if UNROLLLEVEL>=5
        {
          INNERUNROLLBASE6(wksp,0,STOREMUPM2M,LOADM1M2M,SUBEUP,0,pshift-k)
          k+=5;
        }
        if (k <= pshift-3)
        #elif UNROLLLEVEL>=4
        for (;k <= pshift-3;)
        #endif
        #if UNROLLLEVEL>=4
        {
          INNERUNROLLBASE5(wksp,0,STOREMUPM2M,LOADM1M2M,SUBEUP,0,pshift-k)
          k+=4;
        }
        if (k <= pshift-2)
        #elif UNROLLLEVEL>=3
        for (;k <= pshift-2;)
        #endif
        #if UNROLLLEVEL>=3
        {
          INNERUNROLLBASE4(wksp,0,STOREMUPM2M,LOADM1M2M,SUBEUP,0,pshift-k)
          k+=3;
        }
        if (k <= pshift-1)
        #elif UNROLLLEVEL>=2
        for (;k <= pshift-1;)
        #endif
        #if UNROLLLEVEL>=2
        {
          INNERUNROLLBASE3(wksp,0,STOREMUPM2M,LOADM1M2M,SUBEUP,0,pshift-k)
          k+=2;
        }
        if (k <= pshift)
        #else
        for (;k <= pshift;k++)
        #endif
        {
          INNERUNROLLBASE2(wksp,0,STOREMUPM2M,LOADM1M2M,SUBEUP,0,pshift-k)
        }
        //safe calculation of 1/z0
        if(fabs(creal(z0))>fabs(cimag(z0))) {
          itmp=cimag(z0)*(1.0/creal(z0));
          rtmp=creal(z0)+itmp*cimag(z0);
          COMPLEXASSIGN(tmp, 1.0/rtmp, -itmp*(1.0/rtmp));
        }
        else {
          itmp=creal(z0)*(1.0/cimag(z0));
          rtmp=cimag(z0)+itmp*creal(z0);
          COMPLEXASSIGN(tmp, itmp*(1.0/rtmp), -1.0/rtmp);
        }
        COMPLEXASSIGN(invr, creal(tmp)*creal(tmp)-cimag(tmp)*cimag(tmp), creal(tmp)*cimag(tmp)*2);
        if((threadIdx.x&1)==1) { //initiate one thread to z0^-2, the other to z0^-1
          tmp=invr;
        }
        //perform the scaling
        itmp=creal(localcoeff[(threadIdx.x&1)])*cimag(tmp)+cimag(localcoeff[(threadIdx.x&1)])*creal(tmp);
        COMPLEXASSIGN(localcoeff[(threadIdx.x&1)], creal(localcoeff[(threadIdx.x&1)])*creal(tmp)-cimag(localcoeff[(threadIdx.x&1)])*cimag(tmp), itmp);
        for(k=2+(threadIdx.x&1);k<=pshift;k+=2) {
          itmp=creal(invr)*cimag(tmp)+cimag(invr)*creal(tmp); //depending on the definition of COMPLEXASSIGN, a poorly designed macro will overwrite tmp at first step, and use new value in second stage (normally not the case, but this is to prevent future problems)
          COMPLEXASSIGN(tmp, creal(tmp)*creal(invr)-cimag(tmp)*cimag(invr), itmp);
          itmp=creal(localcoeff[k])*cimag(tmp)+cimag(localcoeff[k])*creal(tmp);
          COMPLEXASSIGN(localcoeff[k], creal(localcoeff[k])*creal(tmp)-cimag(localcoeff[k])*cimag(tmp), itmp);
        }
      }
    }
    __syncthreads();
    //write results
    if(boxnr+(threadIdx.x>>3)<(1<<(parentlevel<<1))) { /*safety check. For standard parameters, this should only be an issue for parentlevel=1 and 2, then, it should be a multiple of blockDim/4*/
      //write results, here, each thread has its own box to write results to
      localcoeff=(dcmplx*)&localcoeffs[(threadIdx.x>>1)*(pshift+1)];
      k=((boxnr<<2)+parentstart+(1<<(parentlevel<<1))+threadIdx.x/2)*(pshift+1);
      for(i=(threadIdx.x&1);i<=pshift;i+=2) {
        COMPLEXADD(This[k+i], creal(localcoeff[i]), cimag(localcoeff[i]));
      }
    }
  }
}

//corresponds to shift_p2p in the CPU code, Moves potential expansions by using the Horner scheme with pre and post scaling
__global__ void shift_p2p_cuda_horner_scaled_t(int parentlevel,int parentstart,int pshift,
               dcmplx *This,const SORT_DCMPLX *z,double rminshift  DEBUGVECTORSTRING)
{
//   dcmplx *localcoeff;
  dcmplx tmp, invr;
  double rtmp, itmp;
  SORT_DCMPLX z0;
  __shared__ SORT_DCMPLX z0base[(SHIFT_M2M_MAXTHREADS+3)/4];
  int i, j, k;
  for(int boxnr=blockIdx.x*blockDim.x/8;boxnr<(1<<(parentlevel<<1));boxnr+=gridDim.x*blockDim.x/8) { //loop over all boxes. One block takes care of blockDim.x/2 boxes
//     localcoeff=(dcmplx*)&localcoeffs[(threadIdx.x>>1)*(pshift+1)];
    if(threadIdx.x<blockDim.x/8&&threadIdx.x<(1<<(parentlevel<<1))) {
      z0base[threadIdx.x]=z[boxnr+parentstart+threadIdx.x];
    }
    __syncthreads();
    if(boxnr+(threadIdx.x>>3)<(1<<(parentlevel<<1))) { /*safety check. For standard parameters, this should only be an issue for parentlevel = 1 - 3, then, it should be a multiple of blockDim/4*/
      z0=z[(boxnr<<2)+parentstart+(1<<(parentlevel<<1))+threadIdx.x/2];
      //calculate distances
      COMPLEXASSIGN(z0, creal(z0base[threadIdx.x>>3]) - creal(z0), cimag(z0base[threadIdx.x>>3]) - cimag(z0)); //be careful with this use of COMPLEXASSIGN, remember that it is a macro
      k=(boxnr+(threadIdx.x>>3)+parentstart)*(pshift+1);
      for(i=(threadIdx.x&1);i<pshift+1;i+=2) { //read coefficient parameters
        fetch_double2((double2*)&wksp[i*blockDim.x+(threadIdx.x&(~1))],tcoeff2d,k+i);
      }
//       tmp=z0;
      COMPLEXASSIGN(tmp,creal(z0),cimag(z0));
      COMPLEXASSIGN(invr, creal(tmp)*creal(tmp)-cimag(tmp)*cimag(tmp), creal(tmp)*cimag(tmp)*2);
      if(fabs(creal(z0))<rminshift&&fabs(cimag(z0))<rminshift) { //if overflow is possible, run safe code i.e. traditional horner shift
        if(threadIdx.x&1) {
          for(k=0;k<=pshift;k++) { //perform the Horner scheme
            for(i=pshift-k;i<pshift;i++) {
              wksp[i*blockDim.x+threadIdx.x-1]-=creal(z0)*wksp[(i+1)*blockDim.x+threadIdx.x-1]-cimag(z0)*wksp[(i+1)*blockDim.x+threadIdx.x];
              wksp[i*blockDim.x+threadIdx.x]-=creal(z0)*wksp[(i+1)*blockDim.x+threadIdx.x]+cimag(z0)*wksp[(i+1)*blockDim.x+threadIdx.x-1];
            }
          }
        }
      }
      else { //use scaling
        if((threadIdx.x&1)==1) { //initiate one thread to z0^2, the other to z0
          tmp=invr;
        }
        
        //perform the scaling
        itmp=wksp[(threadIdx.x&1)*blockDim.x+(threadIdx.x&(~1))]*cimag(tmp)+wksp[(threadIdx.x&1)*blockDim.x+(threadIdx.x&(~1))+1]*creal(tmp);
        wksp[(threadIdx.x&1)*blockDim.x+(threadIdx.x&(~1))]=wksp[(threadIdx.x&1)*blockDim.x+(threadIdx.x&(~1))]*creal(tmp)-wksp[(threadIdx.x&1)*blockDim.x+(threadIdx.x&(~1))+1]*cimag(tmp);
        wksp[(threadIdx.x&1)*blockDim.x+(threadIdx.x&(~1))+1]=itmp;
        
        for(k=2+(threadIdx.x&1);k<=pshift;k+=2) {
          itmp=creal(invr)*cimag(tmp)+cimag(invr)*creal(tmp);
          COMPLEXASSIGN(tmp, creal(invr)*creal(tmp)-cimag(invr)*cimag(tmp), itmp);
          itmp=wksp[k*blockDim.x+(threadIdx.x&(~1))]*cimag(tmp)+wksp[k*blockDim.x+(threadIdx.x&(~1))+1]*creal(tmp);
          wksp[k*blockDim.x+(threadIdx.x&(~1))]=wksp[k*blockDim.x+(threadIdx.x&(~1))]*creal(tmp)-wksp[k*blockDim.x+(threadIdx.x&(~1))+1]*cimag(tmp);
          wksp[k*blockDim.x+(threadIdx.x&(~1))+1]=itmp;
        }
//         for(k=0;k<=pshift;k++) { //perform the Horner scheme
//           for(i=pshift-k;i<pshift;i++) {
//             ((double*)&localcoeff[i])[threadIdx.x&1]-=((double*)&localcoeff[i+1])[threadIdx.x&1];
//           }
//         }
        k=1;
        double a[UNROLLLEVEL+1];
        #if UNROLLLEVEL>=7 //loop is performed on different levels depending on UNROLLLEVEL
        for (; k <= pshift-6; k+=7) {
          INNERUNROLLBASE8(wksp,0,STOREMUPM2MT,LOADM1M2MT,SUBEUP,0,pshift-k) //here, imaginary part is taken directly after real part for each operation
        }
        if (k <= pshift-5)
        #elif UNROLLLEVEL>=6 //The pattern repeats itself for each unroll level
        for (;k <= pshift-5;)
        #endif
        #if UNROLLLEVEL>=6
        {
          INNERUNROLLBASE7(wksp,0,STOREMUPM2MT,LOADM1M2MT,SUBEUP,0,pshift-k)
          k+=6;
        }
        if (k <= pshift-4)
        #elif UNROLLLEVEL>=5
        for (;k <= pshift-4;)
        #endif
        #if UNROLLLEVEL>=5
        {
          INNERUNROLLBASE6(wksp,0,STOREMUPM2MT,LOADM1M2MT,SUBEUP,0,pshift-k)
          k+=5;
        }
        if (k <= pshift-3)
        #elif UNROLLLEVEL>=4
        for (;k <= pshift-3;)
        #endif
        #if UNROLLLEVEL>=4
        {
          INNERUNROLLBASE5(wksp,0,STOREMUPM2MT,LOADM1M2MT,SUBEUP,0,pshift-k)
          k+=4;
        }
        if (k <= pshift-2)
        #elif UNROLLLEVEL>=3
        for (;k <= pshift-2;)
        #endif
        #if UNROLLLEVEL>=3
        {
          INNERUNROLLBASE4(wksp,0,STOREMUPM2MT,LOADM1M2MT,SUBEUP,0,pshift-k)
          k+=3;
        }
        if (k <= pshift-1)
        #elif UNROLLLEVEL>=2
        for (;k <= pshift-1;)
        #endif
        #if UNROLLLEVEL>=2
        {
          INNERUNROLLBASE3(wksp,0,STOREMUPM2MT,LOADM1M2MT,SUBEUP,0,pshift-k)
          k+=2;
        }
        if (k <= pshift)
        #else
        for (;k <= pshift;k++)
        #endif
        {
          INNERUNROLLBASE2(wksp,0,STOREMUPM2MT,LOADM1M2MT,SUBEUP,0,pshift-k)
        }
        //safe calculation of 1/z0
        if(fabs(creal(z0))>fabs(cimag(z0))) {
          itmp=cimag(z0)*(1.0/creal(z0));
          rtmp=creal(z0)+itmp*cimag(z0);
          COMPLEXASSIGN(tmp, 1.0/rtmp, -itmp*(1.0/rtmp));
        }
        else {
          itmp=creal(z0)*(1.0/cimag(z0));
          rtmp=cimag(z0)+itmp*creal(z0);
          COMPLEXASSIGN(tmp, itmp*(1.0/rtmp), -1.0/rtmp);
        }
        COMPLEXASSIGN(invr, creal(tmp)*creal(tmp)-cimag(tmp)*cimag(tmp), creal(tmp)*cimag(tmp)*2);
        if((threadIdx.x&1)==1) { //initiate one thread to z0^-2, the other to z0^-1
          tmp=invr;
        }
        //perform the scaling
        itmp=wksp[(threadIdx.x&1)*blockDim.x+(threadIdx.x&(~1))]*cimag(tmp)+wksp[(threadIdx.x&1)*blockDim.x+(threadIdx.x&(~1))+1]*creal(tmp);
        wksp[(threadIdx.x&1)*blockDim.x+(threadIdx.x&(~1))]=wksp[(threadIdx.x&1)*blockDim.x+(threadIdx.x&(~1))]*creal(tmp)-wksp[(threadIdx.x&1)*blockDim.x+(threadIdx.x&(~1))+1]*cimag(tmp);
        wksp[(threadIdx.x&1)*blockDim.x+(threadIdx.x&(~1))+1]=itmp;
        
        for(k=2+(threadIdx.x&1);k<=pshift;k+=2) {
          itmp=creal(invr)*cimag(tmp)+cimag(invr)*creal(tmp);
          COMPLEXASSIGN(tmp, creal(invr)*creal(tmp)-cimag(invr)*cimag(tmp), itmp);
          itmp=wksp[k*blockDim.x+(threadIdx.x&(~1))]*cimag(tmp)+wksp[k*blockDim.x+(threadIdx.x&(~1))+1]*creal(tmp);
          wksp[k*blockDim.x+(threadIdx.x&(~1))]=wksp[k*blockDim.x+(threadIdx.x&(~1))]*creal(tmp)-wksp[k*blockDim.x+(threadIdx.x&(~1))+1]*cimag(tmp);
          wksp[k*blockDim.x+(threadIdx.x&(~1))+1]=itmp;
        }
      }
    }
    __syncthreads();
    //write results
    if(boxnr+(threadIdx.x>>3)<(1<<(parentlevel<<1))) { /*safety check. For standard parameters, this should only be an issue for parentlevel=1 and 2, then, it should be a multiple of blockDim/4*/
      //write results, here, each thread has its own box to write results to
      k=((boxnr<<2)+parentstart+(1<<(parentlevel<<1))+threadIdx.x/2)*(pshift+1);
      for(i=(threadIdx.x&1);i<=pshift;i+=2) {
        COMPLEXADD(This[k+i],  wksp[i*blockDim.x+(threadIdx.x&(~1))], wksp[i*blockDim.x+(threadIdx.x&(~1))+1]);
      }
    }
  }
}

//For speed issues, this one is written as a template with #ifdef checks for if it is logarithmic potential or not
#undef LOGM2PS
__global__ void __launch_bounds__(32) shift_m2ps_cuda_assym(int startlevel, int nlevel, int pshift,
               dcmplx *coeff1, dcmplx *coeff2,SORT_DCMPLX *z0, int** jcptr2, int** kcptr2, int** ir2  DEBUGVECTORSTRING)
{
#include "cudam2psassym2.h"
}
#define LOGM2PS
__global__ void __launch_bounds__(32) shift_m2ps_cuda_assym_log(int startlevel, int nlevel, int pshift,
               dcmplx *coeff1, dcmplx *coeff2,SORT_DCMPLX *z0, int** jcptr2, int** kcptr2, int** ir2  DEBUGVECTORSTRING)
{
#include "cudam2psassym2.h"
}
#endif /*MULTIPOLESHIFTCUDA*/

//reorder functions. Four different versions exists, which takes 1, 2, 3 and for vectors for reordering
__global__ void reorder_cuda_1(double* output1, double* input1, int* index, int count) {
  int localindex;
  for(int i=blockDim.x * blockIdx.x + threadIdx.x;i<count;i+=blockDim.x*gridDim.x) {
    localindex=index[i];
    output1[i]=input1[localindex];
  }
}
__global__ void reorder_cuda_2(double* output1, double* output2, double* input1, double* input2, int* index, int count) {
  int localindex;
  for(int i=blockDim.x * blockIdx.x + threadIdx.x;i<count;i+=blockDim.x*gridDim.x) {
    localindex=index[i];
    output1[i]=input1[localindex];
    output2[i]=input2[localindex];
  }
}
__global__ void reorder_cuda_3(double* output1, double* output2, double* output3, double* input1, double* input2, double* input3, int* index, int count) {
  int localindex;
  for(int i=blockDim.x * blockIdx.x + threadIdx.x;i<count;i+=blockDim.x*gridDim.x) {
    localindex=index[i];
    output1[i]=input1[localindex];
    output2[i]=input2[localindex];
    output3[i]=input3[localindex];
  }
}
__global__ void reorder_cuda_4(double* output1, double* output2, double* output3, double* output4, double* input1, double* input2, double* input3, double* input4, int* index, int count) {
  int localindex;
  for(int i=blockDim.x * blockIdx.x + threadIdx.x;i<count;i+=blockDim.x*gridDim.x) {
    localindex=index[i];
    output1[i]=input1[localindex];
    output2[i]=input2[localindex];
    output3[i]=input3[localindex];
    output4[i]=input4[localindex];
  }
}

//as above, but reorders them in the opposite direction
__global__ void reorder_cuda_1_inv(double* output1, double* input1, int* index, int count) {
  int localindex;
  for(int i=blockDim.x * blockIdx.x + threadIdx.x;i<count;i+=blockDim.x*gridDim.x) {
    localindex=index[i];
    output1[localindex]=input1[i];
  }
}
__global__ void reorder_cuda_2_inv(double* output1, double* output2, double* input1, double* input2, int* index, int count) {
  int localindex;
  for(int i=blockDim.x * blockIdx.x + threadIdx.x;i<count;i+=blockDim.x*gridDim.x) {
    localindex=index[i];
    output1[localindex]=input1[i];
    output2[localindex]=input2[i];
  }
}
__global__ void reorder_cuda_3_inv(double* output1, double* output2, double* output3, double* input1, double* input2, double* input3, int* index, int count) {
  int localindex;
  for(int i=blockDim.x * blockIdx.x + threadIdx.x;i<count;i+=blockDim.x*gridDim.x) {
    localindex=index[i];
    output1[localindex]=input1[i];
    output2[localindex]=input2[i];
    output3[localindex]=input3[i];
  }
}
__global__ void reorder_cuda_4_inv(double* output1, double* output2, double* output3, double* output4, double* input1, double* input2, double* input3, double* input4, int* index, int count) {
  int localindex;
  for(int i=blockDim.x * blockIdx.x + threadIdx.x;i<count;i+=blockDim.x*gridDim.x) {
    localindex=index[i];
    output1[localindex]=input1[i];
    output2[localindex]=input2[i];
    output3[localindex]=input3[i];
    output4[localindex]=input4[i];
  }
}
