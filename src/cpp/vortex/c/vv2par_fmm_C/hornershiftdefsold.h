#ifdef UNROLLHORNEROLD
// #define CROSSORDER
#ifdef SSE2INTRINSIC
#define UPWARDADD(wkspm1,j,base) wkspm1[j]=_mm_add_pd(wkspm1[j],wkspm1[j+1])
#elif defined(CROSSORDER)
#define UPWARDADD(wksp1,j,base) wksp1[2*(j)+0] += wksp1[2*(j+1)+0];  wksp1[2*(j)+1] += wksp1[2*(j+1)+1];
#else
#define UPWARDADD(wksp1,j,base) wksp1[2*(j)+base] += wksp1[2*(j+1)+base]
#endif
#define UPWARDADDU2(wksp1,j,base) UPWARDADD(wksp1,j,base); UPWARDADD(wksp1,j-1,base)
#define UPWARDADDU3(wksp1,j,base) UPWARDADDU2(wksp1,j,base); UPWARDADD(wksp1,j-2,base)
#define UPWARDADDU4(wksp1,j,base) UPWARDADDU3(wksp1,j,base); UPWARDADD(wksp1,j-3,base)
#define UPWARDADDU5(wksp1,j,base) UPWARDADDU4(wksp1,j,base); UPWARDADD(wksp1,j-4,base)
#define UPWARDADDU6(wksp1,j,base) UPWARDADDU5(wksp1,j,base); UPWARDADD(wksp1,j-5,base)
#define UPWARDL0(UPWADD,wksp1,j,base)
#define UPWARDL1(UPWADD,wksp1,j,base) UPWADD(wksp1,j,base)
#define UPWARDL2(UPWADD,wksp1,j,base) UPWARDL1(UPWADD,wksp1,j,base); UPWADD(wksp1,j+1,base)
#define UPWARDL3(UPWADD,wksp1,j,base) UPWARDL2(UPWADD,wksp1,j,base); UPWADD(wksp1,j+2,base)
#define UPWARDL4(UPWADD,wksp1,j,base) UPWARDL3(UPWADD,wksp1,j,base); UPWADD(wksp1,j+3,base)
#define UPWARDL5(UPWADD,wksp1,j,base) UPWARDL4(UPWADD,wksp1,j,base); UPWADD(wksp1,j+4,base)
#define UPWARDL6(UPWADD,wksp1,j,base) UPWARDL5(UPWADD,wksp1,j,base); UPWADD(wksp1,j+5,base)
#define UPWARDL7(UPWADD,wksp1,j,base) UPWARDL6(UPWADD,wksp1,j,base); UPWADD(wksp1,j+6,base)
#define UPWARDL8(UPWADD,wksp1,j,base) UPWARDL7(UPWADD,wksp1,j,base); UPWADD(wksp1,j+7,base)
#define UPWARDEND2(wksp1,j,base) UPWARDADD(wksp1,j,base)
#define UPWARDEND3(wksp1,j,base) UPWARDADDU2(wksp1,j,base); UPWARDEND2(wksp1,j,base)
#define UPWARDEND4(wksp1,j,base) UPWARDADDU3(wksp1,j,base); UPWARDEND3(wksp1,j,base)
#define UPWARDEND5(wksp1,j,base) UPWARDADDU4(wksp1,j,base); UPWARDEND4(wksp1,j,base)
#define UPWARDEND6(wksp1,j,base) UPWARDADDU5(wksp1,j,base); UPWARDEND5(wksp1,j,base)
#define UPSHIFTEND(UPWARDLX,UPWADD,UPWARDENDX,wksp1,j,base) UPWARDLX(UPWADD,wksp1,j,base); UPWARDENDX(wksp1,j,base)

#ifdef SSE2INTRINSIC
#define UPREALIMAG(UPWARDLX,UPWADD,wksp1,j) UPWARDLX(UPWADD,wkspm1,j,0)
#define UPREALIMAGEND(UPWARDLX,UPWADD,UPWARDENDX,wksp1,j,endj) UPWARDLX(UPWADD,wkspm1,j,0); UPWARDENDX(wkspm1,endj,0)
#define UPREALIMAGD(UPWARDLX,UPWADD,wksp1,j) UPWARDLX(UPWADD,wkspm1,j,0); UPWARDLX(UPWADD,wkspm2,j,0)
#define UPREALIMAGENDD(UPWARDLX,UPWADD,UPWARDENDX,wksp1,j,endj) UPWARDLX(UPWADD,wkspm1,j,0); UPWARDENDX(wkspm1,endj,0); UPWARDLX(UPWADD,wkspm2,j,0); UPWARDENDX(wkspm2,endj,0)
#else
#ifdef CROSSORDER
#define UPREALIMAG(UPWARDLX,UPWADD,wksp1,j) UPWARDLX(UPWADD,wksp1,j,0)
#define UPREALIMAGEND(UPWARDLX,UPWADD,UPWARDENDX,wksp1,j,endj) UPWARDLX(UPWADD,wksp1,j,0); UPWARDENDX(wksp1,endj,0)
#else
#define UPREALIMAG(UPWARDLX,UPWADD,wksp1,j) UPWARDLX(UPWADD,wksp1,j,0); UPWARDLX(UPWADD,wksp1,j,1)
#define UPREALIMAGEND(UPWARDLX,UPWADD,UPWARDENDX,wksp1,j,endj) UPWARDLX(UPWADD,wksp1,j,0); UPWARDENDX(wksp1,endj,0); UPWARDLX(UPWADD,wksp1,j,1);  UPWARDENDX(wksp1,endj,1)
#endif
#define UPREALIMAGD(UPWARDLX,UPWADD,wksp1,j) UPREALIMAG(UPWARDLX,UPWADD,wksp1,j); UPREALIMAG(UPWARDLX,UPWADD,wksp2,j)
#define UPREALIMAGENDD(UPWARDLX,UPWADD,UPWARDENDX,wksp1,j,endj) UPREALIMAGEND(UPWARDLX,UPWADD,UPWARDENDX,wksp1,j,endj); UPREALIMAGEND(UPWARDLX,UPWADD,UPWARDENDX,wksp2,j,endj)
#endif /*SSE2INTRINSIC*/
#define INNERUNROLLBASE(UPRI,UPRIEND,UPADD,UPADDEND,wksp1,startj) int j = startj;\
    for (; j < pshift-7; j+=8) {\
      UPRI(UPWARDL8,UPADD,wksp1,j);\
    }\
    if (j < pshift-3) {\
      UPRI(UPWARDL4,UPADD,wksp1,j);\
      j+=4;\
    }\
    if (j < pshift-2) {\
      UPRIEND(UPWARDL3,UPADD,UPADDEND,wksp1,j,j+2);\
    }\
    else if (j < pshift-1) {\
      UPRIEND(UPWARDL2,UPADD,UPADDEND,wksp1,j,j+1);\
    }\
    else if (j < pshift) {\
      UPRIEND(UPWARDL1,UPADD,UPADDEND,wksp1,j,j);\
    }\
    else {\
      j--;\
      UPRIEND(UPWARDL0,UPADD,UPADDEND,wksp1,j,j);\
    }
#define INNERUNROLLUP(UPRI,UPRIEND,UPADD,UPADDEND,wksp1)  INNERUNROLLBASE(UPRI,UPRIEND,UPADD,UPADDEND,wksp1,pshift-k)
#define INNERUNROLL1BASE(UPRI,UPADD,wksp1,startj) int j = startj;\
    for (; j < pshift-7; j+=8) {\
      UPRI(UPWARDL8,UPADD,wksp1,j);\
    }\
    if (j < pshift-3) {\
      UPRI(UPWARDL4,UPADD,wksp1,j);\
      j+=4;\
    }\
    if (j < pshift-2) {\
      UPRI(UPWARDL3,UPADD,wksp1,j);\
    }\
    else if (j < pshift-1) {\
      UPRI(UPWARDL2,UPADD,wksp1,j);\
    }\
    else if (j < pshift) {\
      UPRI(UPWARDL1,UPADD,wksp1,j);\
    }
#define INNERUNROLL1UP(UPRI,wksp1) INNERUNROLL1BASE(UPRI,UPWARDADD,wksp1, pshift-k)
#ifndef UNROLLLEVEL
#define UNROLLLEVEL 6
#endif

/*------------------------------------------------------------------------*/
#ifdef SSE2INTRINSIC
#define DNWARDADD(wkspm1,j,base) wkspm1[j+1]=_mm_add_pd(wkspm1[j+1],wkspm1[j])
#elif defined(CROSSORDER)
#define DNWARDADD(wksp1,j,base) wksp1[2*(j+1)+0] += wksp1[2*(j)+0]; wksp1[2*(j+1)+1] += wksp1[2*(j)+1]
#else
#define DNWARDADD(wksp1,j,base) wksp1[2*(j+1)+base] += wksp1[2*(j)+base]
#endif
#define DNWARDADDU2(wksp1,j,base) DNWARDADD(wksp1,j,base); DNWARDADD(wksp1,j-1,base)
#define DNWARDADDU3(wksp1,j,base) DNWARDADDU2(wksp1,j,base); DNWARDADD(wksp1,j-2,base)
#define DNWARDADDU4(wksp1,j,base) DNWARDADDU3(wksp1,j,base); DNWARDADD(wksp1,j-3,base)
#define DNWARDADDU5(wksp1,j,base) DNWARDADDU4(wksp1,j,base); DNWARDADD(wksp1,j-4,base)
#define DNWARDADDU6(wksp1,j,base) DNWARDADDU5(wksp1,j,base); DNWARDADD(wksp1,j-5,base)
#define DNWARDEND2(wksp1,j,base) DNWARDADD(wksp1,j,base)
#define DNWARDEND3(wksp1,j,base) DNWARDADDU2(wksp1,j,base); DNWARDEND2(wksp1,j,base)
#define DNWARDEND4(wksp1,j,base) DNWARDADDU3(wksp1,j,base); DNWARDEND3(wksp1,j,base)
#define DNWARDEND5(wksp1,j,base) DNWARDADDU4(wksp1,j,base); DNWARDEND4(wksp1,j,base)
#define DNWARDEND6(wksp1,j,base) DNWARDADDU5(wksp1,j,base); DNWARDEND5(wksp1,j,base)
#define INNERUNROLLDN(UPRI,UPRIEND,UPADD,UPADDEND,wksp1)  INNERUNROLLBASE(UPRI,UPRIEND,UPADD,UPADDEND,wksp1,k-1)
#define INNERUNROLL1DN(UPRI,wksp1) INNERUNROLL1BASE(UPRI,DNWARDADD,wksp1, k-1)
#endif