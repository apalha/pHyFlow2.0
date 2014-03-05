#include <math.h>
#include <string.h>
#include "channelpot.h"
#include "channelpotpanel.h"
#include "fmm.h"

#ifdef CHANNELPOT
//Every even coefficient in sum(1/z^k)= 1/2^k*zeta(k), odd coefficients are zero
static double shiftcoeff[]={-8.2246703342411320e-001,
1.3529040421389227e-001,
-3.1791970687014039e-002,
7.8443543452964403e-003,
-1.9550675295465195e-003,
4.8840140944985744e-004,
-1.2207778907898667e-004,
3.0518044502545429e-005,
-7.6294236548863605e-006,
1.9073504523506809e-006,
-4.7683727190518507e-007,
1.1920929665663112e-007,
-2.9802322831796255e-008,
7.4505806246797296e-009,
-1.8626451509656896e-009,
4.6566128741615968e-010,
-1.1641532183371108e-010,
2.9103830457157220e-011,
-7.2759576142098957e-012,
1.8189894035475108e-012,
-4.5474735088656752e-013,
1.1368683772162249e-013,
-2.8421709430404411e-014,
7.1054273576010271e-015,
-1.7763568394002520e-015,
4.4408920985006271e-016,
-1.1102230246251565e-016,
2.7755575615628914e-017,
-6.9388939039072284e-018,
1.7347234759768071e-018,
-4.3368086899420177e-019,
1.0842021724855044e-019,
-2.7105054312137611e-020,
6.7762635780344027e-021,
-1.6940658945086007e-021,
4.2351647362715017e-022,
-1.0587911840678754e-022,
2.6469779601696886e-023,
-6.6174449004242214e-024,
1.6543612251060553e-024,
-4.1359030627651384e-025,
1.0339757656912846e-025,
-2.5849394142282115e-026,
6.4623485355705287e-027,
-1.6155871338926322e-027,
4.0389678347315804e-028,
-1.0097419586828951e-028,
2.5243548967072378e-029,
-6.3108872417680944e-030,
1.5777218104420236e-030,
-3.9443045261050590e-031,
9.8607613152626476e-032,
-2.4651903288156619e-032,
6.1629758220391547e-033,
-1.5407439555097887e-033,
3.8518598887744717e-034,
-9.6296497219361793e-035,
2.4074124304840448e-035,
-6.0185310762101120e-036,
1.5046327690525280e-036,
-3.7615819226313200e-037,
9.4039548065783001e-038,
-2.3509887016445750e-038,
5.8774717541114375e-039};

static double shiftccoeff[]={-4.6740110027233950e-001,
2.9356063208376797e-002,
-2.8941532818842187e-003,
3.1035805059223803e-004,
-3.4082726089650959e-005,
3.7716971662391508e-006,
-4.1848103842300024e-007,
4.6474314758313410e-008,
-5.1628751133195463e-009,
5.7361539491116410e-010,
-6.3733550280887223e-011,
7.0814458878410152e-012,
-7.8682493382889651e-013,
8.7424897184582487e-014,
-9.7138736468101529e-015,
1.0793191413720154e-015,
-1.1992434293277278e-016,
1.3324926748191041e-017,
-1.4805474066921663e-018,
1.6450526701930181e-019,
-1.8278362986507120e-020,
2.0309292200975126e-021,
-2.2565880220803698e-022,
2.5073200244336639e-023,
-2.7859111382195949e-024,
3.0954568202279813e-025,
-3.4393964669135742e-026,
3.8215516299014096e-027,
-4.2461684776672074e-028,
4.7179649751853760e-029,
-5.2421833057613653e-030,
5.8246481175125625e-031,
-6.4718312416805982e-032,
7.1909236018673213e-033,
-7.9899151131859073e-034,
8.8776834590954528e-035,
-9.8640927323282797e-036,
1.0960103035920309e-036,
-1.2177892262133677e-037,
1.3530991402370754e-038,
-1.5034434891523058e-039,
1.6704927657247843e-040,
-1.8561030730275383e-041,
2.0623367478083756e-042,
-2.2914852753426394e-043,
2.5460947503807109e-044,
-2.8289941670896784e-045,
3.1433268523218653e-046,
-3.4925853914687390e-047,
3.8806504349652656e-048,
-4.3118338166280731e-049,
4.7909264629200810e-050,
-5.3232516254667564e-051,
5.9147240282963965e-052,
-6.5719155869959963e-053,
7.3021284299955514e-054,
-8.1134760333283900e-055,
9.0149733703648779e-056,
-1.0016637078183197e-056,
1.1129596753536885e-057,
-1.2366218615040985e-058,
1.3740242905601094e-059,
-1.5266936561778993e-060,
1.6963262846421104e-061};

//local functions
double* setupbinomial(int p);
void shiftm2ld(dcmplx *coeff1, dcmplx *coeff2,double * binomial, double* zpow, int pmax);
void shiftm2ldc(dcmplx *coeff1, dcmplx *coeff2,double * binomial, double* zpow, int pmax);
void shift_m2m_channel_(dcmplx *dcoeff,dcmplx *scoeff, double d, int pmax,int Nb);
void initstreamcoeff(double* ucoeff, double *dcoeff,const double* zr,const double* zi,
                     const double *mr, const double *mi, double z0, double d0, int begin, 
                     int end, int pmaxs,double sigma);
void shiftudstreamcoeff(double* udest, double* ddest, double* usource, double* dsource,
                        double uconst,double dconst,int pmaxs);
void shiftstreamcoeff(double* dest, double* source,double sconst,int pmaxs);
void upwardstreamshift(double* ufcoeff, double* dfcoeff,const double *z0,const double* d0,
                       double sigma,int pmaxs,int nlev,const int* levellist);
void downwardstreamshift(double* ufcoeff, double* dfcoeff,double* ulcoeff, double* dlcoeff,
                         const double *z0,const double* d0,double sigma,int pmaxs,int nlev,const  int* levellist);
void evalstreamcoeff(double* qr,double* qi,const double *zr,const double* zi,
                     const double* ulcoeff,const double* dlcoeff,const double* z0,
                     double d,const int *ixptr,int pmaxs,int Nb,double sigma);

//basic setup function to initialize parameters
void setchannelparams(channelparam* cparams,double width,double tol)
{
//   cparams->H=H;
  cparams->sigma=M_PI/(2*cparams->H);
  double Swidaim=cparams->H/3;
  
  int K=(int)ceil(width/Swidaim); //number of boxes necessary
  if(K<=0) //avoid log2(0)
    K=1;
  cparams->Nlev=2*(int)ceil(0.5*log2(K)); //number of stream expansion levels
  cparams->Nb=1<<cparams->Nlev;                //number of boxes for stream expansions
  cparams->Nt=2*cparams->Nb-1;                 //total number of stream expansion boxes                
  cparams->Swid=(width*1.0000000001)/cparams->Nb; //safety to keep all elements within a box
  if(width<Swidaim)                  //In case of one box, make it atleast as big as Swidaim (avoids divide by zero below)
    cparams->Swid=Swidaim;
  cparams->pmax=(int)ceil(log(1/tol)*Swidaim/cparams->Swid); //number of coefficients for stream expansion
}

//sorts the system using bucket sort for the channel part of the code
void channelsort(void *This0,channelparam* cparams,const double* zr,
                 const double* er,int* ix,int *jx,int* ixptr,int* jxptr,
                 double xmin, int N, int NE)
{
  MPexp *This=(MPexp*)This0;
  int i,j,itmp,itmp2,Nboxlevel,*buckets;
#ifdef RADIALSHRINK
  double dabs;
#endif
  
  cparams->z0=(double*)mxCalloc(cparams->Nt,sizeof(double));
  cparams->d0=(double*)mxCalloc(cparams->Nlev+1,sizeof(double));
  cparams->ixptr=(int*)mxMalloc((cparams->Nb+1)*sizeof(int));
  
  cparams->levellist=(int*)mxMalloc((cparams->Nlev+1)*sizeof(int));
  cparams->levellist[0]=0;
  for(i=1,j=1;i<=cparams->Nlev;i++,j<<=1) {
    cparams->levellist[i]=cparams->levellist[i-1]+j;
  }
  //sort uniform i x-direction using bucketsort
  memset(ixptr,0,(cparams->Nb+1)*sizeof(int));
  for(i=0;i<N;i++) //counting phase
    ixptr[(int)((zr[i]-xmin)/cparams->Swid)+1]++;
  for(i=2;i<=cparams->Nb;i++)//cumulative sum
    ixptr[i]+=ixptr[i-1];
  buckets=(int*)mxCalloc(cparams->Nb,sizeof(int));
  for(i=0;i<N;i++) { //add to buckets
    itmp2=(int)((zr[i]-xmin)/cparams->Swid);
    itmp=ixptr[itmp2]+buckets[itmp2]++;//itmp=ixptr[itmp]+buckets[itmp]++; does not work in visual studio
    ix[itmp]=i;
  }
  memcpy(cparams->ixptr,ixptr,(cparams->Nb+1)*sizeof(int));
  
  if(NE!=0) { //same for evaluation points
    cparams->jxptr=(int*)mxMalloc((cparams->Nb+1)*sizeof(int));
    memset(jxptr,0,(cparams->Nb+1)*sizeof(int));
    for(i=0;i<NE;i++)
      jxptr[(int)((er[i]-xmin)/cparams->Swid)+1]++;
    for(i=2;i<=cparams->Nb;i++)
      jxptr[i]+=jxptr[i-1];
    memset(buckets,0,cparams->Nb*sizeof(int));
    for(i=0;i<NE;i++) {
      itmp2=(int)((er[i]-xmin)/cparams->Swid);
      itmp=jxptr[itmp2]+buckets[itmp2]++;
      jx[itmp]=i;
    }
    memcpy(cparams->jxptr,jxptr,(cparams->Nb+1)*sizeof(int));
  }
  else
    cparams->jxptr=NULL;
  mxFree(buckets);
  
  //setup z and d
  cparams->d0[cparams->Nlev]=cparams->Swid/2;
  for(i=0;i<cparams->Nb;i++) {
    cparams->z0[cparams->levellist[cparams->Nlev]+i]=xmin+(i+0.5)*cparams->Swid;
  }
  for(i=cparams->Nlev-1;i>=0;i--) {
    Nboxlevel=1<<i;
    cparams->d0[i]=2*cparams->d0[i+1];
    for(j=0;j<Nboxlevel;j++)
      cparams->z0[cparams->levellist[i]+j]=0.5*(cparams->z0[cparams->levellist[i+1]+2*j]+cparams->z0[cparams->levellist[i+1]+2*j+1]);
  }
  
  //convert z and d to traditional fmm format
  for(i=0;i<=cparams->Nlev/2;i++) {
    Nboxlevel=1<<(i*2);
#ifdef RADIALSHRINK
    dabs=sqrt(cparams->d0[2*i]*cparams->d0[2*i]+0.25*cparams->H*cparams->H);
#endif
    for(j=0;j<Nboxlevel;j++) {
      COMPLEXASSIGN(This->root[This->lptr[i]+j].z0,cparams->z0[cparams->levellist[2*i]+j],0.5*cparams->H);
      COMPLEXASSIGN(This->root[This->lptr[i]+j].d0,cparams->d0[2*i],0.5*cparams->H);
#ifdef RADIALSHRINK
      This->root[This->lptr[i]+j].absd0=dabs; //not really necessary as no interactions should be performed with these anyway, but to keep compatibilty
#endif
    }
  }
}

//deallocates all variables in channelparam
void cleanupchannelparam(channelparam* cparams)
{
  mxFree(cparams->z0);
  mxFree(cparams->d0);
  mxFree(cparams->levellist);
  mxFree(cparams->ixptr);
  mxFree(cparams->jxptr);
}

//copies z and d to virtual boxes
void replicatezd(void *This0,channelparam* cparams)
{
  MPexp *This=(MPexp*)This0;
  const int Nf = 1 << (This->nlevel << 1);
  const int Nt = (Nf-1)/3+Nf;
  for(int i=0;i<Nt;i++) {
    COMPLEXASSIGN(This->root[Nt+i].z0,creal(This->root[i].z0),-cimag(This->root[i].z0));
    COMPLEXASSIGN(This->root[Nt+i].d0,creal(This->root[i].d0),cimag(This->root[i].d0));
    COMPLEXASSIGN(This->root[2*Nt+i].z0,creal(This->root[i].z0),cimag(This->root[Nt+i].z0)+2*cparams->H);
    COMPLEXASSIGN(This->root[2*Nt+i].d0,creal(This->root[i].d0),cimag(This->root[i].d0));
#ifdef RADIALSHRINK
    This->root[Nt+i].absd0=This->root[i].absd0;
    This->root[2*Nt+i].absd0=This->root[i].absd0;
#endif
    This->root[i+2*Nt].npanel=This->root[i+Nt].npanel=This->root[i].npanel;
  }
}

//copies source positions and strenghts to virtual points and sets up ixptr for it
void replicatepositions(void *This0,channelparam* cparams, void *tmpvar0,int N)
{
  MPexp *This=(MPexp*)This0;
  tmpvariables* tmpvar=(tmpvariables*)tmpvar0;
  const int Nf = 1 << (This->nlevel << 1);
  const int Nt = (Nf-1)/3+Nf;
  for(int i=0;i<N;i++) { //copy positions
    tmpvar->zrnew[i+N]=tmpvar->zrnew[i+N*2]=tmpvar->zrnew[i];
    tmpvar->mrnew[i+N]=tmpvar->mrnew[i+N*2]=tmpvar->mrnew[i];
    tmpvar->zinew[i+N]=-tmpvar->zinew[i]; //first mirror, conjugate
    tmpvar->zinew[i+2*N]=tmpvar->zinew[i+N]+2*cparams->H; //second mirror
  }
  if(tmpvar->minew!=NULL) {
    for(int i=0;i<N;i++) {
      tmpvar->minew[i+N]=tmpvar->minew[i+N*2]=-tmpvar->minew[i];
    }
  }
  for(int i=1;i<=Nf;i++) { //replicate ixptr to point to virtual points
    This->ixptr[i+Nf]=This->ixptr[i]+N;
    This->ixptr[i+2*Nf]=This->ixptr[i]+N*2;
  }
}

//copies coeff1 to virtual boxes. Note that the coefficient is conjugated there
void replicatecoeffs(void *This0,channelparam* cparams)
{
  MPexp *This=(MPexp*)This0;
  const int Nf = 1 << (This->nlevel << 1);
  const int Nt = (Nf-1)/3+Nf;
  const int P = This->pcoeff+1;
  for(int i=0;i<Nt;i++) {
    for(int j=0;j<P;j++) {
      COMPLEXASSIGN(This->root[0].coeff1[(i+Nt)*P+j],creal(This->root[0].coeff1[i*P+j]),-cimag(This->root[0].coeff1[i*P+j]));
      COMPLEXASSIGN(This->root[0].coeff1[(i+2*Nt)*P+j],creal(This->root[0].coeff1[i*P+j]),-cimag(This->root[0].coeff1[i*P+j]));
    }
  }
}
#ifndef SWAPI
#define SWAPI(X,Y) { int tmpi = X; X = Y; Y = tmpi; }
#endif

//copy of MPexp_connect_CPU_ but handles the initial channel sort as well. Only call if channelsort has been called in advance
void MPexp_connect_CPU_channelpot_(void *This0,void *C0,channelparam* cparams,
                                   int *levm2p,int *maxm2p,double cutoff)
/* Once the mesh has been obtained, this function determines the type
   of connection between all multipole boxes. Executed on the CPU. */
{
  // loop over all levels
  MPexp *This=(MPexp*)This0;
  mpSparse *C=(mpSparse*)C0;
  const int Nf = 1 << (This->nlevel << 1);
  const int Nt = (Nf-1)/3+Nf;
  int l=1,startlevel=cparams->Nlev/2;
  for(;l<startlevel;l++) { //lower levels, no connections at all as this will be handled by the stream expansions instead
    const int N0 = This->lptr[l]-This->lptr[l-1];
    const int N1 = This->lptr[l+1]-This->lptr[l];
    C[l].jcptr = (int *)mxCalloc(((2+(l == This->nlevel))*N1+1),sizeof(int));
    C[l].jcptr[0] = 0;
    C[l].kcptr = &C[l].jcptr[N1+1];
    C[l].ir = &C[l].kcptr[(1+(l == This->nlevel))*N1];
  }
  int nc=0;
  if(startlevel!=0){ //last level for stream expansions
    const int N0 = This->lptr[l]-This->lptr[l-1];
    const int N1 = This->lptr[l+1]-This->lptr[l];
    C[l].jcptr = (int *)mxCalloc(((9+(l == This->nlevel))*N1+1),sizeof(int));
    C[l].jcptr[0] = 0;
    C[l].kcptr = &C[l].jcptr[N1+1];
    C[l].ir = &C[l].kcptr[(1+(l == This->nlevel))*N1];
    int i;
    C[l].kcptr[0]=C[l].jcptr[1]=5;
    if(l==This->nlevel)
      C[l].kcptr[N1]=C[l].kcptr[0];
    C[l].ir[0]=1;
    C[l].ir[1]=Nt; //mirror lists
    C[l].ir[2]=1+Nt;
    C[l].ir[3]=2*Nt;
    C[l].ir[4]=1+2*Nt;
    for(i=1;i<N1-1;i++) {
      C[l].kcptr[i]=C[l].jcptr[i+1]=5+7*i;
      if(l==This->nlevel)
        C[l].kcptr[i+N1]=C[l].kcptr[i];
      C[l].ir[C[l].jcptr[i]]=i+1;
      C[l].ir[C[l].jcptr[i]+1]=i-1+Nt; //mirror lists
      C[l].ir[C[l].jcptr[i]+2]=i+Nt;
      C[l].ir[C[l].jcptr[i]+3]=i+1+Nt;
      C[l].ir[C[l].jcptr[i]+4]=i-1+2*Nt;
      C[l].ir[C[l].jcptr[i]+5]=i+2*Nt;
      C[l].ir[C[l].jcptr[i]+6]=i+1+2*Nt;
      nc++;
    }
    C[l].kcptr[i]=C[l].jcptr[i+1]=2+7*i;
    if(l==This->nlevel)
      C[l].kcptr[i+N1]=C[l].kcptr[i];
    C[l].ir[C[l].jcptr[i]]=i-1+Nt; //mirror lists
    C[l].ir[C[l].jcptr[i]+1]=i+Nt;
    C[l].ir[C[l].jcptr[i]+2]=i-1+2*Nt;
    C[l].ir[C[l].jcptr[i]+3]=i+2*Nt;
    nc=N1*7-5;
  }
  else { //special case if channelsort only gives one level
    mxFree(C[0].jcptr);
    C[0].jcptr=(int*)mxMalloc(7*sizeof(int));
    C[0].jcptr[0] = 0;
    C[0].jcptr[1] = 2;
    C[0].kcptr = &C[0].jcptr[2];
    C[0].kcptr[0]=C[0].jcptr[1];
    C[0].kcptr[1]=C[0].jcptr[1];
    C[0].ir = &C[0].kcptr[2];
    C[0].ir[0] = Nt;
    C[0].ir[1] = Nt*2;
    l=0;
    nc=2;
  }
  for (l++; l <= This->nlevel; l++) { //continue as normal
    const int N0 = This->lptr[l]-This->lptr[l-1];
    const int N1 = This->lptr[l+1]-This->lptr[l];
    mpexp *child = &This->root[This->lptr[l]];

    /* allocate connectivity matrix using that all near connections to
       higher index boxes at the previous level implies 4 connections
       for each of 4 children at this level (and 6 additional
       connections to higher indices {0->[1 2 3],1->[2 3],2->[3]}
       between the children of each parent) */
    C[l].jcptr = (int *)mxMalloc(((2+(l == This->nlevel))*N1+1+
                                  (6*N0+4*4*nc))*sizeof(int));
    C[l].jcptr[0] = 0;
    C[l].kcptr = &C[l].jcptr[N1+1];
    C[l].ir = &C[l].kcptr[(1+(l == This->nlevel))*N1];
    

    // setup connections
    nc = 0;  // counter for near connections
    for (int i = 0, j = 0; i < N0; i++)  // parent i
      for (int c = 0; c < 4; c++, j++) { // children j
        // forward- and backward counters in each column
        int fwd = C[l].jcptr[j];
        int bwd = fwd+3-c+4*(C[l-1].kcptr[i]-C[l-1].jcptr[i]);
        C[l].jcptr[j+1] = bwd;
            
        /* near connections to siblings > j: the connection could in
           principle be of the far type (except to the neighbour
           sibling), but this should occur only very rarely */
       
        for (int d = 1; d < 4-c; d++)
#if !defined(CUDASUPPORT) && !defined(M2PSINBOX)
          C[l].ir[fwd++] = j+d;
#else
        if (!mpexp_theta(&child[j], &child[j+d], cutoff))
          // near connection
          C[l].ir[fwd++] = j+d;
        else
          // far connection
          C[l].ir[--bwd] = j+d;
#endif /*CUDASUPPORT*/
            
        // parents (near) connections
        for (int ii = C[l-1].jcptr[i]; ii < C[l-1].kcptr[i]; ii++) {
          /* the children of the near connections of box i become
             either near connections or far connections of box j */
          int k;
          if(C[l-1].ir[ii]<Nt) //handle k for mirrors
            k = C[l-1].ir[ii] << 2; // 0-child of connection
          else if(C[l-1].ir[ii]<2*Nt)
            k = ((C[l-1].ir[ii]-Nt) << 2)+Nt; // 0-child of connection
          else
            k = ((C[l-1].ir[ii]-2*Nt) << 2)+2*Nt; // 0-child of connection
          // (note: minor differences to the C99-code in that
          // interactions through parents are not supported here)

          for (int d = 0; d < 4; d++, k++) {
//             if(j==63)
            // connection according to theta criterion
            if (!mpexp_theta(&child[j],&child[k],cutoff))  {      
              // near connection
              C[l].ir[fwd++] = k;
            }
            else {
              // far connection
              C[l].ir[--bwd] = k;
            }
          }
        }
        C[l].kcptr[j] = bwd;
        /* this splits the column into near connections
           [jcptr[j],kcptr[j]) and far connections
           [kcptr[j],jcptr[j+1]) */
        nc += bwd-C[l].jcptr[j];
        if (C[l].jcptr[j+1]-bwd > *maxm2p)
          *maxm2p = C[l].jcptr[j+1]-bwd;
        if (l < *levm2p && *maxm2p > 0)
          *levm2p = l;
            
        // last level: sort the near connections one step further
        if (l == This->nlevel) {
          int dir1,dir2; // 1: cluster to particle, -1: particle to cluster
          fwd = C[l].jcptr[j]-1;
          bwd = C[l].kcptr[j];
#if defined(CUDASUPPORT) || defined(M2PSINBOX)
          if (j == 0)
            bwd--; // skip box 0
          else {
            if (C[l].ir[fwd+1] == 0)
              fwd++;
#endif
            for ( ; ; ) {
              do
                fwd++;
              while (fwd < bwd &&
                      !(dir1 = mpexp_theta2(&child[j], &child[C[l].ir[fwd]],
                      cutoff)));
              
              do {
                bwd--;
                if (fwd <= bwd &&
                        (dir2 = mpexp_theta2(&child[j], &child[C[l].ir[bwd]],
                        cutoff)))
                  C[l].ir[bwd] *= dir2;
                else
                  break;
              } while (1);
              
              if (fwd >= bwd) break;
              SWAPI(C[l].ir[fwd], C[l].ir[bwd])
              C[l].ir[bwd] *= dir1;
            }
#if defined(CUDASUPPORT) || defined(M2PSINBOX)
          }
#endif
          C[l].kcptr[j+N1] = bwd+1;
          /* this splits the near connections into 'truly near'
             [jcptr[j],kcptr[j+N]) and 'less near'
             [kcptr[j+N],kcptr[j]) connections */
        }
      }
  }
}

//applies the correction due to mirror boxes far away from the channel according
//to chapter 4.2 in "Potential flow in channels" by Greengard (SIAM J. SCI. STAT. COMPUT 1990)
void addmirrorcorrection(void *This0,channelparam* cparams)
{
  MPexp *This=(MPexp*)This0;
  const int P=This->pcoeff+1;
  dcmplx *coeff1s=(dcmplx*)mxCalloc(cparams->Nb*P,sizeof(dcmplx));
  int level=cparams->Nlev/2;
  shift_m2m_channel_(coeff1s,This->root[This->lptr[level]].coeff1,
                     cparams->Swid,This->pcoeff,cparams->Nb);
  
  double *zpow1=(double*)mxMalloc(P*sizeof(double));
  double *zpow2=(double*)mxMalloc(P*sizeof(double));
  double Htmp=cparams->H*cparams->H;
  double *binomial=setupbinomial(This->pcoeff);
  for(int i=0;i<P;i++) {
    zpow1[i]=shiftcoeff[i]/Htmp;
    zpow2[i]=shiftccoeff[i]/Htmp;
    Htmp*=cparams->H*cparams->H;
  }
  for(int m=0;m<cparams->Nb;m++) {
    shiftm2ld(coeff1s+m*P, This->root[This->lptr[level]+m].coeff2,binomial, zpow1, This->pcoeff);
    shiftm2ldc(coeff1s+m*P, This->root[This->lptr[level]+m].coeff2,binomial, zpow2, This->pcoeff);
  }
  mxFree(coeff1s);
  mxFree(zpow1);
  mxFree(zpow2);
  mxFree(binomial);
}

//interactions through stream expansions, see "Potential flow in channels" by Greengard (SIAM J. SCI. STAT. COMPUT 1990)
//note that complex masses are used for potential points instead of pure imaginary as in the article
void streaminteraction(void *This0,channelparam* cparams,double *pr,double* pi,double* qr,
                       double* qi,const double* zr,const double* zi,
                       const double* er,const double* ei,const double* mr,
                       const double* mi,int N,int NE,const void *panels0)
{
  MPexp *This=(MPexp*)This0;
  const panel* panels=(const panel*)panels0;
  const int pmaxs=cparams->pmax;
  const int Nt=cparams->Nt;
  const int Nb=cparams->Nb;
  
  const int nlev=cparams->Nlev;
  const double H=cparams->H;
  const double sigma=cparams->sigma;
  const double *z0=cparams->z0;
  const double *d0=cparams->d0;
  const int *levellist=cparams->levellist;
  const int *ixptr=cparams->ixptr;
  const int *jxptr=cparams->jxptr;

  double *ufcoeff=(double*)mxCalloc(Nt*(pmaxs+1),sizeof(double));
  double *dfcoeff=(double*)mxCalloc(Nt*(pmaxs+1),sizeof(double));
  double *ulcoeff=(double*)mxCalloc(Nt*(pmaxs+1),sizeof(double));
  double *dlcoeff=(double*)mxCalloc(Nt*(pmaxs+1),sizeof(double));
  if(cparams->Nlev) { //no interactions through stream expansions if startlevel=0
    for(int i=0;i<Nb;i++) { //initialize stream expansions
      initstreamcoeff(ufcoeff+(levellist[nlev]+i)*(pmaxs+1), dfcoeff+(levellist[nlev]+i)*(pmaxs+1), 
                      zr, zi, mr, mi, z0[levellist[nlev]+i], d0[nlev], ixptr[i], ixptr[i+1], pmaxs, sigma);
      
      initstreamcoeffpanel(ufcoeff+(levellist[nlev]+i)*(pmaxs+1), dfcoeff+(levellist[nlev]+i)*(pmaxs+1),
                           panels,&This->root[This->lptr[cparams->Nlev/2]+i], 
                           z0[levellist[nlev]+i], d0[nlev],pmaxs,sigma,cparams->H);
    }
    //shift stream expansions
    upwardstreamshift(ufcoeff, dfcoeff,z0,d0,sigma,pmaxs,nlev,levellist);
    
    downwardstreamshift(ufcoeff, dfcoeff,ulcoeff,dlcoeff,z0,d0,sigma,pmaxs,nlev,levellist);
    if(pr!=NULL) { //evaluate stream expansions
      evalstreamcoeff(pr,pi,zr,zi,ulcoeff+(pmaxs+1)*levellist[nlev],dlcoeff+(pmaxs+1)*levellist[nlev],z0+levellist[nlev],d0[nlev],ixptr,pmaxs,Nb,sigma);
    }
    if(NE) {
      evalstreamcoeff(qr,qi,er,ei,ulcoeff+(pmaxs+1)*levellist[nlev],dlcoeff+(pmaxs+1)*levellist[nlev],z0+levellist[nlev],d0[nlev],jxptr,pmaxs,Nb,sigma);
    }
  }
}

//sets up binomial matrix for M2L shift used for the mirror corrections
double* setupbinomial(int p)
{
  double *binomial2 = (double*)mxCalloc(p*(p+1),sizeof(double));
  for (int l = 0; l <= p; l++)
    binomial2[p*l] = -1.0;
  for (int k = 1; k < p; k++) {
    binomial2[k] = -binomial2[k-1];
    for (int l = 1; l <= p; l++)
      binomial2[p*l+k] = binomial2[p*(l-1)+k]-binomial2[p*l+k-1];
  }
  return binomial2;
}

//m2l (m2ps) shift for mirror coefficients
void shiftm2ld(dcmplx *coeff1, dcmplx *coeff2,double * binomial, double* zpow, int pmax)
{
  double tmp;
  for(int m=0;m<=pmax;m++) {
    for(int k=1+((m+1)&1);k<=pmax;k+=2) {
      tmp=binomial[m*pmax+k-1]*zpow[(m+k)/2-1];
      COMPLEXADD(coeff2[m],creal(coeff1[k])*tmp,cimag(coeff1[k])*tmp);
    }
  }
}

//same as above, but with complex conjugate of coeff1
void shiftm2ldc(dcmplx *coeff1, dcmplx *coeff2,double * binomial, double* zpow, int pmax)
{
  double tmp;
  for(int m=0;m<=pmax;m++) {
    for(int k=1+((m+1)&1);k<=pmax;k+=2) {
      tmp=binomial[m*pmax+k-1]*zpow[(m+k)/2-1];
      COMPLEXADD(coeff2[m],creal(coeff1[k])*tmp,-cimag(coeff1[k])*tmp);
    }
  }
}

//special m2m shift that only shifts from the two neighbors to each side and with pure real distance. Not optimized as it will not be called that many times anyway
void shift_m2m_channel_(dcmplx *dcoeff,dcmplx *scoeff, double d, int pmax,int Nb)
/* Multipole to multipole shift from 4 children to 1 parent; This +=
   shift(that[0...p]) {child 0} + shift(that[p+1...2*p]) {child 1} +
   ...  and so on. */
{
  dcmplx *wksp=(dcmplx*)mxMalloc((pmax+1)*sizeof(dcmplx));
  for(int m=0;m<Nb;m++) {
    if(m!=0) { //unless it is the first box, shift the box to the left to this box
      memcpy(wksp, scoeff+(m-1)*(pmax+1), (pmax+1)*sizeof(dcmplx));
      
      for(int k = pmax; k >= 2; k--){
        for(int j = k; j <= pmax; j++) {
          COMPLEXSUB(wksp[j],creal(wksp[j-1])*d,cimag(wksp[j-1])*d);
        }
      }
      for(int k = 0; k <= pmax; k++) {
        COMPLEXADD(dcoeff[m*(pmax+1)+k],creal(wksp[k]),cimag(wksp[k]));
      }
    }
    for(int k = 0; k <= pmax; k++) {
      COMPLEXADD(dcoeff[m*(pmax+1)+k],creal(scoeff[m*(pmax+1)+k]),cimag(scoeff[m*(pmax+1)+k]));
    }
    if(m!=Nb-1) {//unless it is the last box, shift the box to the right to this box
      memcpy(wksp, scoeff+(m+1)*(pmax+1), (pmax+1)*sizeof(dcmplx));
      for(int k = pmax; k >= 2; k--){
        for(int j = k; j <= pmax; j++) {
          COMPLEXADD(wksp[j],creal(wksp[j-1])*d,cimag(wksp[j-1])*d);
        }
      }
      
      for(int k = 0; k <= pmax; k++) {
        COMPLEXADD(dcoeff[m*(pmax+1)+k],creal(wksp[k]),cimag(wksp[k]));
      }
    }
  }
  mxFree(wksp);
}

//initializes the stream coefficients
//ucoeff(k)=2*sigma*(sum(m[i]*exp(2*sigma*z[i]*k)+conj(m[i])*exp(2*sigma*z[i]*k)))
//dcoeff(k)=2*sigma*(sum(m[i]*exp(-2*sigma*z[i]*k)+conj(m[i])*exp(-2*sigma*z[i]*k)))
void initstreamcoeff(double* ucoeff, double *dcoeff,const double* zr,const double* zi,
                     const double *mr, const double *mi, double z0, double d0, int begin, 
                     int end, int pmaxs,double sigma)
{
  double xdist,ydist,mrl,mil,tu,tub,td,tdb,t2,t2b,tmp,stmp,ctmp,stmpb,ctmpb;
  for(int i=begin;i<end;i++) {
    xdist=zr[i]-z0;
    ydist=zi[i];
    mrl=mr[i];
    if(mi!=NULL)
      mil=mi[i];
    else
      mil=0;
    tu=tub=exp(2*sigma*(xdist-d0));//to avoid overflow, use edge of box as expansion center
    td=tdb=exp(-2*sigma*(xdist+d0));//to avoid overflow, use edge of box as expansion center
    t2b=t2=2*sigma*ydist;
    stmp=stmpb=sin(t2);
    ctmp=ctmpb=cos(t2);
    ucoeff[0]+=mrl;
    dcoeff[0]-=mrl;
    ucoeff[1]+=tu*(-mil*stmp+mrl*ctmp);
    dcoeff[1]+=td*(mil*stmp+mrl*ctmp);
    for(int j=2;j<=pmaxs;j++) {
      tu*=tub;
      td*=tdb;
      tmp=ctmp*stmpb+ctmpb*stmp;
      ctmp=ctmp*ctmpb-stmp*stmpb;
      stmp=tmp;
      ucoeff[j]+=tu*(-mil*stmp+mrl*ctmp);
      dcoeff[j]+=td*(mil*stmp+mrl*ctmp);
    }
  }
  ucoeff[0]*=2*sigma;
  dcoeff[0]*=2*sigma;
  for(int j=1;j<=pmaxs;j++) {
    ucoeff[j]*=4*sigma;
    dcoeff[j]*=4*sigma;
  }
}

//shift coefficients to new center. This is given by a[k]=a[k]*exp(2*sigma*z*k) where z is the shifting distance
void shiftudstreamcoeff(double* udest, double* ddest, double* usource, double* dsource,
                        double uconst,double dconst,int pmaxs)
{
  double utmp=uconst,dtmp=dconst;
  udest[0]+=usource[0];
  ddest[0]+=dsource[0];
  udest[1]+=utmp*usource[1];
  ddest[1]+=dtmp*dsource[1];
  for(int i=2;i<=pmaxs;i++) {
    utmp*=uconst;
    dtmp*=dconst;
    udest[i]+=utmp*usource[i];
    ddest[i]+=dtmp*dsource[i];
  }
}

void shiftstreamcoeff(double* dest, double* source,double sconst,int pmaxs)
{
  double tmp=sconst;
  dest[0]+=source[0];
  dest[1]+=tmp*source[1];
  for(int i=2;i<=pmaxs;i++) {
    tmp*=sconst;
    dest[i]+=tmp*source[i];
  }
}

void upwardstreamshift(double* ufcoeff, double* dfcoeff,const double *z0,const double* d0,
                       double sigma,int pmaxs,int nlev,const int* levellist)
{
  int Nboxlevel;
  double uconst,dconst,d;
  for(int i=nlev-1;i>=0;i--) {
    Nboxlevel=1<<i;
    d=2*d0[i+1];
    uconst=exp(-2*sigma*d); //note that the distance is constant for all shifts as it is a uniform code
    for(int j=0;j<Nboxlevel;j++) {
      shiftudstreamcoeff(ufcoeff+(pmaxs+1)*(levellist[i]+j), dfcoeff+(pmaxs+1)*(levellist[i]+j),ufcoeff+(pmaxs+1)*(levellist[i+1]+2*j), dfcoeff+(pmaxs+1)*(levellist[i+1]+2*j),uconst,1,pmaxs);
      shiftudstreamcoeff(ufcoeff+(pmaxs+1)*(levellist[i]+j), dfcoeff+(pmaxs+1)*(levellist[i]+j),ufcoeff+(pmaxs+1)*(levellist[i+1]+2*j+1), dfcoeff+(pmaxs+1)*(levellist[i+1]+2*j+1),1,uconst,pmaxs);
    }
  }
}

void downwardstreamshift(double* ufcoeff, double* dfcoeff,double* ulcoeff, double* dlcoeff,
                         const double *z0,const double* d0,double sigma,int pmaxs,int nlev,const  int* levellist)
{
  int Nboxlevel;
  double uconst,uconst2;
  for(int i=2;i<=nlev;i++) {
    Nboxlevel=1<<i;
    uconst=exp(-4*sigma*d0[i]); //note that the distance is constant for all shifts as it is a uniform code
    uconst2=uconst*uconst; //note that the distance is constant for all shifts as it is a uniform code    for(int j=0;j<Nboxlevel;j++) {
    for(int j=0;j<Nboxlevel;j++) {
      //depending on if a box is the first or second box, it is either two boxes upstream and one downstream or the opposite that are to be shifted
      if(j>=2) {
        shiftstreamcoeff(ulcoeff+(pmaxs+1)*(levellist[i]+j),ufcoeff+(pmaxs+1)*(levellist[i]+j-2),uconst,pmaxs);
        if(j%2) {
          shiftstreamcoeff(ulcoeff+(pmaxs+1)*(levellist[i]+j),ufcoeff+(pmaxs+1)*(levellist[i]+j-3),uconst2,pmaxs);
        }
      }
      if(j<Nboxlevel-2) {
        shiftstreamcoeff(dlcoeff+(pmaxs+1)*(levellist[i]+j),dfcoeff+(pmaxs+1)*(levellist[i]+j+2),uconst,pmaxs);
        if(j%2==0) {
          shiftstreamcoeff(dlcoeff+(pmaxs+1)*(levellist[i]+j),dfcoeff+(pmaxs+1)*(levellist[i]+j+3),uconst2,pmaxs);
        }
      }
      if(j%2)
        shiftudstreamcoeff(ulcoeff+(pmaxs+1)*(levellist[i]+j),dlcoeff+(pmaxs+1)*(levellist[i]+j),ulcoeff+(pmaxs+1)*(levellist[i-1]+j/2),dlcoeff+(pmaxs+1)*(levellist[i-1]+j/2),uconst,1,pmaxs);
      else
        shiftudstreamcoeff(ulcoeff+(pmaxs+1)*(levellist[i]+j),dlcoeff+(pmaxs+1)*(levellist[i]+j),ulcoeff+(pmaxs+1)*(levellist[i-1]+j/2),dlcoeff+(pmaxs+1)*(levellist[i-1]+j/2),1,uconst,pmaxs);
    }
  }
}

//evaluates the velocity from the stream coefficients
void evalstreamcoeff(double* qr,double* qi,const double *zr,const double* zi,
                     const double* ulcoeff,const double* dlcoeff,const double* z0,
                     double d,const int *ixptr,int pmaxs,int Nb,double sigma)
{
  double xdist,ydist,z,xvel,yvel,stmp,ctmp,stmpb,ctmpb,tmp,t2,tu,tub,td,tdb;
  for(int i=0;i<Nb;i++) {
    z=z0[i];
    for(int j=ixptr[i];j<ixptr[i+1];j++) {
      
      xdist=zr[j]-z;
      ydist=zi[j];
      tu=tub=exp(-2*sigma*(xdist+d));
      td=tdb=exp(2*sigma*(xdist-d));
      t2=2*sigma*ydist;
      stmp=stmpb=sin(t2);
      ctmp=ctmpb=cos(t2);
      xvel=-ulcoeff[(pmaxs+1)*i]-dlcoeff[(pmaxs+1)*i]+(-ulcoeff[(pmaxs+1)*i+1]*tu+dlcoeff[(pmaxs+1)*i+1]*td)*ctmp;
      yvel=(ulcoeff[(pmaxs+1)*i+1]*tu+dlcoeff[(pmaxs+1)*i+1]*td)*stmp;
      for(int k=2;k<=pmaxs;k++) {
        tu*=tub;
        td*=tdb;
        tmp=ctmp*stmpb+ctmpb*stmp;
        ctmp=ctmp*ctmpb-stmp*stmpb;
        stmp=tmp;
        xvel+=(-ulcoeff[(pmaxs+1)*i+k]*tu+dlcoeff[(pmaxs+1)*i+k]*td)*ctmp;
        yvel+=(ulcoeff[(pmaxs+1)*i+k]*tu+dlcoeff[(pmaxs+1)*i+k]*td)*stmp;
      }
      qr[j]+=xvel;
      qi[j]+=yvel;
    }
  }
}
#endif /*CHANNELPOT*/