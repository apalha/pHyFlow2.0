#include <string.h>
#include <math.h>
#include <float.h>
#include "channelpotpanel.h"
#include "channelpot.h"
#ifdef CHANNELPOT
void channelsortpanel(MPexp* This,const panel *panels,int* tmpptr1,int* tmpptr2,
                      int** panelptrlist,int endlevel)
{
  int panelsum=This->root[0].npanel,panelptrindex;
  int outcount1,outcount2;
  dcmplx z,d,z0,d0,z1;
  for (int l = 1; l < endlevel; l++) {
    mpexp *parent = &This->root[This->lptr[l-1]];
    mpexp *child = &This->root[This->lptr[l]];
    const int N0 = This->lptr[l]-This->lptr[l-1];
    // maximum number of panels this round
    panelptrlist[l] = (int *)mxMalloc(panelsum*4*3*sizeof(int));
    // the number of panels in the lower level to save memory in the
    // next allocation
    panelsum = 0;
    panelptrindex = 0; // the index in the array

    // loop over all parent boxes
    for (int i = 0, j = 0; i < N0; i++, j += 4) {
      z = parent[i].z0;
      d = parent[i].d0;
      COMPLEXASSIGN(z0,creal(z)-creal(d)/2,cimag(z));
      COMPLEXASSIGN(z1,creal(z)+creal(d)/2,cimag(z));
      COMPLEXASSIGN(d0,creal(d)/2,cimag(d));
      outcount1 = outcount2 = 0;
      if (parent[i].npanel) {
        MPexp_split_panel(z,d,panels,parent[i].panelptr,
                          tmpptr1,tmpptr2,parent[i].npanel,
                          &outcount1,&outcount2,creal(z),1);
      }

      // set up the panel pointers
      child[j].panelptr = panelptrlist[l]+panelptrindex;
      child[j+1].panelptr = panelptrlist[l]+panelptrindex+outcount1;
      panelptrindex += outcount1*2;
      child[j+2].panelptr = panelptrlist[l]+panelptrindex;
      child[j+3].panelptr = panelptrlist[l]+panelptrindex+outcount2;
      panelptrindex += outcount2*2;

      // same procedure for the resulting two boxes...

      if (outcount1)
        MPexp_split_panel(z0,d0,panels,tmpptr1,
                          child[j].panelptr,child[j+1].panelptr,
                          outcount1,
                          &(child[j].npanel),&(child[j+1].npanel),creal(z0),1);
      else
        child[j].npanel = child[j+1].npanel = 0;
      // and again...

      if (outcount2)
        MPexp_split_panel(z1,d0,panels,tmpptr2,
                          child[j+2].panelptr,child[j+3].panelptr,
                          outcount2,
                          &(child[j+2].npanel),&(child[j+3].npanel),creal(z1),1);
      else
        child[j+2].npanel = child[j+3].npanel=0;
      panelsum += child[j].npanel+child[j+1].npanel+
              child[j+2].npanel+child[j+3].npanel;

    }
  }
}
void replicatepanelpositions(panel* panels, int Npanel,double channelheight)
{
  for(int i=0;i<Npanel;i++) {
    memcpy(panels+Npanel+i,panels+i,sizeof(panel));
    panels[i+Npanel].y1=-panels[i+Npanel].y1;
    panels[i+Npanel].y2=-panels[i+Npanel].y2;
    panels[i+Npanel].k1=-panels[i+Npanel].k1;
    panels[i+Npanel].k2=-panels[i+Npanel].k2;
    panels[i+Npanel].theta=-panels[i+Npanel].theta;
    panels[i+Npanel].sintheta=-panels[i+Npanel].sintheta;
    panels[i+Npanel].istrength1=-panels[i+Npanel].istrength1;
    panels[i+Npanel].istrength2=-panels[i+Npanel].istrength2;
    memcpy(panels+2*Npanel+i,panels+Npanel+i,sizeof(panel));
    panels[i+2*Npanel].y1+=2*channelheight;
    panels[i+2*Npanel].y2+=2*channelheight;
  }
}

void replicatepanelindices(MPexp *This,channelparam* cparams, int Npanel)
{
  const int Nf = 1 << (This->nlevel << 1);
  const int Nt = (Nf-1)/3+Nf;
  int *panelptr1=This->root[Nt-1].panelptr+This->root[Nt-1].npanel;
  for(int i=0;i<Nt;i++) { //replicate panelptr to point to virtual points
    This->root[i+Nt].panelptr=panelptr1;
    This->root[i+Nt].npanel= This->root[i].npanel;
    for(int j=0;j<This->root[i].npanel;j++) {
      This->root[i+Nt].panelptr[j]=This->root[i].panelptr[j]+Npanel;
    }
    panelptr1+=This->root[i].npanel;
  }
  for(int i=0;i<Nt;i++) { //replicate panelptr to point to virtual points
    This->root[i+2*Nt].panelptr=panelptr1;
    This->root[i+2*Nt].npanel= This->root[i].npanel;
    for(int j=0;j<This->root[i].npanel;j++) {
      This->root[i+2*Nt].panelptr[j]=This->root[i].panelptr[j]+2*Npanel;
    }
    panelptr1+=This->root[i].npanel;
  }
}

//expansion given by
//a(m)=2*real((((G0*2*sigma*et*m-k)*(exp(2*sigma*(z2-z0)*m)-exp(2*sigma*(z1-z0)*m))+exp(2*sigma*(z2-z0)*m)*k*(2*sigma*m*(z2-z1))))/(2*sigma*et2*m^2))
//b(m)=2*real((((G0*2*sigma*et*m+k)*(exp(2*sigma*(z2-z0)*m)-exp(2*sigma*(z1-z0)*m))+exp(2*sigma*(z2-z0)*m)*k*(2*sigma*m*(z2-z1))))/(2*sigma*et2*m^2))
void initstreamcoeffpanel(double* ucoeff, double *dcoeff, const panel* panels,
                          mpexp *box, double z0, double d0,
                          int pmaxs,double sigma,double channelheight)
{
  panel smallpanel;
  dcmplx z,d;
  COMPLEXASSIGN(z,z0,channelheight/2);
  COMPLEXASSIGN(d,d0,channelheight/2);
  for (int i = 0; i < box->npanel; i++) {
    if (panelinbox(&smallpanel,panels+box->panelptr[i],
                   box->z0,box->d0)) {
      dcmplx z1,z2,G1,k,expit,t1,t2,denom;
      COMPLEXASSIGN(z1,smallpanel.x1-z0,smallpanel.y1);
      COMPLEXASSIGN(z2,smallpanel.x2-z0,smallpanel.y2);
      COMPLEXASSIGN(G1,smallpanel.strength1,smallpanel.istrength1);
      COMPLEXASSIGN(k,(smallpanel.strength2-smallpanel.strength1)*(1/smallpanel.lambda),(smallpanel.istrength2-smallpanel.istrength1)*(1/smallpanel.lambda));
      COMPLEXASSIGN(expit,smallpanel.costheta,smallpanel.sintheta);
      denom=(dcmplx)1/(2*sigma*expit*expit);
      t1=2*sigma*(G1*expit+k*(z2-z1))*denom;
      t2=2*sigma*G1*expit*denom;
      k=k*denom;
      dcmplx exput1,exput2,expdt1,expdt2;
      dcmplx exput10,exput20,expdt10,expdt20;
      double tu,td,t2b,stmp,ctmp;
      tu=exp(2*sigma*(creal(z1)-d0));
      td=exp(-2*sigma*(creal(z1)+d0));
      t2b=2*sigma*cimag(z1);
      stmp=sin(t2b);
      ctmp=cos(t2b);
      COMPLEXASSIGN(exput1,tu*ctmp,tu*stmp);
      COMPLEXASSIGN(expdt1,td*ctmp,-td*stmp);
      COMPLEXASSIGN(exput10,tu*ctmp,tu*stmp);
      COMPLEXASSIGN(expdt10,td*ctmp,-td*stmp);

      tu=exp(2*sigma*(creal(z2)-d0));
      td=exp(-2*sigma*(creal(z2)+d0));
      t2b=2*sigma*cimag(z2);
      stmp=sin(t2b);
      ctmp=cos(t2b);
      COMPLEXASSIGN(exput2,tu*ctmp,tu*stmp);
      COMPLEXASSIGN(expdt2,td*ctmp,-td*stmp);
      COMPLEXASSIGN(exput20,tu*ctmp,tu*stmp);
      COMPLEXASSIGN(expdt20,td*ctmp,-td*stmp);

      ucoeff[0]+=sigma*(smallpanel.strength1+smallpanel.strength2)*smallpanel.lambda;
      dcoeff[0]-=sigma*(smallpanel.strength1+smallpanel.strength2)*smallpanel.lambda;
      ucoeff[1]+=2*creal((t1-k)*exput2-(t2-k)*exput1);
      dcoeff[1]-=2*creal((t1+k)*expdt2-(t2+k)*expdt1);
      for(int j=2;j<=pmaxs;j++) {
        exput1*=exput10;
        expdt1*=expdt10;
        exput2*=exput20;
        expdt2*=expdt20;
        double invj=1/(double)j;
        ucoeff[j]+=2*creal((t1*invj-k*invj*invj)*exput2-(t2*invj-k*invj*invj)*exput1);
        dcoeff[j]-=2*creal((t1*invj+k*invj*invj)*expdt2-(t2*invj+k*invj*invj)*expdt1);
      }
    }
  }
}





//rutines for direct interaction of panels

/*Matlab code for this function (except that matlab code does not check for problems with branch cut
function output=calculatechannelpanelvelocity(evalpoint,panelpositions,panelstrengths,H)
sigma=pi/(2*H);
output=zeros(size(evalpoint));
for k=1:size(panelpositions,2)
  z1mz0=panelpositions(2,k)-panelpositions(1,k);
  theta=atan2(imag(z1mz0),real(z1mz0));
  lambda=hypot(imag(z1mz0),real(z1mz0));
  zmz0=evalpoint-panelpositions(1,k);
  zmz1=evalpoint-panelpositions(2,k);
  lsinh0=log(2*sinh(sigma*(zmz0)));
  lsinh1=log(2*sinh(sigma*(zmz1)));
  zmz0c=evalpoint-conj(panelpositions(1,k));
  zmz1c=evalpoint-conj(panelpositions(2,k));
  lsinh0c=log(2*sinh(sigma*(zmz0c)));
  lsinh1c=log(2*sinh(sigma*(zmz1c)));
  k2=(panelstrengths(2,k)-panelstrengths(1,k))/lambda;
  G0=panelstrengths(1,k);
  G1=panelstrengths(2,k);
  rot=exp(-2i*theta);
  rot0=exp(-1i*theta);
  output=output +rot0.*(G0.*lsinh0-G1.*lsinh1)+0.5.*k2.*rot./sigma.*(sigma.^2*(zmz0.^2-zmz1.^2)+polylog2(exp(-2.*sigma.*zmz0))-polylog2(exp(-2.*sigma.*zmz1)))-...
                (conj(rot0).*(G0.*lsinh0c-G1.*lsinh1c)+0.5.*k2.*conj(rot)./sigma.*(sigma.^2.*(zmz0c.^2-zmz1c.^2)+polylog2(exp(-2.*sigma.*zmz0c))-polylog2(exp(-2.*sigma.*zmz1c))));
end*/
typedef struct {
  double re, im;
} dcx;

//constants for gaussian quadrature
static double x[] = {5.29953250417503070e-003,
2.7712488463383700e-002,
6.7184398806084122e-002,
1.2229779582249845e-001,
1.9106187779867806e-001,
2.7099161117138637e-001,
3.5919822461037054e-001,
4.5249374508118123e-001,
5.4750625491881877e-001,
6.4080177538962946e-001,
7.2900838882861363e-001,
8.0893812220132189e-001,
8.7770220417750155e-001,
9.3281560119391593e-001,
9.7228751153661630e-001,
9.9470046749582497e-001};

static double invx[] = {1.88695889535952260e+002,
3.6084814300285316e+001,
1.4884407954387195e+001,
8.1767622488584184e+000,
5.2339064784744771e+000,
3.6901511293187537e+000,
2.7839781254061595e+000,
2.2099752999250661e+000,
1.8264631518196528e+000,
1.5605449897387780e+000,
1.3717263276034746e+000,
1.2361884952074593e+000,
1.1393385994024070e+000,
1.0720232366612377e+000,
1.0285023597799652e+000,
1.0053277671795176e+000};

static double w[] = {1.35762297058772670e-002,
3.1126761969324106e-002,
4.7579255841246240e-002,
6.2314485627767001e-002,
7.4797994408288479e-002,
8.4578259697501365e-002,
9.1301707522461584e-002,
9.4725305227534154e-002,
9.4725305227534154e-002,
9.1301707522461584e-002,
8.4578259697501323e-002,
7.4797994408288590e-002,
6.2314485627766432e-002,
4.7579255841246240e-002,
3.1126761969324106e-002,
1.3576229705877267e-002};

typedef struct {
  double xdist;
  double ydist;
  dcx G0;
  dcx G1;
  dcx rot;
  dcx rot2;
  double ks;
} panelstruct;

//macros for complex arithmetics
#define DCXASSIGN(x, y, z) (x).re=(y);(x).im=(z);                               //x=(y,z)
#define DCXADD(x, y, z) ((x).re=(y).re+(z).re);((x).im=(y).im+(z).im);          //x=y+z
#define DCXSUB(x, y, z) ((x).re=(y).re-(z).re);((x).im=(y).im-(z).im);          //x=y-z
#define DCXADDC(x, y, z, c) ((x).re=(y).re+c*(z).re);((x).im=(y).im+c*(z).im);  //x=y+c*z, c real
#define DCXSUBC(x, y, z, c) ((x).re=(y).re-c*(z).re);((x).im=(y).im-c*(z).im);  //x=y-c*z, c real
#define DCXABS(x) hypot((x).re, (x).im)                                         //abs(x)
#define DCXNORM(x) (x).re*(x).re+(x).im*(x).im                                  //x*conj(x)
#define DCXMUL(x, y, z, tmp) (tmp)=(y).re*(z).re-(y).im*(z).im;(x).im=(y).re*(z).im+(y).im*(z).re;(x).re=tmp;     //x=y*z
#define DCXMULN(x, y, z, tmp) (tmp)=-(y).re*(z).re+(y).im*(z).im;(x).im=-(y).re*(z).im-(y).im*(z).re;(x).re=tmp;  //x=y*(-z)
#define DCXCMUL(x, y, z, tmp) (tmp)=(y).re*(z).re+(y).im*(z).im;(x).im=(y).re*(z).im-(y).im*(z).re;(x).re=tmp;    //x=conj(y)*z
#define DCXCMULN(x, y, z, tmp) (tmp)=-(y).re*(z).re-(y).im*(z).im;(x).im=-(y).re*(z).im+(y).im*(z).re;(x).re=tmp; //x=conj(y)*(-z)
#define DCXMULC(x, y, z, tmp) (tmp)=(y).re*(z).re+(y).im*(z).im;(x).im=-(y).re*(z).im+(y).im*(z).re;(x).re=tmp;   //x=y*conj(z)
#define DCXSQR(x, y, tmp) (tmp)=(y).re*(y).re-(y).im*(y).im;(x).im=2*(y).re*(y).im;(x).re=tmp;                    //x=y*y
#define DCXLOG(x, y) (x).re=0.5*log((y).re*(y).re+(y).im*(y).im);(x).im=atan2((y).im, (y).re);                    //x=log(y)
#define DCXLOGN(x, y) (x).re=0.5*log((y).re*(y).re+(y).im*(y).im);(x).im=atan2(-(y).im, -(y).re);                 //x=log(-y)
#define DCXRECIPROCAL(x, y, tmp) tmp=(y).re*(y).re+(y).im*(y).im;(x).re=(y).re/tmp;(x).im=-(y).im/tmp;            //x=1/y

dcx polylog2(double re, double im);
void fillpanelstruct(panelstruct *output, const panel* panels, double sigma);
void calculatepanelvelocity(dcx* res0, dcx* res1,dcx* res0c, dcx* res1c, double er, double ei, const panel *panels,panelstruct *panels2,double sigma);
/*------------------------------------------------------------------------*/
void directInteractChannelpotPanel(const panel *panels, double *qr, double *qi,const double* er,const double* ei, const int Npanel,const int N,const double channelheight)
{
  double sigma=M_PI/(2*channelheight);
  panelstruct panels2;
  dcx res0,res1,res0c,res1c;
  for(int i=0;i<Npanel;i++) {
    fillpanelstruct(&panels2, panels+i, sigma);
    for(int j=0;j<N;j++) {
      calculatepanelvelocity(&res0, &res1, &res0c, &res1c, er[j], ei[j], panels, &panels2, sigma);
      qr[j]-=panels2.G0.re*(res0.re+res0c.re)-panels2.G0.im*(res0.im-res0c.im)+panels2.G1.re*(res1.re+res1c.re)-panels2.G1.im*(res1.im-res1c.im);
      qi[j]-=panels2.G0.re*(res0.im+res0c.im)+panels2.G0.im*(res0.re-res0c.re)+panels2.G1.re*(res1.im+res1c.im)+panels2.G1.im*(res1.re-res1c.re);
    }
  }
}
/*------------------------------------------------------------------------*/
//evaluation of the dilogaritm according to "NUMERICAL EVALUATION OF THE DILOGARITHM OF COMPLEX ARGUMENT"
//by Osácar, C., Palacián, J., & Palacios, M in Celestial Mechanics & Dynamical Astronomy, Volume 62, Issue 1, pp.93-98
dcx polylog2(double re, double im) {
  int switchsign, i;
  double k, tmp, dist;
  dcx zk, err, z, y, ctmp;
  DCXASSIGN(y, 0, 0);
  DCXASSIGN(z, re, im);
  dist=hypot(re, im);
  if(dist>1) { //transform to evaluation of abs(z)<1 according to polylog2(z)=-polylog2(1/z)-0.5*log(-z)^2-pi^2/6
    if(im==0) { //special case for real argument, has problems with the branch cut for re>1 otherwise
      if(re<0) {
        tmp=log(-re);
        tmp*=tmp;
        y.re+=0.5*tmp+M_PI*M_PI/6;
        z.re=1/re;
        switchsign=1;
        dist=-z.re;
      }
      else {
        DCXASSIGN(ctmp, log(re), M_PI);//has to choose the other branch here for consistency with log
        DCXMUL(ctmp, ctmp, ctmp, tmp);
        DCXADDC(y, y, ctmp, 0.5);
        y.re+=M_PI*M_PI/6;
        z.re=1/re;
        switchsign=1; //solve -polylog2(z) instead, switch sign at end
        dist=z.re;
      }
    }
    else {
      DCXLOGN(ctmp, z);
      DCXMUL(ctmp, ctmp, ctmp, tmp);
      DCXADDC(y, y, ctmp, 0.5);
      y.re+=M_PI*M_PI/6;
      DCXRECIPROCAL(z, z, tmp); //switch z to 1/z for rest of algorithm
      switchsign=1; //solve -polylog2(z) instead, switch sign at end
      dist=DCXABS(z);
    }
  }
  else
    switchsign=0;
  DCXASSIGN(ctmp, 1-z.re, -z.im); //check if too close to 1, which has bad numerics for gaussian quadrature
  if(DCXNORM(ctmp)<0.25) { //transform according to polylog2(z)=-polylog2(1-z)-log(z)*log(1-z)+pi^2/6
    if(re==1&&im==0) { //special case to avoid NaN
      DCXASSIGN(y, M_PI*M_PI/6, 0);
      return y;
    }
    if(switchsign) { //switch the sign again to solve for -polylog2(z), if first switch has been performed, y needs to change sign here
      DCXASSIGN(y, -y.re, -y.im);
      switchsign=0;
    }
    else
      switchsign=1;
    DCXLOG(zk, ctmp);
    DCXLOG(ctmp, z);
    DCXMUL(ctmp, ctmp, zk, tmp);
    DCXADD(y, y, ctmp);
    y.re-=M_PI*M_PI/6;
    DCXASSIGN(z, 1-z.re, -z.im);
    dist=DCXABS(z);
  }
  if(dist>1.001) //should never happen
    #ifndef C_CODE
    mexErrMsgTxt("Values must be below 1")
    #endif
  if(dist>0.5) { //gaussian quadrature for abs(z)>0.5
    for(i=0;i<16;i++) {
      DCXASSIGN(ctmp, 1-z.re*x[i], -z.im*x[i]);
      DCXLOG(zk, ctmp);
      DCXSUBC(y, y, zk, w[i]*invx[i]);
    }
  }
  else { //evaluate with sum: polylog2(z)=sum{k=1 to infinity}z^k/k^2
    k=1;
    DCXASSIGN(err, 1, 1);
    DCXADD(y, y, z);
    DCXASSIGN(zk, z.re, z.im);
    while(DCXNORM(err)>DBL_EPSILON*DBL_EPSILON) { //convergence criterion
      k=k+1;
      DCXMUL(zk, zk, z, tmp);
      err.re=zk.re/(k*k);
      err.im=zk.im/(k*k);
      DCXADD(y, y, err);
    }
  }
  if(switchsign) {
    y.re=-y.re;
    y.im=-y.im;
  }
  return y;
}

//additional panel parameters only used here
void fillpanelstruct(panelstruct *output, const panel* panels, double sigma)
{
  double tmp;
  output->xdist=panels->x2-panels->x1;
  output->ydist=panels->y2-panels->y1;
  DCXASSIGN(output->G0,panels->strength1,panels->istrength1);
  DCXASSIGN(output->G1,panels->strength2,panels->istrength2);
  DCXASSIGN(output->rot, panels->costheta, -panels->sintheta);
  DCXSQR(output->rot2, output->rot, tmp);
  output->ks=0.5/(panels->lambda*sigma);
}

//calculates the panel velocity giving separate results for starting and ending strength
void calculatepanelvelocity(dcx* res0, dcx* res1,dcx* res0c, dcx* res1c, double er, double ei, const panel *panels,panelstruct *panels2,double sigma)
{
  double  xdist, ydist, xdist2, ydist2, ydistc, ydist2c, rotydist,em2z0, em2z1, tmp;
  dcx zmz0, zmz1, zmz0c, zmz1c, ei2zmz0, ei2zmz1, ei2zmz0c, ei2zmz1c, ctmp, lsinh0,
          lsinh0c, lsinh1, lsinh1c, res, part, part2, localrot;
  int switchsign;
  xdist=er-panels->x1;
  ydist=ei-panels->y1;
  xdist2=er-panels->x2;
  ydist2=ei-panels->y2;
  ydistc=ei+panels->y1;
  ydist2c=ei+panels->y2;
  rotydist=ydist*panels->costheta-xdist*panels->sintheta;
  if(fabs(rotydist)<panels->lambda*1e-12) { //check if this works
    tmp=(2*panels->side-1)*panels->lambda*1e-12*panels->costheta;
    ydist+=tmp;
    ydist2+=tmp;
    tmp=(2*panels->side-1)*panels->lambda*1e-12*panels->sintheta;
    xdist+=tmp;
    xdist2+=tmp;
  }
  switchsign=0;
  DCXASSIGN(localrot, panels2->rot.re, panels2->rot.im);
  if(xdist<0||xdist2<0) {
    //note: this code only checks for branch cut problems for the actual position,
    //not for the conjugate one, meaning that it only gives correct values for the integral inside the channel
    if(ydist<0&&ydist2>0||ydist>0&&ydist2<0) {//integration over branch cut possible
      tmp=(xdist2-xdist)/(ydist2-ydist);
      if(xdist2-ydist2*tmp<0) { //switch sign to avoid problems with branch cut
        xdist=-xdist;
        ydist=-ydist;
        xdist2=-xdist2;
        ydist2=-ydist2;
        ydistc=-ydistc;
        ydist2c=-ydist2c;
        switchsign=1;
        DCXASSIGN(localrot, -panels2->rot.re, -panels2->rot.im);
      }
    }
  }

  DCXASSIGN(zmz0, sigma*xdist, sigma*ydist);
  DCXASSIGN(zmz1, sigma*xdist2, sigma*ydist2);
  DCXASSIGN(zmz0c, zmz0.re, sigma*ydistc);
  DCXASSIGN(zmz1c, zmz1.re, sigma*ydist2c);

  em2z0=exp(-2*zmz0.re);
  DCXASSIGN(ei2zmz0, em2z0*cos(2*zmz0.im), -em2z0*sin(2*zmz0.im));
  DCXASSIGN(ei2zmz0c, em2z0*cos(2*zmz0c.im), -em2z0*sin(2*zmz0c.im));

  DCXASSIGN(ctmp, 1-ei2zmz0.re, -ei2zmz0.im);
  DCXLOG(lsinh0, ctmp)
  DCXADD(lsinh0, lsinh0, zmz0);

  DCXASSIGN(ctmp, 1-ei2zmz0c.re, -ei2zmz0c.im);
  DCXLOG(lsinh0c, ctmp)
  DCXADD(lsinh0c, lsinh0c, zmz0c);

  em2z1=exp(-2*zmz1.re);
  DCXASSIGN(ei2zmz1, em2z1*cos(2*zmz1.im), -em2z1*sin(2*zmz1.im));
  DCXASSIGN(ei2zmz1c, em2z1*cos(2*zmz1c.im), -em2z1*sin(2*zmz1c.im));

  DCXASSIGN(ctmp, 1-ei2zmz1.re, -ei2zmz1.im);
  DCXLOG(lsinh1, ctmp)
  DCXADD(lsinh1, lsinh1, zmz1);

  DCXASSIGN(ctmp, 1-ei2zmz1c.re, -ei2zmz1c.im);
  DCXLOG(lsinh1c, ctmp)
  DCXADD(lsinh1c, lsinh1c, zmz1c);

  DCXMUL(*res0,localrot,lsinh0,tmp);
  DCXMULN(*res1,localrot,lsinh1,tmp);

  DCXCMUL(*res0c,localrot,lsinh0c,tmp);

  DCXCMULN(*res1c,localrot,lsinh1c,tmp);

  part=polylog2(ei2zmz0.re, ei2zmz0.im);
  ctmp=polylog2(ei2zmz1.re, ei2zmz1.im);

  DCXSUB(part, part, ctmp);

  DCXSQR(part2, zmz0, tmp);
  DCXSQR(ctmp, zmz1, tmp);

  DCXSUB(part2, part2, ctmp);
  DCXADD(part, part, part2);

  DCXMUL(ctmp, panels2->rot2, part, tmp);

  DCXADDC(*res1, *res1, ctmp, panels2->ks);
  DCXSUBC(*res0, *res0, ctmp, panels2->ks);

  part=polylog2(ei2zmz0c.re, ei2zmz0c.im);
  ctmp=polylog2(ei2zmz1c.re, ei2zmz1c.im);

  DCXSUB(part, part, ctmp);

  DCXSQR(part2, zmz0c, tmp);
  DCXSQR(ctmp, zmz1c, tmp);
  DCXSUB(part2, part2, ctmp);
  DCXADD(part, part, part2);
  DCXCMUL(ctmp, panels2->rot2, part, tmp);

  DCXADDC(*res1c, *res1c, ctmp, panels2->ks);
  DCXSUBC(*res0c, *res0c, ctmp, panels2->ks);

  if(switchsign) {
    res0->im=-res0->im;
    res0->re=-res0->re;
    res1->im=-res1->im;
    res1->re=-res1->re;
    res0c->im=-res0c->im;
    res0c->re=-res0c->re;
    res1c->im=-res1c->im;
    res1c->re=-res1c->re;
  }
}

#endif /*CHANNELPOT*/
