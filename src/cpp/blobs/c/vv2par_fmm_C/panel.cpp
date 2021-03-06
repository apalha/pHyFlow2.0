#include "panel.h"

#ifndef C_CODE          /*if C_CODE is not defined, we are in MEX file mode*/
#include "mex.h"
#else                   /*if in C file mode, redefine all mex functions to c functions*/
#define mxFree free
#define mxMalloc malloc
#define mexPrintf printf
#define mxLogical int
#define mxArray double
#include <stdio.h>
#endif

/* panel.cpp */
/* A. Goude 2010-01-01 */

int panelinbox(panel* outpanel,const panel* inpanel,dcmplx z0,dcmplx d0)
/*
creates a new panel structure consisting of the panel inside the box.
Returns 1 if panel is in box and 0 in the case the panel is outside the box.
 *outpanel is the new panel and inpanel is the input. z0 and d0 are the box values
*/
{
  double dist, kstrength, kistrength;
  double xmax=creal(z0)+creal(d0);
  double xmin=creal(z0)-creal(d0);
  double ymax=cimag(z0)+cimag(d0);
  double ymin=cimag(z0)-cimag(d0);
  #ifdef DEBUGpanelinbox
  mexPrintf("x1=%f x2=%f y1=%f y2=%f\nxmin=%f xmax=%f ymin=%f ymax=%f\n", inpanel->x1, inpanel->x2, inpanel->y1, inpanel->y2, xmin, xmax, ymin, ymax);
  #endif
  *outpanel=*inpanel;
  if(inpanel->x1>inpanel->x2) {
    #ifdef DEBUGpanelinbox
    mexPrintf("x1 larger than x2\n");
    #endif
    if(inpanel->x2>=xmax) {/*Panel outside box*/
      if((inpanel->x1-xmax)/fmax(creal(d0), fabs(creal(z0)))>1e-14) { /*if it only is slightly outside, This can be due to roundoff errors, therefore, accept it as inside*/
        outpanel->x1=0;
        outpanel->x2=0;
        outpanel->y1=0;
        outpanel->y2=0;
        #ifdef DEBUGpanelinbox
        mexPrintf("x2 larger than xmax,panel outside box\n");
        #endif
        return 0;
      }
    }
    else if(inpanel->x1<=xmin) {/*Panel outside box*/
      if((inpanel->x2-xmin)/fmax(creal(d0), fabs(creal(z0)))<-1e-14) { /*if it only is slightly outside, This can be due to roundoff errors, therefore, accept it as inside*/
        outpanel->x1=0;
        outpanel->x2=0;
        outpanel->y1=0;
        outpanel->y2=0;
        #ifdef DEBUGpanelinbox
        mexPrintf("x1 smaller than xmin,panel outside box\n");
        #endif
        return 0;
      }
    }
    else {
      if(inpanel->x1>xmax) {
        outpanel->y1=inpanel->y1+inpanel->k1*(xmax-inpanel->x1);
        outpanel->x1=xmax;
      }
      if(inpanel->x2<xmin) {
        outpanel->y2=inpanel->y2+inpanel->k1*(xmin-inpanel->x2);
        outpanel->x2=xmin;
      }
    }
  }
  else {
    #ifdef DEBUGpanelinbox
    mexPrintf("x2 larger than x1\n");
    #endif
    if(inpanel->x1>xmax) {/*Panel outside box, Note: if x1==x2=xmax, consider panel inside box*/
      if((inpanel->x2-xmax)/fmax(creal(d0), fabs(creal(z0)))>1e-14) { /*if it only is slightly outside, This can be due to roundoff errors, therefore, accept it as inside*/
        outpanel->x1=0;
        outpanel->x2=0;
        outpanel->y1=0;
        outpanel->y2=0;
        #ifdef DEBUGpanelinbox
        mexPrintf("x1 larger than xmax,panel outside box\n");
        #endif
        return 0;
      }
    }
    else if(inpanel->x2<xmin) {/*Panel outside box*/
      if((inpanel->x1-xmin)/fmax(creal(d0), fabs(creal(z0)))<-1e-14) { /*if it only is slightly outside, This can be due to roundoff errors, therefore, accept it as inside*/
        outpanel->x1=0;
        outpanel->x2=0;
        outpanel->y1=0;
        outpanel->y2=0;
        #ifdef DEBUGpanelinbox
        mexPrintf("x2 smaller than xmin,panel outside box\n");
        #endif
        return 0;
      }
    }
    else {
      if(inpanel->x2>xmax) {
        outpanel->y2=inpanel->y2+inpanel->k1*(xmax-inpanel->x2);
        outpanel->x2=xmax;
        #ifdef DEBUGpanelinbox
        mexPrintf("decreased y2 to %f\n", outpanel->y2);
        #endif
      }
      if(inpanel->x1<xmin) {
        outpanel->y1=inpanel->y1+inpanel->k1*(xmin-inpanel->x1);
        outpanel->x1=xmin;
        #ifdef DEBUGpanelinbox
        mexPrintf("decreased y1 to %f\n", outpanel->y1);
        #endif
      }
    }
  }
  #ifdef DEBUGpanelinbox
  mexPrintf("x1=%f x2=%f y1=%f y2=%f\n", outpanel->x1, outpanel->x2, outpanel->y1, outpanel->y2);
  #endif
  if(outpanel->y1>outpanel->y2) {
    #ifdef DEBUGpanelinbox
    mexPrintf("y1 larger than y2\n");
    #endif
    if(outpanel->y2>=ymax) {/*Panel outside box*/
      if((inpanel->y1-ymax)/fmax(cimag(d0), fabs(cimag(z0)))>1e-14) { /*if it only is slightly outside, This can be due to roundoff errors, therefore, accept it as inside*/
        outpanel->x1=0;
        outpanel->x2=0;
        outpanel->y1=0;
        outpanel->y2=0;
        #ifdef DEBUGpanelinbox
        mexPrintf("y2 larger than ymax,panel outside box\n");
        #endif
        return 0;
      }
    }
    else if(outpanel->y1<=ymin) {/*Panel outside box*/
      if((inpanel->y2-ymin)/fmax(cimag(d0), fabs(cimag(z0)))<-1e-14) { /*if it only is slightly outside, This can be due to roundoff errors, therefore, accept it as inside*/
        outpanel->x1=0;
        outpanel->x2=0;
        outpanel->y1=0;
        outpanel->y2=0;
        #ifdef DEBUGpanelinbox
        mexPrintf("y1 smaller than ymin,panel outside box\n");
        #endif
        return 0;
      }
    }
    else {
      if(outpanel->y1>=ymax) {
        outpanel->x1=outpanel->x1+outpanel->k2*(ymax-outpanel->y1);
        outpanel->y1=ymax;
        #ifdef DEBUGpanelinbox
        mexPrintf("decreased x1 to %f\n", outpanel->x1);
        #endif
      }
      if(outpanel->y2<=ymin) {
        outpanel->x2=outpanel->x2+outpanel->k2*(ymin-outpanel->y2);
        outpanel->y2=ymin;
        #ifdef DEBUGpanelinbox
        mexPrintf("decreased x2 to %f\n", outpanel->x2);
        #endif
      }
    }
  }
  else {
    #ifdef DEBUGpanelinbox
    mexPrintf("y2 larger than y1\n");
    #endif
    if(outpanel->y1>ymax) {/*Panel outside box*/
      if((inpanel->y2-ymax)/fmax(cimag(d0), fabs(cimag(z0)))>1e-14) { /*if it only is slightly outside, This can be due to roundoff errors, therefore, accept it as inside*/
        outpanel->x1=0;
        outpanel->x2=0;
        outpanel->y1=0;
        outpanel->y2=0;
        #ifdef DEBUGpanelinbox
        mexPrintf("y1 larger than ymax,panel outside box\n");
        #endif
        return 0;
      }
    }
    else if(outpanel->y2<ymin) {/*Panel outside box*/
      if((inpanel->y1-ymin)/fmax(cimag(d0), fabs(cimag(z0)))<-1e-14) { /*if it only is slightly outside, This can be due to roundoff errors, therefore, accept it as inside*/
        outpanel->x1=0;
        outpanel->x2=0;
        outpanel->y1=0;
        outpanel->y2=0;
        #ifdef DEBUGpanelinbox
        mexPrintf("y2 smaller than xmax,panel outside box\n");
        #endif
        return 0;
      }
    }
    else {
      if(outpanel->y2>ymax) {
        outpanel->x2=outpanel->x2+outpanel->k2*(ymax-outpanel->y2);
        outpanel->y2=ymax;
        #ifdef DEBUGpanelinbox
        mexPrintf("decreased x2 to %f\n", outpanel->x2);
        #endif
      }
      if(outpanel->y1<ymin) {
        outpanel->x1=outpanel->x1+outpanel->k2*(ymin-outpanel->y1);
        outpanel->y1=ymin;
        #ifdef DEBUGpanelinbox
        mexPrintf("decreased x1 to %f\n", outpanel->x1);
        #endif
      }
    }
  }
  if(inpanel->horizontal) {/*use the direction with smallest k for numerical stability*/
    dist=inpanel->x2-inpanel->x1;
    kstrength=(inpanel->strength2-inpanel->strength1)/dist;
    kistrength=(inpanel->istrength2-inpanel->istrength1)/dist;
    outpanel->strength2=inpanel->strength2+kstrength*(outpanel->x2-inpanel->x2);
    outpanel->strength1=inpanel->strength1+kstrength*(outpanel->x1-inpanel->x1);
    outpanel->istrength2=inpanel->istrength2+kistrength*(outpanel->x2-inpanel->x2);
    outpanel->istrength1=inpanel->istrength1+kistrength*(outpanel->x1-inpanel->x1);
  }
  else {
    dist=inpanel->y2-inpanel->y1;
    kstrength=(inpanel->strength2-inpanel->strength1)/dist;
    kistrength=(inpanel->istrength2-inpanel->istrength1)/dist;
    outpanel->strength2=inpanel->strength2+kstrength*(outpanel->y2-inpanel->y2);
    outpanel->strength1=inpanel->strength1+kstrength*(outpanel->y1-inpanel->y1);
    outpanel->istrength2=inpanel->istrength2+kistrength*(outpanel->y2-inpanel->y2);
    outpanel->istrength1=inpanel->istrength1+kistrength*(outpanel->y1-inpanel->y1);
  }
  outpanel->lambda=hypot(outpanel->x2-outpanel->x1, outpanel->y2-outpanel->y1);
  return 1;
}
void readpanels(panel* panels,const double* realpart,const double* imagpart,const double* strength,const double* istrength,mxLogical* side,int* smoother,double* cutoff,int count)
/*converts the three matrices realpart, imagpart and strength into panel structure. Should work will imagpart==NULL or istrength==NULL
 *The input should be consistent with a 2xN matrix*/
{
  int i;
  double xdist, ydist;
  for(i=0;i<count;i++) {
    panels[i].x1=realpart[2*i];
    panels[i].x2=realpart[2*i+1];
    if(imagpart!=NULL) { //if input was real
      panels[i].y1=imagpart[2*i];
      panels[i].y2=imagpart[2*i+1];
    }
    else {
      panels[i].y1=0;
      panels[i].y2=0;
    }
    panels[i].strength1=strength[2*i];
    panels[i].strength2=strength[2*i+1];
    xdist=panels[i].x2-panels[i].x1;
    ydist=panels[i].y2-panels[i].y1;
    if(xdist==0&&ydist==0) {  //avoid divide by zero
      panels[i].horizontal=1;
      panels[i].k1=0;
      panels[i].k2=0;
      panels[i].theta=0;
    }
    else {
      panels[i].k1=ydist/xdist;
      panels[i].k2=xdist/ydist;
      panels[i].theta=atan2(ydist, xdist);
    }
    panels[i].horizontal=(fabs(xdist)>fabs(ydist));
    panels[i].sintheta=sin(panels[i].theta);
    panels[i].costheta=cos(panels[i].theta);
    panels[i].lambda=hypot(xdist, ydist);
    panels[i].side=side[i];
    panels[i].istrength1=0;
    panels[i].istrength2=0;
    if(istrength!=NULL) { //if input was real
      panels[i].istrength1=istrength[2*i];
      panels[i].istrength2=istrength[2*i+1];
    }
    panels[i].cutoff=cutoff[i];
    panels[i].smoother=smoother[i];
  }
}
//This function is currently not used. Only available for debugging purposes
void writepanels(const panel* panels,double* realpart,double* imagpart,double* strength,double* istrength,mxLogical* side,int count)
/*The inverse of readpanels, converts the panel structure back to realpart,imagpart and strength*/
{
  int i;
  for(i=0;i<count;i++) {
    realpart[2*i]=panels[i].x1;
    realpart[2*i+1]=panels[i].x2;
    if(imagpart!=NULL) {
      imagpart[2*i]=panels[i].y1;
      imagpart[2*i+1]=panels[i].y2;
    }
    strength[2*i]=panels[i].strength1;
    strength[2*i+1]=panels[i].strength2;
    side[i]=panels[i].side;
    if(istrength!=NULL) {
      istrength[2*i]=panels[i].istrength1;
      istrength[2*i+1]=panels[i].istrength2;
    }
  }
}
/*------------------------------------------------------------------------*/
//debugging function only
void printpanels(const panel* panels,int count)
/*prints the panel structure*/
{
  int i;
  for(i=0;i<count;i++) {
    mexPrintf("Panel %d: Start position %f %f, End position %f %f Strength %f+%fi %f+%fi\nk1=%f k2=%f theta=%f lambda=%f horizontal=%d\n", i, panels[i].x1, panels[i].y1, panels[i].x2, panels[i].y2, panels[i].strength1, panels[i].istrength1, panels[i].strength2, panels[i].istrength2, panels[i].k1, panels[i].k2, panels[i].theta, panels[i].lambda, panels[i].horizontal);
  }
}
/*------------------------------------------------------------------------*/
void expandpanel(const panel* panels,dcmplx z0,dcmplx *coeff,int p)
/*Perform multipole expansion of the panels corresponds to mpexp_init
 *panels is the constant panel array, z0 the center of the expansion
 *coeff is where to store the coefficients and p the number of coefficients*/
{
  const dcmplx r1 = (panels->x1+I*panels->y1)-z0;
  const dcmplx r2 = (panels->x2+I*panels->y2)-z0;
  const dcmplx gamma1=(panels->strength1+I*panels->istrength1)*(panels->costheta-I*panels->sintheta);
  const dcmplx gamma2=(panels->strength2+I*panels->istrength2)*(panels->costheta-I*panels->sintheta);
//    const dcmplx gamma1=panels->strength1*(panels->costheta-I*panels->sintheta);
//    const dcmplx gamma2=panels->strength2*(panels->costheta-I*panels->sintheta);
  const dcmplx lambda=panels->x2+I*panels->y2-panels->x1-I*panels->y1;
  const dcmplx m =(gamma2-gamma1)/lambda;
  if(creal(lambda)==0&&cimag(lambda)==0) return; //safety check. Can in rare occations happen due to numerical errors. With lambda==0, panel is unimportant anyway
  /*
    mexPrintf("r1=%f +%fi r2=%f+%fi gamma1=%f+%fi gamma2=%f+%fi lambda=%f+%fi p=%d\n",creal(r1),cimag(r1),creal(r2),cimag(r2),creal(gamma1),cimag(gamma1),creal(gamma2),cimag(gamma2),creal(lambda),cimag(lambda),p);
   */
  dcmplx z1=-m*r1*r1;
  dcmplx z2=m*r2*r2;
  dcmplx zg1=gamma1*r1;
  dcmplx zg2=-gamma2*r2;
  /*
     mexPrintf("z1=%f +%fi z2=%f+%fi zg1=%f+%fi zg2=%f+%fi\n",creal(z1),cimag(z1),creal(z2),cimag(z2),creal(zg1),cimag(zg1),creal(zg2),cimag(zg2));
   */
//    if(panels->strength1!=0||panels->strength2!=0){
//                mexPrintf("Expanding panel\n");
//                printpanels(panels,1);
//            }
  int j;
  for (j = 1; j <= p; j++, z1*=r1, z2*=r2, zg1*=r1, zg2*=r2) {
    /*
        mexPrintf("j=%d,z1=%f +%fi z2=%f+%fi zg1=%f+%fi zg2=%f+%fi\n",j,creal(z1/(j+2)),cimag(z1/(j+2)),creal(z2/(j+2)),cimag(z2/(j+2)),creal(zg1),cimag(zg1),creal(zg2),cimag(zg2));
     */
    coeff[j] += 1.0/j*(z1/(j+1.0)+z2/(j+1.0)+zg1+zg2);
//        mexPrintf("New value=%e+%e,added %e+%ei old value=%e+%ei\n",coeff[j],1.0/j*(z1/(j+1.0)+z2/(j+1.0)+zg1+zg2),coeff[j]-1.0/j*(z1/(j+1.0)+z2/(j+1.0)+zg1+zg2));
    /*
        mexPrintf("New value=%f+%fi Step =%f+%fi\n",creal(1/(j+1)*(z1/(j+2)+z2/(j+2)+zg1+zg2)),cimag(1/(j+1)*(z1/(j+2)+z2/(j+2)+zg1+zg2)),creal(z1/(j+2)+z2/(j+2)+zg1+zg2),cimag(z1/(j+2)+z2/(j+2)+zg1+zg2));
     */
  }
}
void farexpandpanel(const panel* panels,dcmplx z0,dcmplx *coeff,int p)
/*Performs multipole expansion of panel. Corresponds to mpexp_initp
 *panels is the constant panel array, z0 the center of the expansion
 *coeff is where to store the coefficients and p the number of coefficients
 */
{
  const dcmplx r1 = (panels->x1+I*panels->y1)-z0;
  const dcmplx r2 = (panels->x2+I*panels->y2)-z0;
  const dcmplx gamma1=(panels->strength1+I*panels->istrength1)*(panels->costheta-I*panels->sintheta);
  const dcmplx gamma2=(panels->strength2+I*panels->istrength2)*(panels->costheta-I*panels->sintheta);
//    const dcmplx gamma1=panels->strength1*(panels->costheta-I*panels->sintheta);
//    const dcmplx gamma2=panels->strength2*(panels->costheta-I*panels->sintheta);
  const dcmplx invr1=1.0/r1;
  const dcmplx invr2=1.0/r2;
  const dcmplx lambda=panels->x2+I*panels->y2-panels->x1-I*panels->y1;
  const dcmplx m =(gamma2-gamma1)/lambda;

  if(creal(lambda)==0&&cimag(lambda)==0) return;
  dcmplx logdist1, logdist2;
  dcmplx z1=m*invr1;
  dcmplx z2=-m*invr2;
  dcmplx zg1=gamma1*invr1;
  dcmplx zg2=-gamma2*invr2;
  if(creal(r1) > 0 && creal(r2) > 0) { // to prevent that panels cross the branch cut of atan2. Not fully safe, but can only fail for panels belonging to three different quadrants, which should not occur in the far field
    logdist1=I*PI+clog(r1);
    logdist2=I*PI+clog(r2);
  }
  else {
    logdist1=clog(-r1);
    logdist2=clog(-r2);
  }
  if(p<1) //safety
    return;
//    mexPrintf("r1=%f+%fi r2=%f+%fi\ncoeff[0]=%e ",creal(r1),cimag(r1),creal(r2),cimag(r2),coeff[0]);
  coeff[0]+=m*(r1*logdist1-r2*logdist2)+gamma2*logdist2-gamma1*logdist1+gamma2-gamma1;
//    mexPrintf(": %e\n",coeff[0]);
  if(p<2) //can p really be This small???
    return;
//    mexPrintf("coeff[1]=%e ",coeff[1]);
  coeff[1]+=m*(logdist2-logdist1)+zg2+zg1;
//    mexPrintf(": %e\n",coeff[1]);
  zg1*=invr1;
  zg2*=invr2;

  int j;
  for (j = 2; j <= p; j++, z1*=invr1, z2*=invr2, zg1*=invr1, zg2*=invr2) {
//        mexPrintf("coeff[%d]=%e+%ei",j,coeff[j]);
    coeff[j] += 1.0/j*(1.0/(j-1)*(z1+z2)+zg1+zg2);
//        mexPrintf(": %e+%ei added %e+%ei\n",coeff[j],1.0/j*(1.0/(j-1)*(z1+z2)+zg1+zg2));
  }
}

void MPexp_box_panel(MPexp *This,const panel *panels,int panelcount, int* validinput)
//the panel version of MPexp_box_
//This is the box, panels is the panel array, and panelcount its length
{
  for(int i=0;i<panelcount;i++) {
    if(panels[i].x1>This->xmax)
      This->xmax=panels[i].x1;
    if(panels[i].x1<This->xmin)
      This->xmin=panels[i].x1;
    if(panels[i].x2>This->xmax)
      This->xmax=panels[i].x2;
    if(panels[i].x2<This->xmin)
      This->xmin=panels[i].x2;
    if(panels[i].y1>This->ymax)
      This->ymax=panels[i].y1;
    if(panels[i].y1<This->ymin)
      This->ymin=panels[i].y1;
    if(panels[i].y2>This->ymax)
      This->ymax=panels[i].y2;
    if(panels[i].y2<This->ymin)
      This->ymin=panels[i].y2;
#ifdef CHECKNANINPUT
    if(!isfinite(panels[i].x1)||!isfinite(panels[i].x2)||!isfinite(panels[i].y1)||!isfinite(panels[i].y2))
        *validinput=0;
#endif
  }
}
void MPexp_split_panel(dcmplx z0,dcmplx d0,const panel *panels,const int *inptr, int *outptr1,int* outptr2,int incnt,int* outcnt1,int* outcnt2,double splitpoint,int xsplit)
/*Using a defines split point, This function splits the panels of the box.
 *A panel can be includen in both boxes.
 *z0 and d0 are box coordinates
 *panels are the constant panel array
 *inptr contains the indices in panels and incnt the length of inptr
 *outptr1 and outptr2 returns the indices of the panels for each box, and outcnt1 and outcnt2 their length
 *splitpoint the point of the split
 *xsplit==true means split in x-direction, false means y-direction*/
{
  panel smallpanel;
  *outcnt1=0;
  *outcnt2=0;
  if(xsplit) {
    for(int k=0;k<incnt;k++) {
      panelinbox(&smallpanel, panels+inptr[k], z0, d0);
//            printpanels(&smallpanel,1);
      if(smallpanel.x1<splitpoint||smallpanel.x2<splitpoint)
        outptr1[(*outcnt1)++]=inptr[k];
      if(smallpanel.x1>splitpoint||smallpanel.x2>splitpoint)
        outptr2[(*outcnt2)++]=inptr[k];
      if(smallpanel.x1==splitpoint&&smallpanel.x2==splitpoint) /*if both points are on the split point, pick one box*/
        outptr1[(*outcnt1)++]=inptr[k];

    }
  }
  else {
    for(int k=0;k<incnt;k++) {
      panelinbox(&smallpanel, panels+inptr[k], z0, d0);
//            printpanels(&smallpanel,1);
      if(smallpanel.y1<splitpoint||smallpanel.y2<splitpoint)
        outptr1[(*outcnt1)++]=inptr[k];
      if(smallpanel.y1>splitpoint||smallpanel.y2>splitpoint)
        outptr2[(*outcnt2)++]=inptr[k];
      if(smallpanel.y1==splitpoint&&smallpanel.y2==splitpoint) /*if both points are on the split point, pick one box*/
        outptr1[(*outcnt1)++]=inptr[k];
    }
  }
}

#ifndef SWAPI
#define SWAPI(X, Y) { int tmpi = X; X = Y; Y = tmpi; }
#endif
#ifndef SWAPD
#define SWAPD(X, Y) { double tmpd = X; X = Y; Y = tmpd; }
#endif

#ifdef PANELSORT
//debug functions
void printindexlist(const double *z,const int *ix,int begin, int end)
{
  for(int i=begin;i<end;i++) {
    mexPrintf("i=%d ix[i]=%d z[ix[i]]=%.16e\n", i, ix[i], z[ix[i]]);
  }
}
void printdummylist(const double *dummy,int begin,int end)
{
  for(int i=begin;i<end;i++) {
    mexPrintf("i=%d dummy[i]=%.16e\n", i, dummy[i]);
  }
}
// the new version of MPexp_partition. Has an extra input array of dummy vortices
void MPexp_partition_dummy(int begin, int *im, int end, double *z0,
                           int *ix, const double *z,
                           int Nmax,double *dummy,int Ndummy)
/* Partitions the set of points z[ix[begin]] through z[ix[end]]
 * (exclusive) in two sets indicated by [begin,im) and [im,end) such
 * that all points are <= or >= than z0 (which is both input and
 * output). The division is done such that Nlarge <= Nmax, where
 * Nlarge is the larger number of particles in the two lists. It is
 * assumed that This split is possible and hence that Nmax is not too
 * small.
 *
 * The routine is a workhorse and performance is important; for each
 * box which is not a leaf it will be called three times.
 *
 * This modified version includes dummy points in the partitioning
 * to modify the split point. The dummy list will be modified. */
{
  int il, ir, left = begin, right = end;
  int Nsmall, Nlarge;
  int dummyleft=0, dummyright=Ndummy, dummyil, dummyir, dummyNsmall, dummyNlarge, dummysplit;

  // first split according to input z0 (geometric midpoint)
  const double z0_ = *z0; // compiler: This is really a constant!
  il = left-1; ir = right;
  do il++; while (il < right && z[ix[il]] < z0_);
  do ir--; while (left <= ir && z0_ < z[ix[ir]]);
  if (il < ir) {
    SWAPI(ix[il], ix[ir])
    /* after the first swap, there is no need to check for index out
     * of bounds */
    for ( ; ; ) {
      do il++; while (z[ix[il]] < z0_);
      do ir--; while (z0_ < z[ix[ir]]);
      if (il >= ir) break;
      SWAPI(ix[il], ix[ir])
    }
  }

  // Same split for the dummy points
  dummyil = dummyleft-1; dummyir = dummyright;
  do dummyil++; while (dummyil < dummyright && dummy[dummyil] < z0_);
  do dummyir--; while (dummyleft <= dummyir && z0_ < dummy[dummyir]);
  if (dummyil < dummyir) {
    SWAPD(dummy[dummyil], dummy[dummyir])
    /* after the first swap, there is no need to check for index out
     * of bounds */
    for ( ; ; ) {
      do dummyil++; while (dummy[dummyil] < z0_);
      do dummyir--; while (z0_ < dummy[dummyir]);
      if (dummyil >= dummyir) break;
      SWAPD(dummy[dummyil], dummy[dummyir])
    }
  }

  #ifdef DEBUGPANELSORT
  mexPrintf("First sort,z0=%e dummyil=%d dummyir=%d il=%d ir=%d\n", z0_, dummyil, dummyir, il, ir);
//    printindexlist(z,ix,begin,end);
//    printdummylist(dummy,0,Ndummy);
  #endif
  // proceed with largest interval [left,il) or [il,right)
  Nsmall = il-begin;
  Nlarge = end-il;
  //and for the dummies
  dummyNsmall=dummyil;
  dummyNlarge=Ndummy-dummyil;
  #ifdef DEBUGPANELSORT
  mexPrintf("Nsmall=%d Nlarge=%d dummyNsmall=%d dummyNlarge=%d\n", Nsmall, Nlarge, dummyNsmall, dummyNlarge);
  #endif
  if (Nsmall+dummyNsmall*PANELDUMMYFACTOR > Nlarge+dummyNlarge*PANELDUMMYFACTOR) {
    Nlarge = Nsmall;
    *im = right = il;
    dummyNlarge=dummyNsmall;
    dummyright=dummyil;
  }
  else {
    *im = left = il;
    dummyleft=dummyil;
  }
  dummysplit=(dummyright-dummyleft)*PANELDUMMYFACTOR>right-left; //determine if it is best to choose split point from dummies or normal vortices from which part that has most undetermined points
  #ifdef DEBUGPANELSORT
  mexPrintf("vortex interval size=%d, dummy interval size=%d\n", right-left, dummyright-dummyleft);
  #endif
  while (Nlarge +dummyNlarge*PANELDUMMYFACTOR > Nmax) {
    if(dummysplit) { //choose split point from dummies
      #ifdef DEBUGPANELSORT
      mexPrintf("dummy split dummyleft=%d dummyright=%d\n", dummyleft, dummyright);
      #endif
      if (dummyright-dummyleft >= 3) {
        /* "median of 3": find pivot and create boundaries such that
         * left <= left+1 <= right-1 */
        const int mid = (dummyleft+dummyright)>>1;
        SWAPD(dummy[dummyleft+1], dummy[mid])
        if (dummy[dummyleft] > dummy[dummyright-1])
          SWAPD(dummy[dummyleft], dummy[dummyright-1]);
        if (dummy[dummyleft+1] > dummy[dummyright-1])
          SWAPD(dummy[dummyleft+1], dummy[dummyright-1])
        else if (dummy[dummyleft] > dummy[dummyleft+1])
          SWAPD(dummy[dummyleft], dummy[dummyleft+1])
          const double z0_ = *z0 = dummy[dummyleft+1];

        for (dummyil = dummyleft+1, dummyir = dummyright-1; ; ) {
          do dummyil++; while (dummy[dummyil] < z0_);
          do dummyir--; while (z0_ < dummy[dummyir]);
          if (dummyil >= dummyir) break;

          SWAPD(dummy[dummyil], dummy[dummyir])
        }
        SWAPD(dummy[dummyleft+1], dummy[dummyir]) // ir is the correct place for the pivot

        /*with known split point, split the vortices*/
        il = left-1; ir = right;
        do il++; while (il < right && z[ix[il]] < z0_);
        do ir--; while (left <= ir && z0_ < z[ix[ir]]);
        if (il < ir) {
          SWAPI(ix[il], ix[ir])
          /* after the first swap, there is no need to check for index out
           * of bounds */
          for ( ; ; ) {
            do il++; while (z[ix[il]] < z0_);
            do ir--; while (z0_ < z[ix[ir]]);
            if (il >= ir) break;
            SWAPI(ix[il], ix[ir])
          }
        }
        #ifdef DEBUGPANELSORT
        mexPrintf("sort,z0=%e dummyil=%d dummyir=%d il=%d ir=%d\n", z0_, dummyil, dummyir, il, ir);
        printindexlist(z, ix, begin, end);
        printdummylist(dummy, 0, Ndummy);
        #endif
        // the pivot belongs to the shortest interval:
        Nsmall = il-begin; // [begin,ir)
        Nlarge = end-il; // [ir+1,end)
        dummyNsmall = dummyir; // [begin,ir)
        dummyNlarge = Ndummy-dummyir-1; // [ir+1,end)
        #ifdef DEBUGPANELSORT
        mexPrintf("Nsmall=%d Nlarge=%d dummyNsmall=%d dummyNlarge=%d\n", Nsmall, Nlarge, dummyNsmall, dummyNlarge);
        #endif
        if (Nsmall+dummyNsmall*PANELDUMMYFACTOR > Nlarge+dummyNlarge*PANELDUMMYFACTOR) {
          Nlarge = Nsmall;
          *im = right = il;
          dummyNlarge=dummyNsmall;
          dummyright=dummyir;
        }
        else {
          *im = left = il;
          dummyleft=dummyir+1;
        }
        dummysplit=(dummyright-dummyleft)*PANELDUMMYFACTOR>right-left;
        #ifdef DEBUGPANELSORT
        mexPrintf("vortex interval size=%d, dummy interval size=%d\nleft=%d right=%d dummyleft=%d dummyright=%d\n", right-left, dummyright-dummyleft, left, right, dummyleft, dummyright);
        #endif
      }
      else { //if 2 or less dummy points, split with normal points instead
        if(right-left<=0) { //if no more vortices exist, consider it done, otherwise, use the vortices as final split
          *z0 = dummy[dummyleft]>dummy[dummyright-1]?dummy[dummyleft]:dummy[dummyright-1];
          break;
        }
        dummysplit=0;
              }
    }
    if(!dummysplit) {
      #ifdef DEBUGPANELSORT
      mexPrintf("normal split\n");
      #endif
      if (right-left >= 3) {
        /* "median of 3": find pivot and create boundaries such that
         * left <= left+1 <= right-1 */
        const int mid = (left+right)>>1;
        SWAPI(ix[left+1], ix[mid])
        if (z[ix[left]] > z[ix[right-1]])
          SWAPI(ix[left], ix[right-1]);
        if (z[ix[left+1]] > z[ix[right-1]])
          SWAPI(ix[left+1], ix[right-1])
        else if (z[ix[left]] > z[ix[left+1]])
          SWAPI(ix[left], ix[left+1])
          const double z0_ = *z0 = z[ix[left+1]];

        for (il = left+1, ir = right-1; ; ) {
          do il++; while (z[ix[il]] < z0_);
          do ir--; while (z0_ < z[ix[ir]]);
          if (il >= ir) break;

          SWAPI(ix[il], ix[ir])
        }
        SWAPI(ix[left+1], ix[ir]) // ir is the correct place for the pivot

        /*normal split of dummies if vortices are the pivot*/
        dummyil = dummyleft-1; dummyir = dummyright;
        do dummyil++; while (dummyil < dummyright && dummy[dummyil] < z0_);
        do dummyir--; while (dummyleft <= dummyir && z0_ < dummy[dummyir]);
        #ifdef DEBUGPANELSORT
        mexPrintf("before loop dummyil=%d,dummyir=%d\n", dummyil, dummyir);
        #endif
        if (dummyil < dummyir) {
          SWAPD(dummy[dummyil], dummy[dummyir])
          /* after the first swap, there is no need to check for index out
           * of bounds */
          for ( ; ; ) {
            do dummyil++; while (dummy[dummyil] < z0_);
            do dummyir--; while (z0_ < dummy[dummyir]);
            if (dummyil >= dummyir) break;
            SWAPD(dummy[dummyil], dummy[dummyir])
          }
        }
        #ifdef DEBUGPANELSORT
        mexPrintf("sort,z0=%e dummyil=%d dummyir=%d il=%d ir=%d\n", z0_, dummyil, dummyir, il, ir);
//                printindexlist(z,ix,begin,end);
//                printdummylist(dummy,0,Ndummy);
        #endif
        // the pivot belongs to the shortest interval:
        Nsmall = ir-begin; // [begin,ir)
        Nlarge = end-ir-1; // [ir+1,end)
        dummyNsmall=dummyil;
        dummyNlarge=Ndummy-dummyil;
        #ifdef DEBUGPANELSORT
        mexPrintf("Nsmall=%d Nlarge=%d dummyNsmall=%d dummyNlarge=%d\n", Nsmall, Nlarge, dummyNsmall, dummyNlarge);
        #endif
        if (Nsmall+dummyNsmall*PANELDUMMYFACTOR > Nlarge+dummyNlarge*PANELDUMMYFACTOR) {
//                    dummysplit=dummyNsmall*Nlarge>dummyNlarge*Nsmall;
          Nlarge = Nsmall;
          *im = right = ir;
          dummyNlarge=dummyNsmall;
          dummyright=dummyir;
        }
        else {
          /* note: if Nsmall == Nlarge (can only happen if right-left is
           * odd), then the right array will be shorter than the left by
           * one */
          *im = left = ir+1;
//                    dummysplit=dummyNlarge*Nsmall>dummyNsmall*Nlarge;
          dummyleft=dummyil;
        }
        dummysplit=(dummyright-dummyleft)*PANELDUMMYFACTOR>right-left;
        #ifdef DEBUGPANELSORT
        mexPrintf("vortex interval size=%d, dummy interval size=%d\nleft=%d right=%d dummyleft=%d dummyright=%d\n", right-left, dummyright-dummyleft, left, right, dummyleft, dummyright);
        #endif
      }
      else {
        // code above does not work properly for less than three elements
        // the dummies does not matter here
        #ifdef DEBUGPANELSORT
        mexPrintf("Final step for the vortices, left=%d right=%d, left side=%d right side=%d\n", left, right, left-begin, end-right);
        #endif
        if (z[ix[left]] > z[ix[right-1]])
            SWAPI(ix[left], ix[right-1])
            *im = right-1;
        *z0 = z[ix[right-1]];
        break;
      }
    }
  }
  #ifdef DEBUGPANELSORT
  double maxvalue=z[ix[begin]], minvalue=z[ix[begin]];
  for(int k=begin+1;k<end;k++) {
    if(z[ix[k]]<minvalue)
      minvalue=z[ix[k]];
    if(z[ix[k]]>maxvalue)
      maxvalue=z[ix[k]];
  }
  mexPrintf("Split compltete, maxvalue=%e minvalue=%e,splitpoint=%e\n\n", maxvalue, minvalue, *z0);
    #endif
//    if(fabs(*z0-14)<1e-9) {
//        mexPrintf("dummy list split=%.16e\n",*z0);
//        printindexlist(z,ix,begin,end);
//        printdummylist(dummy,0,Ndummy);
//    }
}
#endif
//
//#ifndef BEGIN
//#define BEGIN(begin, i) (begin == -1 ? i+1 : begin)
//#endif
#undef REALDESTINATION
#define REALDESTINATION qr[jx[j]]
#undef IMAGDESTINATION
#define IMAGDESTINATION qi[jx[j]]
#undef EVALPANEL
#define EVALPANEL smallpanel
#undef REALEVAL
#define REALEVAL er[jx[j]]
#undef IMAGEVAL
#define IMAGEVAL ei[jx[j]]
#undef PANELLOOPSTART
#define PANELLOOPSTART begin
#undef PANELLOOPEND
#define PANELLOOPEND end
void mpexp_directInteractPanel(const panel *panels,const mpexp *This,double *qr, double *qi,const double* er,const double* ei,const int *jx,
        int begin, int end)
/*direct interaction for panels. Only implemented for complex panels, but the speed gain of using real panels is most likely small.*/
 {
  panel smallpanel;
  double /*xpaneldist, ypaneldist,*/ xdist, ydist, sqrdist, xterm, yterm, xtmp, ytmp, ax, ay, bx, by, rotxdist, rotydist, logxterm, logyterm,xtmp0,ytmp0;
//    mexPrintf("directInteractPanel panelcount=%d,begin=%d,end=%d\n",This->npanel,begin,end);
  for(int k=0;k<This->npanel;k++) {
    if(!panelinbox(&smallpanel, panels+This->panelptr[k], This->z0, This->d0)) {
//            mexPrintf("Box start=%fx%f, end %fx%f\n",creal(This->z0)+creal(This->d0),cimag(This->z0)+cimag(This->d0),creal(This->z0)-creal(This->d0),cimag(This->z0)-cimag(This->d0));
//            printpanels(&smallpanel,1);
//            printpanels(panels+This->panelptr[k],1);
      continue;
    }
    if(smallpanel.lambda==0){
//            mexPrintf("Panel %d of zero length detected for box at %f+%fi with size %f+%fi\n",This->panelptr[k],creal(This->z0),cimag(This->z0),creal(This->d0),cimag(This->d0));
//            printpanels(&smallpanel,1);
//            mexPrintf("full panel is\n");
//            printpanels(panels+This->panelptr[k],1);
      continue; //panels of zero length have no effect, and causes NANs.
    }
    //xpaneldist=smallpanel.x2-smallpanel.x1;
    //ypaneldist=smallpanel.y2-smallpanel.y1;
    if(smallpanel.istrength1==0&&smallpanel.istrength2==0) {
      #undef COMPLEXPANELS
      #undef PURECOMPLEXPANELS
      #include "paneldirect.h"
    }
    else {
      if(smallpanel.strength1==0&&smallpanel.strength2==0) {
        #ifndef PURECOMPLEXPANELS
        #define PURECOMPLEXPANELS
        #endif
        #include "paneldirect.h"
      }
      else {
        #undef PURECOMPLEXPANELS
        #ifndef COMPLEXPANELS
        #define COMPLEXPANELS
        #endif
        #include "paneldirect.h"
      }
    }
  }
}

#undef REALDESTINATION
#define REALDESTINATION qr[j]
#undef IMAGDESTINATION
#define IMAGDESTINATION qi[j]
#undef EVALPANEL
#define EVALPANEL panels[k]
#undef REALEVAL
#define REALEVAL er[j]
#undef IMAGEVAL
#define IMAGEVAL ei[j]
#undef PANELLOOPSTART
#define PANELLOOPSTART 0
#undef PANELLOOPEND
#define PANELLOOPEND N
void directInteractPanel(const panel *panels, double *qr, double *qi,const double* er,const double* ei, const int panelcount,const int N)
{
  /*direct interaction for panels. Only implemented for complex panels, but the speed gain of using real panels is most likely small.*/
  double /*xpaneldist, ypaneldist,*/ xdist, ydist, sqrdist, xterm, yterm, xtmp, ytmp, ax, ay, bx, by, rotxdist, rotydist, logxterm, logyterm,xtmp0,ytmp0;
//    mexPrintf("directInteractPanel starting\n");
  for(int k=0;k<panelcount;k++) {
    if(panels[k].lambda==0) continue; //panels of zero length have no effect, and causes NANs.
    //xpaneldist=panels[k].x2-panels[k].x1;
    //ypaneldist=panels[k].y2-panels[k].y1;
  if(panels[k].istrength1==0&&panels[k].istrength2==0) {
      #undef COMPLEXPANELS
      #undef PURECOMPLEXPANELS
//            mexPrintf("Real panel\n");
      #include "paneldirect.h"
    }
    else {
      if(panels[k].strength1==0&&panels[k].strength2==0) {
        #ifndef PURECOMPLEXPANELS
        #define PURECOMPLEXPANELS
        #endif
        #include "paneldirect.h"
      }
      else {
        #undef PURECOMPLEXPANELS
        #ifndef COMPLEXPANELS
        #define COMPLEXPANELS
        #endif
        #include "paneldirect.h"
      }
    }
  }
}
void createdummylist(double* dummy,const panel* panels,const int* indices,int N,int xval,dcmplx z,dcmplx d)
/*creates list of dummy vortices to be used by MPexp_partition_dummy*/
{
  panel smallpanel;
  if(xval) {
    for(int i=0;i<N;i++) {
      panelinbox(&smallpanel, panels+indices[i], z, d);
      dummy[2*i]=smallpanel.x1;
      dummy[2*i+1]=smallpanel.x2;
    }
  }
  else {
    for(int i=0;i<N;i++) {
      panelinbox(&smallpanel, panels+indices[i], z, d);
      dummy[2*i]=smallpanel.y1;
      dummy[2*i+1]=smallpanel.y2;
    }
  }
}
#ifdef PANELSHRINKBOX
void panelshrinkbox(mpexp* This,const double *zr,const double *zi,
        const double *er, const double *ei,
        const int *ix, const int *jx,
        int begin1, int end1,
        int begin2, int end2,
        const panel* panels,panel* smallerpanels,double cutoff)

{
  //Start by initializing all values to the center to prevent too many if statements in case zr==NULL, er==NULL etc
  double xmin=creal(This->z0)+creal(This->d0), xmax=creal(This->z0)-creal(This->d0), ymin=cimag(This->z0)+cimag(This->d0), ymax=cimag(This->z0)-cimag(This->d0);
  #ifdef RADIALSHRINK
  double sqrdist, maxsqrdist=0, xtmp, ytmp;
  #else
  panel smallpanel;
  #endif
  for(int i=begin1;i<end1;i++) {
    if(zr[ix[i]]<xmin) xmin=zr[ix[i]];
    if(zr[ix[i]]>xmax) xmax=zr[ix[i]];
    if(zi[ix[i]]<ymin) ymin=zi[ix[i]];
    if(zi[ix[i]]>ymax) ymax=zi[ix[i]];
  }
  for(int i=begin2;i<end2;i++) {
    if(er[jx[i]]<xmin) xmin=er[jx[i]];
    if(er[jx[i]]>xmax) xmax=er[jx[i]];
    if(ei[jx[i]]<ymin) ymin=ei[jx[i]];
    if(ei[jx[i]]>ymax) ymax=ei[jx[i]];
  }
//     xmin=creal(This->z0)-creal(This->d0);
//     xmax=creal(This->z0)+creal(This->d0);
//     ymin=cimag(This->z0)-cimag(This->d0);
//     ymax=cimag(This->z0)+cimag(This->d0);
  for(int i=0;i<This->npanel;i++) {
    #ifdef RADIALSHRINK
    if(panelinbox(smallerpanels+i, panels+This->panelptr[i], This->z0, This->d0)){
      if(smallerpanels[i].x1<xmin) xmin=smallerpanels[i].x1;
      if(smallerpanels[i].x1>xmax) xmax=smallerpanels[i].x1;
      if(smallerpanels[i].x2<xmin) xmin=smallerpanels[i].x2;
      if(smallerpanels[i].x2>xmax) xmax=smallerpanels[i].x2;
      if(smallerpanels[i].y1<ymin) ymin=smallerpanels[i].y1;
      if(smallerpanels[i].y1>ymax) ymax=smallerpanels[i].y1;
      if(smallerpanels[i].y2<ymin) ymin=smallerpanels[i].y2;
      if(smallerpanels[i].y2>ymax) ymax=smallerpanels[i].y2;
    }
    #else
    if(panelinbox(&smallpanel, panels+This->panelptr[i], This->z0, This->d0)){
      if(smallpanel.x1<xmin) xmin=smallpanel.x1;
      if(smallpanel.x1>xmax) xmax=smallpanel.x1;
      if(smallpanel.x2<xmin) xmin=smallpanel.x2;
      if(smallpanel.x2>xmax) xmax=smallpanel.x2;
      if(smallpanel.y1<ymin) ymin=smallpanel.y1;
      if(smallpanel.y1>ymax) ymax=smallpanel.y1;
      if(smallpanel.y2<ymin) ymin=smallpanel.y2;
      if(smallpanel.y2>ymax) ymax=smallpanel.y2;
    }
    #endif

  }
  COMPLEXASSIGN(This->z0, 0.5*(xmin+xmax), 0.5*(ymin+ymax));
//     This->z0=0.5*(xmin+xmax)+0.5*I*(ymin+ymax);
  #ifdef SHRINKCUTOFFCHECK
  if(xmax-xmin<2*cutoff)
    xmin=xmax-2*cutoff;
  if(ymax-ymin<2*cutoff)
    ymin=ymax-2*cutoff;
  #endif
//     This->d0=0.5*(xmax-xmin)+0.5*I*(ymax-ymin);
  COMPLEXASSIGN(This->d0, 0.5*(xmax-xmin), 0.5*(ymax-ymin));
  #ifdef RADIALSHRINK //if radialshrink should be used, make an extra pass and calculate the largest distance
  //uses square distances to avoid extensive calling of sqrt()
  if(end1-begin1+end2-begin2+This->npanel<RADIALSHRINKLIMIT) {
    for(int i=begin1;i<end1;i++) {
      xtmp=zr[ix[i]]-creal(This->z0);
      ytmp=zi[ix[i]]-cimag(This->z0);
      sqrdist=xtmp*xtmp+ytmp*ytmp;
      if(sqrdist>maxsqrdist) maxsqrdist=sqrdist;
    }
    for(int i=begin2;i<end2;i++) {
      xtmp=er[jx[i]]-creal(This->z0);
      ytmp=ei[jx[i]]-cimag(This->z0);
      sqrdist=xtmp*xtmp+ytmp*ytmp;
      if(sqrdist>maxsqrdist) maxsqrdist=sqrdist;
    }
    for(int i=0;i<This->npanel;i++) {
      xtmp=smallerpanels[i].x1-creal(This->z0);
      ytmp=smallerpanels[i].y1-cimag(This->z0);
      sqrdist=xtmp*xtmp+ytmp*ytmp;
      if(sqrdist>maxsqrdist) maxsqrdist=sqrdist;
      xtmp=smallerpanels[i].x2-creal(This->z0);
      ytmp=smallerpanels[i].y2-cimag(This->z0);
      sqrdist=xtmp*xtmp+ytmp*ytmp;
      if(sqrdist>maxsqrdist) maxsqrdist=sqrdist;
    }
    This->absd0=sqrt(maxsqrdist);
  }
  else
    This->absd0=hypot(0.5*(xmax-xmin), 0.5*(ymax-ymin));
  #ifdef SHRINKCUTOFFCHECK
  if(This->absd0<cutoff)
    This->absd0=cutoff;
  #endif
  #endif
//    mexPrintf("Box after shrink: z0=%f+%fi, d0=%f+%fi\n",This->z0,This->d0);
}
#endif
