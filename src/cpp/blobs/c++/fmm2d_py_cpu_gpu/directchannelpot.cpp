#include <stdio.h>
#include "fmm.h"
#include "channelpotpanel.h"
#ifdef CHANNELPOT
void zchannel2(int N,
           const double *zr,const double *zi,
           const double *mr,const double *mi,
           double *pr,double *pi,
           SMOOTHER smooth,double xopt,double cutoff,bool cont,double channelheight);
void rchannel2(int N,
           const double *zr,const double *zi,
           const double *mr,
           double *pr,double *pi,
           SMOOTHER smooth,double xopt,double cutoff,bool cont,double channelheight);
void zchannel(int N,
           const double *zr,const double *zi,
           const double *mr,const double *mi,
           int NE,
           const double *er,const double *ei,
           double *qr,double *qi,
           SMOOTHER smooth,double xopt,double cutoff,bool cont,double channelheight);
void rchannel(int N,
           const double *zr,const double *zi,
           const double *mr,
           int NE,
           const double *er,const double *ei,
           double *qr,double *qi,
           SMOOTHER smooth,double xopt,double cutoff,bool cont,double channelheight);

/*------------------------------------------------------------------------*/
void directInteractchannelpot(int N,
                    const double *zr,const double *zi,
                    const double *mr,const double *mi,
                    int NE,
                    const double *er,const double *ei,
                    double *pr,double *pi,
                    double *qr,double *qi,
                    const panel *panels,int Npanel,double channelheight,
                    SMOOTHER smooth,double xopt,double cutoff,
                    bool cont,double* timing,int printtime)
/* Direct evaluation of the potential from N pointmasses (mr,mi) (with
   mi possibly NULL) at the positions (zr,zi). The evaluation is done
   at the points (zr,zi) if (pr,pi) is not NULL and/or at the points
   (er,ei) if (qr,qi) is not NULL. The result is stored in the vectors
   (pr,pi) and (qr,qi), respectively, which must be allocated and
   cleared prior to call. Logarithmic potential for pot == 0, harmonic
   potential otherwise. The smoothing is done using smoother smooth
   with parameters (xopt,cutoff,cont). */
{
  // evaluation at (zr,zi)
  StartTime(printtime,timing);
  CudaStartTime(printtime,timing);
  if (pr != NULL) {
    // mass is real
    if (mi == NULL)
        rchannel2(N,zr,zi,mr,pr,pi,smooth,xopt,cutoff,cont,channelheight);
    // mass is complex
    else
        zchannel2(N,zr,zi,mr,mi,pr,pi,smooth,xopt,cutoff,cont,channelheight);
    directInteractChannelpotPanel(panels,pr,pi,zr,zi,Npanel,N,channelheight);
  }

  // evaluation at (er,ei)
  if (qr != NULL) {
    // mass is real
    if (mi == NULL)
        rchannel(N,zr,zi,mr,NE,er,ei,qr,qi,smooth,xopt,cutoff,cont,channelheight);
    // mass is complex
    else
        zchannel(N,zr,zi,mr,mi,NE,er,ei,qr,qi,smooth,xopt,cutoff,cont,channelheight);
    directInteractChannelpotPanel(panels,qr,qi,er,ei,Npanel,NE,channelheight);
  }
  StopTime(printtime,timing);
  CudaStopTime(printtime,timing);
  PrintTime("Direct summation System time: %f\n",timing,0,printtime);
  CudaPrintTime("Direct summation Cuda time",timing,1,printtime);

}

void zchannel2(int N,
           const double *restrict zr,const double *restrict zi,
           const double *restrict mr,const double *restrict mi,
           double *restrict pr,double *restrict pi,
           SMOOTHER smooth,double xopt,double cutoff,bool cont,double channelheight)
{
  double shape,scale = 1.0;
  const double sigma=M_PI/(2*channelheight);
  switch (smooth) {
  case DIRAC:
    for (int i = 0; i < N; i++) {
      {
        double im2=2*zi[i];
        const double ad2 = im2*im2;
        double P2,I2;
        im2*=sigma;
        if(ad2==0)
          P2=0;
        else
          P2=sigma/(2-2*cos(2*im2));
        I2=-2*sin(2*im2)*P2;

        pr[i]-=I2*mi[i];
        pi[i]-=I2*mr[i];
      }
      for (int j = i+1; j < N; j++) {
        double re = zr[i]-zr[j],im = zi[i]-zi[j],im2=zi[i]+zi[j];
        const double ad = re*re+im*im, ad2 = re*re+im2*im2;
        double P,P2,R1,R2,I1,I2;

        re*=sigma;
        im*=sigma;
        im2*=sigma;
        P=P2=sigma;

        double e2a=exp(2*re),em2a=exp(-2*re);

        if(ad==0)
          P=0;
        else
          P=sigma/(e2a+em2a-2*cos(2*im));

        if(ad2==0)
          P2=0;
        else
          P2=sigma/(e2a+em2a-2*cos(2*im2));

        R1=(e2a-em2a)*P;
        I1=-2*sin(2*im)*P;
        R2=(e2a-em2a)*P2;
        I2=-2*sin(2*im2)*P2;

        pr[i]-=(R1+R2)*mr[j]-(I1-I2)*mi[j];
        pi[i]-=(I1+I2)*mr[j]+(R1-R2)*mi[j];
        pr[j]-=-(R1+R2)*mr[i]+(I1+I2)*mi[i];
        pi[j]-=(-I1+I2)*mr[i]-(R1-R2)*mi[i];
      }
    }
    break;

  case RANKINE:
    cutoff = xopt*xopt;
    shape = 1.0/cutoff;
    for (int i = 0; i < N; i++) {
      {
        double im2=2*zi[i];
        const double ad2 = im2*im2;
        double P2,I2;
        im2*=sigma;
        if(ad2==0)
          P2=0;
        else {
          if (ad2 < cutoff)
            P2=shape;
          else
            P2=sigma/(2-2*cos(2*im2));
        }
        I2=-2*sin(2*im2)*P2;

        pr[i]-=I2*mi[i];
        pi[i]-=I2*mr[i];
      }
      for (int j = i+1; j < N; j++) {
        double re = zr[i]-zr[j],im = zi[i]-zi[j],im2=zi[i]+zi[j];
        const double ad = re*re+im*im, ad2 = re*re+im2*im2;
        double P,P2,R1,R2,I1,I2;

        re*=sigma;
        im*=sigma;
        im2*=sigma;
        P=P2=sigma;

        double e2a=exp(2*re),em2a=exp(-2*re);

        if(ad==0)
          P=0;
        else {
          if (ad < cutoff)
            P=shape;
          else
            P=sigma/(e2a+em2a-2*cos(2*im));
        }

        if(ad2==0)
          P2=0;
        else {
          if (ad2 < cutoff)
            P2=shape;
          else
            P2=sigma/(e2a+em2a-2*cos(2*im2));
        }

        R1=(e2a-em2a)*P;
        I1=-2*sin(2*im)*P;
        R2=(e2a-em2a)*P2;
        I2=-2*sin(2*im2)*P2;

        pr[i]-=(R1+R2)*mr[j]-(I1-I2)*mi[j];
        pi[i]-=(I1+I2)*mr[j]+(R1-R2)*mi[j];
        pr[j]-=-(R1+R2)*mr[i]+(I1+I2)*mi[i];
        pi[j]-=(-I1+I2)*mr[i]-(R1-R2)*mi[i];
      }
    }
    break;

  case SCULLY:
    cutoff = cutoff*cutoff;
    shape = xopt*xopt;
    if (cont)
      scale = 1.0+shape/cutoff;
    for (int i = 0; i < N; i++) {
      {
        double im2=2*zi[i];
        const double ad2 = im2*im2;
        double P2,I2;
        im2*=sigma;
        if(ad2==0)
          P2=0;
        else {
          if (ad2 < cutoff)
            P2=sigma/(2-2*cos(2*im2)+shape)*scale;
          else
            P2=sigma/(2-2*cos(2*im2));
        }
        I2=-2*sin(2*im2)*P2;

        pr[i]-=I2*mi[i];
        pi[i]-=I2*mr[i];
      }
      for (int j = i+1; j < N; j++) {
        double re = zr[i]-zr[j],im = zi[i]-zi[j],im2=zi[i]+zi[j];
        const double ad = re*re+im*im, ad2 = re*re+im2*im2;
        double P,P2,R1,R2,I1,I2;

        re*=sigma;
        im*=sigma;
        im2*=sigma;
        P=P2=sigma;

        double e2a=exp(2*re),em2a=exp(-2*re);

        if(ad==0)
          P=0;
        else {
          if (ad < cutoff)
            P=sigma/(e2a+em2a-2*cos(2*im)+shape)*scale;
          else
            P=sigma/(e2a+em2a-2*cos(2*im));
        }

        if(ad2==0)
          P2=0;
        else {
          if (ad2 < cutoff)
            P2=sigma/(e2a+em2a-2*cos(2*im2)+shape)*scale;
          else
            P2=sigma/(e2a+em2a-2*cos(2*im2));
        }

        R1=(e2a-em2a)*P;
        I1=-2*sin(2*im)*P;
        R2=(e2a-em2a)*P2;
        I2=-2*sin(2*im2)*P2;

        pr[i]-=(R1+R2)*mr[j]-(I1-I2)*mi[j];
        pi[i]-=(I1+I2)*mr[j]+(R1-R2)*mi[j];
        pr[j]-=-(R1+R2)*mr[i]+(I1+I2)*mi[i];
        pi[j]-=(-I1+I2)*mr[i]-(R1-R2)*mi[i];
      }
    }
    break;

  case OSEEN:
    cutoff = cutoff*cutoff;
    shape = 1.2564312086261696770/(xopt*xopt);
    if (cont)
      scale = -1.0/expm1(-shape*cutoff);
    for (int i = 0; i < N; i++) {
      {
        double im2=2*zi[i];
        const double ad2 = im2*im2;
        double P2,I2;
        im2*=sigma;
        if(ad2==0)
          P2=0;
        else {
          if (ad2 < cutoff)
            P2=-sigma*expm1(-shape*ad2)*scale/(2-2*cos(2*im2));
          else
            P2=sigma/(2-2*cos(2*im2));
        }
        I2=-2*sin(2*im2)*P2;

        pr[i]-=I2*mi[i];
        pi[i]-=I2*mr[i];
      }
      for (int j = i+1; j < N; j++) {
        double re = zr[i]-zr[j],im = zi[i]-zi[j],im2=zi[i]+zi[j];
        const double ad = re*re+im*im, ad2 = re*re+im2*im2;
        double P,P2,R1,R2,I1,I2;

        re*=sigma;
        im*=sigma;
        im2*=sigma;
        P=P2=sigma;

        double e2a=exp(2*re),em2a=exp(-2*re);

        if(ad==0)
          P=0;
        else {
          if (ad < cutoff)
            P=-sigma*expm1(-shape*ad)*scale/(e2a+em2a-2*cos(2*im));
          else
            P=sigma/(e2a+em2a-2*cos(2*im));
        }

        if(ad2==0)
          P2=0;
        else {
          if (ad2 < cutoff)
            P2=-sigma*expm1(-shape*ad2)*scale/(e2a+em2a-2*cos(2*im2));
          else
            P2=sigma/(e2a+em2a-2*cos(2*im2));
        }

        R1=(e2a-em2a)*P;
        I1=-2*sin(2*im)*P;
        R2=(e2a-em2a)*P2;
        I2=-2*sin(2*im2)*P2;

        pr[i]-=(R1+R2)*mr[j]-(I1-I2)*mi[j];
        pi[i]-=(I1+I2)*mr[j]+(R1-R2)*mi[j];
        pr[j]-=-(R1+R2)*mr[i]+(I1+I2)*mi[i];
        pi[j]-=(-I1+I2)*mr[i]-(R1-R2)*mi[i];
      }
    }
    break;
  }
}

void rchannel2(int N,
           const double *restrict zr,const double *restrict zi,
           const double *restrict mr,
           double *restrict pr,double *restrict pi,
           SMOOTHER smooth,double xopt,double cutoff,bool cont,double channelheight)
{
  double shape,scale = 1.0;
  const double sigma=M_PI/(2*channelheight);
  switch (smooth) {
  case DIRAC:
    for (int i = 0; i < N; i++) {
      {
        double im2=2*zi[i];
        const double ad2 = im2*im2;
        double P2,I2;
        im2*=sigma;
        if(ad2==0)
          P2=0;
        else
          P2=sigma/(2-2*cos(2*im2));
        I2=-2*sin(2*im2)*P2;

        pi[i]-=I2*mr[i];
      }
      for (int j = i+1; j < N; j++) {
        double re = zr[i]-zr[j],im = zi[i]-zi[j],im2=zi[i]+zi[j];
        const double ad = re*re+im*im, ad2 = re*re+im2*im2;
        double P,P2,R1,R2,I1,I2;

        re*=sigma;
        im*=sigma;
        im2*=sigma;
        P=P2=sigma;

        double e2a=exp(2*re),em2a=exp(-2*re);

        if(ad==0)
          P=0;
        else
          P=sigma/(e2a+em2a-2*cos(2*im));

        if(ad2==0)
          P2=0;
        else
          P2=sigma/(e2a+em2a-2*cos(2*im2));

        R1=(e2a-em2a)*P;
        I1=-2*sin(2*im)*P;
        R2=(e2a-em2a)*P2;
        I2=-2*sin(2*im2)*P2;

        pr[i]-=(R1+R2)*mr[j];
        pi[i]-=(I1+I2)*mr[j];
        pr[j]-=-(R1+R2)*mr[i];
        pi[j]-=(-I1+I2)*mr[i];
      }
    }
    break;

  case RANKINE:
    cutoff = xopt*xopt;
    shape = 1.0/cutoff;
    for (int i = 0; i < N; i++) {
      {
        double im2=2*zi[i];
        const double ad2 = im2*im2;
        double P2,I2;
        im2*=sigma;
        if(ad2==0)
          P2=0;
        else {
          if (ad2 < cutoff)
            P2=shape;
          else
            P2=sigma/(2-2*cos(2*im2));
        }
        I2=-2*sin(2*im2)*P2;

        pi[i]-=I2*mr[i];
      }
      for (int j = i+1; j < N; j++) {
        double re = zr[i]-zr[j],im = zi[i]-zi[j],im2=zi[i]+zi[j];
        const double ad = re*re+im*im, ad2 = re*re+im2*im2;
        double P,P2,R1,R2,I1,I2;

        re*=sigma;
        im*=sigma;
        im2*=sigma;
        P=P2=sigma;

        double e2a=exp(2*re),em2a=exp(-2*re);

        if(ad==0)
          P=0;
        else {
          if (ad < cutoff)
            P=shape;
          else
            P=sigma/(e2a+em2a-2*cos(2*im));
        }

        if(ad2==0)
          P2=0;
        else {
          if (ad2 < cutoff)
            P2=shape;
          else
            P2=sigma/(e2a+em2a-2*cos(2*im2));
        }

        R1=(e2a-em2a)*P;
        I1=-2*sin(2*im)*P;
        R2=(e2a-em2a)*P2;
        I2=-2*sin(2*im2)*P2;

        pr[i]-=(R1+R2)*mr[j];
        pi[i]-=(I1+I2)*mr[j];
        pr[j]-=-(R1+R2)*mr[i];
        pi[j]-=(-I1+I2)*mr[i];
      }
    }
    break;

  case SCULLY:
    cutoff = cutoff*cutoff;
    shape = xopt*xopt;
    if (cont)
      scale = 1.0+shape/cutoff;
    for (int i = 0; i < N; i++) {
      {
        double im2=2*zi[i];
        const double ad2 = im2*im2;
        double P2,I2;
        im2*=sigma;
        if(ad2==0)
          P2=0;
        else {
          if (ad2 < cutoff)
            P2=sigma/(2-2*cos(2*im2)+shape)*scale;
          else
            P2=sigma/(2-2*cos(2*im2));
        }
        I2=-2*sin(2*im2)*P2;

        pi[i]-=I2*mr[i];
      }
      for (int j = i+1; j < N; j++) {
        double re = zr[i]-zr[j],im = zi[i]-zi[j],im2=zi[i]+zi[j];
        const double ad = re*re+im*im, ad2 = re*re+im2*im2;
        double P,P2,R1,R2,I1,I2;

        re*=sigma;
        im*=sigma;
        im2*=sigma;
        P=P2=sigma;

        double e2a=exp(2*re),em2a=exp(-2*re);

        if(ad==0)
          P=0;
        else {
          if (ad < cutoff)
            P=sigma/(e2a+em2a-2*cos(2*im)+shape)*scale;
          else
            P=sigma/(e2a+em2a-2*cos(2*im));
        }

        if(ad2==0)
          P2=0;
        else {
          if (ad2 < cutoff)
            P2=sigma/(e2a+em2a-2*cos(2*im2)+shape)*scale;
          else
            P2=sigma/(e2a+em2a-2*cos(2*im2));
        }

        R1=(e2a-em2a)*P;
        I1=-2*sin(2*im)*P;
        R2=(e2a-em2a)*P2;
        I2=-2*sin(2*im2)*P2;

        pr[i]-=(R1+R2)*mr[j];
        pi[i]-=(I1+I2)*mr[j];
        pr[j]-=-(R1+R2)*mr[i];
        pi[j]-=(-I1+I2)*mr[i];
      }
    }
    break;

  case OSEEN:
    cutoff = cutoff*cutoff;
    shape = 1.2564312086261696770/(xopt*xopt);
    if (cont)
      scale = -1.0/expm1(-shape*cutoff);
    for (int i = 0; i < N; i++) {
      {
        double im2=2*zi[i];
        const double ad2 = im2*im2;
        double P2,I2;
        im2*=sigma;
        if(ad2==0)
          P2=0;
        else {
          if (ad2 < cutoff)
            P2=-sigma*expm1(-shape*ad2)*scale/(2-2*cos(2*im2));
          else
            P2=sigma/(2-2*cos(2*im2));
        }
        I2=-2*sin(2*im2)*P2;

        pi[i]-=I2*mr[i];
      }
      for (int j = i+1; j < N; j++) {
        double re = zr[i]-zr[j],im = zi[i]-zi[j],im2=zi[i]+zi[j];
        const double ad = re*re+im*im, ad2 = re*re+im2*im2;
        double P,P2,R1,R2,I1,I2;

        re*=sigma;
        im*=sigma;
        im2*=sigma;
        P=P2=sigma;

        double e2a=exp(2*re),em2a=exp(-2*re);

        if(ad==0)
          P=0;
        else {
          if (ad < cutoff)
            P=-sigma*expm1(-shape*ad)*scale/(e2a+em2a-2*cos(2*im));
          else
            P=sigma/(e2a+em2a-2*cos(2*im));
        }

        if(ad2==0)
          P2=0;
        else {
          if (ad2 < cutoff)
            P2=-sigma*expm1(-shape*ad2)*scale/(e2a+em2a-2*cos(2*im2));
          else
            P2=sigma/(e2a+em2a-2*cos(2*im2));
        }

        R1=(e2a-em2a)*P;
        I1=-2*sin(2*im)*P;
        R2=(e2a-em2a)*P2;
        I2=-2*sin(2*im2)*P2;

        pr[i]-=(R1+R2)*mr[j];
        pi[i]-=(I1+I2)*mr[j];
        pr[j]-=-(R1+R2)*mr[i];
        pi[j]-=(-I1+I2)*mr[i];
      }
    }
    break;
  }
}

void zchannel(int N,
           const double *restrict zr,const double *restrict zi,
           const double *restrict mr,const double *restrict mi,
           int NE,
           const double *restrict er,const double *restrict ei,
           double *restrict qr,double *restrict qi,
           SMOOTHER smooth,double xopt,double cutoff,bool cont,double channelheight)
{
  double shape,scale = 1.0;
  const double sigma=M_PI/(2*channelheight);
  switch (smooth) {
  case DIRAC:
    for (int j = 0; j < N; j++)
      for (int i = 0; i < NE; i++) {
        double re = er[i]-zr[j],im = ei[i]-zi[j],im2=ei[i]+zi[j];
        const double ad = re*re+im*im, ad2 = re*re+im2*im2;
        double P,P2,R1,R2,I1,I2;

        re*=sigma;
        im*=sigma;
        im2*=sigma;
        P=P2=sigma;

        double e2a=exp(2*re),em2a=exp(-2*re);

        if(ad==0)
          P=0;
        else
          P=sigma/(e2a+em2a-2*cos(2*im));

        if(ad2==0)
          P2=0;
        else
          P2=sigma/(e2a+em2a-2*cos(2*im2));

        R1=(e2a-em2a)*P;
        I1=-2*sin(2*im)*P;
        R2=(e2a-em2a)*P2;
        I2=-2*sin(2*im2)*P2;

        qr[i]-=(R1+R2)*mr[j]-(I1-I2)*mi[j];
        qi[i]-=(I1+I2)*mr[j]+(R1-R2)*mi[j];
      }
    break;

  case RANKINE:
    cutoff = xopt*xopt;
    shape = 1.0/cutoff;
    for (int j = 0; j < N; j++)
      for (int i = 0; i < NE; i++) {
        double re = er[i]-zr[j],im = ei[i]-zi[j],im2=ei[i]+zi[j];
        const double ad = re*re+im*im, ad2 = re*re+im2*im2;
        double P,P2,R1,R2,I1,I2;

        re*=sigma;
        im*=sigma;
        im2*=sigma;
        P=P2=sigma;

        double e2a=exp(2*re),em2a=exp(-2*re);

        if(ad==0)
          P=0;
        else {
          if (ad < cutoff)
            P=shape;
          else
            P=sigma/(e2a+em2a-2*cos(2*im));
        }

        if(ad2==0)
          P2=0;
        else {
          if (ad2 < cutoff)
            P2=shape;
          else
            P2=sigma/(e2a+em2a-2*cos(2*im2));
        }

        R1=(e2a-em2a)*P;
        I1=-2*sin(2*im)*P;
        R2=(e2a-em2a)*P2;
        I2=-2*sin(2*im2)*P2;

        qr[i]-=(R1+R2)*mr[j]-(I1-I2)*mi[j];
        qi[i]-=(I1+I2)*mr[j]+(R1-R2)*mi[j];
      }
    break;

  case SCULLY:
    cutoff = cutoff*cutoff;
    shape = xopt*xopt;
    if (cont)
      scale = 1.0+shape/cutoff;
    for (int j = 0; j < N; j++)
      for (int i = 0; i < NE; i++) {
        double re = er[i]-zr[j],im = ei[i]-zi[j],im2=ei[i]+zi[j];
        const double ad = re*re+im*im, ad2 = re*re+im2*im2;
        double P,P2,R1,R2,I1,I2;

        re*=sigma;
        im*=sigma;
        im2*=sigma;
        P=P2=sigma;

        double e2a=exp(2*re),em2a=exp(-2*re);

        if(ad==0)
          P=0;
        else {
          if (ad < cutoff)
            P=sigma/(e2a+em2a-2*cos(2*im)+shape)*scale;
          else
            P=sigma/(e2a+em2a-2*cos(2*im));
        }

        if(ad2==0)
          P2=0;
        else {
          if (ad2 < cutoff)
            P2=sigma/(e2a+em2a-2*cos(2*im2)+shape)*scale;
          else
            P2=sigma/(e2a+em2a-2*cos(2*im2));
        }

        R1=(e2a-em2a)*P;
        I1=-2*sin(2*im)*P;
        R2=(e2a-em2a)*P2;
        I2=-2*sin(2*im2)*P2;

        qr[i]-=(R1+R2)*mr[j]-(I1-I2)*mi[j];
        qi[i]-=(I1+I2)*mr[j]+(R1-R2)*mi[j];
      }
    break;

  case OSEEN:
    cutoff = cutoff*cutoff;
    shape = 1.2564312086261696770/(xopt*xopt);
    if (cont)
      scale = -1.0/expm1(-shape*cutoff);
    for (int j = 0; j < N; j++)
      for (int i = 0; i < NE; i++) {
        double re = er[i]-zr[j],im = ei[i]-zi[j],im2=ei[i]+zi[j];
        const double ad = re*re+im*im, ad2 = re*re+im2*im2;
        double P,P2,R1,R2,I1,I2;

        re*=sigma;
        im*=sigma;
        im2*=sigma;
        P=P2=sigma;

        double e2a=exp(2*re),em2a=exp(-2*re);

        if(ad==0)
          P=0;
        else {
          if (ad < cutoff)
            P=-sigma*expm1(-shape*ad)*scale/(e2a+em2a-2*cos(2*im));
          else
            P=sigma/(e2a+em2a-2*cos(2*im));
        }

        if(ad2==0)
          P2=0;
        else {
          if (ad2 < cutoff)
            P2=-sigma*expm1(-shape*ad2)*scale/(e2a+em2a-2*cos(2*im2));
          else
            P2=sigma/(e2a+em2a-2*cos(2*im2));
        }

        R1=(e2a-em2a)*P;
        I1=-2*sin(2*im)*P;
        R2=(e2a-em2a)*P2;
        I2=-2*sin(2*im2)*P2;

        qr[i]-=(R1+R2)*mr[j]-(I1-I2)*mi[j];
        qi[i]-=(I1+I2)*mr[j]+(R1-R2)*mi[j];
      }
    break;
  }
}

void rchannel(int N,
           const double *restrict zr,const double *restrict zi,
           const double *restrict mr,
           int NE,
           const double *restrict er,const double *restrict ei,
           double *restrict qr,double *restrict qi,
           SMOOTHER smooth,double xopt,double cutoff,bool cont,double channelheight)
{
  double shape,scale = 1.0;
  const double sigma=M_PI/(2*channelheight);
  switch (smooth) {
  case DIRAC:
    for (int j = 0; j < N; j++)
      for (int i = 0; i < NE; i++) {
        double re = er[i]-zr[j],im = ei[i]-zi[j],im2=ei[i]+zi[j];
        const double ad = re*re+im*im, ad2 = re*re+im2*im2;
        double P,P2,R1,R2,I1,I2;

        re*=sigma;
        im*=sigma;
        im2*=sigma;
        P=P2=sigma;

        double e2a=exp(2*re),em2a=exp(-2*re);

        if(ad==0)
          P=0;
        else
          P=sigma/(e2a+em2a-2*cos(2*im));

        if(ad2==0)
          P2=0;
        else
          P2=sigma/(e2a+em2a-2*cos(2*im2));

        R1=(e2a-em2a)*P;
        I1=-2*sin(2*im)*P;
        R2=(e2a-em2a)*P2;
        I2=-2*sin(2*im2)*P2;

        qr[i]-=(R1+R2)*mr[j];
        qi[i]-=(I1+I2)*mr[j];
      }
    break;

  case RANKINE:
    cutoff = xopt*xopt;
    shape = 1.0/cutoff;
    for (int j = 0; j < N; j++)
      for (int i = 0; i < NE; i++) {
        double re = er[i]-zr[j],im = ei[i]-zi[j],im2=ei[i]+zi[j];
        const double ad = re*re+im*im, ad2 = re*re+im2*im2;
        double P,P2,R1,R2,I1,I2;

        re*=sigma;
        im*=sigma;
        im2*=sigma;
        P=P2=sigma;

        double e2a=exp(2*re),em2a=exp(-2*re);

        if(ad==0)
          P=0;
        else {
          if (ad < cutoff)
            P=shape;
          else
            P=sigma/(e2a+em2a-2*cos(2*im));
        }

        if(ad2==0)
          P2=0;
        else {
          if (ad2 < cutoff)
            P2=shape;
          else
            P2=sigma/(e2a+em2a-2*cos(2*im2));
        }

        R1=(e2a-em2a)*P;
        I1=-2*sin(2*im)*P;
        R2=(e2a-em2a)*P2;
        I2=-2*sin(2*im2)*P2;

        qr[i]-=(R1+R2)*mr[j];
        qi[i]-=(I1+I2)*mr[j];
      }
    break;

  case SCULLY:
    cutoff = cutoff*cutoff;
    shape = xopt*xopt;
    if (cont)
      scale = 1.0+shape/cutoff;
    for (int j = 0; j < N; j++)
      for (int i = 0; i < NE; i++) {
        double re = er[i]-zr[j],im = ei[i]-zi[j],im2=ei[i]+zi[j];
        const double ad = re*re+im*im, ad2 = re*re+im2*im2;
        double P,P2,R1,R2,I1,I2;

        re*=sigma;
        im*=sigma;
        im2*=sigma;
        P=P2=sigma;

        double e2a=exp(2*re),em2a=exp(-2*re);

        if(ad==0)
          P=0;
        else {
          if (ad < cutoff)
            P=sigma/(e2a+em2a-2*cos(2*im)+shape)*scale;
          else
            P=sigma/(e2a+em2a-2*cos(2*im));
        }

        if(ad2==0)
          P2=0;
        else {
          if (ad2 < cutoff)
            P2=sigma/(e2a+em2a-2*cos(2*im2)+shape)*scale;
          else
            P2=sigma/(e2a+em2a-2*cos(2*im2));
        }

        R1=(e2a-em2a)*P;
        I1=-2*sin(2*im)*P;
        R2=(e2a-em2a)*P2;
        I2=-2*sin(2*im2)*P2;

        qr[i]-=(R1+R2)*mr[j];
        qi[i]-=(I1+I2)*mr[j];
      }
    break;

  case OSEEN:
    cutoff = cutoff*cutoff;
    shape = 1.2564312086261696770/(xopt*xopt);
    if (cont)
      scale = -1.0/expm1(-shape*cutoff);
    for (int j = 0; j < N; j++)
      for (int i = 0; i < NE; i++) {
        double re = er[i]-zr[j],im = ei[i]-zi[j],im2=ei[i]+zi[j];
        const double ad = re*re+im*im, ad2 = re*re+im2*im2;
        double P,P2,R1,R2,I1,I2;

        re*=sigma;
        im*=sigma;
        im2*=sigma;
        P=P2=sigma;

        double e2a=exp(2*re),em2a=exp(-2*re);

        if(ad==0)
          P=0;
        else {
          if (ad < cutoff)
            P=-sigma*expm1(-shape*ad)*scale/(e2a+em2a-2*cos(2*im));
          else
            P=sigma/(e2a+em2a-2*cos(2*im));
        }

        if(ad2==0)
          P2=0;
        else {
          if (ad2 < cutoff)
            P2=-sigma*expm1(-shape*ad2)*scale/(e2a+em2a-2*cos(2*im2));
          else
            P2=sigma/(e2a+em2a-2*cos(2*im2));
        }

        R1=(e2a-em2a)*P;
        I1=-2*sin(2*im)*P;
        R2=(e2a-em2a)*P2;
        I2=-2*sin(2*im2)*P2;

        qr[i]-=(R1+R2)*mr[j];
        qi[i]-=(I1+I2)*mr[j];
      }
    break;
  }
}
#endif /*CHANNELPOT*/
