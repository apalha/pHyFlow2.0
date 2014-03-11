/* paneldirect.h */
/* A. Goude 2010-01-01 */

#ifndef NORMALPANELCODE
#define NORMALPANELCODE(lambda) sqrdist=rotxdist*rotxdist+rotydist*rotydist;\
xterm=1-(lambda)*rotxdist/sqrdist;\
yterm=(lambda)*rotydist/sqrdist;\
logxterm=0.5*log(xterm*xterm+yterm*yterm);\
logyterm=atan2(yterm, xterm);\
xtmp0=rotxdist/(lambda);\
ytmp0=rotydist/(lambda);\
bx=xtmp0*logxterm-ytmp0*logyterm+1;\
by=xtmp0*logyterm+ytmp0*logxterm;\
ax=bx-logxterm;\
ay=by-logyterm;
#endif

#undef NORMALPANELEND
#ifdef PURECOMPLEXPANELS
    #define NORMALPANELEND xtmp=EVALPANEL.istrength1*ax-EVALPANEL.istrength2*bx;\
    ytmp=EVALPANEL.istrength1*ay-EVALPANEL.istrength2*by;\
            REALDESTINATION+=(ytmp*EVALPANEL.costheta-xtmp*EVALPANEL.sintheta);\
            IMAGDESTINATION-=(xtmp*EVALPANEL.costheta+ytmp*EVALPANEL.sintheta);
#else
    #ifdef COMPLEXPANELS
    #define NORMALPANELEND xtmp=EVALPANEL.strength1*ay+EVALPANEL.istrength1*ax-EVALPANEL.strength2*by-EVALPANEL.istrength2*bx;\
    ytmp=-EVALPANEL.strength1*ax+EVALPANEL.istrength1*ay+EVALPANEL.strength2*bx-EVALPANEL.istrength2*by;\
            REALDESTINATION+=(ytmp*EVALPANEL.costheta-xtmp*EVALPANEL.sintheta);\
            IMAGDESTINATION-=(xtmp*EVALPANEL.costheta+ytmp*EVALPANEL.sintheta);
    #else
    #define NORMALPANELEND xtmp=EVALPANEL.strength1*ay-EVALPANEL.strength2*by;\
    ytmp=-EVALPANEL.strength1*ax+EVALPANEL.strength2*bx;\
            REALDESTINATION+=(ytmp*EVALPANEL.costheta-xtmp*EVALPANEL.sintheta);\
            IMAGDESTINATION-=(xtmp*EVALPANEL.costheta+ytmp*EVALPANEL.sintheta);
    #endif
#endif 
#ifndef INV2PI
    #define INV2PI 0.159154943091895
#endif
#ifdef PURECOMPLEXPANELS
    #undef REALPANELS
    #ifndef COMPLEXPANELS
        #define COMPLEXPANELS
    #endif
#else
    #define REALPANELS
#endif

//mexPrintf("panel: %fx%f to %fx%f strength %f+%fi to %f+%fi smoother=%d cutoff=%f\n",EVALPANEL.x1,EVALPANEL.y1,EVALPANEL.x2,EVALPANEL.y2,EVALPANEL.strength1,EVALPANEL.istrength1,EVALPANEL.strength2,EVALPANEL.istrength2,EVALPANEL.smoother,EVALPANEL.cutoff);
if(EVALPANEL.smoother==DIRACSMOOTHER || EVALPANEL.cutoff==0.0) {
    for (int j = PANELLOOPSTART; j < PANELLOOPEND; j++) {
        xdist=REALEVAL-EVALPANEL.x1;
        ydist=IMAGEVAL-EVALPANEL.y1;
        rotxdist=xdist*EVALPANEL.costheta+ydist*EVALPANEL.sintheta;
        rotydist=ydist*EVALPANEL.costheta-xdist*EVALPANEL.sintheta;
        if(fabs(rotydist)<EVALPANEL.lambda*1e-12)
            rotydist=(2*EVALPANEL.side-1)*EVALPANEL.lambda*1e-12;
        NORMALPANELCODE(EVALPANEL.lambda);
        NORMALPANELEND
    }
}
else {
    for (int j = PANELLOOPSTART; j < PANELLOOPEND; j++) { //for speed issues, use a clean inner loop in case of no smoothing
        xdist=REALEVAL-EVALPANEL.x1;
        ydist=IMAGEVAL-EVALPANEL.y1;
        rotxdist=xdist*EVALPANEL.costheta+ydist*EVALPANEL.sintheta;
        rotydist=ydist*EVALPANEL.costheta-xdist*EVALPANEL.sintheta;
        xtmp=0;
        ytmp=0;
        if(rotydist>=EVALPANEL.cutoff||rotydist<=-EVALPANEL.cutoff||rotxdist<=-EVALPANEL.cutoff||rotxdist-EVALPANEL.lambda>=EVALPANEL.cutoff||
                (rotxdist<0&&(rotxdist*rotxdist+rotydist*rotydist)>=EVALPANEL.cutoff*EVALPANEL.cutoff)||
                (rotxdist-EVALPANEL.lambda>0&&((rotxdist-EVALPANEL.lambda)*(rotxdist-EVALPANEL.lambda)+rotydist*rotydist)>=EVALPANEL.cutoff*EVALPANEL.cutoff)) {
            NORMALPANELCODE(EVALPANEL.lambda);
            NORMALPANELEND;
            
        }
        else {
            double invcutoff=1.0/EVALPANEL.cutoff;
            double tmp=sqrt(EVALPANEL.cutoff*EVALPANEL.cutoff-rotydist*rotydist);
            double lowcut=rotxdist-tmp;
            double highcut=rotxdist+tmp;
            
            if(lowcut<=0)
                lowcut=0;
            if(highcut>EVALPANEL.lambda)
                highcut=EVALPANEL.lambda;
            
            #ifdef REALPANELS
            double strengthslope=(EVALPANEL.strength2-EVALPANEL.strength1)/EVALPANEL.lambda;
            double panelstrengthlowcut=EVALPANEL.strength1+strengthslope*lowcut;
            double panelstrengthhighcut=EVALPANEL.strength1+strengthslope*highcut;
//            mexPrintf("strengthslope=%f  panelstrengthlowcut=%f panelstrengthhighcut=%f\n",strengthslope,panelstrengthlowcut,panelstrengthhighcut);
            #endif
            #ifdef COMPLEXPANELS
            double istrengthslope=(EVALPANEL.istrength2-EVALPANEL.istrength1)/EVALPANEL.lambda;
            double istrengthlowcut=EVALPANEL.istrength1+istrengthslope*lowcut;
            double istrengthhighcut=EVALPANEL.istrength1+istrengthslope*highcut;
//            mexPrintf("istrengthslope=%f  ipanelstrengthlowcut=%f ipanelstrengthhighcut=%f\n",istrengthslope,istrengthlowcut,istrengthhighcut);
            #endif
//            mexPrintf("lowcut=%f highcut=%f lambda=%f\n",lowcut,highcut,EVALPANEL.lambda);
            if(lowcut>0){
                NORMALPANELCODE(lowcut);
                #ifdef PURECOMPLEXPANELS
                xtmp=EVALPANEL.istrength1*ax-istrengthlowcut*bx;
                ytmp=EVALPANEL.istrength1*ay-istrengthlowcut*by;
//                mexPrintf("pure complex panel mode xtmp=%f ytmp=%f\n",xtmp,ytmp);
                #else
                #ifdef COMPLEXPANELS
                xtmp=EVALPANEL.strength1*ay-panelstrengthlowcut*by+EVALPANEL.istrength1*ax-istrengthlowcut*bx;
                ytmp=-EVALPANEL.strength1*ax+panelstrengthlowcut*bx+EVALPANEL.istrength1*ay-istrengthlowcut*by;
//                mexPrintf("complex panel mode xtmp=%f ytmp=%f\n",xtmp,ytmp);
                #else
                xtmp=EVALPANEL.strength1*ay-panelstrengthlowcut*by;
                ytmp=-EVALPANEL.strength1*ax+panelstrengthlowcut*bx;
//                mexPrintf("real panel mode xtmp=%f ytmp=%f\n",xtmp,ytmp);
                #endif
                #endif
            }
            double startintegral=-rotxdist+lowcut;
            double stopintegral=-rotxdist+highcut;
            double realbcoeffs,realkcoeffs,imagacoeffs;
            if(EVALPANEL.smoother==HATSMOOTHER) {
                double tmp1=startintegral*startintegral+rotydist*rotydist;
                double stmp1=sqrt(tmp1);
                double tmp2=stopintegral*stopintegral+rotydist*rotydist;
                double stmp2=sqrt(tmp2);
                realbcoeffs=(startintegral*startintegral*startintegral-stopintegral*stopintegral*stopintegral)+
                        0.5*(stmp2*stopintegral*stopintegral*stopintegral-stmp1*startintegral*startintegral*startintegral)*invcutoff;
                realkcoeffs=1.5*(startintegral*startintegral-stopintegral*stopintegral)+
                        8.0/12*(stmp2*tmp2-stmp1*tmp1)*invcutoff;
                imagacoeffs=0;
                if(rotydist!=0) {
                    if(stopintegral+stmp2==0 || startintegral+stmp1==0)
                        tmp=0;
                    else
                        tmp=log((stopintegral+stmp2)/(startintegral+stmp1))*rotydist*rotydist*rotydist;
                    imagacoeffs+=(rotydist*(startintegral*(3*EVALPANEL.cutoff-stmp1)-stopintegral*(3*EVALPANEL.cutoff-stmp2))+tmp)*invcutoff;
                    realbcoeffs+=(0.25*rotydist*rotydist*(stopintegral*stmp2-startintegral*stmp1) - tmp*rotydist/4)*invcutoff;
                }
            }
            else if(EVALPANEL.smoother==RANKINESMOOTHER) {
                
                realkcoeffs=0.5*(startintegral-stopintegral)*(startintegral+stopintegral);
                realbcoeffs=1.0/3*(startintegral*startintegral*startintegral-stopintegral*stopintegral*stopintegral);
                imagacoeffs=(startintegral-stopintegral)*rotydist;
                /*
                 * Vnew=1i.*exp(-1i.*theta(j)).*((m - n).*(3*(m+n).*(a+b.*1i.*z) + 2.*b.*(m.^2 + n.^2 +m.*n) + 6.*a.*1i.*z))./(12*pi*cutoff^2);
                 */
            }
//            mexPrintf("xtmp=%f ytmp=%f\n",xtmp,ytmp);
            #ifdef PURECOMPLEXPANELS
            double itmp=istrengthlowcut+istrengthslope*(rotxdist-lowcut);
            xtmp+=1*invcutoff*invcutoff*(realbcoeffs*istrengthslope+realkcoeffs*itmp);
            ytmp+=1*invcutoff*invcutoff*(imagacoeffs*itmp+rotydist*realkcoeffs*istrengthslope);
            #else
            #ifdef COMPLEXPANELS
            tmp=panelstrengthlowcut+strengthslope*(rotxdist-lowcut);
            double itmp=istrengthlowcut+istrengthslope*(rotxdist-lowcut);
            /*
             * mexPrintf("xtmp=%f ytmp=%f\n",xtmp,ytmp);
             */
            xtmp+=1*invcutoff*invcutoff*(imagacoeffs*tmp+rotydist*realkcoeffs*strengthslope+realbcoeffs*istrengthslope+realkcoeffs*itmp);
            ytmp-=1*invcutoff*invcutoff*(realbcoeffs*strengthslope+realkcoeffs*tmp-imagacoeffs*itmp-rotydist*realkcoeffs*istrengthslope);
//            mexPrintf("tmp=%f itmp=%f,imagacoeffs=%f realkcoeffs=%f realbcoeffs=%f\n",tmp,itmp,imagacoeffs,realkcoeffs,realbcoeffs);
            #else
            tmp=panelstrengthlowcut+strengthslope*(rotxdist-lowcut);
            xtmp+=1*invcutoff*invcutoff*(imagacoeffs*tmp+rotydist*realkcoeffs*strengthslope);
            ytmp-=1*invcutoff*invcutoff*(realbcoeffs*strengthslope+realkcoeffs*tmp);
            #endif
            #endif
//            mexPrintf("v2.1 xtmp=%f ytmp=%f\n",xtmp,ytmp);
            if(highcut<EVALPANEL.lambda) {
                rotxdist-=highcut;
                tmp=EVALPANEL.lambda-highcut;
                NORMALPANELCODE(tmp);
//                mexPrintf("v2.2 xtmp=%f ytmp=%f\n",xtmp,ytmp);
                #ifdef PURECOMPLEXPANELS
                xtmp+=istrengthhighcut*ax-EVALPANEL.istrength2*bx;
                ytmp+=istrengthhighcut*ay-EVALPANEL.istrength2*by;
                #else
                #ifdef COMPLEXPANELS
//                        mexPrintf("v2.3 xtmp=%f ytmp=%f\n",xtmp,ytmp);
                xtmp+=panelstrengthhighcut*ay-EVALPANEL.strength2*by+istrengthhighcut*ax-EVALPANEL.istrength2*bx;
                ytmp+=-panelstrengthhighcut*ax+EVALPANEL.strength2*bx+istrengthhighcut*ay-EVALPANEL.istrength2*by;
//                mexPrintf("panelstrengthhighcut=%f ay=%f panelstrength2=%f by=%f sourcestrengthhighcut=%f ax=%f sourcestrength2=%f bx=%f\n",panelstrengthhighcut,ay,EVALPANEL.strength2,by,istrengthhighcut,ax,EVALPANEL.istrength2,bx);
//                mexPrintf("added %f %f, xtmp=%f ytmp=%f\n",panelstrengthhighcut*ay-EVALPANEL.strength2*by+istrengthhighcut*ax-EVALPANEL.istrength2*bx,-panelstrengthhighcut*ax+EVALPANEL.strength2*bx+istrengthhighcut*ay-EVALPANEL.istrength2*by,xtmp,ytmp);
                #else
                xtmp+=panelstrengthhighcut*ay-EVALPANEL.strength2*by;
                ytmp+=-panelstrengthhighcut*ax+EVALPANEL.strength2*bx;
                #endif
                #endif
//                        mexPrintf("v3 xtmp=%f ytmp=%f\n",xtmp,ytmp);
            }
            REALDESTINATION+=(ytmp*EVALPANEL.costheta-xtmp*EVALPANEL.sintheta);
            IMAGDESTINATION-=(xtmp*EVALPANEL.costheta+ytmp*EVALPANEL.sintheta);
        }
    }
}
