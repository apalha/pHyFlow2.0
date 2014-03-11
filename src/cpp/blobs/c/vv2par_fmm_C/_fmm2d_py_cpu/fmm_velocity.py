import numpy
from _fmm2d_py import fmm2d_py

def fmm_velocity(xBlob,yBlob,wBlob,sigma,k=2.0,xTarget=None,yTarget=None,Ndirect=35,tol=1.0e-6,cutoff=None):

    # convert blob parameters into the parameters of the fmm solver
    alpha = 1.2564312086261696770
    
    ksigmasqr = k*sigma*sigma
    xopt = numpy.sqrt(alpha*ksigmasqr)
    if cutoff==None: cutoff = 5.0*xopt

    if xTarget==None:
        xTarget = xBlob.copy()
        yTarget = yBlob.copy()

    
    vy_fmm_cpu,vx_fmm_cpu = fmm2d_py(xBlob,yBlob,wBlob,xTarget,yTarget,Ndirect,xopt,cutoff,tol)
    
    return -vx_fmm_cpu/(2.0*numpy.pi),-vy_fmm_cpu/(2.0*numpy.pi)
