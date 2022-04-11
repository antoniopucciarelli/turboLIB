import numpy as np 

def Vt2Func(L, rMeanInlet, rMeanOutlet, omega, Vt1):
    '''
    Computation of the exit tangential velocity with respect to:
        L           -- total enthalphy jump 
        rMeanInlet  -- mean inlet radius
        rMeanOutlet -- mean outlet radius 
        omega       -- angular speed 
        Vt1         -- tangential speed at the inlet 
    '''

    # L = U2 * Vt2 - U1 * Vt1 
    Vt2 = ( L + rMeanInlet * omega * Vt1 ) / (rMeanOutlet * omega)

    return Vt2  

def bladeHeightFinder(mFlux, rMean, Va, rho):
    '''
    This function computes the blade height relative to a circular crown inlet.
        inputs:
            mFlux   -- mass flux 
            rMean   -- mean radius 
            Va      -- axial flow velocity 
            rho     -- fluid density 
    '''

    b = mFlux / (2 * np.pi * rMean * Va * rho)

    rHub = (2 * rMean * b - b**2) / (2 * b)

    rTip = rHub + b 

    return b, rHub, rTip 



