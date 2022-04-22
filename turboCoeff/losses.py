# TURBOMACHINERY -- LOSSES LIBRARY
# AUTHOR: antonio pucciarelli 
#
# PROGRAM DESCRIPTION
#   CONTENT: functions that allow to compute the losses in a blade  
#  

# importing libraries 
import numpy as np 
import matplotlib.pyplot as plt
import scipy.integrate as integrate

def Dfactor(W1, W2, beta1, beta2, solidity):
    '''
    This function computes the diffusion factor parameter used for the computation of the aerodynamics losses in the system. 
        inputs:
            W1          -- inlet relative velocity; for the stator W1 == V1 
            W2          -- outlet relative velocity; for the stator W2 == V2
            beta1       -- inlet relative velocity flow angle; for the stator beta1 == alpha1
            beta2       -- outlet relative velocitu flow angle; for the stator beta2 == alpha2
            solidity    -- blade solidity 
    '''

    # check on the flow angles in order to have physical results
    if beta1 < 0:
        beta1 = - beta1 
        beta2 = - beta2 

    # tangential velocity computation
    W1t = W1 * np.sin(np.deg2rad(beta1))
    W2t = W2 * np.sin(np.deg2rad(beta2))

    # D computation 
    D = 1 - W2 / W1 + (W1t - W2t) / (2 * solidity * W1)

    return D 

def DeqFactor(W1=0, W2=0, beta1=0, beta2=0, r1=0, r2=0, Vt1=0, Vt2=0, Va1=0, solidity=0):
    '''
    This function computes the equivalent diffusion factor parameter used for the computation of the aerodynamics losses in the system. 
        inputs:
            W1          -- inlet relative velocity; for the stator W1 == V1 
            W2          -- outlet relative velocity; for the stator W2 == V2
            beta1       -- inlet relative velocity flow angle; for the stator beta1 == alpha1
            beta2       -- outlet relative velocitu flow angle; for the stator beta2 == alpha2
            r1          -- inlet section radius 
            r2          -- outlet section radius 
            Vt1         -- inlet absolute tangential speed 
            Vt2         -- outlet absolute tangential speed
            Va1         -- inlet absolute/relative axial speed 
            solidity    -- blade solidity 
    '''

    # Deq uses the Wmax / W1 value 
    if r1 == 0 or r2 == 0:     
        WmaxW1 = 1.12 + 0.61 * np.cos(np.deg2rad(beta1))**2 / solidity * (np.tan(np.deg2rad(beta1)) - np.tan(np.deg2rad(beta2)))
    else:
        WmaxW1 = 1.12 + 0.61 * np.cos(np.deg2rad(beta1))**2 / solidity * (r1 * Vt1 - r2 * Vt2) / (r1 * Va1)

    # Deq computation 
    Deq = WmaxW1 * W1 / W2 

    return Deq

def lossCoeff(W1=0, W2=0, beta1=0, beta2=0, solidity=1, D=0, r1=0, r2=0, Vt1=0, Vt2=0, Va1=0, kind='equivalent', plot=False, save=False, position='losses.pgf'):
    '''
    This function computes the loss coefficient in a blade section (2D) relating it to a diffusor section (3D). It ONLY applies for design incidence angles.
        inputs:
            W1          -- inlet relative velocity; for the stator W1 == V1 
            W2          -- outlet relative velocity; for the stator W2 == V2
            beta1       -- inlet relative velocity flow angle; for the stator beta1 == alpha1
            beta2       -- outlet relative velocitu flow angle; for the stator beta2 == alpha2
            r1          -- inlet section radius 
            r2          -- outlet section radius 
            Vt1         -- inlet absolute tangential speed 
            Vt2         -- outlet absolute tangential speed
            Va1         -- inlet absolute/relative axial speed 
            solidity    -- blade solidity
            D           -- diffusion factor 
            kind        -- selector for the type of loss coefficient vs diffusion coefficiente relaction 
                        -- values:
                            -- std -> D
                            -- equivalent -> Deq
                        -- chooses between
                            -- D -> 2D diffusion factor based on 3D correlation paremeters
                            -- Deq -> 2D diffusion factor based on 2D correlation parameters 
            plot        -- boolean value for the plotting of the function
            save        -- boolean value for the saving of the plot in vectorail format 
            position    -- path where the figure is saved
    '''

    def lossFunc(D, beta2, solidity, kind):
        '''
        This function computes the loss.
        '''

        if kind == 'std':
            loss = 0.0035 * (1 + 3.5 * D + 37 * D**4) * 2 * solidity / (np.cos(np.deg2rad(beta2)))
        elif kind == 'equivalent':
            loss = 0.004 * (1 + 3.1 * (D - 1)**2 + 0.4 * (D - 1)**8) * 2 * solidity / (np.cos(np.deg2rad(beta2)) * (W1 / W2)**2)

        return loss 

    if W1 != 0 and W2 != 0:
        # diffusion coefficient computation 
        if kind == 'std':
            D = Dfactor(W1, W2, beta1, beta2, solidity)
        elif kind == 'equivalent':
            D = DeqFactor(W1, W2, beta1, beta2, r1, r2, Vt1, Vt2, Va1, solidity)
        # loss computation 
        loss = lossFunc(D, beta2, solidity, kind)
    elif D != 0:
        # loss computation 
        loss = lossFunc(D, beta2, solidity, kind)

    if plot or save:
        if save:
            # setting matplotlib LaTeX export 
            import matplotlib
            matplotlib.use("pgf")
            matplotlib.rcParams.update({
                "pgf.texsystem": "pdflatex",
                'font.family': 'serif',
                'text.usetex': True,
                'pgf.rcfonts': False,
            })

        fig = plt.figure(figsize=(9,9))
        if kind == 'std':
            Dvec = np.linspace(0, 0.6, 300) 
            plt.plot(Dvec, lossFunc(Dvec, beta2, solidity, kind)*np.cos(np.deg2rad(beta2)) / (2*solidity), 'k', linewidth=2)
            plt.ylabel(r'$\omega^{{*}} \ \frac{{cos(\beta_2)}}{{2 \ \sigma}}$')
            plt.xlabel(r'$D^{{*}}$')
        elif kind == 'equivalent':
            Dvec = np.linspace(1, 2, 300)
            plt.plot(Dvec, lossFunc(Dvec, beta2, solidity, kind)*np.cos(np.deg2rad(beta2)) * (W1 / W2)**2 / (2*solidity), 'k', linewidth=2)
            plt.ylabel(r'$\omega^{{*}} \ \frac{{cos(\beta_2)}}{{2 \ \sigma}} \ \frac{{W_1}}{{W_2}}$')
            plt.xlabel(r'$D^{{*}}_{{eq}}$')

        try:
            if kind == 'std':
                plt.plot(D, loss*np.cos(np.deg2rad(beta2)) / (2*solidity), linestyle='', marker='o', color='g', markeredgewidth=1.5, markersize=8, markeredgecolor='k', label=r'$D^{{*}} = {0:.2f}, \omega^{{*}} \ \frac{{cos(\beta_2)}}{{2 \ \sigma}} = {1:.2f}$'.format(D, loss))
            elif kind == 'equivalent':
                plt.plot(D, loss*np.cos(np.deg2rad(beta2)) * (W1 / W2)**2 / (2*solidity), linestyle='', marker='o', color='g', markeredgewidth=1.5, markersize=8, markeredgecolor='k', label=r'$D^{{*}}_{{eq}} = {0:.2f}, \omega^{{*}} \ \frac{{cos(\beta_2)}}{{2 \ \sigma}} \ \frac{{W_1}}{{W_2}} = {1:.2f}$'.format(D, loss))
            plt.legend()
        except:
            pass 
        
        plt.grid(linestyle='--')
        
        if save:
            fig.savefig(position)
        else:  
            plt.show()

    if D != 0 and loss != 0:
        #print(loss)
        return loss, D

def lossHowell(beta1=0, beta2=0, solidity=0, pitch=0, bladeHeight=0, endWall=False):
    '''
    This function computes the section losses due to the 3D phenomena of the blade cascade.
        inputs:
            beta1       -- inlet relative velocity flow angle; for the stator beta1 == alpha1
            beta2       -- outlet relative velocitu flow angle; for the stator beta2 == alpha2
            solidity    -- blade solidity
            pitch       -- cascade pitch 
            bladeHeight -- height of the blade 
            endWall     -- boolean value for identify if the Howel relation is for the end wall losses or for the secondary flow losses 
    ''' 

    # beta angle check -> beta1 > 0 following ASME 
    if beta1 < 0:
        beta1 = - beta1 
        beta2 = - beta2 
    
    # computing average flow angle deflection 
    beta_ = np.arctan( ( np.tan(np.deg2rad(beta1)) + np.tan(np.deg2rad(beta2)) ) / 2)    

    # checking if the section in consideration is close to the wall 
    if endWall:
        # if close to the wall 
        # Cd computation 
        Cd = 0.02 * pitch / bladeHeight
    else:
        # if not close to the wall
        # Cl computation  
        Cl = 2 * np.cos(beta_) * (np.tan(np.deg2rad(beta1)) - np.tan(np.deg2rad(beta2))) / solidity

        # Cd computation with respect to Cl
        Cd = 0.18 * Cl**2 

    # losses computation
    loss = Cd * solidity * np.cos(np.deg2rad(beta1))**2 / np.cos(beta_)**3

    return loss

def machLosses(beta1, beta2, theta, i, W1, M1, Ttr, omegaStar, solidity, Ksh=0.1, R=287.06, gamma=1.4):
    '''
    This funciton computes the mach losses at each blade section.
        inputs:
            beta1       -- inlet relative flow angle 
            beta2       -- outlet relative flow angle 
            theta       -- camber angle
            i           -- design incidence angle 
            W1          -- inlet relative flow speed
            M1          -- inlet relative mach number 
            Ttr         -- inlet relative total temperature 
            omegaStart  -- profile loss without compressibility effect
            solidity    -- section solidity
            Ksh         -- shape factor
    '''
    
    # computing Rc/Rs related to the limit angles 
    Rc = - ( - 9 + (1 - (30 / beta1)**(0.48)) * theta / 8.2 )
    Rs = 10.3 + (2.92 - beta1 / 15.6) * theta / 8.2 

    # computing extremes incidence angles 
    iC = i - Rc / (1 + 0.5 * M1**3)
    iS = i + Rs / (1 + 0.5 * (Ksh * M1)**3)

    # computing mean incidence angle
    iM = iC + (iS - iC) * Rc / (Rc + Rs)

    # critical mach number computation
    WrStar = np.sqrt(2 * gamma / (gamma + 1) * R * Ttr)
    WmaxW1 = 1.12 + 0.61 * np.cos(np.deg2rad(beta1))**2 / solidity * (np.tan(np.deg2rad(beta1)) - np.tan(np.deg2rad(beta2)))
    Wmax   = WmaxW1 * W1
    # critical mach number for the section
    Mc     = M1 * WrStar / Wmax
        
    # computing minimum pressure loss due to compressibility effects 
    # this relation is related to losses without bow shocks for the design conditions
    if M1 < Mc:
        omegaM = omegaStar * (1 + (iM - i)**2 / Rs**2)
    else:
        omegaM = omegaStar * (1 + (iM - i)**2 / Rs**2) + Ksh * ((M1 / Mc - 1) * WrStar / W1)**2

    # computing final loss
    if i < iM:
        omega  = omegaM + omegaM * ((i - iM)/(iC - iM))**2
    else:
        omega  = omegaM + omegaM * ((i - iM)/(iS - iM))**2

    return omega

def tipLosses(r, clearance, VaInMid, VaOutMid, VtInMid, VtOutMid, rhoInMid, rhoOutMid, rhotMid, r1Mid, r2Mid, rHub, rTip, nBlade, chordMid, gammaMid, Nrow, mFlux, deltaPiso):
    '''
    This function computes the rotor clearance losses. The pressure loss distribution is linear from the hub (pressure loss == 0) to the rotor tip (max pressure loss).
        inputs:
            r           -- position of study
            clearance   -- blade clearance 
            VaInMid     -- reference section inlet absolute axial velocity
            VaOutMid    -- reference section outlet absolute axial velocity 
            VtInMid     -- reference section inlet absolute tangential velocity
            VtOutMid    -- reference section outlet absolute tangential velocity 
            rhoInMid    -- reference section inlet density
            rhoOutMid   -- reference section outlet density 
            rhotMid     -- reference section total density 
            r1Mid       -- reference section inlet radius 
            r2Mid       -- reference section outleet radius 
            rHub        -- rotor hub radius 
            rTip        -- rotor tip radius
            nBlade      -- # of rotor blades
            chordMid    -- reference section chord 
            gammaMid    -- reference section stagger angle
            Nrow        -- rotor blade # from 1 to # of turbomachine rotors + stators
            mFlux       -- mass flux 
    '''

    # clearance effect 
    tau = np.pi * clearance * ((r1Mid * rhoInMid * VaInMid) + (r2Mid * rhoOutMid * VaOutMid)) * (r2Mid * VtOutMid - r1Mid * VtInMid)

    # average pressure difference
    deltaP = tau / (nBlade * rTip * clearance * chordMid * np.cos(np.deg2rad(gammaMid)))

    # Uc is used for the reference mass flux 
    Uc = 0.816 * np.sqrt(2 * deltaP / rhotMid) / Nrow**0.2

    # reference mass flux computation
    mFluxC = rhotMid * Uc * nBlade * clearance * chordMid * np.cos(np.deg2rad(gammaMid))

    # total pressure loss alonge the blade due to tip clearance effects
    deltaPtot = deltaP * mFluxC / mFlux

    # following ASME directives -> the pressure loss distribution is linear and starts from 0 at hub 
    # the total pressure loss distribution is such that the integral of pressure loss is equal omegaTip 
    bladeHeight = rTip - rHub

    # computing deltaP at tip 
    # conserving pressure loss: deltaPtot = deltaPtip * bladeHeight / 2 
    deltaPtip = deltaPtot * 2 / bladeHeight
    
    # computing pressure loss at tip 
    lossTip = deltaPtip / deltaPiso

    # computing losses at section r
    # lossTip : bladeHeight = lossR : (r - rHub)
    lossR = lossTip * (r - rHub) / bladeHeight

    return lossR

def shockLosses(theta, tbc, c, gammaStagger, pitch, W1, M1, Ptr1, P1, gamma=1.4):
    '''
    This funciton computes the shock losses at each section due to the presence of a supersonic inlet flow.
        inputs:
            theta   -- blade section geoemtric angle 
            tbc     -- tb / c 
            c       -- section chord 

    '''

    # thetaU computation
    thetau = np.arctan(np.tan(np.deg2rad(theta/4)) + tbc) * 4

    # Ru computation 
    Ru = np.sin(thetau / 2) * 2 / c

    # bu computation 
    bu = np.tan(thetau/4) * c / 2 

    # computing psi 
    psi = 90 - np.rad2deg(thetau)/2 - gammaStagger

    # computing phi 
    phi = np.arctan(pitch * np.cos(np.deg2rad(psi)) / (pitch * np.sin(np.deg2rad(psi)) + Ru))

    # computing sound speed 
    a1 = W1 / M1

    # Prandtl-Meyer expansion 
    PRfunc = lambda W: np.sqrt(W**2/a1**2 - 1) / W

    # generating W vector and computing Ws that will be used later on for the normal shock mach # 
    Wmax = 2.0 * a1
    Wvec = np.linspace(W1, Wmax, 100)

    # allocating storing vector 
    integralVec = np.zeros(len(Wvec))
    
    # computing integral and errors
    for ii in range(len(Wvec)):
        #print(integrate.quad(PRfunc, W1, Wvec[ii])[0])
        integralVec[ii] = np.abs(integrate.quad(PRfunc, W1, Wvec[ii])[0] - phi)

    # getting Ws for the computation of Min 
    WsPos = np.argmin(integralVec)
    Ws = Wvec[WsPos]

    # computing normal mach number that will be used for the computation of the losses 
    Min = np.sqrt(W1/a1 * Ws/a1)

    # total pressure after the normal shock wave 
    Ptr2 = Ptr1 * ((gamma+1)*Min**2/((gamma-1)*Min**2+2))**(gamma/(gamma-1)) * ((gamma+1)/((2*gamma*Min**2)-(gamma-1)))**(1/(gamma-1))

    # computing loss 
    loss = (Ptr1 - Ptr2) / (Ptr1 - P1)

    return loss