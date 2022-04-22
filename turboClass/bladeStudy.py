# TURBOMACHINERY -- BLADE GENERATION 
# AUTHOR: antonio pucciarelli 
#
# PROGRAM DESCRIPTION
#   CONTENT: functions that allow to compute the airfoil sections in a blade 
#  

# importing libraries
import numpy as np
import matplotlib.pyplot as plt  

def bladeStudy(rIn, rOut, omega, rMean, VaMean, VtMeanIn, VtMeanOut, Leu, Tt0, T0, Pt0, P0, eta=1, printout=False, gamma=1.4, R=287.06):
    '''
    This function computes the behaviour of the blade at a radius r with respect to the guideline properties described by the meanline.
        The model used is the FREE VORTEX model with these assumptions: 
            -- inlet enthalpy distribution -> constant along r 
            -- inlet axial speed distribusiton -> constant along r
            -- the mean radius for the inlet and the outlet of the blade is the same 

        inputs:
            rIn         -- blade section inlet radial position
            rOut        -- blade secton outlet radial position 
            omega       -- rotation speed 
            rMean       -- stage mean radius 
            VaMean      -- axial flow speed at mean radius 
            VtMeanIn    -- tangential speed at mean radius inlet 
            VtMeanOut   -- tangential speed at mean radius outlet 
            Leu         -- euler work of the section 
            Tt0         -- inlet total temperature 
            T0          -- inlet static temperature 
            Pt0         -- inlet total pressure 
            P0          -- inlet static pressure 
            eta         -- stage efficiency
            printout    -- boolean value for the printing of the computed quantities 
            gamma       -- specific heat ratio
            R           -- gas constant
    '''

    # cP computation
    cP = gamma / (gamma - 1) * R

    # rotor inlet 
    # Va0 computation 
    Va0 = VaMean

    # rotation speed 
    U0 = omega * rIn

    # phi computation 
    if U0 != 0:
        phi = Va0 / U0
    else: 
        phi = 0

    # psi computation 
    if Leu != 0:
        psi = Leu / U0**2 
    else: 
        psi = 0 

    # lam computation 
    lam = psi * 2 
    
    # Vt0 computation
    Vt0 = VtMeanIn * rMean / rIn 

    # V0 magnitude computation 
    V0 = np.sqrt(Va0**2 + Vt0**2)

    # Wa0 computation 
    Wa0 = Va0 

    # Wt0 computation 
    Wt0 = Vt0 - U0 

    # W0 magnitude computation 
    W0 = np.sqrt(Wa0**2 + Wt0**2)

    # alpha0 angle computation
    alpha0 = np.rad2deg(np.arctan(Vt0/Va0)) 

    # beta0 angle computation 
    beta0 = np.rad2deg(np.arctan(Wt0/Wa0))

    # temperature computation 
    T0 = Tt0 - V0**2 / (2 * cP)

    # speed of sound computation
    a0 = np.sqrt(gamma * R * T0)

    # mach number computation -- absolute
    M0 = V0 / a0

    # mach number computation -- relative 
    Mr0 = W0 / a0 

    # pressure computation
    P0 = Pt0 * (T0/Tt0)**(gamma/(gamma-1))

    # density computation 
    rho0 = P0 / (R * T0)

    # total density computation
    rhot0 = Pt0 / (R * Tt0)
    
    # rotor outlet computation
    # Va1 computation 
    Va1 = VaMean

    # rotation speed 
    U1 = omega * rOut
    
    # Vt1 computation
    Vt1 = VtMeanOut * rMean / rOut 

    # V1 magnitude computation 
    V1 = np.sqrt(Va1**2 + Vt1**2)

    # Wa1 computation
    Wa1 = Va1 

    # Wt1 computation
    Wt1 = Vt1 - U1

    # W1 computation
    W1 = np.sqrt(Wa1**2 + Wt1**2)

    # alpha1 angle computation
    alpha1 = np.rad2deg(np.arctan(Vt1/Va1)) 

    # beta1 angle computation
    beta1 = np.rad2deg(np.arctan(Wt1/Wa1))

    # total temperature computation
    Tt1 = Leu / cP + Tt0

    # ideal temperature computation if the process is completely 
    # isentropic without losses but the work produced is related 
    # to a process that takes into account losses in the stage  
    T1 = Tt1 - V1**2 / (2 * cP)

    # T1 isoentropic computation 
    #   this correction activates only is eta != 1
    T1 = T0 + eta * (T1 - T0)

    # mach number computation 
    a1 = np.sqrt(gamma * R * T1)

    # mach number computation -- absolute
    M1 = V1 / a1

    # mach number computation -- relative 
    Mr1 = W1 / a1 

    # pressure computation
    P1 = P0 * (T1/T0)**(gamma/(gamma-1))

    # total pressure computation
    Pt1 = P1 * (Tt1/T1)**(gamma/(gamma-1))

    # density computation 
    rho1 = P1 / (R * T1)

    # total density computation
    rhot1 = Pt1 / (R * Tt1)

    # reaction degree computation 
    if Leu != 0:
        rD = cP * (T1 - T0) / Leu 
    else: 
        rD = 0

    if printout:
        printLength = 54
        nAdim = int((printLength - len(' ADIMENSIONAL PARAMETERS '))/2)
        print('*' * nAdim + ' ADIMENSIONAL PARAMETERS  ' + '*' * nAdim)
        print('-- rD     = {0:>8.2f}        -- psi    = {1:>8.2f}'.format(rD, psi))
        print('-- lambda = {0:>8.2f}        -- phi    = {1:>8.2f}'.format(lam, phi))
        nWork = int((printLength - len(' WORK '))/2)
        print('*' * nWork + ' WORK ' + '*' * nWork)
        print('-- Leu    = {0:>8.2f} J/kg'.format(Leu))
        nRotation = int((printLength - len(' ROTATION '))/2)
        print('*' * nRotation + ' ROTATION ' + '*' * nRotation)
        print('-- U0  = {0:>8.2f} m/s    -- U1  = {1:>8.2f} rad/s'.format(U0, U1))
        print('*' * printLength)
        nDefinitions = int((printLength - len(' DEFINITIONS '))/2)
        print('\n' + '*' * nDefinitions + ' DEFINITIONS  ' + '*' * nDefinitions)
        print('-- 0 => rotor inlet         -- 1 => rotor outlet')
        nTotTemp = int((printLength - len(' TOTAL TEMPERATURE '))/2)
        print('*' * nTotTemp + ' TOTAL TEMPERATURE  ' + '*' * nTotTemp)
        print('-- Tt0    = {0:>8.2f} K      -- Tt1    = {1:>8.2f} K'.format(Tt0, Tt1))
        nAbsKin = int((printLength - len('KINETICS -- ABSOLUTE '))/2)
        print('*' * nAbsKin + ' KINETICS -- ABSOLUTE ' + '*' * nAbsKin)
        print('-- 0                        -- 1')
        print('-- alpha0 = {0:>8.2f} deg    -- alpha1 = {1:>8.2f} deg'.format(alpha0, alpha1))
        print('-- Va0    = {0:>8.2f} m/s    -- Va1    = {1:>8.2f} m/s'.format(Va0, Va1))
        print('-- Vt0    = {0:>8.2f} m/s    -- Vt1    = {1:>8.2f} m/s'.format(Vt0, Vt1))
        print('-- V0     = {0:>8.2f} m/s    -- V1     = {1:>8.2f} m/s'.format(V0, V1))
        nRelKin = int((printLength - len('KINETICS -- RELATIVE '))/2)
        print('*' * nRelKin + ' KINETICS -- RELATIVE ' + '*' * nRelKin)
        print('-- 0                        -- 1')
        print('-- U0     = {0:>8.2f} m/s    -- U1     = {1:>8.2f} m/s'.format(U0, U1))
        print('-- beta0  = {0:>8.2f} deg    -- beta1  = {1:>8.2f} deg'.format(beta0, beta1))
        print('-- Wa0    = {0:>8.2f} m/s    -- Wa1    = {1:>8.2f} m/s'.format(Wa0, Wa1))
        print('-- Wt0    = {0:>8.2f} m/s    -- Wt1    = {1:>8.2f} m/s'.format(Wt0, Wt1))
        print('-- W0     = {0:>8.2f} m/s    -- W1     = {1:>8.2f} m/s'.format(W0, W1))
        nThermo = int((printLength - len(' THERMODYNAMICS '))/2)
        print('*' * nThermo + ' THERMODYNAMICS ' + '*' * nThermo)
        print('-- 0                        -- 1 ')
        print('-- T0     = {0:>8.2f} K      -- T1     = {1:>8.2f} K    '.format(T0, T1))
        print('-- Tt0    = {0:>8.2f} K      -- Tt1    = {1:>8.2f} K    '.format(Tt0, Tt1))
        print('-- a0     = {0:>8.2f} m/s    -- a1     = {1:>8.2f} m/s  '.format(a0, a1))
        print('-- M0     = {0:>8.2f}        -- M1     = {1:>8.2f}      '.format(M0, M1))
        print('-- Mr0    = {0:>8.2f}        -- Mr1    = {1:>8.2f}      '.format(Mr0, Mr1))
        print('-- P0     = {0:>8.2f} bar    -- P1     = {1:>8.2f} bar  '.format(P0/1e+5, P1/1e+5))
        print('-- Pt0    = {0:>8.2f} bar    -- Pt1    = {1:>8.2f} bar  '.format(Pt0/1e+5, Pt1/1e+5))
        print('-- rho0   = {0:>8.2f} kg/m3  -- rho1   = {1:>8.2f} kg/m3'.format(rho0, rho1))
        print('-- rhot0  = {0:>8.2f} kg/m3  -- rhot1  = {1:>8.2f} kg/m3'.format(rhot0, rhot1))
        nBlade = int((printLength - len(' BLADE DIMENSIONS '))/2)
        print('*' * nBlade + ' BLADE DIMENSIONS ' + '*' * nBlade)
        print('-- 0                        -- 1')
        print('-- rIn    = {0:>8.2f} cm     -- rOut   = {1:>8.2f} cm'.format(rIn*1e+2, rOut*1e+2))
        print('-- rMean0 = {0:>8.2f} cm     -- rMean1 = {0:>8.2f} cm'.format(rMean*1e+2))
        print('*' * printLength) 

    # setting up output vectors 
    adimVec     = [rD, phi, psi, lam]
    rotationVec = [U0, U1]
    abs0Vec     = [Va0, Vt0, V0]
    rel0Vec     = [Wa0, Wt0, W0]
    abs1Vec     = [Va1, Vt1, V1]
    rel1Vec     = [Wa1, Wt1, W1]
    angleVec    = [alpha0, alpha1, beta0, beta1]
    thermo0     = [T0, P0, rho0, Tt0, Pt0, rhot0, M0, Mr0]
    thermo1     = [T1, P1, rho1, Tt1, Pt1, rhot1, M1, Mr1]

    return adimVec, rotationVec, abs0Vec, rel0Vec, abs1Vec, rel1Vec, angleVec, thermo0, thermo1

def optimalPlot(beta1, beta2, error=1e-4):
    '''
    This function computes the optimal incidence for a blade section given flow deflections.
    '''

    # importing libraries 
    from turboCoeff import lieblein
    
    # vector generation 
    thetaVec = np.linspace(20, 70, 10)
    solidityVec = np.linspace(0.8, 1.1, 10)
    tbcVec = np.linspace(0.08, 0.1, 10)

    # figure generation
    fig0, ax0 = plt.subplots(nrows=1, ncols=3, figsize=(8,8))
    fig1, ax1 = plt.subplots(nrows=1, ncols=3, figsize=(8,8))

    # incidence and deviatoric angle 3 fold loop
    for _,theta in enumerate(thetaVec):
        for _,solidity in enumerate(solidityVec):
            for _,tbc in enumerate(tbcVec):
                try: 
                    # incidence and deviatoric angle computation 
                    Ksh = 1 # shape factor 
                    Kti = lieblein.KtiFunc(tbc, solidity, plot=False)
                    i0 = lieblein.i0Func(beta1, solidity, plot=False)
                    n = lieblein.nFunc(beta1, solidity, plot=False)
                    Ktdelta = lieblein.KtdeltaFunc(tbc, solidity, plot=False)
                    delta0 = lieblein.delta0Func(beta1, solidity, plot=False)
                    m = lieblein.mFunc(beta1, solidity, plot=False)

                    # incidence angle 
                    i =  Kti * Ksh * i0 + n * theta
                    # deviation angle 
                    delta = Ktdelta * Ksh * delta0 + m * theta

                    # i = beta1 - k1        -> k1 = beta1 - i     
                    # delta = beta2 - k2    -> k2 = beta2 - delta 
                    # theta = k1 - k2       -> theta = beta1 - i - beta2 + delta 
                    deltaTheta = ((beta1 - i - beta2 + delta) - theta) / theta
                    ax0[0].scatter(theta, deltaTheta*theta)
                    ax0[1].scatter(tbc, deltaTheta*theta)
                    ax0[2].scatter(solidity, deltaTheta*theta)
                    
                    if deltaTheta / theta < error: 
                        ax1[0].scatter(solidity, tbc)
                        ax1[1].scatter(solidity, theta)
                        ax1[2].scatter(tbc, theta)
                except: 
                    pass 

    # ax setup
    ax0[0].grid(linestyle='--')
    ax0[0].set_xlabel(r'$\theta$')
    ax0[0].set_ylabel(r'$\Delta \theta_{{(i^{{*}}, \delta^{{*}})}}$')
    # ax setup
    ax0[1].grid(linestyle='--')
    ax0[1].set_xlabel(r'$\frac{t_b}{c}$')
    ax0[1].set_ylabel(r'$\Delta \theta_{{(i^{{*}}, \delta^{{*}})}}$')
    # ax setup
    ax0[2].grid(linestyle='--')
    ax0[2].set_xlabel(r'$\sigma$')
    ax0[2].set_ylabel(r'$\Delta \theta_{{(i^{{*}}, \delta^{{*}})}}$')

    # ax setup
    fig1.suptitle(r'$\frac{{\Delta \theta}}{{\theta}} < {0:f}$'.format(error))
    ax1[0].grid(linestyle='--')
    ax1[0].set_xlabel(r'$\sigma$')
    ax1[0].set_ylabel(r'$\frac{t_b}{c}$')
    ax1[1].grid(linestyle='--')
    ax1[1].set_xlabel(r'$\sigma$')
    ax1[1].set_ylabel(r'$\theta$')
    ax1[2].grid(linestyle='--')
    ax1[2].set_xlabel(r'$\frac{t_b}{c}$')
    ax1[2].set_ylabel(r'$\theta$')
    
    fig0.tight_layout()
    fig1.tight_layout()
    plt.show()

def iFunc(beta1, tbc, solidity, theta):
    '''
    This function computes the incidence angle following Leiblein model given as inputs:
        inputs:
            beta1       -- inlet relative flow angle 
            tbc         -- tb/c thickness / chord 
            solidity    -- blade solidity
            theta       -- geometric blade deflection 
    '''
    
    # importing libraries 
    from turboCoeff import lieblein

    Ksh = 1                                             # shape correction factor 
    Kti = lieblein.KtiFunc(tbc, solidity, plot=False)   # i thickenss correction factor
    i0  = lieblein.i0Func(beta1, solidity, plot=False)  # i0 computation 
    n   = lieblein.nFunc(beta1, solidity, plot=False)   # n computation 

    return Kti * Ksh * i0 + n * theta

def deltaFunc(beta1, tbc, solidity, theta):
    '''
    This function computes the deflection angle following Leiblein model given as inputs:
        inputs:
            beta1       -- inlet relative flow angle 
            tbc         -- tb/c thickness / chord 
            solidity    -- blade solidity
            theta       -- geometric blade deflection 
    '''

    # importing libraries 
    from turboCoeff import lieblein

    Ksh     = 1                                                 # shape correction factor 
    Ktdelta = lieblein.KtdeltaFunc(tbc, solidity, plot=False)   # delta thickness correction factor
    delta0  = lieblein.delta0Func(beta1, solidity, plot=False)  # delta 0 computation 
    m       = lieblein.mFunc(beta1, solidity, plot=False)       # m computation 

    return Ktdelta * Ksh * delta0 + m * theta

def alphaFunc(ac, solidity, theta, tbc):
    '''
    This function computes the optimal design angle of attach followin the Leiblein approach.
        inputs:
            ac          -- max(x_camber) / chord 
            solidity    -- chord / pitch 
            theta       -- geometrical deflection angle 
            tbc         -- thichness / chord
    '''
    
    # importing libraries
    from turboCoeff import lieblein

    # Ksh set to 0.1
    Ksh = 0.1

    # exponential value 
    e = 0.65 - 0.002 * theta 

    alpha = (3.6 * Ksh * lieblein.KtiFunc(tbc, plot=False) + 0.3532 * theta * ac**0.25) * solidity**e

    return alpha 

def thetaFunc(beta1, beta2, i, delta):
    '''
    This function computes the blade geometric deflection angle given as inputs:
        inputs:
            beta1       -- inlet design flow angle 
            beta2       -- outlet design flow angle 
            i           -- incidence angle 
            delta       -- deviation angle
    '''

    theta = beta1 - i - beta2 + delta 
    
    return theta

def optimalAngles(beta1, beta2, solidity, tbc=0.1, printout=False):
    '''
    This function computes the optimal incidence angle and deviation angle for a blade section given flow deflections.
        inputs: 
            beta1       -- inlet flow angle 
            beta2       -- outlet flow angle 
            solidity    -- section solidity 
            tbc         -- section thickness / chord
            printout    -- boolean value for the printout of the results
    '''

    # importing libraries 
    from scipy import optimize

    def deltaTheta(x):
        '''
        This function computes the error made by the choosing theta, tbc, solidity as blade properties.
            direct inputs:
                x/theta     -- choosed geometrical deflection angle
            global function inputs:
                solidity    -- chord / pitch 
                tbc         -- thickness / chord 
                beta1       -- inlet flow angle 
                beta2       -- outlet flow angle 
        '''

        # variable definition
        theta = x[0] # theta allocation

        # incidence angle computation 
        i = iFunc(beta1, tbc, solidity, theta)

        # deflection angle computation 
        delta = deltaFunc(beta1, tbc, solidity, theta)

        # error computation
        error = np.abs(thetaFunc(beta1, beta2, i, delta) - theta)

        return error 

    # bounds definition 
    epsilon = beta1 - beta2
    bounds  = [(epsilon, None)]
    
    # setting up initial study point 
    x0 = (epsilon*1.05)

    # minimizing system 
    res = optimize.minimize(deltaTheta, x0, bounds=bounds, tol=5e-7)

    # data allocation
    theta = res.x[0]
    i     = iFunc(beta1, tbc, solidity, theta)
    delta = deltaFunc(beta1, tbc, solidity, theta)

    # results printout
    if printout:
        starDim = 34
        bladeDim = np.int16((starDim - len(' ANGLE LOOP '))/2)
        print('*' * bladeDim + ' ANGLE LOOP ' + '*' * bladeDim)
        print('-- converged = ', res.success)
        print('-- error     = {0:>11.2e} deg'.format(deltaTheta(res.x)))
        print('.' * starDim)
        print('-- beta1     = {0:>8.3f}    deg'.format(beta1))
        print('-- beta2     = {0:>8.3f}    deg'.format(beta2))
        print('-- epsilon   = {0:>8.3f}    deg'.format(epsilon))
        print('.' * starDim)
        print('-- i         = {0:>8.3f}    deg'.format(i))
        print('-- delta     = {0:>8.3f}    deg'.format(delta))
        print('-- theta     = {0:>8.3f}    deg'.format(theta))
        print('-- solidity  = {0:>8.3f}'.format(solidity))
        print('-- tb/c      = {0:>8.3f}'.format(tbc))

    return i, delta, theta

def optimalBladeNumber(W1, W2, beta1, beta2, rMean, bladeHeight, r1=0, r2=0, Vt1=0, Vt2=0, Va1=0, bladeInterval=[25,50], ARvec=[1.2, 1.5, 1.8, 2], kind='equivalent', save=False):
    '''
    This function plots the profile losses at the mean line with respect to different aspect ratio and blade number.
        inputs:
            beta1           -- inlet relative flow angle @ mean line 
            beta2           -- outlet relative flow angle @ mean line 
            bladeHeight     -- inlet blade heigth
            bladeInterval   -- array that stores the interval of # of blade study 
            ARvec           -- aspect ratio (bladeHeight / chord) vector of study
            kind            -- enables Leiblein std/equivalent losses
                            -- std => profile losses + secondary flow losses
                            -- equivalent => profile losses 
            save            -- boolean valure for the saving of the plot
    ''' 

    # importing libraries
    from turboCoeff import losses

    # markersize
    markersize=4

    # setting up blade number vector 
    bladeVec = np.arange(bladeInterval[0], bladeInterval[1]+1, 1)

    # figure generation 
    fig, ax = plt.subplots(nrows=2, ncols=1)

    # computing profile losses for different blade number and aspect ratio
    for AR in ARvec:
        # allocating vectors
        lossVec = np.zeros(len(bladeVec))
        Dvec    = np.zeros(len(bladeVec))
        # chord computation
        chord = bladeHeight * AR 
    
        for ii,nBlades in enumerate(bladeVec):
            # pitch computation
            pitch = 2 * np.pi * rMean / nBlades
            # solidity computation
            solidity = chord / pitch 
            
            # loss and diffusion factor computation
            lossVec[ii], Dvec[ii] = losses.profileLosses(W1=W1, W2=W2, beta1=beta1, beta2=beta2, solidity=solidity, D=0, r1=r1, r2=r2, Vt1=Vt1, Vt2=Vt2, Va1=Va1, kind=kind)

        # plotting results
        ax[0].plot(bladeVec, Dvec, linestyle='-', marker='o', markeredgewidth=1.5, markersize=markersize, markeredgecolor='black', label=r'$AR = {0:.2f}$'.format(AR))
        ax[1].plot(bladeVec, lossVec, linestyle='-', marker='o', markeredgewidth=1.5, markersize=markersize, markeredgecolor='black', label=r'$AR = {0:.2f}$'.format(AR))
    
    ax[0].legend(loc='upper left', bbox_to_anchor=[1,1])
    ax[0].set_xlabel('# of blades')
    ax[0].grid(linestyle='--')
    if kind == 'std':
        ax[0].set_ylabel(r'$D$')
    else: 
        ax[0].set_ylabel(r'$D_{eq }$')

    #ax[1].legend()
    ax[1].set_ylabel(r'$\bar{\omega }$')
    ax[1].spines["bottom"].set_position(("axes", 1))
    ax[1].spines["bottom"].set_position(("axes", 0))
    ax[1].xaxis.set_ticks_position('top')
    ax[1].grid(linestyle='--')

    fig.tight_layout()
    plt.show()
