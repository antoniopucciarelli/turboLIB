# TURBOMACHINERY -- ENGINEERING DIMENSION AND COEFFICIENT FUNCTIONS LIBRARY
# AUTHOR: antonio pucciarelli 
#
# PROGRAM DESCRIPTION
#   TURBOMACHINERY: compressor
#   CONTENT: engineering coefficient functions for the computation of turbomachinery performace
#  

# importing libraries
import numpy as np 
import matplotlib.pyplot as plt 
from scipy import interpolate

def efficiency(phi=0, rD=0, plot=False, save=False, position='efficiency.pgf'):
    '''
    This function describes the adimensional parameters for the turbomachinery design.
        inputs: 
            phi         -- flow coefficient
            rD          -- reaction degree
            plot        -- boolean value for the plotting of the Leiblein efficiency chart
            save        -- boolean value for the saving of the charts in vectorial format
            position    -- path where to save the vectorial image (if save == True)
        output:
            eta         -- efficiency with respect Leiblein chart
    '''

    # origin vector 
    originVec = np.array([[0.48, 0.5], [0.5, 0.5],[0.567, 0.5],[0.627, 0.5], [0.75, 0.5], [0.87, 0.5]])

    # radius vector
    radiusVec = np.array([0, 0.2, 0.337, 0.44, 0.609, 0.76])

    # efficiency vector 
    etaVec = np.array([0.926, 0.92, 0.91, 0.90, 0.88, 0.86])

    # radius vs Xorigin interpolation 
    origin2radius = interpolate.interp1d(originVec[:,0], radiusVec, kind='linear')

    # radius vs efficiency intepolation
    radius2eta = interpolate.interp1d(radiusVec, etaVec, kind='linear')

    # computing efficiency from input data 
    if rD != 0 and phi != 0:
        # setting up vector dimension 
        dim = 1000 

        # computing possible point distance from circle center
        deltaX = np.linspace(phi - 0.48, 0, dim)
        deltaY = rD - 0.5 
        
        # computing center position 
        XoriginVec = 0.48 + np.linspace(0, phi - 0.48, dim)
        
        # computing all possible radius 
        deltaRadius = np.sqrt(deltaX**2 + deltaY**2)
        
        # computing radius related to eta circles
        etaRadius = origin2radius(XoriginVec)
        
        # finding minimum error position
        r = np.abs(etaRadius - deltaRadius)
        ii = np.where(r == r.min())

        # computing origin, circle radius and efficiency
        origin = np.array([XoriginVec[ii].item(), 0.5])
        radius = etaRadius[ii].item()
        eta = radius2eta(radius)

    # plotting efficiency
    # plotting vectors
    # theta vector 
    thetaVec = np.linspace(0, 2*np.pi, 1000)
    # label vector 
    labelVec = ['', r'$\eta = 0.92$', r'$\eta = 0.91$', r'$\eta = 0.90$', r'$\eta = 0.88$', r'$\eta = 0.86$']
    # color vector 
    colorVec = ['k', 'b', 'r', 'm', 'c', 'y']
    
    if plot:

        plt.figure(figsize=(9,9))
        plt.title('Lieblein efficiency\nAXIAL COMPRESSOR')

        # plotting efficiency taken from charts 
        plt.plot(0.48, 0.5, 'ok', label=r'$\eta = 0.926$')
        for ii,r in enumerate(radiusVec):
            plt.plot(originVec[ii,0] + np.cos(thetaVec) * r, originVec[ii,1] + np.sin(thetaVec) * r, color=colorVec[ii], label=labelVec[ii])

        # plotting efficiency point computed
        if rD != 0 and phi != 0:
            plt.plot(origin[0] + np.cos(thetaVec) * radius, origin[1] + np.sin(thetaVec) * radius, 'g--', label=r'$\eta = {0:.4f}$'.format(eta))
            plt.plot(phi, rD, linestyle='', marker='s', markeredgewidth=2, markeredgecolor='k', color='g', markersize=8, label=r'$[\phi = {0:.2f}, \chi = {1:.2f}]$'.format(phi, rD))
        
        plt.xlabel(r'$\phi$')
        plt.ylabel(r'$\chi$')
        plt.xlim(0,1.5)
        plt.ylim(0,1.5)
        plt.xticks(np.arange(0,1.6,0.1))
        plt.yticks(np.arange(0,1.6,0.1))
        plt.legend(loc='upper left',bbox_to_anchor=[0,1])  
        plt.grid(linestyle='--')
        plt.tight_layout()
        plt.show()
    
    elif save:
        # setting matplotlib LaTeX export 
        import matplotlib
        matplotlib.use("pgf")
        matplotlib.rcParams.update({
            "pgf.texsystem": "pdflatex",
            'font.family': 'serif',
            'text.usetex': True,
            'pgf.rcfonts': False,
        })

        plt.figure(figsize=(9,9))
        plt.title('Lieblein efficiency\nAXIAL COMPRESSOR')

        # plotting efficiency taken from charts 
        plt.plot(0.48, 0.5, 'ok', label=r'$\eta = 0.926$')
        for ii,r in enumerate(radiusVec):
            plt.plot(originVec[ii,0] + np.cos(thetaVec) * r, originVec[ii,1] + np.sin(thetaVec) * r, color=colorVec[ii], label=labelVec[ii])

        # plotting efficiency point computed
        if rD != 0 and phi != 0:
            plt.plot(origin[0] + np.cos(thetaVec) * radius, origin[1] + np.sin(thetaVec) * radius, 'g--', label=r'$\eta = {0:.4f}$'.format(eta))
            plt.plot(phi, rD, linestyle='', marker='s', markeredgewidth=2, markeredgecolor='k', color='g', markersize=8, label=r'$[\phi = {0:.2f}, \chi = {1:.2f}]$'.format(phi, rD))
        
        plt.xlabel(r'$\phi$')
        plt.ylabel(r'$\chi$')
        plt.xlim(0,1.5)
        plt.ylim(0,1.5)
        plt.xticks(np.arange(0,1.6,0.1))
        plt.yticks(np.arange(0,1.6,0.1))
        plt.legend(loc='upper left',bbox_to_anchor=[0,1])  
        plt.grid(linestyle='--')
        plt.tight_layout()

        # figure saving 
        plt.savefig(position, bbox_inches='tight')

    if rD != 0 and phi != 0:
        return eta

def stagePerf(phi=0, psi=0, perc=1, rD=0.5, phiVec=np.linspace(0,1.5,1000), plot=True, save=False, position='stagePerf.pgf'):
    '''
    This function plots the phi and psi parameter with respect to ASME axial compressor design constraints
        The main constraints for the study of the adimensional parameters are:
            W2/W1 > 0.7     -- to avoid flow separation 
            psiLim = 1      -- for having enough margin from the surge condition 
                            -- ASME -> figure 10.6 psiLim determines line slopes that has to be < 1 in order to avoid surge 
            beta1 < 70 deg  -- for avoiding reduction in cascade performances 
            curve symmetric with respect to reaction degree = 0.5
        inputs:
            phi         -- flow coefficient for the stage 
            psi         -- work coefficient for the stage
            perc        -- percentage of psi with respect to the psi limit 
                        -- this variable activates if:
                            perc != 0 
                            psi == 0 
            rD          -- stage reaction degree 
            phiVec      -- phi vector for the psiLim representation 
            plot        -- boolean value for plotting the chart 
            save        -- boolean value for plotting the chart in vectorial format 
            position    -- path where the pgf file should be saved 
        output:
            psi         -- if perc != 0 and phi != 0 -> computing psi given phi as input 

    '''
    
    # function generation eqn. from ASME 10.28 - 10.30
    Rcap = 0.5 + np.abs(rD - 0.5)
    psiLim = lambda phi: 6 * Rcap / 17 + 0.85 * (0.5/Rcap)**1.18 * phi**(2 + 0.1/Rcap)
    psiBeta = lambda phi: 5 * phi - 2 * rD

    def phiFinder(phiVec, psiTarget, perc=perc, error=5e-3):
        '''
        This function computes the phi value for a defined psi value.
            inputs: 
                phiVec      -- array that stores phi study points
                psiTarget   -- psiLim(phi) = psiTarget
                perc        -- sagety margin with respect to the minimum velocity ratio 
                               such that there is not huge recirculation of flow in the compressor 
        '''

        if perc != 0:
            for _,phi in enumerate(phiVec):
                if np.abs((perc * psiLim(phi) - psiTarget) / (perc * psiLim(phi))) < error:
                    return phi 

    # values computation
    psiLIM = psiLim(phiVec)             # work coefficient limit vs flow coefficient 
    psiBETA = psiBeta(phiVec)           # work coefficient vs flow coefficient for beta = 70deg 
    phi1 = phiFinder(phiVec, 1, perc=1) # phi value such that psiLim == 1

    # computing psi with respect to a percentage and phi 
    if perc != 0 and phi != 0 and psi == 0:      
        # computing psi 
        psi = perc * psiLim(phi)

    # computing phi such that psi == psiTarget
    if phi == 0 and psi != 0:
        phi = phiFinder(phiVec, psi)

    if plot:
        # plotting
        plt.figure(figsize=(8,8))
        plt.plot(phiVec, psiLIM, 'k', label=r'$\psi_{Lim}$')
        plt.plot(phiVec, psiBETA, 'r', label=r'$\beta > 70^{\circ}$')
        plt.plot([phi1, phi1] , [0, np.max(psiLIM)], 'b', label=r'$\phi_{Lim}$')
        if psi != 0 and phi !=0:
            plt.plot(phi, psi, linestyle='', marker='o', markersize=7, markeredgewidth=1, markeredgecolor='k', color='g', label=r'$[\phi = {0:.3f}, \psi = {1:.3f}]$'.format(phi, psi))
        plt.xlim(0,phi1)
        plt.ylim(0,1)
        plt.xlabel(r'$\phi$')
        plt.ylabel(r'$\psi$')
        plt.xticks(np.arange(0,phi1+0.1,0.1))
        plt.yticks(np.arange(0,1.1,0.1))
        if perc == 0:
            plt.title(r'$\frac{W_2}{W_1} > 0.7$' + '    ' + r'$\chi = {0:.2f}$'.format(rD))
        else:
            plt.title(r'$\frac{W_2}{W_1} > 0.7$' + '    ' + r'$\chi = {0:.2f}$'.format(rD) + '    ' + r'$Percentage = {0:d} \%$'.format(int(perc*100)))
        plt.grid(linestyle='--')
        plt.legend(loc='upper left', bbox_to_anchor=(0, 1))
        plt.tight_layout()
        plt.show()

    elif save:
        # setting matplotlib LaTeX export 
        import matplotlib
        matplotlib.use("pgf")
        matplotlib.rcParams.update({
            "pgf.texsystem": "pdflatex",
            'font.family': 'serif',
            'text.usetex': True,
            'pgf.rcfonts': False,
        })
        # plotting 
        plt.figure(figsize=(8,8))
        plt.plot(phiVec, psiLIM, 'k', label=r'$\psi_{Lim}$')
        plt.plot(phiVec, psiBETA, 'r', label=r'$\beta > 70^{\circ}$')
        plt.plot([phi1, phi1] , [0, np.max(psiLIM)], 'b', label=r'$\phi_{Lim}$')
        if psi != 0 and phi !=0:
            plt.plot(phi, psi, 'og', label=r'$[\phi = {0:.3f}, \psi = {1:.3f}]$'.format(phi, psi))
        plt.xlim(0,np.max(phiVec))
        plt.ylim(0,1)
        plt.xlabel(r'$\phi$')
        plt.ylabel(r'$\psi$')
        plt.xticks(np.arange(0,np.max(phiVec)+0.1,0.1))
        plt.yticks(np.arange(0,1.1,0.1))
        if perc == 0:
            plt.title(r'$\frac{W_2}{W_1} > 0.7$' + '    ' + r'$\chi = {0:.2f}$'.format(rD))
        else:
            plt.title(r'$\frac{W_2}{W_1} > 0.7$' + '    ' + r'$\chi = {0:.2f}$'.format(rD) + '    ' + r'$Percentage = {0:d} \%$'.format(int(perc*100)))
        plt.grid(linestyle='--')
        plt.legend(loc='upper left', bbox_to_anchor=(0, 1))
        plt.tight_layout()
        # figure saving 
        plt.savefig(position, bbox_inches='tight')

    if phi != 0 and psi !=0:
        return phi, psi 

def stageProperties(rD, psi, rMean, mFlux, Tt0, Pt0, betaP, T1real=False, printout=False, R=287.06, gamma=1.4):
    '''
    This function allows to compute the properties of a stage following the meanline initial design procedure.
        procedural steps:
            - from rD and psi -> find phi related to W2/W1 > 0.7 condition
            - from phi and rD -> compute eta from Leiblein chart
            - compute Lis and L using [betaP, Tt0, eta]
            - compute rotation parameters 
            - compute Va0, Vt0 and Vt1 
            - compute kinetic and thermodynamic properties at rotor inlet 
            - compute kinetic and thermodynamic properties at rotor outlet
                - efficiency conditions check 
            - compute kinetic and thermodynamic properties at stator outlet

        inputs:
            rD          -- reaction degree at the mean line 
            psi         -- work coefficient at the mean line 
            rMean       -- mean line radius 
            mFlux       -- mass flux 
            Tt0         -- total temperature 
            Pt0         -- tottal pressure 
            betaP       -- pressure ratio 
            T1real      -- boolean value for computation of the real or ideal rotor outlet temperature startig from real work L = Lis / eta
            printout    -- boolean value for the print of the results 
            R           -- gas constant
            gamma       -- specific heat ratio
        
        outputs:
            check return at the bottom
    '''

    # importing libraries  
    from turboCoeff import coeff

    # air properties allocation 
    cP = gamma / (gamma - 1) * R # specific heat ratio @ P cost  [J/kg K]

    # adimensional parameters computation 
    # -- from performance charts 
    phi, psi = stagePerf(psi=psi, rD=rD, plot=False, perc=0.97)
    eta = efficiency(phi=phi, rD=rD, plot=False)

    ######################## WORK ########################
    # ideal compression work 
    Lis = coeff.L_is(Tin=Tt0, beta=betaP, gamma=gamma, kind='compressor')
    # real compression work 
    L = Lis / eta
    # total temperature computation
    Tt1 = L / cP + Tt0

    ######################## ROTATION ######################## 
    # mean section revolution speed from work coefficient
    Umean = np.sqrt(L / psi) 
    # angular speed 
    omega = Umean / rMean
    # rpm 
    n = omega * 60 / (2 * np.pi)

    ######### MAIN VELOCITIES COMPUTATION FROM WORK RESULTS #########
    # Va0, Vt0 computation 
    # axial inlet velocity from flow coefficient
    Va0 = phi * Umean
    # tangential inlet velocity from reaction degree and work coefficient
    lam = 2 * psi 
    Vt0 = (1 - rD - lam/4) * Umean 
    # tangential rotor exit velocity computation
    Vt1 = (1 - rD + lam/4) * Umean 

    ######### ROTOR INLET QUANTITIES ######### 
    # KINETICS 
    # relative speed computation 
    Wa0 = Va0
    Wt0 = Vt0 - Umean
    # veolcity magnitude computation
    V0 = np.sqrt(Va0**2 + Vt0**2)
    W0 = np.sqrt(Wa0**2 + Wt0**2)
    # aerodynamic angles computation 
    alpha0 = np.rad2deg(np.arctan(Vt0/Va0))
    beta0 = np.rad2deg(np.arctan(Wt0/Wa0))
    # THERMODYNAMICS
    # static temperature computation
    T0 = Tt0 - V0**2 / (2*cP)
    # speed of sound computation
    a0 = np.sqrt(gamma * R * T0)
    # mach computation 
    M0 = V0 / a0
    Mr0 = W0 / a0
    # static pressure computation
    P0 = Pt0 / (1 + (gamma - 1)/2 * M0**2)**(gamma/(gamma-1))
    # density computation 
    rhot0 = Pt0 / (R * Tt0)
    rho0 = rhot0 / (1 + (gamma - 1)/2 * M0**2)**(1/(gamma-1))

    ######### ROTOR OUTLET/STATOR INLET QUANTITIES #########  
    # outlet quantities 
    Va1 = Va0 
    Wa1 = Va1 
    Wt1 = Vt1 - Umean
    # velocity magnitude computation
    V1 = np.sqrt(Va1**2 + Vt1**2)
    W1 = np.sqrt(Wa1**2 + Wt1**2)
    # flow angle computation
    alpha1 = np.rad2deg(np.arctan(Vt1/Va1))
    beta1 = np.rad2deg(np.arctan(Wt1/Wa1))

    # THERMODYNAMICS
    # static temperature computation
    T1 = Tt1 - V1**2 / (2*cP)

    ################ EFFICIENCY CORRECTION ################
    # it is assumed that the entropy generation is only on the rotor blade 
    #   this allows to correct thermodynamics only on rotor blade and treating stator blade as isentropic
    #
    # the computed T1 is referred to an isentropic transformation that accounts 
    #   to introduce into the system an amount of energy equal to that necessary
    #   to reach the target betaP with efficiency equal to eta
    #
    # once T1 computed -> if pressure P1 is computed from T1 with the isentropic relation
    #   the pressure obatined (P1) is greater than the real (so that eta < 1) transformation
    # 
    # in order to correct this deficiency:
    # - compute the temperature like the transformation were isentropic 
    # - from this temperature compute the achieved pressure using the isentropic transformation law
    # - compute the density related to the pressure just obtained and the non isentropic temperature 
    #  
    # the changes made by this correction will affect the blade height at the exit 
    #

    if T1real:
        # T1 isoentropic computation 
        T1iso = T0 + eta * (T1 - T0)

        # real pressure computation using the isentropic transformation law and the ideal T1iso temperature
        P1 = Pt0 * (T1iso/Tt0)**(gamma/(gamma-1))
    else: 
        # ideal pressure computation using the isentropic transformation law and the ideal T1 temperature
        P1 = Pt0 * (T1/Tt0)**(gamma/(gamma-1))
    
    # density computation 
    rho1 = P1 / (R * T1)

    # computing other quantities 
    # speed of sound computation 
    a1 = np.sqrt(gamma * R * T1)
    # mach number computation 
    M1 = V1 / a1
    Mr1 = W1 / a1
    # total pressure computation
    Pt1 = P1 * (Tt1/T1)**(gamma/(gamma-1))

    # total density computation
    rhot1 = Pt1 / (R * Tt1)

    ######### STATOR OUTLET QUANTITIES #########
    # main outlet quantities
    # KINETICS 
    Va2 = Va1 
    Wa2 = Va2
    Tt2 = Tt1
    # from thermodynamics
    deltaH = (1 - rD) * L
    T2 = deltaH / cP + T1

    # V2 & Vt2 computation  
    V2 = np.sqrt(2 * cP * (Tt2 - T2))
    # Vt2 computation and numerical correction 
    if np.abs(V2**2 - Va2**2) < 1e-10: 
        Vt2 = 0
    else:
        Vt2 = np.sqrt(V2**2 - Va2**2)

    # W2 & Wt2 computation
    Wt2 = Vt2 - 0 
    W2 = np.sqrt(Wa2**2 + Wt2**2)
    # aerodynamics angles 
    alpha2 = np.rad2deg(np.arctan2(Vt2,Va2))
    beta2 = np.rad2deg(np.arctan(Wt2/Wa2))
    # THERMODYNAMICS
    # sound speed 
    a2 = np.sqrt(gamma * R * T2)
    # mach number computation
    M2 = V2 / a2
    Mr2 = W2 / a2
    # pressure computation
    Pt2 = Pt1
    P2 = P1 * (T2/T1)**(gamma/(gamma-1))
    # density computation 
    rho2 = P2 / (R * T2)
    rhot2 = Pt2 / (R * Tt2)

    ######### BLADE RADIUS #########
    # rotor inlet blade height 
    b0 = mFlux / (rho0 * 2 * np.pi * rMean * Va0)
    # rotor outlet/stator inlet  blade height
    b1 = mFlux / (rho1 * 2 * np.pi * rMean * Va1)
    # stator outlet blade height
    b2 = mFlux / (rho2 * 2 * np.pi * rMean * Va2)

    if printout:
        # print data 
        printLength = 82
        print('*' * printLength)
        if T1real:
            print('\t\tTHE FOLLOWING PAREMETERS ARE REFERRED TO ISENTROPIC\n\t\t  TRANSFORMATION REFERRED TO A WORK L = Lis / eta')
        else:
            print('\t\tTHE FOLLOWING PAREMETERS ARE REFERRED TO REAL\n\t\t  TRANSFORMATION REFERRED TO A WORK L = Lis / eta')
        print('*' * printLength)
        nAdim = int((printLength - len(' ADIMENSIONAL PARAMETERS '))/2)
        print('*' * nAdim + ' ADIMENSIONAL PARAMETERS  ' + '*' * nAdim)
        print('-- rD     = {0:>8.2f}        -- psi    = {1:>8.2f}        -- phi       = {2:>5.2f}'.format(cP * (T1 - T0) / L, psi, phi))
        print('-- eta    = {0:>10.4f}      -- lamdba = {1:>8.2f}        -- Vt0/Umean = {2:>5.2f}'.format(eta, lam, Vt0/Umean))
        nWork = int((printLength - len(' WORK '))/2)
        print('*' * nWork + ' WORK ' + '*' * nWork)
        print('-- Lis    = {0:>8.2f} J/kg   -- L      = {1:>8.2f} J/kg   -- eta    = {2:>10.4f}'.format(Lis, L, eta))
        nRotation = int((printLength - len(' ROTATION '))/2)
        print('*' * nRotation + ' ROTATION ' + '*' * nRotation)
        print('-- Umean  = {0:>8.2f} m/s    -- omega  = {1:>8.2f} rad/s  -- n      = {2:>8.2f} rpm'.format(Umean, omega, n))
        print('*' * printLength)
        nDefinitions = int((printLength - len(' DEFINITIONS '))/2)
        print('\n' + '*' * nDefinitions + ' DEFINITIONS  ' + '*' * nDefinitions)
        print('-- 0 => rotor inlet         -- 1 => rotor outlet        -- 2 => stator outlet')
        nTotTemp = int((printLength - len(' TOTAL TEMPERATURE '))/2)
        print('*' * nTotTemp + ' TOTAL TEMPERATURE  ' + '*' * nTotTemp)
        print('-- Tt0    = {0:>8.2f} K      -- Tt1    = {1:>8.2f} K      -- Tt2    = {2:>8.2f} K'.format(Tt0, Tt1, Tt2))
        nAbsKin = int((printLength - len('KINETICS -- ABSOLUTE '))/2)
        print('*' * nAbsKin + ' KINETICS -- ABSOLUTE ' + '*' * nAbsKin)
        print('-- 0                        -- 1                        -- 2')
        print('-- alpha0 = {0:>8.2f} deg    -- alpha1 = {1:>8.2f} deg    -- alpha2 = {2:>8.2f} deg'.format(alpha0, alpha1, alpha2))
        print('-- Va0    = {0:>8.2f} m/s    -- Va1    = {1:>8.2f} m/s    -- Va2    = {2:>8.2f} m/s'.format(Va0, Va1, Va2))
        print('-- Vt0    = {0:>8.2f} m/s    -- Vt1    = {1:>8.2f} m/s    -- Vt2    = {2:>8.2f} m/s'.format(Vt0, Vt1, Vt2))
        print('-- V0     = {0:>8.2f} m/s    -- V1     = {1:>8.2f} m/s    -- V2     = {2:>8.2f} m/s'.format(V0, V1, V2))
        nRelKin = int((printLength - len('KINETICS -- RELATIVE '))/2)
        print('*' * nRelKin + ' KINETICS -- RELATIVE ' + '*' * nRelKin)
        print('-- 0                        -- 1                        -- 2')
        print('-- U0     = {0:>8.2f} m/s    -- U1     = {1:>8.2f} m/s    -- U2     = {2:>8.2f} m/s'.format(Umean, Umean, 0))
        print('-- beta0  = {0:>8.2f} deg    -- beta1  = {1:>8.2f} deg    -- beta2  = {2:>8.2f} deg'.format(beta0, beta1, beta2))
        print('-- Wa0    = {0:>8.2f} m/s    -- Wa1    = {1:>8.2f} m/s    -- Wa2    = {2:>8.2f} m/s'.format(Wa0, Wa1, Wa2))
        print('-- Wt0    = {0:>8.2f} m/s    -- Wt1    = {1:>8.2f} m/s    -- Wt2    = {2:>8.2f} m/s'.format(Wt0, Wt1, Wt2))
        print('-- W0     = {0:>8.2f} m/s    -- W1     = {1:>8.2f} m/s    -- W2     = {2:>8.2f} m/s'.format(W0, W1, W2))
        nThermo = int((printLength - len(' THERMODYNAMICS '))/2)
        print('*' * nThermo + ' THERMODYNAMICS ' + '*' * nThermo)
        print('-- 0                        -- 1                        -- 2')
        print('-- T0     = {0:>8.2f} K      -- T1     = {1:>8.2f} K      -- T2     = {2:>8.2f} K     '.format(T0, T1, T2))
        print('-- Tt0    = {0:>8.2f} K      -- Tt1    = {1:>8.2f} K      -- Tt2    = {2:>8.2f} K     '.format(Tt0, Tt1, Tt2))
        print('-- a0     = {0:>8.2f} m/s    -- a1     = {1:>8.2f} m/s    -- a2     = {2:>8.2f} m/s   '.format(a0, a1, a2))
        print('-- M0     = {0:>8.2f}        -- M1     = {1:>8.2f}        -- M2     = {2:>8.2f}       '.format(M0, M1, M2))
        print('-- Mr0    = {0:>8.2f}        -- Mr1    = {1:>8.2f}        -- Mr2    = {2:>8.2f}       '.format(Mr0, Mr1, Mr2))
        print('-- P0     = {0:>8.2f} bar    -- P1     = {1:>8.2f} bar    -- P2     = {2:>8.2f} bar   '.format(P0/1e+5, P1/1e+5, P2/1e+5))
        print('-- Pt0    = {0:>8.2f} bar    -- Pt1    = {1:>8.2f} bar    -- Pt2    = {2:>8.2f} bar   '.format(Pt0/1e+5, Pt1/1e+5, Pt2/1e+5))
        print('-- rho0   = {0:>8.2f} kg/m3  -- rho1   = {1:>8.2f} kg/m3  -- rho2   = {2:>8.2f} kg/m3'.format(rho0,rho1,rho2))
        print('-- rhot0  = {0:>8.2f} kg/m3  -- rhot1  = {1:>8.2f} kg/m3  -- rhot2  = {2:>8.2f} kg/m3'.format(rhot0,rhot1,rhot2))
        nBlade = int((printLength - len(' BLADE DIMENSIONS '))/2)
        print('*' * nBlade + ' BLADE DIMENSIONS ' + '*' * nBlade)
        print('-- 0                        -- 1                        -- 2')
        print('-- b0     = {0:>8.2f} cm     -- b1     = {1:>8.2f} cm     -- b2     = {2:>8.2f} cm'.format(b0*1e+2, b1*1e+2, b2*1e+2))
        print('-- rTip0  = {0:>8.2f} cm     -- rTip1  = {1:>8.2f} cm     -- rTip2  = {2:>8.2f} cm'.format((rMean + b0/2)*1e+2, (rMean + b1/2)*1e+2, (rMean + b2/2)*1e+2))
        print('-- rMean0 = {0:>8.2f} cm     -- rMean1 = {0:>8.2f} cm     -- rMean  = {0:>8.2f} cm'.format(rMean))
        print('-- rHub0  = {0:>8.2f} cm     -- rHub1  = {1:>8.2f} cm     -- rHub2  = {2:>8.2f} cm'.format((rMean - b0/2)*1e+2, (rMean - b1/2)*1e+2, (rMean - b2/2)*1e+2))
        print('*' * printLength)

    # return vector values 
    adimVec     = [phi, psi, eta]
    bladeVec    = [b0, b1, b2]
    rotationVec = [Umean, omega, n] 
    V0vec       = [Va0, Vt0, alpha0]
    V1vec       = [Va1, Vt1, alpha1]
    V2vec       = [Va2, Vt2, alpha2]
    W0vec       = [Wa0, Wt0, beta0]
    W1vec       = [Wa1, Wt1, beta1]
    W2vec       = [Wa2, Wt2, beta2] 
    thermo0     = [T0, P0, rho0, Tt0, Pt0, rhot0, M0, Mr0] 
    thermo1     = [T1, P1, rho1, Tt1, Pt1, rhot1, M1, Mr1]
    thermo2     = [T2, P2, rho2, Tt2, Pt2, rhot2, M2, Mr2] 
    work        = [L, Lis]

    return adimVec, bladeVec, rotationVec, V0vec, V1vec, V2vec, W0vec, W1vec, W2vec, thermo0, thermo1, thermo2, work

def stageStudy(mFlux, betaP, rMean, Pt0, Tt0, rDmin=0.5, rDmax=0.75, Vt0UmeanMin=0, Vt0UmeanMax=0.25, R=287.06, gamma=1.4):
    '''
    This function allows to understand the behaviour of the stage with respect to different changes:
        -- reaction degree rD
        -- Vt / Umean 
        -- the rotor properties are computed with respect to the efficiency computed in 'efficiency' function
        -- the fluid properties at the end of the rotor blade are referred to an isentropic transformation 
            this means that the study of the stage is made by L = Lis / eta but the all the quantities are computed 
            following an isentropic process 
                -- the implication of this is that the exit properties of the rotor are not the real ones if 
                    taken into account the efficiency also for the computation of pressure 
                -- the pressure computed in this case is like the stator follows an insentropic process that allows to 
                    have a defined betaP that takes into account the pressure losses in a real stage
                -- to sum up:
                    * work -> referred to real stage -> L = Lis / eta
                    * pressure -> referred to ideal stage such that the work is equal to L 

        inputs:
            mFlux           -- mass flux 
            betaP           -- total pressure ratio
            rMean           -- mean radius 
            Pt0             -- initial total pressure 
            Tt0             -- initial total temperature
            rDmin           -- minimum reaction degree for the rDarray generation 
            rDmax           -- maximum reaction degree for the rDarray generation 
            Vt0UmeanMin     -- minimum Vt0/Umean value for the Vt0UmeanArray generation 
            Vt0UmeanMax     -- maximum Vt0/Umean value for the Vt0UmeanArray generation 
            R               -- gas constant
            gamma           -- specifi heat ratio
        output:
            plots  
    '''

    # vector allocation for the properties study 
    rDarray = np.linspace(rDmin, rDmax, 12)
    Vt0UmeanArray = np.linspace(Vt0UmeanMin, Vt0UmeanMax, 60)

    # setting up figure and subplots 
    fig0, [[ax0, ax1, ax2, ax3], [ax4, ax5, ax6, ax7]] = plt.subplots(figsize=(20,8), nrows=2, ncols=4)

    for _,rD in enumerate(rDarray):

        # vector allocation 
        Vt0UmeanVec = np.array([])
        phiVec      = np.array([])
        psiVec      = np.array([])
        etaVec      = np.array([])
        alpha0Vec   = np.array([])
        UmeanVec    = np.array([])
        Mr0Vec      = np.array([])
        b0Vec       = np.array([])
        b1Vec       = np.array([])
        b2Vec       = np.array([])

        for _,Vt0Umean in enumerate(Vt0UmeanArray):
            # Vt0 = (1 - rD - lam/4) * Umean -> lam = 4 * (1 - rD - Vt0/Umean) 
            # psi = lam / 2 
            lam = (1 - rD - Vt0Umean) * 4 
            psiTarget = lam / 2

            try:                
                # variable needed for the study computation
                adimVec, bladeVec, rotationVec, V0vec, _, _, _, _, _, thermo0, _, _, _ = stageProperties(rD, psiTarget, rMean, mFlux, Tt0, Pt0, betaP, T1real=False, R=R, gamma=gamma)

                # variable allocation -> check stageProperties return  
                phi    = adimVec[0]
                psi    = adimVec[1]
                eta    = adimVec[2]
                alpha0 = V0vec[2]
                Umean  = rotationVec[0]
                Mr0    = thermo0[7]
                b0     = bladeVec[0]
                b1     = bladeVec[1]
                b2     = bladeVec[2]

                # appenda data to vectors 
                Vt0UmeanVec = np.append(Vt0UmeanVec, Vt0Umean)
                phiVec      = np.append(phiVec, phi)
                psiVec      = np.append(psiVec, psi)
                etaVec      = np.append(etaVec, eta)
                alpha0Vec   = np.append(alpha0Vec, alpha0)
                UmeanVec    = np.append(UmeanVec, Umean)
                Mr0Vec      = np.append(Mr0Vec, Mr0)
                b0Vec       = np.append(b0Vec, b0)
                b1Vec       = np.append(b1Vec, b1)
                b2Vec       = np.append(b2Vec, b2)

            except:
                pass
        
        # plotting data
        ax0.plot(Vt0UmeanVec, phiVec,          label=r'$\chi = {0:.2f}$'.format(rD))
        ax1.plot(Vt0UmeanVec, psiVec,          label=r'$\chi = {0:.2f}$'.format(rD))
        ax2.plot(Vt0UmeanVec, etaVec,          label=r'$\chi = {0:.2f}$'.format(rD))
        ax3.plot(Vt0UmeanVec, alpha0Vec,       label=r'$\chi = {0:.2f}$'.format(rD))
        ax4.plot(Vt0UmeanVec, UmeanVec,        label=r'$\chi = {0:.2f}$'.format(rD))
        ax5.plot(Vt0UmeanVec, Mr0Vec,          label=r'$\chi = {0:.2f}$'.format(rD))
        ax6.plot(Vt0UmeanVec, b0Vec,           label=r'$\chi = {0:.2f}$'.format(rD))
        ax7.plot(Vt0UmeanVec, b0Vec/2 + rMean, label=r'$\chi = {0:.2f}$'.format(rD))
    
    # ax0 setup
    ax0.set_title(r'$r_{{mean}} = {0:.3} m$'.format(rMean))
    ax0.set_xlabel(r'$\frac{V_{t0}}{U_{mean}}$')
    ax0.set_ylabel(r'$\phi$')
    ax0.grid(linestyle='--')
    ax0.legend(loc='upper left', bbox_to_anchor=[1,1])

    # ax1 setup 
    ax1.set_title(r'$r_{{mean}} = {0:.3} m$'.format(rMean))
    ax1.set_xlabel(r'$\frac{V_{t0}}{U_{mean}}$')
    ax1.set_ylabel(r'$\psi$')
    ax1.grid(linestyle='--')
    ax1.legend(loc='upper left', bbox_to_anchor=[1,1])

    # ax2 setup
    ax2.set_title(r'$r_{{mean}} = {0:.3} m$'.format(rMean))
    ax2.set_xlabel(r'$\frac{V_{t0}}{U_{mean}}$')
    ax2.set_ylabel(r'$\eta$')
    ax2.grid(linestyle='--')
    ax2.legend(loc='upper left', bbox_to_anchor=[1,1])

    # ax3 setup
    ax3.set_title(r'$r_{{mean}} = {0:.3} m$'.format(rMean))
    ax3.set_xlabel(r'$\frac{V_{t0}}{U_{mean}}$')
    ax3.set_ylabel(r'$\alpha_0$')
    ax3.grid(linestyle='--')
    ax3.legend(loc='upper left', bbox_to_anchor=[1,1])

    # ax4 setup
    ax4.set_title(r'$r_{{mean}} = {0:.3} m$'.format(rMean))
    ax4.set_xlabel(r'$\frac{V_{t0}}{U_{mean}}$')
    ax4.set_ylabel(r'$U_{{mean}}$')
    ax4.grid(linestyle='--')
    ax4.legend(loc='upper left', bbox_to_anchor=[1,1])

    # ax5 setup
    ax5.set_title(r'$r_{{mean}} = {0:.3} m$'.format(rMean))
    ax5.set_xlabel(r'$\frac{V_{t0}}{U_{mean}}$')
    ax5.set_ylabel(r'$M_{r0}$')
    ax5.grid(linestyle='--')
    ax5.legend(loc='upper left', bbox_to_anchor=[1,1])

    # ax6 setup
    ax6.set_title(r'$r_{{mean}} = {0:.3} m$'.format(rMean))
    ax6.set_xlabel(r'$\frac{V_{t0}}{U_{mean}}$')
    ax6.set_ylabel(r'$b_0$')
    ax6.grid(linestyle='--')
    ax6.legend(loc='upper left', bbox_to_anchor=[1,1])

    # ax7 setup
    ax7.set_title(r'$r_{{mean}} = {0:.3} m$'.format(rMean))
    ax7.set_xlabel(r'$\frac{V_{t0}}{U_{mean}}$')
    ax7.set_ylabel(r'$r_{{tip}}$')
    ax7.grid(linestyle='--')
    ax7.legend(loc='upper left', bbox_to_anchor=[1,1])

    # xticks setup 
    ax0.set_xticks(np.arange(Vt0UmeanArray.min(), Vt0UmeanArray.max()+0.05, 0.05))
    ax1.set_xticks(np.arange(Vt0UmeanArray.min(), Vt0UmeanArray.max()+0.05, 0.05))
    ax2.set_xticks(np.arange(Vt0UmeanArray.min(), Vt0UmeanArray.max()+0.05, 0.05))
    ax3.set_xticks(np.arange(Vt0UmeanArray.min(), Vt0UmeanArray.max()+0.05, 0.05))
    ax4.set_xticks(np.arange(Vt0UmeanArray.min(), Vt0UmeanArray.max()+0.05, 0.05))
    ax5.set_xticks(np.arange(Vt0UmeanArray.min(), Vt0UmeanArray.max()+0.05, 0.05))
    ax6.set_xticks(np.arange(Vt0UmeanArray.min(), Vt0UmeanArray.max()+0.05, 0.05))
    ax7.set_xticks(np.arange(Vt0UmeanArray.min(), Vt0UmeanArray.max()+0.05, 0.05))
    
    # tight layout
    fig0.tight_layout()

    # plotting data
    plt.show()

def reactionStudy(mFlux, betaP, rMean, Pt0, Tt0, rDmin=0.5, rDmax=0.73, Vt0UmeanMin=0, Vt0UmeanMax=0.25, save=False, position0='reactionStudy0.pgf', position1='reactionStudy1.pgf'):
    '''
    This function allows the study of the rotor tip and rotor hub given as input some constraints.
        inputs:
            ** constraints ** 
            mFlux       -- mass flux 
            betaP       -- total pressure ratio 
            rMean       -- blade mean radius
            Pt0         -- inlet total pressure 
            Tt0         -- inlet total temperature 
            ** for analysis ** 
            rDmin       -- minimum value of study for the reaction degree
            rDmax       -- maximum value of study for the reaciton degree 
            Vt0UmeanMin -- Vt0/Umean minimum value of study
            Vt0UmeanMax -- Vt0/Umean maximum value of study 
            save        -- boolean value for the saving of the plots 
            position0   -- path where to save the first figure 
            position1   -- path where to save the second figure
    '''

    # importing libraries 
    from turboClass.bladeStudy import bladeStudy 

    # vector allocation for the properties study 
    rDarray = np.linspace(rDmin, rDmax, 10)
    Vt0UmeanArray = np.linspace(Vt0UmeanMin, Vt0UmeanMax, 20)

    # figure allocation
    fig0, [[ax1_11, ax1_12, ax1_13, ax1_14, ax1_15], [ax1_21, ax1_22, ax1_23, ax1_24, ax1_25]] = plt.subplots(figsize=(20,10), nrows=2, ncols=5)
    fig1, [[ax2_11, ax2_12, ax2_13, ax2_14, ax2_15, ax2_16], [ax2_21, ax2_22, ax2_23, ax2_24, ax2_25, ax2_26]] = plt.subplots(figsize=(20,10), nrows=2, ncols=6)
    for _,rD in enumerate(rDarray):

        # vector allocation 
        rDvecHub       = np.array([])
        rDvecTip       = np.array([])
        Vt0UmeanVecHub = np.array([])
        Vt0UmeanVecTip = np.array([])
        alpha0vecHub   = np.array([])
        alpha1vecHub   = np.array([])
        beta0vecHub    = np.array([])
        beta1vecHub    = np.array([])
        alpha0vecTip   = np.array([])
        alpha1vecTip   = np.array([])
        beta0vecTip    = np.array([])
        beta1vecTip    = np.array([])
        M0vecTip       = np.array([])
        M0vecHub       = np.array([])
        M1vecTip       = np.array([])
        M1vecHub       = np.array([])
        Mr0vecTip      = np.array([])
        Mr0vecHub      = np.array([])
        Mr1vecTip      = np.array([])
        Mr1vecHub      = np.array([])

        for _,Vt0Umean in enumerate(Vt0UmeanArray):
            # Vt0 = (1 - rD - lam/4) * Umean -> lam = 4 * (1 - rD - Vt0/Umean) 
            # psi = lam / 2 
            lam = (1 - rD - Vt0Umean) * 4 
            psiTarget = lam / 2

            # plotting charts 
            try:
                ###### MEAN LINE STUDY ###### 
                adimVec, bladeVec, rotationVec, V0vec, V1vec, _, _, _, _, thermo0, _, _, work = stageProperties(rD, psiTarget, rMean, mFlux, Tt0, Pt0, betaP, T1real=False)
                
                # data allocation
                eta   = adimVec[-1]
                b0    = bladeVec[0]
                b1    = bladeVec[1]
                omega = rotationVec[1]
                T0    = thermo0[0]
                P0    = thermo0[1]
                Va0   = V0vec[0]
                Vt0   = V0vec[1]
                Vt1   = V1vec[1]
                L     = work[0]

                ###### BLADE HUB STUDY ######
                # data allocation for 'bladeStudy'
                rIn       = rMean - b0/2
                rOut      = rMean - b1/2
                VaMean    = Va0 
                VtMeanIn  = Vt0 
                VtMeanOut = Vt1 
                Leu       = L 
                
                # computing blade properties 
                adimVec, rotationVec, _, _, _, _, angleVec, thermo0, thermo1 = bladeStudy(rIn, rOut, omega, rMean, VaMean, VtMeanIn, VtMeanOut, Leu, Tt0, T0, Pt0, P0, eta=eta, printout=False)

                # vector update 
                Vt0UmeanVecHub  = np.append(Vt0UmeanVecHub, Vt0Umean)
                rDvecHub        = np.append(rDvecHub, adimVec[0]) 
                alpha0vecHub    = np.append(alpha0vecHub, angleVec[0])
                alpha1vecHub    = np.append(alpha1vecHub, angleVec[1])
                beta0vecHub     = np.append(beta0vecHub, angleVec[2])
                beta1vecHub     = np.append(beta1vecHub, angleVec[3])
                M0vecHub        = np.append(M0vecHub, thermo0[-2])
                M1vecHub        = np.append(M1vecHub, thermo1[-2])
                Mr0vecHub       = np.append(Mr0vecHub, thermo0[-1])
                Mr1vecHub       = np.append(Mr1vecHub, thermo1[-1])

                ###### BLADE TIP STUDY ######
                # data allocation for 'bladeStudy'
                rIn       = rMean + b0/2
                rOut      = rMean + b1/2
                VaMean    = Va0 
                VtMeanIn  = Vt0 
                VtMeanOut = Vt1 
                Leu       = L 
                
                # computing blade properties 
                adimVec, rotationVec, _, _, _, _, angleVec, thermo0, thermo1 = bladeStudy(rIn, rOut, omega, rMean, VaMean, VtMeanIn, VtMeanOut, Leu, Tt0, T0, Pt0, P0, eta=eta, printout=False)

                # vector update
                Vt0UmeanVecTip = np.append(Vt0UmeanVecTip, Vt0Umean)
                rDvecTip       = np.append(rDvecTip, adimVec[0])
                alpha0vecTip   = np.append(alpha0vecTip, angleVec[0])
                alpha1vecTip   = np.append(alpha1vecTip, angleVec[1])
                beta0vecTip    = np.append(beta0vecTip, angleVec[2])
                beta1vecTip    = np.append(beta1vecTip, angleVec[3])
                M0vecTip       = np.append(M0vecTip, thermo0[-2])
                M1vecTip       = np.append(M1vecTip, thermo1[-2])
                Mr0vecTip      = np.append(Mr0vecTip, thermo0[-1])
                Mr1vecTip      = np.append(Mr1vecTip, thermo1[-1])

            except:
                pass 
        
        # plotting results
        # hub plot 
        ax1_11.plot(Vt0UmeanVecHub, rDvecHub,                     )#label=r'$\chi_{{mean}} = {0:.2f}$'.format(rD))
        ax1_12.plot(Vt0UmeanVecHub, M0vecHub,                     )#label=r'$\chi_{{mean}} = {0:.2f}$'.format(rD))
        ax1_13.plot(Vt0UmeanVecHub, M1vecHub,                     )#label=r'$\chi_{{mean}} = {0:.2f}$'.format(rD))
        ax1_14.plot(Vt0UmeanVecHub, Mr0vecHub,                    )#label=r'$\chi_{{mean}} = {0:.2f}$'.format(rD))
        ax1_15.plot(Vt0UmeanVecHub, Mr1vecHub,                    )#label=r'$\chi_{{mean}} = {0:.2f}$'.format(rD))
        ax2_11.plot(Vt0UmeanVecHub, alpha0vecHub,                 )#label=r'$\chi_{{mean}} = {0:.2f}$'.format(rD))
        ax2_12.plot(Vt0UmeanVecHub, alpha1vecHub,                 )#label=r'$\chi_{{mean}} = {0:.2f}$'.format(rD))
        ax2_13.plot(Vt0UmeanVecHub, beta0vecHub,                  )#label=r'$\chi_{{mean}} = {0:.2f}$'.format(rD))
        ax2_14.plot(Vt0UmeanVecHub, beta1vecHub,                  )#label=r'$\chi_{{mean}} = {0:.2f}$'.format(rD))
        ax2_15.plot(Vt0UmeanVecHub, alpha1vecHub - alpha0vecHub,  )#label=r'$\chi_{{mean}} = {0:.2f}$'.format(rD))
        ax2_16.plot(Vt0UmeanVecHub, beta1vecHub - beta0vecHub,    )#label=r'$\chi_{{mean}} = {0:.2f}$'.format(rD))
        # tip plot
        ax1_21.plot(Vt0UmeanVecTip, rDvecTip,                     )#label=r'$\chi_{{mean}} = {0:.2f}$'.format(rD))
        ax1_22.plot(Vt0UmeanVecTip, M0vecTip,                     )#label=r'$\chi_{{mean}} = {0:.2f}$'.format(rD))
        ax1_23.plot(Vt0UmeanVecTip, M1vecTip,                     )#label=r'$\chi_{{mean}} = {0:.2f}$'.format(rD))
        ax1_24.plot(Vt0UmeanVecTip, Mr0vecTip,                    )#label=r'$\chi_{{mean}} = {0:.2f}$'.format(rD))
        ax1_25.plot(Vt0UmeanVecTip, Mr1vecTip,                    )#label=r'$\chi_{{mean}} = {0:.2f}$'.format(rD))
        ax2_21.plot(Vt0UmeanVecTip, alpha0vecTip,                 )#label=r'$\chi_{{mean}} = {0:.2f}$'.format(rD))
        ax2_22.plot(Vt0UmeanVecTip, alpha1vecTip,                 )#label=r'$\chi_{{mean}} = {0:.2f}$'.format(rD))
        ax2_23.plot(Vt0UmeanVecTip, beta0vecTip,                  )#label=r'$\chi_{{mean}} = {0:.2f}$'.format(rD))
        ax2_24.plot(Vt0UmeanVecTip, beta1vecTip,                  )#label=r'$\chi_{{mean}} = {0:.2f}$'.format(rD))
        ax2_25.plot(Vt0UmeanVecTip, alpha1vecTip - alpha0vecTip,  )#label=r'$\chi_{{mean}} = {0:.2f}$'.format(rD))
        ax2_26.plot(Vt0UmeanVecTip, beta1vecTip - beta0vecTip,    )#label=r'$\chi_{{mean}} = {0:.2f}$'.format(rD))

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


    ax1_11.set_xlabel(r'$\frac{V_{t0}}{U_{{mean}}}$')
    ax1_11.set_ylabel(r'$\chi \ @ \ r_{{hub}}$')
    ax1_11.grid(linestyle='--')

    ax1_12.set_xlabel(r'$\frac{V_{t0}}{U_{{mean}}}$')
    ax1_12.set_ylabel(r'$M_{{0}} \ @ \ r_{{hub}}$')
    ax1_12.grid(linestyle='--')

    ax1_13.set_xlabel(r'$\frac{V_{t0}}{U_{{mean}}}$')
    ax1_13.set_ylabel(r'$M_{{1}} \ @ \ r_{{hub}}$')
    ax1_13.grid(linestyle='--')

    ax1_14.set_xlabel(r'$\frac{V_{t0}}{U_{{mean}}}$')
    ax1_14.set_ylabel(r'$M_{{r0}} \ @ \ r_{{hub}}$')
    ax1_14.grid(linestyle='--')

    ax1_15.set_xlabel(r'$\frac{V_{t0}}{U_{{mean}}}$')
    ax1_15.set_ylabel(r'$M_{{r1}} \ @ \ r_{{hub}}$')
    ax1_15.grid(linestyle='--')

    ax2_11.set_xlabel(r'$\frac{V_{t0}}{U_{{mean}}}$')
    ax2_11.set_ylabel(r'$\alpha_{{0}} \ @ \ r_{{hub}}$')
    ax2_11.grid(linestyle='--')

    ax2_12.set_xlabel(r'$\frac{V_{t0}}{U_{{mean}}}$')
    ax2_12.set_ylabel(r'$\alpha_{{1}} \ @ \ r_{{hub}}$')
    ax2_12.grid(linestyle='--')

    ax2_13.set_xlabel(r'$\frac{V_{t0}}{U_{{mean}}}$')
    ax2_13.set_ylabel(r'$\beta_{{0}} \ @ \ r_{{hub}}$')
    ax2_13.grid(linestyle='--')

    ax2_14.set_xlabel(r'$\frac{V_{t0}}{U_{{mean}}}$')
    ax2_14.set_ylabel(r'$\beta_{{1}} \ @ \ r_{{hub}}$')
    ax2_14.grid(linestyle='--')

    ax2_15.set_xlabel(r'$\frac{V_{t0}}{U_{{mean}}}$')
    ax2_15.set_ylabel(r'$\Delta \alpha \ @ \ r_{{hub}}$')
    ax2_15.grid(linestyle='--')

    ax2_16.set_xlabel(r'$\frac{V_{t0}}{U_{{mean}}}$')
    ax2_16.set_ylabel(r'$\Delta \beta \ @ \ r_{{hub}}$')
    ax2_16.grid(linestyle='--')

    ax1_21.set_xlabel(r'$\frac{V_{t0}}{U_{{mean}}}$')
    ax1_21.set_ylabel(r'$\chi \ @ \ r_{{tip}}$')
    ax1_21.grid(linestyle='--')

    ax1_22.set_xlabel(r'$\frac{V_{t0}}{U_{{mean}}}$')
    ax1_22.set_ylabel(r'$M_{{0}} \ @ \ r_{{tip}}$')
    ax1_22.grid(linestyle='--')

    ax1_23.set_xlabel(r'$\frac{V_{t0}}{U_{{mean}}}$')
    ax1_23.set_ylabel(r'$M_{{1}} \ @ \ r_{{tip}}$')
    ax1_23.grid(linestyle='--')

    ax1_24.set_xlabel(r'$\frac{V_{t0}}{U_{{mean}}}$')
    ax1_24.set_ylabel(r'$M_{{r0}} \ @ \ r_{{tip}}$')
    ax1_24.grid(linestyle='--')

    ax1_25.set_xlabel(r'$\frac{V_{t0}}{U_{{mean}}}$')
    ax1_25.set_ylabel(r'$M_{{r1}} \ @ \ r_{{tip}}$')
    ax1_25.grid(linestyle='--')

    ax2_21.set_xlabel(r'$\frac{V_{t0}}{U_{{mean}}}$')
    ax2_21.set_ylabel(r'$\alpha_{{0}} \ @ \ r_{{tip}}$')
    ax2_21.grid(linestyle='--')

    ax2_22.set_xlabel(r'$\frac{V_{t0}}{U_{{mean}}}$')
    ax2_22.set_ylabel(r'$\alpha_{{1}} \ @ \ r_{{tip}}$')
    ax2_22.grid(linestyle='--')

    ax2_23.set_xlabel(r'$\frac{V_{t0}}{U_{{mean}}}$')
    ax2_23.set_ylabel(r'$\beta_{{0}} \ @ \ r_{{tip}}}$')
    ax2_23.grid(linestyle='--')

    ax2_24.set_xlabel(r'$\frac{V_{t0}}{U_{{mean}}}$')
    ax2_24.set_ylabel(r'$\beta_{{1}} \ @ \ r_{{tip}}$')
    ax2_24.grid(linestyle='--')

    ax2_25.set_xlabel(r'$\frac{V_{t0}}{U_{{mean}}}$')
    ax2_25.set_ylabel(r'$\Delta \alpha \ @ \ r_{{tip}}$')
    ax2_25.grid(linestyle='--')

    ax2_26.set_xlabel(r'$\frac{V_{t0}}{U_{{mean}}}$')
    ax2_26.set_ylabel(r'$\Delta \beta \ @ \ r_{{tip}}$')
    ax2_26.grid(linestyle='--')

    #fig0.suptitle(r'$r_{{mean}} = {0:.3f} m$   $\dot{{m}} = {1:.2f} \frac{{kg}}{{s}}$    $\beta_{{T}} = {2:.2f}$'.format(rMean, mFlux, betaP))
    #fig1.suptitle(r'$r_{{mean}} = {0:.3f} m$   $\dot{{m}} = {1:.2f} \frac{{kg}}{{s}}$    $\beta_{{T}} = {2:.2f}$'.format(rMean, mFlux, betaP))
    fig0.suptitle(' ')
    fig1.suptitle(' ')

    leg = []
    for ii in range(len(rDarray)):
        leg.append(r'$\chi = {:.2f}$'.format(rDarray[ii]))

    fig0.legend(leg, loc='upper center', ncol=len(leg))
    fig1.legend(leg, loc='upper center', ncol=len(leg))
    fig0.tight_layout()
    fig1.tight_layout()

    if save:
        # figure saving 
        fig0.savefig(position0, bbox_inches='tight')
        # figure saving 
        fig1.savefig(position1, bbox_inches='tight')

    plt.show()
