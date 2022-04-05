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

def stageStudy(mFlux, betaP, rMean, Pt0, Tt0):
    '''
    This function allows to understand the behaviour of the stage with respect to different changes:
        -- reaction degree rD
        -- Vt / Umean 

        inputs:
            mFlux   -- mass flux 
            betaP   -- total pressure ratio
            rMean   -- mean radius 
            Pt0     -- initial total pressure 
            Tt0     -- initial total temperature 
    '''

    # importing libraries 
    from turboCoeff import coeff 

    # air properties allocation 
    R     = 287.06                  # air gas constant Ru / Mm      [J/kg K]
    gamma = 1.4                     # specific heat ratio           [--]
    cP    = gamma / (gamma - 1) * R # specific heat ratio @ P cost  [J/kg K]
    
    # vector allocation for the properties study 
    rDarray = np.linspace(0.45, 0.65, 12)
    Vt0UmeanArray = np.linspace(0.05, 0.35, 20)

    fig0, [[ax0, ax1, ax2, ax3], [ax4, ax5, ax6, ax7]] = plt.subplots(figsize=(20,8), nrows=2, ncols=4)

    for _,rD in enumerate(rDarray):

        # vector allocation 
        Vt0UmeanVec = np.array([])
        phiVec = np.array([])
        psiVec = np.array([])
        etaVec = np.array([])
        alpha0Vec = np.array([])
        UmeanVec = np.array([])
        rD_etaVec = np.array([])
        Mr0Vec = np.array([])
        b0Vec = np.array([])
        b1Vec = np.array([])
        b2Vec = np.array([])

        for _,Vt0Umean in enumerate(Vt0UmeanArray):
            # Vt0 = (1 - rD - lam/4) * Umean -> lam = 4 * (1 - rD - Vt0/Umean) 
            # psi = lam / 2 
            lam = (1 - rD - Vt0Umean) * 4 
            psiTarget = lam / 2

            # plotting charts 
            try:
                phi, psi = stagePerf(psi=psiTarget, rD=rD, plot=False, perc=0.97)
                eta = efficiency(phi=phi, rD=rD, plot=False)

                # compute flow properties 
                # WORK
                # ideal compression work 
                Lis = coeff.L_is(Tin=Tt0, beta=betaP, gamma=gamma, kind='compressor')
                # real compression work 
                L = Lis / eta
                # total temperature computation
                Tt1 = L / cP + Tt0

                # ROTATION 
                # mean section revolution speed from work coefficient
                Umean = np.sqrt(L / psi) 
                # angular speed 
                omega = Umean / rMean
                # rpm 
                n = omega * 60 / (2 * np.pi)

                # MAIN VELOCITIES COMPUTATION FROM WORK RESULTS
                # Va0, Vt0 computation 
                # axial inlet velocity from flow coefficient
                Va0 = phi * Umean
                # tangential inlet velocity from reaction degree and work coefficient
                Vt0 = (1 - rD - lam/4) * Umean 
                # under Umean assumption & Vt0 assumption
                Vt1 = (1 - rD + lam/4) * Umean 
                Vtinf = (Vt1 + Vt0)/2

                # ROTOR INLET QUANTITIES 
                # main inlet quantities computation
                # kinetics 
                # relative speed computation 
                Wa0 = Va0
                Wt0 = Vt0 - Umean
                # veolcity magnitude computation
                W0 = np.sqrt(Wa0**2 + Wt0**2)
                V0 = np.sqrt(Va0**2 + Vt0**2)
                # aerodynamic angles computation 
                alpha0 = np.rad2deg(np.arctan(Vt0/Va0))
                beta0 = np.rad2deg(np.arctan(Wt0/Wa0))
                # thermodynamics
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

                # ROTOR OUTLET/STATOR INLET QUANTITIES  
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
                # thermodynamics
                # static temperature computation
                T1 = Tt1 - V1**2 / (2*cP)
                # speed of sound computation 
                a1 = np.sqrt(gamma * R * T1)
                # mach number computation 
                M1 = V1 / a1
                Mr1 = W1 / a1
                # pressure computation
                P1 = P0 * (T1/T0)**(gamma/(gamma-1))
                Pt1 = P1 * (Tt1/T1)**(gamma/(gamma-1))
                # density computation
                rhot1 = Pt1 / (R * Tt1)
                rho1 = rhot1 / (1 + (gamma - 1)/2 * M1**2)**(1/(gamma-1))

                # STATOR OUTLET QUANTITIES 
                # main outlet quantities 
                Va2 = Va1 
                Pt2 = Pt1
                Tt2 = Tt1
                deltaH = (1 - rD) * L
                T2 = deltaH / cP + T1 
                V2 = np.sqrt(2 * cP * (Tt2 - T2))
                Vt2 = np.sqrt(V2**2 - Va2**2)
                alpha2 = np.rad2deg(np.arctan2(Vt2,Va2))
                rhot2 = Pt2 / (R * Tt2)
                # thermodynamics
                a2 = np.sqrt(gamma * R * T2)
                M2 = V2 / a2
                P2 = P1 * (T2/T1)**(gamma/(gamma-1))
                rho2 = rhot2 / (1 + (gamma - 1)/2 * M2**2)**(1/(gamma-1))

                # BLADE RADIUS
                # rotor inlet blade height 
                b0 = mFlux / (rho0 * 2 * np.pi * rMean * Va0)
                # rotor outlet/stator inlet  blade height
                b1 = mFlux / (rho1 * 2 * np.pi * rMean * Va1)
                # stator outlet blade height
                b2 = mFlux / (rho2 * 2 * np.pi * rMean * Va2)

                # reaction degree with efficiency 
                rD_eta = cP * (T1 - T0) / L 

                # vector 
                Vt0UmeanVec = np.append(Vt0UmeanVec, Vt0Umean)
                phiVec = np.append(phiVec, phi)
                psiVec = np.append(psiVec, psi)
                etaVec = np.append(etaVec, eta)
                alpha0Vec = np.append(alpha0Vec, alpha0)
                UmeanVec = np.append(UmeanVec, Umean)
                rD_etaVec = np.append(rD_etaVec, rD_eta)
                Mr0Vec = np.append(Mr0Vec, Mr0)
                b0Vec = np.append(b0Vec, b0)
                b1Vec = np.append(b1Vec, b1)
                b2Vec = np.append(b2Vec, b2)

            except:
                pass

        ax0.plot(Vt0UmeanVec, phiVec, label=r'$\chi = {0:.2f}$'.format(rD))
        ax1.plot(Vt0UmeanVec, psiVec, label=r'$\chi = {0:.2f}$'.format(rD))
        ax2.plot(Vt0UmeanVec, etaVec, label=r'$\chi = {0:.2f}$'.format(rD))
        ax3.plot(Vt0UmeanVec, alpha0Vec, label=r'$\chi = {0:.2f}$'.format(rD))
        ax4.plot(Vt0UmeanVec, UmeanVec, label=r'$\chi = {0:.2f}$'.format(rD))
        #ax4.plot(Vt0UmeanVec, rD_etaVec, label=r'$\chi = {0:.2f}$'.format(rD))
        ax5.plot(Vt0UmeanVec, Mr0Vec, label=r'$\chi = {0:.2f}$'.format(rD))
        ax6.plot(Vt0UmeanVec, b0Vec, label=r'$\chi = {0:.2f}$'.format(rD))
        ax7.plot(Vt0UmeanVec, b0Vec/2 + rMean, label=r'$\chi = {0:.2f}$'.format(rD))
    
    # ax0 setup
    ax0.set_title(r'$r_{{mean}} = {0:.2} m$'.format(rMean))
    ax0.set_xlabel(r'$\frac{V_{t0}}{U_{mean}}$')
    ax0.set_ylabel(r'$\phi$')
    ax0.grid(linestyle='--')
    ax0.legend(loc='upper left', bbox_to_anchor=[1,1])

    # ax1 setup 
    ax1.set_title(r'$r_{{mean}} = {0:.2} m$'.format(rMean))
    ax1.set_xlabel(r'$\frac{V_{t0}}{U_{mean}}$')
    ax1.set_ylabel(r'$\psi$')
    ax1.grid(linestyle='--')
    ax1.legend(loc='upper left', bbox_to_anchor=[1,1])

    # ax2 setup
    ax2.set_title(r'$r_{{mean}} = {0:.2} m$'.format(rMean))
    ax2.set_xlabel(r'$\frac{V_{t0}}{U_{mean}}$')
    ax2.set_ylabel(r'$\eta$')
    ax2.grid(linestyle='--')
    ax2.legend(loc='upper left', bbox_to_anchor=[1,1])

    # ax3 setup
    ax3.set_title(r'$r_{{mean}} = {0:.2} m$'.format(rMean))
    ax3.set_xlabel(r'$\frac{V_{t0}}{U_{mean}}$')
    ax3.set_ylabel(r'$\alpha_0$')
    ax3.grid(linestyle='--')
    ax3.legend(loc='upper left', bbox_to_anchor=[1,1])

    # ax4 setup
    ax4.set_title(r'$r_{{mean}} = {0:.2} m$'.format(rMean))
    ax4.set_xlabel(r'$\frac{V_{t0}}{U_{mean}}$')
    ax4.set_ylabel(r'$U_{{mean}}$')
    ax4.grid(linestyle='--')
    ax4.legend(loc='upper left', bbox_to_anchor=[1,1])

    # ax5 setup
    ax5.set_title(r'$r_{{mean}} = {0:.2} m$'.format(rMean))
    ax5.set_xlabel(r'$\frac{V_{t0}}{U_{mean}}$')
    ax5.set_ylabel(r'$M_{r0}$')
    ax5.grid(linestyle='--')
    ax5.legend(loc='upper left', bbox_to_anchor=[1,1])

    # ax6 setup
    ax6.set_title(r'$r_{{mean}} = {0:.2} m$'.format(rMean))
    ax6.set_xlabel(r'$\frac{V_{t0}}{U_{mean}}$')
    ax6.set_ylabel(r'$b_0$')
    ax6.grid(linestyle='--')
    ax6.legend(loc='upper left', bbox_to_anchor=[1,1])

    # ax7 setup
    ax7.set_title(r'$r_{{mean}} = {0:.2} m$'.format(rMean))
    ax7.set_xlabel(r'$\frac{V_{t0}}{U_{mean}}$')
    ax7.set_ylabel(r'$r_{{tip}}$')
    ax7.grid(linestyle='--')
    ax7.legend(loc='upper left', bbox_to_anchor=[1,1])

    # xticks setup 
    ax0.set_xticks(np.arange(0.05, 0.4, 0.05))
    ax1.set_xticks(np.arange(0.05, 0.4, 0.05))
    ax2.set_xticks(np.arange(0.05, 0.4, 0.05))
    ax3.set_xticks(np.arange(0.05, 0.4, 0.05))
    ax4.set_xticks(np.arange(0.05, 0.4, 0.05))
    ax5.set_xticks(np.arange(0.05, 0.4, 0.05))
    ax6.set_xticks(np.arange(0.05, 0.4, 0.05))
    ax7.set_xticks(np.arange(0.05, 0.4, 0.05))
    
    fig0.tight_layout()

    plt.show()