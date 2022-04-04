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
        plt.title('Lieblein efficiency -- axial compressor')

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

        plt.figure(figsize=(8,8))
        plt.title('Lieblein efficiency -- axial compressor')

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

def stagePerf(phi=0, psi=0, perc=0, rD=0.5, phiVec=np.linspace(0,1.1,1000), plot=True, save=False, position='stagePerf.pgf'):
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

    def phiLim(phi, error=5e-3):
        '''
        This function computes the phi value for which psiLim == 1
        '''
        for _,ii in enumerate(phi):
            if np.abs(( psiLim(ii) - 1 ) / psiLim(ii)) < error:
                return ii 

    # values computation
    psiLIM = psiLim(phiVec)    # work coefficient limit vs flow coefficient 
    psiBETA = psiBeta(phiVec)  # work coefficient vs flow coefficient for beta = 70deg 

    # computing psi with respect to a percentage and phi 
    if perc != 0 and phi != 0 and psi == 0: 
        psi = perc * psiLim(phi)

    if plot:
        # plotting
        plt.figure()
        plt.plot(phiVec, psiLIM, 'k', label=r'$\psi_{Lim}$')
        plt.plot(phiVec, psiBETA, 'r', label=r'$\beta > 70^{\circ}$')
        plt.plot([phiLim(phiVec),phiLim(phiVec)] , [0, np.max(psiLIM)], 'b', label=r'$\phi_{Lim}$')
        if psi != 0 and phi !=0:
            plt.plot(phi, psi, 'og', label=r'$[\phi = {0:.3f}, \psi = {1:.3f}]$'.format(phi, psi))
        plt.xlim(0,np.max(phiVec))
        plt.ylim(0,1)
        plt.xlabel(r'$\phi$')
        plt.ylabel(r'$\psi$')
        plt.title(r'$\frac{W_2}{W_1} > 0.7$' + '    ' + r'$\chi = {0:.2f}$'.format(rD))
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
        plt.figure()
        plt.plot(phiVec, psiLIM, 'k', label=r'$\psi_{Lim}$')
        plt.plot(phiVec, psiBETA, 'r', label=r'$\beta > 70^{\circ}$')
        plt.plot([phiLim(phiVec),phiLim(phiVec)] , [0, np.max(psiLIM)], 'b', label=r'$\phi_{Lim}$')
        if psi != 0 and phi !=0:
            plt.plot(phi, psi, 'og', label=r'$[\phi = {0:.3f}, \psi = {1:.3f}]$'.format(phi, psi))
        plt.xlim(0,np.max(phiVec))
        plt.ylim(0,1)
        plt.xlabel(r'$\phi$')
        plt.ylabel(r'$\psi$')
        plt.title(r'$\frac{W_2}{W_1} > 0.7$' + '    ' + r'$\chi = {0:.2f}$'.format(rD))
        plt.grid(linestyle='--')
        plt.legend(loc='upper left', bbox_to_anchor=(0, 1))
        plt.tight_layout()
        # figure saving 
        plt.savefig(position, bbox_inches='tight')

    if perc != 0 and phi != 0 and psi !=0:
        return psi 

