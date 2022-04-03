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

def efficiency(filePath, phi=0, R=0, plot=False, save=False, position='efficiency.pgf'):
    '''
    This function describes the adimensional parameters for the turbomachinery design.
    '''

    # importing data 
    axialEfficiency092 = np.float32(np.loadtxt(filePath + 'axialEfficiency092.txt', comments='#', delimiter=','))
    axialEfficiency091 = np.float32(np.loadtxt(filePath + 'axialEfficiency091.txt', comments='#', delimiter=','))
    axialEfficiency090 = np.float32(np.loadtxt(filePath + 'axialEfficiency090.txt', comments='#', delimiter=','))
    axialEfficiency088 = np.float32(np.loadtxt(filePath + 'axialEfficiency088.txt', comments='#', delimiter=','))
    axialEfficiency086 = np.float32(np.loadtxt(filePath + 'axialEfficiency086.txt', comments='#', delimiter=','))

    # plotting efficiency
    if plot:
        plt.figure(figsize=(8,8))
        plt.title('Lieblein efficiency -- axial compressor')
        plt.plot(0.48, 0.5, 'ok', label=r'$\eta = 0.926$')
        plt.plot(axialEfficiency092[:,0],axialEfficiency092[:,1],'b*',linewidth=3,label=r'$\eta = 0.92$')
        plt.plot(axialEfficiency091[:,0],axialEfficiency091[:,1],'r*',linewidth=3,label=r'$\eta = 0.91$')
        plt.plot(axialEfficiency090[:,0],axialEfficiency090[:,1],'m*',linewidth=3,label=r'$\eta = 0.90$')
        plt.plot(axialEfficiency088[:,0],axialEfficiency088[:,1],'c*',linewidth=3,label=r'$\eta = 0.88$')
        plt.plot(axialEfficiency086[:,0],axialEfficiency086[:,1],'y*',linewidth=3,label=r'$\eta = 0.86$')
        if R != 0 and phi != 0:
            plt.plot(phi, R, linestyle='', marker='s', color='g', markersize=7, label=r'$[\phi = {0:.2f}, \chi = {1:.2f}]$'.format(phi, R))
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
        plt.plot(0.48, 0.5, 'ok', label=r'$\eta = 0.926$')
        plt.plot(axialEfficiency092[:,0],axialEfficiency092[:,1],'b*',linewidth=3,label=r'$\eta = 0.92$')
        plt.plot(axialEfficiency091[:,0],axialEfficiency091[:,1],'r*',linewidth=3,label=r'$\eta = 0.91$')
        plt.plot(axialEfficiency090[:,0],axialEfficiency090[:,1],'m*',linewidth=3,label=r'$\eta = 0.90$')
        plt.plot(axialEfficiency088[:,0],axialEfficiency088[:,1],'c*',linewidth=3,label=r'$\eta = 0.88$')
        plt.plot(axialEfficiency086[:,0],axialEfficiency086[:,1],'y*',linewidth=3,label=r'$\eta = 0.86$')
        if R != 0 and phi != 0:
            plt.plot(phi, R, linestyle='', marker='s', color='g', markersize=7, label=r'$[\phi = {0:.2f}, \chi = {1:.2f}]$'.format(phi, R))
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

def stagePerf(phi=0, psi=0, R=0.5, phiVec=np.linspace(0,1.1,1000), plot=True, save=False, position='stagePerf.pgf'):
    '''
    This function plots the phi and psi parameter with respect to ASME axial compressor design constraints
        The main constraints for the study of the adimensional parameters are:
            W2/W1 < 0.7     -- to avoid flow separation 
            psiLim = 1      -- for having enough margin from the surge condition 
                            -- ASME -> figure 10.6 psiLim determines line slopes that has to be < 1 in order to avoid surge 
            beta1 < 70 deg  -- for avoiding reduction in cascade performances 
            curve symmetric with respect to reaction degree = 0.5
        inputs:
            phi         -- flow coefficient for the stage 
            psi         -- work coefficient for the stage
            R           -- stage reaction degree 
            phiVec      -- phi vector for the psiLim representation 
            plot        -- boolean value for plotting the chart 
            save        -- boolean value for plotting the chart in vectorial format 
            position    -- path where the pgf file should be saved 
    '''
    
    # function generation eqn. from ASME 10.28 - 10.30
    Rcap = 0.5 + np.abs(R - 0.5)
    psiLim = lambda phi: 6 * Rcap / 17 + 0.85 * (0.5/Rcap)**1.18 * phi**(2 + 0.1/Rcap)
    psiBeta = lambda phi: 5 * phi - 2 * R

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

    if plot:
        # plotting
        plt.figure()
        plt.plot(phiVec, psiLIM, 'k', label=r'$\psi_{Lim}$')
        plt.plot(phiVec, psiBETA, 'r', label=r'$\beta > 70^{\circ}$')
        plt.plot([phiLim(phiVec),phiLim(phiVec)] , [0, np.max(psiLIM)], 'b', label=r'$\psi_{Lim} = 1$')
        if psi != 0 and phi !=0:
            plt.plot(phi, psi, 'og', label=r'$[\phi = {0:.2f}, \psi = {1:.2f}]$'.format(phi, psi))
        plt.xlim(0,np.max(phiVec))
        plt.ylim(0,1)
        plt.xlabel(r'$\phi$')
        plt.ylabel(r'$\psi$')
        plt.title(r'$\frac{W_2}{W_1} < 0.7$' + '    ' + r'$\chi = {0:.2f}$'.format(R))
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
        plt.plot([phiLim(phiVec),phiLim(phiVec)] , [0, np.max(psiLIM)], 'b', label=r'$\psi_{Lim} = 1$')
        if psi != 0 and phi !=0:
            plt.plot(phi, psi, 'og', label=r'$[\phi = {0:.2f}, \psi = {1:.2f}]$'.format(phi, psi))
        plt.xlim(0,np.max(phiVec))
        plt.ylim(0,1)
        plt.xlabel(r'$\phi$')
        plt.ylabel(r'$\psi$')
        plt.title(r'$\frac{W_2}{W_1} < 0.7$' + '    ' + r'$\chi = {0:.2f}$'.format(R))
        plt.grid(linestyle='--')
        plt.legend(loc='upper left', bbox_to_anchor=(0, 1))
        plt.tight_layout()
        # figure saving 
        plt.savefig(position, bbox_inches='tight')

