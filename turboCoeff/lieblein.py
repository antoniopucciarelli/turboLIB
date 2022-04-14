# TURBOMACHINERY -- LIEBLEIN BLADE MODELING APPROACH
# AUTHOR: antonio pucciarelli 
#
# PROGRAM DESCRIPTION
#   CONTENT: engineering coefficient functions related to Lieblein blade modeling approach
#  

# importing libraries 
import numpy as np 
import matplotlib.pyplot as plt 

def KtiFunc(tbc=0, tb=0, c=0, plot=False, save=False, position='Kti.pgf'):
    '''
    Shape corrector function from Johnsen and Bullock. 
        inputs:
            tbc         -- tb/c thickness chord ratio
            tb          -- profile thickness 
            c           -- profile chord
            plot        -- boolean value for the data plot 
            save        -- boolean value for the saving on the plot in vectorial format 
            position    -- path for the figure saving  
    '''

    # vector generation 
    tbcVec = np.linspace(0, 0.15, 1000)

    def Ktifunc(tbc):
        '''
        Kti function. It depends on an additional function q(tb/c).
        '''

        # exponent value computation
        q = 0.28 / (0.1 + (tbc)**0.3)

        return (10 * tbc)**q

    # check on inputs
    if tb != 0 and c != 0:
        tbc = tb/c

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
        fig = plt.figure(figsize=(8,8))
        plt.plot(tbcVec, Ktifunc(tbcVec), 'k', linewidth=2)
        plt.grid(linestyle='--')
        plt.title(r'$K_{{t,i}}$')
        plt.xlabel(r'$\frac{t_b}{c}$')
        plt.ylabel(r'$K_{{t,i}}$')
        if tbc != 0:
            Kti = Ktifunc(tbc)
            plt.plot(tbc, Kti, marker='o', markersize=8, markeredgecolor='k', color='g', markeredgewidth=1.5)
        if save:
            fig.savefit(position)
        else:
            plt.show()

    if tbc !=0:
        return Ktifunc(tbc)

def i0Func(beta1=0, solidity=0, plot=False, save=False, position='i0.pgf'):
    '''
    This function computes the i0* value used in the Lieblein modeling approach. 
        inputs:
            beta1       -- inlet relative flow angle 
            solidity    -- c / s -> chord / pitch 
            plot        -- boolean value for the plotting of the figure 
            save        -- boolean value for the saving of the figure 
            position    -- path where the figure is saved  
    '''

    def i0func(beta1, solidity):
        '''
        Computation of i0,10*.
        '''
        
        # exponent computation 
        p = 0.914 + solidity**3 / 160

        # i0* computation 
        i0 = beta1**p / (5 + 46 * np.exp(-2.3*solidity)) - 0.1 * solidity**3 * np.exp((beta1 - 70)/4)

        if i0 == np.Inf:
            print('not ok')
            print(i0)
            print(p)
            print(beta1)
            print(solidity)

        return i0 

    # vector allocation 
    beta1Vec = np.linspace(0, 70, 1000)
    solidityVec = np.linspace(0.4, 2, 5)

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
        fig = plt.figure(figsize=(8,8))
        for _,sol in enumerate(solidityVec):
            plt.plot(beta1Vec, i0func(beta1Vec, sol), linewidth=0.5, label=r'$\sigma = {0:.2f}$'.format(sol))
        if beta1 != 0 and solidity != 0:
            plt.plot(beta1Vec, i0func(beta1Vec, solidity), 'k', linewidth=2, label=r'$\sigma = {0:.2f}$'.format(solidity))
            plt.plot(beta1, i0func(beta1, solidity), linestyle='', marker='o', markersize=8, markeredgecolor='k', color='g', markeredgewidth=1.5, label=r'$\beta_1 = {0:.2f}, i_{{0, 10}}^{{*}} = {1:.2f}$'.format(beta1, i0func(beta1, solidity)))
        plt.xlabel(r'$\beta_1$')
        plt.ylabel(r'$i_{{0,10}}^{*}$')
        plt.legend(loc='upper left')
        plt.grid(linestyle='--')
        if save:
            fig.savefig()
        else:
            plt.show()
    
    if beta1 != 0 and solidity != 0: 
        return i0func(beta1, solidity)

def nFunc(beta1=0, solidity=0, plot=False, save=False, position='n.pgf'):
    '''
    This function computes the n value used in the Lieblein modeling approach. 
        inputs:
            beta1       -- inlet relative flow angle 
            solidity    -- c / s -> chord / pitch 
            plot        -- boolean value for the plotting of the figure 
            save        -- boolean value for the saving of the figure 
            position    -- path where the figure is saved  
    '''

    def nfunc(beta1, solidity):
        '''
        Computation of n. 
        '''

        n = 0.025 * solidity - 0.06 - (beta1 / 90)**(1 + 1.2 * solidity) / (1.5 + 0.43 * solidity)

        return n 

    # vector allocation 
    beta1Vec = np.linspace(0, 70, 1000)
    solidityVec = np.linspace(0.4, 2, 5)

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
        fig = plt.figure(figsize=(8,8))
        for _,sol in enumerate(solidityVec):
            plt.plot(beta1Vec, nfunc(beta1Vec, sol), linewidth=0.5, label=r'$\sigma = {0:.2f}$'.format(sol))
        if beta1 != 0 and solidity != 0:
            plt.plot(beta1Vec, nfunc(beta1Vec, solidity), 'k', linewidth=2, label=r'$\sigma = {0:.2f}$'.format(solidity))
            plt.plot(beta1, nfunc(beta1, solidity), linestyle='', marker='o', markersize=8, markeredgecolor='k', color='g', markeredgewidth=1.5, label=r'$\beta_1 = {0:.2f}, n = {1:.2f}$'.format(beta1, nfunc(beta1, solidity)))
        plt.xlabel(r'$\beta_1$')
        plt.ylabel(r'$n$')
        plt.legend(loc='lower left')
        plt.grid(linestyle='--')
        if save:
            fig.savefig(position)
        else:
            plt.show()
    
    if beta1 != 0 and solidity != 0: 
        return nfunc(beta1, solidity)

def delta0Func(beta1=0, solidity=0, plot=False, save=False, position='delta0.pgf'):
    '''
    This function computes the delta0 value used in the Lieblein modeling approach. 
        inputs:
            beta1       -- inlet relative flow angle 
            solidity    -- c / s -> chord / pitch 
            plot        -- boolean value for the plotting of the figure 
            save        -- boolean value for the saving of the figure 
            position    -- path where the figure is saved  
    '''

    def delta0func(beta1, solidity):
        '''
        Computation of delta0,10*.  
        '''

        delta0 = 0.01 * solidity * beta1 + (0.74 * solidity**1.9 + 3 * solidity) * (beta1 / 90 )**(1.67 + 1.09 * solidity)

        return delta0 

    # vector allocation 
    beta1Vec = np.linspace(0, 70, 1000)
    solidityVec = np.linspace(0.4, 2, 5)

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
    
        fig = plt.figure(figsize=(8,8))
        for _,sol in enumerate(solidityVec):
            plt.plot(beta1Vec, delta0func(beta1Vec, sol), linewidth=0.5, label=r'$\sigma = {0:.2f}$'.format(sol))
        if beta1 != 0 and solidity != 0:
            plt.plot(beta1Vec, delta0func(beta1Vec, solidity), 'k', linewidth=2, label=r'$\sigma = {0:.2f}$'.format(solidity))
            plt.plot(beta1, delta0func(beta1, solidity), linestyle='', marker='o', markersize=8, markeredgecolor='k', color='g', markeredgewidth=1.5, label=r'$\beta_1 = {0:.2f}, \delta_{{0, 10}}^{{*}} = {1:.2f}$'.format(beta1, delta0func(beta1, solidity)))
        plt.xlabel(r'$\beta_1$')
        plt.ylabel(r'$n$')
        plt.legend(loc='upper left')
        plt.grid(linestyle='--')
        if save:
            fig.savefig(position)
        else:
            plt.show()
    
    if beta1 != 0 and solidity != 0: 
        return delta0func(beta1, solidity)    

def mFunc(beta1=0, solidity=0, plot=False, save=False, position='m.pgf'):
    '''
    This function computes the m value used in the Lieblein modeling approach. 
        inputs:
            beta1       -- inlet relative flow angle 
            solidity    -- c / s -> chord / pitch 
            plot        -- boolean value for the plotting of the figure 
            save        -- boolean value for the saving of the figure 
            position    -- path where the figure is saved  
    '''

    def mfunc(beta1, solidity):
        '''
        Computation of m. 
            !!! till now this function work only for NACA-65 airfoils !!!
        '''

        # exponent computation 
        x = beta1 / 100
        b = 0.9625 - 0.17 * x - 0.85 * x**3

        # for NACA-65 airfoils 
        m10 = 0.17 - 0.0333 * x + 0.333 * x**2
        
        # m computation 
        m = m10 / solidity**b 

        return m
    
    # vector allocation 
    beta1Vec = np.linspace(0, 70, 1000)
    solidityVec = np.linspace(0.4, 2, 5)

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

        fig = plt.figure(figsize=(8,8))
        for _,sol in enumerate(solidityVec):
            plt.plot(beta1Vec, mfunc(beta1Vec, sol), linewidth=0.5, label=r'$\sigma = {0:.2f}$'.format(sol))
        if beta1 != 0 and solidity != 0:
            plt.plot(beta1Vec, mfunc(beta1Vec, solidity), 'k', linewidth=2, label=r'$\sigma = {0:.2f}$'.format(solidity))
            plt.plot(beta1, mfunc(beta1, solidity), linestyle='', marker='o', markersize=8, markeredgecolor='k', color='g', markeredgewidth=1.5, label=r'$\beta_1 = {0:.2f}, n = {1:.2f}$'.format(beta1, mfunc(beta1, solidity)))
        plt.xlabel(r'$\beta_1$')
        plt.ylabel(r'$m$')
        plt.legend(loc='best')
        plt.grid(linestyle='--')

        if save:
            fig.savefig(position)
        else:
            plt.show()        
    
    if beta1 != 0 and solidity != 0: 
        return mfunc(beta1, solidity)

def KtdeltaFunc(tbc=0, tb=0, c=0, plot=False, save=False, position='Kti.pgf'):
    '''
    Shape corrector function from Johnsen and Bullock. 
        inputs:
            tbc         -- tb/c thickness chord ratio
            tb          -- profile thickness 
            c           -- profile chord
            plot        -- boolean value for the data plot 
            save        -- boolean value for the saving on the plot in vectorial format 
            position    -- path for the figure saving  
    '''

    # vector generation 
    tbcVec = np.linspace(0, 0.15, 1000)

    def Ktdeltafunc(tbc):
        '''
        Ktdelta function. 
        '''

        Ktdelta = 6.25 * tbc + 37.5 * tbc**2

        return Ktdelta

    # check on inputs
    if tb != 0 and c != 0:
        tbc = tb/c

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

        fig = plt.figure(figsize=(8,8))
        plt.plot(tbcVec, Ktdeltafunc(tbcVec), 'k', linewidth=2)
        plt.grid(linestyle='--')
        plt.title(r'$K_{{t,i}}$')
        plt.xlabel(r'$\frac{t_b}{c}$')
        plt.ylabel(r'$K_{{t,\delta}}$')
        if tbc != 0:
            Kti = Ktdeltafunc(tbc)
            plt.plot(tbc, Kti, marker='o', markersize=8, markeredgecolor='k', color='g', markeredgewidth=1.5)
        
        if save:
            fig.savefig(position)
        else:
            plt.show()

    if tbc !=0:
        return Ktdeltafunc(tbc)

def oFunc(gamma=0, tbc=0.1, solidity=0, pitch=0, Cl=0, plot=False, save=False, position='o.pgf'):
    '''
    This function computes the minimum opening between 2 blades in a cascade row for NACA-65 blade type. 
        inputs:
            gamma       -- stagger angle 
            tbc         -- tb / c 
            solidity    -- section solidity 
            pitch       -- distance between 2 consecutive blades
            Cl          -- lift coefficient of the blade at 0 AOA 
            plot        -- boolean value for the plotting of the results
            save        -- boolean value for the saving of the figure in .pgf format 
            position    -- path where to save the figure 
    '''

    # value used for the opening computation
    def phifunc(gamma, Cl):
        '''
        This function computes phi.
        '''
        
        # converting data to positive values
        gamma = np.abs(gamma)
        Cl = np.abs(Cl)

        return gamma * (1.0 - 0.05 * Cl**1.5) + 5 * Cl**1.5 - 2.0

    # opening computation 
    def ofunc(pitch, solidity, phi, tbc):
        '''
        This function computes the opening, o.
        '''

        return pitch * ((1 - tbc * np.sqrt(solidity)) * np.cos(np.deg2rad(phi)))**(np.sqrt(solidity))

    if Cl != 0:
        o = ofunc(pitch, solidity, phifunc(gamma, Cl), tbc)

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

        gammaVec = np.arange(10, 70, 10)
        ClVec = np.linspace(0.4, 1.4, 300)

        fig = plt.figure(figsize=(9,9))
        for ii in range(len(gammaVec)):
            plt.plot(ClVec, ofunc(pitch, solidity, phifunc(gammaVec[ii], ClVec), tbc)/pitch, linewidth=1, label=r'$\gamma = {0:.2f}$'.format(gammaVec[ii]))
        
        if Cl != 0:
            plt.plot(ClVec, ofunc(pitch, solidity, phifunc(gamma, ClVec), tbc)/pitch, 'k', linewidth=2, label=r'$\gamma = {0:.2f}$'.format(gamma))
            plt.plot(Cl, o/pitch, linestyle='', marker='o', markersize=8, color='g', markeredgecolor='k', markeredgewidth=1.5, label=r'$Cl_0 = {0:.2f}, \gamma = {1:.2f}$'.format(Cl, gamma))

        plt.xlabel(r'$Cl_0$')
        plt.ylabel(r'$\frac{{o}}{{s}}$')
        plt.grid(linestyle='--')
        plt.legend()

        if save:
            fig.savefig(position)
        else:
            plt.show()        
    
    if Cl != 0:
        return o   