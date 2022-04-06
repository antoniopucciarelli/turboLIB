# TURBOMACHINERY -- LIEBLEIN BLADE MODELING APPROACH
# AUTHOR: antonio pucciarelli 
#
# PROGRAM DESCRIPTION
#   CONTENT: engineering coefficient functions related to Lieblein blade modeling approach
#  

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

    if plot:
        fig = plt.figure(figsize=(8,8))
        plt.plot(tbcVec, Ktifunc(tbcVec), 'k', linewidth=2)
        plt.grid(linestyle='--')
        plt.title(r'$K_{{t,i}}$')
        plt.xlabel(r'$\frac{t_b}{c}$')
        plt.ylabel(r'$K_{{t,i}}$')
        if tbc != 0:
            Kti = Ktifunc(tbc)
            plt.plot(tbc, Kti, marker='o', markersize=8, markeredgecolor='k', color='g', markeredgewidth=1.5)
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

        fig = plt.figure(figsize=(8,8))
        plt.plot(tbcVec, Ktifunc(tbcVec), 'k', linewidth=2)
        plt.grid(linestyle='--')
        plt.title(r'$K_{{t,i}}$')
        plt.xlabel(r'$\frac{t_b}{c}$')
        plt.ylabel(r'$K_{{t,i}}$')
        if tbc != 0:
            Kti = Ktifunc(tbc)
            plt.plot(tbc, Kti, marker='o', markersize=8, markeredgecolor='k', markercolor='g', markerwidth=1.5)
        plt.savefig(position)

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
        Computation of i0*. This function uses a n parameter as exponent. 
        '''
        
        # exponent computation 
        p = 0.914 + solidity**3 / 160

        # i0* computation 
        i0 = beta1**p / (5 + 46 * np.exp(-2.3*solidity)) - 0.1 * solidity**3 * np.exp((beta1 - 70)/4)

        return i0 

    # vector allocation 
    beta1Vec = np.linspace(0, 70, 1000)
    solidityVec = np.linspace(0.4, 2, 5)

    if plot: 
        fig = plt.figure(figsize=(8,8))
        for _,sol in enumerate(solidityVec):
            plt.plot(beta1Vec, i0func(beta1Vec, sol), linewidth=0.5, label=r'$\sigma = {0:.2f}$'.format(sol))
        if beta1 != 0 and solidity != 0:
            plt.plot(beta1Vec, i0func(beta1Vec, solidity), 'k', linewidth=2, label=r'$\sigma = {0:.2f}$'.format(solidity))
            plt.plot(beta1, i0func(beta1, solidity), linestyle='', marker='o', markersize=8, markeredgecolor='k', color='g', markeredgewidth=1.5, label=r'$\sigma = {0:.2f}, i_{{0, 10}}^{{*}} = {1:.2f}$'.format(solidity, i0func(beta1, solidity)))
        plt.xlabel(r'$\beta_1$')
        plt.ylabel(r'$i_{{0,10}}^{*}$')
        plt.legend(loc='upper left')
        plt.grid(linestyle='--')
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
        fig = plt.figure(figsize=(8,8))
        for _,sol in enumerate(solidityVec):
            plt.plot(beta1Vec, i0func(beta1Vec, sol), linewidth=0.5, label=r'$\sigma = {0:.2f}$'.format(sol))
        if beta1 != 0 and solidity != 0:
            plt.plot(beta1Vec, i0func(beta1Vec, solidity), 'k', linewidth=2, label=r'$\sigma = {0:.2f}$'.format(solidity))
            plt.plot(beta1, i0func(beta1, solidity), linestyle='', marker='o', markersize=8, markeredgecolor='k', color='g', markeredgewidth=1.5, label=r'$\sigma = {0:.2f}, i_{{0, 10}}^{{*}} = {1:.2f}$'.format(solidity, i0func(beta1, solidity)))
        plt.xlabel(r'$\beta_1$')
        plt.ylabel(r'$i_{{0,10}}^{*}$')
        plt.grid(linestyle='--')
        plt.savefig(position)
    
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
        Computation of n. This function uses a n parameter as exponent. 
        '''

        n = 0.025 * solidity - 0.06 - (beta1 / 90)**(1 + 1.2 * solidity) / (1.5 + 0.43 * solidity)

        return n 

    # vector allocation 
    beta1Vec = np.linspace(0, 70, 1000)
    solidityVec = np.linspace(0.4, 2, 5)

    if plot: 
        fig = plt.figure(figsize=(8,8))
        for _,sol in enumerate(solidityVec):
            plt.plot(beta1Vec, nfunc(beta1Vec, sol), linewidth=0.5, label=r'$\sigma = {0:.2f}$'.format(sol))
        if beta1 != 0 and solidity != 0:
            plt.plot(beta1Vec, nfunc(beta1Vec, solidity), 'k', linewidth=2, label=r'$\sigma = {0:.2f}$'.format(solidity))
            plt.plot(beta1, nfunc(beta1, solidity), linestyle='', marker='o', markersize=8, markeredgecolor='k', color='g', markeredgewidth=1.5, label=r'$\sigma = {0:.2f}, n = {1:.2f}$'.format(solidity, nfunc(beta1, solidity)))
        plt.xlabel(r'$\beta_1$')
        plt.ylabel(r'$n$')
        plt.legend(loc='lower left')
        plt.grid(linestyle='--')
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
        fig = plt.figure(figsize=(8,8))
        for _,sol in enumerate(solidityVec):
            plt.plot(beta1Vec, nfunc(beta1Vec, sol), linewidth=0.5, label=r'$\sigma = {0:.2f}$'.format(sol))
        if beta1 != 0 and solidity != 0:
            plt.plot(beta1Vec, nfunc(beta1Vec, solidity), 'k', linewidth=2, label=r'$\sigma = {0:.2f}$'.format(solidity))
            plt.plot(beta1, nfunc(beta1, solidity), linestyle='', marker='o', markersize=8, markeredgecolor='k', color='g', markeredgewidth=1.5, label=r'$\sigma = {0:.2f}, n = {1:.2f}$'.format(solidity, nfunc(beta1, solidity)))
        plt.xlabel(r'$\beta_1$')
        plt.ylabel(r'$n$')
        plt.legend(loc='lower left')
        plt.grid(linestyle='--')
        plt.savefig(position)
    
    if beta1 != 0 and solidity != 0: 
        return nfunc(beta1, solidity)



