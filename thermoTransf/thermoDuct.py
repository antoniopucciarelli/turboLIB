# TURBOMACHINERY -- LIBRARY FOR THE PLOTTING OF THERMODYNAMIC PROCESSES
# AUTHOR: antonio pucciarelli 
#
# PROGRAM DESCRIPTION
#   NOZZLE PROCESSES: plot and computation 
#       -- nozzle shape 
#       -- nozzle thermodynamics
#       -- nozzle properties used in space propulsion 
#
#   FLUID: air
#  

# importing libraries
import matplotlib.pyplot as plt 
import numpy as np 
import pyromat as pm 
from scipy.optimize import bisect 

# setting up dimensions
unit_pressure = pm.config['unit_pressure'] = 'bar'
unit_temperature = pm.config['unit_temperature'] = 'K'
unit_energy = pm.config['unit_energy'] = 'kJ'
unit_mass = pm.config['unit_mass'] = 'kg'

# air object generation
air = pm.get('ig.air')

def nozzle(Ain=0,At=0,Aout=0,xIn=0,xT=0,xOut=0,dim1=100,dim2=100,color='k',name='nozzle',plot=True):
    '''
    Nozzle shape generation:
        -- linear variation of the nozzle shape 

        function inputs:
            Ain     -- nozzle inlet area 
            At      -- nozzle throat area 
            Aout    -- nozzle exit area
            dim1    -- number of points for convergent description 
                    -- default set to 100
            dim2    -- number of points for divergent description 
                    -- default set to 100 
    '''

    # convergent array description 
    convVec = np.linspace(Ain,At,dim1)
    xConvVec = np.linspace(xIn,xT,dim1)

    # divergent array description 
    divVec = np.linspace(At,Aout,dim2)
    xDivVec = np.linspace(xT,xOut,dim2)

    # complete discretized nozzle 
    nozzleVec = np.append(convVec,divVec)
    xVec = np.append(xConvVec,xDivVec)

    # plotting nozzle 
    if plot:
        plt.plot(xVec, nozzleVec, color=color, label=name)

    return nozzleVec, convVec, divVec, xVec, xConvVec, xDivVec

def nozzleThermo(nozzleVec=[],xVec=[],At=0,xT=0,Tin=0,pIn=0,plot=True):
    '''
    Nozzle thermodynamics computation
        -- thermodynamic quantitties computation for a defined geometry nozzle without friction
    
        function inputs:
            nozzleVec   -- nozzle area geometry
            xVec        -- nozzle longitudinal dimension
            At          -- nozzle throat area 
            xT          -- throat position 
    '''

    # importing air object 
    global air 

    # importing dimensions
    global unit_pressure, unit_temperature 

    # air heat capacity ratio
    gamma = air.gam()

    # nozzleVec dimension 
    dim = len(nozzleVec)

    # area vs Mach thermodynamics -- isentropic transformation 
    def areaMach(M,area,At,gamma):
        '''
        Isentropic transformation 
            -- Mach number vs area

            function inputs:
                M       -- nozzle Mach number at a defined area 
                area    -- nozzle area design
                At      -- nozzle throat area 
                gamma   -- heat capacity ratio
        '''

        return area/At - 1/M * (2/(gamma+1)*(1+(gamma-1)/2*M**2))**((gamma+1)/(2*(gamma-1)))

    # Mach number computation with respect to the defined geometry
    #   -- using numerical approach for the computation of Mach number 
    Mvec = np.zeros(dim)

    for ii, A in enumerate(nozzleVec):
        # bisection process for finding the 0
        # studying for convergent and divergent part wrt x position along the nozzle 
        if xVec[ii] < xT:
            # convergent part 
            extr1 = 1e-10
            extr2 = 1
            Mvec[ii] = bisect(f = areaMach, a = extr1, b = extr2, args=(A,At,gamma))
        else: 
            # divergent part 
            extr1 = 1 + 1e-10
            extr2 = 10
            Mvec[ii] = bisect(f = areaMach, a = extr1, b = extr2, args=(A,At,gamma))

    # temperature and pressure computation 
    # vector allocation 
    Tvec = Tin / (1 + (gamma-1)/2 * Mvec**2);                    
    pVec = pIn / (1 + (gamma-1)/2 * Mvec**2)**(gamma/(gamma-1)); 

    # plotting Mach, temperature and pressure lines
    if plot:
        # figure generation
        fig, ax = plt.subplots(3)
        fig.suptitle('Nozzle analysis')
        # Mach plot
        ax[0].plot(xVec,Mvec,'k')
        ax[0].set_title('Mach')
        ax[0].set_xlabel('x')
        ax[0].set_ylabel('M')
        # temperature plot 
        ax[1].plot(xVec,Tvec,'k')
        ax[1].set_title('Temperature')
        ax[1].set_xlabel('x')
        ax[1].set_ylabel('T [{0:s}]'.format(unit_temperature))
        ax[2].plot(xVec,pVec,'k')
        # pressure plot 
        ax[2].set_title('Pressure')
        ax[2].set_xlabel('x')
        ax[2].set_ylabel('p [{0:s}]'.format(unit_pressure))
        plt.tight_layout()

    return Mvec, Tvec, pVec
