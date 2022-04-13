# TURBOMACHINERY -- LIBRARY FOR THE INITIAL TURBOMACHINERY DESIGN
# AUTHOR: antonio pucciarelli 
#
# PROGRAM DESCRIPTION
#   TURBOMACHINERY DESIGN CLASS:
#       this script sets the rotor/stator class for the turbomachinery desing analisys 
#       

# importing libraries
import numpy as np 
import matplotlib.pyplot as plt 
from geometry import bladeGenerator
import warnings

class section:
    '''
    Blade object, it is used in the stage object.
        AIM:
            --- blade geometry description 
            --- blade thermodynamics 
    '''

    def __init__(self, midpoint, bottom, tip, height):
        '''
        Rotor object declaration: 
            variables:
                midpoint    -- section midpoint
                bottom      -- section lowest point 
                tip         -- section highest point 
                height      -- section height
                velocity    -- section speed 
        '''

        # position 
        self.midpoint = midpoint 
        self.bottom   = bottom 
        self.tip      = tip 
        self.height   = height 

    def allocateDynamics(self, Va, Vt, U):
        '''
        This function allocates the velocity vector for the section object.
            inputs:
                Va  -- axial flow 
                Vt  -- tangential flow 
                U   -- rotation speed 
        '''

        # rotation speed 
        self.U  = U 
        
        # absolute quantities
        self.Va    = Va 
        self.Vt    = Vt
        self.V     = np.sqrt(Va**2 + Vt**2)
        self.alpha = np.rad2deg(np.arctan(Vt/Va))
        
        # relative quantities 
        self.Wa   = Va 
        self.Wt   = Vt - U 
        self.W    = np.sqrt(self.Wa**2 + self.Wt**2)
        self.beta = np.rad2deg(np.arctan(self.Wt/self.Wa))

    def allocateThermodynamics(self, Tt, Pt, T, P, Ttr, Ptr, rho, rhot, s, R=287.06, gamma=1.4):
        '''
        This function computes the thermodynamic properties of a specific point in the blade.
            inputs:
                Tt  -- total temperature
                Pt  -- total pressure 
        '''

        # cP computation
        cP = gamma / (gamma - 1) * R

        # thermodynamic properties computation
        self.Tt = Tt
        self.Ttr = Ttr
        self.T  = T 
        self.Tr = Tt - self.W**2/(2 * cP)
        self.a  = np.sqrt(gamma * R * self.T)
        self.M  = self.V / self.a  
        self.Mr = self.W / self.a 
        self.P  = P 
        self.Pt = Pt 
        self.Ptr = Ptr
        self.rho = rho 
        self.rhot = rhot 
        self.s = s 

    def mFlux(self):
        '''
        This function computes the mass flux in the stream tube.
            -- the axial fluid speed and the density are considered constant in the stream tube.
        '''

        massFlux = np.pi * (self.tip**2 - self.bottom**2) * self.rho * self.Va

        return massFlux 