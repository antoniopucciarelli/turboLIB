# TURBOMACHINERY -- LIBRARY FOR THE INITIAL TURBOMACHINERY DESIGN
# AUTHOR: antonio pucciarelli 
#
# PROGRAM DESCRIPTION
#   TURBOMACHINERY DESIGN CLASS:
#       this script sets the rotor/stator class for the turbomachinery desing analisys 
#       

# importing libraries
import numpy as np 

class section:
    '''
    Blade object, it is used in the stage object.
        AIM:
            --- blade geometry description 
            --- blade thermodynamics 
    '''

    def __init__(self, midpoint, bottom, tip, height, pitch, solidity=1, tbc=0.1):
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
        self.pitch    = pitch
        self.solidity = solidity
        self.tbc      = tbc
        self.Cl       = 0.0
        self.s        = 0.0
        self.rD       = 0.0
        self.theta    = 0
        self.i        = 0

    def allocateKinetics(self, Va, Vt, U):
        '''
        This function allocates the velocity vector for the section object.
            inputs:
                Va  -- axial flow 
                Vt  -- tangential flow 
                U   -- rotation speed 
        '''

        # rotation speed 
        self.U = U 
        
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

    def allocateThermodynamics(self, Tt, Pt, T, P, Ttr, Ptr, rho, rhot, rhotr, s, R=287.06, gamma=1.4):
        '''
        This function computes the thermodynamic properties of a specific point in the blade.
            inputs:
                Tt  -- total temperature
                Pt  -- total pressure 
        '''

        # cP computation
        cP = gamma / (gamma - 1) * R

        # thermodynamic properties computation
        self.s     = s 
        self.T     = T 
        self.Tr    = Tt - self.W**2/(2 * cP)
        self.Tt    = Tt
        self.Ttr   = Ttr
        self.a     = np.sqrt(gamma * R * self.T)
        self.M     = self.V / self.a  
        self.Mr    = self.W / self.a 
        self.P     = P 
        self.Pt    = Pt 
        self.Ptr   = Ptr
        self.rho   = rho 
        self.rhot  = rhot 
        self.rhotr = rhotr

    def allocateQuantities(self, i, delta, solidity, chord, pitch, gamma, Cl, tbc):
        '''
        This function allocates the blade angles.
            inputs:
                i           -- incidence angle 
                delta       -- deviation angle 
                solidity    -- section solidity
                chord       -- section chord 
                pitch       -- section pitch  
                gamma       -- stagger angle
                Cl          -- lift coefficient
                tbc         -- thickness / chord
        '''

        # data allocation 
        self.i          = i 
        self.delta      = delta
        self.solidity   = solidity
        self.chord      = chord 
        self.pitch      = pitch
        self.gamma      = gamma 
        self.Cl         = Cl
        self.tbc        = tbc

    def mFlux(self):
        '''
        This function computes the mass flux in the stream tube.
            -- the axial fluid speed and the density are considered constant in the stream tube.
        '''

        # mass flux computation
        massFlux = np.pi * (self.tip**2 - self.bottom**2) * self.rho * self.Va

        return massFlux 