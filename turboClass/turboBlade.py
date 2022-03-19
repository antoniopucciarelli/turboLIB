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
import warnings

class blade:
    '''
    Blade object, it is used in the stage object.
        AIM:
            --- blade geometry description 
            --- blade thermodynamics 
    '''

    def __init__(self, ID, turboType, turboAim, angle = np.NaN, pitch = np.NaN, chord = np.NaN, absVelocity = np.array([[0,0],[0,0]]), relVelocity = np.array([[0,0],[0,0]]), revVelocity = 0, airfoilName = 'NoName', position=np.NaN):
        '''
        Rotor object declaration: 
            variables:
                ID          --- blade identifier
                turboType   --- blade type: stator/rotor
                turboAim    --- blade aim: compressor/turbine
                angle       --- angle between the chord and the axial direction 
                pitch       --- distance between 2 vanes/blades
                chord       --- vane/blade aerodynamic chord
                airfoilName --- name of the airfoil used in the rotor/stator
                posiiton    --- blade leading edge position 
                absVelocity --- 2x2 vector that stores the absolute velocity at the leading edge and the trailing edge of the blade
                            --- e.g. velocity = | V1x, V1y |
                            ---                 | V2x, V2y |
                relVelocity --- 2x2 vector that store the relative velocity at the leading edge and the trailing edge of the blade 
                            --- e.g. velocity = | W1x, W1y |
                            ---                 | W2x, W2y |
                revVelocity --- revolution velocity only if the blade is a rotor
                            --- e.g. velocity = U
        '''
        
        if turboType != 'stator' or turboType != 'rotor':
            raise ValueError('Invalid blade.turboType -> it should be rotor or stator')
        if turboAim != 'compressor' or turboAim != 'turbine':
            raise ValueError('Invalid blade.turboAim -> it should be compressor or turbine')
        if turboType == 'stator' and revVelocity != 0:
            raise ValueError('Invalid blade.relVelocity -> for a stator the revVelocity == 0')

        self.ID             = ID
        self.turboType      = turboType
        self.angle          = angle
        self.pitch          = pitch
        self.chord          = chord
        self.absVelocity    = absVelocity
        self.relVelocity    = relVelocity
        self.revVelocity    = revVelocity
        self.airfoilName    = airfoilName
        self.position       = position
        self.setSolidity()
    
    def print(self):
        ''' 
        Blade properties printout
        '''
            
        print('\tID:           {0:>10d}'.format(self.ID))
        print('\ttype:         {0:>10s}'.format(self.turboType))
        print('\taim:          {0:>10s}'.format(self.turboAim))
        print('\tairfoil name: {0:>10s}'.format(self.airfoilName))
        print('\tangle:        {0:>10.2f}'.format(self.angle))
        print('\tchord:        {0:>10.2f}'.format(self.chord))
        print('\tpitch:        {0:>10.2f}'.format(self.pitch))
        print('\tLE position:  {0:>10.2f}'.format(self.position))
        print('\tvelocity:')
        print('\trotation')
        print('\t-- U =    {0:4.3f}'.format(self.revVelocity))
        print('\tabsolute')
        print('\t-- V1 = [ {0:4.3f} {1:4.3f}]'.format(self.absVelocity[0,0], self.absVelocity[0,1]))
        print('\t-- V2 = [ {0:4.3f} {1:4.3f}]'.format(self.absVelocity[1,0], self.absVelocity[1,1]))
        print('\trelative')
        print('\t-- W1 = [ {0:4.3f} {1:4.3f}]'.format(self.relVelocity[0,0], self.relVelocity[0,1]))
        print('\t-- W2 = [ {0:4.3f} {1:4.3f}]'.format(self.relVelocity[1,0], self.relVelocity[1,1]))
    
    def setAngle(self, angle):
        '''
        Rotor aerodynamic angle with respect to the axial flow direction
        '''
        self.angle = angle 

    def setSolidity(self):
        '''
        Rotor solidity computation
        '''

        if self.chord != np.NaN and self.pitch != np.NaN:
            self.solidity = self.chord * self.pitch
        else:
            self.solidity = np.NaN
                
    def setVelocity(self, velocity):
        ''' 
        Rotor velocity vector allocation 
        '''

        if velocity.shape == (2,2):
            self.velocity = velocity
        else:
            raise ValueError('Input error: velVec dimension is wrong.')

    def setFoil(self, foilName, plot=False):
        '''
        Rotor blade airfoil shape setting
            foil generation tip: 
                --- "foilName.txt" file should be made with xFoil
                --- each airfoil "foilName.txt" file should be stored in /data/airfoils/ 
        '''

        self.foil = np.loadtxt(foilName, skiprows=1)

        if plot:
            plt.figure()
            plt.plot(self.foil[:,0], self.foil[:,1], 'b')
            plt.title('{0:s}.ID = {1:d}\n{0:s}.chord = {2:2.2f}'.format(self.turboType, self.ID, self.chord))
            plt.axis('equal')
            plt.show()

    def foilRotation(self, unit='deg'):
        '''
        Airfoil rotation function

            function needs:
                self.foil       -- airfoil coordinates data
                                -- expressed in [x, y] numpy array format
                self.angle      -- airfoil angle with respect to the axial direction of the 
                                -- in radiant
                unit            -- defines angle units
                                -- if == deg the foil angle will be converted in radiants
        '''

        # profile rotation due to stagger angle
        gamma = self.angle
        if unit == 'deg':
            gamma = np.deg2rad(gamma)

        rotMatrix = np.matrix([[np.cos(gamma), -np.sin(gamma)],[np.sin(gamma), np.cos(gamma)]])

        # coordinate rotation
        for ii in range(self.foil.shape[0]):
            self.foil[ii,:] = np.matmul(rotMatrix, self.foil[ii,:])

    def plot(self):
        '''
        Rotor airfoil plot
        '''

        plt.figure()
        plt.plot(self.foil[:,0], self.foil[:,1], 'b')
        plt.title('{0:s}.ID = {1:d}\n{0:s}.chord = {2:2.2f}'.format(self.turboType, self.ID, self.chord))
        plt.axis('equal')
        plt.show()

