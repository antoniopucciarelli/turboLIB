# TURBOMACHINERY -- LIBRARY FOR THE INITIAL TURBOMACHINERY DESIGN
# AUTHOR: antonio pucciarelli 
#
# PROGRAM DESCRIPTION
#   TURBOMACHINERY DESIGN CLASS:
#       this script sets the stator class for the turbomachinery initial design 
#       

import numpy as np 
import matplotlib.pyplot as plt 
import warnings

class stator:
    '''
    Stator object, it is used in the stage object.
        AIM:
            --- stator geoemtry description 
            --- stator thermodynamics 
    '''
    def __init__(self, ID, angle = np.NaN, pitch = np.NaN, chord = np.NaN, velocity = np.array([[0,0],[0,0]]), airfoilName = 'NoName'):
        '''
        stator object declaration: 
            variables:
                ID          --- stator name/identifier
                angle       --- angle between the chord and the axial direction 
                pitch       --- distance between 2 vanes/blades in the stator
                chord       --- vane/blade aerodynamic chord
                velocity    --- 2x2 vector that stores the absolute velocity ad the leading edge and the trainling edge of the stator vane
                            --- e.g. velocity = | V1x, V1y |
                            ---                 | V2x, V2y |
                airfoilName --- name of the airfoil used in the stator
        '''
        self.ID             = ID
        self.angle          = angle
        self.pitch          = pitch
        self.chord          = chord
        self.velocity       = velocity
        self.airfoilName    = airfoilName
        self.setSolidity()

    def print(self):
        ''' 
        Stator properties printout
        '''
        
        print('* STATOR *')
        print('-- ID    {0:d}'.format(self.ID))
        print('-- angle {0:2.2f}'.format(self.angle))
        print('-- chord {0:2.2f}'.format(self.chord))
        print('-- pitch {0:2.2f}\n'.format(self.pitch))

    def setAngle(self, angle):
        '''
        Stator metal angle computation
        '''
        self.angle = angle

    def setSolidity(self):
        '''
        Stator solidity computation
        '''

        if self.chord != np.NaN and self.pitch != np.NaN:
            self.solidity = self.chord * self.pitch
        else:
            self.solidity = np.NaN
                
    def setVelocity(self, velocity):
        ''' 
        Stator velocity vector allocation 
        '''

        if velocity.shape == (2,2):
            self.velocity = velocity
        else:
            raise ValueError('Input error: velVec dimension is wrong.')

    def setFoil(self, foilName, plot=False):
        '''
        Stator vane/blade airfoil shape setting
            foil generation tip: 
                --- "foilName.txt" file should be made with xFoil
                --- each airfoil "foilName.txt" file should be stored in /data/airfoils/ 
        '''

        self.foil = np.loadtxt(foilName,skiprows=1)
        self.chord = np.max(self.foil[:,0]) - np.min(self.foil[:,0])

        if plot:
            plt.figure()
            plt.plot(self.foil[:,0], self.foil[:,1], 'b')
            plt.title('stator.ID = %s\nstator.chord = %f'.format(self.ID, self.chord))
            plt.axis('equal')
            plt.show()

    def setScale(self, scaling):
        '''
        Stator vane/blade scaling
        '''

        self.foil = self.foil * scaling

    def foilRotation(self, unit='deg'):
        '''
        Airfoil rotation function

            function needs:
                self.foil       -- airfoil coordinates data
                                -- expressed in [x, y] numpy array format
                self.angle      -- airfoil angle with respect to the axial direction of the 
                                -- in radiant
                unit            -- angle dimension unit descriptor 
                                -- if != deg self.angle is not converted in radiants
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
        Stator airfoil plot
        '''

        plt.figure()
        plt.plot(self.foil[:,0], self.foil[:,1], 'b')
        plt.title('stator.ID = {0:d}\nstator.chord = {1:2.2f}'.format(self.ID, self.chord))
        plt.axis('equal')
        plt.show()
