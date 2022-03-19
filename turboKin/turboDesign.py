# TURBOMACHINERY -- LIBRARY FOR THE INITIAL TURBOMACHINERY DESIGN
# AUTHOR: antonio pucciarelli 
#
# PROGRAM DESCRIPTION
#   TURBOMACHINERY DESIGN: 
#       -- TYPE:
#           -- axial compressor 
#           -- axial turbine 
#       
#       -- FUNCTIONS:
#           -- rotor design 
#           -- stator design 
#           -- blade computation and initial performance 
#       

# importing libraries
import matplotlib.pyplot as plt 
import numpy as np
import warnings

class stator:
    '''
    '''
    def __init__(self, ID, metalAngle = np.NaN, solidity = np.NaN, velVec = np.array([[0,0],[0,0]])):
        self.ID = ID
        self.metalAngle = metalAngle
        self.solidity = solidity
        self.velVec = velVec
        
    def setMetalAngle(self, metalAngle):
        self.metalAngle = metalAngle

    def setSolidity(self, solidity):
        # value check
        if solidity < 0.0:
            raise ValueError('Input error: solidity < 0.')
        else:
            self.solidity = solidity

    def setVelVec(self, velVec):
        if velVec.shape == (2,2):
            self.velVec = velVec
        else:
            raise ValueError('Input error: velVec dimension is wrong.')

    def setFoil(self, foilName, plot=False):
        # foilName.txt file should be made with xFoil
        self.foil = np.loadtxt(foilName,skiprows=1)

        if plot:
            plt.figure()
            plt.plot(self.foil[:,0], self.foil[:,1], 'b')
            plt.title(self.ID)
            plt.axis('equal')
            plt.show()

    def foilRotation(self, dim='deg'):
        '''
        Airfoil rotation function

            function input:
                self
                +.foil          -- airfoil coordinates data
                                -- expressed in [x, y] numpy array format
                +.metalAngle    -- airfoil metal angle
                dim             -- sets the angle dimensions 
                                -- radiants: rad
                                -- degrees: deg
        '''

        # profile rotation due to stagger angle
        gamma = self.metalAngle
        if dim == 'deg':
            # angle conversion
            gamma = np.deg2rad(gamma)
        # matrix computation
        rotMatrix = np.matrix([[np.cos(gamma), -np.sin(gamma)],[np.sin(gamma), np.cos(gamma)]])

        # coordinate rotation
        for ii in range(self.foil.shape[0]):
            self.foil[ii,:] = np.matmul(rotMatrix, self.foil[ii,:])

class rotor:
    '''
    '''
    def __init__(self, ID, metalAngle = np.NaN, solidity = np.NaN, velVec = np.array([[0,0],[0,0]])):
        self.ID = ID
        self.metalAngle = metalAngle 
        self.solidity = solidity 
        self.velVec = velVec

    def setMetalAngle(self, metalAngle):
        self.metalAngle = metalAngle

    def setSolidity(self, solidity):
        # value check
        if solidity < 0.0:
            raise ValueError('Input error: solidity < 0.')
        else:
            self.solidity = solidity

    def setVelVec(self, velVec):
        # velocity array dimension check
        if velVec.shape == (2,2):
            self.velVec = velVec
        else:
            raise ValueError('Input error: velVec dimension is wrong.')
    
    def setFoil(self, foilName, plot=False):
        # foilName.txt file should be made with xFoil
        self.foil = np.loadtxt(foilName,skiprows=1)
        
        if plot:
            plt.figure()
            plt.plot(self.foil[:,0], self.foil[:,1], 'b')
            plt.title(self.ID)
            plt.axis('equal')
            plt.show()

    def foilRotation(self, dim='deg'):
        '''
        Airfoil rotation function

            function input:
                self
                +.foil          -- airfoil coordinates data
                                -- expressed in [x, y] numpy array format
                +.metalAngle    -- airfoil metal angle
                dim             -- sets the angle dimensions 
                                -- radiants: rad
                                -- degrees: deg
        '''

        # profile rotation due to stagger angle
        gamma = self.metalAngle
        if dim == 'deg':
            # angle conversion
            gamma = np.deg2rad(gamma)
        # matrix computation
        rotMatrix = np.matrix([[np.cos(gamma), -np.sin(gamma)],[np.sin(gamma), np.cos(gamma)]])

        # coordinate rotation
        for ii in range(self.foil.shape[0]):
            self.foil[ii,:] = np.matmul(rotMatrix, self.foil[ii,:])

    def plot(self):
        plt.figure()
        plt.plot(self.foil[:,0], self.foil[:,1], 'b')
        plt.title(self.ID)
        plt.axis('equal')
        plt.show()

class stage:
    '''
    '''
    def __init__(self, ID, turboType = 'undefined', rotor = rotor(0), stator = stator(0), interstageGap = np.NaN, intrastageGap = np.NaN, velVec = np.array([[0,0],[0,0],[0,0]])):
        self.ID = ID
        self.turboType = turboType
        self.rotor = rotor
        self.stator = stator
        self.interstageGap = interstageGap
        self.intrastageGap = intrastageGap

        # velocity allocation 
        if velVec.all() == np.zeros([3,2]).all():
            if self.turboType == 'compressor':
                self.velVec = np.array([[rotor.velVec[0,:]],[rotor.velVec[1,:]],[stator.velVec[1,:]]])
            elif self.turboType == 'turbine':
                self.velVec = np.array([[stator.velVec[0,:]],[stator.velVec[1,:]],[rotor.velVec[0,:]]])
            else:
                warnings.warn('self.turboType = undefined ---> setting all velocities to 0.', DeprecationWarning, stacklevel=2)
        else:
            self.velVec = velVec
            # allocating stage velocities to stator and rotor
            if self.turboType == 'compressor':
                self.rotor.setVelVec(velVec[0:1,:])
                self.stator.setVelVec(velVec[1:2,:])
            elif self.turboType == 'turbine':
                self.stator.setVelVec(velVec[0:1,:])
                self.rotor.setVelVec(velVec[1:2,:])
            else:
                warnings.warn('self.rotor.velVec == self.stator.velVec == np.zeros([2,2]) because self.turboType = undefined.', DeprecationWarning, stacklevel=2)
                        
    def setInterstageGap(self, interstageGap):
        self.interstageGap = interstageGap

    def setIntrastageGap(self, intrastageGap):
        self.intrastageGap = intrastageGap

    def setRotor(self, rotor):
        self.rotor = rotor

    def setStator(self, stator):
        self.stator = stator

stage1 = stage(1, 'turbine', rotor = rotor(1), stator = stator(1))
stage1.stator.setFoil('naca1212')
stage1.rotor.setFoil('naca5016', True)
stage1.rotor.setMetalAngle(10)
stage1.rotor.foilRotation()
stage1.rotor.plot()
