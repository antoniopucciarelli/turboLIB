# TURBOMACHINERY -- LIBRARY FOR THE INITIAL TURBOMACHINERY DESIGN
# AUTHOR: antonio pucciarelli 
#
# PROGRAM DESCRIPTION
#   TURBOMACHINERY DESIGN CLASS:
#       this script sets the stage class for the turbomachinery initial design 
#       

import numpy as np
import matplotlib.pyplot as plt 
from turboBlade import blade
import warnings

class stage:
    '''
    Stage object to be used in turbomachine object
    '''
    def __init__(self, ID, turboType = 'undefined', rotor = rotor(0), stator = stator(0), interstageGap = np.NaN, intrastageGap = np.NaN, velocity = np.array([[0,0],[0,0],[0,0]]), origin = np.array([np.NaN, np.NaN]), omega = np.NaN, Rm = np.NaN):
        '''
        Stage object declaration
            stage properties:
                ID              -- object identifier
                turboType       -- stage type   
                                -- avaliable: turbine/compressor
                rotor           -- rotor object
                                -- rotor object data stored in turboRotor.py
                stator          -- stator object 
                                -- stator object datat stored in turboStator.py 
                interstageGap   -- distance between this stage and its following stage
                intrastageGap   -- distance between rotor-stator in the stage
                velocity        -- 3x2 matrix that stores the absolute velocity triangles of the flow
                                -- stator(rotor) inlet -- stator(rotor) outlet -- rotor(stator) outlet
                origin0         -- absolute position of the 1st blade/vane in the turbomachine
                origin1         -- absolute position of the 2nd blade/vane in the turbomachine
                omega           -- revolution speed of the rotor
                Rm              -- average radius
        '''

        self.ID                     = ID
        self.turboType              = turboType
        self.rotor                  = rotor
        self.stator                 = stator
        self.interstageGap          = interstageGap
        self.intrastageGap          = intrastageGap
        self.origin0                = origin[0]
        self.origin1                = origin[1]
        self.omega                  = omega
        self.Rm                     = Rm
        self.setVelocity(velocity, rotor, stator)
        self.setRevSpeed()

    def setRevSpeed(self):
        '''
        Stage velocity rotor speed
        '''
        
        if self.Rm != np.NaN and self.omega != np.Nan:
            self.revSpeed = self.Rm * self.omega
        else:
            self.revSpeed = 0

    def setVelocity(self, velocity, rotor, stator):
        '''
        Stage velocity traingles definition
            -- stage object has self.velocity that stores the absolute velocity triangles in a 3x2 matrix 
            -- stage object has a rotor object and a stator object in it 
                -- this function stores the absolute velocity triangles also in these 2 objects 
            -- the inpute velocity has to be from the inlet to the outlet
                -- the ripartition of the velocity on the rotor and stator will depend on the stage type
                    -- if the stage type is compressor: 1st rotor 2nd stator
                    -- if the stage type is turbine: 1st statro 2nd rotor
            -- the allocation does not take into account the changes of velocity triangles in the intrastage part of the stage
        '''

        # velocity allocation 
        if velocity.all() == np.zeros([3,2]).all():
            if self.turboType == 'compressor':
                self.velocity = np.array([[rotor.velocity[0,:]],[rotor.velocity[1,:]],[stator.velocity[1,:]]])
            elif self.turboType == 'turbine':
                self.velocity = np.array([[stator.velocity[0,:]],[stator.velocity[1,:]],[rotor.velocity[0,:]]])
            else:
                warnings.warn('self.turboType = undefined ---> setting all velocities to 0.', DeprecationWarning, stacklevel=2)
        else:
            self.velocity = velocity
            # allocating stage velocities to stator and rotor
            if self.turboType == 'compressor':
                self.rotor.setVelocity(velocity[0:1,:])
                self.stator.setVelocity(velocity[1:2,:])
            elif self.turboType == 'turbine':
                self.stator.setVelocity(velocity[0:1,:])
                self.rotor.setVelocity(velocity[1:2,:])
            else:
                warnings.warn('self.rotor.velVec == self.stator.velVec == np.zeros([2,2]) because self.turboType = undefined.', DeprecationWarning, stacklevel=2)
    
    def print(self):
        '''
        Stage properties printout
        '''
    
        if self.turboType == 'turbine':
            self.stator.print()
            self.rotor.print()
        else:
            self.rotor.print()
            self.stator.print()

    def setInterstageGap(self, interstageGap):
        '''
        Stage interstage gap setting 
        ''' 

        self.interstageGap = interstageGap

    def setIntrastageGap(self, intrastageGap):
        '''
        Stage intrastage gap setting
        '''

        self.intrastageGap = intrastageGap

    def setRotor(self, rotor):
        '''
        Stage rotor allocation
        '''

        self.rotor = rotor

    def setStator(self, stator):
        ''' 
        Stage stator allocation 
        '''

        self.stator = stator
