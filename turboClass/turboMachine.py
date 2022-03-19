# TURBOMACHINERY -- LIBRARY FOR THE INITIAL TURBOMACHINERY DESIGN
# AUTHOR: antonio pucciarelli 
#
# PROGRAM DESCRIPTION
#   TURBOMACHINERY DESIGN CLASS:
#       this script sets the turbomachine class for the turbomachinery initial design 
#       

import numpy as np
import matplotlib.pyplot as plt 
import warnings
from turboStage import *

class turbomachine:
    def __init__(self, ID, nStage, turboType):
        '''
        Turbomachine object declaration:
            variables:
        '''

        self.ID         = ID
        self.nStage     = nStage
        self.turboType  = turboType
        self.setStage()

    def setStage(self):
        '''
        Turbomachine object stages declaration 
        '''
        
        self.stage = [stage(ii+1, self.turboType, rotor(ii+1), stator(ii+1)) for ii in range(self.nStage)]
    
    def print(self):
        '''
        Turbomachine properties printout
        '''

        print('+++++++++++++++++++')
        print('TURBOMACHINE')
        print('-- ID      {0:d}'.format(self.ID))
        print('-- #stages {0:d}'.format(self.nStage))
        print('-- type    {0:s}'.format(self.turboType))
        print('*******************')
        for ii in range(self.nStage):
            print('---------------')
            print('stage  #{0:1}'.format(self.stage[ii].ID))
            print('---------------')
            self.stage[ii].print()
        print('+++++++++++++++++++')


    def turboPlot(self):
        ''' 
        Turbomacinery stages plot 
        ''' 

        plt.figure()
        
        self.computeGaps()

        for ii in range(self.nStage):
            if self.turboType == 'compressor':
                plt.plot(self.stage[ii].rotor.foil[:,0] + self.stage[ii].origin0, self.stage[ii].rotor.foil[:,1], 'b')
                plt.plot(self.stage[ii].stator.foil[:,0] + self.stage[ii].origin1, self.stage[ii].stator.foil[:,1], 'r')
            else:
                plt.plot(self.stage[ii].stator.foil[:,0] + self.stage[ii].origin0, self.stage[ii].stator.foil[:,1], 'r')
                plt.plot(self.stage[ii].rotor.foil[:,0] + self.stage[ii].origin1, self.stage[ii].rotor.foil[:,1], 'b')

        # important quantities printout
        print('+++++++++++++++++++++++++++++++++++++++')
        print('Stages analysis:')
        for ii in range(self.nStage):
            print('stage #{0}'.format(self.stage[ii].ID))
            print('-- rotor')
            print('\tairfoil name: {0:>10s}'.format(self.stage[ii].rotor.airfoilName))
            print('\tangle:        {0:>10.2f}'.format(self.stage[ii].rotor.angle))
            print('\tchord:        {0:>10.2f}'.format(self.stage[ii].rotor.chord))
            print('\tpitch:        {0:>10.2f}'.format(self.stage[ii].rotor.pitch))
            print('\tLE position:  {0:>10.2f}'.format(self.stage[ii].origin0))
            print('-- stator')
            print('\tairfoil name: {0:>10s}'.format(self.stage[ii].stator.airfoilName))
            print('\tangle:        {0:>10.2f}'.format(self.stage[ii].stator.angle))
            print('\tchord:        {0:>10.2f}'.format(self.stage[ii].stator.chord))
            print('\tpitch:        {0:>10.2f}'.format(self.stage[ii].stator.pitch))
            print('\tLE position:  {0:>10.2f}'.format(self.stage[ii].origin1))
        print('+++++++++++++++++++++++++++++++++++++++')
    
        plt.axis('equal')
        plt.show()
            
    def computeGaps(self):
        '''
        Turbomachinery function that computes all the gaps among stages and rotors-stators
        '''

        # setting up initial conditions
        origin = 0          # turbomachinery axial position origin set to 0

        for ii in range(self.nStage):
            if self.turboType == 'compressor':
                # allocation of the first blade absolute position
                self.stage[ii].origin0 = origin
                origin = origin \
                         + self.stage[ii].rotor.chord * np.cos(np.deg2rad(self.stage[ii].rotor.angle)) \
                         + self.stage[ii].intrastageGap
                # allocation of the second blade absolute position 
                self.stage[ii].origin1 = origin 
                origin = origin \
                         + self.stage[ii].stator.chord * np.cos(np.deg2rad(self.stage[ii].stator.angle)) 
                
                # last stage check 
                if ii != self.nStage-1: 
                   origin = origin + self.stage[ii].interstageGap
                
            else:
                # allocation of the first blade absolute position 
                self.stage[ii].origin0 = origin
                origin = origin \
                         + self.stage[ii].stator.chord * np.cos(np.deg2rad(self.stage[ii].stator.angle)) \
                         + self.stage[ii].intrastageGap
                # allocation of the second blade absolute position 
                self.stage[ii].origin1 = origin
                origin = origin \
                         + self.stage[ii].rotor.chord * np.cos(np.deg2rad(self.stage[ii].rotor.angle)) 
                
                # last stage check 
                if ii != self.nStage-1:
                    origin = origin + self.stage[ii].interstageGap
    
    def velocityVec(self):
        '''
        Velocity triangles plot over a turbomachine blade 
        '''
        
        for ii in range(self.nStage):
            self.stage[ii].setRevSpeed() 
            self.stage[ii].rotor.computeW(self.stage[ii].revSpeed)
                    
        fig, axis = plt.subplots(2)

        for ii in range(self.nStage):
            if self.turboType == 'compressor':
                axis[0].plot(self.stage[ii].rotor.foil[:,0] + self.stage[ii].origin0, self.stage[ii].rotor.foil[:,1], 'b')
                axis[0].plot(self.stage[ii].stator.foil[:,0] + self.stage[ii].origin1, self.stage[ii].stator.foil[:,1], 'r')






