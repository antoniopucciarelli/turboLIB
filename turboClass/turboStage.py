# TURBOMACHINERY -- LIBRARY FOR THE INITIAL TURBOMACHINERY DESIGN
# AUTHOR: antonio pucciarelli 
#
# PROGRAM DESCRIPTION
#   TURBOMACHINERY DESIGN CLASS:
#       this script sets the stage class for the turbomachinery initial design 
#       

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
