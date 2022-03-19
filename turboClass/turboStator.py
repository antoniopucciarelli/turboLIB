# TURBOMACHINERY -- LIBRARY FOR THE INITIAL TURBOMACHINERY DESIGN
# AUTHOR: antonio pucciarelli 
#
# PROGRAM DESCRIPTION
#   TURBOMACHINERY DESIGN CLASS:
#       this script sets the stator class for the turbomachinery initial design 
#       

class stator:
    '''
    Stator object, it is used in the stage object.
        AIM:
            --- stator geoemtry description 
            --- stator thermodynamics 
    '''
    def __init__(self, ID, metalAngle = np.NaN, pitch = np.NaN, chord = np.NaN, velocity = np.array([[0,0],[0,0]])):
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
        '''
        self.ID         = ID
        self.angle      = angle
        self.pitch      = pitch
        self.chord      = chord
        self.velocity   = velocity

        self.setSolidity()
        
    def setMetalAngle(self, angle):
        '''
        Stator metal angle computation
        '''
        pass 

    def setSolidity(self, solidity):
        '''
        Stator solidity computation
        '''

        if self.chord != np.NaN and self.pitch != np.Nan
            self.solidity = self.chord * self.pitch
        else:
            self.solidity = np.NaN
                
    def setVelVec(self, velocity):
        ''' 
        Stator velocity vector allocation 
        '''

        if velVec.shape == (2,2):
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

        if plot:
            plt.figure()
            plt.plot(self.foil[:,0], self.foil[:,1], 'b')
            plt.title('stator.ID = %s\nstator.chord = %f'.format(self.ID, stator.chord))
            plt.axis('equal')
            plt.show()

    def foilRotation(self):
        '''
        Airfoil rotation function

            function needs:
                self.foil       -- airfoil coordinates data
                                -- expressed in [x, y] numpy array format
                self.angle      -- airfoil angle with respect to the axial direction of the 
                                -- in radiant
        '''

        # profile rotation due to stagger angle
        gamma = self.angle
        rotMatrix = np.matrix([[np.cos(gamma), -np.sin(gamma)],[np.sin(gamma), np.cos(gamma)]])

        # coordinate rotation
        for ii in range(self.foil.shape[0]):
            self.foil[ii,:] = np.matmul(rotMatrix, self.foil[ii,:])

stage1 = stage(1, 'turbine', rotor = rotor(1), stator = stator(1))
stage1.stator.setFoil('naca1212')
stage1.rotor.setFoil('naca5016', True)
stage1.rotor.setMetalAngle(10)
stage1.rotor.foilRotation()
stage1.rotor.plot()
