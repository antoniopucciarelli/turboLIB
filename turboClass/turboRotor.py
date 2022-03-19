# TURBOMACHINERY -- LIBRARY FOR THE INITIAL TURBOMACHINERY DESIGN
# AUTHOR: antonio pucciarelli 
#
# PROGRAM DESCRIPTION
#   TURBOMACHINERY DESIGN CLASS:
#       this script sets the rotor class for the turbomachinery initial design 
#       

class rotor:
    '''
    Rotor object, it is used in the stage object.
        AIM:
            --- rotor geoemtry description 
            --- rotor thermodynamics 
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
