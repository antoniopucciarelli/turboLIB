# TURBOMACHINERY -- LIBRARY FOR THE INITIAL TURBOMACHINERY DESIGN
# AUTHOR: antonio pucciarelli 
#
# PROGRAM DESCRIPTION
#   TURBOMACHINERY DESIGN CLASS:
#       this script sets the rotor/stator blade geometry 
#           -- for now the blade geometry is set as a NACA-65 class airfoil modified with respect to the target Cl of the section       
#

# importing libraries
import numpy as np 
import matplotlib.pyplot as plt 

# importing coords
class geometryData:
    def __init__(self, file):
        '''
        Airfoil data object, it is used for blade geometry generation.
        '''
        # data extraction
        data = np.loadtxt(file, comments='#', usecols=(0,1,2))
        # data allocation
        self.x = np.array(data[:,0] / 100)
        self.y = np.array(data[:,1] / 100)
        self.t = np.array(data[:,2] / 100)

    def geometryFitting(self, Cl, chord=1, plot=False):
        '''
        Airfoil geometry fitting, it is used for geometry modification with respect to the given Cl.
            inputs: 
                Cl      -- blade lift coefficient
                chord   -- airfoil chord
                plot    -- showing airfoil coordinates
        '''

        # camber line computation
        yCamber = self.y * Cl * chord
        
        # thickness computation  
        t = self.t * Cl * chord 
        
        # upper line computation 
        yUpper = yCamber + t/2

        # lower line computation 
        yLower = yCamber - t/2  

        # extending chord 
        self.x = self.x * chord 

        if plot:
            plt.figure()
            plt.plot(self.x, yCamber, 'r')
            plt.plot(self.x, yUpper, 'b', linewidth=2)
            plt.plot(self.x, yLower, 'k', linewidth=2)
            plt.axis('equal')
            plt.show()

        # data allocation -- points 
        self.X       = self.x 
        self.Yupper  = yUpper
        self.Ylower  = yLower
        self.T       = t 
        self.Ycamber = yCamber

        # data allocation -- vectors 
        self.upper  = np.array([self.x, yUpper, np.zeros(self.x.shape[0])])
        self.lower  = np.array([self.x, yLower, np.zeros(self.x.shape[0])])
        self.camber = np.array([self.x, yCamber, np.zeros(self.x.shape[0])]) 

    def geometryRotation(self, yaw, pitch=0, plot=False):
        '''
        Airfoil geometry rotation:
            -- yaw      -> [deg] blade to blade plane -- metal angle  
            -- pitch    -> [deg] meridional plane     -- due to annulus diameter variation
            -- plot     -- boolean value for the blade plotting
        '''

        # alpha == roll == 0
        roll = 0
        cA = np.cos(np.deg2rad(roll))
        sA = np.sin(np.deg2rad(roll))
        # beta == pitch --> blade angle in meridional plane
        cB = np.cos(np.deg2rad(pitch))
        sB = np.sin(np.deg2rad(pitch))
        # gamma == yaw --> blade to blade AOA 
        cG = np.cos(np.deg2rad(yaw))
        sG = np.sin(np.deg2rad(yaw))

        # 3D rotation matrix 
        # | cosB * cosG, sinA * sinB * cosG - cosA * sinG, cosA * sinB * cosG + sinA * sinG  |
        # | cosB * sinG, sinA * sinB * sinG + cosA * cosG, cosA * sinB * sinG - sinA * cosG  |
        # |      - sinB,                      sinA * cosB,                      cosA * cosB  |
        rotMatrix = np.matrix([[cB * cG, sA * sB * cG - cA * sG, cA * sB * cG + sA * sG], 
                               [cB * sG, sA * sB * sG + cA * cG, cA * sB * sG - sA * cG],
                               [   - sB,                sA * cB,                cA * cB]])

        # rotating upper surface 
        self.upper = np.matmul(rotMatrix, self.upper)

        # rotation lower surface 
        self.lower = np.matmul(rotMatrix, self.lower)

        # rotation camber line 
        self.camber = np.matmul(rotMatrix, self.camber)

        # variables allocation 
        # -- there was a problem with the 3D plotting (x,y,z) values 
        # -- in order to correct this problem -> raw allocation of data
        x       = np.zeros(self.upper.shape[1],)
        yUpper  = np.zeros(self.upper.shape[1],)
        zUpper  = np.zeros(self.upper.shape[1],)
        yLower  = np.zeros(self.upper.shape[1],)
        zLower  = np.zeros(self.upper.shape[1],)
        yCamber = np.zeros(self.upper.shape[1],)
        zCamber = np.zeros(self.upper.shape[1],)

        for ii in range(self.upper.shape[1]):
            x[ii]       = self.upper[0,ii]
            yUpper[ii]  = self.upper[1,ii]
            zUpper[ii]  = self.upper[2,ii]
            yLower[ii]  = self.lower[1,ii]
            zLower[ii]  = self.lower[2,ii]
            yCamber[ii] = self.camber[1,ii]
            zCamber[ii] = self.camber[2,ii]

        # data reallocation 
        self.upper  = np.array([x, yUpper, zUpper])
        self.lower  = np.array([x, yLower, zLower])
        self.camber = np.array([x, yCamber, zCamber])        

        if plot:
            # plotting blade section 
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.set_box_aspect((np.ptp(x), np.ptp(yUpper), np.ptp(zUpper)))
            plt.plot(self.X, self.Yupper, 0, '--r')
            plt.plot(self.x, self.Ylower, 0, '--r')
            plt.plot(self.x, self.Ycamber, 0, '-*r')
            plt.plot(x, yUpper, zUpper, 'k')
            plt.plot(x, yLower, zLower, 'k')
            plt.plot(x, yCamber, zCamber, '*k')
            plt.title(r'$\alpha_g = $ {0:.2f}'.format(yaw) + '\n' + r'$\zeta_g = $ {0:.2f}'.format(pitch))
            plt.show()

    def middleChord(self):
        '''
        This function computes the middle chord point of the airfoil.
            output:
                middlePoint -- airfoil chord middle point 
        '''

        self.middlePoint = (self.camber[:,0] + self.camber[:,-1]) / 2 

        return self.middlePoint

    def geometryTranslation(self, translationVec, height, plot=False):
        '''
        Airfoil geometry translation:
            inputs:
                translationVec  -- translation vector [x, y, z]
                height          -- airfoil middle chord point height
                plot            -- boolean value for blade plotting 
            !!! the blade should first rotated and then translated !!!
        '''

        # setting up translationVec 
        # the airfoil translation is relative the hub chord centre with the addition of height translation 
        translationVec = translationVec + np.array([0.0, 0.0, height]) - self.middleChord()

        # variables allocation 
        # -- there was a problem with the 3D plotting (x,y,z) values 
        # -- in order to correct this problem -> raw allocation of data
        x       = np.zeros(self.upper.shape[1],)
        yUpper  = np.zeros(self.upper.shape[1],)
        zUpper  = np.zeros(self.upper.shape[1],)
        yLower  = np.zeros(self.upper.shape[1],)
        zLower  = np.zeros(self.upper.shape[1],)
        yCamber = np.zeros(self.upper.shape[1],)
        zCamber = np.zeros(self.upper.shape[1],)

        for ii in range(self.upper.shape[1]):
            x[ii]       = self.upper[0,ii]
            yUpper[ii]  = self.upper[1,ii]
            zUpper[ii]  = self.upper[2,ii]
            yLower[ii]  = self.lower[1,ii]
            zLower[ii]  = self.lower[2,ii]
            yCamber[ii] = self.camber[1,ii]
            zCamber[ii] = self.camber[2,ii]

        if plot:
            # plotting blade section 
            # initial blade + rotated blade
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.set_box_aspect((np.ptp(x), np.ptp(yUpper), np.ptp(zUpper)))
            plt.plot(self.X, self.Yupper, 0, '--r')
            plt.plot(self.x, self.Ylower, 0, '--r')
            plt.plot(self.x, self.Ycamber, 0, '-*r')
            plt.plot(x, yUpper, zUpper, '--b')
            plt.plot(x, yLower, zLower, '--b')
            plt.plot(x, yCamber, zCamber, '-*b')

        # upper surface translation 
        x      = x      + translationVec[0]
        yUpper = yUpper + translationVec[1]
        zUpper = zUpper + translationVec[2]

        # lower surface translation 
        yLower = yLower + translationVec[1]
        zLower = zLower + translationVec[2] 

        # camber line translation
        yCamber = yCamber + translationVec[1]
        zCamber = zCamber + translationVec[2]

        # data reallocation 
        self.upper  = np.array([x, yUpper, zUpper])
        self.lower  = np.array([x, yLower, zLower]) 
        self.camber = np.array([x, yCamber, zCamber])

        if plot:
            # plot the translated blade
            plt.plot(x, yUpper, zUpper, 'k')
            plt.plot(x, yLower, zLower, 'k')
            plt.plot(x, yCamber, zCamber, '--k')
            plt.show()

def writeFacet(file, versor, vec1, vec2, vec3):
    '''
    This function allows writing facet in stl format.
        inputs:
            file    -- file object 
            versor  -- normal versor vector 
            vec1    -- facet triangle (1) point
            vec2    -- facet triangle (2) point
            vec3    -- facet triangle (3) point
    '''
    # writing facet in stl format into file 
    file.write('\tfacet normal {0:f} {1:f} {2:f}\n'.format(versor[0], versor[1], versor[2]))
    file.write('\t\touter loop\n')
    file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec1[0], vec1[1], vec1[2]))
    file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec2[0], vec2[1], vec2[2]))
    file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec3[0], vec3[1], vec3[2]))
    file.write('\t\tendloop\n')
    file.write('\tendfacet\n')

def STLsaving(airfoils, STLname='cad', containerPath='container/'):
    '''
    This function saves the blade in .stl format.
        inputs: 
            airfoils    -- tuple of airfoils objects
            STLname     -- .stl file name
        !!! it is assumed that each airfoil has the same number of description points !!!
        !!! it is assumed that the each airfoil section element is in sequence with respect the hub !!!
    '''
    # stl generation 
    file = open(containerPath + STLname + '.stl', 'w')

    # begin upper + lower surface stl
    file.write('solid bladeSurface\n')

    # loop over all the span positions
    for jj in range(len(airfoils)-1):

        # upper surface stl face generation 
        for ii in range(airfoils[jj].upper.shape[1]-1):
            # versor computation
            vec1 = np.array(airfoils[jj].upper[:,ii])
            vec2 = np.array(airfoils[jj].upper[:,ii+1])
            vec3 = np.array(airfoils[jj+1].upper[:,ii])
            versor = np.cross(vec3 - vec1, vec2 - vec1)
            # direction check -- for the upper surface the versor direction is towards y > 0 
            if versor[1] < 0: 
                versor = - versor 
            # vector normalization
            versor = versor / np.linalg.norm(versor)
            # writing data
            writeFacet(file, versor, vec1, vec2, vec3)

            # versor computation
            vec1 = np.array(airfoils[jj+1].upper[:,ii])
            vec2 = np.array(airfoils[jj+1].upper[:,ii+1])
            vec3 = np.array(airfoils[jj].upper[:,ii+1])
            versor = - np.cross(vec3 - vec1, vec2 - vec1)
            # direction check -- for the upper surface the versor direction is towards y > 0 
            if versor[1] < 0: 
                versor = - versor 
            # vector normalization 
            versor = versor / np.linalg.norm(versor)
            # writing data
            writeFacet(file, versor, vec1, vec2, vec3)

        # lower surface stl face generation 
        for ii in range(airfoils[jj].lower.shape[1]-1):
            # versor computation
            vec1 = np.array(airfoils[jj].lower[:,ii])
            vec2 = np.array(airfoils[jj].lower[:,ii+1])
            vec3 = np.array(airfoils[jj+1].lower[:,ii])
            versor = np.cross(vec3 - vec1, vec2 - vec1)
            # direction check -- for the upper surface the versor direction is towards y < 0 
            if versor[1] > 0: 
                versor = - versor 
            # vectorn normalization
            versor = versor / np.linalg.norm(versor)
            # writing data
            writeFacet(file, versor, vec1, vec2, vec3)

            # versor computation
            vec1 = np.array(airfoils[jj+1].lower[:,ii])
            vec2 = np.array(airfoils[jj+1].lower[:,ii+1])
            vec3 = np.array(airfoils[jj].lower[:,ii+1])
            versor = - np.cross(vec3 - vec1, vec2 - vec1)
            # direction check -- for the upper surface the versor direction is towards y < 0 
            if versor[1] > 0: 
                versor = - versor 
            # vector normalization
            versor = versor / np.linalg.norm(versor)
            # writing data
            writeFacet(file, versor, vec1, vec2, vec3)
    
    # end blade 
    file.write('endsolid bladeSurface\n')

    # begin bottom blade surface -- hub
    file.write('solid bladeBottom\n')

    # closing airfoil [hub]
    for ii in range(airfoils[0].upper.shape[1]-2):
        # versor computation 
        vec1 = np.array(airfoils[0].camber[:,ii+1])
        vec2 = np.array(airfoils[0].upper[:,ii])
        vec3 = np.array(airfoils[0].upper[:,ii+1])
        # versor computation 
        versor = [0, 0, -1]
        # writing data
        writeFacet(file, versor, vec1, vec2, vec3)

    for ii in range(airfoils[0].lower.shape[1]-2):
        # versor computation 
        vec1 = np.array(airfoils[0].camber[:,ii+1])
        vec2 = np.array(airfoils[0].lower[:,ii])
        vec3 = np.array(airfoils[0].lower[:,ii+1])
        # versor computation 
        versor = [0, 0, -1]
        # writing data
        writeFacet(file, versor, vec1, vec2, vec3)

    for ii in range(airfoils[0].upper.shape[1]-2):
        # versor computation 
        vec1 = np.array(airfoils[0].camber[:,ii+2])
        vec2 = np.array(airfoils[0].upper[:,ii+1])
        vec3 = np.array(airfoils[0].camber[:,ii+1])
        # versor computation 
        versor = [0, 0, -1]
        # writing data
        writeFacet(file, versor, vec1, vec2, vec3)
    
    for ii in range(airfoils[0].upper.shape[1]-2):
        # versor computation 
        vec1 = np.array(airfoils[0].camber[:,ii+2])
        vec2 = np.array(airfoils[0].lower[:,ii+1])
        vec3 = np.array(airfoils[0].camber[:,ii+1])
        # versor computation 
        versor = [0, 0, -1]
        # writing data
        writeFacet(file, versor, vec1, vec2, vec3)

    # end blade bottom surface -- hub 
    file.write('endsolid bladeBottom\n')

    # begin top blade surface -- tip
    file.write('solid bladeTop\n')

    # closing airfoil [tip]
    for ii in range(airfoils[-1].upper.shape[1]-2):
        # versor computation 
        vec1 = np.array(airfoils[-1].camber[:,ii+1])
        vec2 = np.array(airfoils[-1].upper[:,ii])
        vec3 = np.array(airfoils[-1].upper[:,ii+1])
        # versor computation 
        versor = [0, 0, 1]
        # writing data
        file.write('\tfacet normal {0:f} {1:f} {2:f}\n'.format(versor[0], versor[1], versor[2]))
        file.write('\t\touter loop\n')
        file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec1[0], vec1[1], vec1[2]))
        file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec2[0], vec2[1], vec2[2]))
        file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec3[0], vec3[1], vec3[2]))
        file.write('\t\tendloop\n')
        file.write('\tendfacet\n')

    for ii in range(airfoils[-1].lower.shape[1]-2):
        # versor computation 
        vec1 = np.array(airfoils[-1].camber[:,ii+1])
        vec2 = np.array(airfoils[-1].lower[:,ii])
        vec3 = np.array(airfoils[-1].lower[:,ii+1])
        # versor computation 
        versor = [0, 0, 1]
        # writing data
        writeFacet(file, versor, vec1, vec2, vec3)

    for ii in range(airfoils[-1].upper.shape[1]-2):
        # versor computation 
        vec1 = np.array(airfoils[-1].camber[:,ii+2])
        vec2 = np.array(airfoils[-1].upper[:,ii+1])
        vec3 = np.array(airfoils[-1].camber[:,ii+1])
        # versor computation 
        versor = [0, 0, 1]
        # writing data
        writeFacet(file, versor, vec1, vec2, vec3)
    
    for ii in range(airfoils[-1].upper.shape[1]-2):
        # versor computation 
        vec1 = np.array(airfoils[-1].camber[:,ii+2])
        vec2 = np.array(airfoils[-1].lower[:,ii+1])
        vec3 = np.array(airfoils[-1].camber[:,ii+1])
        # versor computation 
        versor = [0, 0, 1]
        # writing data
        writeFacet(file, versor, vec1, vec2, vec3)
    
    # end blade top surface -- tip 
    file.write('endsolid bladeTop\n')

    # file closure 
    file.close()




