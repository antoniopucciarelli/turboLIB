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

# importing NACA-65 coords
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

        if plot:
            plt.figure()
            plt.plot(self.x, yCamber, 'r')
            plt.plot(self.x, yUpper, 'b', linewidth=2)
            plt.plot(self.x, yLower, 'k', linewidth=2)
            plt.axis('equal')
            plt.show()

        # data allocation -- points 
        self.X      = self.x 
        self.Yupper = yUpper
        self.Ylower = yLower
        self.T      = t 

        # data allocation -- vectors 
        self.upper = np.array([self.x, yUpper, np.zeros(self.x.shape[0])])
        self.lower = np.array([self.x, yLower, np.zeros(self.x.shape[0])])

        return yCamber, yUpper, yLower, t 

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

        # variables allocation 
        # -- there was a problem with the 3D plotting (x,y,z) values 
        # -- in order to correct this problem -> raw allocation of data
        x  = np.zeros(self.upper.shape[1],)
        y  = np.zeros(self.upper.shape[1],)
        z  = np.zeros(self.upper.shape[1],)
        y1 = np.zeros(self.upper.shape[1],)
        z1 = np.zeros(self.upper.shape[1],)
        
        for ii in range(self.upper.shape[1]):
            x[ii]  = self.upper[0,ii]
            y[ii]  = self.upper[1,ii]
            z[ii]  = self.upper[2,ii]
            y1[ii] = self.lower[1,ii]
            z1[ii] = self.lower[2,ii]

        # data reallocation 
        self.upper = np.array([x, y, z])
        self.lower = np.array([x, y1, z1])        

        if plot:
            # plotting blade section 
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.set_box_aspect((np.ptp(x), np.ptp(y), np.ptp(z)))
            plt.plot(self.X, self.Yupper, 0, '--r')
            plt.plot(self.x, self.Ylower, 0, '--r')
            plt.plot(x, y , z,  'k')
            plt.plot(x, y1, z1, 'k')
            plt.title(r'$\alpha_g = $ {0:.2f}'.format(yaw) + '\n' + r'$\zeta_g = $ {0:.2f}'.format(pitch))
            plt.show()

    def geometryTranslation(self, translationVec, plot=False):
        '''
        Airfoil geometry translation:
            inputs:
                translationVec  -- translation vector 
                plot            -- boolean value for blade plotting 
            !!! the blade should first rotated and then translated !!!
        '''

        # variables allocation 
        # -- there was a problem with the 3D plotting (x,y,z) values 
        # -- in order to correct this problem -> raw allocation of data
        x  = np.zeros(self.upper.shape[1],)
        y  = np.zeros(self.upper.shape[1],)
        z  = np.zeros(self.upper.shape[1],)
        y1 = np.zeros(self.upper.shape[1],)
        z1 = np.zeros(self.upper.shape[1],)
        
        for ii in range(self.upper.shape[1]):
            x[ii]  = self.upper[0,ii]
            y[ii]  = self.upper[1,ii]
            z[ii]  = self.upper[2,ii]
            y1[ii] = self.lower[1,ii]
            z1[ii] = self.lower[2,ii]

        if plot:
            # plotting blade section 
            # initial blade + rotated blade
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.set_box_aspect((np.ptp(x), np.ptp(y), np.ptp(z)))
            plt.plot(self.X, self.Yupper, 0, '--r')
            plt.plot(self.x, self.Ylower, 0, '--r')
            plt.plot(x, y , z,  '--b')
            plt.plot(x, y1, z1, '--b')

        # upper surface translation 
        x = x + translationVec[0]
        y = y + translationVec[1]
        z = z + translationVec[2]

        # lower surface translation 
        y1 = y1 + translationVec[1]
        z1 = z1 + translationVec[2] 

        # data reallocation 
        self.upper = np.array([x, y, z])
        self.lower = np.array([x, y1, z1]) 

        if plot:
            # plot the translated blade
            plt.plot(x, y , z,  'k')
            plt.plot(x, y1, z1, 'k')
            plt.show()

def STLsaving(airfoils, STLname='cad'):
    '''
    This function saves the blade in .stl format.
        inputs: 
            airfoils    -- tuple of airfoils objects
            STLname     -- .stl file name
        !!! it is assumed that each airfoil has the same number of description points !!!
        !!! it is assumed that the each airfoil section element is in sequence with respect the hub !!!
    '''
    # stl generation 
    file = open(STLname + '.stl', 'w')
    file.write('solid blade\n')

    # loop over all the span positions
    for jj in range(len(airfoils)-1):

        # upper surface stl face generation 
        for ii in range(airfoils[jj].upper.shape[1]-1):
            # versor computation
            vec1 = np.array(airfoils[jj].upper[:,ii])
            vec2 = np.array(airfoils[jj].upper[:,ii+1])
            vec3 = np.array(airfoils[jj+1].upper[:,ii])
            versor = np.cross(vec3 - vec1, vec2 - vec1)
            # writing data
            versor = versor / np.linalg.norm(versor)
            file.write('\tfacet normal {0:f} {1:f} {2:f}\n'.format(versor[0], versor[1], versor[2]))
            file.write('\t\touter loop\n')
            file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec1[0], vec1[1], vec1[2]))
            file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec2[0], vec2[1], vec2[2]))
            file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec3[0], vec3[1], vec3[2]))
            file.write('\t\tendloop\n')
            file.write('\tendfacet\n')

            # versor computation
            vec1 = np.array(airfoils[jj+1].upper[:,ii])
            vec2 = np.array(airfoils[jj+1].upper[:,ii+1])
            vec3 = np.array(airfoils[jj].upper[:,ii+1])
            versor = - np.cross(vec3 - vec1, vec2 - vec1)
            # writing data
            versor = versor / np.linalg.norm(versor)
            file.write('\tfacet normal {0:f} {1:f} {2:f}\n'.format(versor[0], versor[1], versor[2]))
            file.write('\t\touter loop\n')
            file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec1[0], vec1[1], vec1[2]))
            file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec2[0], vec2[1], vec2[2]))
            file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec3[0], vec3[1], vec3[2]))
            file.write('\t\tendloop\n')
            file.write('\tendfacet\n')

        # lower surface stl face generation 
        for ii in range(naca65_hub.lower.shape[1]-1):
            # versor computation
            vec1 = np.array(airfoils[jj].lower[:,ii])
            vec2 = np.array(airfoils[jj].lower[:,ii+1])
            vec3 = np.array(airfoils[jj+1].lower[:,ii])
            versor = np.cross(vec3 - vec1, vec2 - vec1)
            # writing data
            versor = versor / np.linalg.norm(versor)
            file.write('\tfacet normal {0:f} {1:f} {2:f}\n'.format(versor[0], versor[1], versor[2]))
            file.write('\t\touter loop\n')
            file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec1[0], vec1[1], vec1[2]))
            file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec2[0], vec2[1], vec2[2]))
            file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec3[0], vec3[1], vec3[2]))
            file.write('\t\tendloop\n')
            file.write('\tendfacet\n')

            # versor computation
            vec1 = np.array(airfoils[jj+1].lower[:,ii])
            vec2 = np.array(airfoils[jj+1].lower[:,ii+1])
            vec3 = np.array(airfoils[jj].lower[:,ii+1])
            versor = - np.cross(vec3 - vec1, vec2 - vec1)
            # writing data
            versor = versor / np.linalg.norm(versor)
            file.write('\tfacet normal {0:f} {1:f} {2:f}\n'.format(versor[0], versor[1], versor[2]))
            file.write('\t\touter loop\n')
            file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec1[0], vec1[1], vec1[2]))
            file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec2[0], vec2[1], vec2[2]))
            file.write('\t\t\tvertex {0:f} {1:f} {2:f}\n'.format(vec3[0], vec3[1], vec3[2]))
            file.write('\t\tendloop\n')
            file.write('\tendfacet\n')

    file.write('endsolid blade\n')
    
    # file closure 
    file.close()

# TIP airfoil setup 
# getting airfoil data
naca65_tip = geometryData('data/airfoils/naca65.txt')

# airfoil shape 
_ = naca65_tip.geometryFitting(Cl=1.1, chord=0.75, plot=False)

# airfoil 3D rotation
naca65_tip.geometryRotation(10, 1, plot=False)

# airfoil 3D translation 
naca65_tip.geometryTranslation([0, 0, 0.2], plot=False)

# HUB airfoil setup 
# getting airfoil data 
naca65_hub = geometryData('data/airfoils/naca65.txt')

# airfoil shape 
_ = naca65_hub.geometryFitting(Cl=1.4, chord=1, plot=False)

# stl generation 
airfoils = [naca65_hub, naca65_tip]
STLsaving(airfoils, STLname='rotor1')
