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
        self.x     = np.array(data[:,0] / 100)
        self.y     = np.array(data[:,1] / 100)
        self.z     = np.zeros(self.x.shape)
        self.t     = np.array(data[:,2] / 100)
        self.chord = np.max(self.x) - np.min(self.x)

    def geometryFitting(self, Cl, chord=1, plot=False):
        '''
        Airfoil geometry fitting, it is used for geometry modification with respect to the given Cl.
            inputs: 
                Cl      -- blade lift coefficient
                chord   -- airfoil chord
                plot    -- showing airfoil coordinates
        '''

        # chord dimension allocation 
        self.chord = chord 

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

        # data allocation -- points 
        self.X       = self.x 
        self.Yupper  = yUpper
        self.Ylower  = yLower
        self.Zupper  = self.z 
        self.Zlower  = self.z 
        self.T       = t 
        self.Ycamber = yCamber
        self.Zcamber = self.z

        # vector allocation 
        self.upper  = np.zeros((self.X.shape[0], 3))
        self.lower  = np.zeros(self.upper.shape)
        self.camber = np.zeros(self.upper.shape)

        # data allocation -- vectors 
        for ii in range(self.upper.shape[0]):
            self.upper[ii,0]  = self.x[ii]
            self.upper[ii,1]  = yUpper[ii]
            self.upper[ii,2]  = 0
            self.lower[ii,0]  = self.x[ii]
            self.lower[ii,1]  = yLower[ii]
            self.lower[ii,2]  = 0 
            self.camber[ii,0] = self.x[ii]
            self.camber[ii,1] = yCamber[ii]
            self.camber[ii,2] = 0

        if plot:
            plt.figure(figsize=(9,9))
            plt.plot(self.camber[:,0], self.camber[:,1], 'r')
            plt.plot(self.upper[:,0], self.upper[:,1], 'b', linewidth=2)
            plt.plot(self.lower[:,0], self.lower[:,1], 'k', linewidth=2)
            plt.axis('equal')
            plt.show()

    def geometryRotation(self, yaw, pitch=0, plot=False):
        '''
        Airfoil geometry rotation:
            -- yaw      -> [deg] blade to blade plane -- metal angle  
            -- pitch    -> [deg] meridional plane     -- due to annulus diameter variation
            -- plot     -- boolean value for the blade plotting
        '''

        # saving rotation properties 
        self.yaw = yaw 
        self.pitch = pitch 

        # alpha == roll == 0
        roll = 0
        cA = np.cos(np.deg2rad(roll))
        sA = np.sin(np.deg2rad(roll))
        # beta == pitch --> blade angle in meridional plane
        cB = np.cos(np.deg2rad(- pitch))
        sB = np.sin(np.deg2rad(- pitch))
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
        for ii in range(self.upper.shape[0]):
            self.upper[ii,:] = np.matmul(rotMatrix, self.upper[ii,:])

        # rotation lower surface 
        for ii in range(self.lower.shape[0]):
            self.lower[ii,:] = np.matmul(rotMatrix, self.lower[ii,:])

        # rotation camber line 
        for ii in range(self.camber.shape[0]):
            self.camber[ii,:] = np.matmul(rotMatrix, self.camber[ii,:])       

        if plot:
            # plotting blade section 
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            if pitch != 0:
                ax.set_box_aspect((np.ptp(self.upper[:,0]), np.ptp(self.upper[:,1]), np.ptp(self.upper[:,2])))
            plt.plot(self.X, self.Yupper, 0, '--r', label='not rotated')
            plt.plot(self.X, self.Ylower, 0, '--r')
            plt.plot(self.X, self.Ycamber, 0, '-*r')
            plt.plot(self.upper[:,0], self.upper[:,1], self.upper[:,2], 'k', label='rotated')
            plt.plot(self.lower[:,0], self.lower[:,1], self.lower[:,2], 'k')
            plt.plot(self.camber[:,0], self.camber[:,1], self.camber[:,2], '*k')
            plt.title(r'$\gamma = {0:.2f}^{{\circ}}$'.format(yaw) + '\n' + r'$\zeta = {0:.2f}^{{\circ}}$'.format(pitch))
            plt.legend(loc='best')
            plt.show()

    def middleChord(self, printout=False):
        '''
        This function computes the middle chord point of the airfoil.
            output:
                middlePoint -- airfoil chord middle point
                printout    -- boolean value for the print of the middle chord position 
        '''

        # middle point computation 
        self.middlePoint = (self.camber[0,:] + self.camber[-1,:]) / 2 

        if printout:
            print('MIDDLE CHORD coordinates:')
            print('x = {0:.2f} m\ty = {1:.2f} m\tz = {2:.2f}'.format(self.middlePoint[0], self.middlePoint[1], self.middlePoint[2]))

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

        if plot:
            # plotting blade section 
            # initial blade + rotated blade
            fig = plt.figure(figsize=(9,9))
            ax = fig.add_subplot(111, projection='3d')
            ax.set_title(r'$v = [ {0:.2f}, {1:.2f}, {2:.2f}]$'.format(translationVec[0], translationVec[1], translationVec[2]))
            try:
                if self.pitch != 0:
                    ax.set_box_aspect((np.ptp(self.upper[:,0]), np.ptp(self.upper[:,1]), np.ptp(self.upper[:,2])))
            except:
                pass
            plt.plot(self.X, self.Yupper, 0, '--r', label='not rotated')
            plt.plot(self.X, self.Ylower, 0, '--r')
            plt.plot(self.X, self.Ycamber, 0, '-*r')
            plt.plot(self.upper[:,0], self.upper[:,1], self.upper[:,2], 'b', label='rotated')
            plt.plot(self.lower[:,0], self.lower[:,1], self.lower[:,2], 'b')
            plt.plot(self.camber[:,0], self.camber[:,1], self.camber[:,2], '*b')

        # upper surface translation 
        for ii in range(self.upper.shape[0]):
            self.upper[ii,0] = self.upper[ii,0] + translationVec[0]
            self.upper[ii,1] = self.upper[ii,1] + translationVec[1]
            self.upper[ii,2] = self.upper[ii,2] + translationVec[2]

        # lower surface translation 
        for ii in range(self.lower.shape[0]):
            self.lower[ii,0] = self.lower[ii,0] + translationVec[0]
            self.lower[ii,1] = self.lower[ii,1] + translationVec[1]
            self.lower[ii,2] = self.lower[ii,2] + translationVec[2]


        # camber line translation
        for ii in range(self.camber.shape[0]):
            self.camber[ii,0] = self.camber[ii,0] + translationVec[0]
            self.camber[ii,1] = self.camber[ii,1] + translationVec[1]
            self.camber[ii,2] = self.camber[ii,2] + translationVec[2]

        if plot:
            # plot the translated blade
            plt.plot(self.upper[:,0], self.upper[:,1], self.upper[:,2], 'k', label='rotated + translated')
            plt.plot(self.lower[:,0], self.lower[:,1], self.lower[:,2], 'k')
            plt.plot(self.camber[:,0], self.camber[:,1], self.camber[:,2], '*k')
            plt.legend()
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

def STLsaving(airfoils, STLname='cad', containerPath='container/', kind='rotor', checkVersor=[False,False]):
    '''
    This function saves the blade in .stl format.
        inputs: 
            airfoils    -- tuple of airfoils objects
            STLname     -- .stl file name
            kind        -- rotor/stator => allows computing the correct direction of the versors
            checkVersor -- boolean array value for the plotting/checking of versor
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
        for ii in range(airfoils[jj].upper.shape[0]-1):
            # versor computation
            vec1 = np.array(airfoils[jj].upper[ii,:])
            vec2 = np.array(airfoils[jj].upper[ii+1,:])
            vec3 = np.array(airfoils[jj+1].upper[ii,:])
            versor = - np.cross(vec3 - vec1, vec2 - vec1)
            # vector normalization
            versor = versor / np.linalg.norm(versor)
            # versor analysis
            if kind == 'stator':
                versor = - versor
            # writing data
            writeFacet(file, versor, vec1, vec2, vec3)

            # versor computation
            vec1 = np.array(airfoils[jj+1].upper[ii,:])
            vec2 = np.array(airfoils[jj+1].upper[ii+1,:])
            vec3 = np.array(airfoils[jj].upper[ii+1,:])
            versor = np.cross(vec3 - vec1, vec2 - vec1)
            # vector normalization 
            versor = versor / np.linalg.norm(versor)
            # versor analysis
            if kind == 'stator':
                versor = - versor
            # writing data
            writeFacet(file, versor, vec1, vec2, vec3)

        if checkVersor[0]:
            plt.plot(airfoils[0].lower[:,0], airfoils[0].lower[:,1], 'r')
            plt.plot(airfoils[0].upper[:,0], airfoils[0].upper[:,1], 'b')
            plt.plot([0,versor[0]], [0,versor[1]], 'b')

        # lower surface stl face generation 
        for ii in range(airfoils[jj].lower.shape[0]-1):
            # versor computation
            vec1 = np.array(airfoils[jj].lower[ii,:])
            vec2 = np.array(airfoils[jj].lower[ii+1,:])
            vec3 = np.array(airfoils[jj+1].lower[ii,:])
            versor = np.cross(vec3 - vec1, vec2 - vec1)
            # direction check -- for the upper surface the versor direction is towards y < 0 
            if versor[1] < 0: 
                versor = - versor 
            # vectorn normalization
            versor = versor / np.linalg.norm(versor)
            # versor analysis
            if kind == 'stator':
                versor = - versor
            # writing data
            writeFacet(file, versor, vec1, vec2, vec3)

            # versor computation
            vec1 = np.array(airfoils[jj+1].lower[ii,:])
            vec2 = np.array(airfoils[jj+1].lower[ii+1,:])
            vec3 = np.array(airfoils[jj].lower[ii+1,:])
            versor = - np.cross(vec3 - vec1, vec2 - vec1)
            # direction check -- for the upper surface the versor direction is towards y < 0 
            if versor[1] < 0: 
                versor = - versor 
            # vector normalization
            versor = versor / np.linalg.norm(versor)
            # versor analysis
            if kind == 'stator':
                versor = - versor
            # writing data
            writeFacet(file, versor, vec1, vec2, vec3)
        
        if checkVersor[1]:
            plt.plot(airfoils[0].lower[:,0], airfoils[0].lower[:,1], 'r')
            plt.plot(airfoils[0].upper[:,0], airfoils[0].upper[:,1], 'b')
            plt.plot([0,versor[0]], [0,versor[1]], 'r')
            plt.show()
    
    # end blade 
    file.write('endsolid bladeSurface\n')

    # begin bottom blade surface -- hub
    file.write('solid bladeBottom\n')

    # closing airfoil [hub]
    for ii in range(airfoils[0].upper.shape[0]-2):
        # versor computation 
        vec1 = np.array(airfoils[0].camber[ii+1,:])
        vec2 = np.array(airfoils[0].upper[ii,:])
        vec3 = np.array(airfoils[0].upper[ii+1,:])
        # versor computation 
        versor = [0, 0, -1]
        # writing data
        writeFacet(file, versor, vec1, vec2, vec3)

    for ii in range(airfoils[0].lower.shape[0]-2):
        # versor computation 
        vec1 = np.array(airfoils[0].camber[ii+1,:])
        vec2 = np.array(airfoils[0].lower[ii,:])
        vec3 = np.array(airfoils[0].lower[ii+1,:])
        # versor computation 
        versor = [0, 0, -1]
        # writing data
        writeFacet(file, versor, vec1, vec2, vec3)

    for ii in range(airfoils[0].upper.shape[0]-2):
        # versor computation 
        vec1 = np.array(airfoils[0].camber[ii+2,:])
        vec2 = np.array(airfoils[0].upper[ii+1,:])
        vec3 = np.array(airfoils[0].camber[ii+1,:])
        # versor computation 
        versor = [0, 0, -1]
        # writing data
        writeFacet(file, versor, vec1, vec2, vec3)
    
    for ii in range(airfoils[0].upper.shape[0]-2):
        # versor computation 
        vec1 = np.array(airfoils[0].camber[ii+2,:])
        vec2 = np.array(airfoils[0].lower[ii+1,:])
        vec3 = np.array(airfoils[0].camber[ii+1,:])
        # versor computation 
        versor = [0, 0, -1]
        # writing data
        writeFacet(file, versor, vec1, vec2, vec3)

    # end blade bottom surface -- hub 
    file.write('endsolid bladeBottom\n')

    # begin top blade surface -- tip
    file.write('solid bladeTop\n')

    # closing airfoil [tip]
    for ii in range(airfoils[-1].upper.shape[0]-2):
        # versor computation 
        vec1 = np.array(airfoils[-1].camber[ii+1,:])
        vec2 = np.array(airfoils[-1].upper[ii,:])
        vec3 = np.array(airfoils[-1].upper[ii+1,:])
        # versor computation 
        versor = [0, 0, 1]
        # writing data
        writeFacet(file, versor, vec1, vec2, vec3)

    for ii in range(airfoils[-1].lower.shape[0]-2):
        # versor computation 
        vec1 = np.array(airfoils[-1].camber[ii+1,:])
        vec2 = np.array(airfoils[-1].lower[ii,:])
        vec3 = np.array(airfoils[-1].lower[ii+1,:])
        # versor computation 
        versor = [0, 0, 1]
        # writing data
        writeFacet(file, versor, vec1, vec2, vec3)

    for ii in range(airfoils[-1].upper.shape[0]-2):
        # versor computation 
        vec1 = np.array(airfoils[-1].camber[ii+2,:])
        vec2 = np.array(airfoils[-1].upper[ii+1,:])
        vec3 = np.array(airfoils[-1].camber[ii+1,:])
        # versor computation 
        versor = [0, 0, 1]
        # writing data
        writeFacet(file, versor, vec1, vec2, vec3)
    
    for ii in range(airfoils[-1].upper.shape[0]-2):
        # versor computation 
        vec1 = np.array(airfoils[-1].camber[ii+2,:])
        vec2 = np.array(airfoils[-1].lower[ii+1,:])
        vec3 = np.array(airfoils[-1].camber[ii+1,:])
        # versor computation 
        versor = [0, 0, 1]
        # writing data
        writeFacet(file, versor, vec1, vec2, vec3)
    
    # end blade top surface -- tip 
    file.write('endsolid bladeTop\n')

    # file closure 
    file.close()

def SCADsaving(nRotorBlades, nStatorBlades, rotorHub, statorHub, rMean, b0, b1, rotorPath='../container/rotor.stl', statorPath='../container/stator.stl', geometryPath='../geometry/'):
    '''
    This function saves the computed results into .scad file that will use for the visualization of the stage.
        inputs:
            nRotorBlades    -- # of blades for the rotor
            nStatorBlades   -- # of blades for the stator
            rotorHub        -- array that stores quantities necessary for the printing    
                            -- [rotor blade hub x0 camber, rotor blade hub y0 camber , rotor blade hub z0 camber]
            statorHub       -- array that stores quantities necessary for the printing    
                            -- [stator blade hub x0 camber, stator blade hub y0 camber , stator blade hub z0 camber]
            rMean           -- mean radius 
                            -- it is the same for rotor and stator due to modeling assumptions
            b0              -- rotor blade inlet blade height
            b1              -- stator blade inlet blade height
            rotorPath       -- path where the rotor .stl file is saved 
            statorPath      -- path wehre the stator .stil file is saved
            geometryPath    -- directory where the file should be stores/saved
    '''

    # stl generation 
    file = open(geometryPath + 'data.scad', 'w')

    # print data for hub generation
    file.write('/*\n')
    file.write('TURBOMACHINERY -- LIBRARY FOR THE INITIAL TURBOMACHINERY DESIGN\n')
    file.write('AUTHOR: antonio pucciarelli\n')
    file.write('PROGRAM DESCRIPTION\n')
    file.write('    TURBOMACHINERY 3D GEOMETRY DATA FILE\n')
    file.write('*/\n\n')
    file.write('// VARIABLES ALLOCATION\n')
    file.write('// rotor properties\n')
    file.write('nRotorBlades         = {0:d};\n'.format(nRotorBlades))
    file.write('rotorHubInletCoords  = [{0}, {1}, {2}];\n'.format(rotorHub[1], rotorHub[2], rotorHub[3] + rMean - b0/2))
    file.write('rotorHubOutletCoords = [{0}, {1}, {2}];\n'.format(rotorHub[4], rotorHub[5], rotorHub[6] + rMean - b0/2))
    file.write('rotorName            = "{0}";\n'.format(rotorPath))
    file.write('// stator properties\n')
    file.write('nStatorBlades         = {0:d};\n'.format(nStatorBlades))
    file.write('statorHubInletCoords  = [{0}, {1}, {2}];\n'.format(statorHub[1], statorHub[2], statorHub[3] + rMean - b1/2))
    file.write('statorHubOutletCoords = [{0}, {1}, {2}];\n'.format(statorHub[4], statorHub[5], statorHub[6] + rMean - b1/2))
    file.write('statorName            = "{0}";\n'.format(statorPath))

    # closing file
    file.close()
