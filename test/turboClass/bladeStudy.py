def bladeGenerator(kind, meanValues, b0, b1, Leu, inletValues, nSections, nBlades=0, STLname='cad', pos='/data/airfoils/naca65.txt', printout=False, plot=False):
    '''
    This function generates the blade parameters at different sections of the blade. 
        inputs:
            kind        -- string value that defines rotor or stator 
                        -- stator/rotor
            meanValues  -- array that stores main mean values.
                        -- [rMean, Umean, VaMean, VtMeanIn, VtMeanOut]
            b0          -- inlet blade height
            b1          -- outlet blade height
            Leu         -- work made at each blade section 
            inletValues -- array that stores the main inlet values for the blade
                        -- [Tt0, Pt0, T0, P0]
            nSections   -- # of section to study for the blade 
            etaVec      -- losses distribution along the blade span 
            hubChord    -- blade chord dimensions at the hub 
            nBlades     -- # of blades on the machine
            STLname     -- name for the .stl file 
            printout    -- boolean value for the printing of the main parameters 
            plot        -- boolean value for the plotting of the blade --> simple line composition plot -> for a better plot open in paraview/openscad the .stl model
    '''

    # importing libraries 
    from geometry import bladeGenerator 
    
    # mean values allocation 
    rMean     = meanValues[0]
    Umean     = meanValues[1]
    VaMean    = meanValues[2]
    VtMeanIn  = meanValues[3]
    VtMeanOut = meanValues[4]
    omega     = Umean / rMean 
    
    # inlet values allocation 
    Tt0 = inletValues[0]
    Pt0 = inletValues[1]
    T0  = inletValues[2]
    P0  = inletValues[3]

    # radial position generation 
    # vector allocation 
    r0Vec = np.zeros(nSections+1)
    r1Vec = np.zeros(nSections+1)

    # elements computation 
    for ii in range(nSections+1):
        r0Vec[ii] = rMean - b0/2 + ii * b0/nSections
        r1Vec[ii] = rMean - b1/2 + ii * b1/nSections

    if plot:
        import random
        ax = plt.axes(projection ='3d')
        number_of_colors = nSections+1
        color = ["#"+''.join([random.choice('0123456789ABCDEF') for _ in range(6)]) for _ in range(number_of_colors)]

    # generation of an object array 
    blade = [] 

    # boolean value for checkin if the studying the hub airfoil section --> needed for the span-wise discretization 
    hubPos = True 

    for ii in range(nSections+1):
        # angle computation from a FREE VORTEX model 
        _, _, _, _, _, _, angleVec, _, _ = bladeStudy(r0Vec[ii], r1Vec[ii], omega, rMean, VaMean, VtMeanIn, VtMeanOut, Leu, Tt0, T0, Pt0, P0, eta=1, printout=False, gamma=1.4, R=287.06)        
        
        # this modeling approach works with beta0 >= 0 and beta1 >= 0
        # -- in order to adopt this model to each blade modeling a boolean 
        #       value for the change is generated and then used later for the 
        #       blade geoemetry generation
        checkAngle = False  

        if kind == 'rotor':
            # angle allocation
            beta0 = angleVec[2]
            beta1 = angleVec[3]
            # chech on angle sign for adopting Lieblein model
            if beta0 < 0:
                beta0 = - beta0
                beta1 = - beta1
                checkAngle = True 

        elif kind == 'stator':
            # angle allocation 
            beta0 = angleVec[0]
            beta1 = angleVec[1]
            # chech on angle sign for adopting Lieblein model
            if beta0 < 0:
                beta0 = - beta0
                beta1 = - beta1
                checkAngle = True 

        # optimal angles computation
        i, delta, theta, solidity, tbc = optimalAngles(beta0, beta1, printout=False)

        # computing optimal alpha angle
        ac = 0.5 # this is valid only for NACA-65
        alpha = alphaFunc(ac, solidity, theta, tbc)

        # computing stagger angle gamma 
        gamma = beta0 - alpha

        # checking angles sign 
        if checkAngle:
            i = - i 
            delta = - delta 
            theta = - theta
            alpha = - alpha 
            gamma = - gamma 

        # computing Cl 
        Cl = ac * np.tan(np.deg2rad(theta)/4) / 0.0551515

        # pitch computation with nBlades 
        pitch = 2 * np.pi * r0Vec[ii] / nBlades 

        # chord computation from solidity
        chord = pitch * solidity 

        # computing blade inclination 
        #zeta = np.rad2deg(np.arctan2(r1Vec[ii]-r0Vec[ii], chord))
        zeta = np.rad2deg(np.arcsin((r1Vec[ii]-r0Vec[ii])/chord))

        # airfoil object generation 
        airfoil = bladeGenerator.geometryData(pos)

        # geometry fitting airfoil -> shape and chord dimensions
        airfoil.geometryFitting(Cl=Cl, chord=chord, plot=False)

        # airfoil 3D rotation
        airfoil.geometryRotation(gamma, zeta, plot=False)

        # translation of the airfoil with respect to the hub position 
        if hubPos: 
            airfoilHub = airfoil 
            hubPos = False 
            translationHub = airfoilHub.middleChord()
            translationHub[2] = 0.0
        
        # translation of profiles
        # translation vector 
        height = r0Vec[ii] - r0Vec[0]
        # airfoil section center
        translationSection = airfoil.middleChord()
        translationSection[0] = 0 
        translationSection[1] = 0
        # translation vector  
        translationVec = translationHub + translationSection
        # translation procedure 
        airfoil.geometryTranslation(translationVec, height, plot=False)

        if plot:
            ax.plot3D(airfoil.upper[:,0], airfoil.upper[:,1], airfoil.upper[:,2], color=color[ii], label=str(ii))
            ax.plot3D(airfoil.lower[:,0], airfoil.lower[:,1], airfoil.lower[:,2], color=color[ii])

        # adding airfoil to blade object array
        blade.append(airfoil)

        # printout 
        if printout:
            starDim = 50
            bladeDim = (starDim - len(' BLADE ANGLES '))/2 
            print('*' * bladeDim + ' BLADE ANGLES ' + '*' * bladeDim)
            print('r             = {0:.3f}'.format(r0Vec[ii]))
            print('i             = {0:.3f}'.format(i))
            print('delta         = {0:.3f}'.format(delta))
            print('alpha         = {0:.3f}'.format(alpha))
            print('beta0 - beta1 = {0:.3f}'.format(beta0 - beta1))
            print('gamma         = {0:.3f}'.format(gamma))
            print('theta         = {0:.3f}'.format(theta))
            print('Cl            = {0:.3f}'.format(np.abs(Cl)))
            print('s             = {0:.3f}'.format(pitch))
            print('c             = {0:.3f}'.format(chord))
            print('zeta          = {0:.3f}\n'.format(zeta))
            print('*' * starDim + '\n')

    if plot:
        if nSections < 10:
            plt.legend()
        plt.title('Blade')
        plt.show()

    # STL file generation 
    bladeGenerator.STLsaving(blade, STLname=STLname)

    # setting up return values 
    inlet = [blade[0].chord, blade[0].camber[0,0], blade[0].camber[0,1], blade[0].camber[0,2], blade[0].camber[-1,0], blade[0].camber[-1,1], blade[0].camber[-1,2]]
    outlet = [blade[-1].chord, blade[-1].upper[0,0], blade[-1].upper[0,1], blade[-1].upper[0,2], blade[-1].upper[-1,0], blade[-1].upper[-1,1], blade[-1].upper[-1,2]]

    return inlet, outlet 