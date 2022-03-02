# TURBOMACHINERY -- LIBRARY FOR THE VELOCITY TRIANGLES ANALYSIS AND STAGES PLOT
# AUTHOR: antonio pucciarelli 
#
# PROGRAM DESCRIPTION
#   TURBOMACHINERY BLADES AND VELOCITY: computation and plot
#       -- airfoil rotation subroutine
#       -- turbomachinery stages plot 
#       -- velocity triangles plot 
#

# importing libraries
import matplotlib.pyplot as plt 
import numpy as np

def airfoilRotation(foil, gamma):
    '''
    Airfoil rotation function

        function inputs:
            foil        -- airfoil coordinates data
                        -- expressed in [x, y] numpy array format
            gamma       -- airfoil rotation angle
                        -- in radiant
    '''

    # profile rotation due to stagger angle
    rotMatrix = np.matrix([[np.cos(gamma), -np.sin(gamma)],[np.sin(gamma), np.cos(gamma)]])

    # coordinate rotation
    for ii in range(foil.shape[0]):
        foil[ii,:] = np.matmul(rotMatrix, foil[ii,:])

    return foil

def turboPlot(nStage=1,rotorName=['naca5016'],rotorAngle=[0],rotorChord=[1],rotorOrig=[0],rotorPitch=[0.5],
                       statorName=['naca1212'],statorAngle=[1],statorChord=[1],statorOrig=[1.5],statorPitch=[0.5],plot=True):
    '''
    Turbomachinery stage(s) plot

        function inputs:
            nStage          -- number of stage for the turbomachinery 
            rotorName       -- airfoil name
                            -- it should be created through xfoil or programs like this 
            rotorAngle      -- alpha: airfoil incidence with respect to the axial line 
                            -- in degree
            rotorChord      -- chord: airfoil scaling through chord dimension
            rotorOrigin     -- origin: airfoil leading edge position 
            rotorPitch      -- s: airfoil pitch
                            -- if you want to print just 1 airfoil instead of 3 set it as 0
            statorName      -- airfoil name
                            -- it should be created through xfoil or programs like this 
            statorAngle     -- alpha: airfoil incidence with respect to the axial line 
                            -- in degree
            statorChord     -- chord: airfoil scaling through chord dimension
            statorOrigin    -- origin: airfoil leading edge position 
            statorPitch     -- s: airfoil pitch
                            -- if you want to print just 1 airfoil instead of 3 set it as 0
    '''

    # angle conversion in radiant 
    rotorAngle = np.deg2rad(rotorAngle)
    statorAngle = np.deg2rad(statorAngle)

    # plotting turbomachinery
    if plot:
        # figure generation
        plt.figure()

        # plotting each stage:
        for ii in range(nStage):
            # rotor
            # loading rotor data
            rotorFoil = np.loadtxt(rotorName[ii],skiprows=1)
            # rotor airfoil scaling
            rotorFoil = rotorFoil * rotorChord[ii]
            # rotor airfoil rotation
            rotorFoil = airfoilRotation(rotorFoil, rotorAngle[ii])
            # rotor airfoil translation 
            rotorFoil[:,0] += rotorOrig[ii]
            
            # stator
            # loading stator data
            statorFoil = np.loadtxt(statorName[ii],skiprows=1)
            # stator airfoil scaling
            statorFoil = statorFoil * statorChord[ii]
            # stator airfoil rotation
            statorFoil = airfoilRotation(statorFoil, statorAngle[ii])
            # stator airfoil translation 
            statorFoil[:,0] += statorOrig[ii]

            # plotting rotor
            plt.plot(rotorFoil[:,0],rotorFoil[:,1],'r',label='Rotor #{0:d}\n{1:s}\n{2:2.2f}$^\circ$'.format(ii+1, rotorName[ii], np.rad2deg(rotorAngle[ii])))
            # plotting stator
            plt.plot(statorFoil[:,0],statorFoil[:,1],'b',label='Stator #{0:d}\n{1:s}\n{2:2.2f}$^\circ$'.format(ii+1, statorName[ii], np.rad2deg(statorAngle[ii])))

            # plotting other blades
            if rotorPitch[ii] != 0:
                plt.plot(rotorFoil[:,0],rotorFoil[:,1]+rotorPitch[ii],'r')
                plt.plot(rotorFoil[:,0],rotorFoil[:,1]-rotorPitch[ii],'r')
            if statorPitch[ii] != 0:
                plt.plot(statorFoil[:,0],statorFoil[:,1]+statorPitch[ii],'b')
                plt.plot(statorFoil[:,0],statorFoil[:,1]-statorPitch[ii],'b')

    # important quantities printout
    print('+++++++++++++++++++++++++++++++++++++++')
    print('Stages analysis:')
    for ii in range(nStage):
        print('stage #{0}'.format(nStage))
        print('-- rotor')
        print('\tairfoil name: {0:>10s}'.format(rotorName[ii]))
        print('\tangle:        {0:>10.2f}'.format(np.rad2deg(rotorAngle[ii])))
        print('\tchord:        {0:>10.2f}'.format(rotorChord[ii]))
        print('\tpitch:        {0:>10.2f}'.format(rotorPitch[ii]))
        print('\tLE position:  {0:>10.2f}'.format(rotorOrig[ii]))
        print('-- stator')
        print('\tairfoil name: {0:>10s}'.format(statorName[ii]))
        print('\tangle:        {0:>10.2f}'.format(np.rad2deg(statorAngle[ii])))
        print('\tchord:        {0:>10.2f}'.format(statorChord[ii]))
        print('\tpitch:        {0:>10.2f}'.format(statorPitch[ii]))
        print('\tLE position:  {0:>10.2f}'.format(statorOrig[ii]))
    print('+++++++++++++++++++++++++++++++++++++++')

    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title='Airfoils')
    plt.axis('equal')
    plt.tight_layout()
    plt.show()            

def velocityVec(airfoilData=['naca1212'],omega=0,D=0,V1=0,V2=0,alpha1=0,alpha2=0,gamma=0,s=0,c=1,velScale=0.5,plot=True):
    '''
    Velocity triangles plot over a turbomachinery blade

        function inputs:
            airfoil     -- airfoil data in .dat format
            omega       -- rotor velocity 
                        -- if it is a stator omega = 0
            D           -- radial distance  
            V1          -- absolute inlet velocity magnitude
            V2          -- absolute outlet velocity magnitude 
            alpha1      -- inlet flow vs blade angle  
                        -- degree 
            alpha2      -- outlet flow vs blade angle
                        -- degree 
            epsilon     -- flow deflection 
            gamma       -- stagger angle 
            s           -- pitch 
                        -- distance between turbomachinery blades 
            c           -- chord 
    '''

    # converting angles into radiants
    alpha1 = np.deg2rad(alpha1)
    alpha2 = np.deg2rad(alpha2)
    gamma = np.deg2rad(gamma)

    # importing geometry
    foil = np.loadtxt(airfoilData, skiprows=1)

    # scaling the airfoil with respect to the input chord
    foil = foil * c

    # airfoil rotation with respect to gamma 
    foil = airfoilRotation(foil, gamma)

    # solidity computation
    if s != 0:
        solidity = c / s
    else: 
        solidity = np.NaN
      
    # normalized and scaled V1 component computation with respect to V1
    V1a = V1 * np.cos(alpha1)  
    V1t = V1 * np.sin(alpha1) 

    # normalized and scaled V2 component computation with respect to V1
    V2a = V2 * np.cos(alpha2)  
    V2t = V2 * np.sin(alpha2) 

    # rotor velocity for an axial non centrifugal/centripetal turbomachine 
    if omega != 0:
        kind = 'rotor'
        U = omega * D/2 
    else:
        kind = 'stator'
        U = 0

    # inlet relative velocity computation 
    W1 = np.sqrt(V1a**2 + (U - V1t)**2)

    # outlet relative velocity computation 
    W2 = np.sqrt(V2a**2 + (U - V2t)**2)

    # important quantities printout
    print('+++++++++++++++++++++++++++++++++++++++')
    print('Blade analysis:')
    print('-- type: {0:s}'.format(kind))
    print('-- solidity = {0:4.2f}'.format(solidity))
    print('-- Inlet:')
    print('\tabsolute   -- V1  = {0:>7.2f} m/s'.format(V1))
    print('\t-- axial   -- V1a = {0:>7.2f} m/s'.format(V1a))
    print('\t-- tangent -- V1t = {0:>7.2f} m/s'.format(V1t))
    print('\thub        -- U1  = {0:>7.2f} m/s'.format(U))
    print('\trelative   -- W1  = {0:>7.2f} m/s'.format(W1))
    print('-- Outlet:')
    print('\tabsolute   -- V2  = {0:>7.2f} m/s'.format(V2))
    print('\t-- axial   -- V2a = {0:>7.2f} m/s'.format(V2a))
    print('\t-- tangent -- V2t = {0:>7.2f} m/s'.format(V2t))
    print('\thub        -- U2  = {0:>7.2f} m/s'.format(U))
    print('\trelative   -- W2  = {0:>7.2f} m/s'.format(W2))
    print('+++++++++++++++++++++++++++++++++++++++')

    # computation of the V1 vector application point
    O1 = ( foil[int(foil.shape[0]/2),:] - np.array([V1a, V1t]) * velScale ) * 1.3
 
    # computation of the V2 vector application point 
    O2 = foil[0,:] * 1.3

    if plot:
        plt.figure()

        # airfoil plot 
        plt.plot(foil[:,0], foil[:,1], 'k')

        # plotting other blades 
        if s != 0: 
            plt.plot(foil[:,0], foil[:,1] + s, 'k')
            plt.plot(foil[:,0], foil[:,1] - s, 'k')

        # inlet plot sclaled with respect to velScale
        # absolute velocity plot 
        abs1 = plt.quiver(O1[0], O1[1], V1a * velScale, V1t * velScale, angles='xy', scale_units='xy', scale=1, color='r', label=r'$V_1$' + ' = {0:>6.2f} '.format(V1) + r'$\frac{m}{s}$')
        # rotor velocity
        rot1 = plt.quiver(O1[0] + V1a * velScale, O1[1] + (V1t - U) * velScale, 0, U * velScale, angles='xy', scale_units='xy', scale=1, color='m', label=r'$U_1$' + ' = {0:>6.2f} '.format(U) + r'$\frac{m}{s}$')
        # relative velocity to the rotor
        rel1 = plt.quiver(O1[0], O1[1], V1a * velScale, (V1t - U) * velScale, angles='xy', scale_units='xy', scale=1, color='g', label=r'$W_1$' + ' = {0:>6.2f} '.format(W1) + r'$\frac{m}{s}$')

        # legend for the inlet velocity 
        leg1 = plt.legend(handles=[abs1, rot1, rel1], loc='upper left', bbox_to_anchor=(1, 1), title='Inlet')

        # outlet plot scaled with respect ot velScale
        # absolute velocity plot 
        abs2 = plt.quiver(O2[0], O2[1], V2a * velScale, V2t * velScale, angles='xy', scale_units='xy', scale=1, color='b', label=r'$V_2$' + ' = {0:>6.2f} '.format(V2) + r'$\frac{m}{s}$')
        # rotor velocity 
        rot2 = plt.quiver(O2[0] + V2a * velScale, O2[1] + (V2t - U) * velScale, 0, U * velScale, angles='xy', scale_units='xy', scale=1, color='m', label=r'$U_2$' + ' = {0:>6.2f} '.format(U) + r'$\frac{m}{s}$')
        # relative velocity to the rotor
        rel2 = plt.quiver(O2[0], O2[1], V2a * velScale, (V2t - U) * velScale, angles='xy', scale_units='xy', scale=1, color='c', label=r'$W_2$' + ' = {0:>6.2f} '.format(W2) + r'$\frac{m}{s}$')

        # chord plot 
        plt.plot([foil[0,0], foil[int(foil.shape[0]/2),0]], [foil[0,1], foil[int(foil.shape[0]/2),1]], '--r')

        # legend for the outlet velocity 
        leg2 = plt.legend(handles=[abs2, rot2, rel2], loc='lower left', bbox_to_anchor=(1, 0), title='Outlet')

        # other turbomachinery properites legend
        if s != 0:
            plt.legend([r'$\frac{c}{s}$' + ' = {:2.2f}'.format(solidity),
                        r'$\alpha_1$' + ' = {:2.2f}'.format(np.rad2deg(alpha1)) + r'$^\circ$',
                        r'$\alpha_2$' + ' = {:2.2f}'.format(np.rad2deg(alpha2)) + r'$^\circ$',
                        r'$\gamma$' + ' = {:2.2f}'.format(np.rad2deg(gamma)) + r'$^\circ$'], loc='center left', bbox_to_anchor=(1, 0.5), title='Other properties')
            plt.gca().add_artist(leg1)
            plt.gca().add_artist(leg2)
        else:
            plt.gca().add_artist(leg1)

        plt.axis('equal')
        plt.xlim(O1[0]*1.05,(O2[0]+V2a*velScale)*1.05)
        plt.tight_layout()
        plt.show()

    return foil, V1a, V1t, V2a, V2t, U, W1, W2
