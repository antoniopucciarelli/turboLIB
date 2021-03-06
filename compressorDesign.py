# importing libraries
from turboClass import turboBlade
from turboCoeff import similarity
from turboCoeff import coeff
from geometry import bladeGenerator
import contextlib

# data
# constraints
mFlux = 100     # mass flux                [kg/s]
betaP = 1.45    # compression ratio        [--]
maxD  = 0.9     # maximum tip diameter     [m]
maxR  = maxD/2  # maximum tip radius       [m]

# inlet values
Pt0 = 1e+5      # inlet total pressure     [Pa]
Tt0 = 300       # inlet total temperature  [K]

# stage hypothesis
# reaction degree
rD = 0.58
# stage mean radius -> radius @ inlet blade midspan
rMean = 0.32
# rotor inlet tangential velocity
Vt0Umean = 0

# Vt0 = (1 - rD - lam/4) * Umean -> lam = 4 * (1 - rD - Vt0/Umean) 
# psi = lam / 2 
lam = (1 - rD - Vt0Umean) * 4 
psiTarget = lam / 2

# declaring blades
# rotor description
nRotorBlades  = 45
ARrotor       = 2.1
# stator description 
nStatorBlades = 35
ARstator      = 2.2
nSection      = 50

# output file generation
file_path = 'compressor_' + str(rD) + '_' + str(rMean) + '_' + str(nRotorBlades) + '_' + str(nStatorBlades) + '.txt'
with open(file_path, "w") as file:
    with contextlib.redirect_stdout(file):
        # generation of mean line properties to be used for the blade assembly 
        adimVec, bladeVec, rotationVec, V0vec, V1vec, V2vec, W0vec, W1vec, W2vec, thermo0, thermo1, _, work = similarity.stageProperties(rD, psiTarget, rMean, mFlux, Tt0, Pt0, betaP, T1real=True, printout=True, save=True, )

        # values allocation 
        Leu                = work[0]
        Lis                = work[1]
        eta                = adimVec[2]
        # rotation 
        omega              = rotationVec[1] * 1.05 # !!! the rotor angular velocity is changed slightly in order to have better rotor outlet velocity distribution
        # blade dimension
        b0                 = bladeVec[0]
        hubRadius          = rMean - b0/2
        # velocity vectors
        rotorVaMeanInlet   = V0vec[0]
        rotorVtMeanInlet   = V0vec[1] 
        rotorVaMeanOutlet  = V1vec[0]
        rotorVtMeanOutlet = Lis/0.82 / (rMean * omega) + rotorVtMeanInlet # !!! after a previous check the eta computed by similarity.stageProperties doesn't fit the actual efficiency of the stage
        statorVtMeanOutlet = rotorVtMeanOutlet 
        statorVtMeanInlet  = V1vec[1]
        statorVaMeanOutlet = V2vec[0]
        statorVtMeanOutlet = V2vec[1]

        # computing mixed vortex design parameter
        WtTarget = 0
        bVal = (WtTarget + omega * hubRadius - rotorVtMeanOutlet * rMean/hubRadius)/(hubRadius - rMean**2/hubRadius)

        # rotor study 
        # rotor object generation 
        print('\n\n-- ROTOR STUDY -- # blades {0:d}'.format(nRotorBlades))
        rotorBlade = turboBlade.blade(ID=1, turboType='rotor', nSection=nSection, inletBladeHeight=b0, outletBladeHeight=b0, inletHubRadius=hubRadius, outletHubRadius=hubRadius, omega=omega, nBlade=nRotorBlades)
        # allocate blade shape 
        rotorBlade.allocateShape(bladeHeight=b0, AR=ARrotor, nBlade=nRotorBlades)
        # blade dimensions allocation -> kinetics inlet
        rotorBlade.allocateKinetics(rMean=rMean, VtMean=rotorVtMeanInlet, VaMean=rotorVaMeanInlet, omega=omega, section='inlet', kind='FV')
        # blade dimensions allocation -> kinetics outlet
        rotorBlade.allocateKinetics(rMean2=rMean, Vt2=rotorVtMeanOutlet, VaMean=rotorVaMeanOutlet, omega=omega, section='outlet', kind='MVD', b=bVal)
        # blade dimensions allocation -> thermodynamics inlet/outlet
        rotorBlade.allocateThermodynamics(Tt0=Tt0, Pt0=Pt0, eta=eta)
        # rotor blade geometry allocation
        rotorBlade.generateGeometry(pos='data/airfoils/naca65.txt', STLname='rotor', plot=False, printout=False)
        # computing the best shape 
        lossVec = rotorBlade.bladeGenerator(mFlux, clearance=1e-3, NISRE=True, STLname='rotor', plot=False, nMaxShape=1)
        # plotting meridional flow 
        rotorBlade.printMeridional(save=True, position0='latex/figures/rotorEntropyFlow.pdf', position1='latex/figures/rotorBetaThermo.pdf')
        # plotting velocity triangles
        rotorBlade.velocityTriangles(sectionNumber=[0, int(nSection/2-1), nSection-1], save=True, position='latex/figures/rotorVelocityTriangle.pdf')
        # computing efficiency
        rotorBlade.computeBladeEfficiency(Va=rotorVaMeanOutlet, lossVec=lossVec)

        # stator study 
        # stator object generation
        print('\n\n-- STATOR STUDY -- # blades {0:d}'.format(nStatorBlades))
        statorBlade = turboBlade.blade(ID=2, turboType='stator', nSection=nSection, inletBladeHeight=b0, outletBladeHeight=b0, inletHubRadius=hubRadius, outletHubRadius=hubRadius, omega=0, nBlade=nStatorBlades)
        # allocate blade shape 
        statorBlade.allocateShape(bladeHeight=b0, AR=ARstator, nBlade=nStatorBlades)
        # blade dimensions allocation 
        # copying data from rotor blade outlet into stator blade inlet
        statorBlade.copySection(blade=rotorBlade, fromSection='outlet', toSection='inlet')
        statorBlade.copySection(blade=rotorBlade, fromSection='outlet', toSection='outlet')

        # setting up velocity distribution at the stato outlet
        VtTarget = 30 # set to this value in order to have convergence of the stator and reduce the stator delta alpha
        def func(r):
            '''
            This is the velocity distribution function.
                VtTarget sets the hub constraint                | -> a second order function is defined 
                statorVtMean = 0 sets the mean line constraint  | -> 

            '''
            a = 4 * VtTarget / b0**2
            Vtnew = a * (r - rMean)**2

            return Vtnew

        # setting up stator outlet kinematics 
        statorBlade.allocateKinetics(rMean=rMean, VtMean=statorVtMeanOutlet, VaMean=statorVaMeanOutlet, omega=0, section='outlet', kind='eqn', func=func)
        # stator blade geometry allocation
        statorBlade.generateGeometry(pos='data/airfoils/naca65.txt', STLname='stator', plot=False, printout=False)
        # computing the best shape 
        lossVec = statorBlade.bladeGenerator(mFlux, clearance=0, NISRE=True, STLname='stator', plot=False, nMaxShape=1, nMaxFlux=100, nMaxS=1)
        # plotting meridional quantities
        statorBlade.printMeridional(save=True, position0='latex/figures/statorEntropyFlow.pdf', position1='latex/figures/statorBetaThermo.pdf')
        # plotting velocity triangle
        statorBlade.velocityTriangles(sectionNumber=[0, int(nSection/2-1), nSection-1], save=True, position='latex/figures/statorVelocityTriangle.pdf')
        # computing efficiency
        statorBlade.computeBladeEfficiency(Va=statorVaMeanOutlet, lossVec=lossVec)

        # .scad file generation 
        nRotorBlades  = rotorBlade.nBlade
        nStatorBlades = statorBlade.nBlade 
        rotorHub      = [rotorBlade.blade[0].chord, rotorBlade.blade[0].camber[0,0], rotorBlade.blade[0].camber[0,1], rotorBlade.blade[0].camber[0,2], rotorBlade.blade[0].camber[-1,0], rotorBlade.blade[0].camber[-1,1], rotorBlade.blade[0].camber[-1,2]]
        statorHub     = [statorBlade.blade[0].chord, statorBlade.blade[0].camber[0,0], statorBlade.blade[0].camber[0,1], statorBlade.blade[0].camber[0,2], statorBlade.blade[0].camber[-1,0], statorBlade.blade[0].camber[-1,1], statorBlade.blade[0].camber[-1,2]]
        b1            = b0 
        bladeGenerator.SCADsaving(nRotorBlades, nStatorBlades, rotorHub, statorHub, rMean, b0, b1, rotorPath='../container/rotor.stl', statorPath='../container/stator.stl', geometryPath='geometry/')

        # computing stage efficiency
        coeff.stageEfficiency(rotorBlade, statorBlade)