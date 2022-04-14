# importing libraries
from turboClass import turboBlade
from turboCoeff import similarity

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
rD = 0.7
# stage mean radius -> radius @ inlet blade midspan
rMean = 0.325
# rotor inlet tangential velocity
Vt0Umean = 0

# Vt0 = (1 - rD - lam/4) * Umean -> lam = 4 * (1 - rD - Vt0/Umean) 
# psi = lam / 2 
lam = (1 - rD - Vt0Umean) * 4 
psiTarget = lam / 2

# generation of mean line properties to be used for the blade assembly 
adimVec, bladeVec, rotationVec, V0vec, V1vec, V2vec, _, _, _, thermo0, _, _, work = similarity.stageProperties(rD, psiTarget, rMean, mFlux, Tt0, Pt0, betaP, T1real=True, printout=True)
 
# values allocation 
nSection     = 50
eta          = adimVec[2]
omega        = rotationVec[1]
b0           = bladeVec[0]
hubRadius    = rMean - b0/2
VaMeanInlet  = V0vec[0]
VtMeanInlet  = V0vec[1] 
VaMeanOutlet = V1vec[0]
VtMeanOutlet = V1vec[1]

# objects generation 
rotorBlade = turboBlade.blade(ID=1, turboType='rotor', nSection=nSection, inletBladeHeight=b0, outletBladeHeight=b0, inletHubRadius=hubRadius, outletHubRadius=hubRadius, omega=omega, nBlade=40)
#statorBlade = turboBlade.blade(ID=1, turboType='stator', nSection=10, inletBladeHeight=b0, outletBladeHeight=b0, inletHubRadius=hubRadius, outletHubRadius=hubRadius, omega=0)

# blade dimensions allocation -> dynamics
rotorBlade.allocateDynamics(rMean=rMean, VtMean=VtMeanInlet, VaMean=VaMeanInlet, omega=omega, section='inlet')
rotorBlade.allocateDynamics(rMean=rMean, VtMean=VtMeanOutlet, VaMean=VaMeanOutlet, omega=omega, section='outlet')
# blade dimensions allocation -> thermodynamics 
rotorBlade.allocateThermodynamics(Tt0=Tt0, Pt0=Pt0, eta=eta)

# blade design through iterative process on radial equilibrium 
rotorBlade.radialEquilibrium(Pt0=Pt0, Tt0=Tt0, mFlux=mFlux ,plot=False)

# rotor blade geometry allocation
rotorBlade.generateShape(pos='data/airfoils/naca65.txt', STLname='cadTest', plot=True, printout=False)

